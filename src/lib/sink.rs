//! The outputs of a run on the worker pool. Workers encode each batch of
//! records into one buffer per destination, and merging appends those buffers
//! in input order. A BGZF output (gzipped FASTX or BAM) is cut into blocks that
//! any worker compresses and its file writes in sequence; uncompressed text and
//! CRAM are written as they merge.

use std::collections::{BTreeMap, HashMap};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex, OnceLock, PoisonError};

use anyhow::{anyhow, Context, Result};
use noodles::bam;

use crate::output::Destination;
use crate::writer::{create_output_file, OutputFormat, OutputRead, OutputWriter};

/// The uncompressed bytes of a full BGZF block.
const BLOCK_SIZE: usize = ::bgzf::BGZF_BLOCK_SIZE;

/// The write buffer of a BGZF output file.
const FILE_BUFFER: usize = 1 << 17;

/// Uncompressed bytes of one destination to compress into one BGZF block, the
/// last of which also ends the file.
pub(crate) struct Block {
    dest: usize,
    seq: u64,
    data: Vec<u8>,
    last: bool,
}

/// A BGZF output written one compressed block at a time, in sequence order
/// whatever order its blocks finish compressing in.
pub(crate) struct BlockFile {
    name: String,
    state: Mutex<InOrder>,
}

/// The file of a [`BlockFile`] and the blocks that finished ahead of their turn.
struct InOrder {
    out: Box<dyn Write + Send>,
    next: u64,
    early: BTreeMap<u64, Vec<u8>>,
}

impl BlockFile {
    fn new(name: String, out: Box<dyn Write + Send>) -> Self {
        Self {
            name,
            state: Mutex::new(InOrder {
                out,
                next: 0,
                early: BTreeMap::new(),
            }),
        }
    }

    /// Write compressed block `seq`, and any blocks after it that wait for it.
    fn write(&self, seq: u64, bytes: Vec<u8>) -> Result<()> {
        let mut state = self.state.lock().unwrap_or_else(PoisonError::into_inner);
        if seq != state.next {
            state.early.insert(seq, bytes);
            return Ok(());
        }
        let mut bytes = Some(bytes);
        while let Some(block) = bytes {
            state
                .out
                .write_all(&block)
                .with_context(|| format!("failed to write output file: {}", self.name))?;
            state.next += 1;
            let next = state.next;
            bytes = state.early.remove(&next);
        }
        Ok(())
    }

    /// Flush the file once its last block is written.
    fn close(&self) -> Result<()> {
        let mut state = self.state.lock().unwrap_or_else(PoisonError::into_inner);
        if !state.early.is_empty() {
            return Err(anyhow!(
                "internal error: blocks of {} were left unwritten",
                self.name
            ));
        }
        state
            .out
            .flush()
            .with_context(|| format!("failed to write output file: {}", self.name))
    }
}

/// Compress `block` with `compressor` and write it to its file in `files`.
pub(crate) fn compress(
    compressor: &mut ::bgzf::Compressor,
    files: &[OnceLock<BlockFile>],
    block: Block,
) -> Result<()> {
    let mut bytes = Vec::new();
    if !block.data.is_empty() {
        compressor
            .compress(&block.data, &mut bytes)
            .map_err(|error| anyhow!("failed to compress a BGZF block: {error}"))?;
    }
    if block.last {
        ::bgzf::Compressor::append_eof(&mut bytes);
    }
    let file = files[block.dest]
        .get()
        .context("internal error: a block has no open file")?;
    file.write(block.seq, bytes)
}

/// A BGZF compressor at `--compression` `level`.
pub(crate) fn compressor(level: u8) -> Result<::bgzf::Compressor> {
    let level = ::bgzf::CompressionLevel::new(level)
        .map_err(|error| anyhow!("invalid --compression {level} for BGZF: {error}"))?;
    Ok(::bgzf::Compressor::new(level))
}

/// What merging writes one destination through.
enum Sink {
    /// Not yet opened: stdout before its first write.
    Closed,
    /// Bytes gathered into blocks for the pool to compress.
    Blocks { pending: Vec<u8>, next: u64 },
    /// Records written as they merge.
    Inline(Box<OutputWriter>),
}

/// Every output of a run, written in input order by merging.
pub(crate) struct Sinks {
    destinations: Arc<[Destination]>,
    files: Arc<[OnceLock<BlockFile>]>,
    sinks: Vec<Sink>,
    compression: u8,
    staged: HashMap<PathBuf, PathBuf>,
}

impl Sinks {
    /// The sinks of `destinations`, each file written at its staged path in
    /// `staged` when it has one.
    pub(crate) fn new(
        destinations: Arc<[Destination]>,
        compression: u8,
        staged: HashMap<PathBuf, PathBuf>,
    ) -> Self {
        let count = destinations.len();
        Self {
            destinations,
            files: (0..count).map(|_| OnceLock::new()).collect(),
            sinks: (0..count).map(|_| Sink::Closed).collect(),
            compression,
            staged,
        }
    }

    /// The BGZF files, for the workers that compress their blocks.
    pub(crate) fn files(&self) -> Arc<[OnceLock<BlockFile>]> {
        Arc::clone(&self.files)
    }

    /// Open every file destination, so each exists even when it receives no
    /// records; stdout opens on its first write.
    pub(crate) fn open_files(&mut self) -> Result<()> {
        for dest in 0..self.sinks.len() {
            if self.destinations[dest].path.is_some() {
                self.open(dest)?;
            }
        }
        Ok(())
    }

    /// Open destination `dest`, writing its header.
    fn open(&mut self, dest: usize) -> Result<()> {
        let Destination {
            path,
            format,
            header,
        } = &self.destinations[dest];
        let file = path
            .as_deref()
            .map(|path| self.staged.get(path).map_or(path, PathBuf::as_path));
        self.sinks[dest] = if format.is_pooled() {
            let (name, out): (String, Box<dyn Write + Send>) = match file {
                Some(file) => {
                    let out = BufWriter::with_capacity(FILE_BUFFER, create_output_file(file)?);
                    let name = path.as_deref().unwrap_or(file).display().to_string();
                    (name, Box::new(out))
                }
                None => {
                    let out = BufWriter::with_capacity(FILE_BUFFER, std::io::stdout());
                    ("<stdout>".to_string(), Box::new(out))
                }
            };
            self.files[dest]
                .set(BlockFile::new(name, out))
                .map_err(|_| anyhow!("internal error: an output opened twice"))?;
            let pending = match format {
                OutputFormat::Bam => {
                    let mut writer = bam::io::Writer::from(Vec::new());
                    writer
                        .write_header(header)
                        .context("failed to write BAM header")?;
                    writer.into_inner()
                }
                _ => Vec::with_capacity(BLOCK_SIZE),
            };
            Sink::Blocks { pending, next: 0 }
        } else {
            let writer =
                OutputWriter::create_with_header(file, *format, self.compression, header.clone())?;
            Sink::Inline(Box::new(writer))
        };
        Ok(())
    }

    /// Append encoded record bytes to destination `dest`, adding the blocks
    /// they fill to `blocks`.
    pub(crate) fn write(
        &mut self,
        dest: usize,
        bytes: &[u8],
        blocks: &mut Vec<Block>,
    ) -> Result<()> {
        if matches!(self.sinks[dest], Sink::Closed) {
            self.open(dest)?;
        }
        match &mut self.sinks[dest] {
            Sink::Blocks { pending, next } => {
                let mut bytes = bytes;
                loop {
                    while pending.len() >= BLOCK_SIZE {
                        blocks.push(full_block(dest, pending, next));
                    }
                    if bytes.is_empty() {
                        return Ok(());
                    }
                    let take = (BLOCK_SIZE - pending.len()).min(bytes.len());
                    pending.extend_from_slice(&bytes[..take]);
                    bytes = &bytes[take..];
                }
            }
            Sink::Inline(writer) => writer.write_encoded(bytes),
            Sink::Closed => Err(anyhow!("internal error: an output did not open")),
        }
    }

    /// Write records to destination `dest`, which encodes them as they merge.
    pub(crate) fn write_records(&mut self, dest: usize, reads: &[OutputRead]) -> Result<()> {
        if matches!(self.sinks[dest], Sink::Closed) {
            self.open(dest)?;
        }
        match &mut self.sinks[dest] {
            Sink::Inline(writer) => writer.write_fragment(reads),
            _ => Err(anyhow!("internal error: records sent to a BGZF output")),
        }
    }

    /// Add the last block of every BGZF output to `blocks`, ending its file.
    pub(crate) fn finish(&mut self, blocks: &mut Vec<Block>) {
        for (dest, sink) in self.sinks.iter_mut().enumerate() {
            if let Sink::Blocks { pending, next } = sink {
                while pending.len() > BLOCK_SIZE {
                    blocks.push(full_block(dest, pending, next));
                }
                blocks.push(Block {
                    dest,
                    seq: *next,
                    data: std::mem::take(pending),
                    last: true,
                });
                *next += 1;
            }
        }
    }

    /// Flush and close every output, finishing each even when one fails, and
    /// return the first error.
    pub(crate) fn close(self) -> Result<()> {
        let mut first_error = None;
        for (dest, sink) in self.sinks.into_iter().enumerate() {
            let closed = match sink {
                Sink::Closed => Ok(()),
                Sink::Blocks { .. } => match self.files[dest].get() {
                    Some(file) => file.close(),
                    None => Ok(()),
                },
                Sink::Inline(writer) => writer.finish(),
            };
            if let Err(error) = closed {
                first_error.get_or_insert(error);
            }
        }
        first_error.map_or(Ok(()), Err)
    }
}

/// Cut the first full block off `pending`.
fn full_block(dest: usize, pending: &mut Vec<u8>, next: &mut u64) -> Block {
    let mut rest = Vec::with_capacity(BLOCK_SIZE);
    rest.extend_from_slice(&pending[BLOCK_SIZE..]);
    pending.truncate(BLOCK_SIZE);
    let block = Block {
        dest,
        seq: *next,
        data: std::mem::replace(pending, rest),
        last: false,
    };
    *next += 1;
    block
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Read;
    use std::path::Path;
    use std::sync::Mutex;

    use crate::input::SniffedFormat;
    use crate::output::provenance_header;

    /// A writer with room for this many more bytes, which fails any write that needs more.
    struct Full(usize);

    impl Write for Full {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            match self.0.checked_sub(buf.len()) {
                Some(room) => {
                    self.0 = room;
                    Ok(buf.len())
                }
                None => Err(std::io::Error::new(
                    std::io::ErrorKind::StorageFull,
                    "disk full",
                )),
            }
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    /// A writer whose bytes can be read back after it is boxed.
    #[derive(Clone, Default)]
    struct Shared(Arc<Mutex<Vec<u8>>>);

    impl Write for Shared {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn test_blocks_are_written_in_sequence_whatever_order_they_arrive() {
        let out = Shared::default();
        let file = BlockFile::new("o/a.fq.gz".to_string(), Box::new(out.clone()));
        for seq in [2, 0, 3, 1] {
            file.write(seq, vec![b'a' + seq as u8]).unwrap();
        }
        file.close().unwrap();
        assert_eq!(*out.0.lock().unwrap(), b"abcd");
    }

    #[test]
    fn test_a_write_error_names_the_file_and_keeps_its_kind() {
        let file = BlockFile::new("o/b.bam".to_string(), Box::new(Full(4)));
        file.write(0, b"fits".to_vec()).unwrap();
        let error = file.write(1, b"does not".to_vec()).unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            "failed to write output file: o/b.bam: disk full"
        );
        let cause = error.downcast_ref::<std::io::Error>().unwrap();
        assert_eq!(cause.kind(), std::io::ErrorKind::StorageFull);
    }

    #[test]
    fn test_an_error_flushing_buffered_bytes_fails_the_close() {
        let buffered = BufWriter::new(Full(0));
        let file = BlockFile::new("o/a.fq.gz".to_string(), Box::new(buffered));
        file.write(0, b"buffered".to_vec()).unwrap();
        let error = file.close().unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            "failed to write output file: o/a.fq.gz: disk full"
        );
    }

    /// Destinations at `paths`, read as FASTQ.
    fn destinations(paths: &[&Path]) -> Arc<[Destination]> {
        let header = provenance_header(None);
        let none = HashMap::new();
        paths
            .iter()
            .map(|path| {
                let path = Some(path.to_path_buf());
                Destination::new(path, &[SniffedFormat::Fastq], &none, &header).unwrap()
            })
            .collect()
    }

    #[test]
    fn test_sinks_fan_bytes_out_to_their_files() {
        let dir = tempfile::tempdir().unwrap();
        let (plain, gzip) = (dir.path().join("a.fq"), dir.path().join("b.fq.gz"));
        let mut sinks = Sinks::new(destinations(&[&plain, &gzip]), 5, HashMap::new());
        sinks.open_files().unwrap();
        let files = sinks.files();
        let mut blocks = Vec::new();
        let big = vec![b'A'; 3 * BLOCK_SIZE + 17];
        sinks.write(0, b"@a\nA\n+\nI\n", &mut blocks).unwrap();
        sinks.write(1, &big[..BLOCK_SIZE - 1], &mut blocks).unwrap();
        sinks.write(1, &big[BLOCK_SIZE - 1..], &mut blocks).unwrap();
        assert_eq!(blocks.len(), 3);
        sinks.finish(&mut blocks);
        let mut compressor = compressor(5).unwrap();
        for block in blocks.into_iter().rev() {
            compress(&mut compressor, &files, block).unwrap();
        }
        sinks.close().unwrap();
        assert_eq!(std::fs::read(&plain).unwrap(), b"@a\nA\n+\nI\n");
        let mut inflated = Vec::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(&gzip).unwrap())
            .read_to_end(&mut inflated)
            .unwrap();
        assert_eq!(inflated, big);
    }
}
