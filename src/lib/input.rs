//! Input reading: identify each input's format by peeking its first bytes,
//! transparently handling a gzip layer, then dispatch to the right reader. The
//! byte-sniffing here uses a single 1 KiB peek to decide FASTA / FASTQ / SAM /
//! BAM / CRAM and whether the stream is gzip-wrapped, without consuming the
//! bytes, so the same reader is handed straight to the matching noodles reader.
//!
//! Above the per-format readers sits the [`FragmentReader`], which opens every
//! `--in` file and yields one [`Fragment`] (a logical record) by pulling a
//! record from each input in lockstep: a fragment's `records[i]` is the record
//! from input file `i`, which is exactly the per-file segment model the matcher
//! and extraction consume. SAM/BAM/CRAM inputs carry their record tags through
//! for the writer to re-emit. A single FASTX input is auto-detected as
//! interleaved, and a single SAM/BAM/CRAM input as paired when its first record
//! has flag `0x1` (two mates per fragment). An input path of `-` reads stdin (at
//! most one). CRAM input decodes against a `--reference` FASTA when given,
//! which must match the CRAM's `@SQ` lines; a CRAM that needs one and has none
//! fails at its first reference-compressed slice.
//!
//! With decode threads to spare, BAM input inflates its BGZF blocks on a worker
//! pool and CRAM input decodes its containers on one. BAM records leave the
//! reader raw ([`PendingRecord::Bam`]) so their decode runs on the matcher
//! workers, and an aligned BAM record passed through to BAM output is written
//! back from those raw bytes without a decode and re-encode, so it keeps the
//! input's tag order and bin.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::sync::LazyLock;

use std::collections::HashMap;

use anyhow::{anyhow, bail, ensure, Context, Result};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::map::{read_group::tag as rgt, ReadGroup};
use noodles::sam::header::record::value::Map;
use noodles::{bam, bgzf, cram, fasta, fastq, sam};
use smallvec::SmallVec;

use crate::cram_raw::Recorder;
use crate::reference::{quietly, read_container, CramDecoder, Reference};

/// A format identified by sniffing the start of an input stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SniffedFormat {
    /// FASTA; first byte `>`.
    Fasta,
    /// FASTQ; first byte `@`, no tab in the first line.
    Fastq,
    /// SAM (text); first byte `@`, first line contains a tab.
    Sam,
    /// BAM; `BAM\x01` after BGZF decompression.
    Bam,
    /// CRAM; `CRAM` magic prefix.
    Cram,
    /// The first bytes match no known signature.
    Unknown,
    /// The input is empty (zero bytes of content) and no format could be
    /// inferred from its file extension; read as zero records rather than
    /// rejected.
    Empty,
}

impl SniffedFormat {
    /// Whether this format carries per-base qualities: FASTA does not,
    /// FASTQ/SAM/BAM/CRAM do. This is the axis inputs must agree on to be
    /// blended; a record assembled across a quality-less and a quality-bearing
    /// segment would have no coherent QUAL.
    pub fn has_qualities(self) -> bool {
        match self {
            SniffedFormat::Fasta | SniffedFormat::Unknown | SniffedFormat::Empty => false,
            SniffedFormat::Fastq
            | SniffedFormat::Sam
            | SniffedFormat::Bam
            | SniffedFormat::Cram => true,
        }
    }
}

/// Classify a slice taken from the start of an input stream (after any gzip
/// layer has been removed). Tolerates short slices and returns `Unknown` when
/// it cannot decide.
pub fn sniff_bytes(bytes: &[u8]) -> SniffedFormat {
    if bytes.starts_with(b"BAM\x01") {
        return SniffedFormat::Bam;
    }
    if bytes.starts_with(b"CRAM") {
        return SniffedFormat::Cram;
    }
    match bytes.first() {
        Some(b'>') => SniffedFormat::Fasta,
        Some(b'@') => sniff_at_prefixed(bytes),
        _ => SniffedFormat::Unknown,
    }
}

/// Input begins with `@`: decide between SAM and FASTQ. A SAM header line is
/// always `@<2 uppercase letters>\t...` (`@HD`/`@SQ`/`@RG`/`@PG`/`@CO`). A FASTQ
/// read-name line is `@<read name>...`; even when it carries tab-separated SAM
/// tags in its comment (the samtools `fastq -T` convention), the token after
/// `@` is a read name, not a two-letter header code. Keying off this start shape
/// (rather than "any tab anywhere in the first line") keeps such a FASTQ from
/// masquerading as SAM, so unmux can read back its own tag-bearing FASTQ output.
fn sniff_at_prefixed(bytes: &[u8]) -> SniffedFormat {
    if bytes.len() >= 4
        && bytes[1].is_ascii_uppercase()
        && bytes[2].is_ascii_uppercase()
        && bytes[3] == b'\t'
    {
        SniffedFormat::Sam
    } else {
        SniffedFormat::Fastq
    }
}

/// Infer a format from a path's extension (a trailing `.gz` is stripped first);
/// the fallback used when an input's content is empty and so cannot be
/// byte-sniffed. Mirrors the extension set of
/// `crate::writer::format_from_extension`. Returns `Unknown` for stdin (`-`) or
/// any unrecognized extension; the caller maps an empty input with no
/// recognized extension to `Empty`.
fn sniff_from_extension(path: &Path) -> SniffedFormat {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    let stem = name.strip_suffix(".gz").unwrap_or(&name);
    if stem.ends_with(".fq") || stem.ends_with(".fastq") {
        SniffedFormat::Fastq
    } else if stem.ends_with(".fa") || stem.ends_with(".fasta") || stem.ends_with(".fna") {
        SniffedFormat::Fasta
    } else if stem.ends_with(".sam") {
        SniffedFormat::Sam
    } else if stem.ends_with(".bam") {
        SniffedFormat::Bam
    } else if stem.ends_with(".cram") {
        SniffedFormat::Cram
    } else {
        SniffedFormat::Unknown
    }
}

/// Open an input source: a file, or stdin when `path` is `-`. The returned
/// reader is a 1 MiB [`BufReader`], large enough to hold any single BGZF block
/// plus headroom for the format peek, and is `Send` so it can move to the
/// read-ahead thread. stdin uses the [`std::io::Stdin`] handle (which is
/// `Send + 'static`, unlike a borrowed `StdinLock`); only one `-` input is
/// allowed (enforced by the grammar), so holding the process-wide stdin is
/// unambiguous.
fn open_source(path: &Path) -> std::io::Result<Box<dyn BufRead + Send>> {
    if path.as_os_str() == "-" {
        Ok(Box::new(BufReader::with_capacity(
            1024 * 1024,
            std::io::stdin(),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(
            1024 * 1024,
            File::open(path)?,
        )))
    }
}

/// Open `path` (or stdin for `-`) and identify the underlying format by peeking
/// its first bytes, transparently decoding a gzip layer for the peek only.
///
/// The returned reader is left at byte 0 (its internal buffer holds the peeked
/// bytes via a non-consuming `fill_buf`, which works over a non-seekable pipe
/// too), so the caller can hand it straight to a noodles text reader, a
/// [`flate2`] decoder, a `bgzf::Reader`, or a `cram::Reader` and start reading
/// from the start. `gzipped` is `true` iff the outer stream starts with `1f
/// 8b`; the caller dispatches by `(format, gzipped)`. A gzip-wrapped CRAM is
/// rejected (it is not a real format).
pub fn sniff_input(path: &Path) -> std::io::Result<(SniffedFormat, bool, Box<dyn BufRead + Send>)> {
    let mut reader = open_source(path)?;

    // Accumulate a prefix big enough to decide the format. A single `fill_buf`
    // can return a SHORT read over a non-seekable pipe (it returns only the
    // first write the kernel delivered, not a full buffer), so a tiny first
    // write would misdetect the format. `take(..).read_to_end` instead loops
    // `read` until the prefix is full or EOF, tolerating short pipe writes. The
    // consumed prefix is replayed in front of the rest of the stream so the
    // caller still reads from byte 0. A regular file's first read already fills
    // the buffer, so this only changes (fixes) the pipe path.
    const SNIFF_PREFIX: u64 = 1024;
    let mut head = Vec::with_capacity(SNIFF_PREFIX as usize);
    (&mut reader).take(SNIFF_PREFIX).read_to_end(&mut head)?;
    let is_gzip = head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b;

    let (format, decoded_empty) = if is_gzip {
        // Decompress a copy of the prefix into scratch space. A decompression
        // failure (truncated member, bad CRC, ...) falls back to `Unknown`
        // rather than propagating: the caller's downstream open surfaces a
        // clearer "failed to read ..." error with the full path.
        let mut scratch = [0u8; 256];
        let mut decoder = flate2::bufread::MultiGzDecoder::new(std::io::Cursor::new(&head));
        let n = decoder.read(&mut scratch).unwrap_or(0);
        let inner = sniff_bytes(&scratch[..n]);
        if inner == SniffedFormat::Cram {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "gzipped CRAM is not a real format; remove the .gz wrapping ({})",
                    path.display()
                ),
            ));
        }
        (inner, n == 0)
    } else {
        (sniff_bytes(&head), head.is_empty())
    };

    // An input with no content cannot be byte-sniffed: fall back to the file
    // extension, so an empty file (e.g. an empty `.fastq.gz` from a run that
    // produced no reads) is read as zero records instead of rejected. A
    // recognized extension picks the concrete format (its reader yields zero
    // records); an unrecognized one (or stdin) becomes `Empty`, also read as
    // zero records. Non-empty unrecognized bytes stay `Unknown`, a real error:
    // the extension is never trusted over real bytes.
    let format = if format == SniffedFormat::Unknown && decoded_empty {
        match sniff_from_extension(path) {
            SniffedFormat::Unknown => SniffedFormat::Empty,
            known => known,
        }
    } else {
        format
    };

    // Replay the consumed prefix in front of the rest of the stream so the
    // caller reads from byte 0.
    let replayed: Box<dyn BufRead + Send> = Box::new(BufReader::with_capacity(
        1024 * 1024,
        std::io::Cursor::new(head).chain(reader),
    ));
    Ok((format, is_gzip, replayed))
}

/// The Phred quality offset. FASTQ stores ASCII Phred+33; the reader decodes it
/// to raw scores so the whole engine works in one quality encoding, and the
/// writer re-applies it for FASTQ/FASTA output.
pub(crate) const PHRED_OFFSET: u8 = 33;

/// One input record for a single segment of a logical record: its name, bases,
/// optional per-base qualities, and (for SAM/BAM/CRAM inputs) the source
/// record's tags carried through to the output.
///
/// Qualities are normalized to raw Phred scores: FASTQ's ASCII Phred+33 is
/// decoded at read time, and SAM/BAM are already raw. The writer re-encodes per
/// output format. FASTA carries no qualities.
#[derive(Debug, Clone, PartialEq)]
pub struct InputRecord {
    /// The read name (the first whitespace-delimited token, the mate-pairing
    /// key).
    pub name: Vec<u8>,
    /// The bases.
    pub bases: Vec<u8>,
    /// The per-base qualities, when the format carries them.
    pub quals: Option<Vec<u8>>,
    /// Carried SAM tags from an alignment input (`None` for FASTX), re-emitted
    /// to the output record.
    pub tags: Option<Data>,
    /// The untouched source record, kept only when the reader passes alignments
    /// through (see [`FragmentReader::keep_alignments`]).
    pub alignment: Option<Alignment>,
}

/// A source alignment kept whole for pass-through output.
#[derive(Debug, Clone, PartialEq)]
pub enum Alignment {
    /// A decoded record, from SAM or CRAM input.
    Record(RecordBuf),
    /// A BAM record's raw bytes, `block_size` prefix included, from BAM input.
    RawBam(Vec<u8>),
}

impl Alignment {
    /// The alignment's `FLAG` field, read without a decode.
    pub fn flags(&self) -> noodles::sam::alignment::record::Flags {
        match self {
            Alignment::Record(record) => record.flags(),
            Alignment::RawBam(raw) => raw_bam_flags(raw),
        }
    }

    /// This alignment as a decoded record, decoding raw BAM bytes on demand.
    pub fn to_record_buf(&self) -> Result<std::borrow::Cow<'_, RecordBuf>> {
        Ok(match self {
            Alignment::Record(record) => std::borrow::Cow::Borrowed(record),
            Alignment::RawBam(raw) => {
                let mut record = RecordBuf::default();
                decode_raw_bam(raw, &mut record)?;
                std::borrow::Cow::Owned(record)
            }
        })
    }

    /// This alignment as an owned decoded record.
    pub fn into_record_buf(self) -> Result<RecordBuf> {
        match self {
            Alignment::Record(record) => Ok(record),
            Alignment::RawBam(raw) => {
                let mut record = RecordBuf::default();
                decode_raw_bam(&raw, &mut record)?;
                Ok(record)
            }
        }
    }
}

/// A placeholder header for BAM record decoding, which never consults it.
static EMPTY_HEADER: LazyLock<sam::Header> = LazyLock::new(sam::Header::default);

/// The size of a BAM record's fixed fields, after its `block_size` prefix.
const BAM_FIXED_FIELDS: usize = 32;

/// Decode raw BAM record bytes (`block_size` prefix included) into `record`.
fn decode_raw_bam(raw: &[u8], record: &mut RecordBuf) -> Result<()> {
    bam::io::Reader::from(raw)
        .read_record_buf(&EMPTY_HEADER, record)
        .context("failed to read BAM record")?;
    Ok(())
}

/// A multithreaded BGZF reader that reports, at end of input, a read error its
/// reader thread hit (a truncated or unreadable frame).
struct CheckedBgzf {
    inner: bgzf::io::MultithreadedReader<Box<dyn BufRead + Send>>,
    finished: bool,
}

impl CheckedBgzf {
    fn new(workers: NonZero<usize>, reader: Box<dyn BufRead + Send>) -> Self {
        Self {
            inner: bgzf::io::MultithreadedReader::with_worker_count(workers, reader),
            finished: false,
        }
    }
}

impl Read for CheckedBgzf {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let available = self.fill_buf()?;
        let n = available.len().min(buf.len());
        buf[..n].copy_from_slice(&available[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for CheckedBgzf {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        if self.finished {
            return Ok(&[]);
        }
        if self.inner.fill_buf()?.is_empty() {
            self.finished = true;
            self.inner.finish()?;
            return Ok(&[]);
        }
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        if !self.finished {
            self.inner.consume(amt);
        }
    }
}

/// Read one raw BAM record, `block_size` prefix included, or `None` at end of
/// input. Rejects a record too short to hold its fixed fields and read name.
fn read_raw_bam(reader: &mut dyn BufRead) -> Result<Option<Vec<u8>>> {
    let mut size = [0u8; 4];
    let mut filled = 0;
    while filled < size.len() {
        match reader.read(&mut size[filled..]) {
            Ok(0) if filled == 0 => return Ok(None),
            Ok(0) => bail!("failed to read BAM record: truncated block size"),
            Ok(n) => filled += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e).context("failed to read BAM record"),
        }
    }
    let block_size = u32::from_le_bytes(size) as usize;
    let mut raw = vec![0u8; 4 + block_size];
    raw[..4].copy_from_slice(&size);
    reader
        .read_exact(&mut raw[4..])
        .context("failed to read BAM record")?;
    let name_end = 4 + BAM_FIXED_FIELDS + raw.get(12).map_or(0, |&n| usize::from(n));
    if block_size < BAM_FIXED_FIELDS || raw.len() < name_end {
        bail!("failed to read BAM record: the record is truncated");
    }
    Ok(Some(raw))
}

/// The `FLAG` field of a raw BAM record.
fn raw_bam_flags(raw: &[u8]) -> noodles::sam::alignment::record::Flags {
    noodles::sam::alignment::record::Flags::from(u16::from_le_bytes([raw[18], raw[19]]))
}

/// The read name of a raw BAM record, empty for a missing (`*`) name.
fn raw_bam_name(raw: &[u8]) -> &[u8] {
    let len = usize::from(raw[12]);
    let name = &raw[4 + BAM_FIXED_FIELDS..4 + BAM_FIXED_FIELDS + len.saturating_sub(1)];
    if name == b"*" {
        &[]
    } else {
        name
    }
}

/// One input record as read, before decoding. A BAM record stays raw so its
/// decode runs on the matcher workers instead of the single read-ahead thread.
#[derive(Debug)]
#[allow(clippy::large_enum_variant)] // boxing would add an allocation per FASTX record
pub enum PendingRecord {
    /// A record its reader already decoded (FASTX, SAM, CRAM).
    Ready(InputRecord),
    /// A raw BAM record, `block_size` prefix included.
    Bam {
        /// The record bytes.
        raw: Vec<u8>,
        /// Whether to keep the record whole for pass-through.
        keep: bool,
    },
}

impl PendingRecord {
    /// The read name, as [`InputRecord::name`] will hold it.
    fn name(&self) -> &[u8] {
        match self {
            PendingRecord::Ready(record) => &record.name,
            PendingRecord::Bam { raw, .. } => raw_bam_name(raw),
        }
    }

    /// Decode this record, reusing `scratch` for a BAM record's fields.
    pub fn decode(self, scratch: &mut RecordBuf) -> Result<InputRecord> {
        match self {
            PendingRecord::Ready(record) => Ok(record),
            PendingRecord::Bam { raw, keep } => {
                decode_raw_bam(&raw, scratch)?;
                let mut input = alignment_to_input(scratch);
                if keep {
                    input.alignment = Some(Alignment::RawBam(raw));
                }
                Ok(input)
            }
        }
    }
}

/// A logical record as read, before its BAM records are decoded: one
/// [`PendingRecord`] per segment, in the order [`Fragment::records`] holds them.
#[derive(Debug)]
pub struct PendingFragment {
    /// The per-segment records.
    pub records: SmallVec<[PendingRecord; 2]>,
}

impl PendingFragment {
    /// Decode every segment into a [`Fragment`], reusing `scratch`.
    pub fn decode(self, scratch: &mut RecordBuf) -> Result<Fragment> {
        let records = self
            .records
            .into_iter()
            .map(|record| record.decode(scratch))
            .collect::<Result<Vec<_>>>()?;
        Ok(Fragment { records })
    }
}

/// A logical record assembled across all inputs: one [`InputRecord`] per
/// 0-based `--in` file, in file-index order. The matcher and extraction see
/// these as the record's per-file segments.
#[derive(Debug, Clone, PartialEq)]
pub struct Fragment {
    /// The per-file records, indexed by input file.
    pub records: Vec<InputRecord>,
}

/// A per-format record reader over one input stream. Holds the concrete noodles
/// reader for the sniffed format; SAM/BAM headers are read once at open and
/// kept for record decoding.
enum FormatReader {
    /// FASTQ (gzip-decoded if needed). The `record` buffer is reused across
    /// reads so noodles refills its internal name/seq/qual Vecs in place rather
    /// than allocating a fresh owned record per call.
    Fastq {
        /// The record reader.
        reader: fastq::io::Reader<Box<dyn BufRead + Send>>,
        /// The reused per-record parse buffer (`read_record` clears and refills
        /// it each call).
        record: fastq::Record,
    },
    /// FASTA (gzip-decoded if needed). The `definition` line buffer is reused
    /// across reads (noodles `read_definition` appends, so it is cleared each
    /// call); the sequence is owned per record.
    Fasta {
        /// The record reader.
        reader: fasta::io::Reader<Box<dyn BufRead + Send>>,
        /// The reused definition-line buffer (cleared before each
        /// `read_definition`).
        definition: String,
    },
    /// Text SAM (gzip-decoded if needed).
    Sam {
        /// The record reader.
        reader: sam::io::Reader<Box<dyn BufRead + Send>>,
        /// The header, needed to decode each record.
        header: Box<sam::Header>,
        /// The reused per-record parse buffer (`read_record_buf` clears and
        /// refills it each call).
        record: RecordBuf,
    },
    /// BAM, its BGZF layer inflated by a worker pool when decode threads allow.
    /// Records are read raw and decoded later by [`PendingRecord::decode`].
    Bam {
        /// The inflated stream, positioned after the header.
        reader: Box<dyn BufRead + Send>,
        /// The header.
        header: Box<sam::Header>,
        /// The reused buffer that skipped secondary and supplementary records
        /// decode into, so a corrupt one fails the read.
        skipped: RecordBuf,
    },
    /// CRAM, read reference-free unless a `--reference` is given. Containers
    /// decode one at a time, or ahead on a worker pool, and their records are
    /// served from `buffered`.
    Cram {
        /// The container source.
        containers: CramContainers,
        /// The header.
        header: Box<sam::Header>,
        /// The decoded records of the current container.
        buffered: std::vec::IntoIter<RecordBuf>,
    },
    /// An empty input (zero records): no underlying reader. `read_next` always
    /// returns `None`. Used when the input has no content and no format could
    /// be inferred from its extension.
    Empty,
}

/// Read the next CRAM container into `container`, keeping its raw bytes; `false` at end of input.
fn read_cram_container(
    reader: &mut cram::io::Reader<Recorder<Box<dyn BufRead + Send>>>,
    container: &mut cram::io::reader::Container,
) -> Result<bool> {
    quietly(|| read_container(reader, container))
        .ok_or_else(|| anyhow!("failed to read CRAM container: it is malformed"))?
        .context("failed to read CRAM container")
}

/// One container's decoded records (`None` at end of input), as a decode worker
/// sends them back.
type ContainerResult = crossbeam_channel::Receiver<Result<Option<Vec<RecordBuf>>>>;

/// The stack for a CRAM decode thread; CRAM 3.1 codecs recurse deeply.
const CRAM_DECODE_STACK: usize = 64 * 1024 * 1024;

/// A CRAM container source: decoded inline, or ahead on a worker pool.
enum CramContainers {
    /// Read and decode each container on the calling thread.
    Inline {
        /// The record reader, keeping each container's raw bytes.
        reader: cram::io::Reader<Recorder<Box<dyn BufRead + Send>>>,
        /// The reusable container buffer.
        container: cram::io::reader::Container,
        /// The header, for decoding.
        header: Box<sam::Header>,
        /// Decodes each container, against the `--reference` FASTA when given.
        decoder: CramDecoder,
        /// Set at the EOF container, whose unread body would misparse as another.
        eof: bool,
    },
    /// Decoded containers in input order: one receiver per container, queued
    /// by a reader thread as it hands the container to a decode worker.
    Pooled {
        /// The per-container results, in input order.
        order: crossbeam_channel::Receiver<ContainerResult>,
    },
}

impl CramContainers {
    /// Start decoding containers from `reader` on `workers` threads.
    fn pooled(
        mut reader: cram::io::Reader<Recorder<Box<dyn BufRead + Send>>>,
        header: sam::Header,
        decoder: &CramDecoder,
        workers: NonZero<usize>,
    ) -> Result<Self> {
        type Job = (
            usize,
            cram::io::reader::Container,
            Vec<u8>,
            crossbeam_channel::Sender<Result<Option<Vec<RecordBuf>>>>,
        );
        let header = std::sync::Arc::new(header);
        let (order_tx, order) = crossbeam_channel::bounded(2 * workers.get());
        let (job_tx, job_rx) = crossbeam_channel::bounded::<Job>(workers.get());
        for _ in 0..workers.get() {
            let job_rx = job_rx.clone();
            let header = std::sync::Arc::clone(&header);
            let decoder = decoder.clone();
            std::thread::Builder::new()
                .stack_size(CRAM_DECODE_STACK)
                .spawn(move || {
                    for (order, container, raw, result_tx) in job_rx {
                        let decoded = decoder.decode_in_order(order, &container, &header, &raw);
                        let _ = result_tx.send(decoded.map(Some));
                    }
                })
                .context("failed to start a CRAM decode thread")?;
        }
        let decoder = decoder.clone();
        std::thread::spawn(move || loop {
            let (result_tx, result_rx) = crossbeam_channel::bounded(1);
            if order_tx.send(result_rx).is_err() {
                return;
            }
            let mut container = cram::io::reader::Container::default();
            match read_cram_container(&mut reader, &mut container) {
                Ok(false) => {
                    let _ = result_tx.send(Ok(None));
                    return;
                }
                Ok(true) => {
                    let job = (
                        decoder.ticket(),
                        container,
                        reader.get_mut().take(),
                        result_tx,
                    );
                    if job_tx.send(job).is_err() {
                        return;
                    }
                }
                Err(e) => {
                    let _ = result_tx.send(Err(e));
                    return;
                }
            }
        });
        Ok(CramContainers::Pooled { order })
    }

    /// The next container's records, or `None` at end of input.
    fn next(&mut self) -> Result<Option<Vec<RecordBuf>>> {
        match self {
            CramContainers::Inline {
                reader,
                container,
                header,
                decoder,
                eof,
            } => {
                if *eof {
                    return Ok(None);
                }
                if !read_cram_container(reader, container)? {
                    *eof = true;
                    return Ok(None);
                }
                decoder
                    .decode(container, header, reader.get_ref().bytes())
                    .map(Some)
            }
            CramContainers::Pooled { order } => {
                let Ok(result) = order.recv() else {
                    return Ok(None);
                };
                result
                    .recv()
                    .map_err(|_| anyhow!("a CRAM decode thread stopped unexpectedly"))?
            }
        }
    }
}

impl FormatReader {
    /// Read the next record from this stream with its `FLAG` field (empty for
    /// FASTX), or `None` at end of input. With `keep` set, secondary and
    /// supplementary records are returned too and each record carries its
    /// untouched source alignment.
    fn read_next(&mut self, keep: bool) -> Result<Option<(InputRecord, Flags)>> {
        match self {
            FormatReader::Fastq { reader, record } => {
                if reader
                    .read_record(record)
                    .context("failed to read FASTQ record")?
                    == 0
                {
                    return Ok(None);
                }
                // A FASTQ record's sequence and quality strings must be the
                // same length. noodles does not enforce this, so a malformed
                // record would otherwise carry the mismatch downstream and
                // panic on an out-of-bounds quality slice during extraction;
                // reject it here instead.
                let (n_bases, n_quals) = (record.sequence().len(), record.quality_scores().len());
                if n_bases != n_quals {
                    bail!(
                        "malformed FASTQ record `{}`: the sequence has {n_bases} base(s) but the quality string has {n_quals} character(s)",
                        String::from_utf8_lossy(&first_token(record.name()))
                    );
                }
                let input = InputRecord {
                    name: first_token(record.name()),
                    bases: record.sequence().to_vec(),
                    quals: Some(
                        record
                            .quality_scores()
                            .iter()
                            .map(|&q| q.saturating_sub(PHRED_OFFSET))
                            .collect(),
                    ),
                    tags: parse_comment_tags(record.description()),
                    alignment: None,
                };
                Ok(Some((input, Flags::empty())))
            }
            FormatReader::Fasta { reader, definition } => {
                // noodles `read_definition` appends to the buffer, so clear the
                // reused line first.
                definition.clear();
                if reader
                    .read_definition(definition)
                    .context("failed to read FASTA definition")?
                    == 0
                {
                    return Ok(None);
                }
                let mut bases = Vec::new();
                reader
                    .read_sequence(&mut bases)
                    .context("failed to read FASTA sequence")?;
                let body = definition.trim_start_matches('>');
                let name = body
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .as_bytes()
                    .to_vec();
                // Lift any SAM tags carried in the description (the definition
                // remainder after the read name), mirroring the FASTQ path.
                let comment = body
                    .trim_start()
                    .split_once([' ', '\t'])
                    .map_or("", |(_, rest)| rest);
                let input = InputRecord {
                    name,
                    bases,
                    quals: None,
                    tags: parse_comment_tags(comment.as_bytes()),
                    alignment: None,
                };
                Ok(Some((input, Flags::empty())))
            }
            FormatReader::Sam {
                reader,
                header,
                record,
            } => loop {
                if reader
                    .read_record_buf(header, record)
                    .context("failed to read SAM record")?
                    == 0
                {
                    return Ok(None);
                }
                let flags = record.flags();
                if keep {
                    return Ok(Some((kept_alignment(std::mem::take(record)), flags)));
                }
                // Secondary/supplementary alignments are not distinct reads;
                // skip to the next primary record so the per-file segment count
                // stays one-per-template.
                if flags.is_secondary() || flags.is_supplementary() {
                    continue;
                }
                return Ok(Some((alignment_to_input(record), flags)));
            },
            FormatReader::Bam { .. } => match self.read_pending(keep)? {
                Some((pending, flags)) => {
                    Ok(Some((pending.decode(&mut RecordBuf::default())?, flags)))
                }
                None => Ok(None),
            },
            FormatReader::Cram {
                containers,
                buffered,
                ..
            } => loop {
                if let Some(record) = buffered.next() {
                    let flags = record.flags();
                    if keep {
                        return Ok(Some((kept_alignment(record), flags)));
                    }
                    if flags.is_secondary() || flags.is_supplementary() {
                        continue;
                    }
                    return Ok(Some((alignment_to_input(&record), flags)));
                }
                match containers.next()? {
                    None => return Ok(None),
                    Some(records) => *buffered = records.into_iter(),
                }
            },
            FormatReader::Empty => Ok(None),
        }
    }

    /// [`Self::read_next`], leaving a BAM record raw for a worker to decode.
    fn read_pending(&mut self, keep: bool) -> Result<Option<(PendingRecord, Flags)>> {
        let FormatReader::Bam {
            reader, skipped, ..
        } = self
        else {
            let read = self.read_next(keep)?;
            return Ok(read.map(|(record, flags)| (PendingRecord::Ready(record), flags)));
        };
        loop {
            let Some(raw) = read_raw_bam(reader.as_mut())? else {
                return Ok(None);
            };
            let flags = raw_bam_flags(&raw);
            if keep || !(flags.is_secondary() || flags.is_supplementary()) {
                return Ok(Some((PendingRecord::Bam { raw, keep }, flags)));
            }
            decode_raw_bam(&raw, skipped)?;
        }
    }

    /// This input's SAM header, or `None` for a FASTX/empty input (which
    /// carries no header).
    fn header(&self) -> Option<&sam::Header> {
        match self {
            FormatReader::Sam { header, .. }
            | FormatReader::Bam { header, .. }
            | FormatReader::Cram { header, .. } => Some(header),
            _ => None,
        }
    }
}

/// The first record of an alignment input, kept whole and read on the CRAM decode stack.
fn probe_first(reader: &mut FormatReader) -> Result<Option<(InputRecord, Flags)>> {
    std::thread::scope(|scope| {
        std::thread::Builder::new()
            .stack_size(CRAM_DECODE_STACK)
            .spawn_scoped(scope, || reader.read_next(true))
            .context("failed to start a thread to read the first record")?
            .join()
            .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
    })
}

/// The first whitespace-delimited token of a name. FASTX names may carry a
/// trailing description; the token before the first space is the mate-pairing
/// key.
fn first_token(name: &[u8]) -> Vec<u8> {
    name.split(|&b| b == b' ' || b == b'\t')
        .next()
        .unwrap_or(name)
        .to_vec()
}

/// Lift well-formed SAM auxiliary tags out of a FASTX read-name comment, the
/// inverse of the writer's `render_tags`. Mirrors `samtools import -T`: the
/// comment is split on tabs and each field that parses as a SAM aux field
/// (`XX:T:VALUE`) is lifted, while a field that is not a well-formed tag (a bare
/// UMI, an Illumina CASAVA string, free text) is skipped and the scan
/// continues. Returns `None` when no tag is lifted, leaving the record tag-free.
fn parse_comment_tags(comment: &[u8]) -> Option<Data> {
    let mut data = Data::default();
    for field in comment.split(|&b| b == b'\t') {
        if let Some((tag, value)) = parse_aux_field(field) {
            data.insert(tag, value);
        }
    }
    (!data.is_empty()).then_some(data)
}

/// Parse one `XX:T:VALUE` SAM aux field, supporting the `A`, `i`, `f`, and `Z`
/// types the writer emits. Returns `None` for any field that does not match the
/// `<alpha><alphanum>:<type>:<value>` shape or whose value is invalid for its
/// type, so a non-tag comment field is skipped rather than mangled. A `Z` value
/// may contain spaces (only tabs separate fields, and the comment is already
/// split on tabs).
fn parse_aux_field(field: &[u8]) -> Option<(Tag, Value)> {
    if field.len() < 5 || field[2] != b':' || field[4] != b':' {
        return None;
    }
    let (b0, b1, ty) = (field[0], field[1], field[3]);
    // SAM tag names are `[A-Za-z][A-Za-z0-9]`.
    if !b0.is_ascii_alphabetic() || !b1.is_ascii_alphanumeric() {
        return None;
    }
    let raw = &field[5..];
    let value = match ty {
        b'A' if raw.len() == 1 => Value::Character(raw[0]),
        b'i' => Value::Int32(std::str::from_utf8(raw).ok()?.parse().ok()?),
        b'f' => Value::Float(std::str::from_utf8(raw).ok()?.parse().ok()?),
        b'Z' => Value::String(raw.to_vec().into()),
        _ => return None,
    };
    Some((Tag::new(b0, b1), value))
}

/// SAM aux tags that describe an alignment to a reference (scores, edit
/// distances, reference coordinates/CIGAR, mate-alignment info, multi-hit
/// bookkeeping). They are stale once a record is reverted to unmapped, so they
/// are stripped on input alongside the alignment flags/positions. The set is
/// the union of `samtools reset`, Picard RevertSam, and the tags bwa mem /
/// bwa-mem2 emit, plus `ms` (`samtools fixmate -m`). Read-group, UMI, barcode,
/// and cell tags are deliberately kept.
pub(crate) const ALIGNMENT_TAGS: &[[u8; 2]] = &[
    *b"AS", *b"NM", *b"MD", *b"SA", *b"MC", *b"MQ", *b"ms", *b"XS", *b"XA", *b"XB", *b"XR", *b"pa",
    *b"UQ", *b"CC", *b"CP", *b"CG", *b"H0", *b"H1", *b"H2", *b"HI", *b"IH", *b"NH", *b"TS",
];

/// Convert an alignment record into an [`InputRecord`], carrying its
/// non-alignment tags through. Alignment-derived tags ([`ALIGNMENT_TAGS`]) are
/// stripped, since the output is unmapped. A read with bases and QUAL `*` has
/// no qualities, while an empty read has an empty quality string.
///
/// A reverse-strand record (flag `0x10`) stores `SEQ` as the reverse complement
/// of the sequenced read and `QUAL` in reversed order; this restores the
/// original read orientation so matching and extraction see the sequenced
/// bases.
fn alignment_to_input(record: &RecordBuf) -> InputRecord {
    let raw_quals = record.quality_scores().as_ref();
    let mut bases = record.sequence().as_ref().to_vec();
    let mut quals = (!raw_quals.is_empty() || bases.is_empty()).then(|| raw_quals.to_vec());
    if record.flags().is_reverse_complemented() {
        bases = bases
            .iter()
            .rev()
            .map(|base| crate::iupac::complement(*base))
            .collect();
        if let Some(quals) = quals.as_mut() {
            quals.reverse();
        }
    }
    let tags: Data = record
        .data()
        .iter()
        .filter(|(tag, _)| {
            let bytes: &[u8; 2] = tag.as_ref();
            !ALIGNMENT_TAGS.contains(bytes)
        })
        .map(|(tag, value)| (tag, value.clone()))
        .collect();
    InputRecord {
        name: record.name().map(|n| n.to_vec()).unwrap_or_default(),
        bases,
        quals,
        tags: Some(tags),
        alignment: None,
    }
}

/// An [`InputRecord`] that also keeps its untouched source alignment, for a run
/// that passes aligned records through unchanged.
fn kept_alignment(record: RecordBuf) -> InputRecord {
    let mut input = alignment_to_input(&record);
    input.alignment = Some(Alignment::Record(record));
    input
}

/// Wrap a sniffed text stream for reading, transparently decoding a plain-gzip
/// layer. BAM is BGZF and decodes inside its own reader, so it never passes
/// through here.
fn text_bufread(reader: Box<dyn BufRead + Send>, gzipped: bool) -> Box<dyn BufRead + Send> {
    if gzipped {
        // Size the decoded-stream buffer the noodles parser refills from to 1
        // MiB (matching the raw file BufReader), not the 8 KiB stdlib-class
        // default: the gzipped path is the hot one, and a small buffer here
        // forces many short refills of the parser over the decode stream.
        Box::new(BufReader::with_capacity(
            1024 * 1024,
            flate2::bufread::MultiGzDecoder::new(reader),
        ))
    } else {
        reader
    }
}

/// Split `workers` across the inputs among `formats` that are `format`, the
/// earlier inputs taking any remainder; other inputs get none.
pub fn decode_shares(
    formats: &[SniffedFormat],
    format: SniffedFormat,
    workers: usize,
) -> Vec<usize> {
    let inputs = formats.iter().filter(|&&f| f == format).count();
    let mut nth = 0;
    formats
        .iter()
        .map(|&f| {
            if f != format {
                return 0;
            }
            let share = workers / inputs + usize::from(nth < workers % inputs);
            nth += 1;
            share
        })
        .collect()
}

/// Open one sniffed input with the matching noodles reader. With
/// `decode_threads`, BAM inflates and CRAM decodes on that many workers.
fn open_format_reader(
    path: &Path,
    (format, gzipped, reader): (SniffedFormat, bool, Box<dyn BufRead + Send>),
    reference: Option<&Reference>,
    decode_threads: Option<NonZero<usize>>,
) -> Result<FormatReader> {
    let format_reader = match format {
        SniffedFormat::Fastq => FormatReader::Fastq {
            reader: fastq::io::Reader::new(text_bufread(reader, gzipped)),
            record: fastq::Record::default(),
        },
        SniffedFormat::Fasta => FormatReader::Fasta {
            reader: fasta::io::Reader::new(text_bufread(reader, gzipped)),
            definition: String::new(),
        },
        SniffedFormat::Sam => {
            let mut reader = sam::io::Reader::new(text_bufread(reader, gzipped));
            let header = reader.read_header().context("failed to read SAM header")?;
            FormatReader::Sam {
                reader,
                header: Box::new(header),
                record: RecordBuf::default(),
            }
        }
        SniffedFormat::Bam => {
            let inflated: Box<dyn BufRead + Send> = match decode_threads {
                Some(workers) => Box::new(CheckedBgzf::new(workers, reader)),
                None => Box::new(bgzf::io::Reader::new(reader)),
            };
            let mut reader = bam::io::Reader::from(inflated);
            let header = reader.read_header().context("failed to read BAM header")?;
            FormatReader::Bam {
                reader: reader.into_inner(),
                header: Box::new(header),
                skipped: RecordBuf::default(),
            }
        }
        SniffedFormat::Cram => {
            // CRAM takes a plain Read (no BGZF/gzip wrapper; gzipped CRAM is
            // already rejected).
            let mut reader = cram::io::Reader::new(Recorder::new(reader));
            let header = quietly(|| reader.read_header())
                .ok_or_else(|| anyhow!("failed to read CRAM header: it is malformed"))?
                .context("failed to read CRAM header")?;
            let decoder = CramDecoder::new(reference.cloned(), path, &header);
            let containers = match decode_threads {
                Some(workers) => CramContainers::pooled(reader, header.clone(), &decoder, workers)?,
                None => CramContainers::Inline {
                    reader,
                    container: cram::io::reader::Container::default(),
                    header: Box::new(header.clone()),
                    decoder,
                    eof: false,
                },
            };
            FormatReader::Cram {
                containers,
                header: Box::new(header),
                buffered: Vec::new().into_iter(),
            }
        }
        SniffedFormat::Empty => FormatReader::Empty,
        SniffedFormat::Unknown => {
            bail!("could not identify the input format of {}", path.display())
        }
    };
    Ok(format_reader)
}

/// Reject inputs that mix formats with and without per-base qualities: FASTA
/// (no qualities) cannot be blended with quality-bearing formats
/// (FASTQ/SAM/BAM/CRAM), since a record assembled across a quality-less and a
/// quality-bearing segment has no coherent QUAL. Quality-bearing formats blend
/// freely (e.g. FASTQ with BAM).
///
/// `Empty` inputs are skipped: an empty input carries no records, so it has no
/// quality axis to conflict with. A genuine empty-vs-non-empty mismatch is not
/// this function's to diagnose; it is caught at read time by the record-count
/// check, which reports the accurate cause ("unequal record counts") rather
/// than a misleading "FASTA cannot be blended" quality error. (Skipping `Empty`
/// is also what keeps an all-empty run readable when one empty input lacks a
/// format-bearing extension.)
fn ensure_compatible_quality(formats: &[SniffedFormat]) -> Result<()> {
    let mut anchor: Option<(usize, SniffedFormat)> = None;
    for (i, &format) in formats.iter().enumerate() {
        if format == SniffedFormat::Empty {
            continue;
        }
        match anchor {
            None => anchor = Some((i, format)),
            Some((j, first)) if first.has_qualities() != format.has_qualities() => {
                bail!(
                    "input files mix formats with and without per-base qualities: file {j} is {first:?} but file {i} is {format:?}; FASTA (no qualities) cannot be blended with quality-bearing formats (FASTQ/SAM/BAM/CRAM)"
                );
            }
            Some(_) => {}
        }
    }
    Ok(())
}

/// One input read group, lifted from an `@RG` header line.
#[derive(Debug, Clone)]
pub struct InputReadGroup {
    /// The RG id (the `RG:Z` value records carry).
    pub id: Vec<u8>,
    /// The `SM` subfield, if present.
    pub sample: Option<String>,
    /// The `LB` subfield, if present.
    pub library: Option<String>,
    /// The original `@RG` map, re-emitted verbatim into output headers.
    pub map: Map<ReadGroup>,
}

/// The read groups (and `@SQ` presence) discovered across all inputs.
#[derive(Debug, Clone, Default)]
pub struct ReadGroupInfo {
    /// The read groups in first-seen order across inputs (ids are unique).
    pub read_groups: Vec<InputReadGroup>,
    /// Whether any input header carries `@SQ` reference sequences.
    pub has_sq: bool,
}

/// Reads logical records across all inputs. With several `--in` files, one
/// record is pulled from each in lockstep; a single FASTX input is
/// auto-detected as interleaved and a single SAM/BAM/CRAM input as paired (two
/// mates per fragment) unless `--per-record` opts out.
pub struct FragmentReader {
    readers: Vec<FormatReader>,
    /// The sniffed format of each input, in `--in` order (blend-compatible by
    /// construction).
    formats: Vec<SniffedFormat>,
    /// A single FASTX input detected as interleaved: pull two mates per
    /// fragment from `readers[0]`.
    interleaved: bool,
    /// Whether a single alignment input is paired, or `None` when pairs are not probed.
    paired: Option<bool>,
    /// Records buffered from `readers[0]` during pair detection, served before
    /// fresh reads. A SAM/BAM/CRAM record is held whole, as
    /// [`Self::keep_alignments`] may still be called.
    buffered: Vec<InputRecord>,
    /// A read error from pair detection, held so the checks before streaming report first.
    probe_error: Option<anyhow::Error>,
    /// Pass SAM/BAM/CRAM records through whole (see [`Self::keep_alignments`]).
    keep_alignments: bool,
    /// The BAM inflate and CRAM decode workers started across all inputs.
    decode_threads: usize,
}

/// What to do about a SAM/BAM/CRAM input whose mates are not next to each other.
const REGROUP_HINT: &str = "sort the input by read name (e.g. `samtools sort -n`), or pass \
                            --per-record to read each record on its own";

impl FragmentReader {
    /// Open every input by sniffing and dispatching, preserving 0-based `--in`
    /// order. With a single FASTX input and `per_record` false, the stream is
    /// probed for interleaving: a leading `/1` mate, or two adjacent records
    /// that share a mate-pairing key, marks it interleaved. A single
    /// SAM/BAM/CRAM input is probed for a first record with flag `0x1`, which
    /// marks it paired.
    pub fn open(paths: &[PathBuf], per_record: bool) -> Result<Self> {
        Self::open_with_reference(paths, per_record, None)
    }

    /// [`Self::open`], decoding CRAM input against `reference` when given; see [`Self::check_reference`].
    pub fn open_with_reference(
        paths: &[PathBuf],
        per_record: bool,
        reference: Option<&Reference>,
    ) -> Result<Self> {
        Self::open_with_threads(paths, per_record, reference, |formats| {
            vec![0; formats.len()]
        })
    }

    /// [`Self::open_with_reference`], inflating BAM and decoding CRAM input on
    /// the workers `decode_threads` allots each input from the sniffed formats;
    /// an input left without one decodes inline.
    pub fn open_with_threads(
        paths: &[PathBuf],
        per_record: bool,
        reference: Option<&Reference>,
        decode_threads: impl FnOnce(&[SniffedFormat]) -> Vec<usize>,
    ) -> Result<Self> {
        let open_error =
            |i: usize, path: &Path| format!("failed to open input file {i} ({})", path.display());
        let sniffed = paths
            .iter()
            .enumerate()
            .map(|(i, path)| {
                sniff_input(path)
                    .with_context(|| format!("failed to read: {}", path.display()))
                    .with_context(|| open_error(i, path))
            })
            .collect::<Result<Vec<_>>>()?;
        let formats: Vec<SniffedFormat> = sniffed.iter().map(|(format, _, _)| *format).collect();
        let shares = decode_threads(&formats);
        ensure!(
            shares.len() == formats.len(),
            "decode workers were allotted for {} inputs, not {}",
            shares.len(),
            formats.len()
        );
        let mut readers = Vec::with_capacity(paths.len());
        for (i, (path, sniffed)) in paths.iter().zip(sniffed).enumerate() {
            let reader = open_format_reader(path, sniffed, reference, NonZero::new(shares[i]))
                .with_context(|| open_error(i, path))?;
            readers.push(reader);
        }
        ensure_compatible_quality(&formats)?;

        let mut buffered = Vec::new();
        let mut interleaved = false;
        let mut paired = None;
        let mut probe_error = None;
        // Only a single input is probed for pairs; several inputs are zipped.
        if !per_record
            && readers.len() == 1
            && matches!(
                readers[0],
                FormatReader::Fastq { .. } | FormatReader::Fasta { .. }
            )
        {
            if let Some((first, _)) = readers[0].read_next(false)? {
                if first.name.ends_with(b"/2") {
                    bail!(
                        "interleaved input starts with a `/2` mate (mates out of order): {}",
                        String::from_utf8_lossy(&first.name)
                    );
                }
                let first_is_r1 = first.name.ends_with(b"/1");
                buffered.push(first);
                if first_is_r1 {
                    interleaved = true;
                } else if let Some((second, _)) = readers[0].read_next(false)? {
                    interleaved = mate_key(&second.name) == mate_key(&buffered[0].name);
                    buffered.push(second);
                }
            }
        } else if !per_record && readers.len() == 1 && readers[0].header().is_some() {
            match probe_first(&mut readers[0]) {
                Ok(Some((first, flags))) => {
                    paired = Some(flags.is_segmented());
                    buffered.push(first);
                }
                Ok(None) => {}
                Err(error) => probe_error = Some(error),
            }
        }

        Ok(Self {
            readers,
            formats,
            interleaved,
            paired,
            buffered,
            probe_error,
            keep_alignments: false,
            decode_threads: shares.iter().sum(),
        })
    }

    /// Pass SAM/BAM/CRAM records through whole: every record, secondary and
    /// supplementary included, is yielded on its own with its untouched source
    /// alignment in [`InputRecord::alignment`], and mates are not paired.
    pub fn keep_alignments(&mut self) {
        self.keep_alignments = true;
    }

    /// The first input's SAM header, or `None` for a FASTX or empty input.
    pub fn alignment_header(&self) -> Option<&sam::Header> {
        self.readers.first().and_then(FormatReader::header)
    }

    /// The sniffed format of each input, in `--in` order. The writer uses this
    /// to mirror the input format to stdout, defaulting a multi-format blend to
    /// uncompressed BAM.
    pub fn formats(&self) -> &[SniffedFormat] {
        &self.formats
    }

    /// The BAM inflate and CRAM decode workers started across all inputs.
    pub fn decode_threads(&self) -> usize {
        self.decode_threads
    }

    /// The read groups declared by the input headers, in first-seen order, plus
    /// whether any input is aligned (`@SQ` present). FASTX inputs contribute
    /// nothing. A repeated RG id whose `SM`/`LB` disagree across inputs errors.
    pub fn input_read_groups(&self) -> Result<ReadGroupInfo> {
        let mut info = ReadGroupInfo::default();
        let mut seen: HashMap<Vec<u8>, usize> = HashMap::new();
        for reader in &self.readers {
            let Some(header) = reader.header() else {
                continue;
            };
            if !header.reference_sequences().is_empty() {
                info.has_sq = true;
            }
            for (id, rg) in header.read_groups() {
                let id: Vec<u8> = AsRef::<[u8]>::as_ref(id).to_vec();
                let sample = rg.other_fields().get(&rgt::SAMPLE).map(|v| v.to_string());
                let library = rg.other_fields().get(&rgt::LIBRARY).map(|v| v.to_string());
                match seen.get(&id) {
                    Some(&i) => {
                        let existing = &info.read_groups[i];
                        if existing.sample != sample || existing.library != library {
                            bail!(
                                "read group `{}` has conflicting SM/LB across inputs",
                                String::from_utf8_lossy(&id)
                            );
                        }
                    }
                    None => {
                        seen.insert(id.clone(), info.read_groups.len());
                        info.read_groups.push(InputReadGroup {
                            id,
                            sample,
                            library,
                            map: rg.clone(),
                        });
                    }
                }
            }
        }
        Ok(info)
    }

    /// Check each CRAM input's `@SQ` lines against `reference`, `paths` being the inputs in `--in` order.
    pub fn check_reference(&self, paths: &[PathBuf], reference: &Reference) -> Result<()> {
        for (i, (reader, path)) in self.readers.iter().zip(paths).enumerate() {
            if let FormatReader::Cram { header, .. } = reader {
                reference.check_header(header, path).with_context(|| {
                    format!(
                        "failed to check input file {i} ({}) against --reference",
                        path.display()
                    )
                })?;
            }
        }
        Ok(())
    }

    /// The number of records every fragment carries: two for a single
    /// auto-detected interleaved or paired input, otherwise one per input file.
    /// This fixes the pass-through body count (one body per record), so output
    /// destinations can be enumerated before streaming.
    pub fn fragment_width(&self) -> usize {
        let paired = self.paired == Some(true) && !self.keep_alignments;
        if self.interleaved || paired {
            2
        } else {
            self.readers.len()
        }
    }

    /// Read the next logical record. A single interleaved or paired input
    /// yields a two-mate fragment; a single plain input yields a one-record
    /// fragment; multiple inputs are zipped one record each. Returns `None`
    /// once every input is exhausted together. A file that ends before the
    /// others, input segments whose read names disagree, a paired record not
    /// followed by its mate, or a mix of paired and unpaired records is a
    /// fail-fast error (each would silently mispair records).
    pub fn next_fragment(&mut self) -> Result<Option<Fragment>> {
        match self.next_pending_fragment()? {
            Some(pending) => pending.decode(&mut RecordBuf::default()).map(Some),
            None => Ok(None),
        }
    }

    /// [`Self::next_fragment`], leaving BAM records raw so the caller can decode
    /// them on another thread with [`PendingFragment::decode`].
    pub fn next_pending_fragment(&mut self) -> Result<Option<PendingFragment>> {
        if let Some(error) = self.probe_error.take() {
            return Err(error);
        }
        if self.readers.len() == 1 {
            let (first, flags) = match self.next_from_first()? {
                Some(read) => read,
                None => return Ok(None),
            };
            if let Some(paired) = self.paired.filter(|_| !self.keep_alignments) {
                return self.flagged_fragment(first, flags, paired).map(Some);
            }
            if !self.interleaved {
                return Ok(Some(PendingFragment {
                    records: smallvec::smallvec![first],
                }));
            }
            let (second, _) = self.next_from_first()?.ok_or_else(|| {
                anyhow!("interleaved input ended on an unpaired mate (odd number of records)")
            })?;
            check_names_concord(&[first.name(), second.name()])?;
            return Ok(Some(PendingFragment {
                records: smallvec::smallvec![first, second],
            }));
        }

        let mut maybe: SmallVec<[Option<PendingRecord>; 4]> =
            SmallVec::with_capacity(self.readers.len());
        for reader in &mut self.readers {
            let read = reader.read_pending(self.keep_alignments)?;
            maybe.push(read.map(|(record, _)| record));
        }
        let present = maybe.iter().filter(|record| record.is_some()).count();
        if present == 0 {
            return Ok(None);
        }
        if present != maybe.len() {
            bail!("input files have unequal record counts (one input ended before the others)");
        }
        let records: SmallVec<[PendingRecord; 2]> = maybe.into_iter().map(Option::unwrap).collect();
        check_names_concord(
            &records
                .iter()
                .map(PendingRecord::name)
                .collect::<SmallVec<[&[u8]; 4]>>(),
        )?;
        Ok(Some(PendingFragment { records }))
    }

    /// The fragment that `first` starts in a single alignment input that is `paired` or not.
    fn flagged_fragment(
        &mut self,
        first: PendingRecord,
        flags: Flags,
        paired: bool,
    ) -> Result<PendingFragment> {
        let name = |record: &PendingRecord| String::from_utf8_lossy(record.name()).into_owned();
        if flags.is_segmented() != paired {
            bail!(
                "input mixes paired and unpaired records (`{}` is {}, unlike the first record); \
                 pass --per-record to read each record on its own",
                name(&first),
                if paired { "unpaired" } else { "paired" }
            );
        }
        if !paired {
            return Ok(PendingFragment {
                records: smallvec::smallvec![first],
            });
        }
        let Some((second, mate_flags)) = self.next_from_first()? else {
            bail!(
                "paired record `{}` ends the input without its mate; {REGROUP_HINT}",
                name(&first)
            );
        };
        if mate_key(second.name()) != mate_key(first.name()) {
            bail!(
                "paired record `{}` is followed by `{}`, not its mate; {REGROUP_HINT}",
                name(&first),
                name(&second)
            );
        }
        let segment = Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::LAST_SEGMENT;
        let is_first = |flags: Flags| flags & segment == Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        let is_last = |flags: Flags| flags & segment == Flags::SEGMENTED | Flags::LAST_SEGMENT;
        let records = if is_first(flags) && is_last(mate_flags) {
            smallvec::smallvec![first, second]
        } else if is_last(flags) && is_first(mate_flags) {
            smallvec::smallvec![second, first]
        } else {
            bail!(
                "adjacent records `{}` have flags {} and {}, not those of one first mate (0x40) \
                 and one last mate (0x80); pass --per-record to read each record on its own",
                name(&first),
                u16::from(flags),
                u16::from(mate_flags)
            );
        };
        Ok(PendingFragment { records })
    }

    /// The next record from the first input with its `FLAG` field, draining
    /// any records buffered during pair detection before reading fresh ones.
    fn next_from_first(&mut self) -> Result<Option<(PendingRecord, Flags)>> {
        while !self.buffered.is_empty() {
            let mut record = self.buffered.remove(0);
            let alignment = record.alignment.as_ref();
            let flags = alignment.map(Alignment::flags).unwrap_or_default();
            if !self.keep_alignments {
                if flags.is_secondary() || flags.is_supplementary() {
                    continue;
                }
                record.alignment = None;
            }
            return Ok(Some((PendingRecord::Ready(record), flags)));
        }
        self.readers[0].read_pending(self.keep_alignments)
    }
}

/// A name's mate-pairing key: the name with a trailing `/1` or `/2` mate suffix
/// removed.
fn mate_key(name: &[u8]) -> &[u8] {
    name.strip_suffix(b"/1")
        .or_else(|| name.strip_suffix(b"/2"))
        .unwrap_or(name)
}

/// Verify every record in a fragment shares one mate-pairing key, catching
/// mis-ordered or desynced inputs that would otherwise silently combine
/// segments from different records.
fn check_names_concord(names: &[&[u8]]) -> Result<()> {
    let Some((first, rest)) = names.split_first() else {
        return Ok(());
    };
    let key = mate_key(first);
    for other in rest {
        if mate_key(other) != key {
            bail!(
                "input segments have mismatched read names (`{}` vs `{}`); the inputs are not in the same record order",
                String::from_utf8_lossy(first),
                String::from_utf8_lossy(other)
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod empty_input_tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;

    fn write(dir: &std::path::Path, name: &str, bytes: &[u8]) -> PathBuf {
        let path = dir.join(name);
        std::fs::File::create(&path)
            .unwrap()
            .write_all(bytes)
            .unwrap();
        path
    }

    /// Open a path and drain it, returning the sniffed formats and the record
    /// count.
    fn drain(path: &Path) -> (Vec<SniffedFormat>, usize) {
        let mut reader = FragmentReader::open(&[path.to_path_buf()], true).unwrap();
        let formats = reader.formats().to_vec();
        let mut count = 0;
        while reader.next_fragment().unwrap().is_some() {
            count += 1;
        }
        (formats, count)
    }

    #[test]
    fn test_sniff_from_extension_maps_known_extensions() {
        assert_eq!(
            sniff_from_extension(Path::new("x.fastq")),
            SniffedFormat::Fastq
        );
        assert_eq!(
            sniff_from_extension(Path::new("x.fq.gz")),
            SniffedFormat::Fastq
        );
        assert_eq!(
            sniff_from_extension(Path::new("x.fasta")),
            SniffedFormat::Fasta
        );
        assert_eq!(
            sniff_from_extension(Path::new("x.fa.gz")),
            SniffedFormat::Fasta
        );
        assert_eq!(
            sniff_from_extension(Path::new("x.fna")),
            SniffedFormat::Fasta
        );
        assert_eq!(sniff_from_extension(Path::new("x.sam")), SniffedFormat::Sam);
        assert_eq!(sniff_from_extension(Path::new("x.bam")), SniffedFormat::Bam);
        assert_eq!(
            sniff_from_extension(Path::new("x.cram")),
            SniffedFormat::Cram
        );
        assert_eq!(
            sniff_from_extension(Path::new("x.txt")),
            SniffedFormat::Unknown
        );
        assert_eq!(sniff_from_extension(Path::new("-")), SniffedFormat::Unknown);
    }

    #[test]
    fn test_empty_fastx_files_read_as_zero_records() {
        let dir = tempdir().unwrap();
        // 0-byte FASTQ and FASTA: sniffed from the extension, drained to zero
        // records.
        let (fmt, n) = drain(&write(dir.path(), "empty.fastq", b""));
        assert_eq!(fmt, vec![SniffedFormat::Fastq]);
        assert_eq!(n, 0);
        let (fmt, n) = drain(&write(dir.path(), "empty.fa", b""));
        assert_eq!(fmt, vec![SniffedFormat::Fasta]);
        assert_eq!(n, 0);
    }

    #[test]
    fn test_empty_gzip_fastq_reads_as_zero_records() {
        let dir = tempdir().unwrap();
        // A valid empty gzip member (decodes to zero bytes), like a sequencer's
        // empty .fastq.gz.
        let encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        let gz = encoder.finish().unwrap();
        let (fmt, n) = drain(&write(dir.path(), "empty.fastq.gz", &gz));
        assert_eq!(fmt, vec![SniffedFormat::Fastq]);
        assert_eq!(n, 0);
    }

    #[test]
    fn test_extensionless_empty_input_is_empty_not_an_error() {
        let dir = tempdir().unwrap();
        // No content and no recognizable extension: resolves to Empty, read as
        // zero records.
        let (fmt, n) = drain(&write(dir.path(), "empty.dat", b""));
        assert_eq!(fmt, vec![SniffedFormat::Empty]);
        assert_eq!(n, 0);
    }

    #[test]
    fn test_header_only_sam_reads_as_zero_records() {
        let dir = tempdir().unwrap();
        let (fmt, n) = drain(&write(dir.path(), "h.sam", b"@HD\tVN:1.6\n"));
        assert_eq!(fmt, vec![SniffedFormat::Sam]);
        assert_eq!(n, 0);
    }

    #[test]
    fn test_header_only_bam_reads_as_zero_records() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("h.bam");
        let header = noodles::sam::Header::default();
        let mut writer = noodles::bam::io::Writer::new(std::fs::File::create(&path).unwrap());
        writer.write_header(&header).unwrap();
        writer.try_finish().unwrap();
        let (fmt, n) = drain(&path);
        assert_eq!(fmt, vec![SniffedFormat::Bam]);
        assert_eq!(n, 0);
    }

    #[test]
    fn test_header_only_cram_reads_as_zero_records() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("h.cram");
        let header = noodles::sam::Header::default();
        let mut writer = noodles::cram::io::Writer::new(std::fs::File::create(&path).unwrap());
        writer.write_header(&header).unwrap();
        writer.try_finish(&header).unwrap();
        let (fmt, n) = drain(&path);
        assert_eq!(fmt, vec![SniffedFormat::Cram]);
        assert_eq!(n, 0);
    }

    #[test]
    fn test_nonempty_unrecognized_bytes_still_error() {
        let dir = tempdir().unwrap();
        let path = write(dir.path(), "junk.dat", b"this is not a sequence file\n");
        assert!(FragmentReader::open(&[path], true).is_err());
    }

    #[test]
    fn test_empty_input_mixed_with_nonempty_errors_on_record_count() {
        let dir = tempdir().unwrap();
        // File 0 is an extensionless empty input (sniffs as Empty); file 1 is a
        // non-empty FASTQ. Open must succeed (an empty input no longer trips
        // the quality-axis guard); the real mismatch surfaces at read time with
        // the accurate record-count message, not a misleading "FASTA cannot be
        // blended" quality error.
        let empty = write(dir.path(), "empty.dat", b"");
        let full = write(dir.path(), "reads.fq", b"@r1\nACGT\n+\nIIII\n");
        let mut reader = FragmentReader::open(&[empty, full], true).unwrap();
        let err = loop {
            match reader.next_fragment() {
                Ok(Some(_)) => continue,
                Ok(None) => {
                    panic!("expected an unequal-record-count error, got clean end of input")
                }
                Err(e) => break e,
            }
        };
        let msg = format!("{err}");
        assert!(
            msg.contains("unequal record counts"),
            "expected the record-count error, got: {msg}"
        );
        assert!(
            !msg.contains("FASTA"),
            "must not surface the misleading quality message: {msg}"
        );
    }

    #[test]
    fn test_all_empty_mixed_extension_and_extensionless_reads_zero_records() {
        let dir = tempdir().unwrap();
        // An extensionless empty input (Empty) plus an empty .fq (Fastq): both
        // carry zero records, so the run must succeed and read as zero
        // fragments rather than erroring on the format mismatch.
        let empty = write(dir.path(), "empty.dat", b"");
        let empty_fq = write(dir.path(), "empty.fq", b"");
        let mut reader = FragmentReader::open(&[empty, empty_fq], true).unwrap();
        assert!(
            reader.next_fragment().unwrap().is_none(),
            "all-empty inputs read as zero records"
        );
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    /// Gzip-compress `bytes` into a fresh buffer.
    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    /// Write `bytes` to a fresh temp file and return the handle (kept alive by
    /// the caller).
    fn temp_with(bytes: &[u8]) -> tempfile::NamedTempFile {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(bytes).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_corrupt_skipped_bam_record_fails_the_read() {
        use noodles::sam::alignment::io::Write as _;
        let text = b"@SQ\tSN:chr1\tLN:100\n\
                     r1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tABCD\n\
                     r1\t256\tchr1\t20\t0\t4M\t*\t0\t0\t*\t*\tXX:Z:broken\n\
                     r2\t2048\tchr1\t30\t60\t4M\t*\t0\t0\tACGT\tABCD\tXX:Z:broken\n\
                     r2\t0\tchr1\t40\t60\t4M\t*\t0\t0\tACGT\tABCD\n";
        let mut reader = sam::io::Reader::new(&text[..]);
        let header = reader.read_header().unwrap();
        let mut writer = bam::io::Writer::from(Vec::new());
        writer.write_header(&header).unwrap();
        for result in reader.record_bufs(&header) {
            writer
                .write_alignment_record(&header, &result.unwrap())
                .unwrap();
        }
        let raw = writer.into_inner();
        let tag = b"XXZbroken\0";
        let corrupted = |nth: usize| {
            let mut raw = raw.clone();
            let at = raw
                .windows(tag.len())
                .enumerate()
                .filter(|(_, window)| *window == tag)
                .nth(nth)
                .unwrap()
                .0;
            raw[at + 2] = b'?';
            let mut bgzf = bgzf::io::Writer::new(Vec::new());
            bgzf.write_all(&raw).unwrap();
            temp_with(&bgzf.finish().unwrap())
        };
        for (nth, threads) in [(0, 0), (1, 0), (0, 2), (1, 2)] {
            let file = corrupted(nth);
            let path = file.path().to_path_buf();
            let mut reader =
                FragmentReader::open_with_threads(&[path], true, None, |_| vec![threads]).unwrap();
            let error = std::iter::from_fn(|| reader.next_fragment().transpose())
                .find_map(Result::err)
                .expect("a corrupt skipped record must fail the read");
            assert!(
                format!("{error:#}").contains("failed to read BAM record"),
                "{error:#}"
            );
        }
    }

    #[test]
    fn test_sniff_bytes_bam_magic() {
        assert_eq!(sniff_bytes(b"BAM\x01rest"), SniffedFormat::Bam);
    }

    #[test]
    fn test_sniff_bytes_cram_magic() {
        assert_eq!(sniff_bytes(b"CRAM\x03"), SniffedFormat::Cram);
    }

    #[test]
    fn test_sniff_bytes_fasta() {
        assert_eq!(sniff_bytes(b">read1\nACGT\n"), SniffedFormat::Fasta);
    }

    #[test]
    fn test_sniff_bytes_fastq() {
        assert_eq!(
            sniff_bytes(b"@read1\nACGT\n+\nIIII\n"),
            SniffedFormat::Fastq
        );
    }

    #[test]
    fn test_sniff_bytes_sam_header_has_tab() {
        assert_eq!(sniff_bytes(b"@HD\tVN:1.6\n@SQ\tSN:1\n"), SniffedFormat::Sam);
    }

    #[test]
    fn test_sniff_bytes_at_prefixed_no_newline_falls_back_to_sam() {
        // `@HD\t...` with no newline in the chunk uses the byte-position rule.
        assert_eq!(sniff_bytes(b"@HD\tVN"), SniffedFormat::Sam);
    }

    #[test]
    fn test_sniff_bytes_at_prefixed_no_newline_falls_back_to_fastq() {
        assert_eq!(sniff_bytes(b"@read_name_only"), SniffedFormat::Fastq);
    }

    #[test]
    fn test_sniff_bytes_unknown_and_empty() {
        assert_eq!(sniff_bytes(b"random"), SniffedFormat::Unknown);
        assert_eq!(sniff_bytes(b""), SniffedFormat::Unknown);
    }

    #[test]
    fn test_sniff_input_plain_fastq() {
        let file = temp_with(b"@r1\nACGT\n+\nIIII\n");
        let (format, gzipped, _reader) = sniff_input(file.path()).unwrap();
        assert_eq!(format, SniffedFormat::Fastq);
        assert!(!gzipped);
    }

    #[test]
    fn test_sniff_input_gzipped_fastq() {
        let file = temp_with(&gzip(b"@r1\nACGT\n+\nIIII\n"));
        let (format, gzipped, _reader) = sniff_input(file.path()).unwrap();
        assert_eq!(format, SniffedFormat::Fastq);
        assert!(gzipped);
    }

    #[test]
    fn test_sniff_input_leaves_reader_at_byte_zero() {
        // The peek must not consume bytes: the caller reads the whole file from
        // the start.
        let payload = b"@r1\nACGT\n+\nIIII\n";
        let file = temp_with(payload);
        let (_format, _gzipped, mut reader) = sniff_input(file.path()).unwrap();
        let mut got = Vec::new();
        reader.read_to_end(&mut got).unwrap();
        assert_eq!(got, payload);
    }

    #[test]
    fn test_sniff_input_gzipped_cram_is_rejected() {
        let file = temp_with(&gzip(b"CRAM\x03somecramdata"));
        // The Ok variant holds a boxed `dyn BufRead` (not `Debug`), so match
        // instead of `unwrap_err`.
        let err = match sniff_input(file.path()) {
            Err(err) => err,
            Ok(_) => panic!("gzipped CRAM must be rejected"),
        };
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("gzipped CRAM"), "{err}");
    }

    /// Read every fragment from a single-input reader over a temp file.
    fn read_all(file: &tempfile::NamedTempFile) -> Vec<Fragment> {
        let mut reader = FragmentReader::open(&[file.path().to_path_buf()], false).unwrap();
        let mut fragments = Vec::new();
        while let Some(fragment) = reader.next_fragment().unwrap() {
            fragments.push(fragment);
        }
        fragments
    }

    #[test]
    fn test_read_fastq_records() {
        let file = temp_with(b"@r1 extra desc\nACGT\n+\nIIII\n@r2\nTTGG\n+\nJJJJ\n");
        let fragments = read_all(&file);
        assert_eq!(fragments.len(), 2);
        assert_eq!(fragments[0].records.len(), 1);
        // The trailing description is dropped; the mate-pairing token is kept.
        assert_eq!(fragments[0].records[0].name, b"r1");
        assert_eq!(fragments[0].records[0].bases, b"ACGT");
        // FASTQ ASCII 'I' (73) decodes to raw Phred score 40.
        assert_eq!(
            fragments[0].records[0].quals.as_deref(),
            Some(&[40u8, 40, 40, 40][..])
        );
        assert!(fragments[0].records[0].tags.is_none());
        assert_eq!(fragments[1].records[0].bases, b"TTGG");
    }

    /// The bytes of a `Z` (string) tag in a record_buf `Data`, for assertions.
    fn z_tag(data: &Data, tag: &[u8; 2]) -> Option<Vec<u8>> {
        use noodles::sam::alignment::record_buf::data::field::Value;
        match data.get(tag) {
            Some(Value::String(s)) => {
                let bytes: &[u8] = s.as_ref();
                Some(bytes.to_vec())
            }
            _ => None,
        }
    }

    #[test]
    fn test_sniff_fastq_with_tab_separated_comment_is_not_sam() {
        // A FASTQ read-name line may carry tab-separated SAM tags in its comment
        // (samtools `fastq -T`). A SAM header line is `@<2 uppercase>\t...`; a
        // read name does not match that shape, even with tabs later in the
        // comment, so it must still sniff as FASTQ.
        assert_eq!(
            sniff_bytes(b"@read1\tRX:Z:ACGT\tBC:Z:GGGG\nACGT\n+\nIIII\n"),
            SniffedFormat::Fastq
        );
        // A real SAM header still sniffs as SAM.
        assert_eq!(sniff_bytes(b"@HD\tVN:1.6\n"), SniffedFormat::Sam);
    }

    #[test]
    fn parse_comment_tags_lifts_valid_sam_tags() {
        let data = parse_comment_tags(b"RX:Z:ACGT\tBC:Z:GGGG").expect("valid tags lift");
        assert_eq!(z_tag(&data, b"RX").as_deref(), Some(&b"ACGT"[..]));
        assert_eq!(z_tag(&data, b"BC").as_deref(), Some(&b"GGGG"[..]));
    }

    #[test]
    fn parse_comment_tags_lifts_integer_tag() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        let data = parse_comment_tags(b"xi:i:42").expect("integer tag lifts");
        assert!(matches!(data.get(b"xi"), Some(Value::Int32(42))));
    }

    #[test]
    fn parse_comment_tags_skips_non_tag_fields_lenient() {
        // A bare-UMI / free-text field is silently skipped; the valid tag in the
        // same comment still lifts (the samtools-import lenient behavior).
        let data =
            parse_comment_tags(b"GATTACA\tRX:Z:ACGT").expect("the valid tag lifts past the junk");
        assert_eq!(z_tag(&data, b"RX").as_deref(), Some(&b"ACGT"[..]));
        assert!(data.get(b"GA").is_none(), "the junk field is not lifted");
    }

    #[test]
    fn parse_comment_tags_returns_none_for_non_tag_comments() {
        // A bare UMI, an Illumina CASAVA string, a malformed type code, and an
        // empty comment all yield no tags.
        assert!(parse_comment_tags(b"ACGTACGT").is_none(), "bare UMI");
        assert!(parse_comment_tags(b"1:N:0:ACGT").is_none(), "CASAVA");
        assert!(parse_comment_tags(b"RX:?:ACGT").is_none(), "bad type code");
        assert!(parse_comment_tags(b"").is_none(), "empty comment");
    }

    #[test]
    fn parse_comment_tags_preserves_spaces_in_z_value() {
        // Only tabs separate fields, so a space inside a Z value is kept verbatim.
        let data = parse_comment_tags(b"CO:Z:hello world").expect("Z tag lifts");
        assert_eq!(z_tag(&data, b"CO").as_deref(), Some(&b"hello world"[..]));
    }

    #[test]
    fn test_read_fastq_lifts_comment_tags() {
        // SAM tags in the FASTQ comment become the record's carried tags.
        let file = temp_with(b"@read1 RX:Z:ACGT\nACGT\n+\nIIII\n");
        let fragments = read_all(&file);
        let data = fragments[0].records[0]
            .tags
            .as_ref()
            .expect("comment tags are lifted");
        assert_eq!(z_tag(data, b"RX").as_deref(), Some(&b"ACGT"[..]));
    }

    #[test]
    fn test_read_fasta_lifts_comment_tags() {
        let file = temp_with(b">read1 BC:Z:GGGG\nACGTACGT\n");
        let fragments = read_all(&file);
        let data = fragments[0].records[0]
            .tags
            .as_ref()
            .expect("comment tags are lifted");
        assert_eq!(z_tag(data, b"BC").as_deref(), Some(&b"GGGG"[..]));
    }

    #[test]
    fn test_read_fastq_gzipped() {
        let file = temp_with(&gzip(b"@r1\nACGT\n+\nIIII\n"));
        let fragments = read_all(&file);
        assert_eq!(fragments.len(), 1);
        assert_eq!(fragments[0].records[0].bases, b"ACGT");
    }

    #[test]
    fn test_read_fasta_has_no_qualities() {
        let file = temp_with(b">r1 desc\nACGTACGT\n>r2\nTTTT\n");
        let fragments = read_all(&file);
        assert_eq!(fragments.len(), 2);
        assert_eq!(fragments[0].records[0].name, b"r1");
        assert_eq!(fragments[0].records[0].bases, b"ACGTACGT");
        assert!(fragments[0].records[0].quals.is_none());
    }

    #[test]
    fn test_multi_file_zip_pairs_by_record_order() {
        let r1 = temp_with(b"@a\nAAAA\n+\nIIII\n@b\nCCCC\n+\nIIII\n");
        let r2 = temp_with(b"@a\nGGGG\n+\nJJJJ\n@b\nTTTT\n+\nJJJJ\n");
        let mut reader =
            FragmentReader::open(&[r1.path().to_path_buf(), r2.path().to_path_buf()], false)
                .unwrap();
        let first = reader.next_fragment().unwrap().unwrap();
        assert_eq!(first.records.len(), 2);
        assert_eq!(first.records[0].bases, b"AAAA"); // input file 0
        assert_eq!(first.records[1].bases, b"GGGG"); // input file 1
        let second = reader.next_fragment().unwrap().unwrap();
        assert_eq!(second.records[0].bases, b"CCCC");
        assert_eq!(second.records[1].bases, b"TTTT");
        assert!(reader.next_fragment().unwrap().is_none());
    }

    #[test]
    fn test_unequal_input_lengths_is_error() {
        let r1 = temp_with(b"@a\nAAAA\n+\nIIII\n@b\nCCCC\n+\nIIII\n"); // 2 records
        let r2 = temp_with(b"@a\nGGGG\n+\nJJJJ\n"); // 1 record
        let mut reader =
            FragmentReader::open(&[r1.path().to_path_buf(), r2.path().to_path_buf()], false)
                .unwrap();
        assert!(reader.next_fragment().unwrap().is_some()); // first pair reads
        let err = reader.next_fragment().unwrap_err(); // r2 exhausted early
        assert!(err.to_string().contains("unequal record counts"), "{err}");
    }

    #[test]
    fn test_read_sam_carries_tags() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        // An unmapped SAM record (flag 4) with an RX:Z tag; the tag is carried
        // through.
        let sam = b"@HD\tVN:1.6\nread1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tRX:Z:AACC\n";
        let file = temp_with(sam);
        let fragments = read_all(&file);
        assert_eq!(fragments.len(), 1);
        let record = &fragments[0].records[0];
        assert_eq!(record.name, b"read1");
        assert_eq!(record.bases, b"ACGT");
        assert_eq!(record.quals.as_ref().map(|q| q.len()), Some(4));
        let data = record
            .tags
            .as_ref()
            .expect("alignment tags carried through");
        match data.get(b"RX") {
            Some(Value::String(seq)) => {
                let value: &[u8] = seq.as_ref();
                assert_eq!(value, &b"AACC"[..]);
            }
            other => panic!("expected RX:Z tag, got {other:?}"),
        }
    }

    #[test]
    fn test_read_sam_empty_read_has_empty_qualities() {
        let sam = b"@HD\tVN:1.6\n\
                    trimmed\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n\
                    unscored\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n";
        let file = temp_with(sam);
        let fragments = read_all(&file);
        assert_eq!(fragments[0].records[0].quals.as_deref(), Some(&[][..]));
        assert_eq!(fragments[1].records[0].quals, None);
    }

    #[test]
    fn test_read_cram_round_trip_reverse_and_secondary() {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
        use noodles::sam::Header;

        // Reference-free CRAM (no @SQ), like unmux's own output and raw demux
        // inputs.
        let header = Header::default();
        let record = |name: &[u8], flags: Flags, seq: &[u8], qual: Vec<u8>| {
            RecordBuf::builder()
                .set_name(name.to_vec())
                .set_flags(flags)
                .set_sequence(Sequence::from(seq.to_vec()))
                .set_quality_scores(QualityScores::from(qual))
                .build()
        };

        // Write a CRAM in memory: forward, reverse-strand, and a secondary
        // record (all unmapped).
        let mut cram_bytes = Vec::new();
        {
            let mut writer = cram::io::Writer::new(&mut cram_bytes);
            writer.write_header(&header).unwrap();
            for rec in [
                record(b"fwd", Flags::UNMAPPED, b"ACGT", vec![10, 11, 12, 13]),
                record(
                    b"rev",
                    Flags::UNMAPPED | Flags::REVERSE_COMPLEMENTED,
                    b"CGTT",
                    vec![15, 16, 17, 18],
                ),
                record(
                    b"sec",
                    Flags::UNMAPPED | Flags::SECONDARY,
                    b"TTTT",
                    vec![20, 20, 20, 20],
                ),
            ] {
                writer.write_alignment_record(&header, &rec).unwrap();
            }
            writer.try_finish(&header).unwrap();
        }

        // Reading exercises the container/slice decode, reverse-strand
        // re-orientation, and secondary skipping (identical to the BAM path,
        // via the shared RecordBuf extraction).
        let file = temp_with(&cram_bytes);
        let fragments = read_all(&file);
        assert_eq!(fragments.len(), 2, "secondary record must be skipped");
        assert_eq!(fragments[0].records[0].name, b"fwd");
        assert_eq!(fragments[0].records[0].bases, b"ACGT");
        assert_eq!(
            fragments[0].records[0].quals.as_deref(),
            Some(&[10u8, 11, 12, 13][..])
        );
        assert_eq!(fragments[1].records[0].name, b"rev");
        assert_eq!(fragments[1].records[0].bases, b"AACG"); // reverse complement of CGTT
    }

    #[test]
    fn test_cram_read_after_eof_stays_at_end() {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};

        let header = sam::Header::default();
        let record = RecordBuf::builder()
            .set_name(b"r1".to_vec())
            .set_flags(Flags::UNMAPPED)
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30; 4]))
            .build();
        let mut cram_bytes = Vec::new();
        {
            let mut writer = cram::io::Writer::new(&mut cram_bytes);
            writer.write_header(&header).unwrap();
            writer.write_alignment_record(&header, &record).unwrap();
            writer.try_finish(&header).unwrap();
        }
        let cram = temp_with(&cram_bytes);
        let mut reader = FragmentReader::open(&[cram.path().to_path_buf()], true).unwrap();
        assert!(reader.next_fragment().unwrap().is_some());
        assert!(reader.next_fragment().unwrap().is_none());
        assert!(reader.next_fragment().unwrap().is_none());
    }

    #[test]
    fn test_cram_unmapped_record_reads_mapq_zero() {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};

        let header = sam::Header::default();
        let record = RecordBuf::builder()
            .set_name(b"r1".to_vec())
            .set_flags(Flags::UNMAPPED)
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30; 4]))
            .build();
        let mut cram_bytes = Vec::new();
        {
            let mut writer = cram::io::Writer::new(&mut cram_bytes);
            writer.write_header(&header).unwrap();
            writer.write_alignment_record(&header, &record).unwrap();
            writer.try_finish(&header).unwrap();
        }
        let cram = temp_with(&cram_bytes);
        let mut reader = FragmentReader::open(&[cram.path().to_path_buf()], true).unwrap();
        reader.keep_alignments();
        let fragment = reader.next_fragment().unwrap().unwrap();
        let Some(Alignment::Record(kept)) = fragment.records[0].alignment.as_ref() else {
            panic!("the source record is kept");
        };
        assert_eq!(kept.mapping_quality().map(u8::from), Some(0));
    }

    #[test]
    fn test_decode_shares_split_workers_across_the_inputs_of_one_format() {
        use SniffedFormat::{Bam, Cram, Fastq};
        assert_eq!(decode_shares(&[Bam], Bam, 2), vec![2]);
        assert_eq!(decode_shares(&[Bam, Cram], Cram, 3), vec![0, 3]);
        assert_eq!(decode_shares(&[Bam, Cram], Bam, 0), vec![0, 0]);
        assert_eq!(decode_shares(&[Fastq, Bam, Bam], Bam, 3), vec![0, 2, 1]);
        assert_eq!(decode_shares(&[Fastq], Bam, 4), vec![0]);
    }

    #[test]
    fn test_pooled_cram_and_bam_decode_match_inline() {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};

        let header = sam::Header::default();
        let records: Vec<RecordBuf> = (0..25_000)
            .map(|i| {
                let flags = if i % 5 == 0 {
                    Flags::UNMAPPED | Flags::SECONDARY
                } else {
                    Flags::UNMAPPED
                };
                RecordBuf::builder()
                    .set_name(format!("r{i}").into_bytes())
                    .set_flags(flags)
                    .set_sequence(Sequence::from(b"ACGTN"[..1 + i % 5].to_vec()))
                    .set_quality_scores(QualityScores::from(vec![30; 1 + i % 5]))
                    .build()
            })
            .collect();
        let mut cram_bytes = Vec::new();
        {
            let mut writer = cram::io::Writer::new(&mut cram_bytes);
            writer.write_header(&header).unwrap();
            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }
            writer.try_finish(&header).unwrap();
        }
        let mut bam_bytes = Vec::new();
        {
            let mut writer = bam::io::Writer::new(&mut bam_bytes);
            writer.write_header(&header).unwrap();
            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }
            writer.try_finish().unwrap();
        }
        for bytes in [cram_bytes, bam_bytes] {
            let file = temp_with(&bytes);
            let read = |threads: usize, keep: bool| {
                let path = file.path().to_path_buf();
                let mut reader =
                    FragmentReader::open_with_threads(&[path], true, None, |_| vec![threads])
                        .unwrap();
                if keep {
                    reader.keep_alignments();
                }
                let mut fragments = Vec::new();
                while let Some(fragment) = reader.next_fragment().unwrap() {
                    fragments.push(fragment);
                }
                fragments
            };
            for keep in [false, true] {
                let inline = read(0, keep);
                assert_eq!(inline.len(), if keep { 25_000 } else { 20_000 });
                assert_eq!(read(3, keep), inline);
            }
        }
    }

    #[test]
    fn test_keep_alignments_yields_every_record_whole() {
        let sam = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
                    p1\t16\tchr1\t10\t60\t4M\t*\t0\t0\tAACC\tABCD\tNM:i:0\n\
                    s1\t256\tchr1\t20\t0\t4M\t*\t0\t0\t*\t*\n\
                    x1\t2048\tchr1\t30\t9\t2H2M\t*\t0\t0\tGG\tII\n";
        let file = temp_with(sam);
        let mut reader = FragmentReader::open(&[file.path().to_path_buf()], false).unwrap();
        reader.keep_alignments();
        assert!(reader.alignment_header().is_some());
        let mut kept = Vec::new();
        while let Some(fragment) = reader.next_fragment().unwrap() {
            kept.push(fragment.records.into_iter().next().unwrap());
        }
        let names: Vec<&[u8]> = kept.iter().map(|r| r.name.as_slice()).collect();
        assert_eq!(names, vec![&b"p1"[..], b"s1", b"x1"]);
        let Some(Alignment::Record(first)) = kept[0].alignment.as_ref() else {
            panic!("the source record is kept");
        };
        assert_eq!(first.alignment_start().map(usize::from), Some(10));
        assert_eq!(first.sequence().as_ref(), b"AACC");
        assert!(first.data().get(&Tag::new(b'N', b'M')).is_some());
        assert_eq!(
            kept[0].bases, b"GGTT",
            "matching still sees the sequenced orientation"
        );
    }

    #[test]
    fn test_default_reader_keeps_no_alignment() {
        let file = temp_with(b"@HD\tVN:1.6\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");
        let fragments = read_all(&file);
        assert!(fragments[0].records[0].alignment.is_none());
    }

    #[test]
    fn test_reference_compressed_cram_without_reference_names_the_sequence() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/aligned.htslib-3.0.cram");
        let mut reader = FragmentReader::open(&[path], true).unwrap();
        let error = loop {
            match reader.next_fragment() {
                Ok(Some(_)) => continue,
                Ok(None) => panic!("a reference-compressed CRAM decoded without a reference"),
                Err(error) => break error,
            }
        };
        assert!(
            format!("{error:#}").contains("needs reference sequence `chr1` to decode"),
            "{error:#}"
        );
    }

    #[test]
    fn test_cram_written_without_a_reference_reads_without_one() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/contigs.no-ref.cram");
        let mut reader = FragmentReader::open(&[path], true).unwrap();
        let mut count = 0;
        while reader.next_fragment().unwrap().is_some() {
            count += 1;
        }
        assert_eq!(count, 26);
    }

    #[test]
    fn test_read_sam_reverse_strand_reorients() {
        // A reverse-strand record (flag 16) stores SEQ reverse-complemented and
        // QUAL reversed; the reader restores the original read orientation.
        // Stored CGTT/0123 -> AACG with reversed quals.
        let sam = b"@HD\tVN:1.6\nrev\t16\t*\t0\t0\t*\t*\t0\t0\tCGTT\t0123\n";
        let file = temp_with(sam);
        let fragments = read_all(&file);
        assert_eq!(fragments.len(), 1);
        assert_eq!(fragments[0].records[0].bases, b"AACG");
        // QUAL "0123" decodes to raw scores [15,16,17,18]; reversed for the
        // 0x10 record.
        assert_eq!(
            fragments[0].records[0].quals.as_deref(),
            Some(&[18u8, 17, 16, 15][..])
        );
    }

    #[test]
    fn test_read_sam_skips_secondary() {
        // The secondary record (flag 256) is not a distinct read and is
        // skipped.
        let sam = b"@HD\tVN:1.6\n\
                    p1\t0\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n\
                    s1\t256\t*\t0\t0\t*\t*\t0\t0\tCCCC\tIIII\n\
                    p2\t0\t*\t0\t0\t*\t*\t0\t0\tGGGG\tIIII\n";
        let file = temp_with(sam);
        let fragments = read_all(&file);
        assert_eq!(fragments.len(), 2);
        assert_eq!(fragments[0].records[0].bases, b"AAAA");
        assert_eq!(fragments[1].records[0].bases, b"GGGG");
    }

    #[test]
    fn test_interleaved_fastq_slash_suffix() {
        // A single FASTQ with /1,/2 mates is auto-detected as interleaved (two
        // segments per read).
        let file = temp_with(
            b"@x/1\nAAAA\n+\nIIII\n@x/2\nCCCC\n+\nIIII\n@y/1\nGGGG\n+\nIIII\n@y/2\nTTTT\n+\nIIII\n",
        );
        let fragments = read_all(&file);
        assert_eq!(fragments.len(), 2);
        assert_eq!(fragments[0].records.len(), 2);
        assert_eq!(fragments[0].records[0].bases, b"AAAA");
        assert_eq!(fragments[0].records[1].bases, b"CCCC");
        assert_eq!(fragments[1].records[0].bases, b"GGGG");
        assert_eq!(fragments[1].records[1].bases, b"TTTT");
    }

    #[test]
    fn test_interleaved_fastq_casava_matching_names() {
        // Casava-style mates share a name before the space; interleaving is
        // detected by that key.
        let file = temp_with(b"@read1 1:N:0:AC\nAAAA\n+\nIIII\n@read1 2:N:0:AC\nCCCC\n+\nIIII\n");
        let fragments = read_all(&file);
        assert_eq!(fragments.len(), 1);
        assert_eq!(fragments[0].records.len(), 2);
        assert_eq!(fragments[0].records[0].bases, b"AAAA");
        assert_eq!(fragments[0].records[1].bases, b"CCCC");
    }

    #[test]
    fn test_per_record_disables_interleaving() {
        // The same /1,/2 stream read with --per-record stays single-end (one
        // record per fragment).
        let file = temp_with(b"@x/1\nAAAA\n+\nIIII\n@x/2\nCCCC\n+\nIIII\n");
        let mut reader = FragmentReader::open(&[file.path().to_path_buf()], true).unwrap();
        let mut count = 0;
        while let Some(fragment) = reader.next_fragment().unwrap() {
            assert_eq!(fragment.records.len(), 1);
            count += 1;
        }
        assert_eq!(count, 2);
    }

    const PAIRED_SAM: &[u8] = b"@HD\tVN:1.6\n\
        a\t77\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n\
        a\t141\t*\t0\t0\t*\t*\t0\t0\tCCCC\tIIII\n\
        b\t77\t*\t0\t0\t*\t*\t0\t0\tGGGG\tIIII\n\
        b\t141\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII\n";

    /// The pairs of [`PAIRED_SAM`] among supplementary and secondary records.
    const PAIRED_AMONG_OTHERS_SAM: &[u8] = b"@HD\tVN:1.6\n\
        a\t2125\t*\t0\t0\t*\t*\t0\t0\tNNNN\tIIII\n\
        a\t77\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n\
        a\t141\t*\t0\t0\t*\t*\t0\t0\tCCCC\tIIII\n\
        b\t77\t*\t0\t0\t*\t*\t0\t0\tGGGG\tIIII\n\
        b\t397\t*\t0\t0\t*\t*\t0\t0\tNNNN\tIIII\n\
        b\t141\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII\n";

    /// The records of SAM text as BAM bytes and as CRAM bytes.
    fn bam_and_cram(sam: &[u8]) -> (Vec<u8>, Vec<u8>) {
        use noodles::sam::alignment::io::Write as _;
        let mut reader = sam::io::Reader::new(sam);
        let header = reader.read_header().unwrap();
        let records: Vec<RecordBuf> = reader
            .record_bufs(&header)
            .collect::<std::io::Result<_>>()
            .unwrap();
        let mut bam_bytes = Vec::new();
        {
            let mut writer = bam::io::Writer::new(&mut bam_bytes);
            writer.write_header(&header).unwrap();
            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }
            writer.try_finish().unwrap();
        }
        let mut cram_bytes = Vec::new();
        {
            let mut writer = cram::io::Writer::new(&mut cram_bytes);
            writer.write_header(&header).unwrap();
            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }
            writer.try_finish(&header).unwrap();
        }
        (bam_bytes, cram_bytes)
    }

    /// The bases of each fragment's records, for one input read on `threads` decode workers.
    fn fragment_bases(bytes: &[u8], per_record: bool, threads: usize) -> Result<Vec<Vec<String>>> {
        let file = temp_with(bytes);
        let path = file.path().to_path_buf();
        let mut reader =
            FragmentReader::open_with_threads(&[path], per_record, None, |_| vec![threads])?;
        let bases = |record: &InputRecord| String::from_utf8_lossy(&record.bases).into_owned();
        let mut fragments = Vec::new();
        while let Some(fragment) = reader.next_fragment()? {
            assert_eq!(fragment.records.len(), reader.fragment_width());
            fragments.push(fragment.records.iter().map(bases).collect());
        }
        Ok(fragments)
    }

    #[test]
    fn test_paired_sam_bam_and_cram_yield_two_record_fragments() {
        for sam in [PAIRED_SAM, PAIRED_AMONG_OTHERS_SAM] {
            let (bam, cram) = bam_and_cram(sam);
            for bytes in [sam.to_vec(), bam, cram] {
                for threads in [0, 3] {
                    assert_eq!(
                        fragment_bases(&bytes, false, threads).unwrap(),
                        vec![vec!["AAAA", "CCCC"], vec!["GGGG", "TTTT"]]
                    );
                }
            }
        }
    }

    #[test]
    fn test_mates_are_ordered_by_their_flags() {
        let sam: &[u8] = b"@HD\tVN:1.6\n\
            a\t141\t*\t0\t0\t*\t*\t0\t0\tCCCC\tIIII\n\
            a\t77\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n\
            b\t77\t*\t0\t0\t*\t*\t0\t0\tGGGG\tIIII\n\
            b\t141\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII\n\
            c\t141\t*\t0\t0\t*\t*\t0\t0\tCGCG\tIIII\n\
            c\t77\t*\t0\t0\t*\t*\t0\t0\tATAT\tIIII\n";
        let (bam, cram) = bam_and_cram(sam);
        for bytes in [sam.to_vec(), bam, cram] {
            for threads in [0, 3] {
                assert_eq!(
                    fragment_bases(&bytes, false, threads).unwrap(),
                    vec![
                        vec!["AAAA", "CCCC"],
                        vec!["GGGG", "TTTT"],
                        vec!["ATAT", "CGCG"]
                    ]
                );
            }
        }
    }

    #[test]
    fn test_pairs_stay_together_across_cram_containers_and_bgzf_blocks() {
        let pairs = 15_000;
        let mut sam = String::from("@HD\tVN:1.6\n");
        for i in 0..pairs {
            sam.push_str(&format!("r{i}\t77\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n"));
            if i % 7 == 0 {
                sam.push_str(&format!("r{i}\t333\t*\t0\t0\t*\t*\t0\t0\tNNNN\tIIII\n"));
            }
            sam.push_str(&format!("r{i}\t141\t*\t0\t0\t*\t*\t0\t0\tCCCC\tIIII\n"));
        }
        let (bam, cram) = bam_and_cram(sam.as_bytes());
        for bytes in [bam, cram] {
            for threads in [0, 3] {
                let fragments = fragment_bases(&bytes, false, threads).unwrap();
                assert_eq!(fragments.len(), pairs);
                assert!(fragments.iter().all(|bases| bases == &["AAAA", "CCCC"]));
            }
        }
    }

    #[test]
    fn test_paired_bam_mates_stay_raw_until_decoded() {
        let (bam, _) = bam_and_cram(PAIRED_SAM);
        let file = temp_with(&bam);
        let mut reader = FragmentReader::open(&[file.path().to_path_buf()], false).unwrap();
        let first = reader.next_pending_fragment().unwrap().unwrap();
        let second = reader.next_pending_fragment().unwrap().unwrap();
        assert!(reader.next_pending_fragment().unwrap().is_none());
        assert!(second
            .records
            .iter()
            .all(|record| matches!(record, PendingRecord::Bam { keep: false, .. })));
        let mut bases: Vec<Vec<Vec<u8>>> = Vec::new();
        for pending in [first, second] {
            let fragment = pending.decode(&mut RecordBuf::default()).unwrap();
            assert!(fragment.records.iter().all(|r| r.alignment.is_none()));
            bases.push(fragment.records.into_iter().map(|r| r.bases).collect());
        }
        let expected: Vec<Vec<Vec<u8>>> = vec![
            vec![b"AAAA".to_vec(), b"CCCC".to_vec()],
            vec![b"GGGG".to_vec(), b"TTTT".to_vec()],
        ];
        assert_eq!(bases, expected);
    }

    #[test]
    fn test_unpairable_alignment_records_are_errors() {
        let mates_error = |flags: &str| {
            format!(
                "adjacent records `a` have flags {flags}, not those of one first mate (0x40) and \
                 one last mate (0x80); pass --per-record"
            )
        };
        let cases: [(&[(&str, u16)], String); 9] = [
            (
                &[("a", 77), ("b", 77), ("a", 141), ("b", 141)],
                "paired record `a` is followed by `b`, not its mate; sort the input".to_string(),
            ),
            (
                &[("a", 77), ("a", 141), ("b", 77)],
                "paired record `b` ends the input without its mate; sort the input".to_string(),
            ),
            (
                &[("a", 77), ("a", 141), ("b", 4)],
                "input mixes paired and unpaired records (`b` is unpaired".to_string(),
            ),
            (
                &[("a", 4), ("b", 77), ("b", 141)],
                "input mixes paired and unpaired records (`b` is paired".to_string(),
            ),
            (&[("a", 77), ("a", 77)], mates_error("77 and 77")),
            (&[("a", 141), ("a", 141)], mates_error("141 and 141")),
            (&[("a", 77), ("a", 4)], mates_error("77 and 4")),
            (&[("a", 13), ("a", 141)], mates_error("13 and 141")),
            (&[("a", 205), ("a", 141)], mates_error("205 and 141")),
        ];
        for (records, expected) in cases {
            let mut sam = String::from("@HD\tVN:1.6\n");
            for (name, flags) in records {
                sam.push_str(&format!("{name}\t{flags}\t*\t0\t0\t*\t*\t0\t0\tA\tI\n"));
            }
            let error = match fragment_bases(sam.as_bytes(), false, 0) {
                Ok(fragments) => panic!("expected `{expected}`, read {fragments:?}"),
                Err(error) => format!("{error:#}"),
            };
            assert!(error.contains(&expected), "{error}");
            assert!(error.contains("--per-record"), "{error}");
        }
    }

    #[test]
    fn test_per_record_reads_paired_alignments_one_per_fragment() {
        let (bam, cram) = bam_and_cram(PAIRED_SAM);
        for bytes in [PAIRED_SAM.to_vec(), bam, cram] {
            assert_eq!(
                fragment_bases(&bytes, true, 0).unwrap(),
                vec![vec!["AAAA"], vec!["CCCC"], vec!["GGGG"], vec!["TTTT"]]
            );
        }
    }

    #[test]
    fn test_keep_alignments_reads_paired_records_one_per_fragment_in_input_order() {
        let (bam, cram) = bam_and_cram(PAIRED_AMONG_OTHERS_SAM);
        for bytes in [PAIRED_AMONG_OTHERS_SAM.to_vec(), bam, cram] {
            let file = temp_with(&bytes);
            let mut reader = FragmentReader::open(&[file.path().to_path_buf()], false).unwrap();
            reader.keep_alignments();
            assert_eq!(reader.fragment_width(), 1);
            let mut flags = Vec::new();
            while let Some(fragment) = reader.next_fragment().unwrap() {
                assert_eq!(fragment.records.len(), 1);
                let alignment = fragment.records[0].alignment.as_ref().unwrap();
                flags.push(u16::from(alignment.flags()));
            }
            assert_eq!(flags, vec![2125, 77, 141, 77, 397, 141]);
        }
    }

    #[test]
    fn test_unreadable_first_alignment_record_fails_the_first_read_not_the_open() {
        let file = temp_with(b"@HD\tVN:1.6\na\t77\t*\tx\t0\t*\t*\t0\t0\tAAAA\tIIII\n");
        let mut reader = FragmentReader::open(&[file.path().to_path_buf()], false).unwrap();
        assert_eq!(reader.fragment_width(), 1);
        let error = reader.next_fragment().unwrap_err();
        assert!(
            format!("{error:#}").contains("failed to read SAM record"),
            "{error:#}"
        );
    }

    #[test]
    fn test_multi_file_name_mismatch_is_error() {
        // Two inputs whose record names disagree would silently mispair; the
        // reader fails fast.
        let r1 = temp_with(b"@a\nAAAA\n+\nIIII\n");
        let r2 = temp_with(b"@b\nGGGG\n+\nIIII\n");
        let mut reader =
            FragmentReader::open(&[r1.path().to_path_buf(), r2.path().to_path_buf()], false)
                .unwrap();
        let err = reader.next_fragment().unwrap_err();
        assert!(err.to_string().contains("mismatched read names"), "{err}");
    }

    #[test]
    fn test_read_bam_round_trip_reverse_and_secondary() {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
        use noodles::sam::Header;

        let header = Header::default();
        let record = |name: &[u8], flags: Flags, seq: &[u8], qual: Vec<u8>| {
            RecordBuf::builder()
                .set_name(name.to_vec())
                .set_flags(flags)
                .set_sequence(Sequence::from(seq.to_vec()))
                .set_quality_scores(QualityScores::from(qual))
                .build()
        };

        // Write a BGZF BAM in memory: forward, reverse-strand, and a secondary
        // record.
        let mut bam_bytes = Vec::new();
        {
            let mut writer = bam::io::Writer::new(&mut bam_bytes);
            writer.write_header(&header).unwrap();
            writer
                .write_alignment_record(
                    &header,
                    &record(b"fwd", Flags::UNMAPPED, b"ACGT", vec![10, 11, 12, 13]),
                )
                .unwrap();
            writer
                .write_alignment_record(
                    &header,
                    &record(
                        b"rev",
                        Flags::REVERSE_COMPLEMENTED,
                        b"CGTT",
                        vec![15, 16, 17, 18],
                    ),
                )
                .unwrap();
            writer
                .write_alignment_record(
                    &header,
                    &record(b"sec", Flags::SECONDARY, b"TTTT", vec![20, 20, 20, 20]),
                )
                .unwrap();
        }

        // Reading exercises the BGZF wiring, reverse-strand re-orientation, and
        // secondary skipping.
        let file = temp_with(&bam_bytes);
        let fragments = read_all(&file);
        assert_eq!(fragments.len(), 2, "secondary record must be skipped");
        assert_eq!(fragments[0].records[0].name, b"fwd");
        assert_eq!(fragments[0].records[0].bases, b"ACGT");
        assert_eq!(
            fragments[0].records[0].quals.as_deref(),
            Some(&[10u8, 11, 12, 13][..])
        );
        assert_eq!(fragments[1].records[0].name, b"rev");
        assert_eq!(fragments[1].records[0].bases, b"AACG"); // reverse complement of CGTT
        assert_eq!(
            fragments[1].records[0].quals.as_deref(),
            Some(&[18u8, 17, 16, 15][..])
        );
    }

    #[test]
    fn test_format_has_qualities() {
        assert!(!SniffedFormat::Fasta.has_qualities());
        assert!(SniffedFormat::Fastq.has_qualities());
        assert!(SniffedFormat::Sam.has_qualities());
        assert!(SniffedFormat::Bam.has_qualities());
        assert!(SniffedFormat::Cram.has_qualities());
    }

    #[test]
    fn test_ensure_compatible_quality() {
        use SniffedFormat::{Bam, Cram, Fasta, Fastq, Sam};
        // All quality-bearing formats blend freely, regardless of FASTX vs
        // alignment.
        assert!(ensure_compatible_quality(&[Fastq, Bam]).is_ok());
        assert!(ensure_compatible_quality(&[Sam, Bam, Cram]).is_ok());
        assert!(ensure_compatible_quality(&[Fastq, Sam, Cram]).is_ok());
        assert!(ensure_compatible_quality(&[Fasta, Fasta]).is_ok()); // both quality-less
        assert!(ensure_compatible_quality(&[Fastq]).is_ok());
        assert!(ensure_compatible_quality(&[]).is_ok());
        // FASTA (no qualities) cannot blend with a quality-bearing format.
        assert!(ensure_compatible_quality(&[Fasta, Fastq]).is_err());
        assert!(ensure_compatible_quality(&[Bam, Fasta]).is_err());
    }

    #[test]
    fn test_fasta_blended_with_quals_format_is_error() {
        // FASTA has no qualities, so blending it with FASTQ is a fail-fast
        // error.
        let fa = temp_with(b">a\nACGT\n");
        let fq = temp_with(b"@a\nGGGG\n+\nIIII\n");
        let err = FragmentReader::open(&[fa.path().to_path_buf(), fq.path().to_path_buf()], false)
            .err()
            .expect("FASTA cannot blend with a quality-bearing format");
        assert!(format!("{err:#}").contains("per-base qualities"), "{err:#}");
    }

    #[test]
    fn test_quals_bearing_blend_ok() {
        // FASTQ + SAM both carry qualities, so the blend is allowed (the SAM
        // segment carries tags).
        let fq = temp_with(b"@a\nACGT\n+\nIIII\n");
        let sam = temp_with(b"@HD\tVN:1.6\na\t4\t*\t0\t0\t*\t*\t0\t0\tGGGG\tIIII\n");
        let mut reader =
            FragmentReader::open(&[fq.path().to_path_buf(), sam.path().to_path_buf()], false)
                .unwrap();
        let fragment = reader.next_fragment().unwrap().unwrap();
        assert_eq!(fragment.records.len(), 2);
        assert_eq!(fragment.records[0].bases, b"ACGT"); // FASTQ
        assert!(fragment.records[0].tags.is_none());
        assert_eq!(fragment.records[1].bases, b"GGGG"); // SAM
        assert!(fragment.records[1].tags.is_some());
    }

    #[test]
    fn test_alignment_tags_are_stripped() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        // Alignment-derived tags (NM/AS/MD/XA/ms/MC) are dropped on read; a UMI
        // tag (RX) is kept.
        let sam = b"@HD\tVN:1.6\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\t\
                    NM:i:3\tAS:i:40\tMD:Z:4\tXA:Z:chr2,+9,4M,1\tms:i:50\tMC:Z:4M\tRX:Z:ACGTACGT\n";
        let file = temp_with(sam);
        let fragments = read_all(&file);
        let data = fragments[0].records[0].tags.as_ref().expect("carried tags");
        for stripped in [b"NM", b"AS", b"MD", b"XA", b"ms", b"MC"] {
            assert!(
                data.get(stripped).is_none(),
                "{} should be stripped",
                String::from_utf8_lossy(stripped)
            );
        }
        match data.get(b"RX") {
            Some(Value::String(value)) => assert_eq!(value.to_string(), "ACGTACGT"),
            other => panic!("RX (a UMI tag) should be kept: {other:?}"),
        }
    }

    #[test]
    fn input_read_groups_reads_sm_lb_and_sq() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("in.sam");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "@HD\tVN:1.6").unwrap();
        writeln!(f, "@SQ\tSN:chr1\tLN:1000").unwrap();
        writeln!(f, "@RG\tID:rg1\tSM:sampleA\tLB:lib1").unwrap();
        writeln!(f, "@RG\tID:rg2\tSM:sampleB").unwrap();
        drop(f);

        let reader = FragmentReader::open(&[path], false).unwrap();
        let info = reader.input_read_groups().unwrap();
        assert!(info.has_sq);
        assert_eq!(info.read_groups.len(), 2);
        assert_eq!(info.read_groups[0].id, b"rg1");
        assert_eq!(info.read_groups[0].sample.as_deref(), Some("sampleA"));
        assert_eq!(info.read_groups[0].library.as_deref(), Some("lib1"));
        assert_eq!(info.read_groups[1].sample.as_deref(), Some("sampleB"));
        assert_eq!(info.read_groups[1].library, None);
    }
}
