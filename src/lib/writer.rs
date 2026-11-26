//! Output writing: pick the output format (from the `--out` extension, or by
//! mirroring the inputs to stdout) and emit assembled output records to a
//! single sink.
//!
//! Qualities arrive as raw Phred scores (the reader normalized them; see
//! [`crate::input`]); this module re-encodes them per output format: ASCII
//! Phred+33 for FASTQ, dropped for FASTA, raw for SAM/BAM. SAM/BAM records
//! carry through the input's tags, the assembled `--tag` fields, and an `RG`
//! read group (matching an `@RG` line in the header passed to
//! [`OutputWriter::create_with_header`]), and a multi-record fragment gets the
//! paired-segment flags. Output is finalized explicitly via
//! [`OutputWriter::finish`] so a failed gzip trailer / BGZF EOF / flush is
//! reported, not silently swallowed. A run's gzipped FASTX and BAM outputs are
//! instead compressed block by block on the worker pool; see `crate::sink`.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

use anyhow::{anyhow, bail, Context, Result};
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::data::field::value::Array;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Data, QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, bgzf, cram, fasta, sam};

use crate::input::{Alignment, SniffedFormat, PHRED_OFFSET};

/// One output record to emit: its name, bases, raw-Phred qualities,
/// carried/derived tags, and the fan-out read group.
pub struct OutputRead<'a> {
    /// The read name (empty becomes a missing name, SAM/BAM `*`).
    pub name: &'a [u8],
    /// The bases.
    pub bases: &'a [u8],
    /// The raw-Phred qualities, when present.
    pub quals: Option<&'a [u8]>,
    /// SAM tags to emit: carried from an alignment input and/or assembled from
    /// `--tag`. On SAM/BAM they become record data fields; on FASTX they render
    /// into the read-name comment.
    pub tags: Option<&'a Data>,
    /// The record's fan-out read group id (`sample[.sub_sample]`), set for an
    /// assigned record on SAM/BAM output (an `RG:Z` data field, matching an
    /// `@RG` header line). FASTX has no read groups, so it is ignored there.
    pub read_group: Option<&'a str>,
    /// Whether the read fails quality control, written as flag `0x200` on SAM/BAM.
    pub qc_fail: bool,
    /// A source alignment to write on SAM/BAM in place of a record built from
    /// the fields above; BAM output keeps a raw BAM record's input bytes. FASTX
    /// ignores it.
    pub alignment: Option<&'a Alignment>,
}

/// A resolved output format. FASTX variants carry whether the sink is
/// gzip-compressed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// FASTQ, optionally gzip-compressed.
    Fastq {
        /// Whether to gzip the output.
        gzip: bool,
    },
    /// FASTA, optionally gzip-compressed.
    Fasta {
        /// Whether to gzip the output.
        gzip: bool,
    },
    /// Text SAM.
    Sam,
    /// BGZF-compressed BAM.
    Bam,
    /// Reference-free unmapped CRAM (its own per-block codecs; `--compression`
    /// does not apply).
    Cram,
}

impl OutputFormat {
    /// Whether a run writes this format as BGZF blocks compressed on the
    /// worker pool: gzipped FASTX (BGZF is valid gzip) and BAM. Text SAM, CRAM
    /// (its own codecs), and uncompressed FASTX are written as they merge.
    pub fn is_pooled(self) -> bool {
        matches!(
            self,
            OutputFormat::Fastq { gzip: true }
                | OutputFormat::Fasta { gzip: true }
                | OutputFormat::Bam
        )
    }
}

/// Resolve the output format. An explicit `--out` extension wins outright
/// (regardless of the input format); a `--out` whose path has no recognized
/// format extension, and stdout (`out = None`), mirror the input format
/// instead, defaulting a multi-format blend to uncompressed BAM. So a user who
/// asks for a format by suffix gets it, while an extensionless or stdout target
/// follows the inputs.
pub fn output_format(out: Option<&Path>, input_formats: &[SniffedFormat]) -> Result<OutputFormat> {
    match out {
        Some(path) => match format_from_extension(path)? {
            Some(format) => Ok(format),
            None => Ok(stdout_format(input_formats)),
        },
        None => Ok(stdout_format(input_formats)),
    }
}

/// Infer the output format from a path's extension (a trailing `.gz` means
/// gzip). Returns `Ok(None)` when no known format extension matches, so the
/// caller can fall back to mirroring the input format (an extensionless or
/// unknown-suffix `--out` then behaves like stdout). A
/// recognized-but-unsupported combination (`.sam.gz`, `.cram.gz`) is a hard
/// error, never a silent fallback.
fn format_from_extension(path: &Path) -> Result<Option<OutputFormat>> {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    let (stem, gzip) = match name.strip_suffix(".gz") {
        Some(stem) => (stem, true),
        None => (name.as_str(), false),
    };
    if stem.ends_with(".fq") || stem.ends_with(".fastq") {
        Ok(Some(OutputFormat::Fastq { gzip }))
    } else if stem.ends_with(".fa") || stem.ends_with(".fasta") || stem.ends_with(".fna") {
        Ok(Some(OutputFormat::Fasta { gzip }))
    } else if stem.ends_with(".sam") {
        if gzip {
            bail!("gzipped SAM output (.sam.gz) is not supported; use .bam for compressed alignment output");
        }
        Ok(Some(OutputFormat::Sam))
    } else if stem.ends_with(".bam") {
        Ok(Some(OutputFormat::Bam))
    } else if stem.ends_with(".cram") {
        if gzip {
            bail!("gzipped CRAM output (.cram.gz) is not supported; CRAM has its own block codecs");
        }
        Ok(Some(OutputFormat::Cram))
    } else {
        Ok(None)
    }
}

/// The stdout output format: mirror the single input format, or default a
/// multi-format blend to uncompressed BAM. stdout is never gzip-wrapped by
/// default.
fn stdout_format(input_formats: &[SniffedFormat]) -> OutputFormat {
    let all_same = input_formats.windows(2).all(|pair| pair[0] == pair[1]);
    if all_same {
        match input_formats.first() {
            Some(SniffedFormat::Fasta) => OutputFormat::Fasta { gzip: false },
            Some(SniffedFormat::Fastq) => OutputFormat::Fastq { gzip: false },
            Some(SniffedFormat::Sam) => OutputFormat::Sam,
            // A bare BAM/CRAM input mirrors to BAM; an empty input list also
            // defaults to BAM.
            _ => OutputFormat::Bam,
        }
    } else {
        OutputFormat::Bam
    }
}

/// The concrete byte sink under a format writer: a plain (buffered) stream, or
/// a gzip encoder over one. Keeping it concrete (not a boxed `dyn Write`) lets
/// [`SinkWriter::finish`] write the gzip trailer and propagate flush errors
/// instead of swallowing them on drop.
enum SinkWriter {
    /// A plain buffered stream (file or stdout); BAM wraps this in its own BGZF
    /// encoder.
    Plain(BufWriter<Box<dyn Write + Send>>),
    /// A gzip encoder over a buffered stream (FASTX.gz output); boxed as it is
    /// much larger than the plain variant.
    Gzip(Box<flate2::write::GzEncoder<BufWriter<Box<dyn Write + Send>>>>),
}

impl Write for SinkWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            SinkWriter::Plain(writer) => writer.write(buf),
            SinkWriter::Gzip(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            SinkWriter::Plain(writer) => writer.flush(),
            SinkWriter::Gzip(writer) => writer.flush(),
        }
    }
}

impl SinkWriter {
    /// Finalize the stream: write the gzip trailer (if any) and flush,
    /// propagating any error.
    fn finish(self) -> Result<()> {
        match self {
            SinkWriter::Plain(mut writer) => writer.flush().context("failed to flush output"),
            SinkWriter::Gzip(writer) => {
                let mut inner = (*writer)
                    .finish()
                    .context("failed to finish the gzip stream")?;
                inner.flush().context("failed to flush output")
            }
        }
    }
}

/// One open output sink, holding the concrete noodles writer for the chosen
/// format.
enum Sink {
    /// FASTQ, written as [`encode_fastq`] encodes it.
    Fastq(SinkWriter),
    /// FASTA.
    Fasta(fasta::io::Writer<SinkWriter>),
    /// Text SAM, plus the header written up front and reused to encode each
    /// record.
    Sam {
        /// The record writer.
        writer: sam::io::Writer<SinkWriter>,
        /// The header.
        header: sam::Header,
    },
    /// BGZF BAM where the bam writer wraps its own BGZF encoder, plus the
    /// header.
    BamInline {
        /// The record writer (its inner `bgzf::io::Writer` does the
        /// compression).
        writer: bam::io::Writer<bgzf::io::Writer<SinkWriter>>,
        /// The header.
        header: sam::Header,
    },
    /// Reference-free unmapped CRAM, plus the header (also needed to finalize
    /// the EOF container).
    Cram {
        /// The record writer.
        writer: cram::io::Writer<SinkWriter>,
        /// The header.
        header: sam::Header,
    },
}

/// Writes assembled output records to a single sink (a file, or stdout). Call
/// [`OutputWriter::finish`] to finalize; dropping without finishing may leave a
/// truncated file (errors unreported).
pub struct OutputWriter {
    sink: Sink,
}

impl OutputWriter {
    /// Create the output sink at `out` (or stdout when `out` is `None`) for the
    /// resolved `format`, creating missing parent directories. `compression`
    /// (0-9) sets the gzip level for FASTX.gz. SAM/BAM get a minimal `@HD`-only
    /// header; use [`OutputWriter::create_with_header`] to supply `@RG` lines
    /// for a demux fan-out.
    pub fn create(out: Option<&Path>, format: OutputFormat, compression: u8) -> Result<Self> {
        Self::create_with_header(out, format, compression, default_header())
    }

    /// Create the output sink with a caller-supplied SAM/BAM `header` (carrying
    /// `@RG` lines for a fan-out target). FASTX output has no header, so
    /// `header` is ignored there.
    pub fn create_with_header(
        out: Option<&Path>,
        format: OutputFormat,
        compression: u8,
        header: sam::Header,
    ) -> Result<Self> {
        let raw: Box<dyn Write + Send> = match out {
            Some(path) => Box::new(create_output_file(path)?),
            None => Box::new(std::io::stdout()),
        };
        let buffered = BufWriter::new(raw);

        let sink = match format {
            OutputFormat::Fastq { gzip } => Sink::Fastq(sink_writer(buffered, gzip, compression)),
            OutputFormat::Fasta { gzip } => Sink::Fasta(fasta::io::Writer::new(sink_writer(
                buffered,
                gzip,
                compression,
            ))),
            OutputFormat::Sam => {
                let mut writer = sam::io::Writer::new(SinkWriter::Plain(buffered));
                writer
                    .write_header(&header)
                    .context("failed to write SAM header")?;
                Sink::Sam { writer, header }
            }
            OutputFormat::Bam => {
                let mut writer =
                    bam::io::Writer::from(bgzf_writer(SinkWriter::Plain(buffered), compression)?);
                writer
                    .write_header(&header)
                    .context("failed to write BAM header")?;
                Sink::BamInline { writer, header }
            }
            OutputFormat::Cram => {
                // Default (empty) reference repository: the output is
                // reference-free unmapped CRAM.
                let mut writer = cram::io::Writer::new(SinkWriter::Plain(buffered));
                writer
                    .write_header(&header)
                    .context("failed to write CRAM header")?;
                Sink::Cram { writer, header }
            }
        };
        Ok(Self { sink })
    }

    /// Write one logical record's output records. A multi-record fragment emits
    /// SAM/BAM paired-segment records (or interleaved FASTX); SAM/BAM allows at
    /// most two (a pair).
    pub fn write_fragment(&mut self, reads: &[OutputRead]) -> Result<()> {
        match &mut self.sink {
            Sink::Fastq(writer) => {
                let mut bytes = Vec::new();
                for read in reads {
                    encode_fastq(&mut bytes, read)?;
                }
                writer
                    .write_all(&bytes)
                    .context("failed to write FASTQ record")?;
            }
            Sink::Fasta(writer) => {
                for read in reads {
                    write_fasta(writer, read)?;
                }
            }
            Sink::Sam { writer, header } => write_alignments(writer, header, reads)?,
            Sink::BamInline { writer, header } => write_alignments(writer, header, reads)?,
            Sink::Cram { writer, header } => write_alignments(writer, header, reads)?,
        }
        Ok(())
    }

    /// Append a fragment's already-encoded record bytes (produced by
    /// [`encode_fragment`] for this sink's format, e.g. on a worker thread),
    /// after the header written at open. The bytes are injected into the sink's
    /// underlying stream so they pass through its compression stage (gzip or
    /// BGZF) exactly as if the records had been encoded here. CRAM is not
    /// supported (its records are not independent per-fragment byte slices);
    /// callers route CRAM through `write_fragment`.
    pub fn write_encoded(&mut self, bytes: &[u8]) -> Result<()> {
        let sink: &mut dyn Write = match &mut self.sink {
            Sink::Fastq(writer) => writer,
            Sink::Fasta(writer) => writer.get_mut(),
            Sink::Sam { writer, .. } => writer.get_mut(),
            // get_mut yields the inner `bgzf::io::Writer`, so the raw record
            // bytes are BGZF-compressed here exactly as the bam writer would
            // have compressed them.
            Sink::BamInline { writer, .. } => writer.get_mut(),
            Sink::Cram { .. } => {
                bail!("internal error: CRAM cannot accept pre-encoded fragment bytes")
            }
        };
        sink.write_all(bytes)
            .context("failed to write pre-encoded record bytes")
    }

    /// Finalize the output: flush buffers, write the gzip trailer / BGZF EOF
    /// block / CRAM EOF container, and propagate any error (so a truncated file
    /// is reported rather than silently produced).
    pub fn finish(self) -> Result<()> {
        match self.sink {
            Sink::Fastq(writer) => writer.finish(),
            Sink::Fasta(writer) => writer.into_inner().finish(),
            Sink::Sam { writer, .. } => writer.into_inner().finish(),
            Sink::BamInline { writer, .. } => {
                let sink = writer
                    .into_inner()
                    .finish()
                    .context("failed to finalize BGZF (BAM EOF block)")?;
                sink.finish()
            }
            Sink::Cram { mut writer, header } => {
                writer
                    .try_finish(&header)
                    .context("failed to finalize the CRAM EOF container")?;
                writer.into_inner().finish()
            }
        }
    }
}

/// A BGZF writer over `inner` at the `--compression` level, which clap validates to 0-9.
fn bgzf_writer<W: Write>(inner: W, compression: u8) -> Result<bgzf::io::Writer<W>> {
    let level = bgzf::io::writer::CompressionLevel::try_from(compression)
        .map_err(|e| anyhow!("invalid --compression {compression} for BGZF: {e}"))?;
    Ok(bgzf::io::writer::Builder::default()
        .set_compression_level(level)
        .build_from_writer(inner))
}

/// Output files written beside their destinations under temporary names, renamed into place by
/// [`StagedOutputs::commit`]. Dropped uncommitted, it removes the temporary files and the
/// directories it created, so a failed run never changes a file that was already there.
#[derive(Debug)]
pub(crate) struct StagedOutputs {
    id: u64,
    files: HashMap<PathBuf, StagedFile>,
    dirs: Vec<PathBuf>,
    committed: bool,
}

impl Default for StagedOutputs {
    fn default() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        Self {
            id: NEXT.fetch_add(1, Ordering::Relaxed),
            files: HashMap::new(),
            dirs: Vec::new(),
            committed: false,
        }
    }
}

/// The temporary files and created directories of every live [`StagedOutputs`], by its id.
type Registry = BTreeMap<u64, (Vec<PathBuf>, Vec<PathBuf>)>;

static REGISTRY: Mutex<Registry> = Mutex::new(BTreeMap::new());

fn registry() -> MutexGuard<'static, Registry> {
    REGISTRY.lock().unwrap_or_else(PoisonError::into_inner)
}

/// Remove `files`, then the directories in `dirs` that are left empty, deepest first.
fn remove_staged(files: &[PathBuf], dirs: &[PathBuf]) {
    for file in files {
        let _ = std::fs::remove_file(file);
    }
    let mut dirs = dirs.to_vec();
    dirs.sort_by_key(|dir| std::cmp::Reverse(dir.components().count()));
    for dir in &dirs {
        let _ = std::fs::remove_dir(dir);
    }
}

/// On SIGINT, SIGTERM or SIGHUP, remove staged outputs and exit 128 plus the signal number.
pub fn clean_up_on_termination() -> io::Result<()> {
    use signal_hook::consts::{SIGHUP, SIGINT, SIGTERM};
    let armed = Arc::new(AtomicBool::new(false));
    for signal in [SIGINT, SIGTERM, SIGHUP] {
        signal_hook::flag::register_conditional_shutdown(signal, 128 + signal, Arc::clone(&armed))?;
        signal_hook::flag::register(signal, Arc::clone(&armed))?;
    }
    let mut signals = signal_hook::iterator::Signals::new([SIGINT, SIGTERM, SIGHUP])?;
    std::thread::Builder::new()
        .name("unmux-signals".to_string())
        .spawn(move || {
            if let Some(signal) = signals.forever().next() {
                let registry = registry();
                for (files, dirs) in registry.values() {
                    remove_staged(files, dirs);
                }
                std::process::exit(128 + signal);
            }
        })?;
    Ok(())
}

/// Where one staged output is written, and the file it replaces on commit.
#[derive(Debug)]
struct StagedFile {
    destination: PathBuf,
    temporary: PathBuf,
}

impl StagedOutputs {
    /// Stage every path in `paths`, noting its missing parent directories. A path that exists but
    /// is not a regular file, such as a device or a FIFO, is written in place; a symbolic link,
    /// dangling or not, is staged beside the file it points to.
    pub(crate) fn stage(&mut self, paths: impl IntoIterator<Item = PathBuf>) {
        for path in paths {
            if self.files.contains_key(&path) {
                continue;
            }
            let target = link_target(&path);
            let destination = match std::fs::metadata(&target) {
                Ok(metadata) if !metadata.is_file() => continue,
                Ok(_) => std::fs::canonicalize(&target).unwrap_or(target),
                Err(_) => target,
            };
            let Some(name) = destination.file_name() else {
                continue;
            };
            let temporary = destination.with_file_name(format!(
                ".{}.unmux-{}.tmp",
                name.to_string_lossy(),
                std::process::id()
            ));
            for dir in destination.ancestors().skip(1) {
                if dir.as_os_str().is_empty() || dir.symlink_metadata().is_ok() {
                    break;
                }
                if !self.dirs.iter().any(|known| known == dir) {
                    self.dirs.push(dir.to_path_buf());
                }
            }
            let file = StagedFile {
                destination,
                temporary,
            };
            self.files.insert(path, file);
        }
        let temporaries = self.files.values().map(|file| file.temporary.clone());
        registry().insert(self.id, (temporaries.collect(), self.dirs.clone()));
    }

    /// The file to write for output `path`: its temporary file when staged, else `path` itself.
    pub(crate) fn path<'a>(&'a self, path: &'a Path) -> &'a Path {
        self.files
            .get(path)
            .map_or(path, |file| file.temporary.as_path())
    }

    /// The temporary file of every staged output, by output path.
    pub(crate) fn temporaries(&self) -> HashMap<PathBuf, PathBuf> {
        self.files
            .iter()
            .map(|(path, file)| (path.clone(), file.temporary.clone()))
            .collect()
    }

    /// Rename every written temporary file over its destination.
    pub(crate) fn commit(&mut self) -> Result<()> {
        let mut registry = registry();
        for file in self.files.values() {
            if file.temporary.symlink_metadata().is_err() {
                continue;
            }
            replace(&file.temporary, &file.destination).with_context(|| {
                format!(
                    "failed to move the finished output into place: {}",
                    file.destination.display()
                )
            })?;
        }
        self.committed = true;
        registry.remove(&self.id);
        Ok(())
    }
}

impl Drop for StagedOutputs {
    fn drop(&mut self) {
        let mut registry = registry();
        if !self.committed {
            let temporaries: Vec<PathBuf> = self
                .files
                .values()
                .map(|file| file.temporary.clone())
                .collect();
            remove_staged(&temporaries, &self.dirs);
        }
        registry.remove(&self.id);
    }
}

/// The file `path` finally points to, following symbolic links even when the last one dangles.
fn link_target(path: &Path) -> PathBuf {
    let mut target = path.to_path_buf();
    for _ in 0..40 {
        match std::fs::read_link(&target) {
            Ok(next) => target = target.parent().unwrap_or(Path::new("")).join(next),
            Err(_) => return target,
        }
    }
    path.to_path_buf()
}

/// Move finished file `temporary` onto `destination`, keeping an existing destination's mode,
/// owner and hard links: renamed when they carry over, else copied into the destination in place.
fn replace(temporary: &Path, destination: &Path) -> io::Result<()> {
    let Ok(existing) = std::fs::metadata(destination) else {
        return std::fs::rename(temporary, destination);
    };
    if carries_over(&existing, temporary) {
        std::fs::set_permissions(temporary, existing.permissions())?;
        return std::fs::rename(temporary, destination);
    }
    let mut out = std::fs::OpenOptions::new()
        .write(true)
        .truncate(true)
        .open(destination)?;
    io::copy(&mut File::open(temporary)?, &mut out)?;
    out.sync_all()?;
    std::fs::remove_file(temporary)
}

/// Whether renaming `temporary` over a file with `existing` metadata keeps its links and owner.
#[cfg(unix)]
fn carries_over(existing: &std::fs::Metadata, temporary: &Path) -> bool {
    use std::os::unix::fs::MetadataExt;
    let Ok(staged) = std::fs::metadata(temporary) else {
        return false;
    };
    existing.nlink() == 1
        && existing.uid() == staged.uid()
        && (existing.gid() == staged.gid()
            || std::os::unix::fs::chown(temporary, None, Some(existing.gid())).is_ok())
}

#[cfg(not(unix))]
fn carries_over(_: &std::fs::Metadata, _: &Path) -> bool {
    true
}

/// A file's device and inode, which two names of one file share.
pub(crate) type FileId = (u64, u64);

/// The device and inode of a file, where the platform reports them.
#[cfg(unix)]
pub(crate) fn file_id(metadata: &std::fs::Metadata) -> Option<FileId> {
    use std::os::unix::fs::MetadataExt;
    Some((metadata.dev(), metadata.ino()))
}

#[cfg(not(unix))]
pub(crate) fn file_id(_: &std::fs::Metadata) -> Option<FileId> {
    None
}

/// Refuse two outputs whose staged files were created as one file, whatever the reason.
pub(crate) fn refuse_shared_files<'a>(
    staged: &StagedOutputs,
    outputs: impl IntoIterator<Item = &'a PathBuf>,
) -> Result<()> {
    let mut seen: HashMap<FileId, &Path> = HashMap::new();
    for path in outputs {
        let Ok(metadata) = std::fs::metadata(staged.path(path)) else {
            continue;
        };
        let Some(id) = file_id(&metadata).filter(|_| metadata.is_file()) else {
            continue;
        };
        match seen.get(&id) {
            Some(first) if *first != path.as_path() => bail!(
                "outputs {first:?} and {path:?} are one file on this file system; give them distinct paths"
            ),
            _ => {
                seen.insert(id, path);
            }
        }
    }
    Ok(())
}

/// Create an output file, making any missing parent directories first. Shared
/// by every output so each creates its file (and dirs) identically.
pub(crate) fn create_output_file(path: &Path) -> Result<File> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).with_context(|| {
                format!("failed to create output directory: {}", parent.display())
            })?;
        }
    }
    File::create(path).with_context(|| format!("failed to create output file: {}", path.display()))
}

/// Wrap a buffered stream in a gzip encoder when `gzip` is set, else leave it
/// plain.
fn sink_writer(
    buffered: BufWriter<Box<dyn Write + Send>>,
    gzip: bool,
    compression: u8,
) -> SinkWriter {
    if gzip {
        SinkWriter::Gzip(Box::new(flate2::write::GzEncoder::new(
            buffered,
            flate2::Compression::new(u32::from(compression)),
        )))
    } else {
        SinkWriter::Plain(buffered)
    }
}

/// A minimal `@HD`-bearing header. noodles only writes the `@HD` line when the
/// header map is present, so a default `sam::Header` (no `@HD`) would produce a
/// headerless, unreadable SAM.
pub(crate) fn default_header() -> sam::Header {
    use noodles::sam::header::record::value::{map, Map};
    sam::Header::builder()
        .set_header(Map::<map::Header>::default())
        .build()
}

/// Append `read` to `out` as one FASTQ record, byte for byte as the noodles
/// writer writes it: the qualities re-encoded to ASCII Phred+33, and input
/// tags carried in the read-name comment (samtools `fastq -T` convention).
fn encode_fastq(out: &mut Vec<u8>, read: &OutputRead) -> Result<()> {
    let quals = read
        .quals
        .context("cannot write FASTQ output for a record with no qualities (a FASTA input, or QUAL `*` in SAM/BAM/CRAM)")?;
    out.push(b'@');
    out.extend_from_slice(read.name);
    if let Some(tags) = read.tags.filter(|tags| !tags.is_empty()) {
        out.push(b' ');
        write_tags(out, tags).context("failed to write FASTQ record")?;
    }
    out.push(b'\n');
    out.extend_from_slice(read.bases);
    out.extend_from_slice(b"\n+\n");
    out.extend(quals.iter().map(|&q| q.saturating_add(PHRED_OFFSET)));
    out.push(b'\n');
    Ok(())
}

/// Write one FASTA record (no qualities), carrying input tags in the
/// description comment.
fn write_fasta<W: Write>(writer: &mut fasta::io::Writer<W>, read: &OutputRead) -> Result<()> {
    let comment = read
        .tags
        .map(render_tags)
        .transpose()
        .context("failed to write FASTA record")?
        .filter(|tags| !tags.is_empty())
        .map(Into::into);
    let record = fasta::Record::new(
        fasta::record::Definition::new(read.name.to_vec(), comment),
        fasta::record::Sequence::from(read.bases.to_vec()),
    );
    writer
        .write_record(&record)
        .context("failed to write FASTA record")
}

/// Render carried SAM tags as a tab-separated `TAG:TYPE:VALUE` comment (the
/// samtools `fastq -T` convention) so FASTX output preserves them and
/// round-trips via `samtools import` / `bwa -C`.
fn render_tags(data: &Data) -> io::Result<Vec<u8>> {
    let mut comment = Vec::new();
    write_tags(&mut comment, data)?;
    Ok(comment)
}

/// Append the comment [`render_tags`] renders for `data` to `out`.
fn write_tags(out: &mut Vec<u8>, data: &Data) -> io::Result<()> {
    for (i, (tag, value)) in data.iter().enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        let bytes: [u8; 2] = tag.into();
        let key = String::from_utf8_lossy(&bytes);
        match value {
            Value::Character(c) => write!(out, "{key}:A:{}", *c as char),
            Value::Int8(n) => write!(out, "{key}:i:{n}"),
            Value::UInt8(n) => write!(out, "{key}:i:{n}"),
            Value::Int16(n) => write!(out, "{key}:i:{n}"),
            Value::UInt16(n) => write!(out, "{key}:i:{n}"),
            Value::Int32(n) => write!(out, "{key}:i:{n}"),
            Value::UInt32(n) => write!(out, "{key}:i:{n}"),
            Value::Float(x) => write!(out, "{key}:f:{x}"),
            Value::String(s) => write!(out, "{key}:Z:{s}"),
            Value::Hex(s) => write!(out, "{key}:H:{s}"),
            Value::Array(array) => write!(out, "{key}:B:{}", render_array(array)),
        }?;
    }
    Ok(())
}

/// Render a `B` array tag value as `subtype,v1,v2,...` (SAM subtype codes
/// c/C/s/S/i/I/f).
fn render_array(array: &Array) -> String {
    fn join<T: std::fmt::Display>(subtype: char, values: &[T]) -> String {
        let mut rendered = subtype.to_string();
        for value in values {
            rendered.push(',');
            rendered.push_str(&value.to_string());
        }
        rendered
    }
    match array {
        Array::Int8(values) => join('c', values),
        Array::UInt8(values) => join('C', values),
        Array::Int16(values) => join('s', values),
        Array::UInt16(values) => join('S', values),
        Array::Int32(values) => join('i', values),
        Array::UInt32(values) => join('I', values),
        Array::Float(values) => join('f', values),
    }
}

/// The quality that BAM and CRAM store per base for a record without qualities.
const MISSING_QUALITY: u8 = 0xFF;

/// The longest read name SAM and BAM allow.
const MAX_NAME_LENGTH: usize = 254;

trait AlignmentSink: AlignmentWrite {
    const NEEDS_U_STORED_AS_T: bool = false;
    /// Adjust `record` so this writer stores it as BAM output does, or refuse what BAM refuses.
    fn prepare(_record: &mut std::borrow::Cow<'_, RecordBuf>) -> io::Result<()> {
        Ok(())
    }
}

impl<W: Write> AlignmentSink for sam::io::Writer<W> {}

/// BAM codes each base in 4 bits, with no code for U.
impl<W: Write> AlignmentSink for bam::io::Writer<W> {
    const NEEDS_U_STORED_AS_T: bool = true;
}

/// The noodles CRAM writer flags every record as storing a quality per base, even one with none.
impl<W: Write> AlignmentSink for cram::io::Writer<W> {
    fn prepare(record: &mut std::borrow::Cow<'_, RecordBuf>) -> io::Result<()> {
        check_name(record)?;
        // noodles ends a missing name with a second stop byte, which readers take as another name.
        if record.name().is_none() {
            *record.to_mut().name_mut() = Some("*".into());
        }
        fill_missing_qualities(record);
        Ok(())
    }
}

/// Refuse a read name that BAM output refuses, with its error; noodles checks none for CRAM.
fn check_name(record: &RecordBuf) -> io::Result<()> {
    let Some(name) = record.name() else {
        return Ok(());
    };
    let reason = if name.len() > MAX_NAME_LENGTH {
        // The length noodles reports counts the NUL that ends a name in BAM.
        let length = name.len() + 1;
        format!("invalid length: expected 1 <= n <= {MAX_NAME_LENGTH}, got {length}")
    } else if name.is_empty()
        || name == "*"
        || name.iter().any(|&b| !b.is_ascii_graphic() || b == b'@')
    {
        format!("invalid name: {name}")
    } else {
        return Ok(());
    };
    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        format!("invalid name: {reason}"),
    ))
}

/// Write a fragment's records as SAM/BAM alignment records, setting
/// paired-segment flags for a multi-record fragment. SAM/BAM holds at most two
/// records per template (a pair).
fn write_alignments<W: AlignmentSink>(
    writer: &mut W,
    header: &sam::Header,
    reads: &[OutputRead],
) -> Result<()> {
    if reads.len() > 2 {
        bail!(
            "SAM/BAM/CRAM output holds at most two reads per record (a pair); got {}",
            reads.len()
        );
    }
    for (index, read) in reads.iter().enumerate() {
        let mut record = match read.alignment {
            Some(alignment) => alignment.to_record_buf()?,
            None => {
                std::borrow::Cow::Owned(build_record(read, record_flags(read, index, reads.len())))
            }
        };
        if W::NEEDS_U_STORED_AS_T {
            store_u_as_t(&mut record);
        }
        W::prepare(&mut record)
            .and_then(|()| writer.write_alignment_record(header, record.as_ref()))
            .context("failed to write SAM/BAM/CRAM record")?;
    }
    Ok(())
}

fn fill_missing_qualities(record: &mut std::borrow::Cow<'_, RecordBuf>) {
    let bases = record.sequence().len();
    if bases > 0 && record.quality_scores().is_empty() {
        *record.to_mut().quality_scores_mut() = QualityScores::from(vec![MISSING_QUALITY; bases]);
    }
}

fn store_u_as_t(record: &mut std::borrow::Cow<'_, RecordBuf>) {
    if memchr::memchr2(b'U', b'u', record.sequence().as_ref()).is_some() {
        for base in record.to_mut().sequence_mut().as_mut() {
            if matches!(*base, b'U' | b'u') {
                *base = b'T';
            }
        }
    }
}

/// Build an unmapped alignment record, carrying tags through, adding the `RG`
/// data field for an assigned record, and omitting an empty name. This is the
/// per-record SAM/BAM encode step (the input to the writer's
/// `write_alignment_record`), exposed so it can be benchmarked in isolation.
pub fn build_record(read: &OutputRead, flags: Flags) -> RecordBuf {
    let mut builder = RecordBuf::builder()
        .set_flags(flags)
        .set_sequence(Sequence::from(read.bases.to_vec()));
    if !read.name.is_empty() {
        builder = builder.set_name(read.name.to_vec());
    }
    if let Some(quals) = read.quals {
        builder = builder.set_quality_scores(QualityScores::from(quals.to_vec()));
    }
    let mut data = read.tags.cloned().unwrap_or_default();
    if let Some(read_group) = read.read_group {
        data.insert(Tag::READ_GROUP, Value::from(read_group.to_string()));
    }
    if !data.is_empty() {
        builder = builder.set_data(data);
    }
    builder.build()
}

/// Encode one fragment's records into the bytes the matching sink writes for
/// `format`, via a throwaway in-memory writer, so they can be produced on a
/// worker thread and appended by [`OutputWriter::write_encoded`] in input order
/// (encode-on-workers). A raw BAM record bound for BAM output keeps its input
/// bytes (tag order and bin included) rather than being re-encoded as
/// [`OutputWriter::write_fragment`] would. `header` is used only for SAM/BAM
/// record encoding and ignored for FASTX. CRAM is NOT supported here (its records are
/// container/slice-coded, not independent per-fragment byte slices); route CRAM
/// through the record path. The throwaway BAM writer is built via `From` (no
/// inner BGZF), so it yields raw record blocks that the destination sink's own
/// BGZF stage then compresses, matching `write_fragment`.
pub fn encode_fragment(
    format: OutputFormat,
    header: &sam::Header,
    reads: &[OutputRead],
) -> Result<Vec<u8>> {
    Ok(match format {
        OutputFormat::Fastq { .. } => {
            let size = reads
                .iter()
                .map(|read| read.name.len() + 2 * read.bases.len());
            let mut bytes = Vec::with_capacity(size.sum::<usize>() + 64 * reads.len());
            for read in reads {
                encode_fastq(&mut bytes, read)?;
            }
            bytes
        }
        OutputFormat::Fasta { .. } => {
            let mut writer = fasta::io::Writer::new(Vec::new());
            for read in reads {
                write_fasta(&mut writer, read)?;
            }
            writer.into_inner()
        }
        // No `write_header`: only the record lines, since the destination sink
        // wrote the header at open.
        OutputFormat::Sam => {
            let mut writer = sam::io::Writer::new(Vec::new());
            write_alignments(&mut writer, header, reads)?;
            writer.into_inner()
        }
        OutputFormat::Bam => match raw_bam_records(header, reads) {
            Some(raw) => raw,
            None => {
                let mut writer = bam::io::Writer::from(Vec::new());
                write_alignments(&mut writer, header, reads)?;
                writer.into_inner()
            }
        },
        OutputFormat::Cram => {
            bail!("internal error: CRAM is not eligible for encode-on-workers")
        }
    })
}

/// The concatenated raw bytes of `reads` when every one is a raw BAM alignment
/// whose reference ids `header` declares, which BAM output writes back verbatim.
fn raw_bam_records(header: &sam::Header, reads: &[OutputRead]) -> Option<Vec<u8>> {
    if reads.len() > 2 {
        return None;
    }
    let references = header.reference_sequences().len();
    let declared = |raw: &[u8], at: usize| {
        let id = raw
            .get(at..at + 4)
            .map(|b| i32::from_le_bytes([b[0], b[1], b[2], b[3]]));
        id == Some(-1)
            || id
                .and_then(|id| usize::try_from(id).ok())
                .is_some_and(|id| id < references)
    };
    let mut bytes = Vec::new();
    for read in reads {
        match read.alignment {
            Some(Alignment::RawBam(raw)) if declared(raw, 4) && declared(raw, 24) => {
                bytes.extend_from_slice(raw);
            }
            _ => return None,
        }
    }
    Some(bytes)
}

/// The flags for `read` at `index` of `total` output records.
fn record_flags(read: &OutputRead, index: usize, total: usize) -> Flags {
    let flags = segment_flags(index, total);
    if read.qc_fail {
        flags | Flags::QC_FAIL
    } else {
        flags
    }
}

/// The flags for output record `index` of `total`: a lone record is just
/// unmapped; a pair gets the segmented + first/last-segment + mate-unmapped
/// bits so downstream tools can re-pair it.
fn segment_flags(index: usize, total: usize) -> Flags {
    if total <= 1 {
        Flags::UNMAPPED
    } else {
        let mut flags = Flags::UNMAPPED | Flags::MATE_UNMAPPED | Flags::SEGMENTED;
        flags |= if index == 0 {
            Flags::FIRST_SEGMENT
        } else {
            Flags::LAST_SEGMENT
        };
        flags
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(unix)]
    #[test]
    fn test_refuse_shared_files_when_staged_outputs_are_one_file() {
        let dir = tempfile::tempdir().unwrap();
        let (a, b, c) = (
            dir.path().join("o/a.fq"),
            dir.path().join("o/b.fq"),
            dir.path().join("o/c.fq"),
        );
        let mut staged = StagedOutputs::default();
        staged.stage([a.clone(), b.clone(), c.clone()]);
        create_output_file(staged.path(&a)).unwrap();
        create_output_file(staged.path(&c)).unwrap();
        assert!(refuse_shared_files(&staged, [&a, &a, &c]).is_ok());
        std::fs::hard_link(staged.path(&a), staged.path(&b)).unwrap();
        let error = refuse_shared_files(&staged, [&a, &b, &c])
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("are one file on this file system"),
            "{error}"
        );
        drop(staged);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn test_staged_outputs_leave_existing_files_untouched_on_failure() {
        let dir = tempfile::tempdir().unwrap();
        let existing = dir.path().join("existing.bam");
        std::fs::write(&existing, b"keep").unwrap();
        let created = dir.path().join("new/nested/out.bam");
        let mut staged = StagedOutputs::default();
        staged.stage([existing.clone(), created.clone()]);
        for path in [&existing, &created] {
            assert_ne!(staged.path(path), path.as_path());
            create_output_file(staged.path(path)).unwrap();
            std::fs::write(staged.path(path), b"partial").unwrap();
        }
        let temporaries = staged.temporaries();
        drop(staged);
        assert_eq!(std::fs::read(&existing).unwrap(), b"keep");
        assert!(!created.exists());
        assert!(!dir.path().join("new").exists());
        assert!(temporaries.values().all(|temporary| !temporary.exists()));
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 1);
    }

    #[test]
    fn test_staged_outputs_replace_their_destinations_on_commit() {
        let dir = tempfile::tempdir().unwrap();
        let existing = dir.path().join("existing.tsv");
        std::fs::write(&existing, b"old").unwrap();
        let link = dir.path().join("link.tsv");
        std::os::unix::fs::symlink(&existing, &link).unwrap();
        let created = dir.path().join("out/out.tsv");
        let unwritten = dir.path().join("unwritten.tsv");
        let mut staged = StagedOutputs::default();
        staged.stage([link.clone(), created.clone(), unwritten.clone()]);
        for (path, text) in [(&link, "new"), (&created, "created")] {
            create_output_file(staged.path(path)).unwrap();
            std::fs::write(staged.path(path), text).unwrap();
        }
        staged.commit().unwrap();
        drop(staged);
        assert!(link.symlink_metadata().unwrap().file_type().is_symlink());
        assert_eq!(std::fs::read_to_string(&existing).unwrap(), "new");
        assert_eq!(std::fs::read_to_string(&created).unwrap(), "created");
        assert!(!unwritten.exists());
        let mut names: Vec<_> = std::fs::read_dir(dir.path())
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .collect();
        names.sort();
        assert_eq!(names, ["existing.tsv", "link.tsv", "out"]);
    }

    /// Stage `path`, write `text` to its temporary file and commit it.
    fn stage_and_commit(path: &Path, text: &str) {
        let mut staged = StagedOutputs::default();
        staged.stage([path.to_path_buf()]);
        assert_ne!(staged.path(path), path);
        create_output_file(staged.path(path)).unwrap();
        std::fs::write(staged.path(path), text).unwrap();
        staged.commit().unwrap();
    }

    #[test]
    fn test_staged_outputs_keep_the_mode_of_the_file_they_replace() {
        use std::os::unix::fs::PermissionsExt;
        let dir = tempfile::tempdir().unwrap();
        let existing = dir.path().join("existing.sam");
        std::fs::write(&existing, b"old").unwrap();
        std::fs::set_permissions(&existing, std::fs::Permissions::from_mode(0o640)).unwrap();
        stage_and_commit(&existing, "new");
        assert_eq!(std::fs::read_to_string(&existing).unwrap(), "new");
        let mode = std::fs::metadata(&existing).unwrap().permissions().mode();
        assert_eq!(mode & 0o777, 0o640);
    }

    #[test]
    fn test_staged_outputs_rewrite_a_hard_linked_file_in_place() {
        use std::os::unix::fs::MetadataExt;
        let dir = tempfile::tempdir().unwrap();
        let existing = dir.path().join("existing.sam");
        let link = dir.path().join("link.sam");
        std::fs::write(&existing, b"old").unwrap();
        std::fs::hard_link(&existing, &link).unwrap();
        let inode = std::fs::metadata(&existing).unwrap().ino();
        stage_and_commit(&existing, "new");
        assert_eq!(std::fs::read_to_string(&link).unwrap(), "new");
        assert_eq!(std::fs::metadata(&existing).unwrap().ino(), inode);
        assert_eq!(std::fs::metadata(&existing).unwrap().nlink(), 2);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 2);
    }

    #[test]
    fn test_staged_outputs_write_through_a_dangling_symbolic_link() {
        let dir = tempfile::tempdir().unwrap();
        let link = dir.path().join("link.sam");
        std::os::unix::fs::symlink("data/target.sam", &link).unwrap();
        stage_and_commit(&link, "new");
        assert!(link.symlink_metadata().unwrap().file_type().is_symlink());
        assert_eq!(
            std::fs::read_to_string(dir.path().join("data/target.sam")).unwrap(),
            "new"
        );
        assert_eq!(std::fs::read_to_string(&link).unwrap(), "new");
    }

    #[test]
    fn test_link_target_follows_a_chain_of_links() {
        let dir = tempfile::tempdir().unwrap();
        let first = dir.path().join("first");
        std::os::unix::fs::symlink("second", &first).unwrap();
        std::os::unix::fs::symlink("sub/third", dir.path().join("second")).unwrap();
        assert_eq!(link_target(&first), dir.path().join("sub/third"));
        let looped = dir.path().join("looped");
        std::os::unix::fs::symlink("looped", &looped).unwrap();
        assert_eq!(link_target(&looped), looped);
    }

    #[test]
    fn test_staged_outputs_write_devices_in_place() {
        let mut staged = StagedOutputs::default();
        staged.stage([PathBuf::from("/dev/null")]);
        assert_eq!(staged.path(Path::new("/dev/null")), Path::new("/dev/null"));
        staged.commit().unwrap();
    }

    #[test]
    fn test_format_from_extension() {
        let f = |p: &str| format_from_extension(Path::new(p)).unwrap().unwrap();
        assert_eq!(f("out.fq"), OutputFormat::Fastq { gzip: false });
        assert_eq!(f("out.fastq.gz"), OutputFormat::Fastq { gzip: true });
        assert_eq!(f("out.fa"), OutputFormat::Fasta { gzip: false });
        assert_eq!(f("out.fasta.gz"), OutputFormat::Fasta { gzip: true });
        assert_eq!(f("out.sam"), OutputFormat::Sam);
        assert_eq!(f("out.bam"), OutputFormat::Bam);
        assert_eq!(f("out.cram"), OutputFormat::Cram);
    }

    #[test]
    fn test_format_from_extension_errors() {
        // .sam.gz / .cram.gz are rejected rather than silently mis-handled.
        assert!(format_from_extension(Path::new("out.sam.gz")).is_err());
        assert!(format_from_extension(Path::new("out.cram.gz")).is_err());
        // An unrecognized extension is not an error: it returns None so the
        // caller can mirror the input format (an extensionless or
        // unknown-suffix --out behaves like stdout).
        assert_eq!(format_from_extension(Path::new("out.txt")).unwrap(), None);
        assert_eq!(format_from_extension(Path::new("results")).unwrap(), None);
    }

    #[test]
    fn test_output_format_extensionless_mirrors_input() {
        use SniffedFormat::{Fasta, Fastq};
        // A concrete extensionless --out mirrors the input format, like stdout,
        // instead of erroring.
        assert_eq!(
            output_format(Some(Path::new("results")), &[Fastq]).unwrap(),
            OutputFormat::Fastq { gzip: false }
        );
        assert_eq!(
            output_format(Some(Path::new("results")), &[Fasta]).unwrap(),
            OutputFormat::Fasta { gzip: false }
        );
        // A mixed-format input defaults to BAM, the same default stdout uses.
        assert_eq!(
            output_format(Some(Path::new("results")), &[Fastq, Fasta]).unwrap(),
            OutputFormat::Bam
        );
        // An explicit extension still wins over the input format.
        assert_eq!(
            output_format(Some(Path::new("out.bam")), &[Fastq]).unwrap(),
            OutputFormat::Bam
        );
        // .sam.gz / .cram.gz stay hard errors (recognized but unsupported),
        // never a silent mirror.
        assert!(output_format(Some(Path::new("out.sam.gz")), &[Fastq]).is_err());
        assert!(output_format(Some(Path::new("out.cram.gz")), &[Fastq]).is_err());
    }

    #[test]
    fn test_stdout_format_mirrors_single_and_defaults_blend() {
        use SniffedFormat::{Bam, Fasta, Fastq, Sam};
        assert_eq!(stdout_format(&[Fastq]), OutputFormat::Fastq { gzip: false });
        assert_eq!(
            stdout_format(&[Fasta, Fasta]),
            OutputFormat::Fasta { gzip: false }
        );
        assert_eq!(stdout_format(&[Sam]), OutputFormat::Sam);
        assert_eq!(stdout_format(&[Bam]), OutputFormat::Bam);
        // A multi-format (quality-bearing) blend defaults to uncompressed BAM.
        assert_eq!(stdout_format(&[Fastq, Bam]), OutputFormat::Bam);
        assert_eq!(stdout_format(&[]), OutputFormat::Bam);
    }

    #[test]
    fn test_record_flags_carry_qc_fail() {
        let flags = |qc_fail: bool, index: usize, total: usize| {
            let read = OutputRead {
                qc_fail,
                ..output_read(b"r1", b"ACGT", None, None)
            };
            u16::from(record_flags(&read, index, total))
        };
        assert_eq!((flags(false, 0, 1), flags(true, 0, 1)), (4, 516));
        assert_eq!((flags(false, 0, 2), flags(true, 0, 2)), (77, 589));
        assert_eq!((flags(false, 1, 2), flags(true, 1, 2)), (141, 653));
    }

    #[test]
    fn test_segment_flags() {
        assert_eq!(segment_flags(0, 1), Flags::UNMAPPED);
        let first = segment_flags(0, 2);
        assert!(first.is_segmented() && first.is_first_segment() && first.is_unmapped());
        let last = segment_flags(1, 2);
        assert!(last.is_segmented() && last.is_last_segment() && last.is_unmapped());
    }

    #[test]
    fn test_bam_carries_at_rg_header_and_record_rg() {
        use std::fs::File;

        use noodles::sam::header::record::value::map::{read_group::tag as rg, Header, ReadGroup};
        use noodles::sam::header::record::value::Map;

        // A header with one @RG (id dna01.lib01, SM dna01, LB lib01) and a read
        // tagged with it.
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("rg.bam");
        let header = sam::Header::builder()
            .set_header(Map::<Header>::default())
            .add_read_group(
                "dna01.lib01",
                Map::<ReadGroup>::builder()
                    .insert(rg::SAMPLE, "dna01")
                    .insert(rg::LIBRARY, "lib01")
                    .build()
                    .unwrap(),
            )
            .build();
        let mut writer =
            OutputWriter::create_with_header(Some(&out), OutputFormat::Bam, 5, header).unwrap();
        writer
            .write_fragment(&[OutputRead {
                name: b"r1",
                bases: b"ACGT",
                quals: Some(&[30, 30, 30, 30]),
                tags: None,
                read_group: Some("dna01.lib01"),
                qc_fail: false,
                alignment: None,
            }])
            .unwrap();
        writer.finish().unwrap();

        let mut reader = bam::io::Reader::new(File::open(&out).unwrap());
        let header = reader.read_header().unwrap();
        assert!(
            header.read_groups().contains_key(&b"dna01.lib01"[..]),
            "@RG line present"
        );
        let record = reader.record_bufs(&header).next().unwrap().unwrap();
        match record.data().get(&Tag::READ_GROUP) {
            Some(Value::String(rg)) => assert_eq!(rg.to_string(), "dna01.lib01"),
            other => panic!("RG missing or wrong type: {other:?}"),
        }
    }

    /// The encode-on-workers invariant: appending `encode_fragment` bytes via
    /// `write_encoded` must produce a byte-for-byte identical file to encoding
    /// inline via `write_fragment`, for every encode-eligible format
    /// (FASTQ/FASTA/SAM/BAM), single-end and paired.
    #[test]
    fn test_write_encoded_matches_write_fragment_byte_for_byte() {
        use std::fs;

        use noodles::sam::header::record::value::map::{
            read_group::tag as rg, Header as HeaderMap, ReadGroup,
        };
        use noodles::sam::header::record::value::Map;

        let header = sam::Header::builder()
            .set_header(Map::<HeaderMap>::default())
            .add_read_group(
                "dna01",
                Map::<ReadGroup>::builder()
                    .insert(rg::SAMPLE, "dna01")
                    .build()
                    .unwrap(),
            )
            .build();

        let q = [40u8; 6];
        // A single-end fragment and a paired fragment (exercises the segment
        // flags), both with quals and a read group.
        let frags: [Vec<OutputRead>; 2] = [
            vec![OutputRead {
                name: b"r1",
                bases: b"ACGTACGT",
                quals: Some(&[30, 31, 32, 33, 34, 35, 36, 37]),
                tags: None,
                read_group: Some("dna01"),
                qc_fail: false,
                alignment: None,
            }],
            vec![
                OutputRead {
                    name: b"r2",
                    bases: b"AAACCC",
                    quals: Some(&q),
                    tags: None,
                    read_group: Some("dna01"),
                    qc_fail: false,
                    alignment: None,
                },
                OutputRead {
                    name: b"r2",
                    bases: b"GGGTTT",
                    quals: Some(&q),
                    tags: None,
                    read_group: Some("dna01"),
                    qc_fail: false,
                    alignment: None,
                },
            ],
        ];

        for format in [
            OutputFormat::Fastq { gzip: false },
            OutputFormat::Fasta { gzip: false },
            OutputFormat::Sam,
            OutputFormat::Bam,
        ] {
            let dir = tempfile::tempdir().unwrap();
            let direct = dir.path().join("direct");
            let encoded = dir.path().join("encoded");

            let mut w =
                OutputWriter::create_with_header(Some(&direct), format, 5, header.clone()).unwrap();
            for f in &frags {
                w.write_fragment(f).unwrap();
            }
            w.finish().unwrap();

            let mut w = OutputWriter::create_with_header(Some(&encoded), format, 5, header.clone())
                .unwrap();
            for f in &frags {
                let bytes = encode_fragment(format, &header, f).unwrap();
                w.write_encoded(&bytes).unwrap();
            }
            w.finish().unwrap();

            assert_eq!(
                fs::read(&direct).unwrap(),
                fs::read(&encoded).unwrap(),
                "encode_fragment+write_encoded diverged from write_fragment for {format:?}"
            );
        }
    }

    fn render_one(tag: [u8; 2], value: Value) -> String {
        let mut data = Data::default();
        data.insert(Tag::from(tag), value);
        String::from_utf8(render_tags(&data).unwrap()).unwrap()
    }

    #[test]
    fn test_encode_fastq_matches_the_noodles_writer() {
        let mut tags = Data::default();
        tags.insert(Tag::from(*b"BC"), Value::String("ACGT-TTGA".into()));
        tags.insert(Tag::from(*b"XI"), Value::Int32(-7));
        tags.insert(Tag::from(*b"XF"), Value::Float(0.25));
        tags.insert(Tag::from(*b"XB"), Value::Array(Array::UInt8(vec![1, 2])));
        let empty = Data::default();
        let quals = [0u8, 2, 30, 41, 93];
        for (name, bases, tags) in [
            (&b"r1"[..], &b"ACGTN"[..], None),
            (b"r2/1", b"ACGTN", Some(&tags)),
            (b"r3", b"ACGTN", Some(&empty)),
            (b"", b"", None),
        ] {
            let read = OutputRead {
                name,
                bases,
                quals: Some(&quals[..bases.len()]),
                tags,
                read_group: None,
                qc_fail: false,
                alignment: None,
            };
            let mut expected = noodles::fastq::io::Writer::new(Vec::new());
            let comment = tags
                .map(|tags| render_tags(tags).unwrap())
                .unwrap_or_default();
            let record = noodles::fastq::Record::new(
                noodles::fastq::record::Definition::new(name.to_vec(), comment),
                bases.to_vec(),
                quals[..bases.len()]
                    .iter()
                    .map(|q| q + PHRED_OFFSET)
                    .collect::<Vec<u8>>(),
            );
            expected.write_record(&record).unwrap();
            let mut encoded = Vec::new();
            encode_fastq(&mut encoded, &read).unwrap();
            assert_eq!(encoded, expected.into_inner());
        }
    }

    #[test]
    fn test_render_tags_every_scalar_value_type() {
        // Each `Value` variant renders to its samtools `fastq -T` field. The
        // signed/unsigned integer widths all share the `i` type code; character,
        // float, string, and hex have their own.
        assert_eq!(render_one(*b"XA", Value::Character(b'Q')), "XA:A:Q");
        assert_eq!(render_one(*b"XB", Value::Int8(-5)), "XB:i:-5");
        assert_eq!(render_one(*b"XC", Value::UInt8(5)), "XC:i:5");
        assert_eq!(render_one(*b"XD", Value::Int16(-300)), "XD:i:-300");
        assert_eq!(render_one(*b"XE", Value::UInt16(300)), "XE:i:300");
        assert_eq!(render_one(*b"XF", Value::Int32(-70000)), "XF:i:-70000");
        assert_eq!(render_one(*b"XG", Value::UInt32(70000)), "XG:i:70000");
        assert_eq!(render_one(*b"XH", Value::Float(1.5)), "XH:f:1.5");
        assert_eq!(
            render_one(*b"XI", Value::String("ACGT".into())),
            "XI:Z:ACGT"
        );
        assert_eq!(render_one(*b"XJ", Value::Hex("CAFE".into())), "XJ:H:CAFE");
    }

    #[test]
    fn test_render_array_every_subtype() {
        // A `B` array renders as `subtype,v1,v2,...` with the SAM subtype code per
        // element width (c/C/s/S/i/I for the integers, f for float).
        assert_eq!(render_array(&Array::Int8(vec![-1, 2])), "c,-1,2");
        assert_eq!(render_array(&Array::UInt8(vec![1, 2])), "C,1,2");
        assert_eq!(render_array(&Array::Int16(vec![-3, 4])), "s,-3,4");
        assert_eq!(render_array(&Array::UInt16(vec![5, 6])), "S,5,6");
        assert_eq!(render_array(&Array::Int32(vec![-7, 8])), "i,-7,8");
        assert_eq!(render_array(&Array::UInt32(vec![9, 10])), "I,9,10");
        assert_eq!(render_array(&Array::Float(vec![1.5])), "f,1.5");
        // And an array value renders through `render_tags` with the `B` code.
        assert_eq!(
            render_one(*b"XK", Value::Array(Array::Int32(vec![1, 2, 3]))),
            "XK:B:i,1,2,3"
        );
    }

    #[test]
    fn test_gzip_fastq_roundtrip_is_valid_gzip() {
        // create() (which supplies the default header) + a gzipped FASTQ sink:
        // exercise the gzip encoder write/finish path and confirm the file is real
        // gzip that decodes back to the record.
        use std::io::Read;
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.fq.gz");
        let mut w =
            OutputWriter::create(Some(&out), OutputFormat::Fastq { gzip: true }, 6).unwrap();
        w.write_fragment(&[OutputRead {
            name: b"r1",
            bases: b"ACGT",
            quals: Some(&[30, 30, 30, 30]),
            tags: None,
            read_group: None,
            qc_fail: false,
            alignment: None,
        }])
        .unwrap();
        w.finish().unwrap();

        let bytes = std::fs::read(&out).unwrap();
        assert_eq!(&bytes[..2], &[0x1f, 0x8b], "gzip magic bytes");
        let mut text = String::new();
        flate2::read::GzDecoder::new(&bytes[..])
            .read_to_string(&mut text)
            .unwrap();
        // qual 30 re-encodes to Phred+33 '?' (30 + 33 = 63).
        assert_eq!(text, "@r1\nACGT\n+\n????\n", "decoded FASTQ: {text:?}");
    }

    #[test]
    fn test_write_encoded_rejects_cram() {
        // CRAM records are not independent per-fragment byte slices, so the
        // encode-on-workers append path refuses them.
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.cram");
        let mut w = OutputWriter::create(Some(&out), OutputFormat::Cram, 0).unwrap();
        let err = w.write_encoded(b"bytes").unwrap_err().to_string();
        assert!(err.contains("CRAM"), "{err}");
    }

    #[test]
    fn test_encode_fragment_rejects_cram() {
        let header = default_header();
        let err = encode_fragment(OutputFormat::Cram, &header, &[])
            .unwrap_err()
            .to_string();
        assert!(err.contains("CRAM"), "{err}");
    }

    /// Write `fragments` to a CRAM and read its records back as SAM text.
    fn cram_round_trip(fragments: &[Vec<OutputRead>]) -> String {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.cram");
        let mut writer = OutputWriter::create(Some(&out), OutputFormat::Cram, 0).unwrap();
        for fragment in fragments {
            writer.write_fragment(fragment).unwrap();
        }
        writer.finish().unwrap();

        let mut reader = cram::io::Reader::new(File::open(&out).unwrap());
        let header = reader.read_header().unwrap();
        let mut text = sam::io::Writer::new(Vec::new());
        for record in reader.records(&header) {
            text.write_alignment_record(&header, &record.unwrap())
                .unwrap();
        }
        String::from_utf8(text.into_inner()).unwrap()
    }

    fn output_read<'a>(
        name: &'a [u8],
        bases: &'a [u8],
        quals: Option<&'a [u8]>,
        tags: Option<&'a Data>,
    ) -> OutputRead<'a> {
        OutputRead {
            name,
            bases,
            quals,
            tags,
            read_group: None,
            qc_fail: false,
            alignment: None,
        }
    }

    #[test]
    fn test_cram_round_trips_records_without_qualities() {
        let mut tags = Data::default();
        tags.insert(Tag::new(b'X', b'Y'), Value::from("hello"));
        let fragments = [
            vec![output_read(b"r1", b"ACGTACGT", None, Some(&tags))],
            vec![output_read(b"r2", b"TTTT", None, None)],
        ];
        assert_eq!(
            cram_round_trip(&fragments),
            "r1\t4\t*\t0\t255\t*\t*\t0\t0\tACGTACGT\t*\tXY:Z:hello\n\
             r2\t4\t*\t0\t255\t*\t*\t0\t0\tTTTT\t*\n"
        );
    }

    #[test]
    fn test_cram_round_trips_records_with_and_without_qualities() {
        let mut tags = Data::default();
        tags.insert(Tag::new(b'X', b'Y'), Value::from("hello"));
        let fragments = [
            vec![output_read(b"r1", b"ACGTACGT", Some(&[30; 8]), Some(&tags))],
            vec![output_read(b"r2", b"GGGG", None, Some(&tags))],
            vec![
                output_read(b"r3", b"TTTTAA", Some(&[0, 1, 2, 3, 4, 5]), None),
                output_read(b"r3", b"CC", None, None),
            ],
            vec![
                output_read(b"r4", b"AAAA", None, None),
                output_read(b"r4", b"CCCC", Some(&[40; 4]), None),
            ],
        ];
        assert_eq!(
            cram_round_trip(&fragments),
            "r1\t4\t*\t0\t255\t*\t*\t0\t0\tACGTACGT\t????????\tXY:Z:hello\n\
             r2\t4\t*\t0\t255\t*\t*\t0\t0\tGGGG\t*\tXY:Z:hello\n\
             r3\t77\t*\t0\t255\t*\t*\t0\t0\tTTTTAA\t!\"#$%&\n\
             r3\t141\t*\t0\t255\t*\t*\t0\t0\tCC\t*\n\
             r4\t77\t*\t0\t255\t*\t*\t0\t0\tAAAA\t*\n\
             r4\t141\t*\t0\t255\t*\t*\t0\t0\tCCCC\tIIII\n"
        );
    }

    #[test]
    fn test_cram_round_trips_a_kept_alignment_without_qualities() {
        let text = b"@HD\tVN:1.6\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\tXY:Z:hello\n";
        let mut reader = sam::io::Reader::new(&text[..]);
        let header = reader.read_header().unwrap();
        let mut record = RecordBuf::default();
        reader.read_record_buf(&header, &mut record).unwrap();
        let alignment = Alignment::Record(record);
        let read = OutputRead {
            alignment: Some(&alignment),
            ..output_read(b"", b"", None, None)
        };
        assert_eq!(
            cram_round_trip(&[vec![read]]),
            "r1\t4\t*\t0\t255\t*\t*\t0\t0\tACGT\t*\tXY:Z:hello\n"
        );
    }

    /// Write `fragments` as BAM or CRAM and read back each record's name (`*` when it has none)
    /// and bases with noodles.
    fn names_and_bases(
        format: OutputFormat,
        fragments: &[Vec<OutputRead>],
    ) -> Vec<(String, String)> {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out");
        let mut writer = OutputWriter::create(Some(&out), format, 0).unwrap();
        for fragment in fragments {
            writer.write_fragment(fragment).unwrap();
        }
        writer.finish().unwrap();

        let file = File::open(&out).unwrap();
        let records: Vec<RecordBuf> = if format == OutputFormat::Cram {
            let mut reader = cram::io::Reader::new(file);
            let header = reader.read_header().unwrap();
            reader.records(&header).collect::<io::Result<_>>().unwrap()
        } else {
            let mut reader = bam::io::Reader::new(file);
            let header = reader.read_header().unwrap();
            reader
                .record_bufs(&header)
                .collect::<io::Result<_>>()
                .unwrap()
        };
        records
            .iter()
            .map(|record| {
                let name = record
                    .name()
                    .map_or("*".to_string(), |name| name.to_string());
                let bases = String::from_utf8(record.sequence().as_ref().to_vec()).unwrap();
                (name, bases)
            })
            .collect()
    }

    #[test]
    fn test_cram_keeps_each_name_with_its_record_around_nameless_records() {
        let text = b"@HD\tVN:1.6\n*\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n";
        let mut reader = sam::io::Reader::new(&text[..]);
        let header = reader.read_header().unwrap();
        let mut record = RecordBuf::default();
        reader.read_record_buf(&header, &mut record).unwrap();
        let kept = Alignment::Record(record);
        let quals = [30; 4];
        let fragments = [
            vec![output_read(b"n1", b"AAAA", Some(&quals), None)],
            vec![output_read(b"", b"CCCC", Some(&quals), None)],
            vec![OutputRead {
                alignment: Some(&kept),
                ..output_read(b"", b"", None, None)
            }],
            vec![output_read(b"n3", b"GGGG", Some(&quals), None)],
            vec![output_read(b"n4", b"TTTT", Some(&quals), None)],
        ];
        let expected = [
            ("n1", "AAAA"),
            ("*", "CCCC"),
            ("*", "ACGT"),
            ("n3", "GGGG"),
            ("n4", "TTTT"),
        ]
        .map(|(name, bases)| (name.to_string(), bases.to_string()));
        assert_eq!(names_and_bases(OutputFormat::Bam, &fragments), expected);
        assert_eq!(names_and_bases(OutputFormat::Cram, &fragments), expected);
    }

    #[test]
    fn test_cram_refuses_the_read_names_bam_refuses_with_the_same_error() {
        let write = |format: OutputFormat, name: &[u8]| {
            let dir = tempfile::tempdir().unwrap();
            let out = dir.path().join("out");
            let mut writer = OutputWriter::create(Some(&out), format, 0).unwrap();
            writer
                .write_fragment(&[output_read(name, b"ACGT", Some(&[30; 4]), None)])
                .map_err(|e| format!("{e:#}"))
        };
        let longest = [b'r'; 254];
        let too_long = [b'r'; 255];
        assert!(write(OutputFormat::Cram, &longest).is_ok());
        assert_eq!(
            write(OutputFormat::Cram, &too_long).unwrap_err(),
            "failed to write SAM/BAM/CRAM record: invalid name: invalid length: expected 1 <= n <= 254, got 256"
        );
        for name in [
            &longest[..],
            &too_long,
            &[b'r'; 300],
            b"r@1",
            b"r 1",
            b"*",
            "r\u{e9}".as_bytes(),
        ] {
            assert_eq!(
                write(OutputFormat::Cram, name),
                write(OutputFormat::Bam, name),
                "{}",
                String::from_utf8_lossy(name)
            );
        }
    }

    #[test]
    fn test_fill_missing_qualities_needs_bases_and_no_qualities() {
        use std::borrow::Cow;
        let record = |bases: &[u8], quals: Option<&[u8]>| {
            build_record(&output_read(b"r1", bases, quals, None), Flags::UNMAPPED)
        };
        let mut filled = Cow::Owned(record(b"ACGT", None));
        fill_missing_qualities(&mut filled);
        assert_eq!(filled.quality_scores().as_ref(), [MISSING_QUALITY; 4]);

        for kept in [
            record(b"ACGT", Some(&[1, 2, 3, 4])),
            record(b"", None),
            record(b"", Some(&[])),
        ] {
            let mut record = Cow::Borrowed(&kept);
            fill_missing_qualities(&mut record);
            assert!(matches!(record, Cow::Borrowed(_)));
        }
    }

    #[test]
    fn test_more_than_two_reads_per_fragment_is_rejected() {
        // SAM/BAM/CRAM holds at most two segments (a pair) per template.
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.sam");
        let mut w = OutputWriter::create(Some(&out), OutputFormat::Sam, 0).unwrap();
        let q = [30u8; 4];
        let reads: Vec<OutputRead> = (0..3)
            .map(|_| OutputRead {
                name: b"r",
                bases: b"ACGT",
                quals: Some(&q),
                tags: None,
                read_group: None,
                qc_fail: false,
                alignment: None,
            })
            .collect();
        let err = w.write_fragment(&reads).unwrap_err().to_string();
        assert!(err.contains("at most two"), "{err}");
    }

    #[test]
    fn test_sink_writer_write_and_flush_delegate() {
        // The `Write` impl delegates to the inner stream for the plain and gzip
        // variants.
        let buffer = || BufWriter::new(Box::new(Vec::new()) as Box<dyn Write + Send>);
        let mut plain = SinkWriter::Plain(buffer());
        assert_eq!(plain.write(b"abc").unwrap(), 3);
        plain.flush().unwrap();
        plain.finish().unwrap();

        let mut gzip = sink_writer(buffer(), true, 6);
        assert_eq!(gzip.write(b"abc").unwrap(), 3);
        gzip.flush().unwrap();
        gzip.finish().unwrap();
    }

    #[test]
    fn test_alignment_is_written_verbatim_to_sam_and_bam() {
        let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
                     r1\t83\tchr1\t10\t60\t2S4M\t=\t40\t-36\tAACCGG\tABCDEF\tMC:Z:6M\tRG:Z:x\n";
        let mut reader = sam::io::Reader::new(&text[..]);
        let header = reader.read_header().unwrap();
        let mut record = RecordBuf::default();
        reader.read_record_buf(&header, &mut record).unwrap();
        let alignment = Alignment::Record(record.clone());
        let read = OutputRead {
            name: b"ignored",
            bases: b"T",
            quals: None,
            tags: None,
            read_group: Some("ignored"),
            qc_fail: false,
            alignment: Some(&alignment),
        };
        let sam_bytes =
            encode_fragment(OutputFormat::Sam, &header, std::slice::from_ref(&read)).unwrap();
        assert_eq!(
            String::from_utf8(sam_bytes).unwrap(),
            "r1\t83\tchr1\t10\t60\t2S4M\t=\t40\t-36\tAACCGG\tABCDEF\tMC:Z:6M\tRG:Z:x\n"
        );

        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.bam");
        let mut w =
            OutputWriter::create_with_header(Some(&out), OutputFormat::Bam, 5, header.clone())
                .unwrap();
        w.write_fragment(&[read]).unwrap();
        w.finish().unwrap();
        let mut reader = bam::io::Reader::new(File::open(&out).unwrap());
        let read_header = reader.read_header().unwrap();
        assert_eq!(read_header.reference_sequences().len(), 1);
        let mut back = RecordBuf::default();
        reader.read_record_buf(&read_header, &mut back).unwrap();
        assert_eq!(back, record);
    }

    #[test]
    fn test_bam_stores_u_bases_as_t_while_sam_keeps_them() {
        let text = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n\
                     r2\t0\tchr1\t10\t60\t4M\t*\t0\t0\tUUgu\tIIII\n";
        let mut reader = sam::io::Reader::new(&text[..]);
        let header = reader.read_header().unwrap();
        let mut record = RecordBuf::default();
        reader.read_record_buf(&header, &mut record).unwrap();
        let alignment = Alignment::Record(record);
        let fragments = [
            vec![output_read(b"r1", b"ACGUacgu", Some(&[30; 8]), None)],
            vec![OutputRead {
                alignment: Some(&alignment),
                ..output_read(b"", b"", None, None)
            }],
        ];

        let dir = tempfile::tempdir().unwrap();
        let direct = dir.path().join("direct.bam");
        let mut writer =
            OutputWriter::create_with_header(Some(&direct), OutputFormat::Bam, 5, header.clone())
                .unwrap();
        for fragment in &fragments {
            writer.write_fragment(fragment).unwrap();
        }
        writer.finish().unwrap();
        let encoded = dir.path().join("encoded.bam");
        let mut writer =
            OutputWriter::create_with_header(Some(&encoded), OutputFormat::Bam, 5, header.clone())
                .unwrap();
        for fragment in &fragments {
            let bytes = encode_fragment(OutputFormat::Bam, &header, fragment).unwrap();
            writer.write_encoded(&bytes).unwrap();
        }
        writer.finish().unwrap();
        for path in [&direct, &encoded] {
            let mut reader = bam::io::Reader::new(File::open(path).unwrap());
            let header = reader.read_header().unwrap();
            let bases: Vec<Vec<u8>> = reader
                .record_bufs(&header)
                .map(|record| record.unwrap().sequence().as_ref().to_vec())
                .collect();
            assert_eq!(bases, [b"ACGTACGT".to_vec(), b"TTGT".to_vec()], "{path:?}");
        }

        let sam: Vec<u8> = fragments
            .iter()
            .flat_map(|fragment| encode_fragment(OutputFormat::Sam, &header, fragment).unwrap())
            .collect();
        let sam = String::from_utf8(sam).unwrap();
        assert!(sam.contains("\tACGUacgu\t"), "{sam}");
        assert!(sam.contains("\tUUgu\t"), "{sam}");
    }

    #[test]
    fn test_raw_bam_with_an_undeclared_reference_id_fails_like_the_encoder() {
        use noodles::sam::alignment::io::Write as _;
        let text = b"@SQ\tSN:chr1\tLN:100\n\
                     r1\t99\tchr1\t10\t60\t4M\t=\t40\t34\tACGT\tABCD\tRG:Z:x\n";
        let mut reader = sam::io::Reader::new(&text[..]);
        let header = reader.read_header().unwrap();
        let mut record = RecordBuf::default();
        reader.read_record_buf(&header, &mut record).unwrap();
        let mut writer = bam::io::Writer::from(Vec::new());
        writer.write_alignment_record(&header, &record).unwrap();
        let raw = writer.into_inner();
        let encode = |raw: Vec<u8>| {
            let alignment = Alignment::RawBam(raw);
            let read = OutputRead {
                name: b"",
                bases: b"",
                quals: None,
                tags: None,
                read_group: None,
                qc_fail: false,
                alignment: Some(&alignment),
            };
            encode_fragment(OutputFormat::Bam, &header, &[read]).map_err(|e| format!("{e:#}"))
        };
        assert_eq!(encode(raw.clone()).unwrap(), raw);
        for (at, id, want) in [
            (4, 7, "invalid reference sequence ID"),
            (24, 9, "invalid mate reference sequence ID"),
            (4, -2, "invalid reference sequence ID"),
        ] {
            let mut bad = raw.clone();
            bad[at..at + 4].copy_from_slice(&i32::to_le_bytes(id));
            let err = encode(bad).unwrap_err();
            assert!(err.contains(want), "{err}");
        }
    }
}
