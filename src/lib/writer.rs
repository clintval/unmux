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
//! reported, not silently swallowed. The per-sample fan-out over these single
//! sinks lives in [`crate::output`]. A gzipped-FASTX or BAM sink can be a
//! `pooled-writer` front-end ([`SinkWriter::Pooled`]) so the many-file
//! fan-out's BGZF compression runs on a shared pool off the consumer thread;
//! the pooled BAM writer is built via `From` so it carries no inner BGZF (the
//! pool is the sole stage).

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::data::field::value::Array;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Data, QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, bgzf, cram, fasta, fastq, sam};
use pooled_writer::PooledWriter;

use crate::input::{SniffedFormat, PHRED_OFFSET};

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
    /// Whether this format *can* be written through the shared BGZF compressor
    /// pool: a pure format property. Gzipped FASTX (BGZF is valid gzip; the
    /// FASTX writer adds no compression of its own) and BAM (built via `From`,
    /// so the bam writer adds no inner BGZF) can be, with the pool as the sole
    /// BGZF stage. Text SAM, CRAM (own codecs), and uncompressed FASTX are not
    /// BGZF. Whether a given destination is actually pooled is the engine's
    /// call (only file destinations are exchanged, and `--out` SAM/BAM/CRAM
    /// stay on the inline single-EOF path).
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
    Plain(BufWriter<Box<dyn Write>>),
    /// A gzip encoder over a buffered stream (FASTX.gz output); boxed as it is
    /// much larger than the plain variant.
    Gzip(Box<flate2::write::GzEncoder<BufWriter<Box<dyn Write>>>>),
    /// A `pooled_writer` front-end that ships blocks to a shared BGZF
    /// compressor pool: the many-file FASTX.gz / BAM fan-out compresses off the
    /// consumer thread. BGZF is valid gzip, so a pooled `.fq.gz`/`.fa.gz` stays
    /// readable by any gzip reader; a pooled BAM gets its BGZF from the pool
    /// (the bam writer adds none). The pool is owned and stopped by the engine
    /// after every writer is finalized.
    Pooled(PooledWriter),
}

impl Write for SinkWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            SinkWriter::Plain(writer) => writer.write(buf),
            SinkWriter::Gzip(writer) => writer.write(buf),
            SinkWriter::Pooled(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            SinkWriter::Plain(writer) => writer.flush(),
            SinkWriter::Gzip(writer) => writer.flush(),
            SinkWriter::Pooled(writer) => writer.flush(),
        }
    }
}

impl SinkWriter {
    /// Finalize the stream: write the gzip trailer (if any) and flush,
    /// propagating any error. For a pooled writer this flushes the final block
    /// and the BGZF EOF into the pool; the pool itself is stopped by the engine
    /// once every writer has been finalized (the mandatory order).
    fn finish(self) -> Result<()> {
        match self {
            SinkWriter::Plain(mut writer) => writer.flush().context("failed to flush output"),
            SinkWriter::Gzip(writer) => {
                let mut inner = (*writer)
                    .finish()
                    .context("failed to finish the gzip stream")?;
                inner.flush().context("failed to flush output")
            }
            // `close` flushes the final data block and the BGZF EOF,
            // propagating any error (the point of calling it over a bare drop).
            // It consumes the writer, so its `Drop` runs the flush once more on
            // the now-empty buffer: the compressor still deflates the empty
            // input into a real (non-EOF, ISIZE=0) block and then appends a
            // second EOF. So the tail is `[...data][EOF][empty block][EOF]`,
            // not two adjacent EOFs. This is harmless: all payload precedes the
            // first EOF, the trailing empty member decodes to nothing, and the
            // last 28 bytes are a canonical EOF, so noodles (and htslib) read
            // every record then stop. A truly single-EOF
            // stream would mean giving up `close`'s error propagation.
            SinkWriter::Pooled(writer) => writer
                .close()
                .context("failed to finish the pooled BGZF stream"),
        }
    }
}

/// One open output sink, holding the concrete noodles writer for the chosen
/// format.
enum Sink {
    /// FASTQ. The `record` buffer is reused across writes so the per-record
    /// SEQ/QUAL Vecs are cleared and refilled in place rather than reallocated
    /// for every output record.
    Fastq {
        /// The record writer.
        writer: fastq::io::Writer<SinkWriter>,
        /// The reused output record (its SEQ/QUAL buffers are cleared and
        /// refilled each write).
        record: fastq::Record,
    },
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
    /// BGZF BAM where the bam writer wraps its own BGZF encoder (stdout and
    /// other non-pooled sinks), plus the header.
    BamInline {
        /// The record writer (its inner `bgzf::io::Writer` does the
        /// compression).
        writer: bam::io::Writer<bgzf::io::Writer<SinkWriter>>,
        /// The header.
        header: sam::Header,
    },
    /// BGZF BAM whose BGZF is done by the `pooled-writer` pool: the bam writer
    /// carries NO inner bgzf layer (built via `From`, not `new`), so the pooled
    /// sink is the sole BGZF stage.
    BamPooled {
        /// The record writer, writing BAM binary straight into the pooled sink.
        writer: bam::io::Writer<SinkWriter>,
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
        let raw: Box<dyn Write> = match out {
            Some(path) => Box::new(create_output_file(path)?),
            None => Box::new(std::io::stdout()),
        };
        let buffered = BufWriter::new(raw);

        let sink = match format {
            OutputFormat::Fastq { gzip } => Sink::Fastq {
                writer: fastq::io::Writer::new(sink_writer(buffered, gzip, compression)),
                record: fastq::Record::default(),
            },
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

    /// Create a writer whose sink is a `pooled_writer` front-end: records are
    /// serialized here and BGZF-compressed on the shared pool, so the many-file
    /// fan-out does not compress on the consumer thread. Poolable formats are
    /// gzipped FASTX (the BGZF output is valid gzip) and BAM (built via `From`,
    /// so the bam writer adds no inner BGZF and the pool is the sole
    /// compressor); every other format stays on [`create_with_header`].
    /// `header` is used for BAM and ignored for FASTX.
    pub fn from_pooled(
        pooled: PooledWriter,
        format: OutputFormat,
        header: sam::Header,
    ) -> Result<Self> {
        let writer = SinkWriter::Pooled(pooled);
        let sink = match format {
            OutputFormat::Fastq { gzip: true } => Sink::Fastq {
                writer: fastq::io::Writer::new(writer),
                record: fastq::Record::default(),
            },
            OutputFormat::Fasta { gzip: true } => Sink::Fasta(fasta::io::Writer::new(writer)),
            OutputFormat::Bam => {
                let mut writer = bam::io::Writer::from(writer);
                writer
                    .write_header(&header)
                    .context("failed to write BAM header")?;
                Sink::BamPooled { writer, header }
            }
            other => bail!("internal error: format {other:?} is not poolable"),
        };
        Ok(Self { sink })
    }

    /// Write one logical record's output records. A multi-record fragment emits
    /// SAM/BAM paired-segment records (or interleaved FASTX); SAM/BAM allows at
    /// most two (a pair).
    pub fn write_fragment(&mut self, reads: &[OutputRead]) -> Result<()> {
        match &mut self.sink {
            Sink::Fastq { writer, record } => {
                for read in reads {
                    write_fastq(writer, record, read)?;
                }
            }
            Sink::Fasta(writer) => {
                for read in reads {
                    write_fasta(writer, read)?;
                }
            }
            Sink::Sam { writer, header } => write_alignments(writer, header, reads)?,
            Sink::BamInline { writer, header } => write_alignments(writer, header, reads)?,
            Sink::BamPooled { writer, header } => write_alignments(writer, header, reads)?,
            Sink::Cram { writer, header } => write_alignments(writer, header, reads)?,
        }
        Ok(())
    }

    /// Append a fragment's already-encoded record bytes (produced by
    /// [`encode_fragment`] for this sink's format, e.g. on a worker thread),
    /// after the header written at open. The bytes are injected into the sink's
    /// underlying stream so they pass through its compression stage (gzip /
    /// inline BGZF / the pooled BGZF) exactly as if the records had been
    /// encoded here, keeping output byte-identical to
    /// [`OutputWriter::write_fragment`]. CRAM is not supported (its records are
    /// not independent per-fragment byte slices); callers route CRAM through
    /// `write_fragment`.
    pub fn write_encoded(&mut self, bytes: &[u8]) -> Result<()> {
        let sink: &mut dyn Write = match &mut self.sink {
            Sink::Fastq { writer, .. } => writer.get_mut(),
            Sink::Fasta(writer) => writer.get_mut(),
            Sink::Sam { writer, .. } => writer.get_mut(),
            // get_mut yields the inner `bgzf::io::Writer`, so the raw record
            // bytes are BGZF-compressed here exactly as the bam writer would
            // have compressed them.
            Sink::BamInline { writer, .. } => writer.get_mut(),
            // get_mut yields the pooled `SinkWriter`; the pool is the sole BGZF
            // stage.
            Sink::BamPooled { writer, .. } => writer.get_mut(),
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
            Sink::Fastq { writer, .. } => writer.into_inner().finish(),
            Sink::Fasta(writer) => writer.into_inner().finish(),
            Sink::Sam { writer, .. } => writer.into_inner().finish(),
            Sink::BamInline { writer, .. } => {
                let sink = writer
                    .into_inner()
                    .finish()
                    .context("failed to finalize BGZF (BAM EOF block)")?;
                sink.finish()
            }
            // The pooled BAM writer has no inner BGZF layer; the pooled sink
            // emits the BGZF EOF on close (do NOT also call the bam writer's
            // bgzf finish - there is none).
            Sink::BamPooled { writer, .. } => writer.into_inner().finish(),
            Sink::Cram { mut writer, header } => {
                writer
                    .try_finish(&header)
                    .context("failed to finalize the CRAM EOF container")?;
                writer.into_inner().finish()
            }
        }
    }
}

/// A BGZF writer over `inner` at the given `--compression` level (0-9). Used by
/// the inline BAM sink so it honors `--compression` (the spec sets the BAM BGZF
/// level from it), matching the pooled BAM bins (which pass the level to the
/// compressor pool). `--compression` is clap-validated to 0-9, so the
/// conversion only errors on a programmer mistake.
fn bgzf_writer<W: Write>(inner: W, compression: u8) -> Result<bgzf::io::Writer<W>> {
    let level = bgzf::io::writer::CompressionLevel::try_from(compression)
        .map_err(|e| anyhow!("invalid --compression {compression} for BGZF: {e}"))?;
    Ok(bgzf::io::writer::Builder::default()
        .set_compression_level(level)
        .build_from_writer(inner))
}

/// Create an output file, making any missing parent directories first. Shared
/// by the inline writers and the pooled-writer exchange so both create files
/// (and dirs) identically.
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
fn sink_writer(buffered: BufWriter<Box<dyn Write>>, gzip: bool, compression: u8) -> SinkWriter {
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

/// Write one FASTQ record (qualities re-encoded to ASCII Phred+33), reusing
/// `record`'s SEQ/QUAL buffers in place to avoid allocating a fresh Vec per
/// output record.
fn write_fastq<W: Write>(
    writer: &mut fastq::io::Writer<W>,
    record: &mut fastq::Record,
    read: &OutputRead,
) -> Result<()> {
    let quals = read
        .quals
        .context("cannot write FASTQ output for a record with no qualities (a FASTA input)")?;
    let seq = record.sequence_mut();
    seq.clear();
    seq.extend_from_slice(read.bases);
    let qual = record.quality_scores_mut();
    qual.clear();
    qual.extend(quals.iter().map(|&q| q.saturating_add(PHRED_OFFSET)));
    // Carry input tags through in the read-name comment (samtools `fastq -T`
    // convention). The definition is small (name + optional comment), so it is
    // rebuilt rather than mutated in place.
    let comment = read.tags.map(render_tags).unwrap_or_default();
    *record.definition_mut() = fastq::record::Definition::new(read.name.to_vec(), comment);
    writer
        .write_record(record)
        .context("failed to write FASTQ record")
}

/// Write one FASTA record (no qualities), carrying input tags in the
/// description comment.
fn write_fasta<W: Write>(writer: &mut fasta::io::Writer<W>, read: &OutputRead) -> Result<()> {
    let comment = read
        .tags
        .map(render_tags)
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
fn render_tags(data: &Data) -> Vec<u8> {
    let mut fields: Vec<String> = Vec::new();
    for (tag, value) in data.iter() {
        let bytes: [u8; 2] = tag.into();
        let key = String::from_utf8_lossy(&bytes);
        let field = match value {
            Value::Character(c) => format!("{key}:A:{}", *c as char),
            Value::Int8(n) => format!("{key}:i:{n}"),
            Value::UInt8(n) => format!("{key}:i:{n}"),
            Value::Int16(n) => format!("{key}:i:{n}"),
            Value::UInt16(n) => format!("{key}:i:{n}"),
            Value::Int32(n) => format!("{key}:i:{n}"),
            Value::UInt32(n) => format!("{key}:i:{n}"),
            Value::Float(x) => format!("{key}:f:{x}"),
            Value::String(s) => format!("{key}:Z:{s}"),
            Value::Hex(s) => format!("{key}:H:{s}"),
            Value::Array(array) => format!("{key}:B:{}", render_array(array)),
        };
        fields.push(field);
    }
    fields.join("\t").into_bytes()
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

/// Write a fragment's records as SAM/BAM alignment records, setting
/// paired-segment flags for a multi-record fragment. SAM/BAM holds at most two
/// records per template (a pair).
fn write_alignments<W: AlignmentWrite>(
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
        let record = build_record(read, segment_flags(index, reads.len()));
        writer
            .write_alignment_record(header, &record)
            .context("failed to write SAM/BAM/CRAM record")?;
    }
    Ok(())
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

/// Encode one fragment's records into the exact bytes the matching sink would
/// write for `format`, via a throwaway in-memory writer. The result is
/// byte-identical to what [`OutputWriter::write_fragment`] emits for the same
/// records (header excluded), so it can be produced on a worker thread and
/// appended verbatim by [`OutputWriter::write_encoded`] in input order
/// (encode-on-workers). `header` is used only for SAM/BAM record encoding and
/// ignored for FASTX. CRAM is NOT supported here (its records are
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
            let mut writer = fastq::io::Writer::new(Vec::new());
            let mut record = fastq::Record::default();
            for read in reads {
                write_fastq(&mut writer, &mut record, read)?;
            }
            writer.into_inner()
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
        OutputFormat::Bam => {
            let mut writer = bam::io::Writer::from(Vec::new());
            write_alignments(&mut writer, header, reads)?;
            writer.into_inner()
        }
        OutputFormat::Cram => {
            bail!("internal error: CRAM is not eligible for encode-on-workers")
        }
    })
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
            }],
            vec![
                OutputRead {
                    name: b"r2",
                    bases: b"AAACCC",
                    quals: Some(&q),
                    tags: None,
                    read_group: Some("dna01"),
                },
                OutputRead {
                    name: b"r2",
                    bases: b"GGGTTT",
                    quals: Some(&q),
                    tags: None,
                    read_group: Some("dna01"),
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
}
