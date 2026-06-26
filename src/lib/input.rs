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
//! interleaved (two mates per fragment). An input path of `-` reads stdin (at
//! most one). CRAM input is read reference-free (unmapped); a
//! reference-compressed CRAM (one carrying `@SQ`) is rejected, since no
//! reference is configured.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, bgzf, cram, fasta, fastq, sam};
use smallvec::SmallVec;

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
    /// BGZF-compressed BAM (the bam reader decodes BGZF internally).
    Bam {
        /// The record reader.
        reader: bam::io::Reader<bgzf::io::Reader<Box<dyn BufRead + Send>>>,
        /// The header, needed to decode each record.
        header: Box<sam::Header>,
        /// The reused per-record parse buffer (`read_record_buf` clears and
        /// refills it each call).
        record: RecordBuf,
    },
    /// Reference-free / unmapped CRAM. The records iterator borrows the reader
    /// and header for its whole lifetime (so it cannot be stored), so we decode
    /// one container at a time into `buffered` and serve from it, mirroring
    /// `cram::io::reader::Slice::records`'s documented public loop.
    Cram {
        /// The record reader.
        reader: cram::io::Reader<Box<dyn BufRead + Send>>,
        /// The header, needed to decode each record.
        header: Box<sam::Header>,
        /// The reusable container buffer (one container of slices is decoded at
        /// a time).
        container: cram::io::reader::Container,
        /// The empty reference repository (reference-free CRAM never consults
        /// it).
        repository: fasta::Repository,
        /// The decoded records of the current container, served one per
        /// `read_next`.
        buffered: std::vec::IntoIter<RecordBuf>,
    },
    /// An empty input (zero records): no underlying reader. `read_next` always
    /// returns `None`. Used when the input has no content and no format could
    /// be inferred from its extension.
    Empty,
}

/// Reject a reference-compressed CRAM (one whose header carries `@SQ` reference
/// sequences). Such a CRAM stores SEQ as deltas against an external reference
/// unmux does not have, and noodles panics (rather than erroring) when decoding
/// it with the empty repository, so it is rejected up front. Reference-free /
/// unmapped CRAM (no `@SQ` - what unmux writes and what raw demux inputs are)
/// reads with the default empty repository. Conservative: an `@SQ`-bearing but
/// all-unmapped CRAM is also rejected, since unmux has no reference repository
/// to decode against.
///
/// Residual (corrupt input only): a CRAM with an EMPTY `@SQ` header but a
/// mapped slice referencing a reference id - invalid per the CRAM spec, so no
/// conforming writer emits it - would still abort via a noodles panic during
/// decode. The slice-level reference context is not exposed by noodles 0.92's
/// public API, so it cannot be pre-checked here; well-formed CRAM
/// (reference-free, or reference-compressed which this rejects) never reaches
/// that path.
fn ensure_reference_free_cram(header: &sam::Header, path: &Path) -> Result<()> {
    let references = header.reference_sequences().len();
    if references > 0 {
        bail!(
            "CRAM input is reference-compressed (carries {references} @SQ reference sequence(s)); \
             unmux demux expects unmapped/reference-free CRAM: {}",
            path.display()
        );
    }
    Ok(())
}

/// Decode the next CRAM container into owned records, or `None` at end of
/// input. Mirrors the public loop documented on
/// `cram::io::reader::Slice::records`: read a container, then for each slice
/// decode its blocks and turn each raw record into an owned [`RecordBuf`]
/// (shared with the BAM/SAM path). Buffers one container at a time (bounded by
/// container size, not the whole file). Pinned to noodles-cram 0.92: the
/// `Container`/`Slice` API and `Slice::records` signature could change on a
/// bump.
fn read_cram_container(
    reader: &mut cram::io::Reader<Box<dyn BufRead + Send>>,
    container: &mut cram::io::reader::Container,
    header: &sam::Header,
    repository: &fasta::Repository,
) -> Result<Option<Vec<RecordBuf>>> {
    if reader
        .read_container(container)
        .context("failed to read CRAM container")?
        == 0
    {
        return Ok(None);
    }
    let compression_header = container
        .compression_header()
        .context("failed to read CRAM compression header")?;
    let mut records = Vec::new();
    for slice in container.slices() {
        let slice = slice.context("failed to read CRAM slice")?;
        let (core, external) = slice
            .decode_blocks()
            .context("failed to decode CRAM slice blocks")?;
        for record in slice
            .records(
                repository.clone(),
                header,
                &compression_header,
                &core,
                &external,
            )
            .context("failed to read CRAM slice records")?
        {
            records.push(
                RecordBuf::try_from_alignment_record(header, &record)
                    .context("failed to convert CRAM record")?,
            );
        }
    }
    Ok(Some(records))
}

impl FormatReader {
    /// Read the next record from this stream, or `None` at end of input.
    fn read_next(&mut self) -> Result<Option<InputRecord>> {
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
                Ok(Some(InputRecord {
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
                }))
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
                Ok(Some(InputRecord {
                    name,
                    bases,
                    quals: None,
                    tags: parse_comment_tags(comment.as_bytes()),
                }))
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
                // Secondary/supplementary alignments are not distinct reads;
                // skip to the next primary record so the per-file segment count
                // stays one-per-template.
                let flags = record.flags();
                if flags.is_secondary() || flags.is_supplementary() {
                    continue;
                }
                return Ok(Some(alignment_to_input(record)));
            },
            FormatReader::Bam {
                reader,
                header,
                record,
            } => loop {
                if reader
                    .read_record_buf(header, record)
                    .context("failed to read BAM record")?
                    == 0
                {
                    return Ok(None);
                }
                let flags = record.flags();
                if flags.is_secondary() || flags.is_supplementary() {
                    continue;
                }
                return Ok(Some(alignment_to_input(record)));
            },
            FormatReader::Cram {
                reader,
                header,
                container,
                repository,
                buffered,
            } => loop {
                if let Some(record) = buffered.next() {
                    let flags = record.flags();
                    if flags.is_secondary() || flags.is_supplementary() {
                        continue;
                    }
                    return Ok(Some(alignment_to_input(&record)));
                }
                match read_cram_container(reader, container, header, repository)? {
                    None => return Ok(None),
                    Some(records) => *buffered = records.into_iter(),
                }
            },
            FormatReader::Empty => Ok(None),
        }
    }
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
const ALIGNMENT_TAGS: &[[u8; 2]] = &[
    *b"AS", *b"NM", *b"MD", *b"SA", *b"MC", *b"MQ", *b"ms", *b"XS", *b"XA", *b"XB", *b"XR", *b"pa",
    *b"UQ", *b"CC", *b"CP", *b"CG", *b"H0", *b"H1", *b"H2", *b"HI", *b"IH", *b"NH", *b"TS",
];

/// Convert an alignment record into an [`InputRecord`], carrying its
/// non-alignment tags through. Alignment-derived tags ([`ALIGNMENT_TAGS`]) are
/// stripped, since the output is unmapped. An empty quality string (SAM `*`)
/// becomes no qualities.
///
/// A reverse-strand record (flag `0x10`) stores `SEQ` as the reverse complement
/// of the sequenced read and `QUAL` in reversed order; this restores the
/// original read orientation so matching and extraction see the sequenced
/// bases.
fn alignment_to_input(record: &RecordBuf) -> InputRecord {
    let raw_quals = record.quality_scores().as_ref();
    let mut bases = record.sequence().as_ref().to_vec();
    let mut quals = (!raw_quals.is_empty()).then(|| raw_quals.to_vec());
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
    }
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

/// Open one input by sniffing its format and dispatching to the matching
/// noodles reader, returning the sniffed format alongside (so the caller can
/// validate the input family and mirror the format to stdout output).
fn open_format_reader(path: &Path) -> Result<(SniffedFormat, FormatReader)> {
    let (format, gzipped, reader) =
        sniff_input(path).with_context(|| format!("failed to read: {}", path.display()))?;
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
            // BAM is BGZF, and the bam reader BGZF-decodes internally, so it
            // takes the raw stream.
            let mut reader = bam::io::Reader::new(reader);
            let header = reader.read_header().context("failed to read BAM header")?;
            FormatReader::Bam {
                reader,
                header: Box::new(header),
                record: RecordBuf::default(),
            }
        }
        SniffedFormat::Cram => {
            // CRAM takes a plain Read (no BGZF/gzip wrapper; gzipped CRAM is
            // already rejected).
            let mut reader = cram::io::Reader::new(reader);
            let header = reader.read_header().context("failed to read CRAM header")?;
            ensure_reference_free_cram(&header, path)?;
            FormatReader::Cram {
                reader,
                header: Box::new(header),
                container: cram::io::reader::Container::default(),
                repository: fasta::Repository::default(),
                buffered: Vec::new().into_iter(),
            }
        }
        SniffedFormat::Empty => FormatReader::Empty,
        SniffedFormat::Unknown => {
            bail!("could not identify the input format of {}", path.display())
        }
    };
    Ok((format, format_reader))
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

/// Reads logical records across all inputs. With several `--in` files, one
/// record is pulled from each in lockstep; a single FASTX input is
/// auto-detected as interleaved (two mates per fragment) unless `--per-record`
/// opts out.
pub struct FragmentReader {
    readers: Vec<FormatReader>,
    /// The sniffed format of each input, in `--in` order (blend-compatible by
    /// construction).
    formats: Vec<SniffedFormat>,
    /// A single FASTX input detected as interleaved: pull two mates per
    /// fragment from `readers[0]`.
    interleaved: bool,
    /// Records buffered from `readers[0]` during interleaving detection, served
    /// before fresh reads.
    buffered: Vec<InputRecord>,
}

impl FragmentReader {
    /// Open every input by sniffing and dispatching, preserving 0-based `--in`
    /// order. With a single FASTX input and `per_record` false, the stream is
    /// probed for interleaving: a leading `/1` mate, or two adjacent records
    /// that share a mate-pairing key, marks it interleaved.
    pub fn open(paths: &[PathBuf], per_record: bool) -> Result<Self> {
        let opened = paths
            .iter()
            .enumerate()
            .map(|(i, path)| {
                open_format_reader(path)
                    .with_context(|| format!("failed to open input file {i} ({})", path.display()))
            })
            .collect::<Result<Vec<(SniffedFormat, FormatReader)>>>()?;

        let formats: Vec<SniffedFormat> = opened.iter().map(|(format, _)| *format).collect();
        ensure_compatible_quality(&formats)?;
        let mut readers: Vec<FormatReader> = opened.into_iter().map(|(_, reader)| reader).collect();

        let mut buffered = Vec::new();
        let mut interleaved = false;
        // Interleaving is only meaningful for a single FASTX input; SAM/BAM
        // segment pairing and multi-file inputs are handled by their own paths.
        if !per_record
            && readers.len() == 1
            && matches!(
                readers[0],
                FormatReader::Fastq { .. } | FormatReader::Fasta { .. }
            )
        {
            if let Some(first) = readers[0].read_next()? {
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
                } else if let Some(second) = readers[0].read_next()? {
                    interleaved = mate_key(&second.name) == mate_key(&buffered[0].name);
                    buffered.push(second);
                }
            }
        }

        Ok(Self {
            readers,
            formats,
            interleaved,
            buffered,
        })
    }

    /// The sniffed format of each input, in `--in` order. The writer uses this
    /// to mirror the input format to stdout, defaulting a multi-format blend to
    /// uncompressed BAM.
    pub fn formats(&self) -> &[SniffedFormat] {
        &self.formats
    }

    /// The number of records every fragment carries: two for a single
    /// auto-detected interleaved input, otherwise one per input file. This
    /// fixes the pass-through body count (one body per record), so output
    /// destinations can be enumerated before streaming.
    pub fn fragment_width(&self) -> usize {
        if self.interleaved {
            2
        } else {
            self.readers.len()
        }
    }

    /// Read the next logical record. A single interleaved input yields a
    /// two-mate fragment; a single plain input yields a one-record fragment;
    /// multiple inputs are zipped one record each. Returns `None` once every
    /// input is exhausted together. A file that ends before the others, or
    /// input segments whose read names disagree, is a fail-fast error (both
    /// would silently mispair records).
    pub fn next_fragment(&mut self) -> Result<Option<Fragment>> {
        if self.readers.len() == 1 {
            let first = match self.next_from_first()? {
                Some(record) => record,
                None => return Ok(None),
            };
            if !self.interleaved {
                return Ok(Some(Fragment {
                    records: vec![first],
                }));
            }
            let second = self.next_from_first()?.ok_or_else(|| {
                anyhow!("interleaved input ended on an unpaired mate (odd number of records)")
            })?;
            check_names_concord(&[&first, &second])?;
            return Ok(Some(Fragment {
                records: vec![first, second],
            }));
        }

        let mut maybe: SmallVec<[Option<InputRecord>; 4]> =
            SmallVec::with_capacity(self.readers.len());
        for reader in &mut self.readers {
            maybe.push(reader.read_next()?);
        }
        let present = maybe.iter().filter(|record| record.is_some()).count();
        if present == 0 {
            return Ok(None);
        }
        if present != maybe.len() {
            bail!("input files have unequal record counts (one input ended before the others)");
        }
        let records: Vec<InputRecord> = maybe.into_iter().map(Option::unwrap).collect();
        check_names_concord(&records.iter().collect::<SmallVec<[&InputRecord; 4]>>())?;
        Ok(Some(Fragment { records }))
    }

    /// The next record from the first input, draining any records buffered
    /// during interleaving detection before reading fresh ones.
    fn next_from_first(&mut self) -> Result<Option<InputRecord>> {
        if !self.buffered.is_empty() {
            return Ok(Some(self.buffered.remove(0)));
        }
        self.readers[0].read_next()
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
fn check_names_concord(records: &[&InputRecord]) -> Result<()> {
    let Some((first, rest)) = records.split_first() else {
        return Ok(());
    };
    let key = mate_key(&first.name);
    for other in rest {
        if mate_key(&other.name) != key {
            bail!(
                "input segments have mismatched read names (`{}` vs `{}`); the inputs are not in the same record order",
                String::from_utf8_lossy(&first.name),
                String::from_utf8_lossy(&other.name)
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
    fn test_cram_reference_compressed_is_rejected() {
        use noodles::sam::header::record::value::{map, Map};
        use noodles::sam::Header;
        use std::num::NonZeroUsize;

        // A header with @SQ reference sequences is (potentially)
        // reference-compressed: rejected up front (unmux has no reference, and
        // noodles would panic decoding it). A reference-free header (no @SQ) is
        // accepted.
        let with_sq = Header::builder()
            .set_header(Map::<map::Header>::default())
            .add_reference_sequence(
                b"chr1".to_vec(),
                Map::<map::ReferenceSequence>::new(NonZeroUsize::new(100).unwrap()),
            )
            .build();
        let err = ensure_reference_free_cram(&with_sq, Path::new("in.cram"))
            .expect_err("reference-compressed CRAM should be rejected");
        assert!(
            format!("{err:#}").contains("reference-compressed"),
            "{err:#}"
        );
        assert!(ensure_reference_free_cram(&Header::default(), Path::new("in.cram")).is_ok());
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
}
