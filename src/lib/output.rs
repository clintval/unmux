//! The demux output layer: resolve `--out`/`--unassigned`/`--remove` path
//! patterns to concrete destinations, assemble the `--tag` SAM fields and `@RG`
//! header for a fan-out, and fan records out to one [`OutputWriter`] per
//! destination.
//!
//! Routing decides *which* target a record belongs to (see [`crate::fanout`]);
//! this module turns that decision into files. A pattern carries placeholders
//! (`%pool`/`%sample`/`%sub_sample`/`%ordinal` for `--out`, `%pool`/`%source`
//! for the raw bins); [`resolve_pattern`] expands them against a
//! [`PathContext`]. The engine enumerates every directed destination and calls
//! [`MultiWriter::create_all`] to create each file up front (a directed file
//! always exists, even empty); [`MultiWriter::write`] then reuses the open
//! writer (with a lazy fallback for any path not pre-created). Gzipped-FASTX
//! and BAM file destinations are compressed by a shared `pooled-writer` BGZF
//! pool, installed via [`MultiWriter::set_pooled`]; text SAM, CRAM,
//! uncompressed FASTX, and stdout use an inline sink.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{bail, Result};
use noodles::sam;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::header::record::value::map::{
    self, program::tag as program_tag, read_group::tag as rg_tag, Program, ReadGroup,
};
use noodles::sam::header::record::value::Map;
use pooled_writer::PooledWriter;

use crate::extract::Extracted;
use crate::fanout::Target;
use crate::grammar::{
    default_qual_tag, OutputPattern, PatternSegment, Placeholder, QualTag, TagBinding,
};
use crate::input::{SniffedFormat, PHRED_OFFSET};
use crate::writer::{output_format, OutputRead, OutputWriter};

/// A non-standard read-group field tag (e.g. `SM`/`LB`/`CN`), as the
/// `Map<ReadGroup>` builder takes.
type ReadGroupTag = map::tag::Other<<ReadGroup as map::Inner>::StandardTag>;

/// The runtime values an output pattern's placeholders expand to. `%sub_sample`
/// with no sub_sample, and any other absent field, expand to the empty string
/// (per the grammar's defaults).
#[derive(Debug, Clone, Copy, Default)]
pub struct PathContext<'a> {
    /// The pool id (`%pool`).
    pub pool: &'a str,
    /// The sample id (`%sample`), when fanning out by sample.
    pub sample: Option<&'a str>,
    /// The sub_sample id (`%sub_sample`); `None` expands to empty.
    pub sub_sample: Option<&'a str>,
    /// The 1-based `--template` ordinal (`%ordinal`), when fanning out per
    /// output record.
    pub ordinal: Option<usize>,
    /// The 0-based input file index (`%source`), for the raw
    /// `--unassigned`/`--remove` bins.
    pub source: Option<usize>,
}

/// Expand a pattern's segments to a path string, pulling placeholder values
/// from `ctx`.
fn render_pattern(pattern: &OutputPattern, ctx: &PathContext) -> String {
    let mut path = String::new();
    for segment in &pattern.segments {
        match segment {
            PatternSegment::Literal(text) => path.push_str(text),
            PatternSegment::Placeholder(placeholder) => {
                path.push_str(&placeholder_value(*placeholder, ctx))
            }
        }
    }
    path
}

/// Resolve an output pattern to a destination: a concrete path, or stdout
/// (`None`) when the pattern is the literal `-` or `/dev/stdout`. Placeholders
/// pull from `ctx`.
pub fn resolve_pattern(pattern: &OutputPattern, ctx: &PathContext) -> Option<PathBuf> {
    let path = render_pattern(pattern, ctx);
    if path == "-" || path == "/dev/stdout" {
        None
    } else {
        Some(PathBuf::from(path))
    }
}

/// Resolve a pool-level metrics path pattern (grammar guarantees only `%pool`
/// is present) to a concrete path. Unlike [`resolve_pattern`], the stdout
/// literals are not special-cased: a metrics TSV is always a real file (the
/// same metrics are also streamed to stderr), so a literal `-` simply names a
/// file `-`.
pub fn resolve_metrics_path(pattern: &OutputPattern, pool: &str) -> PathBuf {
    let ctx = PathContext {
        pool,
        ..Default::default()
    };
    PathBuf::from(render_pattern(pattern, &ctx))
}

/// The string a single placeholder expands to in the current context.
fn placeholder_value(placeholder: Placeholder, ctx: &PathContext) -> String {
    match placeholder {
        Placeholder::Pool => ctx.pool.to_string(),
        Placeholder::Sample => ctx.sample.unwrap_or_default().to_string(),
        Placeholder::SubSample => ctx.sub_sample.unwrap_or_default().to_string(),
        Placeholder::Ordinal => ctx.ordinal.map(|o| o.to_string()).unwrap_or_default(),
        Placeholder::Source => ctx.source.map(|s| s.to_string()).unwrap_or_default(),
    }
}

/// Assemble the `--tag` SAM data fields for one record from its extracted
/// streams. Each binding's sequence tag is the present streams' bases joined by
/// `sep` (a stream absent because its anchoring group did not match is skipped;
/// a tag with no present stream is omitted entirely, so an unmatched group
/// leaves a null tag rather than an empty one). The paired quality tag (from
/// the default seq->qual map, an explicit `qual=`, or off for `qual=none`) is
/// the present streams' qualities re-encoded to ASCII Phred+33 and joined by
/// `qual-sep`; it is emitted only when every present stream carries qualities.
pub fn build_tag_data(bindings: &[TagBinding], streams: &HashMap<&str, Extracted>) -> Data {
    let mut data = Data::default();
    for binding in bindings {
        // Each present stream is paired with its `StreamRef` so a `~`
        // (reverse-complement at this use site) can flip that stream's bases
        // and qualities.
        let present: Vec<_> = binding
            .streams
            .iter()
            .filter_map(|s| streams.get(s.name.as_str()).map(|e| (s, e)))
            .collect();
        if present.is_empty() {
            continue;
        }
        let sep = binding.sep.as_deref().unwrap_or_default().as_bytes();
        // A tag emits the corrected form by default and the observed bases
        // under `raw=true`; a stream with no corrected form (not an
        // error-corrected `@grp` match) yields its observed bases for either. A
        // `~` reverse-complements the chosen form here.
        let bases: Vec<Vec<u8>> = present
            .iter()
            .map(|(s, stream)| {
                let mut bases = stream.tag_bases(binding.raw).to_vec();
                if s.revcomp {
                    crate::iupac::reverse_complement(&mut bases);
                }
                bases
            })
            .collect();
        data.insert(
            two_byte_tag(&binding.tag),
            string_value(&join_with(&bases, sep)),
        );

        if let Some(qual_tag) = qual_tag_name(binding) {
            if present
                .iter()
                .all(|(_, stream)| stream.tag_quals(binding.raw).is_some())
            {
                let qual_sep = binding.qual_sep.as_deref().unwrap_or_default().as_bytes();
                let quals: Vec<Vec<u8>> = present
                    .iter()
                    .map(|(s, stream)| {
                        let mut quals: Vec<u8> = stream
                            .tag_quals(binding.raw)
                            .unwrap()
                            .iter()
                            .map(|q| q.saturating_add(PHRED_OFFSET))
                            .collect();
                        // A `~`-flipped stream's qualities reverse to co-orient
                        // with its rc bases.
                        if s.revcomp {
                            quals.reverse();
                        }
                        quals
                    })
                    .collect();
                data.insert(
                    two_byte_tag(qual_tag),
                    string_value(&join_with(&quals, qual_sep)),
                );
            }
        }
    }
    data
}

/// The quality tag name a binding emits: the explicit `qual=QTAG`, the
/// seq->qual default, or `None` for `qual=none` (and for a sequence tag with no
/// default mapping).
pub(crate) fn qual_tag_name(binding: &TagBinding) -> Option<&str> {
    match &binding.qual {
        QualTag::None => None,
        QualTag::Named(name) => Some(name.as_str()),
        QualTag::Default => default_qual_tag(&binding.tag),
    }
}

/// A SAM data field tag from a validated two-character tag name.
fn two_byte_tag(name: &str) -> Tag {
    let bytes = name.as_bytes();
    Tag::new(bytes[0], bytes[1])
}

/// A `Z` string field value from raw bytes (ASCII bases / Phred+33 qualities).
fn string_value(bytes: &[u8]) -> Value {
    Value::from(String::from_utf8_lossy(bytes).into_owned())
}

/// Insert (or overwrite) a `Z` string SAM tag, e.g. the `--qc-tag`
/// demux-provenance slug. `tag` is a validated two-character name.
pub fn insert_string_tag(data: &mut Data, tag: &str, value: &str) {
    data.insert(two_byte_tag(tag), Value::from(value.to_string()));
}

/// Concatenate `segments` with `sep` between consecutive items (no leading or
/// trailing separator).
fn join_with(segments: &[Vec<u8>], sep: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    for (index, segment) in segments.iter().enumerate() {
        if index > 0 {
            out.extend_from_slice(sep);
        }
        out.extend_from_slice(segment);
    }
    out
}

/// Build the `@PG` program provenance record (`PN:unmux VN:<version>
/// CL:<command line>`); the `ID` is the `add_program` key and `CL` is omitted
/// when no command line was captured.
fn program_record(command_line: Option<&str>) -> Map<Program> {
    let mut builder = Map::<Program>::builder()
        .insert(program_tag::NAME, "unmux")
        .insert(program_tag::VERSION, env!("CARGO_PKG_VERSION"));
    if let Some(command_line) = command_line {
        builder = builder.insert(program_tag::COMMAND_LINE, command_line);
    }
    builder
        .build()
        .expect("a program map has no required fields")
}

/// A SAM header builder seeded with `@HD` and the unmux `@PG` provenance record
/// (mandatory on every SAM/BAM output).
fn provenance_header_builder(command_line: Option<&str>) -> sam::header::Builder {
    sam::Header::builder()
        .set_header(Map::<map::Header>::default())
        .add_program("unmux", program_record(command_line))
}

/// A minimal SAM/BAM header carrying just `@HD` + the `@PG` provenance, for the
/// raw `--unassigned`/`--remove` bins and non-demux output (no `@RG`).
pub fn provenance_header(command_line: Option<&str>) -> sam::Header {
    provenance_header_builder(command_line).build()
}

/// Build a SAM header carrying the `@PG` provenance plus one `@RG` line per
/// fan-out target (`ID` = the target label, `SM` = the sample, `LB` = the
/// sub_sample) and the shared `--rg-tag` fields on every read group. The
/// targets are those whose records land in this output file.
pub fn read_group_header(
    targets: &[&Target],
    rg_tags: &[(String, String)],
    command_line: Option<&str>,
) -> Result<sam::Header> {
    let mut builder = provenance_header_builder(command_line);
    for target in targets {
        let mut read_group =
            Map::<ReadGroup>::builder().insert(rg_tag::SAMPLE, target.sample.as_str());
        if let Some(sub_sample) = &target.sub_sample {
            read_group = read_group.insert(rg_tag::LIBRARY, sub_sample.as_str());
        }
        for (key, value) in rg_tags {
            read_group = read_group.insert(rg_shared_tag(key)?, value.as_str());
        }
        let map = read_group
            .build()
            .expect("a read group map has no required fields");
        builder = builder.add_read_group(target.label(), map);
    }
    Ok(builder.build())
}

/// Build the fallback header for destinations with no per-sample read group:
/// the pass-through `--out` and the raw `--unassigned`/`--remove` bins. It
/// carries the `@PG` provenance plus a single default `@RG` whose
/// `ID`/`SM`/`LB` are the pool id, plus the shared `--rg-tag` fields. This
/// makes every SAM/BAM/CRAM record unmux emits belong to a read group even when
/// no sample was declared (FASTX ignores it). The pool id is always non-empty
/// (it defaults to the input stem, else `pool`).
pub fn default_read_group_header(
    pool: &str,
    rg_tags: &[(String, String)],
    command_line: Option<&str>,
) -> Result<sam::Header> {
    let mut builder = provenance_header_builder(command_line);
    let mut read_group = Map::<ReadGroup>::builder()
        .insert(rg_tag::SAMPLE, pool)
        .insert(rg_tag::LIBRARY, pool);
    for (key, value) in rg_tags {
        read_group = read_group.insert(rg_shared_tag(key)?, value.as_str());
    }
    let map = read_group
        .build()
        .expect("a read group map has no required fields");
    builder = builder.add_read_group(pool.to_string(), map);
    Ok(builder.build())
}

/// Resolve a shared `--rg-tag` key (`CN`/`PL`/`PU`/...) to a non-standard
/// read-group tag. `ID` is the read-group key, not a field, so setting it is a
/// fail-fast error (`SM`/`LB` are already rejected by the grammar).
fn rg_shared_tag(key: &str) -> Result<ReadGroupTag> {
    let bytes = key.as_bytes();
    if bytes.len() != 2 {
        bail!("--rg-tag key `{key}` must be a two-character SAM tag");
    }
    match map::tag::Tag::from([bytes[0], bytes[1]]) {
        map::tag::Tag::Other(other) => Ok(other),
        map::tag::Tag::Standard(_) => {
            bail!("--rg-tag may not set `{key}`; ID is the read-group identifier")
        }
    }
}

/// The destination key for a writer: a concrete path, or the single shared
/// stdout sink.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum Dest {
    /// stdout (the `-` / `/dev/stdout` destination).
    Stdout,
    /// A file path.
    Path(PathBuf),
}

impl Dest {
    /// Map a resolved destination (path or stdout) to a key.
    fn from_path(path: Option<&Path>) -> Self {
        match path {
            Some(path) => Dest::Path(path.to_path_buf()),
            None => Dest::Stdout,
        }
    }
}

/// Fans records out to one [`OutputWriter`] per destination, opening each
/// writer lazily on first write. A SAM/BAM file is created with its precomputed
/// per-sample `@RG` header (looked up by path); a path with no precomputed
/// header (the pass-through `--out` and the raw `--unassigned`/`--remove` bins)
/// gets `default_header`, which carries the `@PG` provenance plus the default
/// pool `@RG`. stdout's format mirrors the inputs.
pub struct MultiWriter {
    writers: HashMap<Dest, OutputWriter>,
    headers: HashMap<Dest, sam::Header>,
    /// Pre-exchanged BGZF pool writers, keyed by destination: a pooled file's
    /// writer is installed from here on open instead of creating its own file.
    /// Populated by [`set_pooled`]; the pool itself is owned and stopped by the
    /// engine after every writer is finalized.
    pooled: HashMap<Dest, PooledWriter>,
    input_formats: Vec<SniffedFormat>,
    compression: u8,
    /// Header for any destination absent from `headers` (pass-through `--out`,
    /// the raw `--unassigned`/`--remove` bins): `@PG` provenance plus the
    /// default pool `@RG`.
    default_header: sam::Header,
}

impl MultiWriter {
    /// Create a fan-out writer. `headers` maps a resolved destination (a path,
    /// or `None` for stdout) to the SAM/BAM header it should be created with
    /// (carrying its per-sample `@RG`); destinations absent from the map get
    /// `default_header` (the `@PG` provenance plus the default pool `@RG`).
    /// `input_formats` picks the stdout mirror format; `compression` is the
    /// gzip/BGZF level.
    pub fn new(
        headers: HashMap<Option<PathBuf>, sam::Header>,
        default_header: sam::Header,
        input_formats: Vec<SniffedFormat>,
        compression: u8,
    ) -> Self {
        let headers = headers
            .into_iter()
            .map(|(path, header)| (Dest::from_path(path.as_deref()), header))
            .collect();
        Self {
            writers: HashMap::new(),
            headers,
            pooled: HashMap::new(),
            input_formats,
            compression,
            default_header,
        }
    }

    /// Install the pre-exchanged BGZF pool writers, keyed by destination path.
    /// Each is consumed when its destination is opened ([`create_all`] /
    /// [`write`]); the file is created by the exchange, so the pooled writer
    /// wraps it directly instead of opening a new file.
    pub fn set_pooled(&mut self, pooled: HashMap<PathBuf, PooledWriter>) {
        self.pooled = pooled
            .into_iter()
            .map(|(path, writer)| (Dest::from_path(Some(path.as_path())), writer))
            .collect();
    }

    /// Eagerly create the writer for each directed file `path`, so a directed
    /// output file always exists even when it receives no records (it is
    /// finalized to a valid empty file). A path that is later written to reuses
    /// the writer opened here; the lazy path in [`MultiWriter::write`] still
    /// covers any destination not listed.
    pub fn create_all(&mut self, paths: &[PathBuf]) -> Result<()> {
        for path in paths {
            let dest = Dest::from_path(Some(path));
            if !self.writers.contains_key(&dest) {
                let writer = self.open(Some(path))?;
                self.writers.insert(dest, writer);
            }
        }
        Ok(())
    }

    /// Write `reads` to the destination `path` (stdout when `None`), opening
    /// its writer on first use.
    pub fn write(&mut self, path: Option<&Path>, reads: &[OutputRead]) -> Result<()> {
        let dest = Dest::from_path(path);
        if !self.writers.contains_key(&dest) {
            let writer = self.open(path)?;
            self.writers.insert(dest.clone(), writer);
        }
        self.writers.get_mut(&dest).unwrap().write_fragment(reads)
    }

    /// Append a fragment's pre-encoded record bytes (produced on a worker by
    /// `encode_fragment` for this destination's format) to `path`, opening its
    /// writer on first use exactly as [`write`] does. The bytes are
    /// byte-identical to what [`write`] would have encoded inline.
    pub fn write_encoded(&mut self, path: Option<&Path>, bytes: &[u8]) -> Result<()> {
        let dest = Dest::from_path(path);
        if !self.writers.contains_key(&dest) {
            let writer = self.open(path)?;
            self.writers.insert(dest.clone(), writer);
        }
        self.writers.get_mut(&dest).unwrap().write_encoded(bytes)
    }

    /// Open the writer for a destination with the right format and header. A
    /// destination with a pre-exchanged pool writer is installed from there
    /// (its file already created by the exchange); every other destination
    /// opens its own file inline.
    fn open(&mut self, path: Option<&Path>) -> Result<OutputWriter> {
        let format = output_format(path, &self.input_formats)?;
        let dest = Dest::from_path(path);
        let header = self
            .headers
            .get(&dest)
            .cloned()
            .unwrap_or_else(|| self.default_header.clone());
        if let Some(pooled) = self.pooled.remove(&dest) {
            // The pooled BAM writer needs the @RG header; the pooled FASTX
            // writers ignore it.
            return OutputWriter::from_pooled(pooled, format, header);
        }
        OutputWriter::create_with_header(path, format, self.compression, header)
    }

    /// Finalize every opened writer (flush, gzip trailers, BGZF EOF), returning
    /// the first error. Every writer is finalized even when one fails, so a
    /// single bad sink does not leave the other fan-out files truncated.
    pub fn finish(self) -> Result<()> {
        let mut first_error = None;
        for (_, writer) in self.writers {
            if let Err(error) = writer.finish() {
                first_error.get_or_insert(error);
            }
        }
        first_error.map_or(Ok(()), Err)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::{OutputPattern, StreamRef};

    /// Parse a pattern through the grammar so tests use real `OutputPattern`
    /// values.
    fn pattern(text: &str, allowed: &[Placeholder]) -> OutputPattern {
        crate::grammar::parse_output_pattern(text, allowed).unwrap()
    }

    #[test]
    fn test_resolve_pattern_sample_and_sub_sample() {
        let pat = pattern(
            "out/%sample.%sub_sample.bam",
            &[Placeholder::Sample, Placeholder::SubSample],
        );
        let ctx = PathContext {
            pool: "lib01",
            sample: Some("dna01"),
            sub_sample: Some("lib01"),
            ..Default::default()
        };
        assert_eq!(
            resolve_pattern(&pat, &ctx),
            Some(PathBuf::from("out/dna01.lib01.bam"))
        );
    }

    #[test]
    fn test_resolve_pattern_missing_sub_sample_is_empty() {
        let pat = pattern(
            "out/%sample.%sub_sample.fq",
            &[Placeholder::Sample, Placeholder::SubSample],
        );
        let ctx = PathContext {
            pool: "p",
            sample: Some("s1"),
            sub_sample: None,
            ..Default::default()
        };
        // %sub_sample with no sub_sample collapses to empty, leaving the
        // doubled separator.
        assert_eq!(
            resolve_pattern(&pat, &ctx),
            Some(PathBuf::from("out/s1..fq"))
        );
    }

    #[test]
    fn test_resolve_pattern_ordinal_and_pool() {
        let pat = pattern(
            "out/%pool.R%ordinal.fq.gz",
            &[Placeholder::Pool, Placeholder::Ordinal],
        );
        let ctx = PathContext {
            pool: "FC1",
            ordinal: Some(2),
            ..Default::default()
        };
        assert_eq!(
            resolve_pattern(&pat, &ctx),
            Some(PathBuf::from("out/FC1.R2.fq.gz"))
        );
    }

    #[test]
    fn test_resolve_pattern_source_for_unassigned() {
        let pat = pattern("unmatched.%source.fq", &[Placeholder::Source]);
        let ctx = PathContext {
            pool: "p",
            source: Some(0),
            ..Default::default()
        };
        assert_eq!(
            resolve_pattern(&pat, &ctx),
            Some(PathBuf::from("unmatched.0.fq"))
        );
    }

    #[test]
    fn test_resolve_metrics_path_expands_pool() {
        // A metrics path is a pool-only pattern: `%pool` expands; the stdout
        // literals are not special-cased (a metrics TSV is always a real file).
        let pat = pattern("metrics/%pool.per_sample.tsv", &[Placeholder::Pool]);
        assert_eq!(
            resolve_metrics_path(&pat, "lib01"),
            PathBuf::from("metrics/lib01.per_sample.tsv")
        );
        let dash = pattern("-", &[]);
        assert_eq!(resolve_metrics_path(&dash, "lib01"), PathBuf::from("-"));
    }

    #[test]
    fn test_resolve_pattern_stdout_literals() {
        for literal in ["-", "/dev/stdout"] {
            let pat = pattern(literal, &[]);
            assert_eq!(resolve_pattern(&pat, &PathContext::default()), None);
        }
    }

    /// One extracted stream from its bases and optional qualities.
    fn extracted(bases: &[u8], quals: Option<&[u8]>) -> Extracted {
        Extracted {
            bases: bases.to_vec(),
            quals: quals.map(<[u8]>::to_vec),
            corrected: None,
        }
    }

    /// The string value of a tag in a Data, for assertions.
    fn tag_str(data: &Data, tag: &str) -> Option<String> {
        match data.get(&two_byte_tag(tag)) {
            Some(Value::String(value)) => Some(value.to_string()),
            _ => None,
        }
    }

    /// A `--tag` binding from its parts (streams given as names, optionally
    /// with a leading `~`).
    fn binding(tag: &str, streams: &[&str], qual: QualTag, sep: Option<&str>) -> TagBinding {
        TagBinding {
            tag: tag.to_string(),
            streams: streams
                .iter()
                .map(|s| match s.strip_prefix('~') {
                    Some(name) => StreamRef {
                        name: name.to_string(),
                        revcomp: true,
                    },
                    None => StreamRef {
                        name: s.to_string(),
                        revcomp: false,
                    },
                })
                .collect(),
            qual,
            sep: sep.map(str::to_string),
            qual_sep: sep.map(|_| " ".to_string()),
            raw: false,
        }
    }

    #[test]
    fn test_build_tag_data_corrected_vs_raw_bases_and_quals() {
        // A stream with a corrected form: a default (raw=false) tag emits the
        // corrected bases+quals; a raw=true tag emits the observed bases+quals.
        // Covers per-tag selection and the qual bytes.
        let mut stream = extracted(b"GAAGGG", Some(&[30, 30, 30, 30, 30, 30]));
        stream.corrected = Some(crate::extract::Corrected {
            bases: b"GAAGAG".to_vec(),
            quals: Some(vec![40, 40, 40, 40, 40, 40]),
        });
        let streams = HashMap::from([("bc", stream)]);

        let corrected = binding("CB", &["bc"], QualTag::Named("CY".to_string()), None);
        let data = build_tag_data(&[corrected], &streams);
        assert_eq!(
            tag_str(&data, "CB").as_deref(),
            Some("GAAGAG"),
            "default emits corrected"
        );
        // raw 40 + 33 = 'I'.
        assert_eq!(
            tag_str(&data, "CY").as_deref(),
            Some("IIIIII"),
            "corrected quals"
        );

        let mut raw = binding("CR", &["bc"], QualTag::Named("CY".to_string()), None);
        raw.raw = true;
        let data = build_tag_data(&[raw], &streams);
        assert_eq!(
            tag_str(&data, "CR").as_deref(),
            Some("GAAGGG"),
            "raw=true emits observed"
        );
        // raw 30 + 33 = '?'.
        assert_eq!(
            tag_str(&data, "CY").as_deref(),
            Some("??????"),
            "observed quals"
        );
    }

    #[test]
    fn test_build_tag_data_single_stream_default_qual() {
        // BC=bc: single stream, default qual tag QT carries the ASCII-encoded
        // qualities.
        let bc = binding("BC", &["bc"], QualTag::Default, None);
        let streams = HashMap::from([("bc", extracted(b"ACGT", Some(&[30, 31, 32, 33])))]);
        let data = build_tag_data(&[bc], &streams);
        assert_eq!(tag_str(&data, "BC").as_deref(), Some("ACGT"));
        // raw 30..=33 + 33 = '?','@','A','B'.
        assert_eq!(tag_str(&data, "QT").as_deref(), Some("?@AB"));
    }

    #[test]
    fn test_build_tag_data_joins_streams_with_sep() {
        // CB=a,b::sep=- joins bases with '-'; qual=none drops the quality tag.
        let cb = binding("CB", &["a", "b"], QualTag::None, Some("-"));
        let streams = HashMap::from([
            ("a", extracted(b"AAAA", Some(&[30; 4]))),
            ("b", extracted(b"TTTT", Some(&[30; 4]))),
        ]);
        let data = build_tag_data(&[cb], &streams);
        assert_eq!(tag_str(&data, "CB").as_deref(), Some("AAAA-TTTT"));
        assert_eq!(tag_str(&data, "CY"), None, "qual=none emits no quality tag");
    }

    #[test]
    fn test_build_tag_data_present_streams_only() {
        // CB=a,b::sep=-: only `a` matched, so CB is just its bases (no trailing
        // separator).
        let cb = binding("CB", &["a", "b"], QualTag::None, Some("-"));
        let streams = HashMap::from([("a", extracted(b"AAAA", Some(&[30; 4])))]);
        let data = build_tag_data(&[cb], &streams);
        assert_eq!(tag_str(&data, "CB").as_deref(), Some("AAAA"));
    }

    #[test]
    fn test_build_tag_data_no_present_stream_omits_tag() {
        let bc = binding("BC", &["bc"], QualTag::Default, None);
        let data = build_tag_data(&[bc], &HashMap::new());
        assert!(data.is_empty(), "an absent stream leaves a null tag");
    }

    #[test]
    fn test_read_group_header_lists_targets() {
        let targets = [
            Target {
                sample: "dna01".to_string(),
                sub_sample: Some("lib01".to_string()),
            },
            Target {
                sample: "dna02".to_string(),
                sub_sample: None,
            },
        ];
        let refs: Vec<&Target> = targets.iter().collect();
        let header = read_group_header(
            &refs,
            &[("PL".to_string(), "ILLUMINA".to_string())],
            Some("unmux --in 0=r.fq"),
        )
        .unwrap();
        // The @PG provenance record is present alongside the @RG lines.
        assert!(header.programs().as_ref().contains_key(&b"unmux"[..]));
        assert!(header.read_groups().contains_key(&b"dna01.lib01"[..]));
        assert!(header.read_groups().contains_key(&b"dna02"[..]));
        let rg = &header.read_groups()[&b"dna01.lib01"[..]];
        assert_eq!(
            rg.other_fields()
                .get(&rg_tag::SAMPLE)
                .map(|v| v.to_string()),
            Some("dna01".to_string())
        );
        assert_eq!(
            rg.other_fields()
                .get(&rg_tag::LIBRARY)
                .map(|v| v.to_string()),
            Some("lib01".to_string())
        );
    }

    #[test]
    fn test_default_read_group_header_uses_pool_for_id_sm_lb() {
        // The fallback header carries one @RG whose ID/SM/LB are all the pool
        // id, alongside the @PG provenance; shared --rg-tag fields validate
        // through the same path as `read_group_header`.
        let header =
            default_read_group_header("lib01", &[("PL".to_string(), "ILLUMINA".to_string())], None)
                .unwrap();
        let rg = header
            .read_groups()
            .get(&b"lib01"[..])
            .expect("@RG ID=pool");
        let field = |tag| rg.other_fields().get(tag).map(|v| v.to_string());
        assert_eq!(field(&rg_tag::SAMPLE).as_deref(), Some("lib01"), "SM=pool");
        assert_eq!(field(&rg_tag::LIBRARY).as_deref(), Some("lib01"), "LB=pool");
    }

    #[test]
    fn test_multiwriter_fans_out_to_separate_files() {
        // Two destinations get their own lazily-opened writer; each holds only
        // its reads.
        let dir = tempfile::tempdir().unwrap();
        let a = dir.path().join("a.fq");
        let b = dir.path().join("b.fq");
        let mut writer = MultiWriter::new(
            HashMap::new(),
            provenance_header(None),
            vec![SniffedFormat::Fastq],
            5,
        );
        let read = |name: &'static [u8]| OutputRead {
            name,
            bases: b"ACGT",
            quals: Some(&[30, 30, 30, 30]),
            tags: None,
            read_group: None,
        };
        writer.write(Some(&a), &[read(b"ra")]).unwrap();
        writer.write(Some(&b), &[read(b"rb")]).unwrap();
        writer.write(Some(&a), &[read(b"ra2")]).unwrap();
        writer.finish().unwrap();

        let names = |path: &std::path::Path| {
            let mut reader =
                crate::input::FragmentReader::open(&[path.to_path_buf()], false).unwrap();
            let mut names = Vec::new();
            while let Some(fragment) = reader.next_fragment().unwrap() {
                names.extend(fragment.records.into_iter().map(|r| r.name));
            }
            names
        };
        assert_eq!(names(&a), vec![b"ra".to_vec(), b"ra2".to_vec()]);
        assert_eq!(names(&b), vec![b"rb".to_vec()]);
    }
}
