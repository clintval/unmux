//! The demux output layer: resolve `--out`/`--unassigned`/`--remove` path
//! patterns to concrete [`Destination`]s, and assemble the `--tag` SAM fields
//! and `@RG` header for a fan-out.
//!
//! Routing decides *which* target a record belongs to (see [`crate::fanout`]);
//! this module turns that decision into files. A pattern carries placeholders
//! (`%pool`/`%sample`/`%sub_sample`/`%ordinal` for `--out`, `%pool`/`%source`
//! for the raw bins); [`resolve_pattern`] expands them against a
//! [`PathContext`]. The engine enumerates every directed destination and
//! creates each file up front, so a directed file always exists, even empty.

use std::collections::HashMap;
use std::path::{Component, Path, PathBuf};

use anyhow::{bail, Result};
use noodles::sam;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::header::record::value::map::{
    self, program::tag as program_tag, read_group::tag as rg_tag, Program, ReadGroup,
};
use noodles::sam::header::record::value::Map;

use crate::extract::Extracted;
use crate::fanout::Target;
use crate::grammar::{
    default_qual_tag, OutputPattern, PatternSegment, Placeholder, QualTag, StreamRef, TagBinding,
};
use crate::input::{SniffedFormat, PHRED_OFFSET};
use crate::writer::{file_id, output_format, FileId, OutputFormat};

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
/// streams. Each binding's sequence tag is its streams' bases joined by `sep`.
/// A stream with fallbacks (`a|b`) takes the first of them that is present.
/// A tag is written only when every one of its streams is present: a stream is
/// absent when its anchoring group did not match, and then the tag and its
/// quality tag are left out, so a partial match never yields a shorter value
/// that could equal another record's. An empty but present stream (a
/// between-anchors region with a missing anchor) keeps its empty slot. The
/// paired quality tag (from the default seq->qual map, an explicit `qual=`, or
/// off for `qual=none`) is the streams' qualities re-encoded to ASCII Phred+33
/// and joined by `qual-sep`; it is emitted only when every stream carries
/// qualities. `CB` and `CR` may both default to `CY`, which a record takes from
/// `CR` where it writes one, as the SAM spec pairs them, and from `CB` elsewhere.
pub fn build_tag_data(bindings: &[TagBinding], streams: &HashMap<&str, Extracted>) -> Data {
    let mut data = Data::default();
    for binding in bindings {
        // Each present stream is paired with its `StreamRef` so a `~`
        // (reverse-complement at this use site) can flip that stream's bases
        // and qualities.
        let Some(present) = chosen_streams(binding, streams) else {
            continue;
        };
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
            let left_to_cr = binding.tag == "CB"
                && binding.qual == QualTag::Default
                && data.get(&two_byte_tag(qual_tag)).is_some();
            if !left_to_cr
                && present
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

/// The stream chosen for each part of a tag, the first present of its choices,
/// or `None` when some part has none present.
fn chosen_streams<'a>(
    binding: &'a TagBinding,
    streams: &'a HashMap<&str, Extracted>,
) -> Option<Vec<(&'a StreamRef, &'a Extracted)>> {
    binding
        .streams
        .iter()
        .zip(&binding.fallbacks)
        .map(|(first, fallbacks)| {
            std::iter::once(first)
                .chain(fallbacks)
                .find_map(|s| streams.get(s.name.as_str()).map(|e| (s, e)))
        })
        .collect()
}

/// Whether a record writes the tag `binding` from a stream without base
/// qualities, such as column text, so no quality tag can describe it.
pub fn lacks_qualities(binding: &TagBinding, streams: &HashMap<&str, Extracted>) -> bool {
    chosen_streams(binding, streams).is_some_and(|chosen| {
        chosen
            .iter()
            .any(|(_, stream)| stream.tag_quals(binding.raw).is_none())
    })
}

/// The quality tag name a binding emits: the explicit `qual=QTAG`, the
/// seq->qual default, or `None` for `qual=none` (and for a sequence tag with no
/// default mapping).
pub(crate) fn qual_tag_name(binding: &TagBinding) -> Option<&str> {
    match &binding.qual {
        QualTag::None | QualTag::Omitted => None,
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
/// Whether a SAM header field can hold `text` as it is.
pub(crate) fn is_header_text(text: &str) -> bool {
    !text.is_empty() && text.chars().all(|c| matches!(c, ' '..='~'))
}

/// `text` with every character that a SAM header field cannot hold escaped.
fn header_text(text: &str) -> String {
    text.chars()
        .map(|c| match c {
            ' '..='~' => c.to_string(),
            _ => c.escape_default().to_string(),
        })
        .collect()
}

fn program_record(command_line: Option<&str>) -> Map<Program> {
    let mut builder = Map::<Program>::builder()
        .insert(program_tag::NAME, "unmux")
        .insert(program_tag::VERSION, env!("CARGO_PKG_VERSION"));
    if let Some(command_line) = command_line {
        builder = builder.insert(program_tag::COMMAND_LINE, header_text(command_line));
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

/// A `Target -> constituent original @RG maps` table: for each fan-out target,
/// the original input `@RG` lines (RG id plus its verbatim `Map<ReadGroup>`)
/// whose records land in that target's output file. Built once per `@RG`-split
/// run and consumed by [`read_group_header`].
pub type ConstituentReadGroups = HashMap<Target, Vec<(Vec<u8>, Map<ReadGroup>)>>;

/// Build a SAM header carrying the `@PG` provenance plus the `@RG` lines for the
/// targets whose records land in this output file.
///
/// With `constituent` `Some` (an `@RG`-split run), each target contributes its
/// *original* input `@RG` line(s) verbatim, deduped by id across targets that
/// share the file, so the split is faithful to the source header (correct
/// `SM`/`LB`/`PL`/`PU`/...); the synthesized per-target `@RG` path is skipped.
/// The shared `--rg-tag` overlay is not applied here: the originals already
/// carry their fields, so applying `--rg-tag` onto them is deferred.
///
/// With `constituent` `None` (normal barcode demux), one synthesized `@RG` is
/// emitted per target (`ID` = the target label, `SM` = the sample, `LB` = the
/// sub_sample) with the shared `--rg-tag` fields on every read group.
pub fn read_group_header(
    targets: &[&Target],
    rg_tags: &[(String, String)],
    command_line: Option<&str>,
    constituent: Option<&ConstituentReadGroups>,
) -> Result<sam::Header> {
    let mut builder = provenance_header_builder(command_line);
    if let Some(constituent) = constituent {
        for (id, map) in file_read_groups(targets, constituent) {
            builder = builder.add_read_group(id, map);
        }
        return Ok(builder.build());
    }
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

/// The original `@RG` lines of the targets sharing one file, deduped and sorted by id.
fn file_read_groups(
    targets: &[&Target],
    constituent: &ConstituentReadGroups,
) -> Vec<(Vec<u8>, Map<ReadGroup>)> {
    let mut seen: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
    let mut deduped: Vec<(Vec<u8>, Map<ReadGroup>)> = Vec::new();
    for target in targets {
        for (id, map) in constituent.get(target).into_iter().flatten() {
            if seen.insert(id.clone()) {
                deduped.push((id.clone(), map.clone()));
            }
        }
    }
    deduped.sort_by(|(a, _), (b, _)| a.cmp(b));
    deduped
}

/// The input header with the unmux `@PG` appended, its `@RG` lines narrowed (in
/// input order) to the targets sharing this file when `targets` is given.
pub fn alignment_header(
    input: &sam::Header,
    targets: Option<(&[&Target], &ConstituentReadGroups)>,
    command_line: Option<&str>,
) -> Result<sam::Header> {
    let mut header = input.clone();
    if let Some((targets, constituent)) = targets {
        let keep: std::collections::HashSet<Vec<u8>> = file_read_groups(targets, constituent)
            .into_iter()
            .map(|(id, _)| id)
            .collect();
        header
            .read_groups_mut()
            .retain(|id, _| keep.contains(AsRef::<[u8]>::as_ref(id)));
    }
    append_program(header.programs_mut(), program_record(command_line));
    Ok(header)
}

/// Append the unmux `@PG` under the first free ID of `unmux`, `unmux.1`, ...,
/// chained by `PP` to the last leaf program unless the existing chains have a
/// dangling `PP` or a cycle, in which case it is left unchained.
fn append_program(programs: &mut sam::header::Programs, mut map: Map<Program>) {
    let inner = programs.as_mut();
    let previous = |id: &[u8]| {
        inner
            .get(id)
            .and_then(|map| map.other_fields().get(&program_tag::PREVIOUS_PROGRAM_ID))
            .map(|pp| AsRef::<[u8]>::as_ref(pp).to_vec())
    };
    let chains_are_sound = inner.keys().all(|start| {
        let mut seen = std::collections::HashSet::new();
        let mut current = AsRef::<[u8]>::as_ref(start).to_vec();
        while let Some(pp) = previous(&current) {
            if !inner.contains_key(pp.as_slice()) || !seen.insert(pp.clone()) {
                return false;
            }
            current = pp;
        }
        true
    });
    if chains_are_sound {
        let referenced: std::collections::HashSet<Vec<u8>> = inner
            .keys()
            .filter_map(|id| previous(id.as_ref()))
            .collect();
        let leaf = inner
            .keys()
            .rev()
            .find(|id| !referenced.contains(AsRef::<[u8]>::as_ref(*id)))
            .cloned();
        if let Some(leaf) = leaf {
            map.other_fields_mut()
                .insert(program_tag::PREVIOUS_PROGRAM_ID, leaf);
        }
    }
    let id = (0..)
        .map(|n| match n {
            0 => "unmux".to_string(),
            n => format!("unmux.{n}"),
        })
        .find(|id| !inner.contains_key(id.as_bytes()))
        .expect("an unused program ID exists");
    inner.insert(id.into(), map);
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

/// One output of a run: a file, or stdout when `path` is `None`, in `format`,
/// opened with `header` when it is SAM, BAM, or CRAM.
#[derive(Debug, Clone)]
pub struct Destination {
    /// The output path, or `None` for stdout.
    pub path: Option<PathBuf>,
    /// The format the path resolves to.
    pub format: OutputFormat,
    /// The header of a SAM, BAM, or CRAM output.
    pub header: sam::Header,
}

impl Destination {
    /// The destination at `path` (stdout when `None`), with the header in
    /// `headers` for that path or else `default_header`.
    pub fn new(
        path: Option<PathBuf>,
        input_formats: &[SniffedFormat],
        headers: &HashMap<Option<PathBuf>, sam::Header>,
        default_header: &sam::Header,
    ) -> Result<Self> {
        let format = output_format(path.as_deref(), input_formats)?;
        let header = headers.get(&path).unwrap_or(default_header).clone();
        Ok(Self {
            path,
            format,
            header,
        })
    }
}

/// What writes one output file, reported when two outputs would be one file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutputOrigin {
    /// An `--out` file of one fan-out target.
    Target(Target),
    /// Records sharing the `--out` writer: `--out` with no target, or a bin given the same pattern.
    Records(String),
    /// A bin such as `--unassigned` or `--remove grp` given a pattern of its own.
    Bin(String),
    /// A metrics TSV written on its own after the run, such as `--metrics-summary`.
    Metrics(&'static str),
}

impl std::fmt::Display for OutputOrigin {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Target(Target {
                sample,
                sub_sample: Some(sub_sample),
            }) => write!(f, "sample {sample}, sub_sample {sub_sample}"),
            Self::Target(target) => write!(f, "sample {}", target.sample),
            Self::Records(option) | Self::Bin(option) => f.write_str(option),
            Self::Metrics(option) => f.write_str(option),
        }
    }
}

/// Refuse outputs that would be one file: paths differing only in case, or two spellings of one.
pub fn check_output_collisions(outputs: &[(PathBuf, OutputOrigin)]) -> Result<()> {
    let mut dirs: HashMap<PathBuf, PathBuf> = HashMap::new();
    let mut by_name: HashMap<String, usize> = HashMap::new();
    let mut by_file: HashMap<FileId, usize> = HashMap::new();
    for (index, (path, origin)) in outputs.iter().enumerate() {
        let Some((name, file)) = file_key(path, &mut dirs) else {
            continue;
        };
        let earlier = by_name
            .get(&name)
            .or_else(|| file.and_then(|file| by_file.get(&file)))
            .copied();
        by_name.entry(name).or_insert(index);
        if let Some(file) = file {
            by_file.entry(file).or_insert(index);
        }
        let Some((first, first_origin)) = earlier.map(|earlier| &outputs[earlier]) else {
            continue;
        };
        if first == path && shares_writer(first_origin, origin) {
            continue;
        }
        let (lexical_a, lexical_b) = (lexical(first), lexical(path));
        let case_only = lexical_a != lexical_b
            && fold(&lexical_a.to_string_lossy()) == fold(&lexical_b.to_string_lossy());
        let advice = match (first_origin, origin) {
            (OutputOrigin::Target(a), OutputOrigin::Target(b)) if case_only => {
                if a.sample == b.sample {
                    "rename one sub_sample"
                } else {
                    "rename one sample"
                }
            }
            _ => "give them distinct paths",
        };
        let reason = if case_only {
            "differ only in case, so they would be one file on a case-insensitive file system"
        } else {
            "are the same file"
        };
        let (a, b) = (first.display(), path.display());
        bail!("outputs \"{a}\" ({first_origin}) and \"{b}\" ({origin}) {reason}; {advice}");
    }
    Ok(())
}

/// Whether two outputs on one path share its writer: both on the `--out` writer, or both bins.
fn shares_writer(a: &OutputOrigin, b: &OutputOrigin) -> bool {
    let out = |o: &OutputOrigin| matches!(o, OutputOrigin::Target(_) | OutputOrigin::Records(_));
    let bin = |o: &OutputOrigin| matches!(o, OutputOrigin::Bin(_));
    (out(a) && out(b)) || (bin(a) && bin(b))
}

/// A file the run reads, reported when an output would replace it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InputOrigin {
    /// A record input, by its 0-based index.
    Reads(usize),
    /// The tag file of a `--group`.
    TagFile(String),
    /// The `--sample-sheet`.
    SampleSheet,
    /// The `--reference` FASTA.
    Reference,
    /// An index of the `--reference` FASTA.
    ReferenceIndex,
}

impl std::fmt::Display for InputOrigin {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Reads(index) => write!(f, "input {index}"),
            Self::TagFile(group) => write!(f, "the tag file of group {group}"),
            Self::SampleSheet => f.write_str("the sample sheet"),
            Self::Reference => f.write_str("the reference"),
            Self::ReferenceIndex => f.write_str("the reference index"),
        }
    }
}

/// Refuse an output that is a file the run reads under any name, which would replace that file.
pub fn check_outputs_over_inputs(
    outputs: &[(PathBuf, OutputOrigin)],
    inputs: &[(PathBuf, InputOrigin)],
) -> Result<()> {
    let files: Vec<_> = inputs
        .iter()
        .filter_map(|(path, origin)| {
            if matches!(origin, InputOrigin::Reads(_)) && path.as_os_str() == "-" {
                return Some(("<stdin>".to_string(), origin, None, Some(stdin_file()?)));
            }
            let (real, file) = read_file(path)?;
            Some((path.display().to_string(), origin, Some(real), file))
        })
        .collect();
    if files.is_empty() {
        return Ok(());
    }
    for (path, origin) in outputs {
        let real = resolve(path, &mut 0);
        let file = match std::fs::metadata(&real) {
            Ok(metadata) if !metadata.is_file() => continue,
            Ok(metadata) => file_id(&metadata),
            Err(_) => None,
        };
        let Some((input, input_origin, ..)) = files.iter().find(|(_, _, input, id)| {
            input.as_ref() == Some(&real) || id.is_some_and(|id| file == Some(id))
        }) else {
            continue;
        };
        let output = path.display();
        bail!(
            "output \"{output}\" ({origin}) is {input_origin} ({input}); give it a distinct path"
        );
    }
    Ok(())
}

/// The resolved path and file of a regular file the run reads, or `None` for a pipe or device.
fn read_file(path: &Path) -> Option<(PathBuf, Option<FileId>)> {
    let metadata = std::fs::metadata(path)
        .ok()
        .filter(std::fs::Metadata::is_file)?;
    // An open descriptor gives the device of a file behind /dev/fd, which stat does not on macOS.
    let opened = std::fs::File::open(path)
        .and_then(|file| file.metadata())
        .ok();
    let file = opened.as_ref().map_or_else(|| file_id(&metadata), file_id);
    Some((resolve(path, &mut 0), file))
}

/// The file standard input reads, when it is a regular file.
#[cfg(unix)]
fn stdin_file() -> Option<FileId> {
    use std::os::fd::AsFd;
    let stdin = std::io::stdin().as_fd().try_clone_to_owned().ok()?;
    let metadata = std::fs::File::from(stdin).metadata().ok()?;
    file_id(&metadata).filter(|_| metadata.is_file())
}

#[cfg(not(unix))]
fn stdin_file() -> Option<FileId> {
    None
}

/// The case-folded absolute path `path` writes and its existing file, or `None` if not regular.
fn file_key(path: &Path, dirs: &mut HashMap<PathBuf, PathBuf>) -> Option<(String, Option<FileId>)> {
    let resolved = match path.file_name() {
        Some(name) => {
            let parent = path.parent().unwrap_or(Path::new(""));
            let candidate = dirs
                .entry(parent.to_path_buf())
                .or_insert_with(|| resolve(parent, &mut 0))
                .join(name);
            match std::fs::symlink_metadata(&candidate) {
                Ok(metadata) if metadata.file_type().is_symlink() => resolve(&candidate, &mut 0),
                _ => candidate,
            }
        }
        None => resolve(path, &mut 0),
    };
    let file = match std::fs::metadata(&resolved) {
        Ok(metadata) if !metadata.is_file() => return None,
        Ok(metadata) => file_id(&metadata),
        Err(_) => None,
    };
    Some((fold(&resolved.to_string_lossy()), file))
}

/// The absolute form of `path`, following every symbolic link in it, even a dangling one.
fn resolve(path: &Path, hops: &mut u32) -> PathBuf {
    let mut real = if path.is_absolute() {
        PathBuf::new()
    } else {
        std::env::current_dir().unwrap_or_default()
    };
    for part in path.components() {
        match part {
            Component::ParentDir => {
                real.pop();
            }
            Component::CurDir => {}
            Component::Normal(name) => {
                let next = real.join(name);
                real = match std::fs::read_link(&next) {
                    Ok(target) if *hops < 40 => {
                        *hops += 1;
                        resolve(&real.join(target), hops)
                    }
                    _ => next,
                };
            }
            root => real.push(root),
        }
    }
    real
}

/// The components of `path` with `.` dropped and each `..` removing the name before it.
fn lexical(path: &Path) -> PathBuf {
    let mut out = PathBuf::new();
    for part in path.components() {
        match part {
            Component::CurDir => {}
            Component::ParentDir
                if matches!(out.components().next_back(), Some(Component::Normal(_))) =>
            {
                out.pop();
            }
            part => out.push(part),
        }
    }
    out
}

/// Fold case per character through upper case, so `ß`/`ss`, `ſ`/`s` and `ς`/`σ` compare equal.
fn fold(text: &str) -> String {
    text.chars()
        .flat_map(char::to_uppercase)
        .flat_map(char::to_lowercase)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::OutputPattern;

    /// Parse a pattern through the grammar so tests use real `OutputPattern`
    /// values.
    fn pattern(text: &str, allowed: &[Placeholder]) -> OutputPattern {
        crate::grammar::parse_output_pattern(text, allowed).unwrap()
    }

    #[test]
    fn test_alignment_header_keeps_input_and_narrows_read_groups() {
        let input = sam::io::Reader::new(
            &b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
               @RG\tID:rg1\tSM:sA\n@RG\tID:rg2\tSM:sB\n@PG\tID:bwa\tPN:bwa\n@CO\tnote\n"[..],
        )
        .read_header()
        .unwrap();
        let target = Target {
            sample: "sA".to_string(),
            sub_sample: None,
        };
        let rg1 = input.read_groups().get(&b"rg1"[..]).unwrap().clone();
        let constituent: ConstituentReadGroups =
            HashMap::from([(target.clone(), vec![(b"rg1".to_vec(), rg1)])]);
        let render = |header: &sam::Header| {
            let mut buf = Vec::new();
            sam::io::Writer::new(&mut buf).write_header(header).unwrap();
            String::from_utf8(buf).unwrap()
        };

        let narrowed =
            alignment_header(&input, Some((&[&target], &constituent)), Some("unmux x")).unwrap();
        let text = render(&narrowed);
        assert!(text.starts_with("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n"));
        assert!(text.contains("@RG\tID:rg1\tSM:sA\n"));
        assert!(!text.contains("ID:rg2"));
        assert!(text.contains("@PG\tID:bwa\tPN:bwa\n"));
        assert!(text.contains("@PG\tID:unmux\tPN:unmux"));
        assert!(text.contains("PP:bwa"));
        assert!(text.contains("@CO\tnote"));

        let full = render(&alignment_header(&input, None, None).unwrap());
        assert!(full.contains("ID:rg1") && full.contains("ID:rg2"));
    }

    fn appended(text: &[u8]) -> Vec<(String, Option<String>)> {
        let input = sam::io::Reader::new(text).read_header().unwrap();
        let header = alignment_header(&input, None, None).unwrap();
        header
            .programs()
            .as_ref()
            .iter()
            .map(|(id, map)| {
                let pp = map
                    .other_fields()
                    .get(&program_tag::PREVIOUS_PROGRAM_ID)
                    .map(|pp| pp.to_string());
                (id.to_string(), pp)
            })
            .collect()
    }

    #[test]
    fn test_unmux_program_chains_to_the_last_leaf() {
        let programs =
            appended(b"@HD\tVN:1.6\n@PG\tID:a\tPN:a\n@PG\tID:b\tPN:b\tPP:a\n@PG\tID:c\tPN:c\n");
        assert_eq!(
            programs.last().unwrap(),
            &("unmux".to_string(), Some("c".to_string()))
        );
        assert_eq!(programs.len(), 4, "one unmux @PG even with two leaves");
    }

    #[test]
    fn test_unmux_program_is_unchained_after_a_dangling_link_or_cycle() {
        let dangling = appended(b"@HD\tVN:1.6\n@PG\tID:a\tPN:a\tPP:missing\n");
        assert_eq!(dangling.last().unwrap(), &("unmux".to_string(), None));
        let cycle = appended(b"@HD\tVN:1.6\n@PG\tID:a\tPN:a\tPP:b\n@PG\tID:b\tPN:b\tPP:a\n");
        assert_eq!(cycle.last().unwrap(), &("unmux".to_string(), None));
    }

    #[test]
    fn test_unmux_program_id_takes_a_counter_suffix() {
        let programs = appended(
            b"@HD\tVN:1.6\n@PG\tID:unmux\tPN:unmux\n@PG\tID:unmux.1\tPN:unmux\tPP:unmux\n",
        );
        assert_eq!(
            programs.last().unwrap(),
            &("unmux.2".to_string(), Some("unmux.1".to_string()))
        );
    }

    #[test]
    fn test_alignment_header_dedups_a_repeated_unmux_program_id() {
        let input = sam::io::Reader::new(&b"@HD\tVN:1.6\n@PG\tID:unmux\tPN:unmux\n"[..])
            .read_header()
            .unwrap();
        let header = alignment_header(&input, None, None).unwrap();
        assert_eq!(header.programs().as_ref().len(), 2);
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
            fallbacks: vec![Vec::new(); streams.len()],
            qual,
            sep: sep.map(str::to_string),
            qual_sep: sep.map(|_| " ".to_string()),
            qual_sep_set: false,
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
    fn test_build_tag_data_one_absent_stream_omits_tag_and_qual_tag() {
        // CB=a+b+c: `b` is absent, so neither CB nor CY is written rather than
        // a CB of `AAAA-TTTT` that a full match of other streams could share.
        let streams = HashMap::from([
            ("a", extracted(b"AAAA", Some(&[30; 4]))),
            ("c", extracted(b"TTTT", Some(&[30; 4]))),
        ]);
        for cb in [
            binding("CB", &["a", "b", "c"], QualTag::Default, Some("-")),
            binding("CB", &["a", "~b", "c"], QualTag::Default, Some("-")),
        ] {
            let data = build_tag_data(&[cb], &streams);
            assert!(data.is_empty(), "{data:?}");
        }
        let mut raw = binding("CR", &["a", "b", "c"], QualTag::Default, Some("-"));
        raw.raw = true;
        assert!(build_tag_data(&[raw], &streams).is_empty());
    }

    #[test]
    fn test_build_tag_data_one_absent_stream_leaves_other_tags() {
        let cb = binding("CB", &["a", "b"], QualTag::Default, Some("-"));
        let rx = binding("RX", &["a"], QualTag::Default, None);
        let streams = HashMap::from([("a", extracted(b"AAAA", Some(&[30; 4])))]);
        let data = build_tag_data(&[cb, rx], &streams);
        assert_eq!(tag_str(&data, "CB"), None);
        assert_eq!(tag_str(&data, "CY"), None);
        assert_eq!(tag_str(&data, "RX").as_deref(), Some("AAAA"));
        assert_eq!(tag_str(&data, "QX").as_deref(), Some("????"));
    }

    #[test]
    fn test_build_tag_data_empty_present_stream_keeps_its_slot() {
        // An empty but present stream is not absent: its empty slot is joined.
        let cb = binding("CB", &["a", "b", "c"], QualTag::Default, Some("-"));
        let streams = HashMap::from([
            ("a", extracted(b"AAAA", Some(&[30; 4]))),
            ("b", extracted(b"", Some(&[]))),
            ("c", extracted(b"TTTT", Some(&[30; 4]))),
        ]);
        let data = build_tag_data(&[cb], &streams);
        assert_eq!(tag_str(&data, "CB").as_deref(), Some("AAAA--TTTT"));
        assert_eq!(tag_str(&data, "CY").as_deref(), Some("????  ????"));
    }

    #[test]
    fn test_build_tag_data_no_present_stream_omits_tag() {
        let bc = binding("BC", &["bc"], QualTag::Default, None);
        let data = build_tag_data(&[bc], &HashMap::new());
        assert!(data.is_empty(), "an absent stream leaves a null tag");
    }

    #[test]
    fn test_build_tag_data_cy_holds_the_cr_qualities_where_cr_writes_them() {
        let cb = binding("CB", &["u"], QualTag::Default, None);
        let mut cr = binding("CR", &["v"], QualTag::Default, None);
        cr.fallbacks = vec![vec![StreamRef {
            name: "n".to_string(),
            revcomp: false,
        }]];
        let tags = |bindings: &[TagBinding], present: &[&str]| {
            let streams: HashMap<&str, Extracted> = [
                ("u", extracted(b"AAAA", Some(&[40; 4]))),
                ("v", extracted(b"CCCC", Some(&[20; 4]))),
                ("n", extracted(b"cell", None)),
            ]
            .into_iter()
            .filter(|(name, _)| present.contains(name))
            .collect();
            build_tag_data(bindings, &streams)
                .iter()
                .map(|(tag, value)| {
                    let tag = String::from_utf8(<[u8; 2]>::from(tag).to_vec()).unwrap();
                    match value {
                        Value::String(text) => format!("{tag}:{text}"),
                        _ => tag,
                    }
                })
                .collect::<Vec<_>>()
                .join(" ")
        };
        let cb_first = [cb.clone(), cr.clone()];
        let cr_first = [cr, cb];
        assert_eq!(tags(&cb_first, &["u", "v"]), "CB:AAAA CY:5555 CR:CCCC");
        assert_eq!(tags(&cr_first, &["u", "v"]), "CR:CCCC CY:5555 CB:AAAA");
        assert_eq!(tags(&cb_first, &["u"]), "CB:AAAA CY:IIII");
        assert_eq!(tags(&cr_first, &["u"]), "CB:AAAA CY:IIII");
        assert_eq!(tags(&cb_first, &["u", "n"]), "CB:AAAA CY:IIII CR:cell");
        assert_eq!(tags(&cr_first, &["u", "n"]), "CR:cell CB:AAAA CY:IIII");
        assert_eq!(tags(&cb_first, &["v"]), "CR:CCCC CY:5555");
        assert_eq!(tags(&cr_first, &["v"]), "CR:CCCC CY:5555");
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
            None,
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
    fn read_group_header_emits_only_constituent_lines() {
        // An @RG-split run: sampleA's file carries the original @RG lines for
        // rg1 and rg3 (the read groups mapping to it), and never sampleB's rg2.
        let sample_a = Target {
            sample: "sampleA".to_string(),
            sub_sample: None,
        };
        let sample_b = Target {
            sample: "sampleB".to_string(),
            sub_sample: None,
        };
        let rg_map = |sm: &str| {
            Map::<ReadGroup>::builder()
                .insert(rg_tag::SAMPLE, sm)
                .build()
                .expect("a read group map has no required fields")
        };
        let mut constituent: ConstituentReadGroups = HashMap::new();
        constituent.insert(
            sample_a.clone(),
            vec![
                (b"rg1".to_vec(), rg_map("sampleA")),
                (b"rg3".to_vec(), rg_map("sampleA")),
            ],
        );
        constituent.insert(sample_b.clone(), vec![(b"rg2".to_vec(), rg_map("sampleB"))]);

        let header = read_group_header(&[&sample_a], &[], None, Some(&constituent)).unwrap();
        let read_groups = header.read_groups();
        assert!(read_groups.contains_key(&b"rg1"[..]), "carries rg1");
        assert!(read_groups.contains_key(&b"rg3"[..]), "carries rg3");
        assert!(
            !read_groups.contains_key(&b"rg2"[..]),
            "never another sample's rg2"
        );
        // The original SM is preserved verbatim (not the synthesized target
        // label).
        let rg1 = &read_groups[&b"rg1"[..]];
        assert_eq!(
            rg1.other_fields()
                .get(&rg_tag::SAMPLE)
                .map(|v| v.to_string()),
            Some("sampleA".to_string())
        );
    }

    #[test]
    fn read_group_header_emits_constituent_lines_in_deterministic_id_order() {
        // The constituent bucket for sampleA is built in reverse id order
        // (rg3 before rg1). Regardless of that bucket order, the emitted
        // @RG lines must come out sorted by id, so the header is
        // deterministic run-to-run even though `build_constituent_rgs`
        // iterates a HashMap upstream.
        let sample_a = Target {
            sample: "sampleA".to_string(),
            sub_sample: None,
        };
        let rg_map = |sm: &str| {
            Map::<ReadGroup>::builder()
                .insert(rg_tag::SAMPLE, sm)
                .build()
                .expect("a read group map has no required fields")
        };
        let mut constituent: ConstituentReadGroups = HashMap::new();
        constituent.insert(
            sample_a.clone(),
            vec![
                (b"rg3".to_vec(), rg_map("sampleA")),
                (b"rg1".to_vec(), rg_map("sampleA")),
            ],
        );

        let header = read_group_header(&[&sample_a], &[], None, Some(&constituent)).unwrap();
        let mut writer = sam::io::Writer::new(Vec::new());
        writer.write_header(&header).unwrap();
        let text = String::from_utf8(writer.into_inner()).unwrap();

        let rg1_pos = text.find("ID:rg1").expect("rg1 present in header text");
        let rg3_pos = text.find("ID:rg3").expect("rg3 present in header text");
        assert!(
            rg1_pos < rg3_pos,
            "expected ID:rg1 before ID:rg3 in header text, got:\n{text}"
        );
    }

    #[test]
    fn test_header_text_escapes_what_a_sam_header_cannot_hold() {
        assert_eq!(header_text("a b ~ \\x"), "a b ~ \\x");
        assert_eq!(header_text("\u{e4}\t\n"), "\\u{e4}\\t\\n");
        assert!(is_header_text("a b~"));
        for text in ["", "\u{e4}", "a\tb", "a\u{7f}"] {
            assert!(!is_header_text(text), "{text:?}");
        }
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

    fn sample(sample: &str, sub_sample: Option<&str>) -> OutputOrigin {
        OutputOrigin::Target(Target {
            sample: sample.to_string(),
            sub_sample: sub_sample.map(str::to_string),
        })
    }

    fn collision(outputs: &[(&str, OutputOrigin)]) -> Option<String> {
        let outputs: Vec<(PathBuf, OutputOrigin)> = outputs
            .iter()
            .map(|(path, origin)| (PathBuf::from(path), origin.clone()))
            .collect();
        check_output_collisions(&outputs)
            .err()
            .map(|e| e.to_string())
    }

    #[test]
    fn test_collision_paths_differing_only_in_case_are_refused() {
        let error = collision(&[
            ("o/s1.fq", sample("s1", None)),
            ("o/S1.fq", sample("S1", None)),
        ]);
        assert_eq!(
            error.as_deref(),
            Some(
                "outputs \"o/s1.fq\" (sample s1) and \"o/S1.fq\" (sample S1) differ only in case, \
                 so they would be one file on a case-insensitive file system; rename one sample"
            )
        );
    }

    #[test]
    fn test_collision_case_folding_covers_directories_and_non_ascii() {
        let error = collision(&[
            ("Run/a.fq", sample("a", None)),
            (
                "run/a.fq",
                OutputOrigin::Records("--unassigned".to_string()),
            ),
        ])
        .unwrap();
        assert!(error.contains("differ only in case"), "{error}");
        assert!(error.contains("give them distinct paths"), "{error}");
        let error = collision(&[
            ("o/\u{c9}.fq", sample("\u{c9}", None)),
            ("o/\u{e9}.fq", sample("\u{e9}", None)),
        ])
        .unwrap();
        assert!(error.contains("differ only in case"), "{error}");
    }

    #[test]
    fn test_collision_sub_samples_differing_only_in_case_are_refused() {
        let error = collision(&[
            ("o/a.x.fq", sample("a", Some("x"))),
            ("o/a.X.fq", sample("a", Some("X"))),
        ])
        .unwrap();
        assert!(
            error.contains("(sample a, sub_sample x) and \"o/a.X.fq\" (sample a, sub_sample X)"),
            "{error}"
        );
        assert!(error.ends_with("rename one sub_sample"), "{error}");
    }

    #[test]
    fn test_collision_identical_record_paths_share_one_writer() {
        assert_eq!(
            collision(&[
                ("all.fq", sample("a", None)),
                ("all.fq", sample("b", None)),
                ("all.fq", OutputOrigin::Records("--unassigned".to_string())),
            ]),
            None
        );
    }

    #[test]
    fn test_collision_bin_given_its_own_pattern_may_not_share_an_out_path() {
        let bin = OutputOrigin::Bin("--unassigned".to_string());
        for out in [
            sample("a", None),
            OutputOrigin::Records("--out".to_string()),
        ] {
            for (first, second) in [(&out, &bin), (&bin, &out)] {
                let error = collision(&[("o/a.fq", first.clone()), ("o/a.fq", second.clone())]);
                assert_eq!(
                    error,
                    Some(format!(
                        "outputs \"o/a.fq\" ({first}) and \"o/a.fq\" ({second}) are the same file; \
                         give them distinct paths"
                    ))
                );
            }
        }
    }

    #[test]
    fn test_collision_bins_given_their_own_patterns_share_one_writer() {
        assert_eq!(
            collision(&[
                ("o/bins.fq", OutputOrigin::Bin("--unassigned".to_string())),
                ("o/bins.fq", OutputOrigin::Bin("--remove grp".to_string())),
            ]),
            None
        );
    }

    #[test]
    fn test_collision_metrics_may_not_share_a_path() {
        let error = collision(&[
            ("o/a.fq", sample("a", None)),
            ("o/a.fq", OutputOrigin::Metrics("--metrics-summary")),
        ]);
        assert_eq!(
            error.as_deref(),
            Some(
                "outputs \"o/a.fq\" (sample a) and \"o/a.fq\" (--metrics-summary) are the same \
                 file; give them distinct paths"
            )
        );
        let error = collision(&[
            ("m.tsv", OutputOrigin::Metrics("--metrics-per-sample")),
            ("M.tsv", OutputOrigin::Metrics("--metrics-summary")),
        ])
        .unwrap();
        assert!(error.contains("differ only in case"), "{error}");
    }

    #[test]
    fn test_collision_distinct_paths_and_directories_pass() {
        assert_eq!(
            collision(&[
                ("o/s1.fq", sample("s1", None)),
                ("o/s2.fq", sample("s2", None)),
                ("p/s1.fq", OutputOrigin::Records("--unassigned".to_string())),
                (
                    "o/s1.fq.gz",
                    OutputOrigin::Records("--remove grp".to_string())
                ),
                ("o/m.tsv", OutputOrigin::Metrics("--metrics-summary")),
            ]),
            None
        );
    }

    #[test]
    fn test_collision_different_spellings_of_one_path_are_refused() {
        for (a, b) in [("o/a.fq", "./o/a.fq"), ("o/a.fq", "o/x/../a.fq")] {
            let error = collision(&[
                (a, sample("a", None)),
                (b, OutputOrigin::Records("--unassigned".to_string())),
            ])
            .unwrap_or_else(|| panic!("{a} and {b} were accepted"));
            assert!(error.contains("are the same file"), "{error}");
        }
    }

    #[cfg(unix)]
    #[test]
    fn test_collision_symbolic_link_to_another_output_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let real = dir.path().join("real");
        std::fs::create_dir(&real).unwrap();
        std::os::unix::fs::symlink(&real, dir.path().join("link")).unwrap();
        std::fs::write(real.join("b.fq"), "").unwrap();
        std::os::unix::fs::symlink(real.join("b.fq"), real.join("to-b.fq")).unwrap();
        for (a, b) in [("real/a.fq", "link/a.fq"), ("real/b.fq", "real/to-b.fq")] {
            let (a, b) = (dir.path().join(a), dir.path().join(b));
            let outputs = [
                (a.clone(), sample("a", None)),
                (b.clone(), sample("b", None)),
            ];
            let error = check_output_collisions(&outputs).unwrap_err().to_string();
            assert!(error.contains("are the same file"), "{error}");
        }
    }

    #[cfg(unix)]
    #[test]
    fn test_collision_existing_non_regular_files_are_skipped() {
        assert_eq!(
            collision(&[
                ("/dev/null", OutputOrigin::Metrics("--metrics-per-sample")),
                ("/dev/null", OutputOrigin::Metrics("--metrics-summary")),
            ]),
            None
        );
    }

    fn over_inputs(outputs: &[PathBuf], inputs: &[PathBuf]) -> Option<String> {
        let outputs: Vec<(PathBuf, OutputOrigin)> = outputs
            .iter()
            .map(|path| (path.clone(), OutputOrigin::Records("--unassigned".into())))
            .collect();
        let inputs: Vec<(PathBuf, InputOrigin)> = inputs
            .iter()
            .enumerate()
            .map(|(index, path)| (path.clone(), InputOrigin::Reads(index)))
            .collect();
        check_outputs_over_inputs(&outputs, &inputs)
            .err()
            .map(|e| e.to_string())
    }

    #[test]
    fn test_over_inputs_reports_which_file_the_run_reads() {
        let dir = tempfile::tempdir().unwrap();
        for (name, origin, shown) in [
            (
                "bc.tsv",
                InputOrigin::TagFile("bc".into()),
                "the tag file of group bc",
            ),
            ("sheet.tsv", InputOrigin::SampleSheet, "the sample sheet"),
            ("ref.fa", InputOrigin::Reference, "the reference"),
            (
                "ref.fa.fai",
                InputOrigin::ReferenceIndex,
                "the reference index",
            ),
        ] {
            let path = dir.path().join(name);
            std::fs::write(&path, "").unwrap();
            let outputs = [(path.clone(), OutputOrigin::Metrics("--metrics-summary"))];
            let error = check_outputs_over_inputs(&outputs, &[(path.clone(), origin)]).unwrap_err();
            let path = path.display();
            assert_eq!(
                error.to_string(),
                format!(
                    "output \"{path}\" (--metrics-summary) is {shown} ({path}); give it a \
                     distinct path"
                )
            );
        }
    }

    #[test]
    fn test_over_inputs_any_spelling_of_an_input_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.fq");
        std::fs::write(&input, "").unwrap();
        for output in ["in.fq", "./in.fq", "missing/../in.fq"] {
            let output = dir.path().join(output);
            assert_eq!(
                over_inputs(std::slice::from_ref(&output), std::slice::from_ref(&input)),
                Some(format!(
                    "output \"{}\" (--unassigned) is input 0 ({}); give it a distinct path",
                    output.display(),
                    input.display()
                ))
            );
        }
    }

    #[cfg(unix)]
    #[test]
    fn test_over_inputs_links_to_an_input_are_refused() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.fq");
        std::fs::write(&input, "").unwrap();
        std::os::unix::fs::symlink("in.fq", dir.path().join("link.fq")).unwrap();
        std::os::unix::fs::symlink(".", dir.path().join("here")).unwrap();
        std::fs::hard_link(&input, dir.path().join("hard.fq")).unwrap();
        for output in ["link.fq", "here/in.fq", "hard.fq", "missing/../hard.fq"] {
            let error = over_inputs(&[dir.path().join(output)], std::slice::from_ref(&input));
            assert!(
                error.is_some_and(|e| e.contains(" is input 0 (")),
                "{output}"
            );
        }
    }

    #[cfg(unix)]
    #[test]
    fn test_over_inputs_a_case_variant_is_refused_only_where_it_is_the_input() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.fq");
        std::fs::write(&input, "").unwrap();
        let folds_case = dir.path().join("IN.FQ").exists();
        for output in ["IN.FQ", "missing/../IN.FQ"] {
            let error = over_inputs(&[dir.path().join(output)], std::slice::from_ref(&input));
            assert_eq!(error.is_some(), folds_case, "{output}: {error:?}");
        }
    }

    #[test]
    fn test_over_inputs_reports_the_0_based_index_of_the_input() {
        let dir = tempfile::tempdir().unwrap();
        let inputs = [dir.path().join("r1.fq"), dir.path().join("r2.fq")];
        for input in &inputs {
            std::fs::write(input, "").unwrap();
        }
        let outputs = [dir.path().join("out.fq"), inputs[1].clone()];
        let error = over_inputs(&outputs, &inputs).unwrap();
        let expected = format!("is input 1 ({})", inputs[1].display());
        assert!(error.contains(&expected), "{error}");
    }

    #[test]
    fn test_over_inputs_outputs_beside_an_input_pass() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.fq");
        std::fs::write(&input, "").unwrap();
        let outputs = ["in.fq.gz", "in", "o/in.fq"].map(|output| dir.path().join(output));
        assert_eq!(over_inputs(&outputs, &[input]), None);
    }

    #[cfg(unix)]
    #[test]
    fn test_over_inputs_devices_are_exempt() {
        let null = [PathBuf::from("/dev/null")];
        assert_eq!(over_inputs(&null, &null), None);
    }

    #[test]
    fn test_collision_rendered_placeholders_differing_only_in_case_are_refused() {
        let out = pattern(
            "out/%pool/%sample.%sub_sample.R%ordinal.fq",
            &[
                Placeholder::Pool,
                Placeholder::Sample,
                Placeholder::SubSample,
                Placeholder::Ordinal,
            ],
        );
        let bin = pattern(
            "out/%pool/UN.%source.R1.fq",
            &[Placeholder::Pool, Placeholder::Source],
        );
        let render = |pool: &str, sample: &str, sub_sample: &str| {
            let ctx = PathContext {
                pool,
                sample: Some(sample),
                sub_sample: Some(sub_sample),
                ordinal: Some(1),
                source: None,
            };
            resolve_pattern(&out, &ctx).unwrap()
        };
        let unassigned = resolve_pattern(
            &bin,
            &PathContext {
                pool: "Lib",
                source: Some(0),
                ..Default::default()
            },
        )
        .unwrap();
        let check = |outputs: &[(PathBuf, OutputOrigin)]| check_output_collisions(outputs).is_err();
        let a = (render("Lib", "a", "x"), sample("a", Some("x")));
        assert!(!check(&[
            a.clone(),
            (render("Lib", "b", "x"), sample("b", Some("x")))
        ]));
        assert!(check(&[
            a.clone(),
            (render("LIB", "a", "x"), sample("a", Some("x")))
        ]));
        assert!(check(&[
            a.clone(),
            (render("Lib", "a", "X"), sample("a", Some("X")))
        ]));
        assert!(check(&[
            (render("Lib", "un", "0"), sample("un", Some("0"))),
            (
                unassigned,
                OutputOrigin::Records("--unassigned".to_string())
            ),
        ]));
    }

    #[cfg(unix)]
    #[test]
    fn test_collision_hard_links_to_one_file_are_refused() {
        let dir = tempfile::tempdir().unwrap();
        let (a, b) = (dir.path().join("a.fq"), dir.path().join("b.fq"));
        std::fs::write(&a, "").unwrap();
        std::fs::hard_link(&a, &b).unwrap();
        let outputs = [(a, sample("a", None)), (b, sample("b", None))];
        let error = check_output_collisions(&outputs).unwrap_err().to_string();
        assert!(
            error.contains("are the same file; give them distinct paths"),
            "{error}"
        );
    }

    #[cfg(unix)]
    #[test]
    fn test_collision_dangling_symbolic_links_are_followed() {
        let dir = tempfile::tempdir().unwrap();
        std::os::unix::fs::symlink("missing.fq", dir.path().join("to-missing.fq")).unwrap();
        std::os::unix::fs::symlink("new", dir.path().join("L")).unwrap();
        for (a, b) in [("missing.fq", "to-missing.fq"), ("new/a.fq", "L/a.fq")] {
            let outputs = [
                (dir.path().join(a), sample("a", None)),
                (
                    dir.path().join(b),
                    OutputOrigin::Records("--unassigned".to_string()),
                ),
            ];
            let error = check_output_collisions(&outputs).unwrap_err().to_string();
            assert!(error.contains("are the same file"), "{a} and {b}: {error}");
        }
    }

    #[test]
    fn test_collision_full_case_folding_is_refused() {
        for (a, b) in [
            ("o/\u{df}.fq", "o/ss.fq"),
            ("o/\u{17f}.fq", "o/s.fq"),
            ("o/\u{fb01}.fq", "o/fi.fq"),
            ("o/x\u{3a3}", "o/x\u{3c3}"),
        ] {
            let error = collision(&[
                (a, OutputOrigin::Records("--unassigned".to_string())),
                (b, OutputOrigin::Records("--remove grp".to_string())),
            ])
            .unwrap_or_else(|| panic!("{a} and {b} were accepted"));
            assert!(error.contains("differ only in case"), "{error}");
        }
    }

    #[test]
    fn test_collision_through_parent_directories_asks_for_distinct_paths() {
        for other in ["b", "A"] {
            let error = collision(&[
                ("o/a/../x.fq", sample("a", None)),
                (&format!("o/{other}/../x.fq"), sample(other, None)),
            ])
            .unwrap();
            assert!(
                error.ends_with("are the same file; give them distinct paths"),
                "{error}"
            );
        }
    }
}
