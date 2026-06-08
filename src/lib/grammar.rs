//! The unmux command grammar: one quote-free language shared by the CLI, config
//! files, and sheets. This module turns the raw, repeatable flag bodies of
//! [`DemuxArgs`] into a structurally validated [`DemuxPlan`] that the engine
//! consumes.
//!
//! Parsing is pure: it never touches the filesystem. Group source files, sample
//! sheets, and the records themselves are opened by later pipeline stages, so
//! the cross-checks done here are the ones that need only the grammar
//! (duplicate-key detection, contiguous input indices, reference existence,
//! attribute conflicts, and the separator rules). The semantic checks that need
//! a tag set on disk (collision, coverage, sequence membership) belong to the
//! fan-out stage.
//!
//! The grammar is built from two structural operators: `name = PRIMARY` (the
//! entity's main binding) and `name :: MORE` (a sub-part or a `key=value`
//! attribute). Distinct keys accumulate across repeated flags; setting the same
//! key twice is a fail-fast error. The full grammar is derivable from `--help`
//! (every flag carries its own examples).

use std::collections::{BTreeMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::demux::DemuxArgs;

/// A 0-based input file index (file `0` is the first input,
/// splitcode-compatible).
pub type FileIndex = usize;

/// An endpoint within a record: an offset from the 5' start, an offset counted
/// back from the 3' end, or the record's full length (the `end` keyword, valid
/// only as an end position).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Endpoint {
    /// `n`: the absolute 0-based offset `n` from the record start.
    FromStart(usize),
    /// `-n`: the position `read_len - n`, counted back from the record end.
    FromEnd(usize),
    /// `end`: the record length.
    ReadEnd,
}

impl Endpoint {
    /// Resolve to a concrete 0-based offset against a record of length `len`,
    /// clamped to `[0, len]`.
    pub fn resolve(self, len: usize) -> usize {
        match self {
            Endpoint::FromStart(n) => n.min(len),
            Endpoint::FromEnd(n) => len.saturating_sub(n),
            Endpoint::ReadEnd => len,
        }
    }
}

/// A `loc=file:start[:end]` search window. An absent `end` means "to the record
/// end".
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Location {
    /// The 0-based input file index the window applies to. There is no "any
    /// file" wildcard; omit `loc` to search every file (a negative file index
    /// is rejected, since it would collide with the from-the-end meaning a
    /// negative number carries in `start`/`end`).
    pub file: FileIndex,
    /// The window start.
    pub start: Endpoint,
    /// The window end; `None` means the record end.
    pub end: Option<Endpoint>,
}

/// An error budget for matching a tag, `dist=mismatch[:indel[:total]]`
/// (splitcode-compatible). An omitted field is `0`; when `total` is unset it
/// becomes `mismatch + indel`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Dist {
    /// Maximum substitutions.
    pub mismatch: usize,
    /// Maximum insertions plus deletions.
    pub indel: usize,
    /// Maximum total edit distance (the value passed to sassy as `k`).
    pub total: usize,
}

/// A truncated-end matching policy, `min:freq` (`partial5`/`partial3`): allow a
/// tag to be cut off at the record's 5'/3' end if at least `min_match` bases
/// align with at most `max_mismatch_freq` mismatches per matched base.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Partial {
    /// Minimum matched bases required.
    pub min_match: usize,
    /// Maximum mismatch frequency over the matched region.
    pub max_mismatch_freq: f64,
}

/// How matches within `dist` are accepted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchMode {
    /// Keep every match within `dist`.
    All,
    /// Require the best match to beat the runner-up by `delta`, else leave the
    /// record unmatched.
    Nearest,
}

/// A `next=GROUP[:lo-hi]` sequential link to a downstream group. An absent
/// window is a plain ordering prerequisite; a window `lo-hi` is relative to
/// this group's matched end and replaces the downstream group's `loc`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NextLink {
    /// The downstream group that must match after this one.
    pub group: String,
    /// An optional relative window `[matchEnd+lo, matchEnd+hi]` for the
    /// downstream group.
    pub window: Option<(usize, usize)>,
}

/// The source of a group's tag set: a file (resolved later) or an inline set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GroupSource {
    /// `NAME=PATH`: a tag-id column plus a sequence column, resolved at load
    /// time.
    File(PathBuf),
    /// `NAME={...}`: tags given inline.
    Inline(Vec<InlineTag>),
}

/// One tag of an inline `{...}` set. When no `id=` is given, splitcode-style,
/// the sequence doubles as the id.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InlineTag {
    /// The tag id.
    pub id: String,
    /// The tag sequence (DNA / IUPAC).
    pub seq: String,
}

/// `anchor=`: which edge of the `loc` window each tag is pinned to. The
/// attribute is optional; when unset (the default) each tag is searched
/// anywhere within the window (it slides).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Anchor {
    /// `5p`: each tag's 5' base is fixed at `loc.start` and matched over
    /// exactly its own length (`[loc.start, loc.start+len(tag))`).
    /// Length-strict; resolves variable-length 5'-anchored groups where a
    /// shorter tag would otherwise slide to a spurious interior offset.
    FivePrime,
    /// `3p`: each tag's 3' base is fixed at `loc.end` and matched over exactly
    /// its own length (`[loc.end-len(tag), loc.end)`). The 3' mirror of
    /// [`Anchor::FivePrime`]; length-strict and (for now) substitution-only.
    ThreePrime,
}

/// The optional `--group NAME::ATTRS`; every field defaults per the grammar's
/// defaults table.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GroupAttrs {
    /// `loc=`: the search window (default: whole record).
    pub loc: Option<Location>,
    /// `dist=`: the error budget (default: exact).
    pub dist: Option<Dist>,
    /// `delta=`: the min best-vs-runner-up gap (only meaningful with
    /// `mode=nearest`).
    pub delta: Option<usize>,
    /// `mode=`: match acceptance (default `all`).
    pub mode: Option<MatchMode>,
    /// `next=`: a downstream sequential link.
    pub next: Option<NextLink>,
    /// `prev=`: an upstream sequential prerequisite.
    pub prev: Option<String>,
    /// `minFindsPerTag=`.
    pub min_finds_per_tag: Option<usize>,
    /// `maxFindsPerTag=`.
    pub max_finds_per_tag: Option<usize>,
    /// `minFindsPerGroup=`.
    pub min_finds_per_group: Option<usize>,
    /// `maxFindsPerGroup=`.
    pub max_finds_per_group: Option<usize>,
    /// `both_strands=` (default `false`).
    pub revcomp: Option<bool>,
    /// `partial5=min:freq`.
    pub partial5: Option<Partial>,
    /// `partial3=min:freq`.
    pub partial3: Option<Partial>,
    /// `match=stream[+stream]`: match the joined `--extract` streams instead of
    /// a record window.
    pub match_streams: Option<Vec<String>>,
    /// `anchor=`: pin tags to a `loc` edge (`5p`/`3p`); unset means the default
    /// sliding window.
    pub anchor: Option<Anchor>,
}

/// A fully assembled tag group: its name, tag set, and attributes.
#[derive(Debug, Clone, PartialEq)]
pub struct Group {
    /// The group name (referenced by `@grp`, selectors, `next`/`prev`).
    pub name: String,
    /// Where the tag set comes from.
    pub source: GroupSource,
    /// The group's attributes.
    pub attrs: GroupAttrs,
}

/// The length of an anchored extraction span past a match end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnchorLen {
    /// A fixed number of bases.
    Bases(usize),
    /// To the record end (`:end`).
    ToEnd,
}

/// The body of an `--extract` span.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpanBody {
    /// `file:start:end` (the trailing number is an end position).
    File {
        /// Input file index.
        file: FileIndex,
        /// Start endpoint.
        start: Endpoint,
        /// End endpoint.
        end: Endpoint,
    },
    /// `@group`: the group's matched span itself (resolves variable-length
    /// barcodes).
    AnchorMatch {
        /// The anchoring group.
        group: String,
    },
    /// `@group+offset:len` (rightward) or `@group-offset:len` (leftward): `len`
    /// bases offset from the group's match. With `+` (`before` = false) the
    /// span starts `offset` bases past the matched END. With `-` (`before` =
    /// true) the span ENDS `offset` bases before the matched START, e.g.
    /// `@grp-0:9` is the 9 bp immediately left of the match (a UMI upstream of
    /// a barcode). The trailing number is a length, not an end; both directions
    /// clamp at the record edge (`:end` reaches the record end for `+`, or the
    /// record start for `-`).
    AnchorOffset {
        /// The anchoring group.
        group: String,
        /// Whether the span lies before the matched start (`-`) rather than
        /// past the matched end (`+`).
        before: bool,
        /// Bases past the match end (`+`) or before the match start (`-`).
        offset: usize,
        /// The span length.
        len: AnchorLen,
    },
    /// `@grpA..@grpB`: the region between two matched anchors.
    Between {
        /// The upstream anchor (region starts at its match end).
        from: String,
        /// The downstream anchor (region ends at its match start).
        to: String,
    },
}

/// An `--extract NAME=SPAN`: a named stream carrying bases and qualities.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Extract {
    /// The stream name (referenced by `--template`/`--tag`).
    pub name: String,
    /// The extraction span.
    pub span: Span,
}

/// A reference to an extract stream in a `--tag` or `--template`. A leading `~`
/// reverse-complements this stream's contribution at this use site (rc the
/// bases, reverse the quals); the stream itself is defined orientation-neutral
/// by its `--extract`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StreamRef {
    /// The extract stream name.
    pub name: String,
    /// Whether a leading `~` reverse-complements this stream's contribution
    /// here.
    pub revcomp: bool,
}

impl std::fmt::Display for StreamRef {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.revcomp {
            write!(f, "~{}", self.name)
        } else {
            f.write_str(&self.name)
        }
    }
}

/// Render a stream list joined by `sep`, each as written (a `~` prefix on a
/// reverse-complemented one).
fn join_streams(streams: &[StreamRef], sep: &str) -> String {
    streams
        .iter()
        .map(StreamRef::to_string)
        .collect::<Vec<_>>()
        .join(sep)
}

/// An extraction span.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Span {
    /// The span body.
    pub body: SpanBody,
}

/// A `--template`: an ordered list of streams assembled into one output record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Template {
    /// The streams to concatenate, in order (each may carry a `~`
    /// reverse-complement).
    pub streams: Vec<StreamRef>,
    /// `raw=true`: emit the observed bases for this body even when a stream has
    /// a corrected form (default `false`: a stream from a matched keeplist
    /// `@grp` is emitted corrected, others observed). Mirrors the `--tag`
    /// `raw=` attribute.
    pub raw: bool,
}

/// The paired quality tag of a SAM sequence tag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QualTag {
    /// Use the built-in seq->qual map (`CB`/`CY`, `CR`/`CY`, `RX`/`QX`,
    /// `BC`/`QT`, `OX`/`BZ`); no quality tag for an unmapped sequence tag.
    Default,
    /// `qual=QTAG`: an explicit quality tag.
    Named(String),
    /// `qual=none`: drop the quality tag.
    None,
}

/// A `--tag` binding: one SAM sequence tag assembled from one or more streams,
/// plus its paired quality tag and join separators.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TagBinding {
    /// The two-character SAM sequence tag (e.g. `CB`).
    pub tag: String,
    /// The streams concatenated into the tag value, in order (each may carry a
    /// `~` reverse-complement).
    pub streams: Vec<StreamRef>,
    /// The paired quality tag policy.
    pub qual: QualTag,
    /// `sep=`: the sequence separator joining two or more streams (defaults to
    /// `-`); `None` for a single-stream tag, which joins nothing.
    pub sep: Option<String>,
    /// `qual-sep=`: the quality separator joining the per-stream qualities of a
    /// multi-stream tag that emits a quality tag (defaults to a single space);
    /// `None` otherwise.
    pub qual_sep: Option<String>,
    /// `raw=true`: emit the observed (as-sequenced) bases for this tag rather
    /// than the corrected form. Default `false`: a tag emits the
    /// error-corrected barcode when its stream came from a matched keeplist
    /// `@grp`. No effect on a stream that was never corrected.
    pub raw: bool,
}

/// An optional finer fan-out level under a sample, written to read-group `LB`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SubSample {
    /// A literal sub_sample id.
    Literal(String),
    /// `%pool`: the pool id.
    Pool,
}

/// A selector over one group: an OR pool of tag tokens, or (when `members` is
/// empty) the whole group (any of its tags).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroupSelector {
    /// The referenced group.
    pub group: String,
    /// Tag tokens (ids or sequences) forming an OR pool; empty means "the whole
    /// group".
    pub members: Vec<String>,
}

/// A `--sample` selector: AND-joined group selectors (joined with `+`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Selector {
    /// The AND terms; a record must satisfy every term.
    pub terms: Vec<GroupSelector>,
}

/// A `--sample SAMPLE[::SUB_SAMPLE]=SELECTOR` fan-out target.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sample {
    /// The sample id (written to read-group `SM`).
    pub sample: String,
    /// The optional sub_sample (written to read-group `LB`).
    pub sub_sample: Option<SubSample>,
    /// The selector that routes records to this target.
    pub selector: Selector,
}

/// Where the fan-out targets come from. These three flags are mutually
/// exclusive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SampleSpec {
    /// No demux: the whole pool is one passthrough target.
    None,
    /// `--sample`: inline targets.
    Inline(Vec<Sample>),
    /// `--sample-sheet`: a table of targets, resolved at load time.
    Sheet(PathBuf),
    /// `--sample-from-group`: 1:1 tag-to-sample over the named group.
    FromGroup(String),
}

/// The selector half of a `--remove SEL[=PATTERN]` rule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RemoveSelector {
    /// A whole group: any record matching any of its tags is removed.
    Group(String),
    /// A specific `group::id` combination.
    GroupTag {
        /// The group.
        group: String,
        /// The tag id.
        id: String,
    },
}

/// A `--remove SEL[=PATTERN]` rule: matching records are removed (tallied as
/// `removed`) and, when a pattern is given, also written there.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemoveRule {
    /// What to remove.
    pub selector: RemoveSelector,
    /// An optional output pattern for the removed records (raw input segments).
    pub pattern: Option<OutputPattern>,
}

/// An output-path placeholder.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Placeholder {
    /// `%pool`: the pool id.
    Pool,
    /// `%sample`: one file per sample.
    Sample,
    /// `%sub_sample`: one file per sub_sample (empty when a sample has none).
    SubSample,
    /// `%ordinal`: the 1-based position in the `--template` list.
    Ordinal,
    /// `%source`: the 0-based input file index (for `--unassigned`/`--remove`).
    Source,
}

impl Placeholder {
    /// The placeholder's name as it appears after `%`.
    pub fn name(self) -> &'static str {
        match self {
            Placeholder::Pool => "pool",
            Placeholder::Sample => "sample",
            Placeholder::SubSample => "sub_sample",
            Placeholder::Ordinal => "ordinal",
            Placeholder::Source => "source",
        }
    }
}

/// One piece of a parsed output pattern.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PatternSegment {
    /// A literal path fragment (percent escapes already decoded).
    Literal(String),
    /// A placeholder to expand at write time.
    Placeholder(Placeholder),
}

/// A parsed `--out`/`--unassigned`/`--remove` output pattern: literal fragments
/// interleaved with placeholders, validated against the placeholders allowed in
/// that flag's context.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputPattern {
    /// The ordered segments.
    pub segments: Vec<PatternSegment>,
}

impl OutputPattern {
    /// Whether the pattern uses the given placeholder.
    pub fn uses(&self, placeholder: Placeholder) -> bool {
        self.segments
            .iter()
            .any(|s| matches!(s, PatternSegment::Placeholder(p) if *p == placeholder))
    }
}

/// A structurally validated demux plan: the parsed form of [`DemuxArgs`] the
/// engine consumes.
#[derive(Debug, Clone, PartialEq)]
pub struct DemuxPlan {
    /// `--pool`: the pool id (fills `%pool`; provenance). `None` defers to the
    /// input-stem default.
    pub pool: Option<String>,
    /// Input record files, indexed by their 0-based file index.
    pub inputs: Vec<PathBuf>,
    /// Tag groups, in declaration order.
    pub groups: Vec<Group>,
    /// Extraction streams, in declaration order.
    pub extracts: Vec<Extract>,
    /// Output record templates (one per output record).
    pub templates: Vec<Template>,
    /// SAM tag bindings.
    pub tags: Vec<TagBinding>,
    /// Shared `@RG` header fields (`SM`/`LB` excluded).
    pub rg_tags: Vec<(String, String)>,
    /// The fan-out source.
    pub samples: SampleSpec,
    /// `--require-samples-explain-all-tags`: require every selector-referenced
    /// group's tags to be claimed by some sample (off by default).
    pub require_samples_explain_all_tags: bool,
    /// `--qc-tag`: the SAM tag for the optional per-record demux-provenance
    /// JSON slug; `None` = off.
    pub qc_tag: Option<String>,
    /// Record-routing skiplist rules.
    pub remove: Vec<RemoveRule>,
    /// `--out`: the output pattern (`None` = stdout in the input format).
    pub out: Option<OutputPattern>,
    /// `--unassigned`: the no-sample-match output pattern.
    pub unassigned: Option<OutputPattern>,
    /// `--metrics-per-sample` TSV path pattern. Pool-level file, so only
    /// `%pool` is valid.
    pub metrics_per_sample: Option<OutputPattern>,
    /// `--metrics-summary` TSV path pattern. Pool-level file, so only `%pool`
    /// is valid.
    pub metrics_summary: Option<OutputPattern>,
    /// `--per-record`: disable auto pair-detection.
    pub per_record: bool,
    /// `--compression` level (0-9).
    pub compression: u8,
    /// `--threads` worker count.
    pub threads: usize,
}

/// Parse and structurally validate the raw flag surface into a [`DemuxPlan`].
pub fn parse_demux(args: &DemuxArgs) -> Result<DemuxPlan> {
    let inputs = parse_inputs(&args.inputs)?;
    let groups = parse_groups(&args.groups)?;
    let extracts = parse_extracts(&args.extracts)?;
    let templates = parse_templates(&args.templates)?;
    let tags = parse_tags(&args.tags)?;
    let rg_tags = parse_rg_tags(&args.rg_tags)?;
    let samples = parse_sample_spec(args)?;
    let remove = parse_removes(&args.remove)?;
    let out = match &args.out {
        Some(p) => Some(parse_output_pattern(
            p,
            &[
                Placeholder::Pool,
                Placeholder::Sample,
                Placeholder::SubSample,
                Placeholder::Ordinal,
            ],
        )?),
        None => None,
    };
    let unassigned = match &args.unassigned {
        Some(p) => Some(parse_output_pattern(
            p,
            &[Placeholder::Pool, Placeholder::Source],
        )?),
        None => None,
    };

    let plan = DemuxPlan {
        pool: args.pool.clone(),
        inputs,
        groups,
        extracts,
        templates,
        tags,
        rg_tags,
        samples,
        require_samples_explain_all_tags: args.require_samples_explain_all_tags,
        qc_tag: match &args.qc_tag {
            Some(tag) => Some(validate_qc_tag(tag)?),
            None => None,
        },
        remove,
        out,
        unassigned,
        metrics_per_sample: parse_metrics_pattern(
            args.metrics_per_sample.as_deref(),
            "--metrics-per-sample",
        )?,
        metrics_summary: parse_metrics_pattern(
            args.metrics_summary.as_deref(),
            "--metrics-summary",
        )?,
        per_record: args.per_record,
        compression: args.compression,
        threads: args.threads,
    };
    validate_plan(&plan)?;
    for warning in plan_warnings(&plan) {
        log::warn!("{warning}");
    }
    Ok(plan)
}

/// Parse `N=PATH` entries into a path list indexed by file index, requiring a
/// contiguous `0..N`.
fn parse_inputs(specs: &[String]) -> Result<Vec<PathBuf>> {
    if specs.is_empty() {
        bail!("no input reads given; provide read files positionally or with `--in N=PATH`");
    }
    let mut by_index: BTreeMap<usize, PathBuf> = BTreeMap::new();
    for spec in specs {
        let (n, path) = spec
            .split_once('=')
            .with_context(|| format!("invalid `--in` entry {spec:?}; expected `N=PATH`"))?;
        let index: usize = n
            .trim()
            .parse()
            .with_context(|| format!("invalid input index {:?} in {spec:?}", n.trim()))?;
        if by_index.insert(index, PathBuf::from(path)).is_some() {
            bail!("input index {index} is assigned more than once");
        }
    }
    let got: Vec<usize> = by_index.keys().copied().collect();
    let expected: Vec<usize> = (0..by_index.len()).collect();
    if got != expected {
        bail!(
            "input indices must be contiguous from 0 (a permutation of 0..{}); got {got:?}",
            by_index.len()
        );
    }
    // A single stdin stream cannot feed two inputs, so at most one `-` is
    // allowed.
    if by_index
        .values()
        .filter(|path| path.as_os_str() == "-")
        .count()
        > 1
    {
        bail!("stdin (`-`) may be used for at most one input");
    }
    Ok(by_index.into_values().collect())
}

/// Whether a `--group`/`--tag` body binds a primary (`=`) or attributes (`::`).
enum Binding<'a> {
    Primary(&'a str),
    Attrs(&'a str),
}

/// Split a `name=...` / `name::...` body into its name and binding kind.
fn split_binding(spec: &str) -> Result<(&str, Binding<'_>)> {
    let name_end = spec.find(|c: char| !is_ident_char(c)).unwrap_or(spec.len());
    let (name, rest) = spec.split_at(name_end);
    if name.is_empty() {
        bail!("missing name in {spec:?}");
    }
    if let Some(attrs) = rest.strip_prefix("::") {
        Ok((name, Binding::Attrs(attrs)))
    } else if let Some(primary) = rest.strip_prefix('=') {
        Ok((name, Binding::Primary(primary)))
    } else {
        bail!("expected `=` or `::` after `{name}` in {spec:?}");
    }
}

/// Builder accumulating a group's source and attributes across repeated flags.
struct GroupBuilder {
    name: String,
    source: Option<GroupSource>,
    attrs: GroupAttrs,
}

/// Fold the `--group` flag bodies into ordered, fully assembled groups.
fn parse_groups(specs: &[String]) -> Result<Vec<Group>> {
    let mut builders: Vec<GroupBuilder> = Vec::new();
    for spec in specs {
        let (name, binding) = split_binding(spec)?;
        validate_ident(name, "group name")?;
        let builder = match builders.iter_mut().find(|b| b.name == name) {
            Some(builder) => builder,
            None => {
                builders.push(GroupBuilder {
                    name: name.to_string(),
                    source: None,
                    attrs: GroupAttrs::default(),
                });
                builders.last_mut().unwrap()
            }
        };
        match binding {
            Binding::Primary(value) => {
                if builder.source.is_some() {
                    bail!("group `{name}` source is defined more than once");
                }
                builder.source = Some(parse_group_source(value, name)?);
            }
            Binding::Attrs(attrs) => {
                apply_group_attrs(&mut builder.attrs, attrs, name)?;
            }
        }
    }
    builders
        .into_iter()
        .map(|b| {
            let source = b.source.with_context(|| {
                format!(
                    "group `{}` has attributes but no tag set (a `{}=SOURCE` is required)",
                    b.name, b.name
                )
            })?;
            Ok(Group {
                name: b.name,
                source,
                attrs: b.attrs,
            })
        })
        .collect()
}

/// Parse a group source: an inline `{...}` set or a file path.
fn parse_group_source(value: &str, group: &str) -> Result<GroupSource> {
    let value = value.trim();
    if let Some(rest) = value.strip_prefix('{') {
        let inner = rest.strip_suffix('}').with_context(|| {
            format!("inline tag set for group `{group}` is missing a closing `}}`: {value:?}")
        })?;
        let mut tags = Vec::new();
        for entry in inner.split(',') {
            let entry = entry.trim();
            if entry.is_empty() {
                continue;
            }
            let tag = match entry.split_once('=') {
                Some((id, seq)) => InlineTag {
                    id: id.trim().to_string(),
                    seq: seq.trim().to_string(),
                },
                // No id given: the sequence doubles as the id (splitcode
                // behavior).
                None => InlineTag {
                    id: entry.to_string(),
                    seq: entry.to_string(),
                },
            };
            tags.push(tag);
        }
        if tags.is_empty() {
            bail!("inline tag set for group `{group}` is empty");
        }
        Ok(GroupSource::Inline(tags))
    } else {
        if value.is_empty() {
            bail!("group `{group}` has an empty source");
        }
        Ok(GroupSource::File(PathBuf::from(value)))
    }
}

/// Apply a comma-list of `key=value` attributes onto a group, erroring on any
/// duplicate key.
fn apply_group_attrs(attrs: &mut GroupAttrs, body: &str, group: &str) -> Result<()> {
    for field in body.split(',') {
        let field = field.trim();
        if field.is_empty() {
            continue;
        }
        // Most attributes are `key=value`; `findOne` is a valueless boolean
        // shorthand (and also accepts `findOne=true|false`).
        let (key, value) = match field.split_once('=') {
            Some((key, value)) => (key.trim(), value.trim()),
            None if field == "findOne" => ("findOne", "true"),
            None => bail!("group `{group}` attribute {field:?} is not `key=value`"),
        };
        macro_rules! set_once {
            ($slot:expr, $parsed:expr) => {{
                if $slot.is_some() {
                    bail!("group `{group}` attribute `{key}` is set more than once");
                }
                $slot = Some($parsed);
            }};
        }
        match key {
            "loc" => set_once!(attrs.loc, parse_location(value)?),
            "dist" => set_once!(attrs.dist, parse_dist(value)?),
            "delta" => set_once!(attrs.delta, parse_usize(value, "delta")?),
            "mode" => set_once!(attrs.mode, parse_mode(value)?),
            "next" => set_once!(attrs.next, parse_next(value)?),
            "prev" => set_once!(attrs.prev, parse_group_ref(value, "prev")?),
            "minFindsPerTag" => {
                set_once!(
                    attrs.min_finds_per_tag,
                    parse_usize(value, "minFindsPerTag")?
                )
            }
            "maxFindsPerTag" => {
                set_once!(
                    attrs.max_finds_per_tag,
                    parse_usize(value, "maxFindsPerTag")?
                )
            }
            "minFindsPerGroup" => set_once!(
                attrs.min_finds_per_group,
                parse_usize(value, "minFindsPerGroup")?
            ),
            "maxFindsPerGroup" => set_once!(
                attrs.max_finds_per_group,
                parse_usize(value, "maxFindsPerGroup")?
            ),
            // `findOne`: shorthand for "this group resolves to exactly one
            // match" -> required (minFindsPerGroup=1), at most one total
            // (maxFindsPerGroup=1, which already implies no tag twice, so
            // maxFindsPerTag=1 is set too for explicitness). Goes through
            // `set_once!`, so it conflicts with an explicit min/max finds cap
            // rather than silently overriding it.
            "findOne" => {
                if parse_bool(value, "findOne")? {
                    set_once!(attrs.max_finds_per_tag, 1);
                    set_once!(attrs.min_finds_per_group, 1);
                    set_once!(attrs.max_finds_per_group, 1);
                }
            }
            "both_strands" => set_once!(attrs.revcomp, parse_bool(value, "both_strands")?),
            "partial5" => set_once!(attrs.partial5, parse_partial(value, "partial5")?),
            "partial3" => set_once!(attrs.partial3, parse_partial(value, "partial3")?),
            "match" => set_once!(attrs.match_streams, parse_match_streams(value)?),
            "anchor" => set_once!(attrs.anchor, parse_anchor(value)?),
            other => bail!("group `{group}` has unknown attribute `{other}`"),
        }
    }
    Ok(())
}

/// Parse `mismatch[:indel[:total]]` with splitcode's defaults and validation.
fn parse_dist(value: &str) -> Result<Dist> {
    let fields: Vec<&str> = value.split(':').collect();
    if fields.len() > 3 {
        bail!("dist must be mismatch[:indel[:total]], got {value:?}");
    }
    let mut nums = [0usize; 3];
    for (i, field) in fields.iter().enumerate() {
        if !field.is_empty() {
            nums[i] = field
                .parse()
                .with_context(|| format!("invalid dist field {field:?} in {value:?}"))?;
        }
    }
    let (mismatch, indel, mut total) = (nums[0], nums[1], nums[2]);
    if total != 0 && (mismatch + indel < total || mismatch > total || indel > total) {
        bail!("dist is inconsistent: {value:?} (total must be >= mismatch, >= indel, and <= mismatch+indel)");
    }
    if total == 0 {
        total = mismatch + indel;
    }
    Ok(Dist {
        mismatch,
        indel,
        total,
    })
}

/// Parse a `loc=file:start[:end]` window.
fn parse_location(value: &str) -> Result<Location> {
    let parts: Vec<&str> = value.split(':').collect();
    if parts.len() < 2 || parts.len() > 3 {
        bail!("loc must be file:start[:end], got {value:?}");
    }
    let file: FileIndex = parts[0].parse().with_context(|| {
        format!(
            "invalid loc file index {:?} in {value:?} (must be a non-negative integer; omit loc to search every file)",
            parts[0]
        )
    })?;
    let start = parse_endpoint(parts[1], false, "loc start")?;
    let end = match parts.get(2) {
        Some(e) => Some(parse_endpoint(e, true, "loc end")?),
        None => None,
    };
    Ok(Location { file, start, end })
}

/// Parse a record endpoint: `end`, a `-n` from-end offset, or an absolute
/// offset.
fn parse_endpoint(token: &str, allow_read_end: bool, ctx: &str) -> Result<Endpoint> {
    if token == "end" {
        if !allow_read_end {
            bail!("`end` is not valid as the {ctx}");
        }
        return Ok(Endpoint::ReadEnd);
    }
    if let Some(rest) = token.strip_prefix('-') {
        let n = rest
            .parse()
            .with_context(|| format!("invalid {ctx} {token:?}"))?;
        return Ok(Endpoint::FromEnd(n));
    }
    let n = token
        .parse()
        .with_context(|| format!("invalid {ctx} {token:?}"))?;
    Ok(Endpoint::FromStart(n))
}

/// Parse `mode=all|nearest`.
fn parse_mode(value: &str) -> Result<MatchMode> {
    match value {
        "all" => Ok(MatchMode::All),
        "nearest" => Ok(MatchMode::Nearest),
        other => bail!("mode must be `all` or `nearest`, got `{other}`"),
    }
}

/// Parse `anchor=5p|3p`. The attribute is omitted for the default (unanchored,
/// sliding) window.
fn parse_anchor(value: &str) -> Result<Anchor> {
    match value {
        "5p" => Ok(Anchor::FivePrime),
        "3p" => Ok(Anchor::ThreePrime),
        other => bail!("anchor must be `5p` or `3p`, got `{other}`"),
    }
}

/// Parse `next=GROUP[:lo-hi]`.
fn parse_next(value: &str) -> Result<NextLink> {
    let (group, window) = match value.split_once(':') {
        Some((g, win)) => {
            let (lo, hi) = win
                .split_once('-')
                .with_context(|| format!("next window must be `lo-hi`, got {win:?}"))?;
            // Bounds are non-negative offsets past the upstream match end; a
            // leading `-` (a negative bound) splits to an empty low, so report
            // the whole window rather than the empty piece.
            let bounds_msg =
                || format!("next window `{win}` bounds must be non-negative integers (`lo-hi`)");
            let lo = lo.trim().parse().with_context(bounds_msg)?;
            let hi = hi.trim().parse().with_context(bounds_msg)?;
            if hi < lo {
                bail!("next window high {hi} is below low {lo}");
            }
            (g, Some((lo, hi)))
        }
        None => (value, None),
    };
    let group = group.trim();
    validate_ident(group, "next group")?;
    Ok(NextLink {
        group: group.to_string(),
        window,
    })
}

/// Parse `partial5`/`partial3` as `min:freq`.
fn parse_partial(value: &str, attr: &str) -> Result<Partial> {
    let (min, freq) = value
        .split_once(':')
        .with_context(|| format!("{attr} must be `min:freq`, got {value:?}"))?;
    let min_match = min
        .trim()
        .parse()
        .with_context(|| format!("invalid {attr} min `{min}`"))?;
    let max_mismatch_freq: f64 = freq
        .trim()
        .parse()
        .with_context(|| format!("invalid {attr} freq `{freq}`"))?;
    if !(0.0..=1.0).contains(&max_mismatch_freq) {
        bail!("{attr} freq must be in [0, 1], got {max_mismatch_freq}");
    }
    Ok(Partial {
        min_match,
        max_mismatch_freq,
    })
}

/// Parse `match=stream[+stream]` into the ordered stream list.
fn parse_match_streams(value: &str) -> Result<Vec<String>> {
    let streams: Vec<String> = value
        .split('+')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if streams.is_empty() {
        bail!("match must list one or more streams, got {value:?}");
    }
    for stream in &streams {
        validate_stream_name(stream)?;
    }
    Ok(streams)
}

/// Parse a bare group reference (`prev=`), validating its identifier.
fn parse_group_ref(value: &str, attr: &str) -> Result<String> {
    let value = value.trim();
    validate_ident(value, &format!("{attr} group"))?;
    Ok(value.to_string())
}

/// Parse the `--extract` flag bodies.
fn parse_extracts(specs: &[String]) -> Result<Vec<Extract>> {
    let mut extracts = Vec::new();
    let mut seen = HashSet::new();
    for spec in specs {
        let extract = parse_extract(spec)?;
        if !seen.insert(extract.name.clone()) {
            bail!(
                "extract stream `{}` is defined more than once",
                extract.name
            );
        }
        extracts.push(extract);
    }
    Ok(extracts)
}

/// Parse one `--extract NAME=SPAN`.
fn parse_extract(spec: &str) -> Result<Extract> {
    let (name, span) = spec
        .split_once('=')
        .with_context(|| format!("invalid `--extract` {spec:?}; expected `NAME=SPAN`"))?;
    let name = name.trim();
    validate_stream_name(name)?;
    let span = parse_span(span.trim())?;
    Ok(Extract {
        name: name.to_string(),
        span,
    })
}

/// Parse an extraction span. Reverse-complement is no longer a span property: a
/// leading `~` here is an error, pointing at the use-site form (`--tag
/// T=~stream` / `--template ~stream`).
pub fn parse_span(value: &str) -> Result<Span> {
    if let Some(rest) = value.strip_prefix('~') {
        bail!(
            "`~` is not allowed on an `--extract` span ({value:?}); reverse-complement at the point of \
             use instead, e.g. `--tag T=~{rest}`-style on the stream name in a `--tag`/`--template`"
        );
    }
    let body = if let Some(rest) = value.strip_prefix('@') {
        parse_anchor_body(rest)?
    } else {
        parse_file_span(value)?
    };
    Ok(Span { body })
}

/// Parse the body after a leading `@` in an extraction span.
fn parse_anchor_body(rest: &str) -> Result<SpanBody> {
    // `@grpA..@grpB`: between two anchors.
    if let Some((from, to)) = rest.split_once("..") {
        let from = from.trim();
        let to = to
            .trim()
            .strip_prefix('@')
            .with_context(|| format!("between-anchor span needs `@grpA..@grpB`, got `@{rest}`"))?
            .trim();
        validate_ident(from, "anchor group")?;
        validate_ident(to, "anchor group")?;
        return Ok(SpanBody::Between {
            from: from.to_string(),
            to: to.to_string(),
        });
    }
    // `@group+offset:len` (past the match end) or `@group-offset:len` (before
    // the match start). Group idents never contain `+`/`-`, so the first such
    // separator splits the group from `offset:len`.
    for (sep, before) in [('+', false), ('-', true)] {
        if let Some((group, tail)) = rest.split_once(sep) {
            let group = group.trim();
            validate_ident(group, "anchor group")?;
            let (offset, len) = tail.split_once(':').with_context(|| {
                format!("anchored span needs `@group{sep}offset:len`, got `@{rest}`")
            })?;
            let offset = offset
                .trim()
                .parse()
                .with_context(|| format!("invalid anchor offset {offset:?}"))?;
            let len = if len.trim() == "end" {
                AnchorLen::ToEnd
            } else {
                AnchorLen::Bases(
                    len.trim()
                        .parse()
                        .with_context(|| format!("invalid anchor length {len:?}"))?,
                )
            };
            return Ok(SpanBody::AnchorOffset {
                group: group.to_string(),
                before,
                offset,
                len,
            });
        }
    }
    // `@group`: the group's matched span itself.
    let group = rest.trim();
    validate_ident(group, "anchor group")?;
    Ok(SpanBody::AnchorMatch {
        group: group.to_string(),
    })
}

/// Parse a `file:start:end` span (the trailing number is an end position).
fn parse_file_span(value: &str) -> Result<SpanBody> {
    let parts: Vec<&str> = value.split(':').collect();
    if parts.len() != 3 {
        bail!("file span must be `file:start:end`, got {value:?}");
    }
    let file = parts[0]
        .parse()
        .with_context(|| format!("invalid span file index {:?} in {value:?}", parts[0]))?;
    let start = parse_endpoint(parts[1], false, "span start")?;
    let end = parse_endpoint(parts[2], true, "span end")?;
    Ok(SpanBody::File { file, start, end })
}

/// Parse the `--template` flag bodies.
fn parse_templates(specs: &[String]) -> Result<Vec<Template>> {
    specs.iter().map(|s| parse_template(s)).collect()
}

/// Parse one `--template stream[+stream][::raw=true|false]`.
fn parse_template(spec: &str) -> Result<Template> {
    // `::ATTRS` (`raw=`) is optional and trails the stream list; stream names
    // cannot contain `:`, so the split is unambiguous.
    let (stream_part, attr_part) = match spec.split_once("::") {
        Some((streams, attrs)) => (streams, Some(attrs)),
        None => (spec, None),
    };
    if stream_part.contains(',') {
        bail!(
            "`--template` joins streams with `+` (e.g. `r1+r2`), not `,`; for separate output reads \
             use `--template r1 r2` or repeat `--template`. Got: {spec:?}"
        );
    }
    let streams: Vec<StreamRef> = stream_part
        .split('+')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(parse_stream_ref)
        .collect();
    if streams.is_empty() {
        bail!("`--template` is empty: {spec:?}");
    }
    for stream in &streams {
        validate_stream_name(&stream.name)?;
    }
    let mut raw = false;
    if let Some(attrs) = attr_part {
        for field in attrs.split(',') {
            let field = field.trim();
            if field.is_empty() {
                continue;
            }
            let (key, value) = field
                .split_once('=')
                .with_context(|| format!("`--template` attribute {field:?} is not `key=value`"))?;
            match key.trim() {
                "raw" => raw = parse_bool(value.trim(), "raw")?,
                other => {
                    bail!("`--template` has unknown attribute `{other}` (only `raw=` is valid)")
                }
            }
        }
    }
    Ok(Template { streams, raw })
}

/// Builder accumulating a tag's streams and attributes across repeated flags.
struct TagBuilder {
    tag: String,
    streams: Option<Vec<StreamRef>>,
    qual: Option<QualTag>,
    sep: Option<String>,
    qual_sep: Option<String>,
    raw: Option<bool>,
}

/// Fold the `--tag` flag bodies into assembled, validated tag bindings.
fn parse_tags(specs: &[String]) -> Result<Vec<TagBinding>> {
    let mut builders: Vec<TagBuilder> = Vec::new();
    for spec in specs {
        let (tag, binding) = split_binding(spec)?;
        validate_sam_tag(tag)?;
        let builder = match builders.iter_mut().find(|b| b.tag == tag) {
            Some(builder) => builder,
            None => {
                builders.push(TagBuilder {
                    tag: tag.to_string(),
                    streams: None,
                    qual: None,
                    sep: None,
                    qual_sep: None,
                    raw: None,
                });
                builders.last_mut().unwrap()
            }
        };
        match binding {
            Binding::Primary(value) => {
                if builder.streams.is_some() {
                    bail!("tag `{tag}` streams are defined more than once");
                }
                if value.contains(',') {
                    bail!(
                        "tag `{tag}` joins streams with `+` (e.g. `{tag}=a+b+c`), not `,`. Got: {value:?}"
                    );
                }
                let streams: Vec<StreamRef> = value
                    .split('+')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(parse_stream_ref)
                    .collect();
                if streams.is_empty() {
                    bail!("tag `{tag}` has no streams: {value:?}");
                }
                for stream in &streams {
                    validate_stream_name(&stream.name)?;
                }
                builder.streams = Some(streams);
            }
            Binding::Attrs(body) => apply_tag_attrs(builder, body, tag)?,
        }
    }
    builders.into_iter().map(finish_tag).collect()
}

/// Apply a comma-list of tag attributes onto a builder, erroring on any
/// duplicate key.
fn apply_tag_attrs(builder: &mut TagBuilder, body: &str, tag: &str) -> Result<()> {
    for field in body.split(',') {
        if field.trim().is_empty() {
            continue;
        }
        // Split on `=` without trimming the value: separator values are taken
        // as written (so a literal single space survives) and `parse_separator`
        // trims only as a fallback. Other attrs trim their own value below.
        let (key, value) = field.split_once('=').with_context(|| {
            format!(
                "tag `{tag}` attribute {:?} is not `key=value`",
                field.trim()
            )
        })?;
        let key = key.trim();
        match key {
            "qual" => {
                let value = value.trim();
                if builder.qual.is_some() {
                    bail!("tag `{tag}` attribute `qual` is set more than once");
                }
                builder.qual = Some(if value == "none" {
                    QualTag::None
                } else {
                    validate_sam_tag(value)?;
                    QualTag::Named(value.to_string())
                });
            }
            "sep" => {
                if builder.sep.is_some() {
                    bail!("tag `{tag}` attribute `sep` is set more than once");
                }
                builder.sep = Some(parse_separator(value, tag, "sep")?);
            }
            "qual-sep" => {
                if builder.qual_sep.is_some() {
                    bail!("tag `{tag}` attribute `qual-sep` is set more than once");
                }
                builder.qual_sep = Some(parse_separator(value, tag, "qual-sep")?);
            }
            "raw" => {
                let value = value.trim();
                if builder.raw.is_some() {
                    bail!("tag `{tag}` attribute `raw` is set more than once");
                }
                builder.raw = Some(parse_bool(value, "raw")?);
            }
            other => bail!("tag `{tag}` has unknown attribute `{other}`"),
        }
    }
    Ok(())
}

/// Parse a tag separator value: a single literal character or a `%XX` escape
/// (after decoding it must be exactly one character; empty and multi-character
/// separators are rejected).
fn parse_separator(value: &str, tag: &str, attr: &str) -> Result<String> {
    // Honor the value as written first, so a literal single space (`sep= `) is
    // a valid separator.
    let decoded = decode_percent(value)?;
    if decoded.chars().count() == 1 {
        return Ok(decoded);
    }
    // Not a single character as written: maybe surrounding whitespace (e.g. a
    // trailing space before a comma in an attr list) padded it, so trim and
    // retry. A genuine single-space separator already matched above, so this
    // never strips an intended whitespace separator.
    let trimmed = decode_percent(value.trim())?;
    if trimmed.chars().count() == 1 {
        return Ok(trimmed);
    }
    bail!(
        "tag `{tag}` `{attr}` must be a single character (a literal char or a %XX escape), got {value:?}"
    );
}

/// The built-in seq->qual tag map; used when a tag binding leaves `qual` at its
/// default.
pub(crate) fn default_qual_tag(tag: &str) -> Option<&'static str> {
    match tag {
        "CB" => Some("CY"),
        // CR is the raw cellular barcode; its quality has no dedicated SAM tag,
        // so it shares CY (the cell-barcode quality) by convention, letting
        // `--tag CR=bc::raw=true` pair a quality with no explicit `qual=`.
        "CR" => Some("CY"),
        "RX" => Some("QX"),
        "BC" => Some("QT"),
        "OX" => Some("BZ"),
        _ => None,
    }
}

/// Finalize a tag builder: require streams and enforce the separator rules.
fn finish_tag(builder: TagBuilder) -> Result<TagBinding> {
    let streams = builder.streams.with_context(|| {
        format!(
            "tag `{}` has attributes but no streams (a `{}=STREAMS` is required)",
            builder.tag, builder.tag
        )
    })?;
    let qual = builder.qual.unwrap_or(QualTag::Default);
    let emits_qual = match &qual {
        QualTag::None => false,
        QualTag::Named(_) => true,
        QualTag::Default => default_qual_tag(&builder.tag).is_some(),
    };
    // A multi-stream tag joins its streams with `sep` (default `-`); when it
    // also emits a quality tag the qualities join with `qual-sep` (default a
    // single space). A single-stream tag joins nothing, so neither separator is
    // needed (any explicitly set value is kept but unused).
    let sep = if streams.len() >= 2 {
        Some(builder.sep.unwrap_or_else(|| "-".to_string()))
    } else {
        builder.sep
    };
    let qual_sep = if streams.len() >= 2 && emits_qual {
        Some(builder.qual_sep.unwrap_or_else(|| " ".to_string()))
    } else {
        builder.qual_sep
    };
    Ok(TagBinding {
        tag: builder.tag,
        streams,
        qual,
        sep,
        qual_sep,
        raw: builder.raw.unwrap_or(false),
    })
}

/// Parse `--rg-tag K=V` lists, rejecting `SM`/`LB` and duplicate keys.
fn parse_rg_tags(specs: &[String]) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    let mut seen = HashSet::new();
    for spec in specs {
        for field in spec.split(',') {
            let field = field.trim();
            if field.is_empty() {
                continue;
            }
            let (key, value) = field
                .split_once('=')
                .with_context(|| format!("invalid `--rg-tag` {field:?}; expected `K=V`"))?;
            let key = key.trim().to_string();
            let value = value.trim().to_string();
            if key == "SM" || key == "LB" {
                bail!("--rg-tag may not set `{key}`; SM and LB come from the sample/sub_sample");
            }
            if !seen.insert(key.clone()) {
                bail!("--rg-tag key `{key}` is set more than once");
            }
            out.push((key, value));
        }
    }
    Ok(out)
}

/// Resolve the mutually exclusive sample flags into a single [`SampleSpec`].
fn parse_sample_spec(args: &DemuxArgs) -> Result<SampleSpec> {
    // The CLI also enforces this via clap, but a config file or a direct
    // library call could supply more than one; reject that here rather than
    // silently resolving by precedence.
    let set = [
        args.sample_from_group.is_some(),
        args.sample_sheet.is_some(),
        !args.samples.is_empty(),
    ]
    .into_iter()
    .filter(|x| *x)
    .count();
    if set > 1 {
        bail!("--sample, --sample-sheet, and --sample-from-group are mutually exclusive; set only one");
    }
    if let Some(group) = &args.sample_from_group {
        let group = group.trim();
        validate_ident(group, "sample-from-group")?;
        return Ok(SampleSpec::FromGroup(group.to_string()));
    }
    if let Some(sheet) = &args.sample_sheet {
        return Ok(SampleSpec::Sheet(sheet.clone()));
    }
    if !args.samples.is_empty() {
        let mut samples = Vec::new();
        let mut seen = HashSet::new();
        for spec in &args.samples {
            let sample = parse_sample(spec)?;
            let key = (sample.sample.clone(), sample.sub_sample.clone());
            if !seen.insert(key) {
                bail!(
                    "sample `{}`{} is defined more than once (combine its members in one selector)",
                    sample.sample,
                    match &sample.sub_sample {
                        Some(SubSample::Literal(s)) => format!("::{s}"),
                        Some(SubSample::Pool) => "::%pool".to_string(),
                        None => String::new(),
                    }
                );
            }
            samples.push(sample);
        }
        return Ok(SampleSpec::Inline(samples));
    }
    Ok(SampleSpec::None)
}

/// Parse one `--sample SAMPLE[::SUB_SAMPLE]=SELECTOR`.
fn parse_sample(spec: &str) -> Result<Sample> {
    let (left, selector) = spec.split_once('=').with_context(|| {
        format!("invalid `--sample` {spec:?}; expected `SAMPLE[::SUB_SAMPLE]=SELECTOR`")
    })?;
    let (sample, sub_sample) = match left.split_once("::") {
        Some((s, sub)) => {
            let sub = sub.trim();
            let sub_sample = if sub.is_empty() {
                None
            } else if sub == "%pool" {
                Some(SubSample::Pool)
            } else {
                Some(SubSample::Literal(sub.to_string()))
            };
            (s.trim(), sub_sample)
        }
        None => (left.trim(), None),
    };
    validate_ident(sample, "sample name")?;
    Ok(Sample {
        sample: sample.to_string(),
        sub_sample,
        selector: parse_selector(selector.trim())?,
    })
}

/// Parse a `--sample` selector: `+`-joined AND terms, each `group::members` or
/// a bare `group`.
fn parse_selector(value: &str) -> Result<Selector> {
    let mut terms = Vec::new();
    for term in value.split('+') {
        let term = term.trim();
        if term.is_empty() {
            bail!("empty selector term in {value:?}");
        }
        let selector = match term.split_once("::") {
            Some((group, members)) => {
                let group = group.trim();
                validate_ident(group, "selector group")?;
                let members: Vec<String> = members
                    .split(',')
                    .map(|m| m.trim().to_string())
                    .filter(|m| !m.is_empty())
                    .collect();
                if members.is_empty() {
                    bail!("group selector `{group}` lists no members");
                }
                GroupSelector {
                    group: group.to_string(),
                    members,
                }
            }
            None => {
                validate_ident(term, "selector group")?;
                GroupSelector {
                    group: term.to_string(),
                    members: Vec::new(),
                }
            }
        };
        terms.push(selector);
    }
    Ok(Selector { terms })
}

/// Parse the `--remove` flag bodies.
fn parse_removes(specs: &[String]) -> Result<Vec<RemoveRule>> {
    specs.iter().map(|s| parse_remove(s)).collect()
}

/// Parse one `--remove SEL[=PATTERN]`.
fn parse_remove(spec: &str) -> Result<RemoveRule> {
    let (sel, pattern) = match spec.split_once('=') {
        Some((sel, pattern)) => (sel.trim(), Some(pattern.trim())),
        None => (spec.trim(), None),
    };
    if sel.is_empty() {
        bail!("--remove requires a selector before `=`: {spec:?}");
    }
    let selector = match sel.split_once("::") {
        Some((group, id)) => {
            let group = group.trim();
            let id = id.trim();
            validate_ident(group, "remove group")?;
            if id.is_empty() {
                bail!("--remove `{sel}` has an empty tag id after `::`");
            }
            RemoveSelector::GroupTag {
                group: group.to_string(),
                id: id.to_string(),
            }
        }
        None => {
            validate_ident(sel, "remove group")?;
            RemoveSelector::Group(sel.to_string())
        }
    };
    let pattern = match pattern {
        Some(p) => Some(parse_output_pattern(
            p,
            &[Placeholder::Pool, Placeholder::Source],
        )?),
        None => None,
    };
    Ok(RemoveRule { selector, pattern })
}

/// Parse an output pattern, expanding `%placeholder` tokens and `%XX` escapes,
/// and rejecting any placeholder not allowed in this flag's context.
pub(crate) fn parse_output_pattern(raw: &str, allowed: &[Placeholder]) -> Result<OutputPattern> {
    const PLACEHOLDERS: &[Placeholder] = &[
        // Listed longest-first so a longer placeholder name is never
        // pre-empted by a shorter one that prefixes it.
        Placeholder::SubSample,
        Placeholder::Sample,
        Placeholder::Source,
        Placeholder::Ordinal,
        Placeholder::Pool,
    ];
    let mut segments: Vec<PatternSegment> = Vec::new();
    let mut literal = String::new();
    let mut i = 0;
    while i < raw.len() {
        let rest = &raw[i..];
        if let Some(after) = rest.strip_prefix('%') {
            if let Some(placeholder) = PLACEHOLDERS
                .iter()
                .find(|p| after.starts_with(p.name()))
                .copied()
            {
                if !allowed.contains(&placeholder) {
                    bail!(
                        "placeholder `%{}` is not allowed in this pattern: {raw:?}",
                        placeholder.name()
                    );
                }
                if !literal.is_empty() {
                    segments.push(PatternSegment::Literal(std::mem::take(&mut literal)));
                }
                segments.push(PatternSegment::Placeholder(placeholder));
                i += 1 + placeholder.name().len();
                continue;
            }
            let bytes = after.as_bytes();
            if bytes.len() >= 2 && bytes[0].is_ascii_hexdigit() && bytes[1].is_ascii_hexdigit() {
                let byte = u8::from_str_radix(&after[..2], 16).unwrap();
                literal.push(byte as char);
                i += 3;
                continue;
            }
            bail!("unrecognized `%` escape in pattern {raw:?} (use a known placeholder or a %XX hex escape)");
        }
        let c = rest.chars().next().unwrap();
        literal.push(c);
        i += c.len_utf8();
    }
    if !literal.is_empty() {
        segments.push(PatternSegment::Literal(literal));
    }
    if segments.is_empty() {
        bail!("output pattern is empty");
    }
    Ok(OutputPattern { segments })
}

/// Parse a `--metrics-*` path as a pool-level output pattern. These TSVs are
/// one row (or one block) per pool, so only `%pool` is meaningful;
/// `%sample`/`%sub_sample`/`%ordinal`/`%source` cannot resolve for a metrics
/// file and are rejected fail-fast (via [`parse_output_pattern`]) rather than
/// passed through literally. `None` (the flag was omitted) stays `None`.
fn parse_metrics_pattern(path: Option<&Path>, flag: &str) -> Result<Option<OutputPattern>> {
    let Some(path) = path else { return Ok(None) };
    let Some(raw) = path.to_str() else {
        bail!("{flag} path is not valid UTF-8: {path:?}");
    };
    let pattern = parse_output_pattern(raw, &[Placeholder::Pool])
        .with_context(|| format!("{flag}: only `%pool` is valid in a metrics path"))?;
    Ok(Some(pattern))
}

/// Decode `%XX` percent escapes in a value (e.g. `%20` -> a space). A lone `%`
/// not starting a valid escape is kept literal.
fn decode_percent(value: &str) -> Result<String> {
    let bytes = value.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%'
            && i + 2 < bytes.len()
            && bytes[i + 1].is_ascii_hexdigit()
            && bytes[i + 2].is_ascii_hexdigit()
        {
            out.push(u8::from_str_radix(&value[i + 1..i + 3], 16).unwrap());
            i += 3;
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    String::from_utf8(out)
        .with_context(|| format!("value is not valid UTF-8 after percent-decoding: {value:?}"))
}

/// Parse a non-negative integer attribute value.
fn parse_usize(value: &str, attr: &str) -> Result<usize> {
    value.parse().with_context(|| {
        format!("invalid `{attr}` value {value:?} (expected a non-negative integer)")
    })
}

/// Parse a boolean attribute value (`true`/`false`; never `0`/`1`).
fn parse_bool(value: &str, attr: &str) -> Result<bool> {
    match value {
        "true" => Ok(true),
        "false" => Ok(false),
        other => bail!("`{attr}` must be `true` or `false`, got `{other}`"),
    }
}

/// Whether a character is valid in an identifier (group/tag/sample name).
fn is_ident_char(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_'
}

/// Validate a bare identifier (group/sample name): non-empty, `[A-Za-z0-9_]`.
fn validate_ident(name: &str, kind: &str) -> Result<()> {
    if name.is_empty() {
        bail!("{kind} is empty");
    }
    if !name.chars().all(is_ident_char) {
        bail!("{kind} `{name}` has invalid characters (allowed: letters, digits, `_`)");
    }
    Ok(())
}

/// Parse one `+`-joined stream reference: an optional leading `~`
/// (reverse-complement this stream's contribution at this use site) followed by
/// the stream name.
fn parse_stream_ref(token: &str) -> StreamRef {
    match token.strip_prefix('~') {
        Some(name) => StreamRef {
            name: name.trim().to_string(),
            revcomp: true,
        },
        None => StreamRef {
            name: token.to_string(),
            revcomp: false,
        },
    }
}

/// Validate a stream name (`--extract`/`--template`/`--tag` streams): like an
/// identifier but `.` is allowed (e.g. `umi.r1`).
fn validate_stream_name(name: &str) -> Result<()> {
    if name.is_empty() {
        bail!("stream name is empty");
    }
    if !name.chars().all(|c| is_ident_char(c) || c == '.') {
        bail!("stream name `{name}` has invalid characters (allowed: letters, digits, `_`, `.`)");
    }
    Ok(())
}

/// Validate a two-character SAM tag (`[A-Za-z][A-Za-z0-9]`).
fn validate_sam_tag(tag: &str) -> Result<()> {
    let bytes = tag.as_bytes();
    let ok = bytes.len() == 2 && bytes[0].is_ascii_alphabetic() && bytes[1].is_ascii_alphanumeric();
    if !ok {
        bail!("`{tag}` is not a valid SAM tag (must be a letter then a letter or digit)");
    }
    Ok(())
}

/// Validate the `--qc-tag` tag: a valid SAM tag AND in the spec's local-use
/// namespace (starts with `X`/`Y`/`Z`, or contains a lowercase letter).
/// Two-uppercase tags are reserved by the spec for future standardization and
/// must not be squatted.
fn validate_qc_tag(tag: &str) -> Result<String> {
    validate_sam_tag(tag)?;
    let local = tag.starts_with(['X', 'Y', 'Z']) || tag.chars().any(|c| c.is_ascii_lowercase());
    if !local {
        bail!(
            "--qc-tag `{tag}` must be a SAM local-use tag (start with X/Y/Z or contain a lowercase letter); two-uppercase tags are reserved by the SAM spec"
        );
    }
    Ok(tag.to_string())
}

/// Validate cross-entity references and attribute conflicts on the assembled
/// plan.
fn validate_plan(plan: &DemuxPlan) -> Result<()> {
    let group_names: HashSet<&str> = plan.groups.iter().map(|g| g.name.as_str()).collect();
    let extract_names: HashSet<&str> = plan.extracts.iter().map(|e| e.name.as_str()).collect();
    let groups_with_match: HashSet<&str> = plan
        .groups
        .iter()
        .filter(|g| g.attrs.match_streams.is_some())
        .map(|g| g.name.as_str())
        .collect();

    // Group-level attribute conflicts and reference existence.
    for group in &plan.groups {
        let attrs = &group.attrs;
        if attrs.match_streams.is_some() && attrs.loc.is_some() {
            bail!(
                "group `{}` sets both `match` and `loc`; a match= group has no read window",
                group.name
            );
        }
        if let Some(anchor) = attrs.anchor {
            // An anchored placement fixes each tag's anchored edge at a loc
            // boundary and matches over the tag's own length via a direct
            // kernel. Reject the combinations that kernel does not (yet) model,
            // fail-fast, rather than silently mis-matching.
            let (edge, boundary) = match anchor {
                Anchor::ThreePrime => ("3'", "loc.end"),
                Anchor::FivePrime => ("5'", "loc.start"),
            };
            if attrs.match_streams.is_some() {
                bail!(
                    "group `{}` cannot combine anchor= with match= (a match= group has no read window to anchor in)",
                    group.name
                );
            }
            // An anchored group needs its boundary: either its own `loc`, or a
            // relative `next=GROUP:lo-hi` window (which the engine supplies at
            // match time).
            let is_relative_next_target = plan.groups.iter().any(|g| {
                g.attrs
                    .next
                    .as_ref()
                    .is_some_and(|n| n.group == group.name && n.window.is_some())
            });
            if attrs.loc.is_none() && !is_relative_next_target {
                bail!(
                    "group `{}` is anchored but has no {edge} anchor ({boundary}): give it a `loc`, or make it a relative `next=GROUP:lo-hi` window target",
                    group.name
                );
            }
            if attrs.partial5.is_some() || attrs.partial3.is_some() {
                bail!(
                    "group `{}` cannot combine anchor= with partial5/partial3 (anchoring fixes the {edge} edge)",
                    group.name
                );
            }
            if attrs.revcomp == Some(true) {
                bail!(
                    "group `{}` anchor= supports the forward strand only (both_strands is not supported)",
                    group.name
                );
            }
            // anchor=3p uses a substitution-only kernel for now (no 3'-pinned
            // indel DP yet), so reject an indel budget rather than silently
            // ignoring it.
            if anchor == Anchor::ThreePrime && attrs.dist.as_ref().is_some_and(|d| d.indel > 0) {
                bail!(
                    "group `{}` anchor=3p is substitution-only; set indel=0 in dist",
                    group.name
                );
            }
        }
        if let Some(next) = &attrs.next {
            require_group(
                &group_names,
                &next.group,
                &format!("group `{}` next=", group.name),
            )?;
            // A relative `next` window replaces the DOWNSTREAM group's `loc`,
            // so that group may not also pin an absolute `loc`. The upstream
            // group keeps its own `loc` (where it searches to find its match
            // end); only the named target is constrained.
            if next.window.is_some() {
                // The window is computed from the upstream group's match
                // coordinates, but a `match=` group matches a synthetic
                // joined-stream text with no read coordinates, so it can only
                // be a bare ordering link, never the source of a relative
                // window.
                if attrs.match_streams.is_some() {
                    bail!(
                        "group `{}` is a match= group and cannot drive a relative `next` window (it has no read coordinates); use a bare `next={}`",
                        group.name,
                        next.group
                    );
                }
                if let Some(target) = plan.groups.iter().find(|g| g.name == next.group) {
                    if target.attrs.loc.is_some() {
                        bail!(
                            "group `{}` has an absolute `loc` but is the target of a relative `next` window from `{}` (the window replaces its `loc`)",
                            next.group,
                            group.name
                        );
                    }
                    // The window replaces the target's `loc`, but a `match=`
                    // group matches a synthetic joined-stream text and has no
                    // read window for the window to replace, so aiming a
                    // relative window at one is meaningless (symmetric to the
                    // match=-source rejection).
                    if target.attrs.match_streams.is_some() {
                        bail!(
                            "group `{}` aims a relative `next` window at match= group `{}`, which has no read window; use a bare `next={}`",
                            group.name,
                            next.group,
                            next.group
                        );
                    }
                }
            }
        }
        if let Some(prev) = &attrs.prev {
            require_group(&group_names, prev, &format!("group `{}` prev=", group.name))?;
        }
        if let Some(streams) = &attrs.match_streams {
            for stream in streams {
                require_extract(
                    &extract_names,
                    stream,
                    &format!("group `{}` match=", group.name),
                )?;
            }
        }
    }

    // A downstream group may have only one upstream `next` link; two groups
    // both pointing `next=` at it is an ambiguous ordering (whose window
    // wins?), so fail fast naming the conflict.
    let mut next_targets: HashSet<&str> = HashSet::new();
    for group in &plan.groups {
        if let Some(next) = &group.attrs.next {
            if !next_targets.insert(next.group.as_str()) {
                bail!(
                    "group `{}` is the target of more than one `next=` link; a group may follow only one upstream group",
                    next.group
                );
            }
        }
    }

    // Extract anchor references.
    for extract in &plan.extracts {
        for group in anchor_groups(&extract.span.body) {
            require_group(&group_names, group, &format!("extract `{}`", extract.name))?;
            if groups_with_match.contains(group) {
                bail!(
                    "extract `{}` anchors on `@{group}`, but `{group}` is a match= group (it has no read span to anchor)",
                    extract.name
                );
            }
        }
    }

    // Template and tag stream references.
    for (i, template) in plan.templates.iter().enumerate() {
        for stream in &template.streams {
            require_extract(
                &extract_names,
                &stream.name,
                &format!("template #{}", i + 1),
            )?;
        }
    }
    for tag in &plan.tags {
        for stream in &tag.streams {
            require_extract(&extract_names, &stream.name, &format!("tag `{}`", tag.tag))?;
        }
    }

    // Selector group references (samples + removes).
    match &plan.samples {
        SampleSpec::Inline(samples) => {
            for sample in samples {
                // A `+` selector AND-joins terms across DIFFERENT groups;
                // naming the same group twice asks one group's single best
                // match to be two tags at once, which no read can satisfy (the
                // route would be dead). Reject it and point at the OR-pool
                // form.
                let mut seen_groups = HashSet::new();
                for term in &sample.selector.terms {
                    require_group(
                        &group_names,
                        &term.group,
                        &format!("sample `{}` selector", sample.sample),
                    )?;
                    if !seen_groups.insert(term.group.as_str()) {
                        bail!(
                            "sample `{}` selector names group `{}` more than once; combine its tags as an OR pool (`{}::a,b`), not an AND (`+`)",
                            sample.sample,
                            term.group,
                            term.group
                        );
                    }
                }
            }
        }
        SampleSpec::FromGroup(group) => {
            require_group(&group_names, group, "--sample-from-group")?;
        }
        SampleSpec::Sheet(_) | SampleSpec::None => {}
    }
    for rule in &plan.remove {
        let group = match &rule.selector {
            RemoveSelector::Group(g) => g,
            RemoveSelector::GroupTag { group, .. } => group,
        };
        require_group(&group_names, group, "--remove selector")?;
    }

    // Rule 3 (hard error): `%sample`/`%sub_sample` in `--out` requires sample
    // targets. With no sample spec the placeholder can only ever expand to the
    // empty string (one file), which is almost always a forgotten `--sample`.
    // (A `--sample-sheet`/`--sample-from-group` always resolves to >=1 target,
    // or errors at load.)
    if matches!(plan.samples, SampleSpec::None) {
        if let Some(out) = &plan.out {
            if out.uses(Placeholder::Sample) || out.uses(Placeholder::SubSample) {
                bail!(
                    "`--out` uses `%sample`/`%sub_sample` but no samples are defined; add `--sample`, `--sample-sheet`, or `--sample-from-group` (or remove the placeholder)"
                );
            }
        }
    }

    Ok(())
}

/// The groups an extraction span body anchors on.
pub(crate) fn anchor_groups(body: &SpanBody) -> Vec<&str> {
    match body {
        SpanBody::File { .. } => Vec::new(),
        SpanBody::AnchorMatch { group } => vec![group.as_str()],
        SpanBody::AnchorOffset { group, .. } => vec![group.as_str()],
        SpanBody::Between { from, to } => vec![from.as_str(), to.as_str()],
    }
}

/// Error if `name` is not a declared group.
fn require_group(groups: &HashSet<&str>, name: &str, referrer: &str) -> Result<()> {
    if !groups.contains(name) {
        bail!("{referrer} references undefined group `{name}`");
    }
    Ok(())
}

/// Error if `name` is not a declared extract stream.
fn require_extract(extracts: &HashSet<&str>, name: &str, referrer: &str) -> Result<()> {
    if !extracts.contains(name) {
        bail!("{referrer} references undefined extract stream `{name}`");
    }
    Ok(())
}

/// Non-fatal plan warnings, surfaced by [`parse_demux`] via `log::warn!`.
/// Returned as strings so the detection logic is unit-testable without a
/// logger. Two cases: a no-sample run that still defines groups (rule 2: it
/// does no demultiplexing), and any single group nothing consumes (rule 1).
fn plan_warnings(plan: &DemuxPlan) -> Vec<String> {
    let mut warnings = Vec::new();

    // Rule 2: groups defined but no sample target -> no demultiplexing happens
    // (everything passes through to `--out`). Structurally indistinguishable
    // from an intentional annotate-only run, so this can only be a warning,
    // never an error.
    if matches!(plan.samples, SampleSpec::None) && !plan.groups.is_empty() {
        warnings.push(
            "groups are defined but no `--sample`, `--sample-sheet`, or `--sample-from-group` was \
             given; reads will not be demultiplexed (they pass through to `--out` undemultiplexed). \
             Add a sample target to name a single sample or split by sample."
                .to_string(),
        );
    }

    // Rule 1: a group that nothing consumes has no effect on output. "Consumed"
    // is read broadly (a find constraint, a next/prev link, or a match= all do
    // real routing work) so the warning never fires on a load-bearing group;
    // only a genuinely inert one is flagged. A `--sample-sheet` may reference
    // any group (resolved later, at load), so suppress the warning entirely
    // under a sheet.
    if !matches!(plan.samples, SampleSpec::Sheet(_)) {
        let emitted: HashSet<&str> = plan
            .templates
            .iter()
            .flat_map(|t| t.streams.iter())
            .chain(plan.tags.iter().flat_map(|t| t.streams.iter()))
            .map(|s| s.name.as_str())
            .collect();
        let match_refs: HashSet<&str> = plan
            .groups
            .iter()
            .filter_map(|g| g.attrs.match_streams.as_ref())
            .flatten()
            .map(String::as_str)
            .collect();
        let mut linked: HashSet<&str> = HashSet::new();
        for group in &plan.groups {
            if let Some(next) = &group.attrs.next {
                linked.insert(group.name.as_str());
                linked.insert(next.group.as_str());
            }
            if let Some(prev) = &group.attrs.prev {
                linked.insert(group.name.as_str());
                linked.insert(prev.as_str());
            }
        }
        for group in &plan.groups {
            if group_is_consumed(plan, group, &emitted, &match_refs, &linked) {
                continue;
            }
            warnings.push(format!(
                "group `{}` is defined but nothing consumes its matches (no sample selector, \
                 `--remove`, find constraint, `next`/`prev` link, `match=`, or extract emitted by a \
                 `--tag`/`--template`); it has no effect on output.",
                group.name
            ));
        }
    }

    // An extract is consumed only when a `--tag`/`--template` emits its stream
    // or a `match=` group matches on it. One that nothing references is
    // computed and then dropped (the body falls back to pass-through), with no
    // effect on output, so flag it rather than silently doing nothing.
    let consumed_streams: HashSet<&str> = plan
        .templates
        .iter()
        .flat_map(|t| t.streams.iter())
        .chain(plan.tags.iter().flat_map(|t| t.streams.iter()))
        .map(|s| s.name.as_str())
        .chain(
            plan.groups
                .iter()
                .filter_map(|g| g.attrs.match_streams.as_ref())
                .flatten()
                .map(String::as_str),
        )
        .collect();
    for extract in &plan.extracts {
        if !consumed_streams.contains(extract.name.as_str()) {
            warnings.push(format!(
                "extract `{}` is defined but unused (no `--tag`/`--template`/`match=` consumes it); \
                 it has no effect on output.",
                extract.name
            ));
        }
    }

    // A multi-stream `--template` concatenates its streams into ONE output read
    // (`+` combines streams, as in `--tag`). Separate output reads (e.g. an
    // R1/R2 pair) each need their own `--template`, so a `+` here may be an
    // unintended concat where a list of reads was meant; warn and show the
    // per-read fix.
    for template in &plan.templates {
        if template.streams.len() > 1 {
            let per_read = template
                .streams
                .iter()
                .map(|s| format!("--template {s}"))
                .collect::<Vec<_>>()
                .join(" ");
            warnings.push(format!(
                "`--template {}` assembles {} streams into ONE output read (`+` concatenates \
                 streams, as in `--tag`). For separate output reads (e.g. an R1/R2 pair), use one \
                 `--template` per read: `{}`.",
                join_streams(&template.streams, "+"),
                template.streams.len(),
                per_read,
            ));
        }
    }

    // `raw=true` differs from the default only when a stream has a corrected
    // form, i.e. it is an `@grp` self-span match (a keeplist group; `match=`
    // groups cannot be `@grp`-anchored). Flag a `raw=true` tag or template that
    // can never differ, so a user does not believe they are emitting an
    // observed sequence that is in fact identical to the corrected output.
    let correctable: HashSet<&str> = plan
        .extracts
        .iter()
        .filter(|extract| matches!(extract.span.body, SpanBody::AnchorMatch { .. }))
        .map(|extract| extract.name.as_str())
        .collect();
    let warn_raw = |kind: &str, label: String, streams: &[StreamRef]| -> Option<String> {
        if streams
            .iter()
            .any(|s| correctable.contains(s.name.as_str()))
        {
            return None;
        }
        Some(format!(
            "{kind} `{label}` sets `raw=true` but none of its streams ({}) are error-corrected, so \
             it emits the same observed bases as the default; `raw` differs from corrected only for \
             a matched `@grp` self-span stream.",
            join_streams(streams, ",")
        ))
    };
    for tag in &plan.tags {
        if tag.raw {
            warnings.extend(warn_raw("tag", tag.tag.clone(), &tag.streams));
        }
    }
    for template in &plan.templates {
        if template.raw {
            warnings.extend(warn_raw(
                "--template",
                join_streams(&template.streams, ","),
                &template.streams,
            ));
        }
    }

    warnings
}

/// Whether anything in the plan consumes `group`'s matches (see
/// [`plan_warnings`] rule 1). `emitted` is every extract-stream name a
/// `--template`/`--tag` emits, `match_refs` every stream a `match=` group
/// consumes, and `linked` every group on either end of a `next=`/`prev=` link.
fn group_is_consumed(
    plan: &DemuxPlan,
    group: &Group,
    emitted: &HashSet<&str>,
    match_refs: &HashSet<&str>,
    linked: &HashSet<&str>,
) -> bool {
    let name = group.name.as_str();
    let attrs = &group.attrs;

    // A sample selector / `--sample-from-group` names it.
    match &plan.samples {
        SampleSpec::Inline(samples)
            if samples
                .iter()
                .any(|s| s.selector.terms.iter().any(|t| t.group == name)) =>
        {
            return true;
        }
        SampleSpec::FromGroup(g) if g == name => return true,
        _ => {}
    }
    // A `--remove` rule names it.
    if plan.remove.iter().any(|r| {
        let g = match &r.selector {
            RemoveSelector::Group(g) => g,
            RemoveSelector::GroupTag { group, .. } => group,
        };
        g == name
    }) {
        return true;
    }
    // A find constraint can reject reads on this group (the documented keeplist
    // idiom), so it gates routing globally even when no sample/tag names it.
    if attrs.min_finds_per_group.is_some_and(|v| v >= 1)
        || attrs.max_finds_per_group.is_some()
        || attrs.min_finds_per_tag.is_some_and(|v| v >= 1)
        || attrs.max_finds_per_tag.is_some()
    {
        return true;
    }
    // It orders/gates another group (next/prev), or is a match= group doing
    // cross-stream work.
    if linked.contains(name) || attrs.match_streams.is_some() {
        return true;
    }
    // An `--extract` anchored to it feeds a `--tag`/`--template`
    // (annotate/transform) or a `match=`.
    plan.extracts.iter().any(|e| {
        anchor_groups(&e.span.body).contains(&name)
            && (emitted.contains(e.name.as_str()) || match_refs.contains(e.name.as_str()))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal [`DemuxArgs`] with a single input, overridable per test.
    fn args() -> DemuxArgs {
        DemuxArgs {
            pool: None,
            inputs: vec!["0=reads.fq".to_string()],
            groups: Vec::new(),
            extracts: Vec::new(),
            templates: Vec::new(),
            tags: Vec::new(),
            rg_tags: Vec::new(),
            samples: Vec::new(),
            sample_sheet: None,
            sample_from_group: None,
            require_samples_explain_all_tags: false,
            remove: Vec::new(),
            out: None,
            unassigned: None,
            metrics_per_sample: None,
            metrics_summary: None,
            qc_tag: None,
            per_record: false,
            compression: 5,
            threads: 8,
            command_line: None,
        }
    }

    #[test]
    fn test_parse_inputs_contiguous_ok() {
        let inputs = parse_inputs(&[
            "0=a.fq".to_string(),
            "2=c.fq".to_string(),
            "1=b.fq".to_string(),
        ])
        .unwrap();
        assert_eq!(
            inputs,
            vec![
                PathBuf::from("a.fq"),
                PathBuf::from("b.fq"),
                PathBuf::from("c.fq")
            ]
        );
    }

    #[test]
    fn test_parse_inputs_gap_is_error() {
        let err = parse_inputs(&["0=a.fq".to_string(), "2=c.fq".to_string()]).unwrap_err();
        assert!(err.to_string().contains("contiguous"), "{err}");
    }

    #[test]
    fn test_parse_inputs_duplicate_index_is_error() {
        let err = parse_inputs(&["0=a.fq".to_string(), "0=b.fq".to_string()]).unwrap_err();
        assert!(err.to_string().contains("more than once"), "{err}");
    }

    #[test]
    fn test_parse_inputs_empty_is_error() {
        let err = parse_inputs(&[]).unwrap_err();
        assert!(err.to_string().contains("no input reads"), "{err}");
    }

    #[test]
    fn test_qc_tag_accepts_local_use_and_rejects_reserved() {
        // The default ZS and any lowercase-bearing local-use tag are accepted;
        // a two-uppercase tag (reserved by the SAM spec for future standard
        // use) is rejected fail-fast.
        for tag in ["ZS", "Xx", "ZQ"] {
            let mut args = args();
            args.qc_tag = Some(tag.to_string());
            assert_eq!(parse_demux(&args).unwrap().qc_tag.as_deref(), Some(tag));
        }
        let mut args = args();
        args.qc_tag = Some("SC".to_string());
        let err = parse_demux(&args).unwrap_err().to_string();
        assert!(err.contains("local-use"), "{err}");
    }

    #[test]
    fn test_out_sample_placeholder_without_samples_errors() {
        // Rule 3: `%sample` with no sample targets can only expand to "" ->
        // hard error.
        let mut a = args();
        a.out = Some("out/%sample.fq".to_string());
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("no samples are defined"), "{err}");
    }

    #[test]
    fn test_out_subsample_placeholder_without_samples_errors() {
        let mut a = args();
        a.out = Some("out/%sample.%sub_sample.fq".to_string());
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("no samples are defined"), "{err}");
    }

    #[test]
    fn test_out_sample_placeholder_with_sample_ok() {
        let mut a = args();
        a.groups = vec!["g={ACGT}".to_string()];
        a.samples = vec!["s=g::ACGT".to_string()];
        a.out = Some("out/%sample.fq".to_string());
        assert!(parse_demux(&a).is_ok());
    }

    #[test]
    fn test_warn_multi_stream_template_concatenates() {
        // A `+` `--template` concatenates into ONE read; warn and point at the
        // per-read form.
        let mut a = args();
        a.extracts = vec!["x=0:0:4".to_string(), "y=0:4:8".to_string()];
        a.templates = vec!["x+y".to_string()];
        a.out = Some("out.fq".to_string());
        let plan = parse_demux(&a).unwrap();
        let w = plan_warnings(&plan);
        assert!(
            w.iter()
                .any(|m| m.contains("ONE output read") && m.contains("--template x --template y")),
            "{w:?}"
        );
    }

    #[test]
    fn test_warn_groups_without_sample_is_undemuxed() {
        // Rule 2: groups but no sample spec -> a (non-fatal) "no
        // demultiplexing" warning.
        let mut a = args();
        a.groups = vec!["g={ACGT}".to_string()];
        let plan = parse_demux(&a).unwrap();
        let w = plan_warnings(&plan);
        assert!(
            w.iter().any(|m| m.contains("will not be demultiplexed")),
            "{w:?}"
        );
    }

    #[test]
    fn test_no_warning_for_pure_passthrough() {
        // Pure pass-through (no groups, no sample) is a valid mode and must not
        // warn.
        let mut a = args();
        a.out = Some("out.fq".to_string());
        let plan = parse_demux(&a).unwrap();
        assert!(
            plan_warnings(&plan).is_empty(),
            "{:?}",
            plan_warnings(&plan)
        );
    }

    #[test]
    fn test_warn_unused_extract() {
        // An extract that no --tag/--template/match= consumes is computed and
        // dropped (pass-through), so flag it. With `--template r` it is
        // consumed and the warning must not fire.
        let mut a = args();
        a.extracts = vec!["r=0:9:end".to_string()];
        a.out = Some("out.fq".to_string());
        let plan = parse_demux(&a).unwrap();
        assert!(
            plan_warnings(&plan)
                .iter()
                .any(|m| m.contains("extract `r` is defined but unused")),
            "{:?}",
            plan_warnings(&plan)
        );

        let mut consumed = args();
        consumed.extracts = vec!["r=0:9:end".to_string()];
        consumed.templates = vec!["r".to_string()];
        consumed.out = Some("out.fq".to_string());
        let plan = parse_demux(&consumed).unwrap();
        assert!(
            !plan_warnings(&plan).iter().any(|m| m.contains("unused")),
            "a consumed extract must not warn: {:?}",
            plan_warnings(&plan)
        );
    }

    #[test]
    fn test_warn_inert_group_with_no_consumer() {
        // Rule 1: `g` is defined but consumed by nothing (a sample names a
        // different group).
        let mut a = args();
        a.groups = vec![
            "g={ACGT}".to_string(),
            "g::loc=0:0,dist=0".to_string(),
            "g2={AAAA}".to_string(),
        ];
        a.samples = vec!["s=g2::AAAA".to_string()];
        a.out = Some("out/%sample.fq".to_string());
        let plan = parse_demux(&a).unwrap();
        let w = plan_warnings(&plan);
        assert!(
            w.iter()
                .any(|m| m.contains("group `g`") && m.contains("nothing consumes")),
            "{w:?}"
        );
        // The consumed group `g2` is not flagged.
        assert!(!w.iter().any(|m| m.contains("group `g2`")), "{w:?}");
    }

    #[test]
    fn test_warn_raw_tag_on_uncorrected_stream() {
        // raw=true on a positional (never-corrected) stream can never differ
        // from the default, so it is flagged; raw=true on a corrected `@grp`
        // self-span is not.
        let mut a = args();
        a.groups = vec!["g={ACGT}".to_string(), "g::loc=0:0,dist=1".to_string()];
        a.extracts = vec!["pos=0:0:4".to_string(), "bc=@g".to_string()];
        a.tags = vec![
            "OX=pos".to_string(),
            "OX::raw=true".to_string(),
            "CR=bc".to_string(),
            "CR::raw=true".to_string(),
        ];
        a.samples = vec!["s=g::ACGT".to_string()];
        a.out = Some("out/%sample.bam".to_string());
        let plan = parse_demux(&a).unwrap();
        let w = plan_warnings(&plan);
        assert!(
            w.iter()
                .any(|m| m.contains("tag `OX`") && m.contains("raw=true")),
            "raw on the positional stream is flagged: {w:?}"
        );
        assert!(
            !w.iter().any(|m| m.contains("tag `CR`")),
            "raw on the corrected @grp stream is not flagged: {w:?}"
        );
    }

    #[test]
    fn test_minfinds_group_is_consumed_no_inert_warning() {
        // A `minFindsPerGroup` keeplist gate is load-bearing -> not inert (no
        // rule-1 warning).
        let mut a = args();
        a.groups = vec!["g={ACGT}".to_string(), "g::minFindsPerGroup=1".to_string()];
        a.out = Some("out.fq".to_string());
        let plan = parse_demux(&a).unwrap();
        assert!(
            !plan_warnings(&plan)
                .iter()
                .any(|m| m.contains("nothing consumes")),
            "{:?}",
            plan_warnings(&plan)
        );
    }

    #[test]
    fn test_group_consumed_via_tag_extract_no_inert_warning() {
        // Annotate-only: group -> @g extract -> CB tag. `g` is consumed even
        // without a sample (this is the footgun's structural twin, so it warns
        // about no-demux but NOT about an inert group).
        let mut a = args();
        a.groups = vec!["g={ACGT}".to_string(), "g::loc=0:0".to_string()];
        a.extracts = vec!["bc=@g".to_string()];
        a.tags = vec!["CB=bc".to_string()];
        a.out = Some("out.fq".to_string());
        let plan = parse_demux(&a).unwrap();
        let w = plan_warnings(&plan);
        assert!(!w.iter().any(|m| m.contains("nothing consumes")), "{w:?}");
        assert!(
            w.iter().any(|m| m.contains("will not be demultiplexed")),
            "{w:?}"
        );
    }

    #[test]
    fn test_parse_dist_single_field_sets_total() {
        assert_eq!(
            parse_dist("1").unwrap(),
            Dist {
                mismatch: 1,
                indel: 0,
                total: 1
            }
        );
    }

    #[test]
    fn test_parse_dist_two_fields_total_is_sum() {
        assert_eq!(
            parse_dist("1:1").unwrap(),
            Dist {
                mismatch: 1,
                indel: 1,
                total: 2
            }
        );
    }

    #[test]
    fn test_parse_dist_three_fields_explicit() {
        assert_eq!(
            parse_dist("2:1:2").unwrap(),
            Dist {
                mismatch: 2,
                indel: 1,
                total: 2
            }
        );
    }

    #[test]
    fn test_parse_dist_inconsistent_total_is_error() {
        // total 3 exceeds mismatch+indel (1+1).
        let err = parse_dist("1:1:3").unwrap_err();
        assert!(err.to_string().contains("inconsistent"), "{err}");
    }

    #[test]
    fn test_parse_dist_too_many_fields_is_error() {
        assert!(parse_dist("1:1:1:1").is_err());
    }

    #[test]
    fn test_parse_location_file_start() {
        assert_eq!(
            parse_location("1:0").unwrap(),
            Location {
                file: 1,
                start: Endpoint::FromStart(0),
                end: None
            }
        );
    }

    #[test]
    fn test_parse_location_negative_file_index_is_error() {
        // A negative file index (formerly `-1` = "any file") is rejected: omit
        // `loc` to search every file. This also avoids overloading the sign,
        // which means "from the read end" in start/end.
        let err = parse_location("-1:5:end").unwrap_err().to_string();
        assert!(err.contains("file index"), "names the bad field: {err}");
    }

    #[test]
    fn test_parse_location_negative_start() {
        assert_eq!(
            parse_location("0:-10:end").unwrap(),
            Location {
                file: 0,
                start: Endpoint::FromEnd(10),
                end: Some(Endpoint::ReadEnd)
            }
        );
    }

    #[test]
    fn test_parse_span_file_to_end() {
        let span = parse_span("0:9:end").unwrap();
        assert_eq!(
            span.body,
            SpanBody::File {
                file: 0,
                start: Endpoint::FromStart(9),
                end: Endpoint::ReadEnd
            }
        );
    }

    #[test]
    fn test_parse_span_negative_end() {
        let span = parse_span("0:5:-2").unwrap();
        assert_eq!(
            span.body,
            SpanBody::File {
                file: 0,
                start: Endpoint::FromStart(5),
                end: Endpoint::FromEnd(2)
            }
        );
    }

    #[test]
    fn test_parse_span_rejects_tilde() {
        // `~` is no longer a span property; it moved to the point of use
        // (`--tag T=~stream`).
        let err = parse_span("~0:0:9").unwrap_err().to_string();
        assert!(err.contains("not allowed on an `--extract` span"), "{err}");
    }

    #[test]
    fn test_parse_span_anchor_match() {
        let span = parse_span("@grp_cbt").unwrap();
        assert_eq!(
            span.body,
            SpanBody::AnchorMatch {
                group: "grp_cbt".to_string()
            }
        );
    }

    #[test]
    fn test_parse_span_anchor_offset() {
        let span = parse_span("@grp_cbt+19:9").unwrap();
        assert_eq!(
            span.body,
            SpanBody::AnchorOffset {
                group: "grp_cbt".to_string(),
                before: false,
                offset: 19,
                len: AnchorLen::Bases(9)
            }
        );
    }

    #[test]
    fn test_parse_span_anchor_offset_before() {
        // `@grp-0:9`: the 9 bp immediately to the left of the match (a UMI
        // upstream of a barcode).
        let span = parse_span("@grp_cbt-0:9").unwrap();
        assert_eq!(
            span.body,
            SpanBody::AnchorOffset {
                group: "grp_cbt".to_string(),
                before: true,
                offset: 0,
                len: AnchorLen::Bases(9)
            }
        );
    }

    #[test]
    fn test_parse_span_anchor_offset_to_end() {
        let span = parse_span("@grp_cbt+28:end").unwrap();
        assert_eq!(
            span.body,
            SpanBody::AnchorOffset {
                group: "grp_cbt".to_string(),
                before: false,
                offset: 28,
                len: AnchorLen::ToEnd
            }
        );
    }

    #[test]
    fn test_parse_span_between_anchors() {
        let span = parse_span("@grpA..@grpB").unwrap();
        assert_eq!(
            span.body,
            SpanBody::Between {
                from: "grpA".to_string(),
                to: "grpB".to_string()
            }
        );
    }

    #[test]
    fn test_parse_partial_min_freq() {
        let partial = parse_partial("3:0.1", "partial5").unwrap();
        assert_eq!(partial.min_match, 3);
        assert!((partial.max_mismatch_freq - 0.1).abs() < 1e-9);
    }

    #[test]
    fn test_parse_partial_freq_out_of_range_is_error() {
        assert!(parse_partial("3:1.5", "partial5").is_err());
    }

    #[test]
    fn test_parse_next_plain_link() {
        assert_eq!(
            parse_next("grp_cb3").unwrap(),
            NextLink {
                group: "grp_cb3".to_string(),
                window: None
            }
        );
    }

    #[test]
    fn test_parse_next_relative_window() {
        assert_eq!(
            parse_next("grp_cb3:3-5").unwrap(),
            NextLink {
                group: "grp_cb3".to_string(),
                window: Some((3, 5))
            }
        );
    }

    #[test]
    fn test_parse_mode_values() {
        assert_eq!(parse_mode("all").unwrap(), MatchMode::All);
        assert_eq!(parse_mode("nearest").unwrap(), MatchMode::Nearest);
        assert!(parse_mode("closest").is_err());
    }

    #[test]
    fn test_parse_bool_rejects_numeric() {
        assert!(parse_bool("1", "both_strands").is_err());
        assert!(parse_bool("true", "both_strands").unwrap());
    }

    #[test]
    fn test_parse_groups_accumulate_source_then_attrs() {
        let groups = parse_groups(&[
            "grp_cbt=cbt.tsv".to_string(),
            "grp_cbt::loc=1:0,dist=1,maxFindsPerTag=1".to_string(),
        ])
        .unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(
            groups[0].source,
            GroupSource::File(PathBuf::from("cbt.tsv"))
        );
        assert_eq!(
            groups[0].attrs.dist,
            Some(Dist {
                mismatch: 1,
                indel: 0,
                total: 1
            })
        );
        assert_eq!(groups[0].attrs.max_finds_per_tag, Some(1));
    }

    #[test]
    fn test_parse_groups_inline_set_seq_as_id() {
        let groups = parse_groups(&["g={AACGT, sci_8=ACAGGCG}".to_string()]).unwrap();
        let GroupSource::Inline(tags) = &groups[0].source else {
            panic!("expected inline source");
        };
        assert_eq!(
            tags[0],
            InlineTag {
                id: "AACGT".to_string(),
                seq: "AACGT".to_string()
            }
        );
        assert_eq!(
            tags[1],
            InlineTag {
                id: "sci_8".to_string(),
                seq: "ACAGGCG".to_string()
            }
        );
    }

    #[test]
    fn test_parse_groups_duplicate_attr_is_error() {
        let err =
            parse_groups(&["g=f.tsv".to_string(), "g::dist=1,dist=2".to_string()]).unwrap_err();
        assert!(err.to_string().contains("set more than once"), "{err}");
    }

    #[test]
    fn test_parse_groups_duplicate_source_is_error() {
        let err = parse_groups(&["g=a.tsv".to_string(), "g=b.tsv".to_string()]).unwrap_err();
        assert!(err.to_string().contains("defined more than once"), "{err}");
    }

    #[test]
    fn test_parse_groups_attrs_without_source_is_error() {
        let err = parse_groups(&["g::dist=1".to_string()]).unwrap_err();
        assert!(err.to_string().contains("no tag set"), "{err}");
    }

    #[test]
    fn test_parse_groups_unknown_attr_is_error() {
        let err = parse_groups(&["g=f.tsv".to_string(), "g::wat=1".to_string()]).unwrap_err();
        assert!(err.to_string().contains("unknown attribute"), "{err}");
    }

    #[test]
    fn test_findone_shorthand() {
        // `findOne` (bare or =true) sets the three finds caps; =false is a
        // no-op; a conflict with an explicit finds cap is a fail-fast.
        let attrs = |spec: &str| {
            parse_groups(&["g={ACGT}".to_string(), format!("g::{spec}")]).unwrap()[0]
                .attrs
                .clone()
        };
        let bare = attrs("findOne");
        assert_eq!(bare.max_finds_per_tag, Some(1));
        assert_eq!(bare.min_finds_per_group, Some(1));
        assert_eq!(bare.max_finds_per_group, Some(1));
        assert_eq!(attrs("findOne=true").max_finds_per_group, Some(1));
        assert_eq!(attrs("findOne=false").max_finds_per_group, None);
        let err = parse_groups(&[
            "g={ACGT}".to_string(),
            "g::findOne,maxFindsPerGroup=2".to_string(),
        ])
        .unwrap_err();
        assert!(err.to_string().contains("set more than once"), "{err}");
    }

    #[test]
    fn test_parse_tags_single_stream_needs_no_sep() {
        let tags = parse_tags(&["BC=cbt".to_string()]).unwrap();
        assert_eq!(tags[0].streams[0].name, "cbt");
        assert!(!tags[0].streams[0].revcomp);
        assert_eq!(tags[0].qual, QualTag::Default);
        assert!(!tags[0].raw, "raw defaults to false (corrected)");
    }

    #[test]
    fn test_parse_tags_raw_attribute() {
        let tags = parse_tags(&["CR=bc".to_string(), "CR::raw=true".to_string()]).unwrap();
        assert!(tags[0].raw, "raw=true sets the raw flag");
        // raw= takes a boolean; a non-boolean is rejected, and a duplicate is
        // rejected.
        assert!(parse_tags(&["CR=bc".to_string(), "CR::raw=yes".to_string()]).is_err());
        let dup = parse_tags(&[
            "CR=bc".to_string(),
            "CR::raw=true".to_string(),
            "CR::raw=false".to_string(),
        ])
        .unwrap_err();
        assert!(dup.to_string().contains("set more than once"), "{dup}");
    }

    #[test]
    fn test_default_qual_tag_includes_raw_cell() {
        // CR (raw cellular barcode) pairs CY, like the corrected CB; OX (raw
        // UMI) pairs BZ.
        assert_eq!(default_qual_tag("CR"), Some("CY"));
        assert_eq!(default_qual_tag("OX"), Some("BZ"));
    }

    #[test]
    fn test_parse_tags_multi_stream_defaults_sep_and_qual_sep() {
        // A multi-stream tag needs no explicit separators: `sep` defaults to
        // `-`, and (because RX has a default quality tag QX) `qual-sep`
        // defaults to a single space.
        let tags = parse_tags(&["RX=a+b".to_string()]).unwrap();
        assert_eq!(tags[0].sep.as_deref(), Some("-"));
        assert_eq!(tags[0].qual_sep.as_deref(), Some(" "));
    }

    #[test]
    fn test_parse_tags_single_stream_has_no_default_separators() {
        // A single-stream tag joins nothing, so neither separator is defaulted.
        let tags = parse_tags(&["RX=a".to_string()]).unwrap();
        assert_eq!(tags[0].sep, None);
        assert_eq!(tags[0].qual_sep, None);
    }

    #[test]
    fn test_parse_tags_qual_none_skips_qual_sep_default() {
        // CB::qual=none drops the quality tag, so qual-sep is not defaulted;
        // sep still defaults to `-`.
        let tags = parse_tags(&["CB=a+b+c".to_string(), "CB::qual=none".to_string()]).unwrap();
        assert_eq!(tags[0].qual, QualTag::None);
        assert_eq!(tags[0].sep.as_deref(), Some("-"));
        assert_eq!(tags[0].qual_sep, None);
    }

    #[test]
    fn test_parse_separator_accepts_literal_space_and_percent_escape() {
        // A literal single space is honored as written; the %20 escape decodes
        // to the same.
        let lit = parse_tags(&["CB=a+b".to_string(), "CB::sep= ".to_string()]).unwrap();
        assert_eq!(lit[0].sep.as_deref(), Some(" "));
        let esc = parse_tags(&["CB=a+b".to_string(), "CB::sep=%20".to_string()]).unwrap();
        assert_eq!(esc[0].sep.as_deref(), Some(" "));
    }

    #[test]
    fn test_parse_separator_trims_only_to_salvage_a_single_char() {
        // A value that is not a single char as written is trimmed and retried
        // (` - ` -> `-`)...
        let ok = parse_tags(&["CB=a+b".to_string(), "CB::sep= - ".to_string()]).unwrap();
        assert_eq!(ok[0].sep.as_deref(), Some("-"));
        // ...but a genuinely multi-character value still errors.
        let err = parse_tags(&["CB=a+b".to_string(), "CB::sep=ab".to_string()]).unwrap_err();
        assert!(err.to_string().contains("single character"), "{err}");
    }

    #[test]
    fn test_parse_tags_qual_sep_percent_decoded() {
        let tags = parse_tags(&[
            "RX=a+b".to_string(),
            "RX::sep=-".to_string(),
            "RX::qual=BZ".to_string(),
            "RX::qual-sep=%20".to_string(),
        ])
        .unwrap();
        assert_eq!(tags[0].qual, QualTag::Named("BZ".to_string()));
        assert_eq!(tags[0].qual_sep.as_deref(), Some(" "));
    }

    #[test]
    fn test_parse_tags_invalid_sam_tag_is_error() {
        assert!(parse_tags(&["C=a".to_string()]).is_err());
        assert!(parse_tags(&["123=a".to_string()]).is_err());
    }

    #[test]
    fn test_parse_tags_empty_separator_is_error() {
        let err = parse_tags(&["CB=a+b".to_string(), "CB::sep=".to_string()]).unwrap_err();
        assert!(err.to_string().contains("single character"), "{err}");
    }

    #[test]
    fn test_parse_tags_multi_char_separator_is_error() {
        let err = parse_tags(&["CB=a+b".to_string(), "CB::sep=abc".to_string()]).unwrap_err();
        assert!(err.to_string().contains("single character"), "{err}");
    }

    #[test]
    fn test_parse_rg_tags_rejects_sm() {
        let err = parse_rg_tags(&["SM=foo".to_string()]).unwrap_err();
        assert!(err.to_string().contains("may not set"), "{err}");
    }

    #[test]
    fn test_parse_rg_tags_comma_list() {
        let tags = parse_rg_tags(&["CN=ACME,PL=ILLUMINA".to_string()]).unwrap();
        assert_eq!(
            tags,
            vec![
                ("CN".to_string(), "ACME".to_string()),
                ("PL".to_string(), "ILLUMINA".to_string())
            ]
        );
    }

    #[test]
    fn test_parse_sample_with_sub_sample_pool() {
        let sample = parse_sample("dna01::%pool=grp_cbt::cbt01,cbt02").unwrap();
        assert_eq!(sample.sample, "dna01");
        assert_eq!(sample.sub_sample, Some(SubSample::Pool));
        assert_eq!(sample.selector.terms.len(), 1);
        assert_eq!(sample.selector.terms[0].group, "grp_cbt");
        assert_eq!(sample.selector.terms[0].members, vec!["cbt01", "cbt02"]);
    }

    #[test]
    fn test_parse_sample_bare_group_and_and_terms() {
        let sample = parse_sample("s1=grp_cbt::a,b+grp_plate").unwrap();
        assert_eq!(sample.sub_sample, None);
        assert_eq!(sample.selector.terms.len(), 2);
        assert!(sample.selector.terms[1].members.is_empty());
        assert_eq!(sample.selector.terms[1].group, "grp_plate");
    }

    #[test]
    fn test_parse_sample_spec_duplicate_key_is_error() {
        let mut a = args();
        a.samples = vec!["s1=g::a".to_string(), "s1=g::b".to_string()];
        let err = parse_sample_spec(&a).unwrap_err();
        assert!(err.to_string().contains("more than once"), "{err}");
    }

    #[test]
    fn test_parse_sample_spec_distinct_sub_samples_ok() {
        let mut a = args();
        a.samples = vec!["s1=g::a".to_string(), "s1::lib2=g::b".to_string()];
        assert!(matches!(parse_sample_spec(&a).unwrap(), SampleSpec::Inline(v) if v.len() == 2));
    }

    #[test]
    fn test_parse_sample_spec_conflicting_flags_is_error() {
        let mut a = args();
        a.samples = vec!["s1=g::a".to_string()];
        a.sample_from_group = Some("g".to_string());
        let err = parse_sample_spec(&a).unwrap_err();
        assert!(err.to_string().contains("mutually exclusive"), "{err}");
    }

    #[test]
    fn test_validate_sample_selector_repeats_group_is_error() {
        // `+` AND-joins distinct groups; naming one group twice is an
        // unsatisfiable (dead) route, since a group has a single best match per
        // read. Reject it at validation.
        let mut a = args();
        a.groups = vec!["g=g.tsv".to_string()];
        a.samples = vec!["s1=g::a+g::b".to_string()];
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("more than once"), "{err}");
    }

    #[test]
    fn test_parse_remove_group_only() {
        let rule = parse_remove("phix").unwrap();
        assert_eq!(rule.selector, RemoveSelector::Group("phix".to_string()));
        assert!(rule.pattern.is_none());
    }

    #[test]
    fn test_parse_remove_group_tag_with_pattern() {
        let rule = parse_remove("grp::bad01=out/removed.%source.fq.gz").unwrap();
        assert_eq!(
            rule.selector,
            RemoveSelector::GroupTag {
                group: "grp".to_string(),
                id: "bad01".to_string()
            }
        );
        assert!(rule.pattern.as_ref().unwrap().uses(Placeholder::Source));
    }

    #[test]
    fn test_parse_remove_empty_selector_is_error() {
        let err = parse_remove("=out/x.fq").unwrap_err();
        assert!(err.to_string().contains("requires a selector"), "{err}");
    }

    #[test]
    fn test_parse_output_pattern_placeholders() {
        let pattern = parse_output_pattern(
            "output/%sample.%sub_sample.bam",
            &[Placeholder::Sample, Placeholder::SubSample],
        )
        .unwrap();
        assert_eq!(
            pattern.segments,
            vec![
                PatternSegment::Literal("output/".to_string()),
                PatternSegment::Placeholder(Placeholder::Sample),
                PatternSegment::Literal(".".to_string()),
                PatternSegment::Placeholder(Placeholder::SubSample),
                PatternSegment::Literal(".bam".to_string()),
            ]
        );
    }

    #[test]
    fn test_parse_output_pattern_disallowed_placeholder_is_error() {
        let err = parse_output_pattern("out/%source.fq", &[Placeholder::Sample]).unwrap_err();
        assert!(err.to_string().contains("not allowed"), "{err}");
    }

    #[test]
    fn test_parse_output_pattern_unknown_escape_is_error() {
        let err = parse_output_pattern("out/%nope.fq", &[Placeholder::Pool]).unwrap_err();
        assert!(err.to_string().contains("unrecognized"), "{err}");
    }

    #[test]
    fn test_parse_output_pattern_hex_escape() {
        let pattern = parse_output_pattern("a%20b", &[]).unwrap();
        assert_eq!(
            pattern.segments,
            vec![PatternSegment::Literal("a b".to_string())]
        );
    }

    #[test]
    fn test_decode_percent_space_and_tab() {
        assert_eq!(decode_percent("%20").unwrap(), " ");
        assert_eq!(decode_percent("a%09b").unwrap(), "a\tb");
        assert_eq!(decode_percent("plain").unwrap(), "plain");
        assert_eq!(decode_percent("100%").unwrap(), "100%");
    }

    #[test]
    fn test_validate_plan_undefined_anchor_group_is_error() {
        let mut a = args();
        a.extracts = vec!["cbt=@missing".to_string()];
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("undefined group"), "{err}");
    }

    #[test]
    fn test_validate_plan_template_unknown_stream_is_error() {
        let mut a = args();
        a.templates = vec!["r1+r2".to_string()];
        let err = parse_demux(&a).unwrap_err();
        assert!(
            err.to_string().contains("undefined extract stream"),
            "{err}"
        );
    }

    #[test]
    fn test_parse_template_raw_attribute() {
        // `--template cdna::raw=true` parses the streams plus the raw flag;
        // bare templates default to raw=false; an unknown template attribute is
        // rejected.
        let names = |t: &Template| t.streams.iter().map(|s| s.name.clone()).collect::<Vec<_>>();
        let plain = parse_templates(&["cdna".to_string()]).unwrap();
        assert_eq!(names(&plain[0]), vec!["cdna".to_string()]);
        assert!(!plain[0].raw);
        let raw = parse_templates(&["cdna::raw=true".to_string()]).unwrap();
        assert!(raw[0].raw);
        let multi = parse_templates(&["a+b::raw=true".to_string()]).unwrap();
        assert_eq!(names(&multi[0]), vec!["a".to_string(), "b".to_string()]);
        assert!(multi[0].raw);
        let err = parse_templates(&["cdna::nope=1".to_string()]).unwrap_err();
        assert!(err.to_string().contains("unknown attribute"), "{err}");
    }

    #[test]
    fn test_validate_plan_match_with_loc_is_error() {
        let mut a = args();
        a.extracts = vec!["i7=1:0:8".to_string()];
        a.groups = vec![
            "sb=meta.tsv".to_string(),
            "sb::match=i7,loc=0:0".to_string(),
        ];
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("both `match` and `loc`"), "{err}");
    }

    #[test]
    fn test_parse_group_anchored() {
        let mut a = args();
        a.groups = vec!["g={ACGT}".to_string(), "g::loc=0:0:4,anchor=5p".to_string()];
        let plan = parse_demux(&a).unwrap();
        let g = plan.groups.iter().find(|g| g.name == "g").unwrap();
        assert_eq!(g.attrs.anchor, Some(Anchor::FivePrime));
        // An unknown anchor value is rejected.
        let mut bad = args();
        bad.groups = vec![
            "g={ACGT}".to_string(),
            "g::loc=0:0:4,anchor=nope".to_string(),
        ];
        assert!(parse_demux(&bad)
            .unwrap_err()
            .to_string()
            .contains("anchor must be"));
    }

    #[test]
    fn test_parse_group_anchor3p() {
        let mut a = args();
        a.groups = vec!["g={ACGT}".to_string(), "g::loc=0:0:4,anchor=3p".to_string()];
        let plan = parse_demux(&a).unwrap();
        let g = plan.groups.iter().find(|g| g.name == "g").unwrap();
        assert_eq!(g.attrs.anchor, Some(Anchor::ThreePrime));
    }

    #[test]
    fn test_parse_group_both_strands() {
        let mut a = args();
        a.groups = vec![
            "g={ACGT}".to_string(),
            "g::loc=0:0:4,both_strands=true".to_string(),
        ];
        let plan = parse_demux(&a).unwrap();
        let g = plan.groups.iter().find(|g| g.name == "g").unwrap();
        assert_eq!(g.attrs.revcomp, Some(true));
    }

    #[test]
    fn test_validate_anchor3p_combos() {
        let err_for = |attr: &str| {
            let mut a = args();
            a.groups = vec!["g={ACGT}".to_string(), format!("g::{attr}")];
            parse_demux(&a).unwrap_err().to_string()
        };
        // Mirrors anchor5p's restrictions: needs a 3' anchor, forward strand
        // only.
        assert!(err_for("anchor=3p").contains("no 3' anchor"));
        assert!(err_for("loc=0:0:4,anchor=3p,both_strands=true").contains("forward strand"));
        // anchor3p is substitution-only: an indel budget is rejected, a
        // mismatch-only budget is ok.
        assert!(err_for("loc=0:0:6,anchor=3p,dist=1:1").contains("substitution-only"));
        let mut ok = args();
        ok.groups = vec![
            "g={ACGT}".to_string(),
            "g::loc=0:0:6,anchor=3p,dist=1".to_string(),
        ];
        assert!(parse_demux(&ok).is_ok());
    }

    #[test]
    fn test_validate_anchored_rejects_unsupported_combos() {
        let err_for = |attr: &str| {
            let mut a = args();
            a.groups = vec!["g={ACGT}".to_string(), format!("g::{attr}")];
            parse_demux(&a).unwrap_err().to_string()
        };
        assert!(err_for("anchor=5p").contains("no 5' anchor"));
        assert!(err_for("loc=0:0:4,anchor=5p,partial5=2:0.1").contains("partial"));
        assert!(err_for("loc=0:0:4,anchor=5p,both_strands=true").contains("forward strand"));
        // An indel budget is now supported under anchor=5p (no longer
        // rejected).
        let mut ok = args();
        ok.groups = vec![
            "g={ACGT}".to_string(),
            "g::loc=0:0:6,anchor=5p,dist=1:1".to_string(),
        ];
        assert!(parse_demux(&ok).is_ok());
        // anchor=5p + match= (the match= group has no loc, so the match= guard
        // fires first).
        let mut m = args();
        m.extracts = vec!["i7=1:0:8".to_string()];
        m.groups = vec![
            "g=meta.tsv".to_string(),
            "g::match=i7,anchor=5p".to_string(),
        ];
        assert!(parse_demux(&m)
            .unwrap_err()
            .to_string()
            .contains("anchor= with match="));
    }

    #[test]
    fn test_validate_plan_anchor_on_match_group_is_error() {
        let mut a = args();
        a.extracts = vec!["i7=1:0:8".to_string(), "x=@sb".to_string()];
        a.groups = vec!["sb=meta.tsv".to_string(), "sb::match=i7".to_string()];
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("match= group"), "{err}");
    }

    #[test]
    fn test_validate_plan_downstream_loc_with_relative_next_is_error() {
        // The relative `next` window replaces g2's loc; g2 may not also pin an
        // absolute loc.
        let mut a = args();
        a.groups = vec![
            "g1=a.tsv".to_string(),
            "g2=b.tsv".to_string(),
            "g1::next=g2:3-5".to_string(),
            "g2::loc=2:0".to_string(),
        ];
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("relative `next` window"), "{err}");
    }

    #[test]
    fn test_validate_plan_upstream_loc_with_relative_next_is_ok() {
        // g1 keeps its own loc (where g1 searches); only the target g2 is
        // constrained.
        let mut a = args();
        a.groups = vec![
            "g1=a.tsv".to_string(),
            "g2=b.tsv".to_string(),
            "g1::next=g2:3-5,loc=0:0".to_string(),
        ];
        assert!(parse_demux(&a).is_ok());
    }

    #[test]
    fn test_validate_plan_match_group_relative_next_is_error() {
        // A match= group has no read coordinates, so it cannot drive a relative
        // `next` window.
        let mut a = args();
        a.extracts = vec!["i7=0:0:8".to_string()];
        a.groups = vec![
            "sb=m.tsv".to_string(),
            "sb::match=i7,next=g2:3-5".to_string(),
            "g2=b.tsv".to_string(),
        ];
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("cannot drive a relative"), "{err}");
    }

    #[test]
    fn test_validate_plan_relative_next_at_match_target_is_error() {
        // The symmetric case: a relative `next` window may not TARGET a match=
        // group either (it has no read window for the window to replace).
        let mut a = args();
        a.extracts = vec!["i7=0:0:8".to_string()];
        a.groups = vec![
            "g1=a.tsv".to_string(),
            "g1::loc=0:0:4,next=sb:3-5".to_string(),
            "sb=m.tsv".to_string(),
            "sb::match=i7".to_string(),
        ];
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("aims a relative"), "{err}");
    }

    #[test]
    fn test_validate_plan_negative_next_window_is_clear_error() {
        // A negative window bound splits to an empty low; the error names the
        // whole window and its non-negative requirement rather than complaining
        // about an "empty string".
        let mut a = args();
        a.groups = vec![
            "g1=a.tsv".to_string(),
            "g2=b.tsv".to_string(),
            "g1::next=g2:-1-3".to_string(),
        ];
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("non-negative"), "{err}");
    }

    #[test]
    fn test_validate_plan_duplicate_next_target_is_error() {
        // Two groups pointing `next=` at g3 is an ambiguous ordering and is
        // rejected.
        let mut a = args();
        a.groups = vec![
            "g1=a.tsv".to_string(),
            "g2=b.tsv".to_string(),
            "g3=c.tsv".to_string(),
            "g1::next=g3".to_string(),
            "g2::next=g3".to_string(),
        ];
        let err = parse_demux(&a).unwrap_err();
        assert!(err.to_string().contains("more than one `next=`"), "{err}");
    }

    #[test]
    fn test_parse_demux_splitseq_roundtrip() {
        let mut a = args();
        a.pool = Some("lib01".to_string());
        a.inputs = vec![
            "0=lib01.R1.fastq.gz".to_string(),
            "1=lib01.R2.fastq.gz".to_string(),
            "2=lib01.I1.fastq.gz".to_string(),
        ];
        a.groups = vec![
            "grp_cbt=assets/cbt.tags.tsv".to_string(),
            "grp_cbt::loc=1:0,dist=1,maxFindsPerTag=1,minFindsPerGroup=1,maxFindsPerGroup=1".to_string(),
            "grp_cb2=assets/cb2.tags.tsv".to_string(),
            "grp_cb2::loc=2:0,dist=1,next=grp_cb3,maxFindsPerTag=1,minFindsPerGroup=1,maxFindsPerGroup=1".to_string(),
            "grp_cb3=assets/cb3.tags.tsv".to_string(),
            "grp_cb3::loc=2:38,dist=1,maxFindsPerTag=1,minFindsPerGroup=1,maxFindsPerGroup=1".to_string(),
        ];
        a.extracts = vec![
            "umi.r1=0:0:9".to_string(),
            "r1=0:9:end".to_string(),
            "cbt=@grp_cbt".to_string(),
            "umi.r2=@grp_cbt+19:9".to_string(),
            "r2=@grp_cbt+28:end".to_string(),
            "cb2=@grp_cb2".to_string(),
            "cb3=@grp_cb3".to_string(),
        ];
        a.templates = vec!["r1".to_string(), "r2".to_string()];
        a.tags = vec![
            "RX=umi.r1+umi.r2".to_string(),
            "RX::sep=-".to_string(),
            "RX::qual=BZ".to_string(),
            "RX::qual-sep=%20".to_string(),
            "CB=cbt+cb2+cb3".to_string(),
            "CB::sep=-".to_string(),
            "CB::qual=none".to_string(),
            "BC=cbt".to_string(),
        ];
        a.rg_tags = vec!["CN=ACME".to_string(), "PL=ILLUMINA".to_string()];
        a.samples = vec![
            "dna01::%pool=grp_cbt::cbt01,cbt02,cbt03,cbt04,cbt05,cbt06".to_string(),
            "dna02::%pool=grp_cbt::cbt07,cbt08,cbt09".to_string(),
        ];
        a.out = Some("output/%sample.%sub_sample.raw.unmapped.bam".to_string());

        let plan = parse_demux(&a).expect("split-seq-like command should parse");
        assert_eq!(plan.inputs.len(), 3);
        assert_eq!(plan.groups.len(), 3);
        assert_eq!(plan.extracts.len(), 7);
        assert_eq!(plan.tags.len(), 3);
        // grp_cb2 -> grp_cb3 sequential link.
        let cb2 = plan.groups.iter().find(|g| g.name == "grp_cb2").unwrap();
        assert_eq!(cb2.attrs.next.as_ref().unwrap().group, "grp_cb3");
        // CB has three streams, dash sep, no quality tag.
        let cb = plan.tags.iter().find(|t| t.tag == "CB").unwrap();
        assert_eq!(cb.streams.len(), 3);
        assert_eq!(cb.qual, QualTag::None);
        match &plan.samples {
            SampleSpec::Inline(v) => assert_eq!(v.len(), 2),
            other => panic!("expected inline samples, got {other:?}"),
        }
        assert!(plan.out.as_ref().unwrap().uses(Placeholder::Sample));
    }

    #[test]
    fn test_parse_demux_dual_index_roundtrip() {
        let mut a = args();
        a.pool = Some("FC1.L001".to_string());
        a.inputs = vec![
            "0=r1.fq.gz".to_string(),
            "1=i1.fq.gz".to_string(),
            "2=i2.fq.gz".to_string(),
            "3=r2.fq.gz".to_string(),
        ];
        a.extracts = vec![
            "i7=1:0:8".to_string(),
            "i5=2:0:8".to_string(),
            "t.r1=0:0:end".to_string(),
            "t.r2=3:0:end".to_string(),
        ];
        a.groups = vec![
            "sample_bc=metadata.tsv".to_string(),
            "sample_bc::match=i7+i5,mode=nearest,dist=1,delta=2".to_string(),
        ];
        a.templates = vec!["t.r1".to_string(), "t.r2".to_string()];
        a.rg_tags = vec!["PL=ILLUMINA".to_string()];
        a.sample_from_group = Some("sample_bc".to_string());
        a.unassigned = Some("out/unmatched.%source.fq.gz".to_string());
        a.out = Some("out/%sample.R%ordinal.fq.gz".to_string());

        let plan = parse_demux(&a).expect("dual-index command should parse");
        let sb = plan.groups.iter().find(|g| g.name == "sample_bc").unwrap();
        assert_eq!(
            sb.attrs.match_streams,
            Some(vec!["i7".to_string(), "i5".to_string()])
        );
        assert_eq!(sb.attrs.mode, Some(MatchMode::Nearest));
        assert_eq!(sb.attrs.delta, Some(2));
        assert_eq!(plan.samples, SampleSpec::FromGroup("sample_bc".to_string()));
        assert!(plan.unassigned.as_ref().unwrap().uses(Placeholder::Source));
        assert!(plan.out.as_ref().unwrap().uses(Placeholder::Ordinal));
    }
}
