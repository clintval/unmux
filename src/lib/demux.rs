//! The demux engine: parse records, extract tags/UMIs (carrying qualities), and
//! fan out to per-sample FASTX/SAM/BAM/CRAM (SAM/BAM/CRAM written unmapped,
//! except an aligned `@RG` split, which keeps its alignments), splitcode-style,
//! on the `sassy` matcher.
//!
//! [`run_demux`] parses the quote-free flag surface into a
//! [`crate::DemuxPlan`], then drives the engine: load each `--group` tag set,
//! compile the matcher groups and the [`crate::fanout::Routing`] table, and,
//! per fragment, match the groups, route the record, and fan it out. An
//! assigned record is written to its per-sample `--out` (carrying `--tag` SAM
//! fields and an `@RG` read group on SAM/BAM); a removed or unmatched record
//! goes to the raw `--remove` / `--unassigned` bins (with `%source` fan-out) or
//! is dropped.
//! The sequential group attributes (next/prev/match/partial/correct), `@PG`
//! provenance + per-sample `@RG`, and the metrics TSVs are all implemented.
//! Every stage runs on one pool of workers (see `crate::pool`): batches of
//! fragments are read, matched, and encoded in parallel, merged into the
//! outputs in input order, and compressed block by block.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, OnceLock, PoisonError};

use anyhow::{anyhow, bail, ensure, Context, Result};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::Header;
use smallvec::SmallVec;
use thousands::Separable;

use crate::chunks::{Chunk, FastqSource};
use crate::extract::{
    extract, span_range, Corrected, Extracted, GroupHit, GroupHits, MatchSpan, Segment,
};
use crate::fanout::{
    compile_routing, expand_from_group, load_sample_sheet_columns, Disposition, RemoveTarget,
    Route, Routing, SheetColumn, SheetColumns, Target,
};
use crate::grammar::{
    column_tag_qualities, default_qual_tag, Anchor, DemuxPlan, Endpoint, GroupSource, MatchMode,
    OutputPattern, Placeholder, QualTag, SampleSpec, SpanBody, TagBinding, Template,
};
use crate::input::{
    Alignment, ChunkFragments, Fragment, FragmentReader, InputRecord, PendingFragment,
    PendingRecord, SniffedFormat,
};
use crate::matcher::{
    match_verdict, match_verdict_after, search_patterns, CompiledGroup, GroupVerdict, MatchStrand,
    Prereq, RelativeWindow, Scratch,
};
use crate::metrics::Metrics;
use crate::output::{
    alignment_header, build_tag_data, check_output_collisions, default_read_group_header,
    insert_string_tag, is_header_text, lacks_qualities, qual_tag_name, read_group_header,
    resolve_metrics_path, resolve_pattern, ConstituentReadGroups, Destination, OutputOrigin,
    PathContext,
};
use crate::pool::{self, End, Readout};
use crate::qc::{removed_slug, routed_slug, unassigned_slug, Unassigned};
use crate::sink::{Block, BlockFile, Sinks};
use crate::tags::{is_column_value, load_tag_file_columns, TagColumns, TagSet};
use crate::writer::{encode_fragment, output_format, OutputFormat, OutputRead, StagedOutputs};

/// Resolved arguments for the `unmux` demux command.
#[derive(Debug, Clone)]
pub struct DemuxArgs {
    /// `--pool`: identifier for the input pool (fills `%pool`; provenance).
    pub pool: Option<String>,
    /// `--in N=PATH`: input record assignments by 0-based file index
    /// (FASTX/SAM/BAM/CRAM; gzip and format auto-detected; `PATH` may be `-`
    /// for stdin, at most once).
    pub inputs: Vec<String>,
    /// `--group`: tag-group definitions / attributes (`NAME=SOURCE` or
    /// `NAME::ATTRS`).
    pub groups: Vec<String>,
    /// `--extract`: extraction spans (`NAME=SPAN`).
    pub extracts: Vec<String>,
    /// `--template`: stream lists assembled into record `SEQ`/`QUAL`.
    pub templates: Vec<String>,
    /// `--tag`: SAM tag bindings / attributes (`TAG=STREAMS` or `TAG::ATTRS`).
    pub tags: Vec<String>,
    /// `--rg-tag`: shared `@RG` header fields (`K=V`: `CN`, `PL`, `PU`, `DS`,
    /// ...) applied to all output read groups; `SM`/`LB` come from the
    /// sample/sub_sample.
    pub rg_tags: Vec<String>,
    /// `--sample`: fan-out targets (`SAMPLE[::SUB_SAMPLE]=SELECTOR`).
    pub samples: Vec<String>,
    /// `--sample-sheet`: table form of `--sample`.
    pub sample_sheet: Option<PathBuf>,
    /// `--sample-sheet-ignore`: sample-sheet columns to skip unread (each must
    /// be in the header and be neither reserved nor a declared group).
    pub sample_sheet_ignore: Vec<String>,
    /// `--sample-from-group`: treat each matched tag of a group as its own
    /// sample (1:1).
    pub sample_from_group: Option<String>,
    /// `--require-samples-explain-all-tags`: error unless every tag of every
    /// group named by a `--sample`/`--sample-sheet` selector is claimed by
    /// some sample. Off by default.
    pub require_samples_explain_all_tags: bool,
    /// `--remove`: record-routing skiplist entries `SEL[=PATTERN]`. A record
    /// matching SEL (a group name, or a `group::id` combination) is removed and
    /// tallied as `removed`; with `=PATTERN` it is also written there (raw
    /// input segments), else dropped. SEL is required.
    pub remove: Vec<String>,
    /// `--out`: output path or path pattern; placeholders
    /// (`%pool`/`%sample`/`%sub_sample`/`%ordinal`) drive the layout
    /// (`%source` applies only to `--unassigned`/`--remove`). `-` or
    /// `/dev/stdout` writes stdout in the input format. Format is otherwise
    /// inferred from the extension (FASTX/SAM/BAM/CRAM); FASTX carries `--tag`
    /// tags in the read-name comment.
    pub out: Option<String>,
    /// `--unassigned`: output path/pattern for records matching no sample
    /// (written as raw input segments; shares the `--out` placeholders).
    pub unassigned: Option<String>,
    /// `--metrics-per-sample`: per-sample metrics TSV
    /// (`pool`/`sample`/`sub_sample` + wide metric columns; headered,
    /// concatenable, joinable to the summary on `pool`). A pool-level file, so
    /// the path may use `%pool` (the only valid placeholder;
    /// `%sample`/`%sub_sample`/`%ordinal`/`%source` are rejected fail-fast).
    pub metrics_per_sample: Option<PathBuf>,
    /// `--metrics-summary`: pool-level summary metrics TSV (`pool` + wide
    /// metric columns; one row per pool; headered, concatenable). The path may
    /// use `%pool` (the only valid placeholder).
    pub metrics_summary: Option<PathBuf>,
    /// `--metrics-unassigned`: unassigned record counts by reason; only `%pool` is valid.
    pub metrics_unassigned: Option<PathBuf>,
    /// `--qc-tag`: emit a per-record JSON demux-provenance slug under this SAM
    /// tag (a local-use tag, default `ZS`). `None` = off (the default).
    pub qc_tag: Option<String>,
    /// `--per-record`: disable auto pair-detection (interleaving and
    /// SAM/BAM/CRAM pairs are auto-detected otherwise).
    pub per_record: bool,
    /// `--compression`: level (0-9) for BGZF (BAM) and gzip (FASTX.gz) outputs;
    /// CRAM uses its own per-block codecs. Output format is inferred from the
    /// file extension.
    pub compression: u8,
    /// `--threads`: worker thread count.
    pub threads: usize,
    /// `--reference`: an indexed FASTA for reading reference-compressed CRAM.
    pub reference: Option<PathBuf>,
    /// `--skip-reference-check`: do not check `--reference` against CRAM `@SQ`
    /// lines.
    pub skip_reference_check: bool,
    /// The full command line, recorded in the output `@PG CL` provenance field.
    /// `None` (e.g. a library/test caller) omits the `CL` field.
    pub command_line: Option<String>,
}

/// Run the `demux` subcommand.
///
/// The raw flag surface is first parsed and structurally validated into a
/// [`crate::DemuxPlan`] (the grammar parser), which surfaces user errors as
/// fail-fast messages. The engine then matches tag groups, extracts streams,
/// assembles `--template` bodies and `--tag` fields, and routes each record to
/// a per-sample fan-out target (with an `@RG` read group on SAM/BAM) or to the
/// `--remove`/`--unassigned` bins, with `@PG` provenance in every header and
/// the metrics TSVs emitted at the end.
pub fn run_demux(args: DemuxArgs) -> Result<()> {
    let plan = crate::grammar::parse_demux(&args)?;
    log::debug!("parsed plan: {plan:?}");
    // The run-shape summary (input/output counts + per-item enumeration) is
    // logged from `run_engine` via `log_run_shape`, once both the inputs
    // (opened + format-sniffed) and the output destinations (resolved from the
    // plan + routing) are known.
    run_engine(&plan, args.command_line.as_deref())
}

/// A short label for a sniffed input format.
fn sniffed_label(format: SniffedFormat) -> &'static str {
    match format {
        SniffedFormat::Fasta => "FASTA",
        SniffedFormat::Fastq => "FASTQ",
        SniffedFormat::Sam => "SAM",
        SniffedFormat::Bam => "BAM",
        SniffedFormat::Cram => "CRAM",
        SniffedFormat::Unknown => "unknown",
        SniffedFormat::Empty => "empty",
    }
}

/// A short label for a resolved output format.
fn output_label(format: OutputFormat) -> &'static str {
    match format {
        OutputFormat::Fastq { gzip: true } => "FASTQ.gz",
        OutputFormat::Fastq { gzip: false } => "FASTQ",
        OutputFormat::Fasta { gzip: true } => "FASTA.gz",
        OutputFormat::Fasta { gzip: false } => "FASTA",
        OutputFormat::Sam => "SAM",
        OutputFormat::Bam => "BAM",
        OutputFormat::Cram => "CRAM",
    }
}

/// Log the run shape once inputs and outputs are both resolved: a summary count
/// line, then one line per indexed input and one per indexed output (never
/// crammed onto one line), so a run's full I/O is auditable from the log alone.
/// stdin shows as `<stdin>`, the implicit stdout sink (no `--out`) as
/// `<stdout>`; each line carries the sniffed/resolved format. `dests` is the
/// directed-file set from [`enumerate_dests`] (per-sample `--out`,
/// `--unassigned`/`--remove` bins); the stdout sink is added when assigned
/// records have no `--out` path of their own. The SAM-tag count dedupes each
/// `--tag` sequence tag, its paired quality tag, and the `--qc-tag` slug by
/// name.
fn log_run_shape(plan: &DemuxPlan, input_formats: &[SniffedFormat], dests: &[PathBuf]) {
    let mut sam_tags: HashSet<&str> = HashSet::new();
    for binding in &plan.tags {
        sam_tags.insert(binding.tag.as_str());
        if let Some(qual_tag) = qual_tag_name(binding) {
            sam_tags.insert(qual_tag);
        }
    }
    sam_tags.extend(plan.qc_tag.as_deref());

    let stdout_sink = plan.out.is_none();
    let n_outputs = dests.len() + usize::from(stdout_sink);
    log::info!(
        "pool={:?}: {} input(s), {} output(s), {} group(s), {} extract(s), {} template(s), {} SAM tag(s)",
        plan.pool,
        plan.inputs.len(),
        n_outputs,
        plan.groups.len(),
        plan.extracts.len(),
        plan.templates.len(),
        sam_tags.len(),
    );

    for (index, path) in plan.inputs.iter().enumerate() {
        let label = if path.as_os_str() == "-" {
            "<stdin>".to_string()
        } else {
            path.display().to_string()
        };
        match input_formats.get(index) {
            Some(&format) => log::info!("  input  {index}: {label} ({})", sniffed_label(format)),
            None => log::info!("  input  {index}: {label}"),
        }
    }

    let mut index = 0;
    if stdout_sink {
        match output_format(None, input_formats) {
            Ok(format) => log::info!("  output {index}: <stdout> ({})", output_label(format)),
            Err(_) => log::info!("  output {index}: <stdout>"),
        }
        index += 1;
    }
    for path in dests {
        let label = if path.as_os_str() == "-" || path.as_os_str() == "/dev/stdout" {
            "<stdout>".to_string()
        } else {
            path.display().to_string()
        };
        match output_format(Some(path), input_formats) {
            Ok(format) => log::info!("  output {index}: {label} ({})", output_label(format)),
            Err(_) => log::info!("  output {index}: {label}"),
        }
        index += 1;
    }
}

/// Load every `--group`'s tag set once, keyed by group name. The matcher (via
/// [`compile_groups`]) and the router (via [`build_routing`]) share these, so a
/// file is read and an inline set built only once. Also returns, by group name,
/// the values of the columns that `@group::COLUMN` extracts read.
fn load_tag_sets(
    plan: &DemuxPlan,
    resolved: &HashMap<String, crate::ResolvedReadGroup>,
) -> Result<(HashMap<String, TagSet>, HashMap<String, TagColumns>)> {
    let mut tag_sets = HashMap::with_capacity(plan.groups.len());
    let mut tag_columns = HashMap::new();
    for group in &plan.groups {
        let wanted: BTreeSet<&str> = plan
            .extracts
            .iter()
            .filter_map(|extract| match &extract.span.body {
                SpanBody::Column { group: of, column } if *of == group.name => {
                    Some(column.as_str())
                }
                _ => None,
            })
            .collect();
        let wanted: Vec<&str> = wanted.into_iter().collect();
        let (tag_set, columns) = match &group.source {
            GroupSource::File(path) => load_tag_file_columns(path, &wanted)
                .with_context(|| format!("failed to load tag set for group `{}`", group.name))?,
            GroupSource::Inline(tags) => {
                let ids: Vec<String> = tags
                    .iter()
                    .map(|tag| if tag.id.is_empty() { &tag.seq } else { &tag.id })
                    .cloned()
                    .collect();
                let columns = wanted
                    .iter()
                    .map(|column| (column.to_string(), ids.clone()))
                    .collect();
                (TagSet::from_inline(tags), columns)
            }
            GroupSource::ReadGroup(_) => (resolved[&group.name].tag_set.clone(), TagColumns::new()),
        };
        tag_sets.insert(group.name.clone(), tag_set);
        tag_columns.insert(group.name.clone(), columns);
    }
    Ok((tag_sets, tag_columns))
}

/// A column stream's lookup: its group's canonical index and one value per tag.
struct ColumnStream<'a> {
    group: usize,
    values: &'a [String],
}

/// Every column stream's lookup, keyed by stream name.
fn resolve_column_streams<'a>(
    plan: &'a DemuxPlan,
    tag_columns: &'a HashMap<String, TagColumns>,
    group_index: &HashMap<String, usize>,
) -> HashMap<&'a str, ColumnStream<'a>> {
    plan.extracts
        .iter()
        .filter_map(|extract| match &extract.span.body {
            SpanBody::Column { group, column } => Some((
                extract.name.as_str(),
                ColumnStream {
                    group: group_index[group],
                    values: &tag_columns[group][column],
                },
            )),
            _ => None,
        })
        .collect()
}

/// Resolve every `@RG` group in the plan against the input header's read
/// groups, once, keyed by group name. Built right after the reader opens (so
/// the header is available) and shared by the tag-set loader, the group
/// compiler, and the routing builder ([`load_tag_sets`], `compile_groups`,
/// `build_routing`).
fn resolve_read_group_groups(
    plan: &DemuxPlan,
    rg_info: &crate::input::ReadGroupInfo,
) -> Result<HashMap<String, crate::ResolvedReadGroup>> {
    let mut resolved = HashMap::new();
    for group in &plan.groups {
        if let GroupSource::ReadGroup(key) = group.source {
            resolved.insert(
                group.name.clone(),
                crate::resolve_read_group(&group.name, rg_info, key)?,
            );
        }
    }
    Ok(resolved)
}

/// A group's sequential prerequisite by upstream NAME, the form built before
/// the final group order (and thus the canonical indices) is known.
/// [`compile_groups`] orders the groups, then resolves each of these to the
/// index-based [`Prereq`] the matching path uses.
struct NamePrereq {
    /// Every upstream group that must have matched first, by name (deduped).
    upstreams: Vec<String>,
    /// The relative window's upstream group name, plus its `lo`/`hi`, or `None`
    /// for a bare/`prev` prerequisite.
    window: Option<(String, usize, usize)>,
}

/// A warning when an `anchor=5p` group has an explicit, concrete `loc` end and
/// a tag whose anchored span (`loc.start + len(tag)`) runs past it: only
/// `loc.start` anchors, so the `loc` end declares the window's intended extent
/// but does not bound an anchored tag (its match runs to the tag's full
/// length). A `:end` / from-end `loc` end is record-length-relative, so no
/// static bound is checked. `None` when there is nothing to warn about.
fn anchor5p_window_overrun_warning(group: &CompiledGroup) -> Option<String> {
    if group.anchor != Some(Anchor::FivePrime) {
        return None;
    }
    let loc = group.loc.as_ref()?;
    let (Endpoint::FromStart(start), Some(Endpoint::FromStart(end))) = (loc.start, loc.end) else {
        return None;
    };
    let width = end.checked_sub(start)?;
    let tag = group.tags.iter().find(|tag| tag.len() > width)?;
    Some(format!(
        "group `{}`: anchor=5p tag of length {} anchored at loc start {} runs past the loc end {} (window width {}); the loc end does not bound an anchored tag, so its match runs to the tag's full length",
        group.name,
        tag.len(),
        start,
        end,
        width
    ))
}

/// The widest window a group searches for its tags, over every record length,
/// or `None` when nothing bounds it: an `anchor=5p` tag runs past its window
/// end, and a `match=` group has no window.
fn widest_window(group: &crate::grammar::Group, groups: &[crate::grammar::Group]) -> Option<usize> {
    if group.attrs.anchor == Some(Anchor::FivePrime) || group.attrs.match_streams.is_some() {
        return None;
    }
    let relative = groups.iter().find_map(|upstream| {
        let next = upstream.attrs.next.as_ref()?;
        (next.group == group.name).then_some(next.window).flatten()
    });
    if let Some((lo, hi)) = relative {
        return Some(hi.saturating_sub(lo));
    }
    let loc = group.attrs.loc.as_ref()?;
    match (loc.start, loc.end.unwrap_or(Endpoint::ReadEnd)) {
        (Endpoint::FromStart(start), Endpoint::FromStart(end)) => Some(end.saturating_sub(start)),
        (Endpoint::FromEnd(start), Endpoint::FromEnd(end)) => Some(start.saturating_sub(end)),
        (Endpoint::FromEnd(start), Endpoint::FromStart(end)) => Some(start.min(end)),
        (Endpoint::FromEnd(start), Endpoint::ReadEnd) => Some(start),
        (Endpoint::ReadEnd, _) => Some(0),
        (Endpoint::FromStart(_), _) => None,
    }
}

/// Compile each `--group` into a matcher [`CompiledGroup`] over its loaded tag
/// set, return them in dependency order (each group after the upstream group
/// its `next`/`prev` link depends on), and the `name -> canonical index` map
/// for that final order. The sequential prerequisites are resolved to canonical
/// indices once the final order (and thus the indices) is known.
fn compile_groups(
    plan: &DemuxPlan,
    tag_sets: &HashMap<String, TagSet>,
    resolved: &HashMap<String, crate::ResolvedReadGroup>,
) -> Result<(Vec<CompiledGroup>, HashMap<String, usize>)> {
    let mut compiled = Vec::with_capacity(plan.groups.len());
    let mut name_prereqs = Vec::with_capacity(plan.groups.len());
    for group in &plan.groups {
        let is_read_group = matches!(group.source, GroupSource::ReadGroup(_));
        let attrs = &group.attrs;
        let dist = attrs.dist.unwrap_or_default();
        let mode = attrs.mode.unwrap_or(MatchMode::All);
        // An `@RG` group's tags are routing tokens, not DNA, so they are not validated.
        if !is_read_group {
            let search = crate::iupac::Search {
                dist,
                mode: attrs.mode,
                anchor: attrs.anchor,
                both_strands: attrs.revcomp.unwrap_or(false),
                window: widest_window(group, &plan.groups),
            };
            for warning in
                crate::iupac::validate_group(&group.name, &tag_sets[&group.name].seqs(), &search)?
            {
                log::warn!("{warning}");
            }
        }
        let tags = tag_sets[&group.name].seqs();
        let mut compiled_group = CompiledGroup {
            name: group.name.clone(),
            patterns: search_patterns(&tags),
            tags,
            dist,
            loc: attrs.loc.clone(),
            revcomp: attrs.revcomp.unwrap_or(false),
            mode,
            delta: attrs.delta.unwrap_or(1),
            min_finds_per_tag: attrs.min_finds_per_tag,
            max_finds_per_tag: attrs.max_finds_per_tag,
            min_finds_per_group: attrs.min_finds_per_group,
            max_finds_per_group: attrs.max_finds_per_group,
            // Resolved to indices below, once the final group order is known.
            prereq: None,
            match_streams: attrs.match_streams.clone(),
            partial5: attrs.partial5,
            partial3: attrs.partial3,
            anchor: attrs.anchor,
            encoded: Vec::new(),
            batched_tags: Default::default(),
            read_group: resolved.get(&group.name).map(|r| r.rg_to_idx.clone()),
        };
        // Precompute the sassy v2 batched search batches (and their covered tag
        // set) for large equal-length tag buckets. Skipped for an `@RG` group:
        // its tags are routing tokens, not DNA, so a sassy search over them is
        // meaningless.
        if !is_read_group {
            compiled_group.encode();
            if let Some(warning) = anchor5p_window_overrun_warning(&compiled_group) {
                log::warn!("{warning}");
            }
        }
        compiled.push(compiled_group);
        name_prereqs.push(prereq_for(group, &plan.groups, &plan.extracts));
    }

    let (ordered, ordered_prereqs) = order_by_prereq(compiled, name_prereqs)?;

    // The canonical index of each group is its position in this final,
    // never-mutated order; group names are unique. Resolve the name-based
    // prerequisites to those indices once, here.
    let group_index: HashMap<String, usize> = ordered
        .iter()
        .enumerate()
        .map(|(idx, group)| (group.name.clone(), idx))
        .collect();
    let mut ordered = ordered;
    for (group, name_prereq) in ordered.iter_mut().zip(ordered_prereqs) {
        group.prereq = resolve_prereq(name_prereq, &group_index);
    }
    Ok((ordered, group_index))
}

/// The estimated share of reads with a barcode error past which shifted tags are warned about.
const WIDE_WINDOW_SHARE: f64 = 0.003;

/// A warning when, often enough to matter, another tag can match reads at a shifted position.
fn wide_window_warning(
    index: usize,
    groups: &[CompiledGroup],
    plan: &DemuxPlan,
    routing: &Routing,
) -> Option<String> {
    let group = &groups[index];
    let dist = group.dist;
    let edits = dist.total > 0 && dist.mismatch + dist.indel > 0;
    if !edits || group.anchor.is_some() || group.read_group.is_some() {
        return None;
    }
    let (width, window) = search_width(group, plan)?;
    let mut lengths = std::collections::BTreeMap::new();
    for tag in group.tags.iter().filter(|tag| !tag.is_empty()) {
        *lengths.entry(tag.len()).or_insert(0usize) += 1;
    }
    let (&shortest, &longest) = (lengths.keys().next()?, lengths.keys().next_back()?);
    let wider = width.checked_sub(shortest).filter(|&wider| wider > 0)?;
    let share = shifted_share(&lengths, width, &dist, group.revcomp);
    let nearest = group.mode == MatchMode::Nearest;
    let best = longest + 2 * dist.indel;
    if share <= WIDE_WINDOW_SHARE
        || (nearest && dist.indel > 0 && width <= best)
        || (!nearest && only_removed(index, groups, plan, routing))
    {
        return None;
    }
    let narrow = if width > longest {
        Some(format!("narrow {window} to {longest} bp"))
    } else if group.match_streams.is_none() {
        Some("pin the tags with anchor=5p".to_string())
    } else {
        None
    };
    let advice = match (nearest, dist.indel > 0, narrow) {
        (false, false, Some(narrow)) => {
            format!("use mode=nearest, or {narrow}, losing reads whose barcode an indel upstream shifts")
        }
        (false, true, _) if width > best => format!(
            "use mode=nearest with {window} narrowed to {best} bp, which measured best with an indel budget"
        ),
        (false, _, _) => "use mode=nearest".to_string(),
        (true, true, _) => {
            format!("narrow {window} to {best} bp, which measured best with an indel budget")
        }
        (true, false, Some(narrow)) => format!("{narrow} if no indel upstream shifts the barcode"),
        (true, false, None) => return None,
    };
    let tags = if shortest == longest {
        "the tags"
    } else {
        "the shortest tag"
    };
    let outcome = if nearest {
        "mode=nearest leaves those reads unmatched"
    } else {
        "the leftmost match wins"
    };
    Some(format!(
        "group `{}`: the {window} window is {wider} bp wider than {tags}, so for an estimated {:.1}% of reads with an error in the barcode another tag matches at a shifted position as closely as the true tag, and {outcome}; {advice}",
        group.name,
        100.0 * share
    ))
}

/// The width a group's tags slide over when it is the same for every record, and its name.
fn search_width(group: &CompiledGroup, plan: &DemuxPlan) -> Option<(usize, String)> {
    let relative = group
        .prereq
        .as_ref()
        .and_then(|prereq| prereq.window.as_ref());
    if let Some(relative) = relative {
        let name = format!("next={}:{}-{}", group.name, relative.lo, relative.hi);
        return Some((relative.hi.saturating_sub(relative.lo), name));
    }
    if let Some(loc) = &group.loc {
        let width = fixed_width(loc.start, loc.end.unwrap_or(Endpoint::ReadEnd))?;
        return Some((width, "loc".to_string()));
    }
    let streams = group.match_streams.as_ref()?;
    let mut width = 0;
    for stream in streams {
        let extract = plan
            .extracts
            .iter()
            .find(|extract| &extract.name == stream)?;
        width += match &extract.span.body {
            SpanBody::File { start, end, .. } => fixed_width(*start, *end)?,
            SpanBody::AnchorOffset {
                len: crate::grammar::AnchorLen::Bases(len),
                ..
            } => *len,
            _ => return None,
        };
    }
    Some((width, format!("match={}", streams.join("+"))))
}

/// The width from `start` to `end` when it is the same for every record long enough to hold it.
fn fixed_width(start: Endpoint, end: Endpoint) -> Option<usize> {
    match (start, end) {
        (Endpoint::FromStart(start), Endpoint::FromStart(end)) => Some(end.saturating_sub(start)),
        (Endpoint::FromEnd(start), Endpoint::FromEnd(end)) => Some(start.saturating_sub(end)),
        (Endpoint::FromEnd(start), Endpoint::ReadEnd) => Some(start),
        _ => None,
    }
}

/// The estimated share of reads with one barcode error that a shifted tag matches as closely.
fn shifted_share(
    lengths: &std::collections::BTreeMap<usize, usize>,
    width: usize,
    dist: &crate::grammar::Dist,
    both_strands: bool,
) -> f64 {
    // A tag's chance of matching random bases as is or by one substitution, insertion or deletion.
    let chance = |len: usize| {
        let bases = len as f64;
        let subs = if dist.mismatch > 0 { 3.0 * bases } else { 0.0 };
        let indels = if dist.indel > 0 {
            4.0 * bases + 2.0
        } else {
            0.0
        };
        (1.0 + subs + indels) / 4f64.powf(bases)
    };
    let claims = |own: usize| -> f64 {
        lengths
            .iter()
            .map(|(&len, &count)| {
                // A tag no longer than the barcode reaches past it only from a shifted start.
                let (others, starts) = if len <= own {
                    (count - usize::from(len == own), width - own)
                } else {
                    (count, (width + 1).saturating_sub(len))
                };
                (others * starts) as f64 * chance(len)
            })
            .sum()
    };
    let worst = lengths
        .keys()
        .filter(|&&own| own <= width)
        .map(|&own| claims(own))
        .fold(0.0, f64::max);
    let strands = if both_strands { 2.0 } else { 1.0 };
    1.0 - (-strands * worst).exp()
}

/// Whether whole-group `--remove` rules alone read a group's match, which any of its tags satisfies.
fn only_removed(
    index: usize,
    groups: &[CompiledGroup],
    plan: &DemuxPlan,
    routing: &Routing,
) -> bool {
    let name = groups[index].name.as_str();
    let removes: Vec<_> = routing
        .removes()
        .iter()
        .filter(|remove| remove.group == index)
        .collect();
    !removes.is_empty()
        && removes.iter().all(|remove| remove.tag.is_none())
        && !routing.routes_on(index)
        && !groups.iter().any(|group| {
            group
                .prereq
                .as_ref()
                .is_some_and(|prereq| prereq.upstreams.contains(&index))
        })
        && !plan
            .extracts
            .iter()
            .any(|extract| crate::grammar::anchor_groups(&extract.span.body).contains(&name))
}

/// Resolve a name-based [`NamePrereq`] into the index-based [`Prereq`] the
/// matching path consumes, using the canonical `name -> index` map. Every
/// referenced name is a declared group, so the lookups always succeed.
fn resolve_prereq(
    prereq: Option<NamePrereq>,
    group_index: &HashMap<String, usize>,
) -> Option<Prereq> {
    let prereq = prereq?;
    let upstreams = prereq
        .upstreams
        .iter()
        .map(|name| group_index[name])
        .collect();
    let window = prereq.window.map(|(name, lo, hi)| RelativeWindow {
        upstream: group_index[&name],
        lo,
        hi,
    });
    Some(Prereq { upstreams, window })
}

/// The sequential prerequisite of `group` (by upstream name): an upstream group
/// whose `next=` link targets it (carrying the link's optional relative
/// window), the upstream group named by this group's own `prev=`, or, for a
/// `match=` group, every group that anchors one of its referenced `--extract`
/// streams (those anchors must have matched before this group's joined value
/// exists).
fn prereq_for(
    group: &crate::grammar::Group,
    groups: &[crate::grammar::Group],
    extracts: &[crate::grammar::Extract],
) -> Option<NamePrereq> {
    let mut upstreams: Vec<String> = Vec::new();
    let mut window = None;
    // Any group whose `next=` targets this one is an upstream prerequisite (the
    // grammar already rejects two `next=` links to the same group, so at most
    // one carries a relative window).
    for upstream in groups {
        if let Some(next) = &upstream.attrs.next {
            if next.group == group.name {
                upstreams.push(upstream.name.clone());
                if let Some((lo, hi)) = next.window {
                    window = Some((upstream.name.clone(), lo, hi));
                }
            }
        }
    }
    // This group's own `prev=` is an additional prerequisite (deduped against a
    // `next=` source that names the same group, since `next`/`prev` are
    // symmetric expressions of one edge).
    if let Some(prev) = &group.attrs.prev {
        if !upstreams.contains(prev) {
            upstreams.push(prev.clone());
        }
    }
    // A `match=` group matches the joined value of named `--extract` streams,
    // so every group those streams anchor on must have matched first. Without
    // this, an anchor group declared after the `match=` group would not yet
    // have a hit when the join is assembled, and the match would route to
    // unassigned. A stream anchored on this same group contributes no
    // dependency (and would be a cycle anyway).
    if let Some(streams) = &group.attrs.match_streams {
        for stream in streams {
            if let Some(spec) = extracts.iter().find(|e| &e.name == stream) {
                for anchor in crate::grammar::anchor_groups(&spec.span.body) {
                    if anchor != group.name.as_str()
                        && !upstreams.iter().any(|u| u.as_str() == anchor)
                    {
                        upstreams.push(anchor.to_string());
                    }
                }
            }
        }
    }
    if upstreams.is_empty() {
        None
    } else {
        Some(NamePrereq { upstreams, window })
    }
}

/// Order compiled groups so every group follows the upstream group it depends
/// on, by repeatedly emitting groups whose prerequisite is already emitted
/// (independent groups keep declaration order). The name-based prerequisites
/// travel alongside their group so the caller can resolve them to canonical
/// indices once this final order is fixed. A `next`/`prev` cycle makes no
/// progress and is a fail-fast error.
fn order_by_prereq(
    groups: Vec<CompiledGroup>,
    prereqs: Vec<Option<NamePrereq>>,
) -> Result<(Vec<CompiledGroup>, Vec<Option<NamePrereq>>)> {
    let mut emitted: Vec<CompiledGroup> = Vec::with_capacity(groups.len());
    let mut emitted_prereqs: Vec<Option<NamePrereq>> = Vec::with_capacity(groups.len());
    let mut emitted_names: HashSet<String> = HashSet::new();
    let mut remaining: Vec<(CompiledGroup, Option<NamePrereq>)> =
        groups.into_iter().zip(prereqs).collect();
    while !remaining.is_empty() {
        let before = remaining.len();
        let mut still = Vec::new();
        for (group, prereq) in remaining {
            let ready = prereq
                .as_ref()
                .is_none_or(|p| p.upstreams.iter().all(|u| emitted_names.contains(u)));
            if ready {
                emitted_names.insert(group.name.clone());
                emitted.push(group);
                emitted_prereqs.push(prereq);
            } else {
                still.push((group, prereq));
            }
        }
        remaining = still;
        if remaining.len() == before {
            let cycle: Vec<&str> = remaining.iter().map(|(g, _)| g.name.as_str()).collect();
            bail!(
                "unmux: cyclic next/prev group dependency among: {}",
                cycle.join(", ")
            );
        }
    }
    Ok((emitted, emitted_prereqs))
}

/// The pool id filling `%pool`, including in `%pool` sub_samples: `--pool` when
/// given, else the longest common prefix of the input file names without their
/// format extensions (trimmed of trailing separators), else `pool`.
fn pool_id(plan: &DemuxPlan) -> String {
    if let Some(pool) = &plan.pool {
        return pool.clone();
    }
    let names: Vec<String> = plan
        .inputs
        .iter()
        .filter_map(|path| path.file_name())
        .map(|name| without_format_extension(&name.to_string_lossy()).to_string())
        .collect();
    let prefix = longest_common_prefix(&names);
    let trimmed = prefix.trim_end_matches(['.', '_', '-', ' ']);
    if trimmed.is_empty() {
        "pool".to_string()
    } else {
        trimmed.to_string()
    }
}

/// The file name without the extension of a format unmux reads, or the name itself.
fn without_format_extension(name: &str) -> &str {
    let without = |name: &str, extension: &str| {
        let stem = name.len().checked_sub(extension.len() + 1)?;
        let dotted = name.get(stem..)?;
        (dotted.starts_with('.') && dotted[1..].eq_ignore_ascii_case(extension)).then_some(stem)
    };
    let inner = &name[..without(name, "gz").unwrap_or(name.len())];
    ["fq", "fastq", "fa", "fasta", "fna", "sam", "bam", "cram"]
        .iter()
        .find_map(|extension| without(inner, extension))
        .filter(|&stem| stem > 0)
        .map_or(name, |stem| &name[..stem])
}

/// The longest common prefix of a set of strings, in whole characters.
fn longest_common_prefix(strings: &[String]) -> String {
    let Some(first) = strings.first() else {
        return String::new();
    };
    let mut end = first.len();
    for other in &strings[1..] {
        end = first
            .bytes()
            .zip(other.bytes())
            .take(end)
            .take_while(|(a, b)| a == b)
            .count();
    }
    while !first.is_char_boundary(end) {
        end -= 1;
    }
    first[..end].to_string()
}

/// The read group id of a record with no sample: the pool, or `pool` when a SAM header cannot hold it.
fn pool_read_group<'a>(
    plan: &DemuxPlan,
    pool: &'a str,
    writes_alignments: bool,
) -> Result<&'a str> {
    if is_header_text(pool) {
        return Ok(pool);
    }
    if writes_alignments && plan.pool.is_some() {
        bail!(
            "--pool `{}` cannot be written to a SAM/BAM/CRAM header; use printable ASCII",
            pool.escape_debug()
        );
    }
    if writes_alignments {
        log::warn!(
            "the input name `{}` cannot be a read group in a SAM/BAM/CRAM header; using `pool` \
             instead, set --pool to choose one",
            pool.escape_debug()
        );
    }
    Ok("pool")
}

/// Build the routing table from the plan's fan-out source (`--sample` inline, a
/// `--sample-sheet`, or `--sample-from-group`) plus the `--remove` rules,
/// sharing the already-loaded tag sets. Also returns the sheet columns that
/// `%sample::COLUMN` extracts read.
fn build_routing(
    plan: &DemuxPlan,
    tag_sets: &HashMap<String, TagSet>,
    group_index: &HashMap<String, usize>,
    pool: &str,
    resolved: &HashMap<String, crate::ResolvedReadGroup>,
) -> Result<(Routing, SheetColumns)> {
    let mut sheet_columns = SheetColumns::new();
    let (samples, demux) = match &plan.samples {
        SampleSpec::None => (Vec::new(), false),
        SampleSpec::Inline(samples) => (samples.clone(), true),
        SampleSpec::Sheet(path, ignore) => {
            let declared: HashSet<&str> = plan.groups.iter().map(|g| g.name.as_str()).collect();
            let ignored: HashSet<&str> = ignore.iter().map(String::as_str).collect();
            let reads: Vec<(&str, &str)> = plan
                .extracts
                .iter()
                .filter_map(|extract| match &extract.span.body {
                    SpanBody::SampleColumn { column } => {
                        Some((extract.name.as_str(), column.as_str()))
                    }
                    _ => None,
                })
                .collect();
            let (samples, columns) = load_sample_sheet_columns(path, &declared, &ignored, &reads)?;
            sheet_columns = columns;
            (samples, true)
        }
        SampleSpec::FromGroup(group) => match resolved.get(group) {
            Some(r) => (r.samples.clone(), true),
            None => {
                let tag_set = tag_sets.get(group).with_context(|| {
                    format!("no loaded tag set for --sample-from-group `{group}`")
                })?;
                (expand_from_group(group, tag_set), true)
            }
        },
    };
    let routing = compile_routing(
        &samples,
        &plan.remove,
        tag_sets,
        group_index,
        pool,
        demux,
        plan.require_samples_explain_all_tags,
    )?;
    Ok((routing, sheet_columns))
}

/// Store each `%sample::COLUMN` stream's value per route on the routes, where
/// `None` is a blank cell. Fails on a value that cannot be a SAM tag or on two
/// sheet rows that route one record to one target with different values.
fn attach_sample_columns(
    plan: &DemuxPlan,
    routing: &mut Routing,
    sheet_columns: &SheetColumns,
) -> Result<()> {
    let is_sheet = matches!(plan.samples, SampleSpec::Sheet(..));
    let mut names = Vec::new();
    let mut columns = Vec::new();
    let mut has_cells = false;
    for extract in &plan.extracts {
        let SpanBody::SampleColumn { column } = &extract.span.body else {
            continue;
        };
        let read = match sheet_columns.get(column) {
            Some(read) => read.clone(),
            None if column.eq_ignore_ascii_case("sub_sample") => SheetColumn::SubSample,
            None => SheetColumn::Sample,
        };
        let values: Vec<Option<String>> = match read {
            SheetColumn::Sample => routing
                .route_targets()
                .map(|target| Some(target.sample.clone()))
                .collect(),
            SheetColumn::SubSample => routing
                .route_targets()
                .map(|target| target.sub_sample.clone())
                .collect(),
            SheetColumn::Cells(values) => {
                has_cells = true;
                values
            }
        };
        for (route, value) in values.iter().enumerate() {
            let Some(value) = value.as_deref().filter(|value| !is_column_value(value)) else {
                continue;
            };
            if is_sheet {
                bail!(
                    "sample sheet row {}: value {value:?} in column {column:?} may hold only \
                     printable ASCII and spaces",
                    route + 1
                );
            }
            bail!(
                "extract `{}` reads {column:?}, but value {value:?} may hold only printable ASCII \
                 and spaces",
                extract.name
            );
        }
        log_blank_sample_column(plan, &extract.name, column, &values);
        names.push(column.as_str());
        columns.push(values);
    }
    routing.set_columns(columns)?;
    if has_cells {
        routing.check_column_conflicts(&names)?;
    }
    Ok(())
}

/// Log how many routes have no value for a `%sample::COLUMN` stream, and the
/// tags their records go without.
fn log_blank_sample_column(
    plan: &DemuxPlan,
    stream: &str,
    column: &str,
    values: &[Option<String>],
) {
    let blank = values.iter().filter(|value| value.is_none()).count();
    let (fallen_back, left_out): (Vec<&TagBinding>, Vec<&TagBinding>) = plan
        .tags
        .iter()
        .filter(|tag| tag.all_streams().any(|s| s.name == stream))
        .partition(|tag| !tag.is_last_choice(stream));
    if blank == 0 || (fallen_back.is_empty() && left_out.is_empty()) {
        return;
    }
    let names = |tags: &[&TagBinding]| {
        tags.iter()
            .map(|tag| tag.tag.as_str())
            .collect::<Vec<_>>()
            .join(" or ")
    };
    let rows = match &plan.samples {
        SampleSpec::Sheet(..) => "sheet rows".to_string(),
        SampleSpec::FromGroup(group) => format!("tags of group `{group}`"),
        _ => "samples".to_string(),
    };
    let mut outcomes = Vec::new();
    if !left_out.is_empty() {
        outcomes.push(format!("get no {} tag", names(&left_out)));
    }
    if !fallen_back.is_empty() {
        outcomes.push(format!(
            "fall back to the next stream for {}",
            names(&fallen_back)
        ));
    }
    log::info!(
        "column {column:?} is blank in {blank} of {} {rows}; records they route {}",
        values.len(),
        outcomes.join(" and ")
    );
}

/// Build the static `Target -> constituent original @RG maps` table for an
/// `@RG`-split run: for each RG id, the target its records route to paired with
/// that read group's original `@RG` map (looked up by id). One target
/// accumulates every original read group whose records land in its file, so the
/// output header can list them verbatim (correct SM/LB/PL/PU/...) instead of a
/// synthesized per-target `@RG`. Empty for any non-`@RG` run.
fn build_constituent_rgs(
    plan: &DemuxPlan,
    resolved: &HashMap<String, crate::ResolvedReadGroup>,
    rg_info: &crate::input::ReadGroupInfo,
) -> ConstituentReadGroups {
    let mut out: ConstituentReadGroups = HashMap::new();
    let SampleSpec::FromGroup(group) = &plan.samples else {
        return out;
    };
    let Some(r) = resolved.get(group) else {
        return out;
    };
    let pool = pool_id(plan);
    let by_id: HashMap<&[u8], &crate::input::InputReadGroup> = rg_info
        .read_groups
        .iter()
        .map(|rg| (rg.id.as_slice(), rg))
        .collect();
    for (rg_id, &tag_idx) in &r.rg_to_idx {
        let target = crate::fanout::resolve_target(&r.samples[tag_idx], &pool);
        if let Some(rg) = by_id.get(rg_id.as_slice()) {
            out.entry(target)
                .or_default()
                .push((rg_id.clone(), rg.map.clone()));
        }
    }
    out
}

/// A SAM/BAM `--out` destination (a path, or `None` for stdout), the fan-out
/// targets whose records land there, and its format.
type OutDestination = (Option<PathBuf>, Vec<Target>, OutputFormat);

/// The SAM/BAM/CRAM `--out` destinations and the fan-out targets whose records
/// land in each (a single-file `--out` collects every target; a `%sample`
/// pattern splits them; no `--out` is stdout), in first-seen order. FASTX
/// destinations have no header, so they are excluded.
fn out_destinations(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
) -> Result<Vec<OutDestination>> {
    let mut by_dest: Vec<OutDestination> = Vec::new();
    let mut index: HashMap<Option<PathBuf>, usize> = HashMap::new();
    for target in routing.targets() {
        let ctx = PathContext {
            pool,
            sample: Some(&target.sample),
            sub_sample: target.sub_sample.as_deref(),
            ordinal: None,
            source: None,
        };
        let path = plan
            .out
            .as_ref()
            .and_then(|pattern| resolve_pattern(pattern, &ctx));
        let format = output_format(path.as_deref(), input_formats)?;
        if matches!(
            format,
            OutputFormat::Sam | OutputFormat::Bam | OutputFormat::Cram
        ) {
            match index.get(&path) {
                Some(&i) => by_dest[i].1.push(target.clone()),
                None => {
                    index.insert(path.clone(), by_dest.len());
                    by_dest.push((path, vec![target.clone()], format));
                }
            }
        }
    }
    Ok(by_dest)
}

/// Every directed output FILE path the engine could write, computed from the
/// plan and routing before streaming so each can be created up front: a
/// directed file always exists, even when it receives no records (an empty but
/// valid file). stdout (a `None` path) is excluded; it cannot be pre-created
/// and is a single sink. [`destination_table`] resolves the same paths for
/// the writes of each record.
fn enumerate_dests(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
    fragment_width: usize,
) -> Result<Vec<PathBuf>> {
    let mut seen: HashSet<PathBuf> = HashSet::new();
    Ok(
        labeled_dests(plan, routing, pool, input_formats, fragment_width)?
            .into_iter()
            .map(|(path, _)| path)
            .filter(|path| seen.insert(path.clone()))
            .collect(),
    )
}

/// Every directed output file path, before deduplication, with what writes it.
fn labeled_dests(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
    fragment_width: usize,
) -> Result<Vec<(PathBuf, OutputOrigin)>> {
    use crate::grammar::RemoveSelector;

    let mut paths: Vec<(PathBuf, OutputOrigin)> = Vec::new();

    // `--out`: one set of paths per assigned target (a `%sample` pattern fans
    // out; a plain pattern dedups to one). A pure pass-through run (no targets)
    // still writes to `--out` with no sample.
    if let Some(pattern) = &plan.out {
        let body_count = if plan.templates.is_empty() {
            fragment_width
        } else {
            plan.templates.len()
        };
        let targets = routing.targets();
        if targets.is_empty() {
            for path in out_paths(pattern, pool, None, body_count, input_formats)? {
                paths.push((path, OutputOrigin::Records("--out".to_string())));
            }
        }
        for target in targets {
            for path in out_paths(pattern, pool, Some(target), body_count, input_formats)? {
                paths.push((path, OutputOrigin::Target(target.clone())));
            }
        }
    }
    // The `--unassigned` and `--remove` bins, which share the `--out` writer when given its pattern.
    let could_share = plan.out.as_ref().is_some_and(|out| {
        !(out.uses(Placeholder::Sample)
            || out.uses(Placeholder::SubSample)
            || out.uses(Placeholder::Ordinal))
    });
    let bin = |pattern: &OutputPattern, option: String| {
        if plan.out.as_ref() == Some(pattern) {
            OutputOrigin::Records(option)
        } else {
            OutputOrigin::Bin {
                option,
                could_share,
            }
        }
    };
    if let Some(pattern) = &plan.unassigned {
        for path in bin_paths(pattern, pool, plan.inputs.len()) {
            paths.push((path, bin(pattern, "--unassigned".to_string())));
        }
    }
    for rule in &plan.remove {
        if let Some(pattern) = &rule.pattern {
            let selector = match &rule.selector {
                RemoveSelector::Group(group) => group.clone(),
                RemoveSelector::GroupTag { group, id } => format!("{group}::{id}"),
            };
            for path in bin_paths(pattern, pool, plan.inputs.len()) {
                paths.push((path, bin(pattern, format!("--remove {selector}"))));
            }
        }
    }
    Ok(paths)
}

/// Refuse outputs that would be one file with each other or with an input, before any is created.
fn check_outputs(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
    fragment_width: usize,
) -> Result<()> {
    let mut outputs = labeled_dests(plan, routing, pool, input_formats, fragment_width)?;
    for (pattern, option) in [
        (&plan.metrics_per_sample, "--metrics-per-sample"),
        (&plan.metrics_summary, "--metrics-summary"),
        (&plan.metrics_unassigned, "--metrics-unassigned"),
    ] {
        if let Some(pattern) = pattern {
            outputs.push((
                resolve_metrics_path(pattern, pool),
                OutputOrigin::Metrics(option),
            ));
        }
    }
    crate::output::check_outputs_over_inputs(&outputs, &read_files(plan))?;
    check_output_collisions(&outputs)
}

/// Every file the run reads: the inputs, the tag files, the sample sheet and the reference.
fn read_files(plan: &DemuxPlan) -> Vec<(PathBuf, crate::output::InputOrigin)> {
    use crate::output::InputOrigin;

    let mut files: Vec<(PathBuf, InputOrigin)> = plan
        .inputs
        .iter()
        .enumerate()
        .map(|(index, path)| (path.clone(), InputOrigin::Reads(index)))
        .collect();
    for group in &plan.groups {
        if let GroupSource::File(path) = &group.source {
            files.push((path.clone(), InputOrigin::TagFile(group.name.clone())));
        }
    }
    if let SampleSpec::Sheet(path, _) = &plan.samples {
        files.push((path.clone(), InputOrigin::SampleSheet));
    }
    if let Some(reference) = &plan.reference {
        let (fai, gzi) = crate::reference::index_paths(reference);
        files.push((reference.clone(), InputOrigin::Reference));
        for index in std::iter::once(fai).chain(gzi) {
            files.push((index, InputOrigin::ReferenceIndex));
        }
    }
    files
}

/// The `--out` file path(s) for one target (or `None` for a pure pass-through),
/// mirroring the destination resolution in `write_body`: a non-alignment
/// template set with more than one body fans out over `%ordinal`, everything
/// else is a single destination.
fn out_paths(
    pattern: &OutputPattern,
    pool: &str,
    target: Option<&Target>,
    body_count: usize,
    input_formats: &[SniffedFormat],
) -> Result<Vec<PathBuf>> {
    let sample = target.map(|t| t.sample.as_str());
    let sub_sample = target.and_then(|t| t.sub_sample.as_deref());
    // The format is fixed by the extension (probe with ordinal 1), independent
    // of the ordinal.
    let probe = PathContext {
        pool,
        sample,
        sub_sample,
        ordinal: Some(1),
        source: None,
    };
    let probe_path = resolve_pattern(pattern, &probe);
    let format = output_format(probe_path.as_deref(), input_formats)?;
    let is_alignment = matches!(
        format,
        OutputFormat::Sam | OutputFormat::Bam | OutputFormat::Cram
    );
    let mut paths = Vec::new();
    if !is_alignment && body_count > 1 {
        for index in 0..body_count {
            let ctx = PathContext {
                pool,
                sample,
                sub_sample,
                ordinal: Some(index + 1),
                source: None,
            };
            if let Some(path) = resolve_pattern(pattern, &ctx) {
                paths.push(path);
            }
        }
    } else {
        let ctx = PathContext {
            pool,
            sample,
            sub_sample,
            ordinal: (body_count == 1).then_some(1),
            source: None,
        };
        if let Some(path) = resolve_pattern(pattern, &ctx) {
            paths.push(path);
        }
    }
    Ok(paths)
}

/// The bin file path(s) for a `--unassigned` / `--remove` pattern, mirroring
/// `write_raw_bin`: with more than one input source and a `%source` pattern the
/// segments fan out one file per source, otherwise everything lands in a single
/// file.
fn bin_paths(pattern: &OutputPattern, pool: &str, num_sources: usize) -> Vec<PathBuf> {
    let mut paths = Vec::new();
    if num_sources > 1 && pattern.uses(Placeholder::Source) {
        for source in 0..num_sources {
            let ctx = PathContext {
                pool,
                source: Some(source),
                ..Default::default()
            };
            if let Some(path) = resolve_pattern(pattern, &ctx) {
                paths.push(path);
            }
        }
    } else {
        let ctx = PathContext {
            pool,
            source: Some(0),
            ..Default::default()
        };
        if let Some(path) = resolve_pattern(pattern, &ctx) {
            paths.push(path);
        }
    }
    paths
}

/// The destinations of a run, by index, and where each kind of write goes.
struct DestinationTable {
    /// Every destination, stdout included when a write can go there.
    list: Arc<[Destination]>,
    /// The `--out` layout of each fan-out target, by [`Route::target_index`].
    targets: Vec<OutLayout>,
    /// The `--out` layout of a record passed through with no fan-out target.
    pass_through: Option<OutLayout>,
    /// The `--unassigned` bin.
    unassigned: Option<BinLayout>,
    /// The bin of each `--remove` rule, in rule order.
    removes: Vec<Option<BinLayout>>,
}

/// Where the bodies of a record bound for `--out` go.
struct OutLayout {
    /// Whether the bodies are written as SAM, BAM, or CRAM records.
    is_alignment: bool,
    /// The destinations, or why a record cannot be written to `--out`.
    dests: std::result::Result<OutDests, String>,
    /// The `@RG` id of a SAM, BAM, or CRAM record.
    read_group: Option<String>,
}

/// The `--out` destination of every body, or of each body by `%ordinal`.
enum OutDests {
    Single(usize),
    PerOrdinal(Vec<usize>),
}

/// The bin destination of every raw segment, or of each by `%source`.
enum BinLayout {
    Single(usize),
    PerSource(Vec<usize>),
}

/// Destinations gathered by path, each listed once in the order first seen.
struct DestinationList<'a> {
    list: Vec<Destination>,
    index: HashMap<Option<PathBuf>, usize>,
    input_formats: &'a [SniffedFormat],
    headers: &'a HashMap<Option<PathBuf>, Header>,
    default_header: &'a Header,
    keeps_alignments: bool,
}

impl DestinationList<'_> {
    /// The index of the destination at `path` (stdout when `None`), added if new.
    fn add(&mut self, path: Option<PathBuf>) -> Result<usize> {
        if let Some(&at) = self.index.get(&path) {
            return Ok(at);
        }
        let dest = Destination::new(
            path.clone(),
            self.input_formats,
            self.headers,
            self.default_header,
        )?;
        self.list.push(dest);
        self.index.insert(path, self.list.len() - 1);
        Ok(self.list.len() - 1)
    }

    /// Add `option`'s bin destination like [`Self::add`], declaring its records' read groups.
    fn add_bin(&mut self, path: Option<PathBuf>, option: &str) -> Result<usize> {
        let at = self.add(path)?;
        let dest = &mut self.list[at];
        if !is_alignment_format(dest.format) || is_dev_null(dest.path.as_deref()) {
            return Ok(at);
        }
        let samples = self.headers.get(&dest.path).map(Header::read_groups);
        for (id, read_group) in self.default_header.read_groups() {
            let sample = samples.and_then(|samples| samples.get(id));
            // Aligned records keep their input read groups, which the samples' header shares.
            if sample.is_some_and(|sample| !self.keeps_alignments || sample != read_group) {
                let place = match &dest.path {
                    Some(path) => format!("in \"{}\"", path.display()),
                    None => "on stdout".to_string(),
                };
                bail!(
                    "{option} records would share read group \"{id}\" with a sample {place}; set \
                     --pool to another name"
                );
            }
            if !dest.header.read_groups().contains_key(id) {
                dest.header
                    .read_groups_mut()
                    .insert(id.clone(), read_group.clone());
            }
        }
        Ok(at)
    }
}

/// Resolve every destination a run writes and the layout of each kind of
/// write, in the order [`labeled_dests`] lists the files. A record with no
/// sample carries the `read_group` id.
fn destination_table(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    read_group: &str,
    body_count: usize,
    list: &mut DestinationList,
) -> Result<DestinationTable> {
    let targets = routing.targets();
    let pass_through = match targets.is_empty() {
        true => Some(out_layout(plan, pool, read_group, None, body_count, list)?),
        false => None,
    };
    let targets = targets
        .iter()
        .map(|target| out_layout(plan, pool, read_group, Some(target), body_count, list))
        .collect::<Result<_>>()?;
    let unassigned = plan
        .unassigned
        .as_ref()
        .map(|pattern| bin_layout(pattern, "--unassigned", pool, plan.inputs.len(), list))
        .transpose()?;
    let removes = routing
        .removes()
        .iter()
        .map(|remove| {
            let option = format!("--remove {}", remove.selector);
            let pattern = remove.pattern.as_ref();
            pattern
                .map(|pattern| bin_layout(pattern, &option, pool, plan.inputs.len(), list))
                .transpose()
        })
        .collect::<Result<_>>()?;
    Ok(DestinationTable {
        list: std::mem::take(&mut list.list).into(),
        targets,
        pass_through,
        unassigned,
        removes,
    })
}

/// The `--out` layout of the bodies of a record routed to `target` (`None`
/// for a pass-through record, which carries the `read_group` id), mirroring
/// [`out_paths`].
fn out_layout(
    plan: &DemuxPlan,
    pool: &str,
    read_group: &str,
    target: Option<&Target>,
    body_count: usize,
    list: &mut DestinationList,
) -> Result<OutLayout> {
    let sample = target.map(|t| t.sample.as_str());
    let sub_sample = target.and_then(|t| t.sub_sample.as_deref());
    let path = |ordinal: Option<usize>| {
        let ctx = PathContext {
            pool,
            sample,
            sub_sample,
            ordinal,
            source: None,
        };
        plan.out.as_ref().and_then(|p| resolve_pattern(p, &ctx))
    };
    let format = output_format(path(Some(1)).as_deref(), list.input_formats)?;
    let is_alignment = matches!(
        format,
        OutputFormat::Sam | OutputFormat::Bam | OutputFormat::Cram
    );
    let multi_fastx = plan
        .out
        .as_ref()
        .filter(|_| !is_alignment && body_count > 1);
    let dests = match multi_fastx {
        Some(pattern) if !pattern.uses(Placeholder::Ordinal) => Err(format!(
            "unmux: multi-read FASTX output needs the %ordinal placeholder in --out to \
             separate the {body_count} template reads"
        )),
        Some(_) => Ok(OutDests::PerOrdinal(
            (1..=body_count)
                .map(|ordinal| list.add(path(Some(ordinal))))
                .collect::<Result<_>>()?,
        )),
        None => Ok(OutDests::Single(
            list.add(path((body_count == 1).then_some(1)))?,
        )),
    };
    Ok(OutLayout {
        is_alignment,
        dests,
        read_group: target
            .map(Target::label)
            .or_else(|| Some(read_group.to_string())),
    })
}

/// The layout of a `--unassigned` or `--remove` bin, mirroring [`bin_paths`].
fn bin_layout(
    pattern: &OutputPattern,
    option: &str,
    pool: &str,
    num_sources: usize,
    list: &mut DestinationList,
) -> Result<BinLayout> {
    let path = |source: usize| {
        let ctx = PathContext {
            pool,
            source: Some(source),
            ..Default::default()
        };
        resolve_pattern(pattern, &ctx)
    };
    if num_sources > 1 && pattern.uses(Placeholder::Source) {
        let dests = (0..num_sources).map(|source| list.add_bin(path(source), option));
        return Ok(BinLayout::PerSource(dests.collect::<Result<_>>()?));
    }
    Ok(BinLayout::Single(list.add_bin(path(0), option)?))
}

/// Precompute the `@RG`-bearing SAM/BAM header for each distinct `--out`
/// destination, so a per-sample file lists its own read groups and a
/// single-file `--out` lists every target's. FASTX destinations and the raw
/// bins are absent from the map and get the default `@PG`-only header.
/// With `input_header` set (aligned records pass through), each header is the
/// input header narrowed to the file's own `@RG` lines instead.
fn out_headers(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
    command_line: Option<&str>,
    constituent: Option<&ConstituentReadGroups>,
    input_header: Option<&Header>,
) -> Result<HashMap<Option<PathBuf>, noodles::sam::Header>> {
    let mut headers = HashMap::new();
    for (dest, targets, _) in out_destinations(plan, routing, pool, input_formats)? {
        let refs: Vec<&Target> = targets.iter().collect();
        let header = match (input_header, constituent) {
            (Some(input), Some(constituent)) => {
                alignment_header(input, Some((&refs, constituent)), command_line)?
            }
            _ => read_group_header(&refs, &plan.rg_tags, command_line, constituent)?,
        };
        headers.insert(dest, header);
    }
    Ok(headers)
}

/// The SAM/BAM/CRAM formats among every destination a run can write: `--out`
/// (a file pattern or stdout) and the `--unassigned`/`--remove` bins.
fn alignment_formats(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
) -> Result<Vec<OutputFormat>> {
    let mut formats: Vec<OutputFormat> = out_destinations(plan, routing, pool, input_formats)?
        .into_iter()
        .filter(|(path, _, _)| !is_dev_null(path.as_deref()))
        .map(|(_, _, format)| format)
        .collect();
    if plan.out.is_none() {
        formats.push(output_format(None, input_formats)?);
    }
    let bins = plan
        .unassigned
        .iter()
        .chain(plan.remove.iter().filter_map(|rule| rule.pattern.as_ref()));
    for pattern in bins {
        let ctx = PathContext {
            pool,
            source: Some(0),
            ..Default::default()
        };
        let path = resolve_pattern(pattern, &ctx);
        if !is_dev_null(path.as_deref()) {
            formats.push(output_format(path.as_deref(), input_formats)?);
        }
    }
    formats.retain(|format| {
        matches!(
            format,
            OutputFormat::Sam | OutputFormat::Bam | OutputFormat::Cram
        )
    });
    Ok(formats)
}

/// Whether any destination of the run is SAM/BAM/CRAM, with or without samples.
fn writes_alignments(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
) -> Result<bool> {
    if !alignment_formats(plan, routing, pool, input_formats)?.is_empty() {
        return Ok(true);
    }
    let Some(pattern) = &plan.out else {
        return Ok(false);
    };
    let ctx = PathContext {
        pool,
        source: Some(0),
        ..Default::default()
    };
    let path = resolve_pattern(pattern, &ctx);
    Ok(matches!(
        output_format(path.as_deref(), input_formats)?,
        OutputFormat::Sam | OutputFormat::Bam | OutputFormat::Cram
    ))
}

/// Whether an `@RG` split passes its aligned records through unchanged: the
/// input is aligned and some destination is SAM/BAM/CRAM. Fails fast on a run
/// that would have to drop alignments to proceed, or that cannot encode them.
fn passes_alignments_through(
    plan: &DemuxPlan,
    aligned_split: bool,
    formats: &[OutputFormat],
) -> Result<bool> {
    if !aligned_split || formats.is_empty() {
        return Ok(false);
    }
    if !plan.templates.is_empty() {
        bail!(
            "the input is aligned (@SQ present), so this @RG split keeps each record whole, but \
             --template would rewrite its sequence and qualities; remove --template, write FASTX \
             instead, or unmap the input first (e.g. `samtools reset`)"
        );
    }
    if plan.inputs.len() > 1 {
        bail!(
            "an @RG split that keeps alignments takes a single input, got {}; merge the inputs \
             first (e.g. `samtools merge`)",
            plan.inputs.len()
        );
    }
    if formats.contains(&OutputFormat::Cram) {
        bail!(
            "writing aligned records to CRAM is not supported; write BAM and convert it (e.g. \
             `samtools view -C -T ref.fa`)"
        );
    }
    let claimed = plan
        .tags
        .iter()
        .flat_map(|binding| std::iter::once(binding.tag.as_str()).chain(qual_tag_name(binding)))
        .chain(plan.qc_tag.as_deref());
    for tag in claimed {
        if crate::input::ALIGNMENT_TAGS.contains(&tag_bytes(tag)) {
            bail!(
                "--tag/--qc-tag `{tag}` would overwrite an alignment tag of records kept whole by \
                 this aligned @RG split; choose another tag"
            );
        }
    }
    let by_sequence = plan
        .groups
        .iter()
        .find(|group| {
            !matches!(group.source, GroupSource::ReadGroup(_))
                && crate::grammar::group_has_effect(plan, group)
        })
        .map(|group| format!("group `{}` matches", group.name))
        .or_else(|| {
            plan.tags
                .first()
                .map(|binding| format!("--tag `{}` is drawn from", binding.tag))
        });
    if let Some(what) = by_sequence {
        bail!(
            "this aligned @RG split keeps secondary and supplementary records, but {what} record \
             sequences, which those records may lack (SEQ `*`) or hard-clip, so they could be \
             routed or tagged apart from their primary; drop it, write FASTX, or unmap the input \
             first (e.g. `samtools reset`)"
        );
    }
    Ok(true)
}

/// Whether `path` discards what is written to it, so its format does not matter.
fn is_dev_null(path: Option<&Path>) -> bool {
    path == Some(Path::new("/dev/null"))
}

/// Refuse FASTX outputs in a run that keeps aligned records whole for its SAM/BAM outputs.
fn refuse_mixed_outputs(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
) -> Result<()> {
    use crate::grammar::RemoveSelector;

    let mut outputs = Vec::new();
    match &plan.out {
        Some(pattern) => {
            for target in routing.targets() {
                let ctx = PathContext {
                    pool,
                    sample: Some(&target.sample),
                    sub_sample: target.sub_sample.as_deref(),
                    ordinal: Some(1),
                    source: None,
                };
                outputs.push(("--out".to_string(), resolve_pattern(pattern, &ctx)));
            }
        }
        None => outputs.push(("--out".to_string(), None)),
    }
    let bin = PathContext {
        pool,
        source: Some(0),
        ..Default::default()
    };
    if let Some(pattern) = &plan.unassigned {
        outputs.push(("--unassigned".to_string(), resolve_pattern(pattern, &bin)));
    }
    for rule in &plan.remove {
        if let Some(pattern) = &rule.pattern {
            let option = match &rule.selector {
                RemoveSelector::Group(group) => format!("--remove {group}"),
                RemoveSelector::GroupTag { group, id } => format!("--remove {group}::{id}"),
            };
            outputs.push((option, resolve_pattern(pattern, &bin)));
        }
    }
    let mut seen = HashSet::new();
    let (mut alignments, mut fastx) = (Vec::new(), Vec::new());
    for (option, path) in outputs {
        if is_dev_null(path.as_deref()) {
            continue;
        }
        let keeps_whole = is_alignment_format(output_format(path.as_deref(), input_formats)?);
        if !seen.insert((keeps_whole, option.clone())) {
            continue;
        }
        let output = match path {
            Some(path) => format!("\"{}\" ({option})", path.display()),
            None => "stdout".to_string(),
        };
        if keeps_whole {
            alignments.push(output);
        } else {
            fastx.push(output);
        }
    }
    if !fastx.is_empty() {
        bail!(
            "this aligned @RG split cannot write both SAM/BAM, which keeps each record whole, and \
             FASTX, which pairs mates: it writes SAM/BAM to {} and FASTX to {}; make every output \
             SAM/BAM, or every output FASTX",
            alignments.join(", "),
            fastx.join(", ")
        );
    }
    Ok(())
}

/// A validated two-character SAM tag name as bytes.
fn tag_bytes(tag: &str) -> [u8; 2] {
    let bytes = tag.as_bytes();
    [bytes[0], bytes[1]]
}

/// The most threads that inflate BAM or decode CRAM input.
pub const MAX_DECODE_THREADS: usize = 4;

/// The BAM inflate / CRAM decode workers for each input of `formats` at
/// `threads` total threads. None below 5 threads, where the worker reading an
/// input decodes it inline. The CRAM inputs share a quarter of the threads
/// first, as their decode is the heaviest stage; the BAM inputs share an eighth
/// from 16 threads, below which the reading worker keeps up with their inflate,
/// out of what the CRAM inputs leave of [`MAX_DECODE_THREADS`].
pub fn decode_workers(threads: usize, formats: &[SniffedFormat]) -> Vec<usize> {
    use SniffedFormat::{Bam, Cram};
    let cram = if threads >= 5 && formats.contains(&Cram) {
        (threads / 4).min(MAX_DECODE_THREADS)
    } else {
        0
    };
    let bam = if threads >= 16 && formats.contains(&Bam) {
        (threads / 8).min(MAX_DECODE_THREADS - cram)
    } else {
        0
    };
    let cram = crate::input::decode_shares(formats, Cram, cram);
    let bam = crate::input::decode_shares(formats, Bam, bam);
    cram.iter().zip(bam).map(|(cram, bam)| cram + bam).collect()
}

/// The worker pool size for `--threads N` once `decode_threads` of them
/// inflate or decode input; at least one worker.
pub fn pool_workers(threads: usize, decode_threads: usize) -> usize {
    threads.saturating_sub(decode_threads).max(1)
}

/// The static per-run engine context shared across fragments: the parsed plan,
/// the compiled groups (in dependency order), the canonical `name -> index` map
/// (for resolving `@grp` anchor names to positional slots), the routing table,
/// the read group of a record with no sample, and the destinations.
struct Engine<'a> {
    plan: &'a DemuxPlan,
    groups: &'a [CompiledGroup],
    group_index: &'a HashMap<String, usize>,
    /// The lookup of every `@group::COLUMN` extract stream, by stream name.
    column_streams: &'a HashMap<&'a str, ColumnStream<'a>>,
    /// Every tag that has a column stream and a default quality tag, with that
    /// quality tag.
    column_tag_qualities: &'a [(&'a str, &'a str)],
    routing: &'a Routing,
    /// The read group id of a record with no sample (see [`pool_read_group`]).
    pool_read_group: &'a str,
    /// Extract-stream names referenced exactly once across all `--template`s
    /// and `--tag`s. A single-stream template whose stream is movable is the
    /// sole consumer of that extract, so the assembled body can MOVE the
    /// extracted bases/quals rather than copy them (the common case for
    /// whole-record template records, which are large; avoids a per-record copy
    /// + alloc).
    movable_streams: &'a HashSet<&'a str>,
    /// Extract-stream names drawn on by any `--template`/`--tag`, the
    /// denominator basis for the `frac_bases_unextracted` metric.
    /// Plan-invariant, so it is built once here rather than rebuilt for every
    /// retained record in [`frac_bases`].
    routed_streams: &'a HashSet<&'a str>,
    /// Every destination with its format and header, shared read-only so a
    /// worker encodes SAM/BAM records with the exact header their output was
    /// opened with, and where each kind of write goes.
    destinations: &'a DestinationTable,
    /// Whether this is an `@RG`-split run (`--sample-from-group` over an `@RG`
    /// source). When set, an assigned record keeps its original `RG:Z` (its
    /// `@RG` line is carried verbatim by the header, see [`build_constituent_rgs`])
    /// rather than being rewritten to a synthesized per-target read group.
    is_read_group_split: bool,
}

/// Run the demux engine: match tag groups per record, route each record with
/// the compiled [`Routing`], fan it out to per-sample `--out` (carrying `--tag`
/// fields and an `@RG` read group) or to the `--remove` / `--unassigned` bins,
/// tally metrics, and emit the metrics TSVs / stderr log.
fn run_engine(plan: &DemuxPlan, command_line: Option<&str>) -> Result<()> {
    let reference = plan
        .reference
        .as_deref()
        .map(|path| {
            crate::reference::Reference::open(path, !plan.skip_reference_check, plan.threads)
        })
        .transpose()?;
    let mut reader = FragmentReader::open_with_threads(
        &plan.inputs,
        plan.per_record,
        reference.as_ref(),
        |formats| decode_workers(plan.threads.max(1), formats),
    )?;
    let input_formats = reader.formats().to_vec();
    if plan.reference.is_some() && !input_formats.contains(&SniffedFormat::Cram) {
        log::warn!("--reference is unused: no input is CRAM");
    }
    let decode_threads = reader.decode_threads();
    // Only read (and conflict-check) the input `@RG` header when an `@RG`
    // group is actually present: `input_read_groups` bails on a repeated RG id
    // with conflicting SM/LB across inputs, which is a legitimate outcome for
    // a non-`@RG` run that never consults SM/LB. Gating keeps that path
    // behavior-identical to before this feature existed.
    let has_read_group_source = plan
        .groups
        .iter()
        .any(|g| matches!(g.source, GroupSource::ReadGroup(_)));
    let rg_info = if has_read_group_source {
        reader.input_read_groups()?
    } else {
        crate::input::ReadGroupInfo::default()
    };
    let resolved = resolve_read_group_groups(plan, &rg_info)?;

    let (tag_sets, tag_columns) = load_tag_sets(plan, &resolved)?;
    let (groups, group_index) = compile_groups(plan, &tag_sets, &resolved)?;
    let column_streams = resolve_column_streams(plan, &tag_columns, &group_index);
    let column_tag_qualities = column_tag_qualities(plan);
    let pool = pool_id(plan);
    let (mut routing, sheet_columns) =
        build_routing(plan, &tag_sets, &group_index, &pool, &resolved)?;
    for index in 0..groups.len() {
        if let Some(warning) = wide_window_warning(index, &groups, plan, &routing) {
            log::warn!("{warning}");
        }
    }
    attach_sample_columns(plan, &mut routing, &sheet_columns)?;

    // An `@RG`-split run (`--sample-from-group` over an `@RG` source) is faithful
    // to the input header: each output file lists its constituent original `@RG`
    // lines and each assigned record keeps its original `RG:Z`. The
    // `Target -> constituent @RG maps` table is statically known before reading
    // any record, so it is built once here.
    let is_read_group_split =
        matches!(&plan.samples, SampleSpec::FromGroup(g) if resolved.contains_key(g));
    let constituent = if is_read_group_split {
        Some(build_constituent_rgs(plan, &resolved, &rg_info))
    } else {
        None
    };

    let alignment_formats = alignment_formats(plan, &routing, &pool, &input_formats)?;
    let read_group = pool_read_group(
        plan,
        &pool,
        writes_alignments(plan, &routing, &pool, &input_formats)?,
    )?;
    let keep_alignments = passes_alignments_through(
        plan,
        is_read_group_split && rg_info.has_sq,
        &alignment_formats,
    )?;
    if keep_alignments {
        refuse_mixed_outputs(plan, &routing, &pool, &input_formats)?;
        log::info!(
            "input is aligned (@SQ present); records pass through unchanged with their alignments"
        );
        reader.keep_alignments();
    }
    let input_header = if keep_alignments {
        reader.alignment_header().cloned()
    } else {
        None
    };

    if is_read_group_split {
        if let Some(pattern) = &plan.out {
            if pattern.uses(Placeholder::Sample) || pattern.uses(Placeholder::SubSample) {
                let mut by_path: HashMap<PathBuf, Vec<String>> = HashMap::new();
                for target in routing.targets() {
                    for path in out_paths(pattern, &pool, Some(target), 1, &input_formats)? {
                        by_path.entry(path).or_default().push(target.label());
                    }
                }
                for (path, labels) in &by_path {
                    if labels.len() > 1 {
                        log::warn!(
                            "targets {} merge into one file `{}`; add a placeholder such as \
                             %sub_sample to separate them",
                            labels.join(", "),
                            path.display()
                        );
                    }
                }
            }
        }
    }

    let headers = out_headers(
        plan,
        &routing,
        &pool,
        &input_formats,
        command_line,
        constituent.as_ref(),
        input_header.as_ref(),
    )?;
    // The fallback header for any destination absent from `headers`
    // (pass-through `--out`, the raw `--unassigned`/`--remove` bins): `@PG`
    // provenance plus a default `@RG` whose ID/SM/LB are the pool id, so every
    // SAM/BAM/CRAM record carries a read group even with no `--sample`. A
    // worker encodes records with the header their output opens with, so the
    // bytes stay identical.
    let default_header = match &input_header {
        Some(input) => alignment_header(input, None, command_line)?,
        None => default_read_group_header(read_group, &plan.rg_tags, command_line)?,
    };
    // Declared first so it drops after every writer, removing what a failed run wrote.
    let mut staged = StagedOutputs::default();
    // Every directed output file, created up front so a directed file always
    // exists, even empty.
    let dests = enumerate_dests(
        plan,
        &routing,
        &pool,
        &input_formats,
        reader.fragment_width(),
    )?;
    check_outputs(
        plan,
        &routing,
        &pool,
        &input_formats,
        reader.fragment_width(),
    )?;
    let body_count = match plan.templates.len() {
        0 => reader.fragment_width(),
        templates => templates,
    };
    let mut list = DestinationList {
        list: Vec::new(),
        index: HashMap::new(),
        input_formats: &input_formats,
        headers: &headers,
        default_header: &default_header,
        keeps_alignments: keep_alignments,
    };
    let destinations = destination_table(plan, &routing, &pool, read_group, body_count, &mut list)?;
    log_run_shape(plan, &input_formats, &dests);
    if let Some(reference) = &reference {
        reader.check_reference(&plan.inputs, reference)?;
    }
    let metrics_paths: Vec<PathBuf> = [
        &plan.metrics_per_sample,
        &plan.metrics_summary,
        &plan.metrics_unassigned,
    ]
    .into_iter()
    .flatten()
    .map(|pattern| resolve_metrics_path(pattern, &pool))
    .collect();
    staged.stage(dests.iter().chain(&metrics_paths).cloned());
    let mut sinks = Sinks::new(
        Arc::clone(&destinations.list),
        plan.compression,
        staged.temporaries(),
    );
    crate::sink::compressor(plan.compression)?;

    // Defensive max(1): clap already guarantees --threads >= 1, but guard for
    // library callers.
    let threads = plan.threads.max(1);
    let workers = pool_workers(threads, decode_threads);
    log::info!(
        "threads={threads}: {workers} worker(s) read, match, and compress; {decode_threads} decode \
         BAM or CRAM input"
    );

    // Streams referenced exactly once across templates + tags can be MOVED into
    // their (sole) output body instead of copied (see Engine::movable_streams).
    let mut stream_refs: HashMap<&str, usize> = HashMap::new();
    for stream in plan
        .templates
        .iter()
        .flat_map(|t| t.streams.iter())
        .chain(plan.tags.iter().flat_map(|t| t.all_streams()))
    {
        *stream_refs.entry(stream.name.as_str()).or_insert(0) += 1;
    }
    let movable_streams: HashSet<&str> = stream_refs
        .iter()
        .filter(|(_, &count)| count == 1)
        .map(|(&name, _)| name)
        .collect();
    // Every stream any template/tag draws on (the keys of `stream_refs`): the
    // plan-invariant basis for `frac_bases`, hoisted here so it is not rebuilt
    // per retained read.
    let routed_streams: HashSet<&str> = stream_refs.keys().copied().collect();

    let engine = Engine {
        plan,
        groups: &groups,
        group_index: &group_index,
        column_streams: &column_streams,
        column_tag_qualities: &column_tag_qualities,
        routing: &routing,
        pool_read_group: read_group,
        movable_streams: &movable_streams,
        routed_streams: &routed_streams,
        destinations: &destinations,
        is_read_group_split,
    };
    let metrics = Metrics::new(pool.clone(), routing.targets());

    // Create the output files and run every stage on the worker pool. The
    // outputs are closed whatever the outcome, and the run error comes first.
    let opened = (|| -> Result<()> {
        sinks.open_files()?;
        for path in &metrics_paths {
            let temporary = staged.path(path);
            if temporary != path.as_path() {
                crate::writer::create_output_file(temporary)?;
            }
        }
        crate::writer::refuse_shared_files(&staged, dests.iter().chain(&metrics_paths))
    })();
    let run = DemuxRun {
        engine: &engine,
        source: Source::new(reader),
        files: sinks.files(),
        merger: Mutex::new(Merger {
            sinks,
            metrics,
            progress: Progress::new(),
            dropped: PerBaseTags::default(),
        }),
        targets: routing.targets().len(),
        compression: plan.compression,
    };
    let run_result = opened.and_then(|()| pool::run(&run, workers));
    let Merger {
        sinks,
        metrics,
        progress,
        dropped,
    } = run
        .merger
        .into_inner()
        .unwrap_or_else(PoisonError::into_inner);
    drop(progress);
    let closed = sinks.close();
    run_result?;
    closed?;

    // Provenance lives in the header written up front: `@PG`
    // (program/version/command line) and the per-sample `@RG`. The demux stats
    // go to stderr and the `--metrics-*` TSVs below (they are not embedded back
    // into the output headers, which would force a full post-pass rewrite of
    // each file).
    metrics.log();
    if is_read_group_split && plan.unassigned.is_none() && metrics.unassigned() > 0 {
        log::warn!(
            "{} records carried no RG:Z, or one the header does not declare, and were dropped; \
             keep them with --unassigned",
            metrics.unassigned().separate_with_commas()
        );
    }
    if !dropped.is_empty() {
        log::warn!(
            "--template changed the bases of some records, so their per-base tags were left out: {}",
            dropped.names()
        );
    }
    if let Some(pattern) = &plan.metrics_per_sample {
        metrics.write_per_sample_tsv(staged.path(&resolve_metrics_path(pattern, &pool)))?;
    }
    if let Some(pattern) = &plan.metrics_summary {
        metrics.write_summary_tsv(staged.path(&resolve_metrics_path(pattern, &pool)))?;
    }
    if let Some(pattern) = &plan.metrics_unassigned {
        metrics
            .write_unassigned_tsv(staged.path(&resolve_metrics_path(pattern, &pool)), &groups)?;
    }
    staged.commit()
}

/// The result of matching every group against one record: either a
/// find-constraint failure (the record is unassigned) or the per-group hits,
/// indexed positionally by canonical group index, that drive routing and body
/// assembly. Owned (one `Vec<Option<GroupHit>>`, no per-record heap maps or
/// `String` clones).
enum FragmentMatch {
    /// A group's find constraint (e.g. an unmet `minFindsPerGroup`) was not
    /// met; carries the reason for the optional `--qc-tag` slug, plus the
    /// partial hits accumulated before the failure (the upstream groups that
    /// matched), so the slug can report what did match.
    Unassigned {
        /// The unassigned reason (the offending group's index + the specific
        /// find failure).
        reason: Unassigned,
        /// The hits filled before the failing group; the failing group and
        /// later slots are `None`.
        hits: Vec<Option<GroupHit>>,
    },
    /// Every group matched or was skipped; the routing inputs are ready.
    Routed {
        /// One slot per group in canonical order; `None` means the group did
        /// not match.
        hits: Vec<Option<GroupHit>>,
    },
}

/// One fully-prepared output record: owned bytes ready for the writer, with no
/// borrows from the engine or the input fragment, so it is `Send` and can be
/// built on a worker and handed to merging. `tags` is `None` when the record
/// carries no SAM data (an empty tag set emits no field).
struct PreparedRead {
    name: Vec<u8>,
    bases: Vec<u8>,
    quals: Option<Vec<u8>>,
    tags: Option<Data>,
    read_group: Option<String>,
    qc_fail: bool,
    alignment: Option<Alignment>,
}

impl PreparedRead {
    /// Borrow this prepared record as the writer's [`OutputRead`], just before the write.
    fn as_output(&self) -> OutputRead<'_> {
        OutputRead {
            name: &self.name,
            bases: &self.bases,
            quals: self.quals.as_deref(),
            tags: self.tags.as_ref(),
            read_group: self.read_group.as_deref(),
            qc_fail: self.qc_fail,
            alignment: self.alignment.as_ref(),
        }
    }
}

/// What merging writes for one [`PreparedWrite`]: bytes already encoded for
/// the destination's format on the worker (appended verbatim, the common path),
/// or the structured records for a format that cannot be pre-encoded to an
/// independent per-fragment byte slice (CRAM).
enum WritePayload {
    /// Pre-encoded record bytes (from [`encode_fragment`]) for FASTX, SAM, and
    /// BAM.
    Encoded(Vec<u8>),
    /// Structured records, encoded as they merge: CRAM, which is not
    /// pre-encodable per fragment.
    Records(Vec<PreparedRead>),
}

/// One write merging must issue: the index of its destination and its
/// payload, in the order the record's writes are issued.
struct PreparedWrite {
    dest: usize,
    payload: WritePayload,
}

/// The metrics tally for one fragment, owned so it crosses the worker
/// boundary. Commutative integer sums, added in input order (order is
/// irrelevant to the sums).
enum Tally {
    Unassigned(Unassigned),
    Removed,
    PassThrough {
        denom: u64,
        unext: u64,
    },
    Assigned {
        target: usize,
        denom: u64,
        unext: u64,
    },
    NonPrimary,
}

/// The fully-prepared outcome of one fragment: the per-record base lengths (for
/// the record-length metric), the metric tally, and the writes to issue in
/// order. Owned + `Send`, so the whole route/extract/assemble/encode pipeline
/// runs on a worker, leaving merging to only tally and issue the (inherently
/// serial, input-ordered) writes.
struct Prepared {
    lengths: SmallVec<[usize; 4]>,
    tally: Tally,
    writes: Vec<PreparedWrite>,
    /// The per-base tags its `--template` bodies left out.
    dropped: PerBaseTags,
}

/// The result of assembling a routed record's body: the writes to issue with the
/// per-base tags its `--template` bodies left out, or the name of a `--template`
/// stream whose anchoring group did not match (so the record is rerouted to the
/// `--unassigned` bin).
enum BodyOutcome {
    Assembled(Vec<PreparedWrite>, PerBaseTags),
    MissingStream(String),
}

/// Match every group against one record in dependency order: the pure,
/// parallelizable half of per-record work. No writer, no metrics; it borrows
/// the record and the shared engine and returns an owned [`FragmentMatch`], so
/// it runs on a worker with a thread-local [`Scratch`].
fn match_fragment(
    engine: &Engine,
    scratch: &mut Scratch,
    fragment: &Fragment,
) -> Result<FragmentMatch> {
    let base_segments: SmallVec<[&[u8]; 4]> = fragment
        .records
        .iter()
        .map(|r| r.bases.as_slice())
        .collect();
    let segments: SmallVec<[Segment; 4]> = fragment
        .records
        .iter()
        .map(|record| Segment {
            bases: &record.bases,
            quals: record.quals.as_deref(),
        })
        .collect();
    // One positional slot per group, in canonical order; `None` until (and
    // unless) the group hits. The presence of `Some` IS the match predicate,
    // preserving today's `contains_key`/`get`->`None` semantics for routing,
    // anchoring, and prerequisite checks.
    let mut hits: Vec<Option<GroupHit>> = vec![None; engine.groups.len()];
    // `groups` is in dependency order, so an upstream group's hit (and its
    // span) is already recorded when its downstream group is matched (a
    // relative `next` link windows the downstream search; a bare `next`/`prev`
    // link skips the downstream group when the upstream did not match).
    for (idx, group) in engine.groups.iter().enumerate() {
        if let Some(rg_to_idx) = &group.read_group {
            if let Some(tag_idx) = match_read_group(rg_to_idx, fragment) {
                hits[idx] = Some(GroupHit {
                    tag_idx,
                    span: None,
                    subs: 0,
                    indels: 0,
                    revcomp: false,
                });
            }
            continue;
        }
        let verdict = match &group.match_streams {
            // A `match=` group matches the joined `--extract` streams, not a
            // read window.
            Some(streams) => match_joined(engine, group, scratch, streams, &segments, &hits)?,
            None => match_group_for(group, scratch, &base_segments, &hits),
        };
        // A failed find constraint (e.g. an unmet minFindsPerGroup) leaves the
        // read unassigned; carry the specific failure and the partial hits so
        // far for the optional QC slug.
        if let Some(failure) = verdict.failure {
            return Ok(FragmentMatch::Unassigned {
                reason: Unassigned::FindConstraint {
                    group: idx,
                    failure,
                },
                hits,
            });
        }
        if let Some(best) = &verdict.best {
            // A `match=` group matches a synthetic joined text, not a read
            // span, so it anchors nothing (the grammar already forbids `@grp`
            // on a match= group): its slot records the tag but no span, which
            // the @grp consumer treats like a missing anchor.
            let span = group.match_streams.is_none().then_some(MatchSpan {
                file: best.file,
                start: best.start,
                end: best.end,
            });
            hits[idx] = Some(GroupHit {
                tag_idx: best.tag_idx,
                span,
                subs: best.subs,
                indels: best.indels,
                revcomp: best.strand == MatchStrand::ReverseComplement,
            });
        }
    }

    Ok(FragmentMatch::Routed { hits })
}

/// The full per-record worker stage: match every group, then
/// route/extract/assemble/encode the record into an owned, `Send` [`Prepared`].
/// This is everything except the metric tally and the inherently-serial,
/// input-ordered writes, so the whole pipeline runs on a worker with a
/// thread-local [`Scratch`].
fn process_fragment(
    engine: &Engine,
    scratch: &mut Scratch,
    fragment: &Fragment,
) -> Result<Prepared> {
    let fmatch = match_fragment(engine, scratch, fragment)?;
    let mut prepared = prepare(engine, fragment, fmatch)?;
    if is_non_primary(fragment) {
        prepared.tally = Tally::NonPrimary;
        let destinations = &engine.destinations.list;
        prepared
            .writes
            .retain(|write| is_alignment_format(destinations[write.dest].format));
    }
    encode_writes(engine, &mut prepared.writes)?;
    Ok(prepared)
}

/// Whether a fragment is a kept secondary or supplementary alignment: it goes
/// only to SAM/BAM destinations and is left out of the metrics.
fn is_non_primary(fragment: &Fragment) -> bool {
    fragment.records.iter().any(|record| {
        record.alignment.as_ref().is_some_and(|alignment| {
            let flags = alignment.flags();
            flags.is_secondary() || flags.is_supplementary()
        })
    })
}

/// Whether `format` is SAM, BAM, or CRAM.
fn is_alignment_format(format: OutputFormat) -> bool {
    matches!(
        format,
        OutputFormat::Sam | OutputFormat::Bam | OutputFormat::Cram
    )
}

/// Encode each write's records into the destination's on-disk bytes here on the
/// worker, so serial merging only appends them in input order rather than
/// running the per-record SAM/BAM/FASTX encode itself. CRAM stays structured
/// (its records are container/slice-coded, not independent per-fragment byte
/// slices) and is encoded as it merges. The format and header are the
/// destination's, which its output opens with, so the appended bytes match an
/// encode at merging, except that a raw BAM record keeps its input bytes.
fn encode_writes(engine: &Engine, writes: &mut [PreparedWrite]) -> Result<()> {
    for write in writes {
        let destination = &engine.destinations.list[write.dest];
        // CRAM is container/slice-coded, so it cannot be pre-encoded per fragment.
        if matches!(destination.format, OutputFormat::Cram) {
            continue;
        }
        let bytes = {
            let WritePayload::Records(reads) = &write.payload else {
                continue;
            };
            let outputs: SmallVec<[OutputRead; 2]> =
                reads.iter().map(PreparedRead::as_output).collect();
            encode_fragment(destination.format, &destination.header, &outputs)?
        };
        write.payload = WritePayload::Encoded(bytes);
    }
    Ok(())
}

/// Prepare one matched record into an owned, `Send` [`Prepared`] (the
/// per-record lengths, the metric tally, and the writes to issue): route it,
/// compute `frac_bases`, and assemble its body / bins, but touch neither the
/// writer nor the metrics, so it runs on a worker (alongside matching) while
/// only the tally and the ordered writes stay serial. The assigned/pass-through
/// decision reflects the post-assembly outcome (a missing `--template` stream
/// reroutes to unassigned).
fn prepare(engine: &Engine, fragment: &Fragment, fmatch: FragmentMatch) -> Result<Prepared> {
    let lengths: SmallVec<[usize; 4]> = fragment.records.iter().map(|r| r.bases.len()).collect();
    let hits = match fmatch {
        FragmentMatch::Unassigned { reason, hits } => {
            let writes = prepare_unassigned(engine, fragment, &reason, &hits)?;
            return Ok(Prepared {
                lengths,
                tally: Tally::Unassigned(reason),
                writes,
                dropped: PerBaseTags::default(),
            });
        }
        FragmentMatch::Routed { hits } => hits,
    };
    let segments: SmallVec<[Segment; 4]> = fragment
        .records
        .iter()
        .map(|record| Segment {
            bases: &record.bases,
            quals: record.quals.as_deref(),
        })
        .collect();

    let mut dropped = PerBaseTags::default();
    let (tally, writes) = match engine.routing.route(&hits) {
        Disposition::Removed(remove) => (
            Tally::Removed,
            prepare_removed(engine, fragment, remove, &hits)?,
        ),
        Disposition::Unassigned => unassigned(engine, fragment, Unassigned::NoSample, &hits)?,
        Disposition::PassThrough => {
            let (denom, unext) = frac_bases(engine, &segments, &hits)?;
            match prepare_body(engine, None, &hits, &segments, fragment)? {
                BodyOutcome::Assembled(writes, left_out) => {
                    dropped = left_out;
                    (Tally::PassThrough { denom, unext }, writes)
                }
                BodyOutcome::MissingStream(missing) => {
                    unassigned(engine, fragment, Unassigned::MissingStream(missing), &hits)?
                }
            }
        }
        Disposition::Assigned(route) => {
            let (denom, unext) = frac_bases(engine, &segments, &hits)?;
            match prepare_body(engine, Some(route), &hits, &segments, fragment)? {
                BodyOutcome::Assembled(writes, left_out) => {
                    dropped = left_out;
                    (
                        Tally::Assigned {
                            target: route.target_index(),
                            denom,
                            unext,
                        },
                        writes,
                    )
                }
                BodyOutcome::MissingStream(missing) => {
                    unassigned(engine, fragment, Unassigned::MissingStream(missing), &hits)?
                }
            }
        }
    };
    Ok(Prepared {
        lengths,
        tally,
        writes,
        dropped,
    })
}

/// The tally and writes of a record routed to the `--unassigned` bin.
fn unassigned(
    engine: &Engine,
    fragment: &Fragment,
    reason: Unassigned,
    hits: &GroupHits,
) -> Result<(Tally, Vec<PreparedWrite>)> {
    let writes = prepare_unassigned(engine, fragment, &reason, hits)?;
    Ok((Tally::Unassigned(reason), writes))
}

/// Count one [`Prepared`] fragment in `metrics`.
fn tally(metrics: &mut Metrics, prepared: &Prepared) {
    if !matches!(prepared.tally, Tally::NonPrimary) {
        metrics.record_processed();
        metrics.record_read_lengths(prepared.lengths.iter().copied());
    }
    match &prepared.tally {
        Tally::Unassigned(reason) => metrics.record_unassigned(reason),
        Tally::Removed => metrics.record_removed(),
        Tally::PassThrough { denom, unext } => metrics.record_pass_through(*denom, *unext),
        Tally::Assigned {
            target,
            denom,
            unext,
        } => metrics.record_assigned(*target, *denom, *unext),
        Tally::NonPrimary => metrics.record_non_primary(),
    }
}

/// The most fragments in one batch, the unit of work a worker prepares at once.
const READ_CHUNK: usize = 2048;

/// The most bytes of records in one batch, so a batch of long records stays small.
const BATCH_BYTES: usize = 4 << 20;

/// Records between progress log lines; a final total is logged on completion.
const PROGRESS_UNIT: u64 = 500_000;

/// A running record counter that logs a comma-formatted total every
/// [`PROGRESS_UNIT`] records and a final total on drop (so completion and the
/// error path both report). Logging through this module's own `log` target
/// scopes the line as `unmux::demux`, like the rest of the engine's logging.
struct Progress {
    seen: u64,
}

impl Progress {
    fn new() -> Self {
        Self { seen: 0 }
    }

    /// Count `records` more, logging the running total at each
    /// [`PROGRESS_UNIT`] boundary they pass.
    fn add(&mut self, records: u64) {
        let passed = self.seen / PROGRESS_UNIT;
        self.seen += records;
        for unit in passed + 1..=self.seen / PROGRESS_UNIT {
            let total = unit * PROGRESS_UNIT;
            log::info!("Demultiplexed {} records", total.separate_with_commas());
        }
    }
}

impl Drop for Progress {
    /// Log the final total unless the last record already landed on a boundary
    /// (which would duplicate the line); a zero count stays silent.
    fn drop(&mut self) {
        if !self.seen.is_multiple_of(PROGRESS_UNIT) {
            log::info!("Demultiplexed {} records", self.seen.separate_with_commas());
        }
    }
}

/// The next batch of fragments from `reader`, or `None` at the end of the
/// input. A read error propagates.
fn next_batch(reader: &mut FragmentReader) -> Result<Option<Vec<PendingFragment>>> {
    let mut batch = Vec::new();
    let mut bytes = 0;
    while batch.len() < READ_CHUNK && bytes < BATCH_BYTES {
        match reader.next_pending_fragment()? {
            Some(fragment) => {
                bytes += fragment_bytes(&fragment);
                batch.push(fragment);
            }
            None => break,
        }
    }
    Ok((!batch.is_empty()).then_some(batch))
}

/// About how many bytes the records of a fragment hold.
fn fragment_bytes(fragment: &PendingFragment) -> usize {
    let size = |record: &PendingRecord| match record {
        PendingRecord::Ready(record) => record.name.len() + 2 * record.bases.len(),
        PendingRecord::Bam { raw, .. } => raw.len(),
    };
    fragment.records.iter().map(size).sum()
}

/// Decode a fragment's raw records, reusing `record`, then run [`process_fragment`].
fn process_pending(
    engine: &Engine,
    scratch: &mut Scratch,
    record: &mut RecordBuf,
    pending: PendingFragment,
) -> Result<Prepared> {
    let fragment = pending.decode(record)?;
    process_fragment(engine, scratch, &fragment)
}

/// A batch of fragments prepared on a worker: what merging writes and counts.
struct Processed {
    /// The fragments prepared without error.
    fragments: u64,
    /// Their metrics.
    tally: Metrics,
    /// Their encoded bytes, one buffer per destination.
    bytes: Vec<(usize, Vec<u8>)>,
    /// Their records for destinations that encode them as they merge, in order.
    records: Vec<(usize, Vec<PreparedRead>)>,
    /// The per-base tags their `--template` bodies left out.
    dropped: PerBaseTags,
    /// The error that stopped the batch after the fragments above.
    error: Option<anyhow::Error>,
}

impl Processed {
    /// Add one prepared fragment, finding the buffer of each destination
    /// through `slots`, by destination the index of its buffer in `bytes`.
    fn add(&mut self, prepared: Prepared, slots: &mut [usize]) {
        tally(&mut self.tally, &prepared);
        self.dropped |= prepared.dropped;
        for write in prepared.writes {
            match write.payload {
                WritePayload::Encoded(bytes) => {
                    let slot = &mut slots[write.dest];
                    if *slot == usize::MAX {
                        *slot = self.bytes.len();
                        self.bytes.push((write.dest, Vec::new()));
                    }
                    self.bytes[*slot].1.extend_from_slice(&bytes);
                }
                WritePayload::Records(reads) => self.records.push((write.dest, reads)),
            }
        }
        self.fragments += 1;
    }
}

/// The state of merging: the outputs, the metrics, the record count, and the
/// per-base tags left out of `--template` bodies.
struct Merger {
    sinks: Sinks,
    metrics: Metrics,
    progress: Progress,
    dropped: PerBaseTags,
}

/// What each worker keeps across the jobs it runs.
struct DemuxWorker {
    scratch: Scratch,
    record: RecordBuf,
    fastq: noodles::fastq::Record,
    compressor: Option<::bgzf::Compressor>,
    slots: Vec<usize>,
}

/// The fragments read while opening FASTQ inputs, and the error that ended
/// them, if one did.
type Head = (Vec<PendingFragment>, Option<anyhow::Error>);

/// Where the batches of a run come from.
enum Source {
    /// Every input read a fragment at a time, by one worker at a time.
    Fragments {
        reader: Mutex<(FragmentReader, u64)>,
        exhausted: AtomicBool,
    },
    /// FASTQ inputs read in chunks of records, each input by one worker at a
    /// time, after the fragments read while opening them, which are batch 0.
    Fastq {
        source: FastqSource,
        head: Mutex<Option<Head>>,
        interleaved: bool,
        /// The input paths, when several inputs may not hold both mates of a pair.
        pairs: Option<Vec<PathBuf>>,
    },
}

impl Source {
    /// Read every input of `reader` in chunks when they are all FASTQ.
    fn new(reader: FragmentReader) -> Self {
        match reader.into_fastq() {
            Ok(inputs) => {
                let step = if inputs.interleaved { 2 } else { 1 };
                let source = FastqSource::new(inputs.streams, step, READ_CHUNK * step, BATCH_BYTES);
                Source::Fastq {
                    source,
                    head: Mutex::new(Some(inputs.head)),
                    interleaved: inputs.interleaved,
                    pairs: inputs.pairs,
                }
            }
            Err(reader) => Source::Fragments {
                reader: Mutex::new((*reader, 0)),
                exhausted: AtomicBool::new(false),
            },
        }
    }
}

/// A batch of fragments, or the chunks of FASTQ records that hold them.
enum Batch {
    Fragments(Vec<PendingFragment>),
    Chunks(Vec<Chunk>),
}

/// A demux run on the worker pool: inputs read in batches, workers that
/// prepare and encode the fragments of each, and outputs merged in input order.
struct DemuxRun<'a> {
    engine: &'a Engine<'a>,
    source: Source,
    merger: Mutex<Merger>,
    files: Arc<[OnceLock<BlockFile>]>,
    targets: usize,
    compression: u8,
}

impl DemuxRun<'_> {
    fn merger(&self) -> std::sync::MutexGuard<'_, Merger> {
        self.merger.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Prepare `pending` into `processed`, or keep its error there; whether it
    /// was prepared.
    fn prepare(
        &self,
        worker: &mut DemuxWorker,
        pending: Result<PendingFragment>,
        processed: &mut Processed,
    ) -> bool {
        let (scratch, record) = (&mut worker.scratch, &mut worker.record);
        match pending.and_then(|pending| process_pending(self.engine, scratch, record, pending)) {
            Ok(prepared) => {
                processed.add(prepared, &mut worker.slots);
                true
            }
            Err(error) => {
                processed.error = Some(error);
                false
            }
        }
    }
}

impl pool::Stages for DemuxRun<'_> {
    type Batch = Batch;
    type Done = Processed;
    type Block = Block;
    type Worker = DemuxWorker;

    fn streams(&self) -> usize {
        match &self.source {
            Source::Fragments { .. } => 1,
            Source::Fastq { source, .. } => source.inputs(),
        }
    }

    fn worker(&self) -> DemuxWorker {
        DemuxWorker {
            scratch: Scratch::new(),
            record: RecordBuf::default(),
            fastq: noodles::fastq::Record::default(),
            compressor: None,
            slots: vec![usize::MAX; self.files.len()],
        }
    }

    fn can_read(&self, stream: usize) -> bool {
        match &self.source {
            Source::Fragments { exhausted, .. } => !exhausted.load(Ordering::Acquire),
            Source::Fastq { source, .. } => source.can_read(stream),
        }
    }

    fn read(&self, _: &mut DemuxWorker, stream: usize) -> Readout<Batch> {
        let mut readout = Readout::nothing();
        match &self.source {
            Source::Fragments { reader, exhausted } => {
                let mut reader = reader.lock().unwrap_or_else(PoisonError::into_inner);
                let (fragments, next) = &mut *reader;
                let index = *next;
                match next_batch(fragments) {
                    Ok(Some(batch)) => {
                        *next += 1;
                        readout.batches.push((index, Batch::Fragments(batch)));
                    }
                    Ok(None) => readout.end = Some(End::Exhausted(index)),
                    Err(error) => readout.end = Some(End::Failed(index, error)),
                }
                if readout.end.is_some() {
                    exhausted.store(true, Ordering::Release);
                }
            }
            Source::Fastq { source, head, .. } => {
                let head = head.lock().unwrap_or_else(PoisonError::into_inner).take();
                if let Some((fragments, error)) = head {
                    readout.batches.push((0, Batch::Fragments(fragments)));
                    if let Some(error) = error {
                        source.stop();
                        readout.end = Some(End::Failed(1, error));
                        return readout;
                    }
                }
                let framed = source.read(stream);
                let chunks = framed.chunks.into_iter();
                let batches = chunks.map(|(index, chunks)| (index + 1, Batch::Chunks(chunks)));
                readout.batches.extend(batches);
                readout.end = framed.total.map(|total| End::Exhausted(total + 1));
            }
        }
        readout
    }

    fn process(&self, worker: &mut DemuxWorker, batch: Batch) -> Processed {
        let mut processed = Processed {
            fragments: 0,
            tally: Metrics::tally(self.targets),
            bytes: Vec::new(),
            records: Vec::new(),
            dropped: PerBaseTags::default(),
            error: None,
        };
        match batch {
            Batch::Fragments(fragments) => {
                for pending in fragments {
                    if !self.prepare(worker, Ok(pending), &mut processed) {
                        break;
                    }
                }
            }
            Batch::Chunks(chunks) => {
                let (interleaved, pairs) = match &self.source {
                    Source::Fastq {
                        interleaved, pairs, ..
                    } => (*interleaved, pairs.as_deref()),
                    Source::Fragments { .. } => (false, None),
                };
                let mut fragments = ChunkFragments::new(chunks, interleaved, pairs);
                loop {
                    let pending = match fragments.next(&mut worker.fastq) {
                        Ok(Some(pending)) => Ok(pending),
                        Ok(None) => break,
                        Err(error) => Err(error),
                    };
                    if !self.prepare(worker, pending, &mut processed) {
                        break;
                    }
                }
            }
        }
        for &(dest, _) in &processed.bytes {
            worker.slots[dest] = usize::MAX;
        }
        processed
    }

    fn merge(&self, done: Processed, blocks: &mut Vec<Block>) -> Result<()> {
        let mut merger = self.merger();
        let Merger {
            sinks,
            metrics,
            progress,
            dropped,
        } = &mut *merger;
        metrics.absorb(&done.tally);
        *dropped |= done.dropped;
        for (dest, bytes) in &done.bytes {
            sinks.write(*dest, bytes, blocks)?;
        }
        for (dest, reads) in &done.records {
            let outputs: SmallVec<[OutputRead; 2]> =
                reads.iter().map(PreparedRead::as_output).collect();
            sinks.write_records(*dest, &outputs)?;
        }
        progress.add(done.fragments);
        done.error.map_or(Ok(()), Err)
    }

    fn finish(&self, blocks: &mut Vec<Block>) -> Result<()> {
        self.merger().sinks.finish(blocks);
        Ok(())
    }

    fn compress(&self, worker: &mut DemuxWorker, block: Block) -> Result<()> {
        let compressor = match &mut worker.compressor {
            Some(compressor) => compressor,
            empty => empty.insert(crate::sink::compressor(self.compression)?),
        };
        crate::sink::compress(compressor, &self.files, block)
    }
}

/// The `frac_bases_unextracted` numerator and denominator for one retained
/// record: the denominator is the total bases of input files a
/// `--template`/`--tag` stream draws from, and the numerator is the bases of
/// those files not covered by any such stream's source range (a union per file,
/// so overlaps are not double-counted). A run with no `--template` is a
/// pass-through whose whole record is the output, so nothing is unextracted.
///
/// Known limitation: a routed `@a..@b` between-anchors stream that resolves to
/// an empty stream (an anchor unmatched, or the anchors bounding no region)
/// contributes no range, so its file is not added to the denominator. This is a
/// rare diagnostic-metric edge; file spans and matched `@grp` anchors (the
/// common cases) are exact.
fn frac_bases(engine: &Engine, segments: &[Segment], hits: &GroupHits) -> Result<(u64, u64)> {
    let plan = engine.plan;
    if plan.templates.is_empty() {
        let total: u64 = segments.iter().map(|s| s.bases.len() as u64).sum();
        return Ok((total, 0));
    }
    let routed = engine.routed_streams;
    let mut covered: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for spec in &plan.extracts {
        if routed.contains(spec.name.as_str()) {
            if let Some((file, start, end)) =
                span_range(&spec.span, segments, hits, engine.group_index)?
            {
                if end > start {
                    covered.entry(file).or_default().push((start, end));
                }
            }
        }
    }
    let mut denominator = 0u64;
    let mut unextracted = 0u64;
    for (&file, intervals) in &covered {
        let file_len = segments.get(file).map_or(0, |s| s.bases.len());
        denominator += file_len as u64;
        unextracted += file_len.saturating_sub(union_len(intervals)) as u64;
    }
    Ok((denominator, unextracted))
}

/// The total length covered by a set of `[start, end)` intervals, merging
/// overlaps (so a base covered by two streams counts once).
fn union_len(intervals: &[(usize, usize)]) -> usize {
    let mut sorted = intervals.to_vec();
    sorted.sort_unstable();
    let mut covered = 0;
    let mut cursor = 0;
    for (start, end) in sorted {
        let start = start.max(cursor);
        if end > start {
            covered += end - start;
            cursor = end;
        }
    }
    covered
}

/// Route a fragment by its carried `RG:Z` for an `@RG` group: the tag index of
/// the read group the record belongs to, or `None` (no RG tag, or an id absent
/// from the header).
fn match_read_group(rg_to_idx: &HashMap<Vec<u8>, usize>, fragment: &Fragment) -> Option<usize> {
    for record in &fragment.records {
        if let Some(data) = &record.tags {
            if let Some(Value::String(rg)) = data.get(&Tag::READ_GROUP) {
                return rg_to_idx.get(AsRef::<[u8]>::as_ref(rg)).copied();
            }
        }
    }
    None
}

/// Match a `match=` group against the concatenation of its named `--extract`
/// streams (no delimiter), each as [`join_bases`] takes it. A referenced stream
/// that is absent (its anchoring group did not match, or it resolved to
/// nothing) makes the group a non-match. The tags are searched within the
/// joined text as a single synthetic segment, so
/// `dist`/`mode`/`delta`/`bothStrands` apply over the concatenation.
fn match_joined(
    engine: &Engine,
    group: &CompiledGroup,
    scratch: &mut Scratch,
    streams: &[String],
    segments: &[Segment],
    hits: &GroupHits,
) -> Result<GroupVerdict> {
    let mut joined = Vec::new();
    for name in streams {
        let Some(spec) = engine.plan.extracts.iter().find(|e| &e.name == name) else {
            return Ok(GroupVerdict::unmatched(group));
        };
        match extract_stream(engine, spec, segments, hits)? {
            Some(extracted) => joined.extend_from_slice(join_bases(&extracted)),
            None => return Ok(GroupVerdict::unmatched(group)),
        }
    }
    if joined.is_empty() {
        return Ok(GroupVerdict::unmatched(group));
    }
    Ok(match_verdict(group, scratch, &[joined.as_slice()], None))
}

/// A stream's bases in a `match=` join: its corrected form unless the tag has an ambiguity code.
fn join_bases(extracted: &Extracted) -> &[u8] {
    match &extracted.corrected {
        Some(corrected) if crate::iupac::degeneracy(&corrected.bases) == 1 => &corrected.bases,
        _ => &extracted.bases,
    }
}

/// Match one group, honoring its sequential prerequisite: a group with no
/// prerequisite searches at its own `loc`; a downstream group whose upstream
/// did not match is skipped (unmatched); a relative `next=GROUP:lo-hi` link
/// windows the search past the upstream match, in that match's orientation.
fn match_group_for(
    group: &CompiledGroup,
    scratch: &mut Scratch,
    segments: &[&[u8]],
    hits: &GroupHits,
) -> GroupVerdict {
    let Some(prereq) = &group.prereq else {
        return match_verdict(group, scratch, segments, group.loc.as_ref());
    };
    // Skip the group when any prerequisite upstream did not match (its slot is
    // `None`). Upstreams are earlier in dependency order, so their slots are
    // already filled.
    if prereq
        .upstreams
        .iter()
        .any(|&upstream| hits[upstream].is_none())
    {
        return GroupVerdict::unmatched(group);
    }
    // A relative `next` link windows the search relative to its upstream's
    // match (present, since the upstream is one of the now-satisfied upstreams;
    // an upstream is always a non-`match=` group, so its span is `Some`);
    // otherwise the group searches at its own `loc`.
    match &prereq.window {
        Some(window) => {
            let (span, reverse) = hits[window.upstream]
                .as_ref()
                .and_then(|hit| Some((hit.span?, hit.revcomp)))
                .expect("a relative-window upstream is a matched non-match= group with a span");
            match_verdict_after(group, scratch, segments, &span, reverse, window)
        }
        None => match_verdict(group, scratch, segments, group.loc.as_ref()),
    }
}

/// Route a removed record to its `--remove` bin (raw input segments), or drop
/// it when the rule has no pattern. With `--qc-tag`, each written segment
/// carries the `removed` slug naming the rule that fired.
fn prepare_removed(
    engine: &Engine,
    fragment: &Fragment,
    remove: &RemoveTarget,
    hits: &GroupHits,
) -> Result<Vec<PreparedWrite>> {
    let rule = engine
        .routing
        .removes()
        .iter()
        .position(|rule| std::ptr::eq(rule, remove));
    let layout = rule.and_then(|rule| engine.destinations.removes[rule].as_ref());
    let Some(layout) = layout else {
        return Ok(Vec::new());
    };
    let slug = engine.plan.qc_tag.as_ref().map(|tag| {
        let base_segments: Vec<&[u8]> = fragment
            .records
            .iter()
            .map(|r| r.bases.as_slice())
            .collect();
        (
            tag.as_str(),
            removed_slug(
                engine.groups,
                hits,
                &base_segments,
                &remove.selector,
                remove.group,
                remove.tag,
            ),
        )
    });
    prepare_raw_bin(
        layout,
        engine.pool_read_group,
        fragment,
        slug.as_ref().map(|(tag, value)| (*tag, value.as_str())),
    )
}

/// Route an unassigned record to the `--unassigned` bin (raw input segments),
/// or drop it when no bin is configured. With `--qc-tag`, each written segment
/// carries the unassigned-reason slug.
fn prepare_unassigned(
    engine: &Engine,
    fragment: &Fragment,
    reason: &Unassigned,
    hits: &GroupHits,
) -> Result<Vec<PreparedWrite>> {
    match &engine.destinations.unassigned {
        Some(layout) => {
            let slug = engine.plan.qc_tag.as_ref().map(|tag| {
                let base_segments: Vec<&[u8]> = fragment
                    .records
                    .iter()
                    .map(|r| r.bases.as_slice())
                    .collect();
                (
                    tag.as_str(),
                    unassigned_slug(engine.groups, hits, &base_segments, reason),
                )
            });
            prepare_raw_bin(
                layout,
                engine.pool_read_group,
                fragment,
                slug.as_ref().map(|(tag, value)| (*tag, value.as_str())),
            )
        }
        None => Ok(Vec::new()),
    }
}

/// Prepare a record's raw input segments for a `--unassigned` / `--remove` bin.
/// With `%source` over multiple `--in` files the segments fan out one file per
/// input file (record `i` came from file `i`); a single input (whether plain or
/// auto-detected interleaved) is one `%source` (0), so its mates stay together.
/// Without `%source`, every segment goes to one file.
fn prepare_raw_bin(
    layout: &BinLayout,
    read_group: &str,
    fragment: &Fragment,
    slug: Option<(&str, &str)>,
) -> Result<Vec<PreparedWrite>> {
    match layout {
        BinLayout::PerSource(dests) => fragment
            .records
            .iter()
            .zip(dests)
            .map(|(record, &dest)| {
                Ok(PreparedWrite {
                    dest,
                    payload: WritePayload::Records(vec![prepared_raw_read(
                        record, read_group, slug,
                    )?]),
                })
            })
            .collect(),
        &BinLayout::Single(dest) => Ok(vec![PreparedWrite {
            dest,
            payload: WritePayload::Records(
                fragment
                    .records
                    .iter()
                    .map(|record| prepared_raw_read(record, read_group, slug))
                    .collect::<Result<_>>()?,
            ),
        }]),
    }
}

/// One raw input segment as an owned prepared record for a
/// `--unassigned`/`--remove` bin: the record's bases/quals/name and its carried
/// input tags, optionally augmented with the `--qc-tag` slug. Carries the
/// default pool read group so a SAM/BAM/CRAM bin matches the default `@RG` in
/// its header; FASTX output ignores the read group.
fn prepared_raw_read(
    record: &InputRecord,
    read_group: &str,
    slug: Option<(&str, &str)>,
) -> Result<PreparedRead> {
    let tags = match slug {
        Some((tag, value)) => {
            let mut data = record.tags.clone().unwrap_or_default();
            insert_string_tag(&mut data, tag, value);
            Some(data)
        }
        None => record.tags.clone(),
    };
    let alignment = match (record.alignment.clone(), slug) {
        (Some(alignment), Some((tag, value))) => {
            let mut record = alignment.into_record_buf()?;
            insert_string_tag(record.data_mut(), tag, value);
            Some(Alignment::Record(record))
        }
        (alignment, _) => alignment,
    };
    Ok(PreparedRead {
        name: record.name.clone(),
        bases: record.bases.clone(),
        quals: record.quals.clone(),
        tags,
        read_group: Some(read_group.to_string()),
        qc_fail: record.qc_fail,
        alignment,
    })
}

/// The canonical (corrected) tag sequence for a matched `@grp` self-span, or
/// `None` when the group did not match. Resolved lazily from the group's
/// positional hit (its matched tag index) against the engine's group slice (the
/// tag's canonical bases), so no per-record correction buffer is carried across
/// the matcher boundary.
fn canonical_correction<'a>(engine: &'a Engine, hits: &GroupHits, group: &str) -> Option<&'a [u8]> {
    let idx = *engine.group_index.get(group)?;
    let compiled = &engine.groups[idx];
    let hit = hits.get(idx)?.as_ref()?;
    compiled.tags.get(hit.tag_idx).map(Vec::as_slice)
}

/// Resolve a stream, carrying the matched tag as the corrected form of an `@grp` self-span.
fn extract_stream(
    engine: &Engine,
    spec: &crate::grammar::Extract,
    segments: &[Segment],
    hits: &GroupHits,
) -> Result<Option<Extracted>> {
    let Some(mut extracted) = extract(&spec.span, segments, hits, engine.group_index)? else {
        return Ok(None);
    };
    if let SpanBody::AnchorMatch { group } = &spec.span.body {
        if let Some(canonical) = canonical_correction(engine, hits, group) {
            extracted.corrected = Some(corrected_form(&extracted, canonical));
        }
    }
    Ok(Some(extracted))
}

/// The matched tag's value for a column stream, or `None` when its group did not
/// match.
fn column_value<'a>(engine: &'a Engine, hits: &GroupHits, stream: &str) -> Option<&'a str> {
    let column = engine.column_streams.get(stream)?;
    let hit = hits.get(column.group)?.as_ref()?;
    column.values.get(hit.tag_idx).map(String::as_str)
}

/// Assemble the output body for an assigned (`route = Some`) or pass-through
/// (`route = None`) record into owned
/// [`PreparedWrite`]s: resolve its extraction streams, build the `--template`
/// bodies and the `--tag` data (leaving out a tag with a blank `%sample::COLUMN`
/// stream), attach the `@RG` read group on SAM/BAM, and resolve the `--out`
/// destination(s). Returns [`BodyOutcome::MissingStream`] when a `--template`
/// references a stream whose group did not match (the caller reroutes it to the
/// `--unassigned` bin and tallies it unassigned). Pure (no writer, no metrics),
/// so it runs on the worker.
fn prepare_body(
    engine: &Engine,
    route: Option<&Route>,
    hits: &GroupHits,
    segments: &[Segment],
    fragment: &Fragment,
) -> Result<BodyOutcome> {
    let plan = engine.plan;
    let target = route.map(Route::target);
    let mut streams: HashMap<&str, Extracted> = HashMap::new();
    let mut slot = 0;
    for spec in &plan.extracts {
        if matches!(spec.span.body, SpanBody::SampleColumn { .. }) {
            let value = route.and_then(|route| route.column(slot));
            slot += 1;
            if let Some(value) = value {
                streams.insert(spec.name.as_str(), Extracted::from_text(value));
            }
            continue;
        }
        if matches!(spec.span.body, SpanBody::Column { .. }) {
            if let Some(value) = column_value(engine, hits, &spec.name) {
                streams.insert(spec.name.as_str(), Extracted::from_text(value));
            }
            continue;
        }
        if let Some(extracted) = extract_stream(engine, spec, segments, hits)? {
            streams.insert(spec.name.as_str(), extracted);
        }
    }

    let sample = target.map(|t| t.sample.as_str());
    let sub_sample = target.and_then(|t| t.sub_sample.as_deref());
    let layout = match route {
        Some(route) => &engine.destinations.targets[route.target_index()],
        None => engine
            .destinations
            .pass_through
            .as_ref()
            .context("internal error: no --out layout for a record without a sample")?,
    };
    let is_alignment = layout.is_alignment;

    let (bodies, dropped) = match assemble_bodies(
        plan,
        fragment,
        &mut streams,
        engine.movable_streams,
        is_alignment,
    ) {
        Ok(assembled) => assembled,
        // A templated stream is absent (its anchoring group did not match): the
        // caller reroutes to unassigned, naming the missing stream in the QC
        // slug.
        Err(missing) => return Ok(BodyOutcome::MissingStream(missing)),
    };

    let tag_data = build_tag_data(&plan.tags, &streams);
    let left_out = left_out_tags(&plan.tags, &tag_data);
    let text_qualities: Vec<(&str, &str)> = engine
        .column_tag_qualities
        .iter()
        .copied()
        .filter(|(tag, _)| {
            plan.tags
                .iter()
                .find(|binding| binding.tag == *tag)
                .is_some_and(|binding| lacks_qualities(binding, &streams))
        })
        .collect();
    // An assigned read carries its sample's read group; a pass-through read (no
    // target) carries the default pool read group, so a SAM/BAM/CRAM `--out`
    // matches the default `@RG` in its header. On an `@RG`-split run the record
    // keeps its carried original `RG:Z` instead (already in `read.tags` via
    // `merge_data`); `None` leaves it untouched, since `writer::build_record`
    // only inserts `RG:Z` when a read group is set.
    let read_group = match engine.is_read_group_split {
        true => None,
        false => layout.read_group.clone(),
    };

    let mut datas: Vec<Data> = bodies
        .iter()
        .map(|body| merge_data(body.carried.as_ref(), &tag_data, &text_qualities, &left_out))
        .collect();
    // With `--qc-tag`, every output body of this record carries the
    // routing-outcome slug (assigned with its sample, or pass-through). Built
    // once from the read's group hits and observed bases.
    let slug = plan.qc_tag.as_ref().map(|tag| {
        let base_segments: Vec<&[u8]> = segments.iter().map(|segment| segment.bases).collect();
        (
            tag.as_str(),
            routed_slug(engine.groups, hits, &base_segments, sample, sub_sample),
        )
    });
    if let Some((tag, slug)) = &slug {
        for data in &mut datas {
            insert_string_tag(data, tag, slug);
        }
    }

    // Multi-read FASTX to a file fans out one file per `%ordinal`; SAM/BAM (a
    // pair in one file), single-read FASTX, and stdout all write every read to
    // one destination.
    let dests = layout
        .dests
        .as_ref()
        .map_err(|message| anyhow!("{message}"))?;
    let writes = match dests {
        OutDests::PerOrdinal(dests) => {
            ensure!(
                dests.len() == bodies.len(),
                "internal error: {} bodies for {} --out ordinals",
                bodies.len(),
                dests.len()
            );
            bodies
                .into_iter()
                .zip(datas)
                .zip(dests)
                .map(|((body, data), &dest)| PreparedWrite {
                    dest,
                    payload: WritePayload::Records(vec![prepared_body_read(body, data, None)]),
                })
                .collect()
        }
        &OutDests::Single(dest) => {
            let reads = bodies
                .into_iter()
                .zip(datas)
                .map(|(mut body, data)| match body.alignment.take() {
                    Some(record) => passed_through_read(
                        record,
                        &tag_data,
                        slug.as_ref().map(|(tag, value)| (*tag, value.as_str())),
                    ),
                    None => {
                        let read_group = if is_alignment {
                            read_group.clone()
                        } else {
                            None
                        };
                        Ok(prepared_body_read(body, data, read_group))
                    }
                })
                .collect::<Result<_>>()?;
            vec![PreparedWrite {
                dest,
                payload: WritePayload::Records(reads),
            }]
        }
    };
    Ok(BodyOutcome::Assembled(writes, dropped))
}

/// One assembled body as an owned prepared record: the body's bytes are moved
/// (no copy), the merged tag `Data` becomes the record's tags (omitted when
/// empty, so an unmatched record leaves a null tag), and `read_group` is the
/// SAM/BAM `@RG` id (already gated to alignment output by the caller).
fn prepared_body_read(body: Body, data: Data, read_group: Option<String>) -> PreparedRead {
    PreparedRead {
        name: body.name,
        bases: body.bases,
        quals: body.quals,
        tags: (!data.is_empty()).then_some(data),
        read_group,
        qc_fail: body.qc_fail,
        alignment: None,
    }
}

/// An aligned record passed through whole, gaining only `--tag` and `--qc-tag`.
/// Without either it stays as read, so raw BAM bytes are written back verbatim.
fn passed_through_read(
    alignment: Alignment,
    tag_data: &Data,
    slug: Option<(&str, &str)>,
) -> Result<PreparedRead> {
    let inserted = tag_data
        .iter()
        .any(|(_, value)| !matches!(value, Value::String(text) if text.is_empty()));
    let alignment = if !inserted && slug.is_none() {
        alignment
    } else {
        let mut record = alignment.into_record_buf()?;
        let data = record.data_mut();
        for (tag, value) in tag_data.iter() {
            if !matches!(value, Value::String(text) if text.is_empty()) {
                data.insert(tag, value.clone());
            }
        }
        if let Some((tag, value)) = slug {
            insert_string_tag(data, tag, value);
        }
        Alignment::Record(record)
    };
    Ok(PreparedRead {
        name: Vec::new(),
        bases: Vec::new(),
        quals: None,
        tags: None,
        read_group: None,
        qc_fail: false,
        alignment: Some(alignment),
    })
}

/// One assembled output record: its name, bases, qualities, and the carried
/// input tags it should keep.
struct Body {
    name: Vec<u8>,
    bases: Vec<u8>,
    quals: Option<Vec<u8>>,
    carried: Option<Data>,
    qc_fail: bool,
    alignment: Option<Alignment>,
}

/// SAM tags tied to a read's exact bases: modifications, per-base qualities, kinetics, moves.
const PER_BASE_TAGS: [[u8; 2]; 25] = [
    *b"MM", *b"ML", *b"MN", *b"OQ", *b"BQ", *b"E2", *b"U2", *b"CQ", *b"CS", *b"FZ", *b"PT", *b"BI",
    *b"BD", *b"ip", *b"pw", *b"fi", *b"fp", *b"ri", *b"rp", *b"dq", *b"dt", *b"iq", *b"mq", *b"sq",
    *b"mv",
];

/// A set of [`PER_BASE_TAGS`], one bit per tag.
#[derive(Clone, Copy, Default)]
struct PerBaseTags(u32);

impl PerBaseTags {
    fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// The tags of the set, comma-separated in [`PER_BASE_TAGS`] order.
    fn names(self) -> String {
        PER_BASE_TAGS
            .iter()
            .enumerate()
            .filter(|&(index, _)| self.0 & 1 << index != 0)
            .map(|(_, tag)| String::from_utf8_lossy(tag))
            .collect::<Vec<_>>()
            .join(", ")
    }
}

impl std::ops::BitOrAssign for PerBaseTags {
    fn bitor_assign(&mut self, other: Self) {
        self.0 |= other.0;
    }
}

/// The tags of `data` without its per-base ones, and which per-base tags it held.
fn without_per_base_tags(data: &Data) -> (Data, PerBaseTags) {
    let mut dropped = PerBaseTags::default();
    let kept = data
        .iter()
        .filter(|(tag, _)| {
            let per_base = PER_BASE_TAGS.iter().position(|bytes| tag.as_ref() == bytes);
            if let Some(index) = per_base {
                dropped.0 |= 1 << index;
            }
            per_base.is_none()
        })
        .map(|(tag, value)| (tag, value.clone()))
        .collect();
    (kept, dropped)
}

/// Assemble a record's output bodies: the raw input segments (each carrying its
/// tags, and its source alignment when `keep_alignment` is set) when there is no
/// `--template`, else one body per `--template`. A template body carries the
/// tags of the read that supplies its first segment, without the per-base ones
/// unless the body is that read unchanged; the per-base tags left out are
/// returned with the bodies. Returns `Err(missing_stream_name)` when a template
/// references a stream that was not extracted (its anchoring group did not
/// match), naming the absent stream for the `--qc-tag` slug.
fn assemble_bodies(
    plan: &DemuxPlan,
    fragment: &Fragment,
    streams: &mut HashMap<&str, Extracted>,
    movable: &HashSet<&str>,
    keep_alignment: bool,
) -> std::result::Result<(Vec<Body>, PerBaseTags), String> {
    if plan.templates.is_empty() {
        let bodies = fragment
            .records
            .iter()
            .map(
                |record| match record.alignment.as_ref().filter(|_| keep_alignment) {
                    Some(alignment) => Body {
                        name: Vec::new(),
                        bases: Vec::new(),
                        quals: None,
                        carried: None,
                        qc_fail: false,
                        alignment: Some(alignment.clone()),
                    },
                    None => Body {
                        name: record.name.clone(),
                        bases: record.bases.clone(),
                        quals: record.quals.clone(),
                        carried: record.tags.clone(),
                        qc_fail: record.qc_fail,
                        alignment: None,
                    },
                },
            )
            .collect();
        return Ok((bodies, PerBaseTags::default()));
    }
    let name = fragment
        .records
        .first()
        .map(|record| record.name.clone())
        .unwrap_or_default();
    let qc_fail = fragment.records.iter().any(|record| record.qc_fail);
    let mut bodies = Vec::with_capacity(plan.templates.len());
    let mut dropped = PerBaseTags::default();
    for template in &plan.templates {
        let mut sources = template
            .streams
            .iter()
            .filter_map(|stream| streams.get(stream.name.as_str())?.source);
        let first = sources.next();
        let joined = sources.any(|source| Some(source) != first);
        // A single-stream template whose stream is referenced nowhere else is
        // that extract's sole consumer, so move its (possibly large)
        // bases/quals into the body instead of copying. The moved stream is in
        // no `--tag`, so the later `build_tag_data` does not need it. A
        // `raw=false` template of a corrected stream emits the corrected form;
        // otherwise the observed bases move.
        let (bases, quals) = match template.streams.as_slice() {
            [only] if movable.contains(only.name.as_str()) => {
                let mut extracted = streams
                    .remove(only.name.as_str())
                    .ok_or_else(|| only.name.clone())?;
                let (mut bases, mut quals) = match extracted.corrected.take() {
                    Some(corrected) if !template.raw => (corrected.bases, corrected.quals),
                    _ => (extracted.bases, extracted.quals),
                };
                // A `~` on the (sole, moved) stream reverse-complements its
                // contribution in place.
                if only.revcomp {
                    crate::iupac::reverse_complement(&mut bases);
                    if let Some(quals) = quals.as_mut() {
                        quals.reverse();
                    }
                }
                (bases, quals)
            }
            _ => assemble_template(template, streams)?,
        };
        let source = match first {
            Some(file) => fragment.records.get(file),
            None => fragment.records.first(),
        };
        let carried = source.and_then(|record| {
            let tags = record.tags.as_ref()?;
            let unchanged = !joined
                && template.streams.iter().all(|stream| !stream.revcomp)
                && bases == record.bases
                && quals == record.quals;
            if unchanged {
                return Some(tags.clone());
            }
            let (kept, left_out) = without_per_base_tags(tags);
            dropped |= left_out;
            Some(kept)
        });
        bodies.push(Body {
            name: name.clone(),
            bases,
            quals,
            carried,
            qc_fail,
            alignment: None,
        });
    }
    Ok((bodies, dropped))
}

/// The `--tag` bindings that `tag_data` leaves out because a stream is absent,
/// each with its explicit `qual=` tag, if any.
fn left_out_tags<'a>(
    bindings: &'a [TagBinding],
    tag_data: &Data,
) -> Vec<(&'a str, Option<&'a str>)> {
    bindings
        .iter()
        .filter(|binding| tag_data.get(&tag_bytes(&binding.tag)).is_none())
        .map(|binding| {
            let quality = match &binding.qual {
                QualTag::Named(name) => Some(name.as_str()),
                _ => None,
            };
            (binding.tag.as_str(), quality)
        })
        .collect()
}

/// Merge the carried input tags (if any) with the assembled `--tag` fields,
/// with `--tag` winning a key collision. `column_qualities` pairs each tag that
/// has a column stream with its default quality tag, which is dropped from the
/// carried tags when that tag is rewritten without qualities. A `left_out` tag
/// is dropped from the carried tags, along with its explicit and default
/// quality tags unless another carried tag that this record does not rewrite
/// keeps one (a `CY` with a `CR`).
fn merge_data(
    carried: Option<&Data>,
    tag_data: &Data,
    column_qualities: &[(&str, &str)],
    left_out: &[(&str, Option<&str>)],
) -> Data {
    let mut data = carried.cloned().unwrap_or_default();
    for (tag, quality) in column_qualities {
        if is_stale_quality(&data, tag_data, tag, quality) {
            let quality = Tag::from(tag_bytes(quality));
            data = data
                .iter()
                .filter(|(tag, _)| *tag != quality)
                .map(|(tag, value)| (tag, value.clone()))
                .collect();
        }
    }
    let is_left_out = |tag: &[u8; 2]| left_out.iter().any(|(name, _)| tag_bytes(name) == *tag);
    let keeps_quality = |quality: &str| {
        data.iter().any(|(tag, _)| {
            let kept = tag_data.get(&tag).is_none();
            let tag = <[u8; 2]>::from(tag);
            kept && !is_left_out(&tag)
                && std::str::from_utf8(&tag).is_ok_and(|tag| default_qual_tag(tag) == Some(quality))
        })
    };
    let dropped: Vec<[u8; 2]> = left_out
        .iter()
        .flat_map(|(tag, quality)| {
            let qualities = quality.iter().copied().chain(default_qual_tag(tag));
            let qualities = qualities.filter(|quality| !keeps_quality(quality));
            std::iter::once(*tag).chain(qualities).map(tag_bytes)
        })
        .collect();
    if !dropped.is_empty() {
        data = data
            .iter()
            .filter(|(tag, _)| !dropped.contains(&<[u8; 2]>::from(*tag)))
            .map(|(tag, value)| (tag, value.clone()))
            .collect();
    }
    for (tag, value) in tag_data.iter() {
        data.insert(tag, value.clone());
    }
    data
}

/// Whether the carried `quality` tag described the old value of `tag`, which
/// `tag_data` rewrites without qualities. A `CY` is the quality of `CR`, so one
/// carried with a `CR` that is kept outlives a rewritten `CB`.
fn is_stale_quality(carried: &Data, tag_data: &Data, tag: &str, quality: &str) -> bool {
    let has = |data: &Data, tag: &str| data.get(&tag_bytes(tag)).is_some();
    let keeps_cr = tag == "CB" && has(carried, "CR") && !has(tag_data, "CR");
    has(carried, quality) && has(tag_data, tag) && !has(tag_data, quality) && !keeps_cr
}

/// Build the corrected form of a matched `@grp` stream: always the declared
/// canonical sequence, in the orientation the user supplied it (the read strand
/// is a use-site choice via `~`, not baked in here). Degenerate / IUPAC
/// positions are emitted verbatim. The observed qualities, already in the
/// tag's orientation, are copied and resized to the corrected length.
/// Non-destructive: the observed bases/quals on `extracted` are left intact for
/// the record body and for a `raw=true` tag/template.
fn corrected_form(extracted: &Extracted, canonical: &[u8]) -> Corrected {
    let bases = canonical.to_vec();
    let quals = extracted
        .quals
        .as_ref()
        .map(|quals| fit_quals(quals.clone(), bases.len()));
    Corrected { bases, quals }
}

/// Resize observed qualities to `len`: identical for a same-length
/// (substitution) correction; truncated when longer; padded with the last
/// observed score (or a neutral 30) when shorter.
fn fit_quals(mut quals: Vec<u8>, len: usize) -> Vec<u8> {
    let pad = quals.last().copied().unwrap_or(30);
    quals.resize(len, pad);
    quals
}

/// Concatenate a template's streams into one output record's bases and
/// qualities; qualities are kept only when every stream carries them. Returns
/// `Err(stream_name)` when a referenced stream is absent (its anchoring group
/// did not match), which makes the record unassigned.
fn assemble_template(
    template: &Template,
    streams: &HashMap<&str, Extracted>,
) -> std::result::Result<(Vec<u8>, Option<Vec<u8>>), String> {
    let mut bases = Vec::new();
    let mut quals = Some(Vec::new());
    for stream in &template.streams {
        let extracted = streams
            .get(stream.name.as_str())
            .ok_or_else(|| stream.name.clone())?;
        let stream_bases = extracted.tag_bases(template.raw);
        let stream_quals = extracted.tag_quals(template.raw);
        if stream.revcomp {
            // Reverse-complement this stream's contribution at the use site.
            let mut rc = stream_bases.to_vec();
            crate::iupac::reverse_complement(&mut rc);
            bases.extend_from_slice(&rc);
            match (quals.as_mut(), stream_quals) {
                (Some(accumulated), Some(stream_quals)) => {
                    accumulated.extend(stream_quals.iter().rev());
                }
                _ => quals = None,
            }
        } else {
            bases.extend_from_slice(stream_bases);
            match (quals.as_mut(), stream_quals) {
                (Some(accumulated), Some(stream_quals)) => {
                    accumulated.extend_from_slice(stream_quals)
                }
                _ => quals = None,
            }
        }
    }
    Ok((bases, quals))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use std::path::Path;

    /// A minimal [`DemuxArgs`] with everything off, overridden per test.
    fn base_args() -> DemuxArgs {
        DemuxArgs {
            pool: None,
            inputs: Vec::new(),
            groups: Vec::new(),
            extracts: Vec::new(),
            templates: Vec::new(),
            tags: Vec::new(),
            rg_tags: Vec::new(),
            samples: Vec::new(),
            sample_sheet: None,
            sample_sheet_ignore: Vec::new(),
            sample_from_group: None,
            require_samples_explain_all_tags: false,
            remove: Vec::new(),
            out: None,
            unassigned: None,
            metrics_per_sample: None,
            metrics_summary: None,
            metrics_unassigned: None,
            qc_tag: None,
            per_record: false,
            compression: 5,
            threads: 8,
            reference: None,
            skip_reference_check: false,
            command_line: None,
        }
    }

    #[test]
    fn test_decode_workers_fit_inside_the_thread_budget() {
        use SniffedFormat::{Bam, Cram, Fastq};
        let mixed = [&[Bam, Cram][..], &[Cram, Bam], &[Bam, Bam, Cram]];
        let single = [&[Bam][..], &[Bam, Bam], &[Cram], &[Fastq, Bam], &[Fastq]];
        for formats in single.into_iter().chain(mixed) {
            for threads in 1..=40 {
                let decode: usize = decode_workers(threads, formats).iter().sum();
                assert!(decode <= MAX_DECODE_THREADS);
                let workers = pool_workers(threads, decode);
                assert!(workers >= 1);
                assert_eq!(
                    workers + decode,
                    threads,
                    "{formats:?} at {threads} threads"
                );
            }
        }
        assert_eq!(decode_workers(8, &[Bam]), [0]);
        assert_eq!(decode_workers(16, &[Bam]), [2]);
        assert_eq!(decode_workers(8, &[Cram]), [2]);
        assert_eq!(decode_workers(4, &[Cram]), [0]);
        assert_eq!(decode_workers(64, &[Fastq]), [0]);
    }

    #[test]
    fn test_decode_workers_go_to_cram_inputs_first() {
        use SniffedFormat::{Bam, Cram};
        for threads in 5..=7 {
            assert_eq!(decode_workers(threads, &[Bam, Cram]), [0, 1]);
        }
        assert_eq!(decode_workers(8, &[Bam, Cram]), [0, 2]);
        assert_eq!(decode_workers(12, &[Cram, Bam]), [3, 0]);
        assert_eq!(decode_workers(16, &[Bam, Cram]), [0, 4]);
        assert_eq!(decode_workers(32, &[Bam, Bam, Cram]), [0, 0, 4]);
        assert_eq!(decode_workers(32, &[Bam, Bam]), [2, 2]);
    }

    #[test]
    fn test_pool_workers_keep_one_worker() {
        assert_eq!(pool_workers(1, 0), 1);
        assert_eq!(pool_workers(8, 0), 8);
        assert_eq!(pool_workers(8, 2), 6);
        assert_eq!(pool_workers(2, 4), 1);
    }

    #[test]
    fn test_low_thread_output_matches_high_thread() {
        // --threads 1 and 2 must produce the same
        // demultiplexed records as a high-thread run. Plain FASTQ is
        // byte-identical across every thread count. Compressed output is
        // compared as decoded records via read_back rather than raw bytes.
        let dir = tempfile::tempdir().unwrap();
        let n = READ_CHUNK + 137; // spans a chunk boundary plus a partial tail
        let mut fq = Vec::new();
        for i in 0..n {
            let bc = match i % 3 {
                0 => "AAAAAAAA",
                1 => "TTTTTTTT",
                _ => "GGGGGGGG",
            };
            fq.extend_from_slice(format!("@r{i}\n{bc}CGCG\n+\nIIIIIIIIIIII\n").as_bytes());
        }
        let input = write_file(dir.path(), "in.fq", &fq);

        let run = |threads: usize, out: &Path| {
            let mut args = base_args();
            args.threads = threads;
            args.inputs = vec![format!("0={}", input.display())];
            args.groups = vec![
                "grp={D1=AAAAAAAA,D2=TTTTTTTT}".to_string(),
                "grp::loc=0:0,dist=0".to_string(),
            ];
            args.out = Some(out.display().to_string());
            run_demux(args).unwrap();
        };

        // Plain FASTQ: byte-identical across all thread counts.
        let ref_fq = dir.path().join("ref.fq");
        run(16, &ref_fq);
        for threads in [1usize, 2, 3] {
            let out = dir.path().join(format!("t{threads}.fq"));
            run(threads, &out);
            assert_eq!(
                std::fs::read(&out).unwrap(),
                std::fs::read(&ref_fq).unwrap(),
                "plain FASTQ at --threads {threads} must match the 16-thread run byte-for-byte"
            );
        }

        // Compressed BAM: compare decoded records.
        let ref_bam = dir.path().join("ref.bam");
        run(16, &ref_bam);
        let ref_records = read_back(&ref_bam);
        assert!(
            !ref_records.is_empty(),
            "the reference run produced records"
        );
        for threads in [1usize, 2, 3] {
            let out = dir.path().join(format!("t{threads}.bam"));
            run(threads, &out);
            assert_eq!(
                read_back(&out),
                ref_records,
                "BAM records at --threads {threads} must match the 16-thread run"
            );
        }
    }

    /// A BAM file's decompressed header bytes and its raw records, each with its
    /// `block_size` prefix.
    fn raw_bam(path: &Path) -> (Vec<u8>, Vec<Vec<u8>>) {
        use std::io::Read;
        let mut bytes = Vec::new();
        noodles::bgzf::io::Reader::new(File::open(path).unwrap())
            .read_to_end(&mut bytes)
            .unwrap();
        let int = |at: usize| i32::from_le_bytes(bytes[at..at + 4].try_into().unwrap()) as usize;
        let mut at = 8 + int(4);
        let references = int(at);
        at += 4;
        for _ in 0..references {
            at += 4 + int(at) + 4;
        }
        let header = bytes[..at].to_vec();
        let mut records = Vec::new();
        while at < bytes.len() {
            let end = at + 4 + int(at);
            records.push(bytes[at..end].to_vec());
            at = end;
        }
        (header, records)
    }

    #[test]
    fn test_bam_input_output_is_identical_across_threads() {
        use noodles::sam::alignment::io::Write as _;
        let dir = tempfile::tempdir().unwrap();
        let mut sam = String::from(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000000\n@SQ\tSN:chr2\tLN:10000000\n\
             @RG\tID:rg0\tSM:s0\n@RG\tID:rg1\tSM:s1\n@RG\tID:rg2\tSM:s2\n",
        );
        let n = READ_CHUNK + 4321;
        for i in 0..n {
            let chrom = if i < n / 2 { "chr1" } else { "chr2" };
            let pos = 1 + (i % (n / 2)) * 3;
            let flag = match i % 7 {
                0 => 99,
                1 => 147,
                2 => 256,
                3 => 2048,
                4 => 16,
                5 => 73,
                _ => 0,
            };
            let len = 20 + i % 3;
            let seq: String = (0..len)
                .map(|j| ["A", "C", "G", "T", "N"][(i + j) % 5])
                .collect();
            let qual: String = (0..len)
                .map(|j| char::from(b'!' + ((i + j) % 40) as u8))
                .collect();
            sam.push_str(&format!(
                "r{i}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{len}M\t=\t{pos}\t0\t{seq}\t{qual}\t\
                 RG:Z:rg{rg}\tNM:i:{nm}\tMC:Z:{len}M\tXB:B:s,1,-2,{xb}\n",
                mapq = i % 61,
                rg = i % 3,
                xb = i % 30000,
                nm = i % 4,
            ));
        }
        let mut reader = noodles::sam::io::Reader::new(sam.as_bytes());
        let header = reader.read_header().unwrap();
        let input = dir.path().join("in.bam");
        let mut writer = noodles::bam::io::Writer::new(File::create(&input).unwrap());
        writer.write_header(&header).unwrap();
        let mut record = RecordBuf::default();
        while reader.read_record_buf(&header, &mut record).unwrap() > 0 {
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish().unwrap();
        let (_, input_records) = raw_bam(&input);

        let run = |threads: usize, ext: &str| -> PathBuf {
            let out = dir.path().join(format!("t{threads}-{ext}"));
            let mut args = base_args();
            args.threads = threads;
            args.inputs = vec![format!("0={}", input.display())];
            args.groups = vec!["rg=@RG".to_string()];
            args.sample_from_group = Some("rg".to_string());
            args.out = Some(format!("{}/%sample.{ext}", out.display()));
            // FASTQ output is not a pass-through, and these paired records have no mates.
            args.per_record = ext == "fq";
            run_demux(args).unwrap();
            out
        };

        let reference = run(1, "bam");
        for rg in 0..3 {
            let (_, records) = raw_bam(&reference.join(format!("rg{rg}.bam")));
            let expected: Vec<Vec<u8>> = input_records
                .iter()
                .enumerate()
                .filter(|(i, _)| i % 3 == rg)
                .map(|(_, raw)| raw.clone())
                .collect();
            assert_eq!(
                records, expected,
                "rg{rg} holds its input records byte for byte"
            );
        }
        for threads in [2usize, 4, 8] {
            let out = run(threads, "bam");
            for rg in 0..3 {
                let name = format!("rg{rg}.bam");
                assert_eq!(
                    raw_bam(&out.join(&name)),
                    raw_bam(&reference.join(&name)),
                    "{name} at --threads {threads} matches --threads 1"
                );
            }
        }

        let reference = run(1, "fq");
        for threads in [2usize, 4, 8] {
            let out = run(threads, "fq");
            for rg in 0..3 {
                let name = format!("rg{rg}.fq");
                assert_eq!(
                    std::fs::read(out.join(&name)).unwrap(),
                    std::fs::read(reference.join(&name)).unwrap(),
                    "{name} at --threads {threads} matches --threads 1"
                );
            }
        }
    }

    /// Write `bytes` to `dir/name` and return the path.
    fn write_file(dir: &Path, name: &str, bytes: &[u8]) -> PathBuf {
        let path = dir.join(name);
        let mut file = File::create(&path).unwrap();
        file.write_all(bytes).unwrap();
        file.flush().unwrap();
        path
    }

    /// Read every record back out of an output file (sniffed), for round-trip
    /// assertions. A directed output file now always exists; a zero-record one
    /// is an empty (0-byte) file, which has no format to sniff, so it reads
    /// back as no records.
    fn read_back(path: &Path) -> Vec<crate::input::InputRecord> {
        if std::fs::metadata(path)
            .map(|m| m.len() == 0)
            .unwrap_or(true)
        {
            return Vec::new();
        }
        let mut reader = FragmentReader::open(&[path.to_path_buf()], false).unwrap();
        let mut records = Vec::new();
        while let Some(fragment) = reader.next_fragment().unwrap() {
            records.extend(fragment.records);
        }
        records
    }

    /// Read the first record out of a BAM, for SAM-tag (corrected/raw)
    /// assertions.
    fn first_bam_record(path: &Path) -> noodles::sam::alignment::RecordBuf {
        let mut reader = noodles::bam::io::Reader::new(File::open(path).unwrap());
        let header = reader.read_header().unwrap();
        reader.record_bufs(&header).next().unwrap().unwrap()
    }

    /// The string value of a two-character SAM tag on a record, or `None` if
    /// absent/not a string.
    fn string_tag(record: &noodles::sam::alignment::RecordBuf, tag: [u8; 2]) -> Option<String> {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        match record.data().get(&Tag::new(tag[0], tag[1])) {
            Some(Value::String(value)) => Some(value.to_string()),
            _ => None,
        }
    }

    #[test]
    fn test_passthrough_fastq_to_fastq() {
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nTTTT\n+\nJJJJ\n",
        );
        let out = dir.path().join("out.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let records = read_back(&out);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, b"r1");
        assert_eq!(records[0].bases, b"ACGTACGT");
        // 'I' (raw 40) round-trips ASCII -> raw -> ASCII -> raw.
        assert_eq!(records[0].quals.as_deref(), Some(&[40u8; 8][..]));
        assert_eq!(records[1].bases, b"TTTT");
    }

    #[test]
    fn test_passthrough_fastq_to_bam() {
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r1\nACGTACGT\n+\nIIIIIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let records = read_back(&out);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].bases, b"ACGTACGT");
        assert_eq!(records[0].quals.as_deref(), Some(&[40u8; 8][..]));
    }

    #[test]
    fn test_carve_with_extract_and_template() {
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r1\nAAACCCGGGTTT\n+\n012345678901\n");
        let out = dir.path().join("out.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.extracts = vec!["body=0:9:end".to_string()];
        args.templates = vec!["body".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let records = read_back(&out);
        assert_eq!(records.len(), 1);
        // read[9:end] of the 12 bp read keeps the last 3 bases.
        assert_eq!(records[0].bases, b"TTT");
    }

    #[test]
    fn test_tag_emits_into_fastq_comment() {
        // --tag BC=bc emits the barcode (and its default QT quality tag) into
        // the FASTQ comment.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r1\nACGTACGT\n+\nIIIIIIII\n");
        let out = dir.path().join("out.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.extracts = vec!["bc=0:0:4".to_string()];
        args.tags = vec!["BC=bc".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.contains("BC:Z:ACGT"), "BC carried in comment: {text}");
        assert!(
            text.contains("QT:Z:"),
            "default QT quality tag present: {text}"
        );
    }

    #[test]
    fn test_group_match_drives_anchored_extraction() {
        // A matched group anchors `@g`-relative extraction; the read after the
        // barcode is carved out.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@r1\nTTTTACGTACGTGGGGGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n",
        );
        let out = dir.path().join("out.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec!["g={ACGTACGT}".to_string(), "g::loc=0:0".to_string()];
        args.extracts = vec!["rest=@g+0:end".to_string()];
        args.templates = vec!["rest".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let records = read_back(&out);
        assert_eq!(records.len(), 1);
        // The barcode matches at [4, 12); `@g+0:end` keeps the bases after it.
        assert_eq!(records[0].bases, b"GGGGGGGG");
    }

    #[test]
    fn test_required_group_routes_assigned_and_unassigned() {
        // minFindsPerGroup=1 makes the group required: a matching read goes to
        // --out, a non-matching read goes to --unassigned.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@m\nACGTACGT\n+\nIIIIIIII\n@u\nTTTTTTTT\n+\nIIIIIIII\n",
        );
        let out = dir.path().join("out.fq");
        let unassigned = dir.path().join("un.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec!["g={ACGT}".to_string(), "g::minFindsPerGroup=1".to_string()];
        args.out = Some(out.display().to_string());
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        let assigned = read_back(&out);
        assert_eq!(assigned.len(), 1);
        assert_eq!(assigned[0].name, b"m");
        let unmatched = read_back(&unassigned);
        assert_eq!(unmatched.len(), 1);
        assert_eq!(unmatched[0].name, b"u");
    }

    #[test]
    fn test_qc_tag_emits_routing_slug_for_assigned_and_unassigned() {
        // --qc-tag ZS writes a per-record JSON demux-provenance slug. An
        // assigned read carries its sample, matched group,
        // observed-vs-corrected bases, and edit breakdown; an unassigned read
        // carries the failure reason and offending group. Off by default; here
        // opted in via qc_tag.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            // r_ok: observed GAAGGG (1 mismatch from canonical GAAGAG) matches
            // and routes to dna01. r_no: TTTTTT matches nothing, fails
            // minFindsPerGroup=1, lands unassigned.
            b"@r_ok\nGAAGGGAAAA\n+\nIIIIIIIIII\n@r_no\nTTTTTTAAAA\n+\nIIIIIIIIII\n",
        );
        let out = dir.path().join("out.bam");
        let unassigned = dir.path().join("un.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={dna01=GAAGAG}".to_string(),
            "grp::loc=0:0,dist=1,minFindsPerGroup=1".to_string(),
        ];
        args.samples = vec!["dna01=grp::dna01".to_string()];
        args.out = Some(out.display().to_string());
        args.unassigned = Some(unassigned.display().to_string());
        args.qc_tag = Some("ZS".to_string());
        run_demux(args).unwrap();

        let assigned = first_bam_record(&out);
        assert_eq!(
            string_tag(&assigned, *b"ZS").as_deref(),
            Some(
                r#"{"v":1,"outcome":"assigned","sample":"dna01","sub_sample":null,"groups":[{"g":"grp","tag":"GAAGAG","loc":"0:0:6","obs":"GAAGGG","sub":1,"ind":0}]}"#
            )
        );

        let unmatched = first_bam_record(&unassigned);
        assert_eq!(
            string_tag(&unmatched, *b"ZS").as_deref(),
            Some(
                r#"{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"grp","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}"#
            )
        );
    }

    #[test]
    fn test_bare_next_chain_searches_own_loc() {
        // grp_a::next=grp_b is an ordering link only: grp_b (required) still
        // searches at its own loc. A read carrying both passes; one missing
        // grp_b's tag fails grp_b and is unassigned.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@ok\nACGTGGGG\n+\nIIIIIIII\n@no\nACGTAAAA\n+\nIIIIIIII\n",
        );
        let out = dir.path().join("out.fq");
        let unassigned = dir.path().join("un.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp_a={ACGT}".to_string(),
            "grp_a::loc=0:0:4,next=grp_b".to_string(),
            "grp_b={GGGG}".to_string(),
            "grp_b::loc=0:4:8,minFindsPerGroup=1".to_string(),
        ];
        args.out = Some(out.display().to_string());
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        assert_eq!(read_back(&out)[0].name, b"ok");
        assert_eq!(read_back(&unassigned)[0].name, b"no");
    }

    #[test]
    fn test_prev_prerequisite_skips_downstream() {
        // grp_y::prev=grp_x: when grp_x does not match, grp_y is skipped (not
        // searched) even though the read carries grp_y's tag, so the read is
        // unassigned (feature: prev prerequisite).
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@ok\nACGTGGGG\n+\nIIIIIIII\n@skip\nTTTTGGGG\n+\nIIIIIIII\n",
        );
        let out = dir.path().join("out.fq");
        let unassigned = dir.path().join("un.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp_x={ACGT}".to_string(),
            "grp_x::loc=0:0:4".to_string(),
            "grp_y={GGGG}".to_string(),
            "grp_y::loc=0:4:8,prev=grp_x,minFindsPerGroup=1".to_string(),
        ];
        args.out = Some(out.display().to_string());
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        // `ok` has grp_x then grp_y; `skip` has grp_y's GGGG but no grp_x, so
        // grp_y is never searched.
        assert_eq!(read_back(&out)[0].name, b"ok");
        let unmatched = read_back(&unassigned);
        assert_eq!(unmatched.len(), 1);
        assert_eq!(unmatched[0].name, b"skip");
    }

    #[test]
    fn test_relative_next_windows_downstream() {
        // grp_a::next=grp_b:2-6 windows grp_b at [grp_a.end+2, grp_a.end+6].
        // grp_b (no loc of its own) matches only when its tag sits in that
        // window.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@in\nAAAACCGGGG\n+\nIIIIIIIIII\n@out\nAAAAGGGGCC\n+\nIIIIIIIIII\n",
        );
        let out = dir.path().join("out.fq");
        let unassigned = dir.path().join("un.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp_a={AAAA}".to_string(),
            "grp_a::loc=0:0:4,next=grp_b:2-6".to_string(),
            "grp_b={GGGG}".to_string(),
            "grp_b::minFindsPerGroup=1".to_string(),
        ];
        args.out = Some(out.display().to_string());
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        // `in`: AAAA at [0,4], window [6,10] holds GGGG -> matches. `out`: GGGG
        // at [4,8] is outside the [6,10] window -> grp_b fails -> unassigned.
        assert_eq!(read_back(&out)[0].name, b"in");
        assert_eq!(read_back(&unassigned)[0].name, b"out");
    }

    #[test]
    fn test_match_over_joined_streams_demuxes_dual_index() {
        // the sample barcode is the i7+i5 concatenation. `match=i7+i5` matches
        // the joined 16 bp value against the group, which no per-read window
        // could find (each input is 8 bp).
        let dir = tempfile::tempdir().unwrap();
        let i7 = write_file(
            dir.path(),
            "i7.fq",
            b"@p0\nAAAAAAAA\n+\nIIIIIIII\n@p1\nAAAAAAAA\n+\nIIIIIIII\n",
        );
        let i5 = write_file(
            dir.path(),
            "i5.fq",
            b"@p0\nCCCCCCCC\n+\nIIIIIIII\n@p1\nGGGGGGGG\n+\nIIIIIIII\n",
        );
        let unassigned = dir.path().join("un.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", i7.display()), format!("1={}", i5.display())];
        args.extracts = vec![
            "i7=0:0:8".to_string(),
            "i5=1:0:8".to_string(),
            "body=0:0:end".to_string(),
        ];
        args.templates = vec!["body".to_string()];
        args.groups = vec![
            "sample_bc={s1=AAAAAAAACCCCCCCC}".to_string(),
            "sample_bc::match=i7+i5,dist=0".to_string(),
        ];
        args.sample_from_group = Some("sample_bc".to_string());
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        // p0's joined barcode equals s1 -> routed to sample s1; p1's
        // (AAAA..GGGG..) matches nothing.
        let s1 = read_back(&dir.path().join("s1.fq"));
        assert_eq!(s1.len(), 1);
        assert_eq!(s1[0].name, b"p0");
        let unmatched = read_back(&unassigned);
        assert!(unmatched.iter().any(|r| r.name == b"p1"), "p1 unassigned");
    }

    #[test]
    fn test_match_group_ordered_after_anchor_declared_later() {
        // A `match=` group whose stream anchors on a group declared AFTER it.
        // The match join only exists once the anchor group has matched, so the
        // anchor must be ordered first. Without the match=-anchored-extract
        // prerequisite, declaration order would leave the anchor unmatched at
        // join time and the read would route to unassigned instead of its
        // sample.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nAAAACCCCGGGG\n+\nIIIIIIIIIIII\n");
        let unassigned = dir.path().join("un.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        // `sample_bc` is declared first, but it matches the joined value of
        // `umi`, which is anchored on `grpZ` (declared second). The
        // prerequisite must reorder `grpZ` ahead of `sample_bc`.
        args.groups = vec![
            "sample_bc={s1=GGGG}".to_string(),
            "sample_bc::match=umi,dist=0".to_string(),
            "grpZ={GGGG}".to_string(),
            "grpZ::loc=0:8:12".to_string(),
        ];
        args.extracts = vec!["umi=@grpZ".to_string(), "body=0:0:end".to_string()];
        args.templates = vec!["body".to_string()];
        args.sample_from_group = Some("sample_bc".to_string());
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        // grpZ matches GGGG, umi resolves to GGGG, sample_bc's join equals s1
        // -> routed to s1.
        let s1 = read_back(&dir.path().join("s1.fq"));
        assert_eq!(s1.len(), 1, "the read routed to sample s1");
        assert_eq!(s1[0].name, b"r");
        assert!(read_back(&unassigned).is_empty(), "nothing unassigned");
    }

    #[test]
    fn test_match_joins_the_corrected_form_of_a_grp_stream() {
        // The observed ACGTACGA is corrected to ACGTACGT, which the exact cell join matches.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nACGTACGAGGGG\n+\nIIIIIIIIIIII\n");
        let unassigned = dir.path().join("un.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "bc={x=ACGTACGT}".to_string(),
            "bc::loc=0:0:8,dist=1".to_string(),
            "cell={c1=ACGTACGT}".to_string(),
            "cell::match=b".to_string(),
        ];
        args.extracts = vec!["b=@bc".to_string()];
        args.sample_from_group = Some("cell".to_string());
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        let c1 = read_back(&dir.path().join("c1.fq"));
        assert_eq!(c1.len(), 1, "the read routed to cell c1");
        assert_eq!(c1[0].name, b"r");
        assert!(read_back(&unassigned).is_empty(), "nothing unassigned");
    }

    #[test]
    fn test_corrected_by_default_in_tag_and_body() {
        // A matched keeplist `@grp` stream is error-corrected by default
        // (raw=false) everywhere it is emitted: both the SAM tag and the
        // templated read body carry the canonical sequence.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nGAAGGGTTTT\n+\nIIIIIIIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec!["g={GAAGAG}".to_string(), "g::loc=0:0,dist=1".to_string()];
        args.extracts = vec!["bc=@g".to_string()];
        args.templates = vec!["bc".to_string()];
        args.tags = vec!["CB=bc".to_string(), "CB::qual=none".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let record = first_bam_record(&out);
        // Observed GAAGGG (1 mismatch) corrects to GAAGAG in both CB and the
        // body.
        assert_eq!(string_tag(&record, *b"CB").as_deref(), Some("GAAGAG"));
        let body: &[u8] = record.sequence().as_ref();
        assert_eq!(body, b"GAAGAG", "the body is corrected by default");
    }

    #[test]
    fn test_raw_true_tag_emits_observed() {
        // `raw=true` opts a tag back to the observed bases: CB carries the
        // canonical sequence while CR (raw=true) carries the observed barcode
        // from the same match.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nGAAGGGTTTT\n+\nIIIIIIIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec!["g={GAAGAG}".to_string(), "g::loc=0:0,dist=1".to_string()];
        args.extracts = vec!["bc=@g".to_string()];
        args.templates = vec!["bc".to_string()];
        args.tags = vec![
            "CB=bc".to_string(),
            "CB::qual=none".to_string(),
            "CR=bc".to_string(),
            "CR::raw=true".to_string(),
            "CR::qual=none".to_string(),
        ];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let record = first_bam_record(&out);
        assert_eq!(
            string_tag(&record, *b"CB").as_deref(),
            Some("GAAGAG"),
            "CB is corrected"
        );
        assert_eq!(
            string_tag(&record, *b"CR").as_deref(),
            Some("GAAGGG"),
            "CR (raw=true) is the observed barcode"
        );
    }

    #[test]
    fn test_template_raw_reverts_to_observed_body() {
        // `--template bc::raw=true` reverts the body to the observed bases
        // while a default tag still corrects: the body keeps GGGTCC and CB
        // carries the canonical GGATCC.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nGGGTCCTTTT\n+\nIIIIIIIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec!["g={GGATCC}".to_string(), "g::loc=0:0,dist=1".to_string()];
        args.extracts = vec!["bc=@g".to_string()];
        args.templates = vec!["bc::raw=true".to_string()];
        args.tags = vec!["CB=bc".to_string(), "CB::qual=none".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let record = first_bam_record(&out);
        let body: &[u8] = record.sequence().as_ref();
        assert_eq!(
            body, b"GGGTCC",
            "raw=true template reverts the body to observed"
        );
        assert_eq!(
            string_tag(&record, *b"CB").as_deref(),
            Some("GGATCC"),
            "the tag still corrects by default"
        );
    }

    #[test]
    fn test_three_group_sequential_chain() {
        // A bare-next chain grp_a -> grp_b -> grp_c: a read matching all three
        // passes; one missing the third group's tag fails it and is unassigned.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@all\nAAAACCCCGGGG\n+\nIIIIIIIIIIII\n@no_c\nAAAACCCCTTTT\n+\nIIIIIIIIIIII\n",
        );
        let out = dir.path().join("out.fq");
        let unassigned = dir.path().join("un.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp_a={AAAA}".to_string(),
            "grp_a::loc=0:0:4,next=grp_b,minFindsPerGroup=1".to_string(),
            "grp_b={CCCC}".to_string(),
            "grp_b::loc=0:4:8,next=grp_c,minFindsPerGroup=1".to_string(),
            "grp_c={GGGG}".to_string(),
            "grp_c::loc=0:8:12,minFindsPerGroup=1".to_string(),
        ];
        args.out = Some(out.display().to_string());
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        assert_eq!(read_back(&out)[0].name, b"all");
        assert_eq!(read_back(&unassigned)[0].name, b"no_c");
    }

    #[test]
    fn test_group_with_next_source_and_prev_honors_both() {
        // grp_b is BOTH the target of grp_a::next=grp_b AND declares
        // prev=grp_c, so it has two prerequisites. If grp_c does not match,
        // grp_b must be skipped even though grp_a matched and grp_b's own tag
        // is present (neither prerequisite may be silently dropped).
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nAAAACCCCTTTT\n+\nIIIIIIIIIIII\n");
        let out = dir.path().join("out.fq");
        let unassigned = dir.path().join("un.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp_a={AAAA}".to_string(),
            "grp_a::loc=0:0:4,next=grp_b".to_string(),
            "grp_c={GGGG}".to_string(),
            "grp_c::loc=0:8:12".to_string(),
            "grp_b={CCCC}".to_string(),
            "grp_b::loc=0:4:8,prev=grp_c,minFindsPerGroup=1".to_string(),
        ];
        args.out = Some(out.display().to_string());
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        // grp_c (GGGG) is absent at [8,12] (the read has TTTT), so grp_b's
        // prev=grp_c is unmet and grp_b is skipped -> the read is unassigned,
        // not passed through.
        let reached_out = out.exists() && !read_back(&out).is_empty();
        assert!(
            !reached_out,
            "read must not reach --out: grp_b's prev=grp_c is unmet"
        );
        let unmatched = read_back(&unassigned);
        assert_eq!(unmatched.len(), 1);
        assert_eq!(unmatched[0].name, b"r");
    }

    #[test]
    fn test_cyclic_next_links_fail_fast() {
        // grp_a -> grp_b -> grp_a is a dependency cycle; compilation fails
        // fast.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nAAAACCCC\n+\nIIIIIIII\n");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp_a={AAAA}".to_string(),
            "grp_a::next=grp_b".to_string(),
            "grp_b={CCCC}".to_string(),
            "grp_b::next=grp_a".to_string(),
        ];
        args.out = Some(dir.path().join("out.fq").display().to_string());
        let err = run_demux(args).unwrap_err();
        assert!(err.to_string().contains("cyclic"), "{err}");
    }

    #[test]
    fn test_prev_induced_cycle_fails_fast() {
        // A contradiction routed through `prev`: grp_b is grp_a's next AND
        // requires grp_c, while grp_c requires grp_b, so grp_b and grp_c each
        // depend on the other. Because both prerequisites of grp_b are kept
        // (not just the next= source), the cycle is detected.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nAAAACCCCGGGG\n+\nIIIIIIIIIIII\n");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp_a={AAAA}".to_string(),
            "grp_a::next=grp_b".to_string(),
            "grp_b={CCCC}".to_string(),
            "grp_b::prev=grp_c".to_string(),
            "grp_c={GGGG}".to_string(),
            "grp_c::prev=grp_b".to_string(),
        ];
        args.out = Some(dir.path().join("out.fq").display().to_string());
        let err = run_demux(args).unwrap_err();
        assert!(err.to_string().contains("cyclic"), "{err}");
    }

    #[test]
    fn test_use_site_tilde_flips_corrected_and_observed() {
        // `~` at the point of use reverse-complements that stream's
        // contribution there: `--tag CB=~bc` flips the corrected canonical, and
        // `--template ~bc::raw=true` flips the observed bases. The extract
        // itself keeps the orientation of the match.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nAGCCTTTT\n+\nIIIIIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec!["g={AACC}".to_string(), "g::loc=0:0:4,dist=1".to_string()];
        args.extracts = vec!["bc=@g".to_string()];
        args.templates = vec!["~bc::raw=true".to_string()];
        args.tags = vec!["CB=~bc".to_string(), "CB::qual=none".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let record = first_bam_record(&out);
        // CB: observed AGCC corrects to canonical AACC, then `~` -> GGTT (rc of
        // the canonical).
        assert_eq!(string_tag(&record, *b"CB").as_deref(), Some("GGTT"));
        // Body (`~` + raw=true): rc of the observed AGCC -> GGCT.
        let body: &[u8] = record.sequence().as_ref();
        assert_eq!(
            body, b"GGCT",
            "the ~ raw=true body is the rc of the observed bases"
        );
    }

    #[test]
    fn test_both_strands_group_emits_declared_canonical() {
        // A `bothStrands=true` group matches antisense, but the corrected tag
        // is always the declared canonical, and the raw CR is the observed
        // bases in the same (the tag's) orientation.
        let dir = tempfile::tempdir().unwrap();
        // [10:16] = CACGAT; tag AACGTG matches its rc (CACGTT) at dist=1.
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@r\nGGGGGGGGGGCACGATAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n",
        );
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "g={AACGTG}".to_string(),
            "g::loc=0:10:16,bothStrands=true,dist=1".to_string(),
        ];
        args.extracts = vec!["bc=@g".to_string()];
        args.tags = vec![
            "CB=bc".to_string(),
            "CB::qual=none".to_string(),
            "CR=bc".to_string(),
            "CR::raw=true".to_string(),
            "CR::qual=none".to_string(),
        ];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let record = first_bam_record(&out);
        // CB: the declared canonical AACGTG, regardless of the antisense match
        // (use `~` for read strand).
        assert_eq!(string_tag(&record, *b"CB").as_deref(), Some("AACGTG"));
        // CR (raw): the observed read bases CACGAT reverse-complemented into
        // the tag's orientation, one substitution from AACGTG.
        assert_eq!(string_tag(&record, *b"CR").as_deref(), Some("ATCGTG"));
    }

    #[test]
    fn test_both_strands_corrected_quals_co_orient_with_declared() {
        // The corrected qualities must co-orient with the DECLARED canonical
        // bases. For an antisense match the read's quals run opposite to the
        // declared orientation, so they are reversed into it.
        // Distinct per-base quals make the orientation observable.
        let dir = tempfile::tempdir().unwrap();
        // [10:16] = CACGAT with quals raw 20..=25 (Phred+33 = "56789:"); tag
        // AACGTG matches its rc.
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@r\nGGGGGGGGGGCACGATAAAA\n+\nIIIIIIIIII56789:IIII\n",
        );
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "g={AACGTG}".to_string(),
            "g::loc=0:10:16,bothStrands=true,dist=1".to_string(),
        ];
        args.extracts = vec!["bc=@g".to_string()];
        // CB defaults to the CY quality tag, carrying the corrected qualities.
        args.tags = vec!["CB=bc".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let record = first_bam_record(&out);
        assert_eq!(string_tag(&record, *b"CB").as_deref(), Some("AACGTG"));
        // CY reversed to ":98765" so it lines up with the declared AACGTG (a
        // co-orientation bug would leave the read-strand "56789:").
        assert_eq!(string_tag(&record, *b"CY").as_deref(), Some(":98765"));
    }

    #[test]
    fn test_both_strands_group_forward_match_quals_not_reversed() {
        // A `bothStrands=true` group also searches the forward strand. On a
        // forward match the observed quals already co-orient with the declared
        // canonical, so they must NOT be reversed: the per-match strand governs
        // the correction, not the group's `bothStrands` attribute.
        let dir = tempfile::tempdir().unwrap();
        // [10:16] = AACGTG (the declared tag, forward) with quals raw 20..=25
        // (Phred+33 = "56789:").
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@r\nGGGGGGGGGGAACGTGAAAA\n+\nIIIIIIIIII56789:IIII\n",
        );
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "g={AACGTG}".to_string(),
            "g::loc=0:10:16,bothStrands=true,dist=1".to_string(),
        ];
        args.extracts = vec!["bc=@g".to_string()];
        // CB defaults to the CY quality tag, carrying the corrected qualities.
        args.tags = vec!["CB=bc".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let record = first_bam_record(&out);
        assert_eq!(string_tag(&record, *b"CB").as_deref(), Some("AACGTG"));
        // CY stays "56789:" (read order); reversing it to ":98765" would be the
        // bug of using the group's `bothStrands` config instead of the match
        // strand.
        assert_eq!(string_tag(&record, *b"CY").as_deref(), Some("56789:"));
    }

    #[test]
    fn test_a_read_and_its_reverse_complement_carry_identical_streams() {
        struct Lcg(u64);
        impl Lcg {
            fn below(&mut self, n: usize) -> usize {
                self.0 = self
                    .0
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                (self.0 >> 33) as usize % n
            }
            fn bases(&mut self, n: usize) -> Vec<u8> {
                (0..n).map(|_| b"ACGT"[self.below(4)]).collect()
            }
        }
        let barcodes: [&[u8]; 4] = [b"ACGTTGCA", b"TTGCAACG", b"GATCCTAG", b"CCATGGTA"];
        let mut rng = Lcg(0x5eed);
        let mut fastq = Vec::new();
        for i in 0..300 {
            let mut barcode = barcodes[rng.below(4)].to_vec();
            if rng.below(3) == 0 {
                let at = rng.below(barcode.len());
                barcode[at] = b"ACGT"[rng.below(4)];
            }
            let lengths = [rng.below(12), rng.below(4), rng.below(16), rng.below(12)];
            let read = [
                rng.bases(lengths[0]),
                b"CTACACGACGCTCTTCCGATCT".to_vec(),
                rng.bases(lengths[1]),
                barcode,
                rng.bases(10 + lengths[2]),
                b"GGTACCTTAGCA".to_vec(),
                rng.bases(lengths[3]),
            ]
            .concat();
            let quals: Vec<u8> = (0..read.len())
                .map(|_| b'#' + rng.below(40) as u8)
                .collect();
            let mut flipped = read.clone();
            crate::iupac::reverse_complement(&mut flipped);
            let flipped_quals: Vec<u8> = quals.iter().rev().copied().collect();
            for (name, bases, quals) in [("f", &read, &quals), ("r", &flipped, &flipped_quals)] {
                fastq.extend(format!("@{name}{i}\n").bytes());
                fastq.extend(
                    bases
                        .iter()
                        .chain(b"\n+\n")
                        .chain(quals.iter())
                        .chain(b"\n"),
                );
            }
        }
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", &fastq);
        let out = dir.path().join("out.sam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "ad={CTACACGACGCTCTTCCGATCT}".to_string(),
            "ad::bothStrands=true,next=bc:0-12".to_string(),
            "bc={a=ACGTTGCA,b=TTGCAACG,c=GATCCTAG,d=CCATGGTA}".to_string(),
            "bc::dist=1".to_string(),
            "tr={GGTACCTTAGCA}".to_string(),
            "tr::bothStrands=true".to_string(),
        ];
        args.extracts = vec![
            "cb=@bc".to_string(),
            "umi=@bc+0:10".to_string(),
            "up=@ad-0:6".to_string(),
            "gap=@ad..@bc".to_string(),
            "ins=@bc..@tr".to_string(),
            "tail=@tr+0:end".to_string(),
        ];
        args.tags = [
            "CB=cb",
            "CR=cb",
            "CR::raw=true",
            "CR::qual=XC",
            "RX=umi",
            "XU=up",
            "XU::qual=XV",
            "XG=gap",
            "XG::qual=XH",
            "XI=ins",
            "XI::qual=XJ",
            "XT=tail",
            "XT::qual=XW",
        ]
        .map(String::from)
        .to_vec();
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let text = std::fs::read_to_string(&out).unwrap();
        let records: Vec<Vec<&str>> = text
            .lines()
            .filter(|line| !line.starts_with('@'))
            .map(|line| line.split('\t').skip(11).collect())
            .collect();
        assert_eq!(records.len(), 600);
        for pair in records.chunks(2) {
            assert_eq!(pair[0], pair[1]);
        }
        let barcoded = records
            .iter()
            .filter(|tags| tags.iter().any(|tag| tag.starts_with("CB:Z:")))
            .count();
        assert!(barcoded > 500, "{barcoded} of 600 records carry a CB");
    }

    #[test]
    fn test_partial_corrects_to_full_canonical() {
        // A partial (truncated) match emits the full canonical tag in CB
        // (padding the corrected qualities to its length); a raw=true template
        // keeps only the observed bases that matched.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nCGATCGTTTT\n+\nIIIIIIIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "g={ATCGATCG}".to_string(),
            "g::loc=0:0,partial5=4:0.1".to_string(),
        ];
        args.extracts = vec!["bc=@g".to_string()];
        args.templates = vec!["bc::raw=true".to_string()];
        // CB defaults to the CY quality tag, which carries the padded corrected
        // qualities.
        args.tags = vec!["CB=bc".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let record = first_bam_record(&out);
        // CB: only CGATCG (6 nt) matched, but correction emits the full
        // canonical ATCGATCG (8 nt).
        assert_eq!(string_tag(&record, *b"CB").as_deref(), Some("ATCGATCG"));
        assert_eq!(
            string_tag(&record, *b"CY").map(|q| q.len()),
            Some(8),
            "corrected qualities padded to the canonical length"
        );
        // Body (raw=true): the observed 6 nt that matched.
        let body: &[u8] = record.sequence().as_ref();
        assert_eq!(
            body, b"CGATCG",
            "the raw=true body keeps the observed matched bases"
        );
    }

    #[test]
    fn test_movable_template_honors_corrected_and_raw() {
        // A corrected `@grp` stream consumed by exactly one `--template` and NO
        // `--tag` is "movable" (its bases are moved into the body instead of
        // copied). That fast-path must still honor the corrected form by
        // default and revert to observed under raw=true.
        let run = |template: &str, name: &str| -> Vec<u8> {
            let dir = tempfile::tempdir().unwrap();
            let input = write_file(dir.path(), "in.fq", b"@r\nGAAGGGTTTT\n+\nIIIIIIIIII\n");
            let out = dir.path().join(name);
            let mut args = base_args();
            args.inputs = vec![format!("0={}", input.display())];
            args.groups = vec!["g={GAAGAG}".to_string(), "g::loc=0:0,dist=1".to_string()];
            args.extracts = vec!["bc=@g".to_string()];
            // Exactly one consumer (this template), no --tag, so `bc` is
            // movable.
            args.templates = vec![template.to_string()];
            args.out = Some(out.display().to_string());
            run_demux(args).unwrap();
            read_back(&out)[0].bases.clone()
        };
        // Default (raw=false): the movable body is the corrected canonical
        // GAAGAG.
        assert_eq!(run("bc", "corrected.fq"), b"GAAGAG");
        // raw=true: the movable body reverts to the observed GAAGGG.
        assert_eq!(run("bc::raw=true", "raw.fq"), b"GAAGGG");
    }

    #[test]
    fn test_anchor_offset_not_corrected() {
        // Correction rewrites only the @grp self-span; an @grp+offset region
        // has no corrected form, so its body keeps the observed bases.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nAGCCGGGGTT\n+\nIIIIIIIIII\n");
        let out = dir.path().join("out.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec!["g={AACC}".to_string(), "g::loc=0:0:4,dist=1".to_string()];
        // @g matches AGCC at [0,4]; @g+0:4 is [4,8] = GGGG, which correction
        // must leave untouched.
        args.extracts = vec!["off=@g+0:4".to_string()];
        args.templates = vec!["off".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        assert_eq!(read_back(&out)[0].bases, b"GGGG");
    }

    #[test]
    fn test_anchor5p_window_overrun_warning() {
        use crate::grammar::{Anchor, Endpoint, Location};
        use crate::matcher::CompiledGroup;
        let mut g = CompiledGroup::new("grp", vec![vec![b'A'; 9]]);
        g.anchor = Some(Anchor::FivePrime);
        g.loc = Some(Location {
            file: 0,
            start: Endpoint::FromStart(0),
            end: Some(Endpoint::FromStart(8)),
        });
        // A 9 nt tag in an 8-wide anchored window overruns the declared end ->
        // warn.
        assert!(anchor5p_window_overrun_warning(&g)
            .is_some_and(|w| w.contains("runs past the loc end")));
        // Tags within the window, or no concrete end, or an unanchored window
        // -> no warning.
        g.tags = vec![vec![b'A'; 8], vec![b'A'; 6]];
        assert!(anchor5p_window_overrun_warning(&g).is_none());
        g.tags = vec![vec![b'A'; 9]];
        g.loc = Some(Location {
            file: 0,
            start: Endpoint::FromStart(0),
            end: None,
        });
        assert!(
            anchor5p_window_overrun_warning(&g).is_none(),
            "no concrete end"
        );
        g.anchor = None;
        assert!(
            anchor5p_window_overrun_warning(&g).is_none(),
            "window is unaffected"
        );
    }

    #[test]
    fn compile_groups_marks_read_group() {
        // An @RG group's compiled tags are synthetic routing tokens (e.g.
        // `rg1#0`), not DNA, so compile_groups must skip the IUPAC
        // validation/sassy encode() sequence machinery for it and instead
        // carry the RG-id -> tag-idx map on `read_group`. A normal group's
        // `read_group` stays `None`.
        let mut args = base_args();
        args.inputs = vec!["0=/dev/null".to_string()];
        args.groups = vec!["norm={s1=AAAA}".to_string(), "rg=@RG".to_string()];
        let plan = crate::grammar::parse_demux(&args).unwrap();

        let mut rg_to_idx = HashMap::new();
        rg_to_idx.insert(b"rg1".to_vec(), 0usize);
        let resolved: HashMap<String, crate::ResolvedReadGroup> = HashMap::from([(
            "rg".to_string(),
            crate::ResolvedReadGroup {
                tag_set: crate::tags::TagSet {
                    entries: vec![crate::tags::TagEntry {
                        id: "rg1#0".to_string(),
                        seq: "rg1#0".to_string(),
                        sub_sample: None,
                    }],
                },
                samples: Vec::new(),
                rg_to_idx,
            },
        )]);

        let (tag_sets, _) = load_tag_sets(&plan, &resolved).unwrap();
        let (compiled, _index) = compile_groups(&plan, &tag_sets, &resolved).unwrap();

        let norm = compiled.iter().find(|g| g.name == "norm").unwrap();
        assert!(
            norm.read_group.is_none(),
            "a normal group's read_group is None"
        );

        let rg = compiled.iter().find(|g| g.name == "rg").unwrap();
        let read_group = rg
            .read_group
            .as_ref()
            .expect("an @RG group's read_group is Some");
        assert_eq!(read_group[b"rg1".as_slice()], 0);
        // The sequence machinery is skipped: no batched-search encoding, and
        // the rest of the record-window attributes stay at their bare
        // defaults (the grammar rejects any attrs on an `@RG` group).
        assert!(
            rg.encoded.is_empty(),
            "an @RG group skips sassy encode() (its tags are routing tokens, not DNA)"
        );
        assert_eq!(rg.dist, crate::grammar::Dist::default());
        assert!(rg.loc.is_none());
    }

    #[test]
    fn test_widest_window_is_the_largest_over_every_record_length() {
        for (attrs, widest) in [
            ("g::loc=0:2:10", Some(8)),
            ("g::loc=0:10:2", Some(0)),
            ("g::loc=0:-10", Some(10)),
            ("g::loc=0:-10:-4", Some(6)),
            ("g::loc=0:-10:6", Some(6)),
            ("g::loc=0:2", None),
            ("g::loc=0:2:-3", None),
            ("g::loc=0:2:10,anchor=5p", None),
            ("g::loc=0:2:10,anchor=3p", Some(8)),
            ("g::dist=1", None),
        ] {
            let mut args = base_args();
            args.inputs = vec!["0=/dev/null".to_string()];
            args.groups = vec!["g={a=ACGT}".to_string(), attrs.to_string()];
            let plan = crate::grammar::parse_demux(&args).unwrap();
            assert_eq!(
                widest_window(&plan.groups[0], &plan.groups),
                widest,
                "{attrs}"
            );
        }
        let mut args = base_args();
        args.inputs = vec!["0=/dev/null".to_string()];
        args.groups = ["u={a=AAAA}", "u::loc=0:0:4,next=g:2-9", "g={a=ACGT}"]
            .map(str::to_string)
            .to_vec();
        let plan = crate::grammar::parse_demux(&args).unwrap();
        let target = plan.groups.iter().find(|group| group.name == "g").unwrap();
        assert_eq!(widest_window(target, &plan.groups), Some(7));
    }

    #[test]
    fn test_compile_groups_defaults_nearest_delta_to_one() {
        let mut args = base_args();
        args.inputs = vec!["0=/dev/null".to_string()];
        args.groups = vec![
            "g={s1=TGAGC,s2=TGAGA}".to_string(),
            "g::loc=0:0:5,dist=1,mode=nearest".to_string(),
        ];
        let plan = crate::grammar::parse_demux(&args).unwrap();
        let resolved = HashMap::new();
        let (tag_sets, _) = load_tag_sets(&plan, &resolved).unwrap();
        let (compiled, _index) = compile_groups(&plan, &tag_sets, &resolved).unwrap();
        let outcome =
            |read: &[u8]| crate::matcher::match_group(&compiled[0], &mut Scratch::new(), &[read]);
        let matched_tags = |result: &crate::matcher::GroupOutcome| -> BTreeSet<usize> {
            result.accepted.iter().map(|hit| hit.tag_idx).collect()
        };

        let winner = outcome(b"TGAGCCCC");
        assert_eq!(matched_tags(&winner), BTreeSet::from([0, 1]));
        assert_eq!(
            winner.best.map(|hit| (hit.tag_idx, hit.cost)),
            Some((0, 0)),
            "a one-cost gap over the runner-up meets the default delta"
        );

        let tie = outcome(b"TGAGTCCC");
        assert_eq!(matched_tags(&tie), BTreeSet::from([0, 1]));
        assert_eq!(tie.best, None, "a tie is left unmatched with no delta");
    }

    #[test]
    fn test_anchor5p_relative_next_window() {
        // anchor5p as a relative next= window target (no own loc): grp_a
        // matches AAAA at [0,4], its next=grp_b:0-4 window starts at grp_a's
        // match end (4), and grp_b (anchor5p) anchors CCCC there -> BC carries
        // the [4,8) barcode. Validates anchor5p + a relative next= window.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nAAAACCCCGGGG\n+\nIIIIIIIIIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp_a={AAAA}".to_string(),
            "grp_a::loc=0:0:4,next=grp_b:0-4,minFindsPerGroup=1".to_string(),
            "grp_b={b=CCCC}".to_string(),
            "grp_b::anchor=5p,minFindsPerGroup=1".to_string(),
        ];
        args.extracts = vec!["bc=@grp_b".to_string()];
        args.tags = vec!["BC=bc".to_string(), "BC::qual=none".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let record = first_bam_record(&out);
        assert_eq!(
            string_tag(&record, *b"BC").as_deref(),
            Some("CCCC"),
            "grp_b anchored at the relative window start (grp_a end = 4)"
        );
    }

    #[test]
    fn test_anchored_resolves_variable_length_group() {
        // End-to-end split-seq-like grp_cbt fix: variable-length 5'-anchored
        // tags routed by sample. Read AGTACTCT matches the 7nt tag dna01 at
        // [0,7); under an unanchored window the 6nt tag dna02 slides to a
        // spurious offset (2 hits -> maxFindsPerGroup=1 drops the read to
        // unassigned), but anchor=5p pins each tag at the start so only dna01
        // matches and the read routes to dna01.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nAGTACTCTGGGG\n+\nIIIIIIIIIIII\n");
        let un = dir.path().join("un.fq");
        let run = |place_attr: &str| -> (usize, usize) {
            let mut args = base_args();
            args.inputs = vec![format!("0={}", input.display())];
            args.groups = vec![
                "grp={dna01=AGTACTC,dna02=TACTCA}".to_string(),
                format!("grp::loc=0:0:8,dist=1,maxFindsPerGroup=1{place_attr}"),
            ];
            args.sample_from_group = Some("grp".to_string());
            args.out = Some(format!("{}/%sample.fq", dir.path().display()));
            args.unassigned = Some(un.display().to_string());
            run_demux(args).unwrap();
            (
                read_back(&dir.path().join("dna01.fq")).len(),
                read_back(&un).len(),
            )
        };
        // an unanchored window: the spurious 2nd hit trips maxFindsPerGroup=1
        // -> unassigned.
        assert_eq!(run(""), (0, 1), "window slides -> dropped");
        // anchor=5p: a single hit -> routed to dna01.
        assert_eq!(run(",anchor=5p"), (1, 0), "anchored -> dna01");
    }

    #[test]
    fn test_sam_input_tags_carry_through_to_sam() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        // Pass-through SAM -> SAM preserves arbitrary, non-alignment, non-demux
        // input tags. (bwa's XA/XB would be stripped as alignment-derived, so
        // this uses custom local tags ZX/ZY.)
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.sam",
            b"@HD\tVN:1.6\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tZX:Z:keepme\tZY:i:42\n",
        );
        let out = dir.path().join("out.sam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let records = read_back(&out);
        assert_eq!(records.len(), 1);
        let data = records[0].tags.as_ref().expect("tags carried through");
        match data.get(b"ZX") {
            Some(Value::String(value)) => assert_eq!(value.to_string(), "keepme"),
            other => panic!("ZX missing or wrong type: {other:?}"),
        }
        assert!(data.get(b"ZY").is_some(), "ZY carried through");
    }

    #[test]
    fn test_sam_input_tags_carry_through_to_fastq_comment() {
        // SAM -> FASTQ carries the tags into the read-name comment (samtools
        // `fastq -T` style). Custom local tags (ZX/ZY); bwa XA/XB would be
        // stripped as alignment-derived.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.sam",
            b"@HD\tVN:1.6\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tZX:Z:keepme\tZY:i:42\n",
        );
        let out = dir.path().join("out.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let text = std::fs::read_to_string(&out).unwrap();
        assert!(
            text.contains("ZX:Z:keepme"),
            "FASTQ comment carries ZX: {text}"
        );
        assert!(text.contains("ZY:i:42"), "FASTQ comment carries ZY: {text}");
    }

    #[test]
    fn test_paired_templates_to_bam_set_pair_flags() {
        // Two templates to a BAM emit a properly-flagged read pair (first/last
        // segment).
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r1\nAAAAACCCCC\n+\nIIIIIIIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.extracts = vec!["a=0:0:5".to_string(), "b=0:5:end".to_string()];
        args.templates = vec!["a".to_string(), "b".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let mut reader = noodles::bam::io::Reader::new(File::open(&out).unwrap());
        let header = reader.read_header().unwrap();
        let records: Vec<_> = reader
            .record_bufs(&header)
            .collect::<std::io::Result<_>>()
            .unwrap();
        assert_eq!(records.len(), 2);
        assert!(records[0].flags().is_segmented() && records[0].flags().is_first_segment());
        assert!(records[1].flags().is_segmented() && records[1].flags().is_last_segment());
        let first: &[u8] = records[0].sequence().as_ref();
        let second: &[u8] = records[1].sequence().as_ref();
        assert_eq!(first.to_vec(), b"AAAAA".to_vec());
        assert_eq!(second.to_vec(), b"CCCCC".to_vec());
    }

    #[test]
    fn test_metrics_tsvs_report_fan_out_and_unextracted() {
        // Per-sample + summary metrics TSVs are written with the right counts
        // and a non-zero frac_bases_unextracted when a --template drops part of
        // the read.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            // s1: AAAA-prefixed (8 bp); s2: TTTT-prefixed; one unmatched.
            b"@a\nAAAACCCC\n+\nIIIIIIII\n@b\nTTTTGGGG\n+\nIIIIIIII\n@u\nGGGGGGGG\n+\nIIIIIIII\n",
        );
        let per_sample = dir.path().join("per_sample.tsv");
        let summary = dir.path().join("summary.tsv");
        let mut args = base_args();
        args.pool = Some("pool1".to_string());
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={s1=AAAA,s2=TTTT}".to_string(),
            "grp::loc=0:0:4,dist=0".to_string(),
        ];
        // Carve only the first 4 bp into the body, dropping the trailing 4 (so
        // 4/8 unextracted).
        args.extracts = vec!["bc=0:0:4".to_string()];
        args.templates = vec!["bc".to_string()];
        args.samples = vec!["s1=grp::s1".to_string(), "s2=grp::s2".to_string()];
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));
        args.metrics_per_sample = Some(per_sample.clone());
        args.metrics_summary = Some(summary.clone());
        run_demux(args).unwrap();

        let per = std::fs::read_to_string(&per_sample).unwrap();
        // s1 and s2 each get one read; half the bases (4 of 8) are unextracted.
        assert!(per.contains("pool1\ts1\t\t1\t"), "s1 row: {per}");
        assert!(
            per.contains("0.500000"),
            "frac_bases_unextracted 4/8: {per}"
        );
        let sum = std::fs::read_to_string(&summary).unwrap();
        // total 3, assigned 2, unassigned 1.
        assert!(
            sum.lines()
                .nth(1)
                .unwrap()
                .starts_with("pool1\t3\t2\t0\t1\t0\t"),
            "summary: {sum}"
        );
    }

    #[test]
    fn test_without_format_extension() {
        for (name, stem) in [
            ("lib1.fq", "lib1"),
            ("lib1.fastq.gz", "lib1"),
            ("LIB.BAM", "LIB"),
            ("a.b.cram", "a.b"),
            ("reads.txt", "reads.txt"),
            ("reads.txt.gz", "reads.txt.gz"),
            ("x.gz", "x.gz"),
            ("fq", "fq"),
            (".fq", ".fq"),
            (".fq.gz", ".fq.gz"),
            ("\u{e9}.fq", "\u{e9}"),
            ("a\u{e9}", "a\u{e9}"),
        ] {
            assert_eq!(without_format_extension(name), stem, "{name}");
        }
    }

    #[test]
    fn test_longest_common_prefix_ends_at_a_whole_character() {
        let prefix = |names: &[&str]| {
            longest_common_prefix(
                &names
                    .iter()
                    .map(|name| name.to_string())
                    .collect::<Vec<_>>(),
            )
        };
        assert_eq!(prefix(&["\u{e9}1", "\u{e8}1"]), "");
        assert_eq!(prefix(&["a\u{e9}1", "a\u{e8}1", "a\u{ea}"]), "a");
        assert_eq!(prefix(&["lib_R1", "lib_R2"]), "lib_R");
        assert_eq!(prefix(&["solo"]), "solo");
        assert_eq!(prefix(&[]), "");
    }

    #[test]
    fn test_metrics_paths_expand_pool_placeholder() {
        // `%pool` in a `--metrics-*` path expands to the pool id (the same
        // value `--out` uses), so the files land at the resolved names rather
        // than at literal `%pool` paths.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@a\nAAAACCCC\n+\nIIIIIIII\n");
        let mut args = base_args();
        args.pool = Some("lib01".to_string());
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={s1=AAAA}".to_string(),
            "grp::loc=0:0:4,dist=0".to_string(),
        ];
        args.extracts = vec!["bc=0:0:4".to_string()];
        args.templates = vec!["bc".to_string()];
        args.samples = vec!["s1=grp::s1".to_string()];
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));
        args.metrics_per_sample = Some(dir.path().join("%pool.per_sample.tsv"));
        args.metrics_summary = Some(dir.path().join("%pool.summary.tsv"));
        run_demux(args).unwrap();

        // The literal `%pool` path is never created; the expanded `lib01.*`
        // paths are.
        assert!(!dir.path().join("%pool.per_sample.tsv").exists());
        assert!(!dir.path().join("%pool.summary.tsv").exists());
        let per = std::fs::read_to_string(dir.path().join("lib01.per_sample.tsv")).unwrap();
        assert!(per.contains("lib01\ts1\t"), "per-sample row: {per}");
        let sum = std::fs::read_to_string(dir.path().join("lib01.summary.tsv")).unwrap();
        assert!(sum.lines().nth(1).unwrap().starts_with("lib01\t"), "{sum}");
    }

    #[test]
    fn test_metrics_path_rejects_non_pool_placeholder() {
        // A metrics file is pool-level, so a `%sample` placeholder cannot
        // resolve: it is a fail-fast parse error (clear message), not a
        // silently-literal path.
        let mut args = base_args();
        args.inputs = vec!["0=/dev/null".to_string()];
        args.metrics_per_sample = Some(PathBuf::from("out/%sample.tsv"));
        let err = crate::grammar::parse_demux(&args).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("--metrics-per-sample"),
            "names the flag: {msg}"
        );
        assert!(
            msg.contains("%pool"),
            "points at the valid placeholder: {msg}"
        );
        assert!(msg.contains("not allowed"), "rejects %sample: {msg}");
    }

    #[test]
    fn test_require_samples_explain_all_tags_flags_unclaimed_tag() {
        // A group with two barcodes but only one claimed by a sample:
        // `--require-samples-explain-all-tags` fails fast and names the orphan;
        // off (the default) the same spec runs fine.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@a\nAAAACCCC\n+\nIIIIIIII\n");
        let groups = vec![
            "grp={s1=AAAA,s2=TTTT}".to_string(),
            "grp::loc=0:0:4,dist=0".to_string(),
        ];

        // Off by default: claiming only s1 is allowed.
        let mut ok = base_args();
        ok.inputs = vec![format!("0={}", input.display())];
        ok.groups = groups.clone();
        ok.extracts = vec!["bc=0:0:4".to_string()];
        ok.templates = vec!["bc".to_string()];
        ok.samples = vec!["a=grp::s1".to_string()];
        ok.out = Some(format!("{}/%sample.fq", dir.path().display()));
        run_demux(ok).unwrap();

        // On: the unclaimed s2 is a fail-fast error that names the flag and the
        // orphan tag.
        let mut strict = base_args();
        strict.inputs = vec![format!("0={}", input.display())];
        strict.groups = groups;
        strict.extracts = vec!["bc=0:0:4".to_string()];
        strict.templates = vec!["bc".to_string()];
        strict.samples = vec!["a=grp::s1".to_string()];
        strict.require_samples_explain_all_tags = true;
        strict.out = Some(format!("{}/%sample.fq", dir.path().display()));
        let err = run_demux(strict).unwrap_err().to_string();
        assert!(
            err.contains("require-samples-explain-all-tags"),
            "names the flag: {err}"
        );
        assert!(err.contains("grp::s2"), "names the unclaimed tag: {err}");
    }

    /// The first bytes of a BGZF block: gzip magic with the FEXTRA flag set,
    /// then the `BC` subfield. flate2's plain gzip has no FEXTRA, so this
    /// distinguishes a pool-compressed file from the inline gzip path.
    fn is_bgzf(path: &Path) -> bool {
        let bytes = std::fs::read(path).unwrap();
        bytes.len() >= 14 && bytes[0..4] == [0x1f, 0x8b, 0x08, 0x04] && &bytes[12..14] == b"BC"
    }

    #[test]
    fn test_fastx_gz_is_bgzf_and_round_trips() {
        // Each .fq.gz of a gzipped-FASTX fan-out is BGZF-framed (valid gzip), round-trips its reads, and a
        // zero-read sample is an empty but valid BGZF file (the EOF block
        // only).
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@a\nAAAACCCC\n+\nIIIIIIII\n@b\nTTTTGGGG\n+\nIIIIIIII\n",
        );
        let mut args = base_args();
        args.threads = 8;
        args.inputs = vec![format!("0={}", input.display())];
        // s3 is declared but matches nothing in this input, so its .fq.gz is
        // empty.
        args.groups = vec![
            "grp={s1=AAAA,s2=TTTT,s3=CCCC}".to_string(),
            "grp::loc=0:0:4,dist=0".to_string(),
        ];
        args.samples = vec![
            "s1=grp::s1".to_string(),
            "s2=grp::s2".to_string(),
            "s3=grp::s3".to_string(),
        ];
        args.out = Some(format!("{}/%sample.fq.gz", dir.path().display()));
        run_demux(args).unwrap();

        let s1 = dir.path().join("s1.fq.gz");
        assert!(is_bgzf(&s1), "the pool BGZF-framed the .fq.gz");
        let s1_recs = read_back(&s1);
        assert_eq!(s1_recs.len(), 1);
        assert_eq!(s1_recs[0].bases, b"AAAACCCC");
        assert_eq!(
            read_back(&dir.path().join("s2.fq.gz"))[0].bases,
            b"TTTTGGGG"
        );

        // The zero-read sample's file exists and is a valid (BGZF EOF only)
        // gzip.
        let s3 = dir.path().join("s3.fq.gz");
        assert!(s3.exists(), "a directed zero-read .fq.gz still exists");
        assert!(is_bgzf(&s3), "and is a valid BGZF stream");
    }

    #[test]
    fn test_bam_output_honors_compression_level() {
        // unmux maps --compression to the BAM BGZF level. The same data at
        // level 0 (stored) must produce a clearly larger file than at level 9
        // (max), and both must round-trip every record.
        let dir = tempfile::tempdir().unwrap();
        let bases = *b"ACGT";
        let n = 1000u32;
        let mut fq = Vec::new();
        for i in 0..n {
            fq.extend_from_slice(format!("@r{i}\n").as_bytes());
            for j in 0..40u32 {
                let h = i
                    .wrapping_mul(2_654_435_761)
                    .wrapping_add(j.wrapping_mul(40_503));
                fq.push(bases[((h >> 13) & 3) as usize]);
            }
            fq.extend_from_slice(b"\n+\n");
            fq.extend_from_slice(&[b'I'; 40]);
            fq.push(b'\n');
        }
        let input = write_file(dir.path(), "in.fq", &fq);

        let run = |level: u8, name: &str| -> PathBuf {
            let out = dir.path().join(name);
            let mut args = base_args();
            args.compression = level;
            args.inputs = vec![format!("0={}", input.display())];
            args.out = Some(out.display().to_string());
            run_demux(args).unwrap();
            out
        };
        let stored = run(0, "c0.bam");
        let packed = run(9, "c9.bam");

        let size = |p: &Path| std::fs::metadata(p).unwrap().len();
        assert!(
            size(&stored) > size(&packed),
            "level 0 ({}) must be larger than level 9 ({})",
            size(&stored),
            size(&packed)
        );
        // Both round-trip every record.
        for path in [&stored, &packed] {
            let mut reader = noodles::bam::io::Reader::new(File::open(path).unwrap());
            let header = reader.read_header().unwrap();
            let count = reader.record_bufs(&header).count();
            assert_eq!(
                count,
                n as usize,
                "all records survive in {}",
                path.display()
            );
        }
    }

    #[test]
    fn test_bam_unassigned_bin_round_trips() {
        // A BAM `--unassigned` bin must read back, with one BGZF layer and a sound EOF.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@a\nGGGG\n+\nIIII\n"); // matches no sample
        let mut args = base_args();
        args.threads = 4;
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={s1=AAAA}".to_string(),
            "grp::loc=0:0:4,dist=0".to_string(),
        ];
        args.samples = vec!["s1=grp::s1".to_string()];
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));
        args.unassigned = Some(format!("{}/un.bam", dir.path().display()));
        run_demux(args).unwrap();

        let un = dir.path().join("un.bam");
        let mut reader = noodles::bam::io::Reader::new(File::open(&un).unwrap());
        let header = reader.read_header().unwrap();
        let records: Vec<_> = reader
            .record_bufs(&header)
            .collect::<std::result::Result<_, _>>()
            .unwrap();
        assert_eq!(
            records.len(),
            1,
            "the unassigned read is a valid BAM record"
        );
    }

    #[test]
    fn test_bgzf_error_path_returns_clean_error_not_panic() {
        // A write error while BGZF outputs are open must surface a clean Err,
        // not a panic. The run errors because a 3-input pass-through emits 3
        // reads per record, which alignment output cannot hold, beside a
        // `.fq.gz` unassigned bin.
        let dir = tempfile::tempdir().unwrap();
        let in0 = write_file(dir.path(), "0.fq", b"@r\nAAAA\n+\nIIII\n");
        let in1 = write_file(dir.path(), "1.fq", b"@r\nCCCC\n+\nIIII\n");
        let in2 = write_file(dir.path(), "2.fq", b"@r\nGGGG\n+\nIIII\n");
        let mut args = base_args();
        args.threads = 4;
        args.inputs = vec![
            format!("0={}", in0.display()),
            format!("1={}", in1.display()),
            format!("2={}", in2.display()),
        ];
        args.out = Some(format!("{}/p.bam", dir.path().display()));
        args.unassigned = Some(format!("{}/unmatched.fq.gz", dir.path().display()));

        let err = run_demux(args).unwrap_err();
        assert!(
            err.to_string().contains("at most two reads"),
            "the actionable write error surfaces, not a pool panic: {err}"
        );
    }

    #[test]
    fn test_fastx_gz_deterministic_across_threads() {
        // BGZF output is byte-identical regardless of --threads: merging cuts
        // blocks in input order and BGZF block compression is deterministic at
        // a fixed level.
        let dir = tempfile::tempdir().unwrap();
        let mut fq = Vec::new();
        for i in 0..500 {
            let bc = if i % 2 == 0 { "AAAA" } else { "TTTT" };
            fq.extend_from_slice(format!("@r{i}\n{bc}CGCG\n+\nIIIIIIII\n").as_bytes());
        }
        let input = write_file(dir.path(), "in.fq", &fq);
        let run = |threads: usize, tag: &str| {
            let mut args = base_args();
            args.threads = threads;
            args.inputs = vec![format!("0={}", input.display())];
            args.groups = vec![
                "grp={s1=AAAA,s2=TTTT}".to_string(),
                "grp::loc=0:0:4,dist=0".to_string(),
            ];
            args.samples = vec!["s1=grp::s1".to_string(), "s2=grp::s2".to_string()];
            args.out = Some(format!("{}/{tag}_%sample.fq.gz", dir.path().display()));
            run_demux(args).unwrap();
        };
        run(4, "a");
        run(16, "b");
        for sample in ["s1", "s2"] {
            assert_eq!(
                std::fs::read(dir.path().join(format!("a_{sample}.fq.gz"))).unwrap(),
                std::fs::read(dir.path().join(format!("b_{sample}.fq.gz"))).unwrap(),
                "pooled {sample}.fq.gz identical across thread counts"
            );
        }
    }

    #[test]
    fn test_sample_fanout_to_per_sample_fastq() {
        // Two samples partition a barcode group; --out %sample writes one file
        // per sample.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@a\nAAAAAAAA\n+\nIIIIIIII\n@b\nTTTTTTTT\n+\nIIIIIIII\n",
        );
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={D1=AAAAAAAA,D2=TTTTTTTT}".to_string(),
            "grp::loc=0:0,dist=0".to_string(),
        ];
        args.samples = vec!["s1=grp::D1".to_string(), "s2=grp::D2".to_string()];
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));
        run_demux(args).unwrap();

        let s1 = read_back(&dir.path().join("s1.fq"));
        assert_eq!(s1.len(), 1);
        assert_eq!(s1[0].name, b"a");
        let s2 = read_back(&dir.path().join("s2.fq"));
        assert_eq!(s2.len(), 1);
        assert_eq!(s2[0].name, b"b");
    }

    #[test]
    fn test_parallel_output_is_deterministic_across_thread_counts() {
        // Records are written in input order, so a single-file --out is
        // byte-identical regardless of --threads. The input spans several READ_CHUNK chunks so
        // cross-chunk ordering is exercised too.
        let dir = tempfile::tempdir().unwrap();
        let n = READ_CHUNK * 2 + 137; // > 2 chunks plus a partial tail
        let mut fq = Vec::new();
        for i in 0..n {
            let bc = match i % 3 {
                0 => "AAAAAAAA", // matches D1
                1 => "TTTTTTTT", // matches D2
                _ => "GGGGGGGG", // matches neither -> still passes through
            };
            fq.extend_from_slice(format!("@r{i}\n{bc}CGCG\n+\nIIIIIIIIIIII\n").as_bytes());
        }
        let input = write_file(dir.path(), "in.fq", &fq);

        let run = |threads: usize, out: &Path| {
            let mut args = base_args();
            args.threads = threads;
            args.inputs = vec![format!("0={}", input.display())];
            args.groups = vec![
                "grp={D1=AAAAAAAA,D2=TTTTTTTT}".to_string(),
                "grp::loc=0:0,dist=0".to_string(),
            ];
            args.out = Some(out.display().to_string());
            run_demux(args).unwrap();
        };

        let out_min = dir.path().join("min.fq");
        let out_many = dir.path().join("many.fq");
        run(3, &out_min);
        run(16, &out_many);
        assert_eq!(
            std::fs::read(&out_min).unwrap(),
            std::fs::read(&out_many).unwrap(),
            "single-file output must be identical across thread counts"
        );
    }

    #[test]
    fn test_parallel_bam_output_deterministic_across_thread_counts() {
        // Records are encoded on the workers and written in input order, so a
        // single-file BAM `--out` is byte-identical regardless of --threads.
        // The input spans several READ_CHUNK chunks so cross-chunk ordering and
        // multi-job parallel encode (with_min_len coarsens each chunk into
        // several jobs) are both exercised; the headers match because the test
        // sets no command line (no `@PG CL`).
        let dir = tempfile::tempdir().unwrap();
        let n = READ_CHUNK * 2 + 137; // > 2 chunks plus a partial tail
        let mut fq = Vec::new();
        for i in 0..n {
            let bc = match i % 3 {
                0 => "AAAAAAAA",
                1 => "TTTTTTTT",
                _ => "GGGGGGGG",
            };
            fq.extend_from_slice(format!("@r{i}\n{bc}CGCG\n+\nIIIIIIIIIIII\n").as_bytes());
        }
        let input = write_file(dir.path(), "in.fq", &fq);

        let run = |threads: usize, out: &Path| {
            let mut args = base_args();
            args.threads = threads;
            args.inputs = vec![format!("0={}", input.display())];
            args.groups = vec![
                "grp={D1=AAAAAAAA,D2=TTTTTTTT}".to_string(),
                "grp::loc=0:0,dist=0".to_string(),
            ];
            args.out = Some(out.display().to_string());
            run_demux(args).unwrap();
        };

        // Both thread counts must give identical bytes.
        let out_min = dir.path().join("min.bam");
        let out_many = dir.path().join("many.bam");
        run(4, &out_min);
        run(16, &out_many);
        assert_eq!(
            std::fs::read(&out_min).unwrap(),
            std::fs::read(&out_many).unwrap(),
            "single-file BAM output (encode-on-workers) must be identical across thread counts"
        );
    }

    #[test]
    fn test_sample_fanout_bam_has_read_group() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::header::record::value::map::read_group::tag as rgt;

        // A per-sample BAM lists its @RG (ID=label, SM, LB, plus shared
        // --rg-tag fields) and tags each record with that RG; %pool sets LB to
        // the pool id.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@a\nAAAAAAAA\n+\nIIIIIIII\n");
        let mut args = base_args();
        args.pool = Some("lib01".to_string());
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec!["grp={D1=AAAAAAAA}".to_string(), "grp::loc=0:0".to_string()];
        args.samples = vec!["dna01::%pool=grp::D1".to_string()];
        args.rg_tags = vec!["PL=ILLUMINA".to_string()];
        args.out = Some(format!("{}/%sample.%sub_sample.bam", dir.path().display()));
        run_demux(args).unwrap();

        let path = dir.path().join("dna01.lib01.bam");
        let mut reader = noodles::bam::io::Reader::new(File::open(&path).unwrap());
        let header = reader.read_header().unwrap();
        let read_group = &header.read_groups()[&b"dna01.lib01"[..]];
        let field = |tag| read_group.other_fields().get(tag).map(|v| v.to_string());
        assert_eq!(field(&rgt::SAMPLE).as_deref(), Some("dna01"));
        assert_eq!(field(&rgt::LIBRARY).as_deref(), Some("lib01"));
        assert_eq!(field(&rgt::PLATFORM).as_deref(), Some("ILLUMINA"));
        let record = reader.record_bufs(&header).next().unwrap().unwrap();
        match record.data().get(&Tag::READ_GROUP) {
            Some(Value::String(rg)) => assert_eq!(rg.to_string(), "dna01.lib01"),
            other => panic!("RG missing or wrong type: {other:?}"),
        }
    }

    #[test]
    fn test_remove_drops_and_optionally_writes() {
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@a\nAAAA\n+\nIIII\n@c\nCCCC\n+\nIIII\n",
        );
        let groups = vec![
            "grp={A=AAAA,C=CCCC}".to_string(),
            "grp::loc=0:0,dist=0".to_string(),
        ];

        // Bare `--remove grp::C` drops the C reads; A passes through.
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = groups.clone();
        args.remove = vec!["grp::C".to_string()];
        let kept = dir.path().join("kept.fq");
        args.out = Some(kept.display().to_string());
        run_demux(args).unwrap();
        let kept = read_back(&kept);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].name, b"a");

        // `--remove grp::C=PATTERN` also writes the removed reads to the
        // pattern.
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = groups.clone();
        let removed = dir.path().join("removed.fq");
        args.remove = vec![format!("grp::C={}", removed.display())];
        args.out = Some(dir.path().join("kept2.fq").display().to_string());
        run_demux(args).unwrap();
        let removed = read_back(&removed);
        assert_eq!(removed.len(), 1);
        assert_eq!(removed[0].name, b"c");

        // With --qc-tag, a removed record in a BAM bin carries the `removed`
        // slug naming the rule.
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = groups;
        let removed_bam = dir.path().join("removed.bam");
        args.remove = vec![format!("grp::C={}", removed_bam.display())];
        args.out = Some(dir.path().join("kept3.fq").display().to_string());
        args.qc_tag = Some("ZS".to_string());
        run_demux(args).unwrap();
        let record = first_bam_record(&removed_bam);
        assert_eq!(
            string_tag(&record, *b"ZS").as_deref(),
            Some(
                r#"{"v":1,"outcome":"removed","rule":"grp::C","group":"grp","tag":"CCCC","groups":[{"g":"grp","tag":"CCCC","loc":"0:0:4","obs":"CCCC","sub":0,"ind":0}]}"#
            )
        );
    }

    #[test]
    fn test_unassigned_fans_out_by_source() {
        // A required group that does not match leaves the read unassigned;
        // %source fans its two input segments into one file per input.
        let dir = tempfile::tempdir().unwrap();
        let in0 = write_file(dir.path(), "r1.fq", b"@r\nAAAA\n+\nIIII\n");
        let in1 = write_file(dir.path(), "r2.fq", b"@r\nCCCC\n+\nIIII\n");
        let mut args = base_args();
        args.inputs = vec![
            format!("0={}", in0.display()),
            format!("1={}", in1.display()),
        ];
        args.groups = vec![
            "g={GGGG}".to_string(),
            "g::loc=0:0,minFindsPerGroup=1".to_string(),
        ];
        args.unassigned = Some(format!("{}/un.%source.fq", dir.path().display()));
        args.out = Some(format!("{}/out.fq", dir.path().display()));
        run_demux(args).unwrap();

        let s0 = read_back(&dir.path().join("un.0.fq"));
        assert_eq!(s0.len(), 1);
        assert_eq!(s0[0].bases, b"AAAA");
        let s1 = read_back(&dir.path().join("un.1.fq"));
        assert_eq!(s1.len(), 1);
        assert_eq!(s1[0].bases, b"CCCC");
    }

    #[test]
    fn test_tag_cb_written_to_bam() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;

        // --tag CB=bc writes the carved barcode as a CB:Z field; the body is
        // the carved template.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@a\nAAAACCCC\n+\nIIIIIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.extracts = vec!["bc=0:0:4".to_string(), "body=0:4:end".to_string()];
        args.templates = vec!["body".to_string()];
        args.tags = vec!["CB=bc".to_string(), "CB::qual=none".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let mut reader = noodles::bam::io::Reader::new(File::open(&out).unwrap());
        let header = reader.read_header().unwrap();
        let record = reader.record_bufs(&header).next().unwrap().unwrap();
        let sequence: &[u8] = record.sequence().as_ref();
        assert_eq!(sequence, b"CCCC", "body is the carved template");
        match record.data().get(&Tag::new(b'C', b'B')) {
            Some(Value::String(cb)) => assert_eq!(cb.to_string(), "AAAA"),
            other => panic!("CB missing or wrong type: {other:?}"),
        }
    }

    #[test]
    fn test_sample_from_group_fans_out_one_per_tag() {
        // --sample-from-group makes one sample per tag id (a 1:1 shape), no
        // --sample listing.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@a\nAAAA\n+\nIIII\n@b\nCCCC\n+\nIIII\n",
        );
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "sample_bc={s1=AAAA,s2=CCCC}".to_string(),
            "sample_bc::loc=0:0,dist=0".to_string(),
        ];
        args.sample_from_group = Some("sample_bc".to_string());
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));
        run_demux(args).unwrap();

        let s1 = read_back(&dir.path().join("s1.fq"));
        assert_eq!(s1.len(), 1);
        assert_eq!(s1[0].name, b"a");
        let s2 = read_back(&dir.path().join("s2.fq"));
        assert_eq!(s2.len(), 1);
        assert_eq!(s2[0].name, b"b");
    }

    #[test]
    fn test_sample_sheet_fans_out_through_engine() {
        // A --sample-sheet resolves each barcode cell against its group and
        // fans out per sample.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@a\nAAAA\n+\nIIII\n@b\nCCCC\n+\nIIII\n",
        );
        let sheet = write_file(
            dir.path(),
            "sheet.tsv",
            b"sample\tgrp\ndna01\tA\ndna02\tB\n",
        );
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={A=AAAA,B=CCCC}".to_string(),
            "grp::loc=0:0,dist=0".to_string(),
        ];
        args.sample_sheet = Some(sheet);
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));
        run_demux(args).unwrap();

        assert_eq!(read_back(&dir.path().join("dna01.fq"))[0].name, b"a");
        assert_eq!(read_back(&dir.path().join("dna02.fq"))[0].name, b"b");
    }

    #[test]
    fn test_sample_sheet_ignored_column_fans_out_through_engine() {
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@a\nAAAA\n+\nIIII\n@b\nCCCC\n+\nIIII\n",
        );
        let sheet = write_file(
            dir.path(),
            "sheet.tsv",
            b"sample\tassay\tgrp\ndna01\tscATAC\tA\ndna02\tanything at all\tB\n",
        );
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={A=AAAA,B=CCCC}".to_string(),
            "grp::loc=0:0,dist=0".to_string(),
        ];
        args.sample_sheet = Some(sheet);
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));

        let err = run_demux(args.clone()).unwrap_err();
        assert!(format!("{err:#}").contains("`assay`"), "{err:#}");

        args.sample_sheet_ignore = vec!["assay".to_string()];
        run_demux(args).unwrap();
        assert_eq!(read_back(&dir.path().join("dna01.fq"))[0].name, b"a");
        assert_eq!(read_back(&dir.path().join("dna02.fq"))[0].name, b"b");
    }

    #[test]
    fn test_multi_read_fastx_fans_out_by_ordinal() {
        // Two templates to FASTX need %ordinal; each output read goes to its
        // R%ordinal file.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@a\nAAAACCCC\n+\nIIIIIIII\n");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.extracts = vec!["r1=0:0:4".to_string(), "r2=0:4:end".to_string()];
        args.templates = vec!["r1".to_string(), "r2".to_string()];
        args.out = Some(format!("{}/out.R%ordinal.fq", dir.path().display()));
        run_demux(args).unwrap();

        let r1 = read_back(&dir.path().join("out.R1.fq"));
        assert_eq!(r1.len(), 1);
        assert_eq!(r1[0].bases, b"AAAA");
        let r2 = read_back(&dir.path().join("out.R2.fq"));
        assert_eq!(r2.len(), 1);
        assert_eq!(r2[0].bases, b"CCCC");
    }

    #[test]
    fn test_multi_read_fastx_without_ordinal_is_error() {
        // Two templates to a single FASTX file with no %ordinal cannot be
        // separated; fail fast.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@a\nAAAACCCC\n+\nIIIIIIII\n");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.extracts = vec!["r1=0:0:4".to_string(), "r2=0:4:end".to_string()];
        args.templates = vec!["r1".to_string(), "r2".to_string()];
        args.out = Some(dir.path().join("out.fq").display().to_string());
        let err = run_demux(args).unwrap_err();
        assert!(err.to_string().contains("%ordinal"), "{err}");
    }

    #[test]
    fn test_unassigned_interleaved_single_input_is_one_source() {
        // A single auto-detected interleaved input is one %source (0); its
        // mates stay together in un.0.fq rather than being split across
        // un.0/un.1.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@r/1\nAAAA\n+\nIIII\n@r/2\nTTTT\n+\nIIII\n",
        );
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "g={GGGG}".to_string(),
            "g::loc=0:0,minFindsPerGroup=1".to_string(),
        ];
        args.unassigned = Some(format!("{}/un.%source.fq", dir.path().display()));
        run_demux(args).unwrap();

        let both = read_back(&dir.path().join("un.0.fq"));
        assert_eq!(both.len(), 2, "both mates land in the single source-0 file");
        assert_eq!(both[0].bases, b"AAAA");
        assert_eq!(both[1].bases, b"TTTT");
        assert!(
            !dir.path().join("un.1.fq").exists(),
            "an interleaved single input has no source 1"
        );
    }

    #[test]
    fn test_single_file_bam_lists_every_read_group() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;

        // A single-file --out (no %sample) holds every assigned read, separable
        // by @RG: the header lists one @RG per sample and each record carries
        // its own RG.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@a\nAAAA\n+\nIIII\n@b\nTTTT\n+\nIIII\n",
        );
        let out = dir.path().join("pool.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={D1=AAAA,D2=TTTT}".to_string(),
            "grp::loc=0:0,dist=0".to_string(),
        ];
        args.samples = vec!["s1=grp::D1".to_string(), "s2=grp::D2".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let mut reader = noodles::bam::io::Reader::new(File::open(&out).unwrap());
        let header = reader.read_header().unwrap();
        assert!(header.read_groups().contains_key(&b"s1"[..]), "@RG s1");
        assert!(header.read_groups().contains_key(&b"s2"[..]), "@RG s2");
        let records: Vec<_> = reader
            .record_bufs(&header)
            .collect::<std::io::Result<_>>()
            .unwrap();
        let rg =
            |record: &noodles::sam::alignment::RecordBuf| match record.data().get(&Tag::READ_GROUP)
            {
                Some(Value::String(value)) => value.to_string(),
                _ => panic!("record has no RG"),
            };
        assert_eq!(rg(&records[0]), "s1");
        assert_eq!(rg(&records[1]), "s2");
    }

    #[test]
    fn test_passthrough_bam_carries_default_pool_read_group() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::header::record::value::map::read_group::tag as rgt;

        // A pass-through (no --sample) BAM still carries a read group so the
        // uBAM is valid downstream: a single default @RG whose ID/SM/LB are the
        // pool id, and every record references it via RG:Z:<pool>.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@a\nACGT\n+\nIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.pool = Some("lib01".to_string());
        args.inputs = vec![format!("0={}", input.display())];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let mut reader = noodles::bam::io::Reader::new(File::open(&out).unwrap());
        let header = reader.read_header().unwrap();
        let group = header
            .read_groups()
            .get(&b"lib01"[..])
            .expect("@RG ID=pool");
        let field = |tag| group.other_fields().get(tag).map(|v| v.to_string());
        assert_eq!(field(&rgt::SAMPLE).as_deref(), Some("lib01"), "SM=pool");
        assert_eq!(field(&rgt::LIBRARY).as_deref(), Some("lib01"), "LB=pool");
        let record = reader.record_bufs(&header).next().unwrap().unwrap();
        match record.data().get(&Tag::READ_GROUP) {
            Some(Value::String(value)) => {
                assert_eq!(value.to_string(), "lib01", "read RG:Z:<pool>")
            }
            _ => panic!("pass-through read has no RG"),
        }
    }

    #[test]
    fn test_bam_output_carries_pg_provenance() {
        use noodles::sam::header::record::value::map::program::tag as pg;

        // Every BAM gets a single @PG with ID/PN unmux, the running version,
        // and the command line.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nACGT\n+\nIIII\n");
        let out = dir.path().join("out.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.out = Some(out.display().to_string());
        args.command_line = Some("unmux --in 0=in.fq --out out.bam".to_string());
        run_demux(args).unwrap();

        let mut reader = noodles::bam::io::Reader::new(File::open(&out).unwrap());
        let header = reader.read_header().unwrap();
        let programs = header.programs().as_ref();
        assert_eq!(programs.len(), 1, "exactly one @PG");
        let program = programs.get(&b"unmux"[..]).expect("@PG unmux");
        let field = |tag| program.other_fields().get(tag).map(|v| v.to_string());
        assert_eq!(field(&pg::NAME).as_deref(), Some("unmux"));
        assert_eq!(
            field(&pg::VERSION).as_deref(),
            Some(env!("CARGO_PKG_VERSION"))
        );
        assert!(
            field(&pg::COMMAND_LINE).unwrap().contains("--out out.bam"),
            "CL carries the command line"
        );
    }

    #[test]
    fn test_cram_output_round_trips_with_provenance() {
        // Reference-free CRAM fan-out: each per-sample CRAM is readable and
        // carries @PG and its @RG, with the record sequence intact.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@a\nAAAA\n+\nIIII\n@b\nTTTT\n+\nIIII\n",
        );
        let mut args = base_args();
        args.pool = Some("pool1".to_string());
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={s1=AAAA,s2=TTTT}".to_string(),
            "grp::loc=0:0:4,dist=0".to_string(),
        ];
        args.samples = vec!["s1=grp::s1".to_string(), "s2=grp::s2".to_string()];
        args.out = Some(format!("{}/%sample.cram", dir.path().display()));
        args.command_line = Some("unmux demux --out %sample.cram".to_string());
        run_demux(args).unwrap();

        let path = dir.path().join("s1.cram");
        let mut reader = noodles::cram::io::Reader::new(File::open(&path).unwrap());
        let header = reader.read_header().unwrap();
        assert!(
            header.programs().as_ref().contains_key(&b"unmux"[..]),
            "@PG present"
        );
        assert!(header.read_groups().contains_key(&b"s1"[..]), "@RG s1");
        let records: Vec<_> = reader
            .records(&header)
            .collect::<std::io::Result<_>>()
            .unwrap();
        assert_eq!(records.len(), 1);
        let sequence: &[u8] = records[0].sequence().as_ref();
        assert_eq!(sequence, b"AAAA");
    }

    #[test]
    fn test_cram_paired_with_tag_and_rg_round_trips() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;

        // A read pair (two templates) carved to CRAM with a CB tag and a
        // sample: the pair flags, RG, CB, sequences, and qualities are all
        // preserved on round-trip.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@r\nAAAACCCCGGGG\n+\nIIIIIIIIIIII\n");
        let out = dir.path().join("out.cram");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={s1=AAAA}".to_string(),
            "grp::loc=0:0:4,dist=0".to_string(),
        ];
        args.extracts = vec![
            "bc=0:0:4".to_string(),
            "r1=0:4:8".to_string(),
            "r2=0:8:12".to_string(),
        ];
        args.templates = vec!["r1".to_string(), "r2".to_string()];
        args.tags = vec!["CB=bc".to_string(), "CB::qual=none".to_string()];
        args.samples = vec!["s1=grp::s1".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let mut reader = noodles::cram::io::Reader::new(File::open(&out).unwrap());
        let header = reader.read_header().unwrap();
        let records: Vec<_> = reader
            .records(&header)
            .collect::<std::io::Result<_>>()
            .unwrap();
        assert_eq!(records.len(), 2, "a read pair");
        assert!(records[0].flags().is_segmented() && records[0].flags().is_first_segment());
        assert!(records[1].flags().is_segmented() && records[1].flags().is_last_segment());
        let seq0: &[u8] = records[0].sequence().as_ref();
        let seq1: &[u8] = records[1].sequence().as_ref();
        assert_eq!(seq0, b"CCCC");
        assert_eq!(seq1, b"GGGG");
        assert_eq!(
            records[0].quality_scores().as_ref().len(),
            4,
            "qualities preserved"
        );
        let data = records[0].data();
        match data.get(&Tag::READ_GROUP) {
            Some(Value::String(rg)) => assert_eq!(rg.to_string(), "s1"),
            other => panic!("RG missing: {other:?}"),
        }
        match data.get(&Tag::new(b'C', b'B')) {
            Some(Value::String(cb)) => assert_eq!(cb.to_string(), "AAAA"),
            other => panic!("CB missing: {other:?}"),
        }
    }

    #[test]
    fn test_directed_files_exist_even_when_empty() {
        // Every directed output FILE exists even with no reads: a declared
        // `%sample.fq` that gets no reads is an empty (0-byte) file, and a
        // directed `--unassigned` bin that gets no reads exists too. (stdout is
        // the only directed sink not pre-created.)
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@a\nAAAA\n+\nIIII\n");
        let unassigned = dir.path().join("un.fq");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={s1=AAAA,s2=TTTT}".to_string(),
            "grp::loc=0:0:4,dist=0".to_string(),
        ];
        args.samples = vec!["s1=grp::s1".to_string(), "s2=grp::s2".to_string()];
        args.out = Some(format!("{}/%sample.fq", dir.path().display()));
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        // s1 received the only read; s2 received none but its file still
        // exists, empty.
        assert_eq!(read_back(&dir.path().join("s1.fq")).len(), 1);
        let s2 = dir.path().join("s2.fq");
        assert!(s2.exists(), "a directed zero-read sample file exists");
        assert_eq!(std::fs::metadata(&s2).unwrap().len(), 0, "and is empty");
        // No read was unassigned, but the directed bin still exists (empty).
        assert!(
            unassigned.exists(),
            "a directed --unassigned bin exists even when empty"
        );
        assert!(read_back(&unassigned).is_empty());
    }

    #[test]
    fn test_zero_read_sample_creates_empty_bam() {
        // A declared sample that receives no reads still gets its directed
        // output file: an empty but valid BAM (header + EOF, zero records).
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(dir.path(), "in.fq", b"@a\nAAAA\n+\nIIII\n");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={s1=AAAA,s2=TTTT}".to_string(),
            "grp::loc=0:0:4,dist=0".to_string(),
        ];
        // s2 matches no read in this input.
        args.samples = vec!["s1=grp::s1".to_string(), "s2=grp::s2".to_string()];
        args.out = Some(format!("{}/%sample.bam", dir.path().display()));
        run_demux(args).unwrap();

        assert!(dir.path().join("s1.bam").exists());
        let s2 = dir.path().join("s2.bam");
        assert!(
            s2.exists(),
            "a directed but zero-read sample still gets a file"
        );

        // It is a valid, empty BAM: @PG header present, zero records.
        let mut reader = noodles::bam::io::Reader::new(File::open(&s2).unwrap());
        let header = reader.read_header().unwrap();
        assert!(
            header.programs().as_ref().contains_key(&b"unmux"[..]),
            "@PG present"
        );
        assert_eq!(reader.record_bufs(&header).count(), 0, "no records");
    }

    #[test]
    fn test_metrics_counts_removed_reads() {
        // --remove reads are tallied in the summary `reads_removed` column,
        // distinct from unassigned.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.fq",
            b"@a\nAAAA\n+\nIIII\n@c\nCCCC\n+\nIIII\n",
        );
        let summary = dir.path().join("summary.tsv");
        let mut args = base_args();
        args.pool = Some("p".to_string());
        args.inputs = vec![format!("0={}", input.display())];
        args.groups = vec![
            "grp={A=AAAA,C=CCCC}".to_string(),
            "grp::loc=0:0:4,dist=0".to_string(),
        ];
        args.remove = vec!["grp::C".to_string()];
        args.out = Some(dir.path().join("out.fq").display().to_string());
        args.metrics_summary = Some(summary.clone());
        run_demux(args).unwrap();

        let sum = std::fs::read_to_string(&summary).unwrap();
        // total 2; A passes through (assigned 0, pass_through 1); C removed
        // (removed 1); unassigned 0.
        let row = sum.lines().nth(1).unwrap();
        assert!(
            row.starts_with("p\t2\t0\t1\t0\t1\t"),
            "removed tallied: {sum}"
        );
    }

    #[test]
    fn test_passthrough_merges_input_tag_with_demux_tag() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;

        // A pass-through SAM read keeps its carried input tag (ZX) and also
        // gains the assembled --tag field (BC); --tag fields and carried tags
        // coexist on the record.
        let dir = tempfile::tempdir().unwrap();
        let input = write_file(
            dir.path(),
            "in.sam",
            b"@HD\tVN:1.6\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tZX:Z:keep\n",
        );
        let out = dir.path().join("out.sam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", input.display())];
        args.extracts = vec!["bc=0:0:4".to_string()];
        args.tags = vec!["BC=bc".to_string(), "BC::qual=none".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let mut reader =
            noodles::sam::io::Reader::new(std::io::BufReader::new(File::open(&out).unwrap()));
        let header = reader.read_header().unwrap();
        let record = reader.record_bufs(&header).next().unwrap().unwrap();
        let value = |tag| match record.data().get(&tag) {
            Some(Value::String(v)) => Some(v.to_string()),
            _ => None,
        };
        assert_eq!(value(Tag::new(b'Z', b'X')).as_deref(), Some("keep"));
        assert_eq!(value(Tag::new(b'B', b'C')).as_deref(), Some("ACGT"));
    }

    #[test]
    fn test_zero_read_run_creates_all_directed_files() {
        // The sequencer produced no data: both inputs are empty. Every directed
        // output (the per-sample fan-out file and the unassigned bin) must
        // still be created as a valid empty BAM (header + EOF), so a no-data
        // run never fails a pipeline. read_back decoding them to zero records,
        // plus a non-zero file size, proves they are valid containers and not
        // 0-byte stubs.
        let dir = tempfile::tempdir().unwrap();
        let i7 = write_file(dir.path(), "i7.fq", b""); // 0-byte FASTQ, sniffed by extension
        let i5 = write_file(dir.path(), "i5.fq", b"");
        let unassigned = dir.path().join("un.bam");
        let mut args = base_args();
        args.inputs = vec![format!("0={}", i7.display()), format!("1={}", i5.display())];
        args.extracts = vec![
            "i7=0:0:8".to_string(),
            "i5=1:0:8".to_string(),
            "body=0:0:end".to_string(),
        ];
        args.templates = vec!["body".to_string()];
        args.groups = vec![
            "sample_bc={s1=AAAAAAAACCCCCCCC}".to_string(),
            "sample_bc::match=i7+i5,dist=0".to_string(),
        ];
        args.sample_from_group = Some("sample_bc".to_string());
        args.out = Some(format!("{}/%sample.bam", dir.path().display()));
        args.unassigned = Some(unassigned.display().to_string());
        run_demux(args).unwrap();

        let s1 = dir.path().join("s1.bam");
        assert!(
            s1.exists(),
            "per-sample file is created even with zero input reads"
        );
        assert!(
            unassigned.exists(),
            "unassigned bin is created even with zero input reads"
        );
        // Valid empty BAM: non-zero (header + EOF) and decodes to zero records.
        assert!(
            std::fs::metadata(&s1).unwrap().len() > 0,
            "s1.bam has a header + EOF, not 0 bytes"
        );
        assert!(read_back(&s1).is_empty(), "s1.bam holds zero records");
        assert!(
            std::fs::metadata(&unassigned).unwrap().len() > 0,
            "un.bam has a header + EOF"
        );
        assert!(
            read_back(&unassigned).is_empty(),
            "un.bam holds zero records"
        );
    }

    #[test]
    fn match_read_group_routes_by_rg_tag() {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;

        let rg_to_idx: HashMap<Vec<u8>, usize> = [(b"rg1".to_vec(), 0), (b"rg2".to_vec(), 1)]
            .into_iter()
            .collect();

        let with_rg = |rg: &[u8]| InputRecord {
            name: b"r".to_vec(),
            bases: b"ACGT".to_vec(),
            quals: None,
            tags: Some(
                [(Tag::READ_GROUP, Value::String(rg.to_vec().into()))]
                    .into_iter()
                    .collect(),
            ),
            qc_fail: false,
            alignment: None,
        };

        let fragment = Fragment {
            records: vec![with_rg(b"rg2")],
        };
        assert_eq!(
            match_read_group(&rg_to_idx, &fragment),
            Some(1),
            "routes to the tag index of the record's RG:Z"
        );

        let no_rg = InputRecord {
            name: b"r".to_vec(),
            bases: b"ACGT".to_vec(),
            quals: None,
            tags: None,
            qc_fail: false,
            alignment: None,
        };
        let fragment = Fragment {
            records: vec![no_rg],
        };
        assert_eq!(
            match_read_group(&rg_to_idx, &fragment),
            None,
            "no RG tag at all means no route"
        );

        let fragment = Fragment {
            records: vec![with_rg(b"rg3")],
        };
        assert_eq!(
            match_read_group(&rg_to_idx, &fragment),
            None,
            "an RG id absent from the map means no route"
        );
    }

    #[test]
    fn test_passes_alignments_through_only_for_an_aligned_split_into_alignments() {
        let mut args = base_args();
        args.inputs = vec!["0=in.bam".to_string()];
        let plan = crate::grammar::parse_demux(&args).unwrap();
        let bam = [OutputFormat::Bam];
        assert!(passes_alignments_through(&plan, true, &bam).unwrap());
        assert!(!passes_alignments_through(&plan, false, &bam).unwrap());
        assert!(!passes_alignments_through(&plan, true, &[]).unwrap());
    }

    #[test]
    fn test_passes_alignments_through_refuses_a_template() {
        let mut args = base_args();
        args.inputs = vec!["0=in.bam".to_string()];
        args.extracts = vec!["head=0:0:4".to_string()];
        args.templates = vec!["head".to_string()];
        let plan = crate::grammar::parse_demux(&args).unwrap();
        let err = passes_alignments_through(&plan, true, &[OutputFormat::Sam])
            .unwrap_err()
            .to_string();
        assert!(err.contains("--template"), "{err}");
        assert!(!passes_alignments_through(&plan, true, &[]).unwrap());
    }

    #[test]
    fn test_passes_alignments_through_refuses_several_inputs() {
        let mut args = base_args();
        args.inputs = vec!["0=a.bam".to_string(), "1=b.bam".to_string()];
        let plan = crate::grammar::parse_demux(&args).unwrap();
        let err = passes_alignments_through(&plan, true, &[OutputFormat::Bam])
            .unwrap_err()
            .to_string();
        assert!(err.contains("single input"), "{err}");
    }

    #[test]
    fn test_passes_alignments_through_refuses_cram() {
        let mut args = base_args();
        args.inputs = vec!["0=in.bam".to_string()];
        args.reference = Some(PathBuf::from("ref.fa"));
        let plan = crate::grammar::parse_demux(&args).unwrap();
        let err = passes_alignments_through(&plan, true, &[OutputFormat::Bam, OutputFormat::Cram])
            .unwrap_err()
            .to_string();
        assert!(err.contains("CRAM is not supported"), "{err}");
    }

    #[test]
    fn test_passes_alignments_through_refuses_tags_it_would_overwrite() {
        for tag in ["NM", "MD", "MC", "SA"] {
            let mut args = base_args();
            args.inputs = vec!["0=in.bam".to_string()];
            args.extracts = vec!["head=0:0:4".to_string()];
            args.tags = vec![format!("{tag}=head")];
            let plan = crate::grammar::parse_demux(&args).unwrap();
            let err = passes_alignments_through(&plan, true, &[OutputFormat::Bam])
                .unwrap_err()
                .to_string();
            assert!(
                err.contains(&format!("`{tag}` would overwrite an alignment tag")),
                "{err}"
            );
        }
    }

    #[test]
    fn test_passes_alignments_through_refuses_a_tag_drawn_from_sequence() {
        let mut args = base_args();
        args.inputs = vec!["0=in.bam".to_string()];
        args.extracts = vec!["head=0:0:4".to_string()];
        args.tags = vec!["XH=head".to_string()];
        let plan = crate::grammar::parse_demux(&args).unwrap();
        let err = passes_alignments_through(&plan, true, &[OutputFormat::Bam])
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("--tag `XH` is drawn from record sequences"),
            "{err}"
        );
        assert!(!passes_alignments_through(&plan, true, &[]).unwrap());
    }

    #[test]
    fn test_passes_alignments_through_refuses_a_sequence_group() {
        let mut args = base_args();
        args.inputs = vec!["0=in.bam".to_string()];
        args.groups = vec![
            "rg=@RG::SM".to_string(),
            "ad={X=AACCGGTT}".to_string(),
            "ad::loc=0:0:8".to_string(),
        ];
        args.sample_from_group = Some("rg".to_string());
        args.remove = vec!["ad".to_string()];
        let plan = crate::grammar::parse_demux(&args).unwrap();
        let err = passes_alignments_through(&plan, true, &[OutputFormat::Sam])
            .unwrap_err()
            .to_string();
        assert!(err.contains("group `ad` matches record sequences"), "{err}");
        assert!(!passes_alignments_through(&plan, true, &[]).unwrap());
    }

    #[test]
    fn test_non_primary_fragments_are_detected() {
        use noodles::sam::alignment::record::Flags;
        let fragment = |flags: Flags| Fragment {
            records: vec![InputRecord {
                name: b"r".to_vec(),
                bases: Vec::new(),
                quals: None,
                tags: None,
                qc_fail: false,
                alignment: Some(Alignment::Record(
                    RecordBuf::builder().set_flags(flags).build(),
                )),
            }],
        };
        assert!(is_non_primary(&fragment(Flags::SECONDARY)));
        assert!(is_non_primary(&fragment(Flags::SUPPLEMENTARY)));
        assert!(!is_non_primary(&fragment(Flags::empty())));
        let mut unkept = fragment(Flags::SECONDARY);
        unkept.records[0].alignment = None;
        assert!(!is_non_primary(&unkept));
    }

    #[test]
    fn test_column_streams_are_identical_across_threads() {
        let dir = tempfile::tempdir().unwrap();
        let tags = write_file(
            dir.path(),
            "tags.tsv",
            b"id\tseq\tcell\nw1\tAACCGG\tA+B\nw2\tTTGGCC\t\n",
        );
        let mut fq = Vec::new();
        for i in 0..3000 {
            let bc = ["AACCGG", "AACCGT", "TTGGCC", "GAGAGA"][i % 4];
            fq.extend_from_slice(format!("@r{i}\n{bc}ACGTACGT\n+\nABCDEFGHIJKLMN\n").as_bytes());
        }
        let input = write_file(dir.path(), "in.fq", &fq);
        let run = |threads: usize, ext: &str| {
            let out = dir.path().join(format!("{threads}.{ext}"));
            let mut args = base_args();
            args.threads = threads;
            args.inputs = vec![format!("0={}", input.display())];
            args.groups = vec![
                format!("bc={}", tags.display()),
                "bc::loc=0:0:6,dist=1".to_string(),
            ];
            args.extracts = [
                "b=@bc",
                "cell=@bc::cell",
                "i=@bc::id",
                "u=0:6:10",
                "r=0:10:end",
            ]
            .map(String::from)
            .to_vec();
            args.templates = vec!["r".to_string()];
            args.tags = ["BC=b", "CB=cell+u", "xi=i"].map(String::from).to_vec();
            args.pool = Some("p".to_string());
            args.out = Some(out.display().to_string());
            run_demux(args).unwrap();
            out
        };
        for ext in ["sam", "bam", "cram", "fa", "fq", "fq.gz"] {
            let serial = read_back(&run(1, ext));
            assert_eq!(serial.len(), 3000);
            for threads in [2, 3, 8] {
                assert!(
                    serial == read_back(&run(threads, ext)),
                    "{ext} with {threads} threads"
                );
            }
        }
        let records = read_back(&dir.path().join("1.fq"));
        let text = |record: &InputRecord, tag: [u8; 2]| match record.tags.as_ref()?.get(&tag) {
            Some(Value::String(value)) => Some(value.to_string()),
            _ => None,
        };
        for (i, record) in records.iter().enumerate() {
            let (cell, id) = [
                (Some("A+B-ACGT"), Some("w1")),
                (Some("A+B-ACGT"), Some("w1")),
                (Some("TTGGCC-ACGT"), Some("w2")),
                (None, None),
            ][i % 4];
            assert_eq!(text(record, *b"CB").as_deref(), cell, "record {i}");
            assert_eq!(text(record, *b"xi").as_deref(), id, "record {i}");
            assert_eq!(text(record, *b"CY"), None, "record {i}");
        }
    }

    #[test]
    fn test_sample_columns_are_per_row_and_identical_across_threads() {
        let dir = tempfile::tempdir().unwrap();
        let sheet = write_file(
            dir.path(),
            "sheet.tsv",
            b"sample\ta\tb\tnucleus\ns1\ta1\tb1\tN1\ns1\ta1\tb2\tN2+N3\ns2\ta2\tb1\t\ns2\ta2\tb2\t007\n",
        );
        let mut fq = Vec::new();
        for i in 0..3000 {
            let bc = ["AAAAGGGG", "AAAATTTT", "CCCCGGGG", "CCCCTTTT", "GAGAGGGG"][i % 5];
            fq.extend_from_slice(format!("@r{i}\n{bc}ACGTACGT\n+\nABCDEFGHIJKLMNOP\n").as_bytes());
        }
        let input = write_file(dir.path(), "in.fq", &fq);
        let run = |threads: usize, ext: &str| {
            let out = dir.path().join(format!("{threads}.{ext}"));
            let mut args = base_args();
            args.threads = threads;
            args.inputs = vec![format!("0={}", input.display())];
            args.groups = [
                "a={a1=AAAA,a2=CCCC}",
                "a::loc=0:0:4",
                "b={b1=GGGG,b2=TTTT}",
                "b::loc=0:4:8",
            ]
            .map(String::from)
            .to_vec();
            args.extracts = [
                "nuc=%sample::nucleus",
                "sm=%sample::sample",
                "ba=@a",
                "bb=@b",
                "r=0:8:end",
            ]
            .map(String::from)
            .to_vec();
            args.templates = vec!["r".to_string()];
            args.tags = ["nu=nuc+ba+bb", "xn=nuc", "xs=sm", "CB=ba+bb"]
                .map(String::from)
                .to_vec();
            args.sample_sheet = Some(sheet.clone());
            args.pool = Some("p".to_string());
            args.out = Some(out.display().to_string());
            run_demux(args).unwrap();
            out
        };
        for ext in ["sam", "bam", "cram", "fq", "fq.gz"] {
            let serial = read_back(&run(1, ext));
            assert_eq!(serial.len(), 2400);
            for threads in [2, 8] {
                assert!(
                    serial == read_back(&run(threads, ext)),
                    "{ext} with {threads} threads"
                );
            }
        }
        let records = read_back(&dir.path().join("1.sam"));
        let text = |record: &InputRecord, tag: [u8; 2]| match record.tags.as_ref()?.get(&tag) {
            Some(Value::String(value)) => Some(value.to_string()),
            _ => None,
        };
        for (i, record) in records.iter().enumerate() {
            let (nu, xn, xs, cb) = [
                (Some("N1-AAAA-GGGG"), Some("N1"), "s1", "AAAA-GGGG"),
                (Some("N2+N3-AAAA-TTTT"), Some("N2+N3"), "s1", "AAAA-TTTT"),
                (None, None, "s2", "CCCC-GGGG"),
                (Some("007-CCCC-TTTT"), Some("007"), "s2", "CCCC-TTTT"),
            ][i % 4];
            assert_eq!(text(record, *b"nu").as_deref(), nu, "record {i}");
            assert_eq!(text(record, *b"xn").as_deref(), xn, "record {i}");
            assert_eq!(text(record, *b"xs").as_deref(), Some(xs), "record {i}");
            assert_eq!(text(record, *b"CB").as_deref(), Some(cb), "record {i}");
            assert_eq!(
                text(record, *b"CY").as_deref(),
                Some("ABCD EFGH"),
                "record {i}"
            );
        }
    }

    #[test]
    fn test_merge_data_drops_only_a_stale_default_quality_tag() {
        let data = |tags: &[&str]| -> Data {
            tags.iter()
                .map(|tag| (Tag::from(tag_bytes(tag)), Value::from(tag.to_string())))
                .collect()
        };
        let merged = |carried: &[&str], written: &[&str], columns: &[(&str, &str)]| {
            merge_data(Some(&data(carried)), &data(written), columns, &[])
                .iter()
                .map(|(tag, _)| String::from_utf8(<[u8; 2]>::from(tag).to_vec()).unwrap())
                .collect::<Vec<_>>()
                .join(" ")
        };
        let cb = &[("CB", "CY")];
        assert_eq!(merged(&["CB", "CY", "NH", "XX"], &["CB"], cb), "CB NH XX");
        assert_eq!(merged(&["CB", "CY"], &["CB"], &[]), "CB CY");
        assert_eq!(merged(&["CB", "CY"], &[], cb), "CB CY");
        assert_eq!(merged(&["CB", "CY"], &["CB", "CY"], cb), "CB CY");
        assert_eq!(merged(&["CY"], &["CB"], cb), "CB");
        assert_eq!(merged(&["CR", "CB", "CY"], &["CB"], cb), "CR CB CY");
        assert_eq!(merged(&["CR", "CB", "CY"], &["CB", "CR"], cb), "CR CB");
        let cr = &[("CR", "CY")];
        assert_eq!(merged(&["CB", "CY"], &["CR"], cr), "CB CR");
        assert_eq!(merged(&["CR", "CY"], &["CR"], cr), "CR");
        assert_eq!(merged(&["CR", "CB", "CY"], &["CR"], cr), "CR CB");
        assert_eq!(merged(&["CR", "CB", "CY"], &["CB"], cr), "CR CB CY");
        let both = &[("CB", "CY"), ("CR", "CY")];
        assert_eq!(merged(&["CR", "CB", "CY"], &["CB", "CR"], both), "CR CB");
        assert_eq!(merged(&["CR", "CB", "CY"], &["CB"], both), "CR CB CY");
        assert_eq!(merged(&["CB", "CY", "QX"], &["CB"], cb), "CB QX");
    }

    #[test]
    fn test_merge_data_drops_a_carried_tag_whose_binding_is_left_out() {
        let data = |tags: &[&str]| -> Data {
            tags.iter()
                .map(|tag| (Tag::from(tag_bytes(tag)), Value::from(tag.to_string())))
                .collect()
        };
        let merged = |carried: &[&str], written: &[&str], left_out: &[(&str, Option<&str>)]| {
            merge_data(Some(&data(carried)), &data(written), &[], left_out)
                .iter()
                .map(|(tag, _)| String::from_utf8(<[u8; 2]>::from(tag).to_vec()).unwrap())
                .collect::<Vec<_>>()
                .join(" ")
        };
        let cb = &[("CB", None)];
        assert_eq!(merged(&["CB", "CY", "RX", "QX"], &[], cb), "RX QX");
        assert_eq!(
            merged(&["CB", "CY", "RX", "QX"], &["RX", "QX"], cb),
            "RX QX"
        );
        assert_eq!(merged(&["RX", "QX"], &[], cb), "RX QX");
        assert_eq!(merged(&["CR", "CB", "CY"], &[], cb), "CR CY");
        assert_eq!(
            merged(&["CR", "CB", "CY"], &[], &[("CB", None), ("CR", None)]),
            ""
        );
        assert_eq!(merged(&["CB", "CY", "XY"], &[], &[("CB", Some("XY"))]), "");
        assert_eq!(
            merged(&["xc", "BZ", "XX"], &[], &[("xc", Some("BZ"))]),
            "XX"
        );
    }

    #[test]
    fn test_merge_data_keeps_no_carried_quality_tag_for_a_rewritten_tag() {
        let data = |tags: &[&str]| -> Data {
            tags.iter()
                .map(|tag| (Tag::from(tag_bytes(tag)), Value::from(tag.to_string())))
                .collect()
        };
        let merged = |carried: &[&str], written: &[&str], left_out: &[(&str, Option<&str>)]| {
            merge_data(Some(&data(carried)), &data(written), &[], left_out)
                .iter()
                .map(|(tag, _)| String::from_utf8(<[u8; 2]>::from(tag).to_vec()).unwrap())
                .collect::<Vec<_>>()
                .join(" ")
        };
        let cr = &[("CR", None)];
        assert_eq!(merged(&["CB", "CR", "CY"], &["CB"], cr), "CB");
        assert_eq!(merged(&["CB", "CY", "CR"], &["CB"], cr), "CB");
        assert_eq!(merged(&["CR", "CB", "CY"], &["CR"], &[("CB", None)]), "CR");
        assert_eq!(merged(&["CB", "CR", "CY"], &["CB", "CY"], cr), "CB CY");
        assert_eq!(merged(&["CB", "CR", "CY"], &[], cr), "CB CY");
        assert_eq!(merged(&["CR", "CB", "CY"], &[], &[("CB", None)]), "CR CY");
    }
}
