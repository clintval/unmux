//! The demux engine: parse records, extract tags/UMIs (carrying qualities), and
//! fan out to per-sample FASTX/SAM/BAM/CRAM (SAM/BAM/CRAM written unmapped),
//! splitcode-style, on the `sassy` matcher.
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
//! provenance + per-sample `@RG`, the metrics TSVs, CRAM output, and the rayon
//! record-parallel + pooled-BGZF engine are all implemented; matching is
//! parallelized while the writer stays on the consumer thread.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use std::fs::File;
use std::io::BufWriter;

use anyhow::{anyhow, bail, Context, Result};
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::Header;
use pooled_writer::{bgzf::BgzfCompressor, Pool, PoolBuilder, PooledWriter};
use rayon::prelude::*;
use smallvec::SmallVec;
use thousands::Separable;

use crate::extract::{
    extract, span_range, Corrected, Extracted, GroupHit, GroupHits, MatchSpan, Segment,
};
use crate::fanout::{
    compile_routing, expand_from_group, load_sample_sheet, Disposition, RemoveTarget, Routing,
    Target,
};
use crate::grammar::{
    Anchor, DemuxPlan, Endpoint, GroupSource, MatchMode, OutputPattern, Placeholder, SampleSpec,
    SpanBody, Template,
};
use crate::input::{Fragment, FragmentReader, InputRecord, SniffedFormat};
use crate::matcher::{
    match_group, match_group_at, window_after, CompiledGroup, GroupOutcome, MatchStrand, Prereq,
    RelativeWindow, Scratch,
};
use crate::metrics::Metrics;
use crate::output::{
    build_tag_data, default_read_group_header, insert_string_tag, qual_tag_name, read_group_header,
    resolve_metrics_path, resolve_pattern, MultiWriter, PathContext,
};
use crate::qc::{removed_slug, routed_slug, unassigned_slug, Unassigned};
use crate::tags::{load_tag_file, TagSet};
use crate::writer::{encode_fragment, output_format, OutputFormat, OutputRead};

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
    /// `--qc-tag`: emit a per-record JSON demux-provenance slug under this SAM
    /// tag (a local-use tag, default `ZS`). `None` = off (the default).
    pub qc_tag: Option<String>,
    /// `--per-record`: disable auto pair-detection (interleaving is
    /// auto-detected otherwise).
    pub per_record: bool,
    /// `--compression`: level (0-9) for BGZF (BAM) and gzip (FASTX.gz) outputs;
    /// CRAM uses its own per-block codecs. Output format is inferred from the
    /// file extension.
    pub compression: u8,
    /// `--threads`: worker thread count.
    pub threads: usize,
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
/// file is read and an inline set built only once.
fn load_tag_sets(plan: &DemuxPlan) -> Result<HashMap<String, TagSet>> {
    let mut tag_sets = HashMap::with_capacity(plan.groups.len());
    for group in &plan.groups {
        let tag_set = match &group.source {
            GroupSource::File(path) => load_tag_file(path)
                .with_context(|| format!("failed to load tag set for group `{}`", group.name))?,
            GroupSource::Inline(tags) => TagSet::from_inline(tags),
        };
        tag_sets.insert(group.name.clone(), tag_set);
    }
    Ok(tag_sets)
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

/// Compile each `--group` into a matcher [`CompiledGroup`] over its loaded tag
/// set, return them in dependency order (each group after the upstream group
/// its `next`/`prev` link depends on), and the `name -> canonical index` map
/// for that final order. The sequential prerequisites are resolved to canonical
/// indices once the final order (and thus the indices) is known.
fn compile_groups(
    plan: &DemuxPlan,
    tag_sets: &HashMap<String, TagSet>,
) -> Result<(Vec<CompiledGroup>, HashMap<String, usize>)> {
    let mut compiled = Vec::with_capacity(plan.groups.len());
    let mut name_prereqs = Vec::with_capacity(plan.groups.len());
    for group in &plan.groups {
        let attrs = &group.attrs;
        let dist = attrs.dist.unwrap_or_default();
        let mode = attrs.mode.unwrap_or(MatchMode::All);
        // Guard the IUPAC tag set up front: reject a combinatorial explosion or
        // an indistinguishable tag pair, and warn on a within-`dist` overlap
        // that `mode=nearest` would only disambiguate.
        for warning in
            crate::iupac::validate_group(&group.name, &tag_sets[&group.name].seqs(), &dist, mode)?
        {
            log::warn!("{warning}");
        }
        let mut compiled_group = CompiledGroup {
            name: group.name.clone(),
            tags: tag_sets[&group.name].seqs(),
            dist,
            loc: attrs.loc.clone(),
            revcomp: attrs.revcomp.unwrap_or(false),
            mode,
            delta: attrs.delta.unwrap_or(0),
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
        };
        // Precompute the sassy v2 batched search batches (and their covered tag
        // set) for large equal-length tag buckets.
        compiled_group.encode();
        if let Some(warning) = anchor5p_window_overrun_warning(&compiled_group) {
            log::warn!("{warning}");
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
/// given, else the longest common prefix of the input file names (trimmed of
/// trailing separators), else `pool`.
fn pool_id(plan: &DemuxPlan) -> String {
    if let Some(pool) = &plan.pool {
        return pool.clone();
    }
    let names: Vec<String> = plan
        .inputs
        .iter()
        .filter_map(|path| path.file_name().map(|n| n.to_string_lossy().into_owned()))
        .collect();
    let prefix = longest_common_prefix(&names);
    let trimmed = prefix.trim_end_matches(['.', '_', '-', ' ']);
    if trimmed.is_empty() {
        "pool".to_string()
    } else {
        trimmed.to_string()
    }
}

/// The longest common prefix (by bytes) of a set of strings.
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
    first[..end].to_string()
}

/// Build the routing table from the plan's fan-out source (`--sample` inline, a
/// `--sample-sheet`, or `--sample-from-group`) plus the `--remove` rules,
/// sharing the already-loaded tag sets.
fn build_routing(
    plan: &DemuxPlan,
    tag_sets: &HashMap<String, TagSet>,
    group_index: &HashMap<String, usize>,
    pool: &str,
) -> Result<Routing> {
    let (samples, demux) = match &plan.samples {
        SampleSpec::None => (Vec::new(), false),
        SampleSpec::Inline(samples) => (samples.clone(), true),
        SampleSpec::Sheet(path) => {
            let declared: HashSet<&str> = plan.groups.iter().map(|g| g.name.as_str()).collect();
            (load_sample_sheet(path, &declared)?, true)
        }
        SampleSpec::FromGroup(group) => {
            let tag_set = tag_sets
                .get(group)
                .with_context(|| format!("no loaded tag set for --sample-from-group `{group}`"))?;
            (expand_from_group(group, tag_set), true)
        }
    };
    compile_routing(
        &samples,
        &plan.remove,
        tag_sets,
        group_index,
        pool,
        demux,
        plan.require_samples_explain_all_tags,
    )
}

/// A SAM/BAM `--out` destination (a path, or `None` for stdout), the fan-out
/// targets whose records land there, and its format.
type OutDestination = (Option<PathBuf>, Vec<Target>, OutputFormat);

/// The SAM/BAM/CRAM `--out` destinations and the fan-out targets whose records
/// land in each (a single-file `--out` collects every target; a `%sample`
/// pattern splits them), in first-seen order. FASTX destinations have no
/// header, so they are excluded. Shared by the `@RG` header precompute and the
/// pooling exclusion (the SAM/CRAM `--out` files among these stay on the inline
/// single-EOF path; BAM `--out` is pooled).
fn out_destinations(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
) -> Result<Vec<OutDestination>> {
    let mut by_dest: Vec<OutDestination> = Vec::new();
    let mut index: HashMap<Option<PathBuf>, usize> = HashMap::new();
    if let Some(pattern) = &plan.out {
        for target in routing.targets() {
            let ctx = PathContext {
                pool,
                sample: Some(&target.sample),
                sub_sample: target.sub_sample.as_deref(),
                ordinal: None,
                source: None,
            };
            let path = resolve_pattern(pattern, &ctx);
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
    }
    Ok(by_dest)
}

/// Every directed output FILE path the engine could write, computed from the
/// plan and routing before streaming so each can be created up front: a
/// directed file always exists, even when it receives no records (an empty but
/// valid file). stdout (a `None` path) is excluded; it cannot be pre-created
/// and is a single sink. This mirrors the path resolution in `write_body` /
/// `write_raw_bin`; the lazy fallback in `MultiWriter::write` covers any path
/// this misses, so a mismatch never drops a record.
fn enumerate_dests(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
    fragment_width: usize,
) -> Result<Vec<PathBuf>> {
    let mut paths: Vec<PathBuf> = Vec::new();
    let mut seen: HashSet<PathBuf> = HashSet::new();
    let add = |path: Option<PathBuf>, seen: &mut HashSet<PathBuf>, paths: &mut Vec<PathBuf>| {
        if let Some(path) = path {
            if seen.insert(path.clone()) {
                paths.push(path);
            }
        }
    };

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
                add(Some(path), &mut seen, &mut paths);
            }
        }
        for target in targets {
            for path in out_paths(pattern, pool, Some(target), body_count, input_formats)? {
                add(Some(path), &mut seen, &mut paths);
            }
        }
    }
    // The `--unassigned` and `--remove` bins (raw segments; `%source` fans out
    // over the input files).
    if let Some(pattern) = &plan.unassigned {
        for path in bin_paths(pattern, pool, plan.inputs.len()) {
            add(Some(path), &mut seen, &mut paths);
        }
    }
    for rule in &plan.remove {
        if let Some(pattern) = &rule.pattern {
            for path in bin_paths(pattern, pool, plan.inputs.len()) {
                add(Some(path), &mut seen, &mut paths);
            }
        }
    }
    Ok(paths)
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

/// Precompute the `@RG`-bearing SAM/BAM header for each distinct `--out`
/// destination, so a per-sample file lists its own read groups and a
/// single-file `--out` lists every target's. FASTX destinations and the raw
/// bins are absent from the map and get the default `@PG`-only header.
fn out_headers(
    plan: &DemuxPlan,
    routing: &Routing,
    pool: &str,
    input_formats: &[SniffedFormat],
    command_line: Option<&str>,
) -> Result<HashMap<Option<PathBuf>, noodles::sam::Header>> {
    let mut headers = HashMap::new();
    for (dest, targets, _) in out_destinations(plan, routing, pool, input_formats)? {
        let refs: Vec<&Target> = targets.iter().collect();
        headers.insert(dest, read_group_header(&refs, &plan.rg_tags, command_line)?);
    }
    Ok(headers)
}

/// The most BGZF compressor threads the writer pool uses, regardless of
/// `--threads`. Profiling the 2M dual-index `.fq.gz` demux showed compression
/// is ~half the CPU (not the minor stage originally assumed), and the warm
/// match cache leaves matcher threads idle, so the pool is allowed a meaningful
/// share; the cap bounds it at high `--threads` where compression parallelism
/// has diminishing returns and BGZF pool overhead grows.
pub const MAX_WRITER_THREADS: usize = 6;

/// How one `--threads N` is split across the engine's roles. Our input is a
/// single lockstep cursor, so one read-ahead thread suffices; the main/consumer
/// thread owns the writer and applies outcomes in input order; `writer` threads
/// run the BGZF compressor pool (zero when there is no pooled output); every
/// remaining thread joins the rayon matcher pool. For compressed output the
/// spare threads (after the reader + consumer) are split evenly between the
/// BGZF pool and the matchers (compression is ~half the CPU and the warm match
/// cache idles matchers); with no pooled output every spare thread is a
/// matcher.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ThreadBudget {
    /// The read-ahead/decompress thread.
    pub reader: usize,
    /// The orchestrator/consumer thread (owns the writer; applies outcomes in
    /// order).
    pub main: usize,
    /// The BGZF compressor pool size; zero when nothing is pooled.
    pub writer: usize,
    /// The rayon matcher pool size (the parallel CPU stage).
    pub matchers: usize,
}

/// The BGZF compressor pool size for `threads` total threads, given whether the
/// run has any poolable (compressed) output. Zero with nothing to compress, and
/// zero when there are too few threads to staff a pool: a pool needs reader +
/// main + at least one writer + at least one matcher, i.e. `threads >= 4`;
/// below that, compression runs inline on the consumer (`SinkWriter::Gzip` /
/// `Sink::BamInline`). At `threads >= 4` the spare threads (after reader +
/// main) split evenly between the BGZF pool and the matchers (compression is
/// ~half the CPU and the warm match cache idles matchers), capped at
/// [`MAX_WRITER_THREADS`]; the even split keeps `matchers >= writer`.
fn writer_pool_size(threads: usize, has_poolable: bool) -> usize {
    if !has_poolable || threads < 4 {
        return 0;
    }
    (threads.saturating_sub(2) / 2).clamp(1, MAX_WRITER_THREADS)
}

/// Partition `--threads N` into the engine's [`ThreadBudget`] after reserving
/// `writer_threads` for the BGZF pool. Assumes the caller has checked `threads`
/// leaves at least one matcher; the `max(1)` is a defensive floor.
pub fn partition_threads(threads: usize, writer_threads: usize) -> ThreadBudget {
    let reader = 1;
    let main = 1;
    let matchers = threads
        .saturating_sub(reader + main + writer_threads)
        .max(1);
    ThreadBudget {
        reader,
        main,
        writer: writer_threads,
        matchers,
    }
}

/// The static per-run engine context shared across fragments: the parsed plan,
/// the compiled groups (in dependency order), the canonical `name -> index` map
/// (for resolving `@grp` anchor names to positional slots), the routing table,
/// the resolved pool id, and the input formats.
struct Engine<'a> {
    plan: &'a DemuxPlan,
    groups: &'a [CompiledGroup],
    group_index: &'a HashMap<String, usize>,
    routing: &'a Routing,
    pool: &'a str,
    input_formats: &'a [SniffedFormat],
    /// Extract-stream names referenced exactly once across all `--template`s
    /// and `--tag`s. A single-stream template whose stream is movable is the
    /// sole consumer of that extract, so the assembled body can MOVE the
    /// extracted bases/quals rather than copy them (the common case for
    /// whole-record template records, which are large; avoids a per-record copy
    /// + alloc on the consumer).
    movable_streams: &'a HashSet<&'a str>,
    /// Extract-stream names drawn on by any `--template`/`--tag`, the
    /// denominator basis for the `frac_bases_unextracted` metric.
    /// Plan-invariant, so it is built once here rather than rebuilt for every
    /// retained record in [`frac_bases`].
    routed_streams: &'a HashSet<&'a str>,
    /// Per-destination SAM/BAM headers (`@RG` lines), keyed by output path
    /// (`None` = stdout), shared read-only so a worker can encode SAM/BAM
    /// records (encode-on-workers) with the exact header the consumer's writer
    /// was opened with. FASTX ignores it.
    headers: &'a HashMap<Option<PathBuf>, Header>,
    /// The fallback header for a destination absent from `headers` (the raw
    /// `--unassigned`/`--remove` bins, or a non-demux run): the minimal `@HD` +
    /// `@PG`-provenance header, matching `MultiWriter::open`'s fallback so
    /// worker-encoded bytes stay byte-identical.
    default_header: &'a Header,
}

/// Run the demux engine: match tag groups per record, route each record with
/// the compiled [`Routing`], fan it out to per-sample `--out` (carrying `--tag`
/// fields and an `@RG` read group) or to the `--remove` / `--unassigned` bins,
/// tally metrics, and emit the metrics TSVs / stderr log.
fn run_engine(plan: &DemuxPlan, command_line: Option<&str>) -> Result<()> {
    let tag_sets = load_tag_sets(plan)?;
    let (groups, group_index) = compile_groups(plan, &tag_sets)?;
    let pool = pool_id(plan);
    let routing = build_routing(plan, &tag_sets, &group_index, &pool)?;

    let reader = FragmentReader::open(&plan.inputs, plan.per_record)?;
    let input_formats = reader.formats().to_vec();
    let headers = out_headers(plan, &routing, &pool, &input_formats, command_line)?;
    // The fallback header for any destination absent from `headers`
    // (pass-through `--out`, the raw `--unassigned`/`--remove` bins): `@PG`
    // provenance plus a default `@RG` whose ID/SM/LB are the pool id, so every
    // SAM/BAM/CRAM record carries a read group even with no `--sample`. A
    // worker uses it for encode-on-workers exactly where `MultiWriter::open`
    // would, so the bytes stay identical.
    let default_header = default_read_group_header(&pool, &plan.rg_tags, command_line)?;
    // Declared before `writer` so that, on any unwind, `writer` (and its
    // PooledWriters) drops before the pool: a PooledWriter dropped after the
    // pool is stopped would panic. The explicit shutdown below also enforces
    // this order; the declaration order is the belt-and-suspenders guarantee.
    let mut writer_pool: Option<Pool> = None;
    let mut writer = MultiWriter::new(
        headers.clone(),
        default_header.clone(),
        input_formats.clone(),
        plan.compression,
    );
    // Every directed output file, created up front so a directed file always
    // exists, even empty.
    let dests = enumerate_dests(
        plan,
        &routing,
        &pool,
        &input_formats,
        reader.fragment_width(),
    )?;
    log_run_shape(plan, &input_formats, &dests);

    // Route poolable destinations (gzipped FASTX and BAM) through a shared BGZF
    // compressor pool, so compression runs off the consumer thread (the
    // single-file BAM `--out` of a combinatorial run is otherwise
    // serial-compression-bound). The pooled BAM sink emits the BGZF EOF on
    // close (writer.rs `Sink::BamPooled` finish), so a user-facing BAM `--out`
    // stays valid (asserted by the BAM round-trip tests + `samtools
    // quickcheck`). Only CRAM/SAM `--out` stay inline: CRAM has its own
    // per-block codec (not BGZF) and SAM is uncompressed text, so neither is
    // poolable (and `OutputFormat::is_pooled` already excludes them; this keeps
    // them inline regardless).
    let out_alignment: HashSet<PathBuf> = out_destinations(plan, &routing, &pool, &input_formats)?
        .into_iter()
        .filter(|(_, _, format)| !matches!(format, OutputFormat::Bam))
        .filter_map(|(dest, _, _)| dest)
        .collect();
    let poolable: Vec<PathBuf> = dests
        .iter()
        .filter(|path| !out_alignment.contains(*path))
        .filter(|path| {
            output_format(Some(path), &input_formats)
                .map(OutputFormat::is_pooled)
                .unwrap_or(false)
        })
        .cloned()
        .collect();
    // Defensive max(1): clap already guarantees --threads >= 1, but guard for
    // library callers.
    let threads = plan.threads.max(1);
    let writer_threads = writer_pool_size(threads, !poolable.is_empty());
    let budget = partition_threads(threads, writer_threads);
    match threads {
        1 => log::info!("threads=1: serial (read + match + write inline, no pools)"),
        2 => log::info!("threads=2: reader=1, match + write inline (no rayon, no BGZF pool)"),
        _ => log::info!(
            "threads={}: reader={}, main={}, matchers={}, writer={} (rayon + BGZF pool)",
            threads,
            budget.reader,
            budget.main,
            budget.matchers,
            budget.writer,
        ),
    }

    if writer_threads > 0 {
        let mut builder = PoolBuilder::<BufWriter<File>, BgzfCompressor>::new()
            .threads(budget.writer)
            .compression_level(plan.compression)
            .map_err(|e| anyhow!("invalid --compression {} for BGZF: {e}", plan.compression))?;
        let mut exchanged: HashMap<PathBuf, PooledWriter> = HashMap::new();
        for path in &poolable {
            let file = crate::writer::create_output_file(path)?;
            exchanged.insert(path.clone(), builder.exchange(BufWriter::new(file)));
        }
        writer_pool = Some(
            builder
                .build()
                .map_err(|e| anyhow!("failed to start the BGZF writer pool: {e}"))?,
        );
        writer.set_pooled(exchanged);
    }

    // Streams referenced exactly once across templates + tags can be MOVED into
    // their (sole) output body instead of copied (see Engine::movable_streams).
    let mut stream_refs: HashMap<&str, usize> = HashMap::new();
    for streams in plan
        .templates
        .iter()
        .map(|t| &t.streams)
        .chain(plan.tags.iter().map(|t| &t.streams))
    {
        for stream in streams {
            *stream_refs.entry(stream.name.as_str()).or_insert(0) += 1;
        }
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
        routing: &routing,
        pool: &pool,
        input_formats: &input_formats,
        movable_streams: &movable_streams,
        routed_streams: &routed_streams,
        headers: &headers,
        default_header: &default_header,
    };
    let mut metrics = Metrics::new(pool.clone(), routing.targets());

    // Create the output files and run the parallel engine, capturing the result
    // WITHOUT `?` so the mandatory shutdown always runs in order even on error:
    // `finish` closes every PooledWriter (flushing its BGZF EOF), and only then
    // is the pool stopped (joining the compressor threads). Reversing that
    // order would hang, and a PooledWriter dropped after `stop_pool` would
    // panic. `MultiWriter::finish` finalizes every writer even when one errors,
    // so the pool is always safe to stop once it returns; we stop it
    // unconditionally (its threads never leak) and surface the most actionable
    // error first: the run error, then a finalize error, then a pool-stop
    // error.
    let run_result = (|| -> Result<()> {
        writer.create_all(&dests)?;
        match threads {
            1 => drive_serial(&engine, reader, &mut writer, &mut metrics),
            2 => drive_pipelined(&engine, reader, &mut writer, &mut metrics),
            _ => drive_parallel(&engine, budget, reader, &mut writer, &mut metrics),
        }
    })();
    let finished = writer.finish();
    let stopped = match writer_pool {
        Some(mut writer_pool) => writer_pool
            .stop_pool()
            .map_err(|e| anyhow!("failed to stop the BGZF writer pool: {e}")),
        None => Ok(()),
    };
    run_result?;
    finished?;
    stopped?;

    // Provenance lives in the header written up front: `@PG`
    // (program/version/command line) and the per-sample `@RG`. The demux stats
    // go to stderr and the `--metrics-*` TSVs below (they are not embedded back
    // into the output headers, which would force a full post-pass rewrite of
    // each file).
    metrics.log();
    if let Some(pattern) = &plan.metrics_per_sample {
        metrics.write_per_sample_tsv(&resolve_metrics_path(pattern, &pool))?;
    }
    if let Some(pattern) = &plan.metrics_summary {
        metrics.write_summary_tsv(&resolve_metrics_path(pattern, &pool))?;
    }
    Ok(())
}

/// The result of matching every group against one record: either a
/// find-constraint failure (the record is unassigned) or the per-group hits,
/// indexed positionally by canonical group index, that drive routing and body
/// assembly. Owned and `Send` (one `Vec<Option<GroupHit>>`, no per-record heap
/// maps or `String` clones), so it crosses from a rayon matcher worker to the
/// serial consumer.
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
/// built on a rayon worker and handed to the serial writer. `tags` is `None`
/// when the record carries no SAM data (an empty tag set emits no field).
struct PreparedRead {
    name: Vec<u8>,
    bases: Vec<u8>,
    quals: Option<Vec<u8>>,
    tags: Option<Data>,
    read_group: Option<String>,
}

impl PreparedRead {
    /// Borrow this prepared record as the writer's [`OutputRead`], just before
    /// the (serial) write.
    fn as_output(&self) -> OutputRead<'_> {
        OutputRead {
            name: &self.name,
            bases: &self.bases,
            quals: self.quals.as_deref(),
            tags: self.tags.as_ref(),
            read_group: self.read_group.as_deref(),
        }
    }
}

/// What the consumer writes for one [`PreparedWrite`]: bytes already encoded
/// for the destination's format on the worker (encode-on-workers; appended
/// verbatim, the common path), or the structured records for a format that
/// cannot be pre-encoded to an independent per-fragment byte slice (CRAM).
enum WritePayload {
    /// Pre-encoded record bytes (from [`encode_fragment`]) for the alignment
    /// formats (SAM/BAM), whose heavier encode is worth moving off the
    /// consumer.
    Encoded(Vec<u8>),
    /// Structured records, encoded on the consumer: CRAM (not pre-encodable per
    /// fragment) and the FASTX formats (light text encode that does not pay for
    /// the worker handoff).
    Records(Vec<PreparedRead>),
}

/// One write the consumer must issue: a destination (`None` is stdout) and its
/// payload, in the exact order the inline `write_body` / `write_raw_bin`
/// emitted them.
struct PreparedWrite {
    path: Option<PathBuf>,
    payload: WritePayload,
}

/// The metrics tally for one fragment, owned (`Target` is `Clone`) so it
/// crosses the worker boundary. Commutative integer sums, applied on the
/// consumer in input order (order is irrelevant to the sums).
enum Tally {
    Unassigned,
    Removed,
    PassThrough {
        denom: u64,
        unext: u64,
    },
    Assigned {
        target: Target,
        denom: u64,
        unext: u64,
    },
}

/// The fully-prepared outcome of one fragment: the per-record base lengths (for
/// the record-length metric), the metric tally, and the writes to issue in
/// order. Owned + `Send`, so the whole route/extract/assemble/encode pipeline
/// can run on a rayon worker, leaving the consumer to only tally and issue the
/// (inherently serial, input-ordered) writes.
struct Prepared {
    lengths: SmallVec<[usize; 4]>,
    tally: Tally,
    writes: Vec<PreparedWrite>,
}

/// The result of assembling a routed record's body: the writes to issue, or the
/// name of a `--template` stream whose anchoring group did not match (so the
/// record is rerouted to the `--unassigned` bin).
enum BodyOutcome {
    Assembled(Vec<PreparedWrite>),
    MissingStream(String),
}

/// Match every group against one record in dependency order: the pure,
/// parallelizable half of per-record work. No writer, no metrics; it borrows
/// the record and the shared engine and returns an owned [`FragmentMatch`], so
/// it runs on a rayon worker with a thread-local [`Scratch`] while the writer
/// and metrics stay on the serial consumer.
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
        let outcome = match &group.match_streams {
            // A `match=` group matches the joined `--extract` streams, not a
            // read window.
            Some(streams) => match_joined(
                group,
                scratch,
                streams,
                &engine.plan.extracts,
                &segments,
                &hits,
                engine.group_index,
            )?,
            None => match_group_for(group, scratch, &base_segments, &hits),
        };
        // A failed find constraint (e.g. an unmet minFindsPerGroup) leaves the
        // read unassigned; carry the specific failure and the partial hits so
        // far for the optional QC slug.
        if let Some(failure) = group.find_failure(&outcome) {
            return Ok(FragmentMatch::Unassigned {
                reason: Unassigned::FindConstraint {
                    group: idx,
                    failure,
                },
                hits,
            });
        }
        if let Some(best) = &outcome.best {
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
/// input-ordered writes, so the whole pipeline runs on a rayon worker with a
/// thread-local [`Scratch`] while only [`apply_prepared`] stays on the
/// consumer.
fn process_fragment(
    engine: &Engine,
    scratch: &mut Scratch,
    fragment: &Fragment,
) -> Result<Prepared> {
    let fmatch = match_fragment(engine, scratch, fragment)?;
    let mut prepared = prepare(engine, fragment, fmatch)?;
    encode_writes(engine, &mut prepared.writes)?;
    Ok(prepared)
}

/// Encode each write's records into the destination's on-disk bytes here on the
/// worker (encode-on-workers), so the serial consumer only appends them in
/// input order rather than running the per-record SAM/BAM/FASTX encode itself.
/// CRAM stays structured (its records are container/slice-coded, not
/// independent per-fragment byte slices) and is encoded by the consumer. The
/// format and header are resolved exactly as the consumer's `MultiWriter::open`
/// does (per path, falling back to the provenance default), so the appended
/// bytes are byte-identical.
fn encode_writes(engine: &Engine, writes: &mut [PreparedWrite]) -> Result<()> {
    for write in writes {
        let format = output_format(write.path.as_deref(), engine.input_formats)?;
        // Encode only the alignment formats here. BAM/SAM go through the
        // heavier build_record + serialize path, so moving it to the worker is
        // a measured ~1.8x win on BAM `--out`; FASTQ/FASTA text formatting is
        // light enough that the per-fragment byte handoff makes
        // encode-on-workers a slight net loss (measured), and CRAM is
        // container/slice-coded (not pre-encodable per fragment). Those stay on
        // the serial consumer path.
        if !matches!(format, OutputFormat::Bam | OutputFormat::Sam) {
            continue;
        }
        let bytes = {
            let WritePayload::Records(reads) = &write.payload else {
                continue;
            };
            let header = engine
                .headers
                .get(&write.path)
                .unwrap_or(engine.default_header);
            let outputs: SmallVec<[OutputRead; 2]> =
                reads.iter().map(PreparedRead::as_output).collect();
            encode_fragment(format, header, &outputs)?
        };
        write.payload = WritePayload::Encoded(bytes);
    }
    Ok(())
}

/// Prepare one matched record into an owned, `Send` [`Prepared`] (the
/// per-record lengths, the metric tally, and the writes to issue): route it,
/// compute `frac_bases`, and assemble its body / bins, but touch neither the
/// writer nor the metrics. This is the pure half of the old `apply_match`, so
/// it can run on a rayon worker (alongside matching) while only the tally + the
/// ordered writes stay serial. The assigned/pass-through decision reflects the
/// post-assembly outcome (a missing `--template` stream reroutes to
/// unassigned), exactly as the inline path did.
fn prepare(engine: &Engine, fragment: &Fragment, fmatch: FragmentMatch) -> Result<Prepared> {
    let lengths: SmallVec<[usize; 4]> = fragment.records.iter().map(|r| r.bases.len()).collect();
    let hits = match fmatch {
        FragmentMatch::Unassigned { reason, hits } => {
            let writes = prepare_unassigned(engine, fragment, &reason, &hits);
            return Ok(Prepared {
                lengths,
                tally: Tally::Unassigned,
                writes,
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

    let (tally, writes) = match engine.routing.route(&hits) {
        Disposition::Removed(remove) => (
            Tally::Removed,
            prepare_removed(engine, fragment, remove, &hits),
        ),
        Disposition::Unassigned => (
            Tally::Unassigned,
            prepare_unassigned(engine, fragment, &Unassigned::NoSample, &hits),
        ),
        Disposition::PassThrough => {
            let (denom, unext) = frac_bases(engine, &segments, &hits)?;
            match prepare_body(engine, None, &hits, &segments, fragment)? {
                BodyOutcome::Assembled(writes) => (Tally::PassThrough { denom, unext }, writes),
                BodyOutcome::MissingStream(missing) => (
                    Tally::Unassigned,
                    prepare_unassigned(
                        engine,
                        fragment,
                        &Unassigned::MissingStream(missing),
                        &hits,
                    ),
                ),
            }
        }
        Disposition::Assigned(target) => {
            let (denom, unext) = frac_bases(engine, &segments, &hits)?;
            match prepare_body(engine, Some(target), &hits, &segments, fragment)? {
                BodyOutcome::Assembled(writes) => (
                    Tally::Assigned {
                        target: target.clone(),
                        denom,
                        unext,
                    },
                    writes,
                ),
                BodyOutcome::MissingStream(missing) => (
                    Tally::Unassigned,
                    prepare_unassigned(
                        engine,
                        fragment,
                        &Unassigned::MissingStream(missing),
                        &hits,
                    ),
                ),
            }
        }
    };
    Ok(Prepared {
        lengths,
        tally,
        writes,
    })
}

/// Apply one [`Prepared`] on the serial consumer: tally it and issue its writes
/// in order. The writer and metrics are single-threaded; everything else now
/// happens on the worker in [`prepare`].
fn apply_prepared(
    prepared: Prepared,
    writer: &mut MultiWriter,
    metrics: &mut Metrics,
) -> Result<()> {
    metrics.record_processed();
    metrics.record_read_lengths(prepared.lengths.iter().copied());
    match &prepared.tally {
        Tally::Unassigned => metrics.record_unassigned(),
        Tally::Removed => metrics.record_removed(),
        Tally::PassThrough { denom, unext } => metrics.record_pass_through(*denom, *unext),
        Tally::Assigned {
            target,
            denom,
            unext,
        } => metrics.record_assigned(target, *denom, *unext),
    }
    for write in &prepared.writes {
        match &write.payload {
            // Encoded on the worker: the consumer only appends the bytes (the
            // merge-point lever).
            WritePayload::Encoded(bytes) => writer.write_encoded(write.path.as_deref(), bytes)?,
            // CRAM fallback: the consumer encodes the structured reads.
            WritePayload::Records(reads) => {
                let outputs: SmallVec<[OutputRead; 2]> =
                    reads.iter().map(PreparedRead::as_output).collect();
                writer.write(write.path.as_deref(), &outputs)?;
            }
        }
    }
    Ok(())
}

/// Fragments per chunk handed from the read-ahead thread to the matcher pool.
/// Large enough that the per-worker [`Scratch`] (match cache + searcher),
/// rebuilt per chunk by `map_init`, amortizes its warm-up across many records
/// (see [`MATCHER_MIN_JOB`]); bounded by [`READ_AHEAD_CHUNKS`] so memory stays
/// capped.
const READ_CHUNK: usize = 32768;
/// Minimum records a single rayon job processes with one `map_init`
/// [`Scratch`]. Coarsening the split (vs rayon's default fine granularity) lets
/// the per-job match cache warm up and the searcher be reused across thousands
/// of records rather than ~one, which is the whole point of the cache.
const MATCHER_MIN_JOB: usize = 2048;
/// Bounded read-ahead: at most this many chunks queued ahead of the consumer,
/// so a slow consumer applies back-pressure to the reader instead of growing
/// memory.
const READ_AHEAD_CHUNKS: usize = 2;

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

    /// Count one record, logging the running total at each [`PROGRESS_UNIT`]
    /// boundary.
    fn record(&mut self) {
        self.seen += 1;
        if self.seen.is_multiple_of(PROGRESS_UNIT) {
            log::info!("Demultiplexed {} records", self.seen.separate_with_commas());
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

/// Pull up to `chunk_size` fragments from the reader into one owned chunk, or
/// `Ok(None)` at end of input. Shared by the read-ahead thread (`read_chunks`)
/// and the serial driver (`drive_serial`), which both need the same bounded
/// pull. A read error propagates.
fn next_chunk(reader: &mut FragmentReader, chunk_size: usize) -> Result<Option<Vec<Fragment>>> {
    let mut chunk = Vec::with_capacity(chunk_size);
    for _ in 0..chunk_size {
        match reader.next_fragment()? {
            Some(fragment) => chunk.push(fragment),
            None => break,
        }
    }
    if chunk.is_empty() {
        Ok(None)
    } else {
        Ok(Some(chunk))
    }
}

/// Match every fragment in a chunk serially with one warm `Scratch`, preserving
/// positional order. The serial counterpart to the rayon `par_iter` in
/// `drive_parallel`; used by the 1- and 2-thread drivers, where the consumer
/// does the matching itself.
fn match_chunk_serial(
    engine: &Engine,
    scratch: &mut Scratch,
    chunk: &[Fragment],
) -> Vec<Result<Prepared>> {
    chunk
        .iter()
        .map(|fragment| process_fragment(engine, scratch, fragment))
        .collect()
}

/// Apply one chunk's prepared outcomes to the writer in input order, tallying
/// metrics and counting progress, stopping at (and returning) the first error.
/// Shared by all three drivers; only the serial, inherently-ordered writes
/// happen here.
fn apply_chunk(
    prepared: Vec<Result<Prepared>>,
    writer: &mut MultiWriter,
    metrics: &mut Metrics,
    progress: &mut Progress,
) -> Result<()> {
    for outcome in prepared {
        apply_prepared(outcome?, writer, metrics)?;
        progress.record();
    }
    Ok(())
}

/// The read-ahead thread: pull fragments into bounded chunks and hand each to
/// the matcher pool over `tx`. Returns when the input is exhausted or the
/// consumer hangs up; a read error propagates as `Err`, surfaced when the
/// consumer joins this thread. Dropping `tx` on return closes the channel,
/// signalling end-of-input to the consumer.
fn read_chunks(
    mut reader: FragmentReader,
    tx: std::sync::mpsc::SyncSender<Vec<Fragment>>,
    chunk_size: usize,
) -> Result<()> {
    while let Some(chunk) = next_chunk(&mut reader, chunk_size)? {
        if tx.send(chunk).is_err() {
            return Ok(()); // the consumer hung up (it bailed); stop reading.
        }
    }
    Ok(())
}

/// Drive the engine on a single thread (`--threads 1`): no read-ahead thread,
/// no rayon pool, no BGZF pool. Read a chunk, match it serially with one warm
/// `Scratch`, then apply it, all on this thread; compression, if any, runs
/// inline on the sink writers. Uses exactly one core, so it never
/// oversubscribes a single-core runner. Output is input-ordered, like the other
/// drivers.
fn drive_serial(
    engine: &Engine,
    mut reader: FragmentReader,
    writer: &mut MultiWriter,
    metrics: &mut Metrics,
) -> Result<()> {
    let mut progress = Progress::new();
    let mut scratch = Scratch::new();
    while let Some(chunk) = next_chunk(&mut reader, READ_CHUNK)? {
        let prepared = match_chunk_serial(engine, &mut scratch, &chunk);
        apply_chunk(prepared, writer, metrics, &mut progress)?;
    }
    Ok(())
}

/// Drive the engine on two threads (`--threads 2`): a read-ahead thread streams
/// chunks while this (consumer) thread matches each chunk serially with one
/// warm `Scratch` and applies it. No rayon matcher pool and no BGZF pool
/// (compression is inline), so the run uses exactly two cores: one
/// reading/decompressing, one matching + writing. A consumer error takes
/// precedence over a read error or a reader panic, matching `drive_parallel`.
fn drive_pipelined(
    engine: &Engine,
    reader: FragmentReader,
    writer: &mut MultiWriter,
    metrics: &mut Metrics,
) -> Result<()> {
    let mut progress = Progress::new();
    let mut scratch = Scratch::new();
    let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<Fragment>>(READ_AHEAD_CHUNKS);
    let reader_handle = std::thread::spawn(move || read_chunks(reader, tx, READ_CHUNK));

    let mut outcome: Result<()> = Ok(());
    for chunk in &rx {
        let prepared = match_chunk_serial(engine, &mut scratch, &chunk);
        if let Err(error) = apply_chunk(prepared, writer, metrics, &mut progress) {
            outcome = Err(error);
            break;
        }
    }
    drop(rx);
    let reader_result = reader_handle.join();
    outcome?;
    match reader_result {
        Ok(read_outcome) => read_outcome,
        Err(_) => Err(anyhow!("the read-ahead thread panicked")),
    }
}

/// Drive the engine in parallel: a read-ahead thread streams bounded chunks of
/// fragments to a rayon worker pool, and the consumer applies each chunk's
/// outcomes to the single writer in input order. The full per-record pipeline
/// ([`process_fragment`]: match + route + extract + assemble + encode) runs on
/// the pool with a thread-local [`Scratch`] per worker, returning owned
/// [`Prepared`] outcomes; only the metric tally and the inherently-serial
/// writes stay on the consumer (main) thread, while the now-`Send` reader is
/// moved to the read-ahead thread. Output stays in input order (the serial
/// sinks and round-trip tests observe it, and it makes runs deterministic
/// regardless of `--threads`).
fn drive_parallel(
    engine: &Engine,
    budget: ThreadBudget,
    reader: FragmentReader,
    writer: &mut MultiWriter,
    metrics: &mut Metrics,
) -> Result<()> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(budget.matchers)
        .build()
        .context("failed to build the matcher thread pool")?;

    // Progress: a line every PROGRESS_UNIT records and a final total when this
    // counter drops (on return, including the error path), so a long run
    // reports steadily and on completion.
    let mut progress = Progress::new();

    let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<Fragment>>(READ_AHEAD_CHUNKS);
    let reader_handle = std::thread::spawn(move || read_chunks(reader, tx, READ_CHUNK));

    let mut outcome: Result<()> = Ok(());
    for chunk in &rx {
        // Match + route + extract + assemble + encode every read in the chunk
        // in parallel, preserving positional order in `prepared`; only the
        // metric tally and the (serial, ordered) writes stay on the consumer.
        let prepared: Vec<Result<Prepared>> = pool.install(|| {
            chunk
                .par_iter()
                .with_min_len(MATCHER_MIN_JOB)
                .map_init(Scratch::new, |scratch, fragment| {
                    process_fragment(engine, scratch, fragment)
                })
                .collect()
        });
        if let Err(error) = apply_chunk(prepared, writer, metrics, &mut progress) {
            outcome = Err(error);
            break;
        }
    }
    // Unblock the reader if we bailed early (closing `rx` makes its next `send`
    // fail), then join it. A consumer error is the most actionable, so it takes
    // precedence over both a read error and a reader-thread panic; absent one,
    // surface a read error or the panic.
    drop(rx);
    let reader_result = reader_handle.join();
    outcome?;
    match reader_result {
        Ok(read_outcome) => read_outcome,
        Err(_) => Err(anyhow!("the read-ahead thread panicked")),
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
/// an empty stream (one anchor unmatched) contributes no range, so its file is
/// not added to the denominator. This is a rare diagnostic-metric edge; file
/// spans and matched `@grp` anchors (the common cases) are exact.
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

/// Match a `match=` group against the concatenation of its named `--extract`
/// streams (no delimiter). A referenced stream that is absent (its anchoring
/// group did not match, or it resolved to nothing) makes the group a non-match.
/// The tags are searched within the joined text as a single synthetic segment,
/// so `dist`/`mode`/`delta`/`both_strands` apply over the concatenation.
fn match_joined(
    group: &CompiledGroup,
    scratch: &mut Scratch,
    streams: &[String],
    extracts: &[crate::grammar::Extract],
    segments: &[Segment],
    hits: &GroupHits,
    group_index: &HashMap<String, usize>,
) -> Result<GroupOutcome> {
    let mut joined = Vec::new();
    for name in streams {
        let Some(spec) = extracts.iter().find(|e| &e.name == name) else {
            return Ok(GroupOutcome::unmatched());
        };
        match extract(&spec.span, segments, hits, group_index)? {
            Some(extracted) => joined.extend_from_slice(&extracted.bases),
            None => return Ok(GroupOutcome::unmatched()),
        }
    }
    if joined.is_empty() {
        return Ok(GroupOutcome::unmatched());
    }
    Ok(match_group(group, scratch, &[joined.as_slice()]))
}

/// Match one group, honoring its sequential prerequisite: a group with no
/// prerequisite searches at its own `loc`; a downstream group whose upstream
/// did not match is skipped (unmatched); a relative `next=GROUP:lo-hi` link
/// windows the search relative to the upstream match end.
fn match_group_for(
    group: &CompiledGroup,
    scratch: &mut Scratch,
    segments: &[&[u8]],
    hits: &GroupHits,
) -> GroupOutcome {
    let Some(prereq) = &group.prereq else {
        return match_group(group, scratch, segments);
    };
    // Skip the group when any prerequisite upstream did not match (its slot is
    // `None`). Upstreams are earlier in dependency order, so their slots are
    // already filled.
    if prereq
        .upstreams
        .iter()
        .any(|&upstream| hits[upstream].is_none())
    {
        return GroupOutcome::unmatched();
    }
    // A relative `next` link windows the search relative to its upstream's
    // match (present, since the upstream is one of the now-satisfied upstreams;
    // an upstream is always a non-`match=` group, so its span is `Some`);
    // otherwise the group searches at its own `loc`.
    match &prereq.window {
        Some(window) => {
            let span = hits[window.upstream]
                .as_ref()
                .and_then(|hit| hit.span)
                .expect("a relative-window upstream is a matched non-match= group with a span");
            match_group_at(
                group,
                scratch,
                segments,
                Some(&window_after(&span, window.lo, window.hi)),
            )
        }
        None => match_group(group, scratch, segments),
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
) -> Vec<PreparedWrite> {
    let Some(pattern) = &remove.pattern else {
        return Vec::new();
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
        pattern,
        engine.pool,
        engine.plan.inputs.len(),
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
) -> Vec<PreparedWrite> {
    match &engine.plan.unassigned {
        Some(pattern) => {
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
                pattern,
                engine.pool,
                engine.plan.inputs.len(),
                fragment,
                slug.as_ref().map(|(tag, value)| (*tag, value.as_str())),
            )
        }
        None => Vec::new(),
    }
}

/// Prepare a record's raw input segments for a `--unassigned` / `--remove` bin.
/// With `%source` over multiple `--in` files the segments fan out one file per
/// input file (record `i` came from file `i`); a single input (whether plain or
/// auto-detected interleaved) is one `%source` (0), so its mates stay together.
/// Without `%source`, every segment lands in one file.
fn prepare_raw_bin(
    pattern: &OutputPattern,
    pool: &str,
    num_sources: usize,
    fragment: &Fragment,
    slug: Option<(&str, &str)>,
) -> Vec<PreparedWrite> {
    if num_sources > 1 && pattern.uses(Placeholder::Source) {
        fragment
            .records
            .iter()
            .enumerate()
            .map(|(source, record)| {
                let ctx = PathContext {
                    pool,
                    source: Some(source),
                    ..Default::default()
                };
                PreparedWrite {
                    path: resolve_pattern(pattern, &ctx),
                    payload: WritePayload::Records(vec![prepared_raw_read(record, pool, slug)]),
                }
            })
            .collect()
    } else {
        // A single source: `%source` (if present) is always 0, and every
        // segment lands together.
        let ctx = PathContext {
            pool,
            source: Some(0),
            ..Default::default()
        };
        vec![PreparedWrite {
            path: resolve_pattern(pattern, &ctx),
            payload: WritePayload::Records(
                fragment
                    .records
                    .iter()
                    .map(|record| prepared_raw_read(record, pool, slug))
                    .collect(),
            ),
        }]
    }
}

/// One raw input segment as an owned prepared record for a
/// `--unassigned`/`--remove` bin: the record's bases/quals/name and its carried
/// input tags, optionally augmented with the `--qc-tag` slug. Carries the
/// default pool read group (`RG:Z:<pool>`) so a SAM/BAM/CRAM bin matches the
/// default `@RG` in its header; FASTX output ignores the read group.
fn prepared_raw_read(record: &InputRecord, pool: &str, slug: Option<(&str, &str)>) -> PreparedRead {
    let tags = match slug {
        Some((tag, value)) => {
            let mut data = record.tags.clone().unwrap_or_default();
            insert_string_tag(&mut data, tag, value);
            Some(data)
        }
        None => record.tags.clone(),
    };
    PreparedRead {
        name: record.name.clone(),
        bases: record.bases.clone(),
        quals: record.quals.clone(),
        tags,
        read_group: Some(pool.to_string()),
    }
}

/// The canonical (corrected) tag sequence for a matched `@grp` self-span, plus
/// whether this match landed on the reverse-complement strand, or `None` when
/// the group did not match. Resolved lazily from the group's positional hit
/// (its matched tag index and strand) against the engine's group slice (the
/// tag's canonical bases), so no per-record correction buffer is carried across
/// the matcher boundary. The strand is the per-match outcome, not the group's
/// `both_strands` attribute: a `both_strands=true` group also searches the
/// forward strand, and a forward match must not have its observed quals
/// reversed.
fn canonical_correction<'a>(
    engine: &'a Engine,
    hits: &GroupHits,
    group: &str,
) -> Option<(&'a [u8], bool)> {
    let idx = *engine.group_index.get(group)?;
    let compiled = &engine.groups[idx];
    let hit = hits.get(idx)?.as_ref()?;
    let canonical = compiled.tags.get(hit.tag_idx).map(Vec::as_slice)?;
    Some((canonical, hit.revcomp))
}

/// Assemble the output body for an assigned (`target = Some`) or pass-through
/// (`target = None`) record into owned [`PreparedWrite`]s: resolve its
/// extraction streams, build the `--template` bodies and the `--tag` data,
/// attach the `@RG` read group on SAM/BAM, and resolve the `--out`
/// destination(s). Returns [`BodyOutcome::MissingStream`] when a `--template`
/// references a stream whose group did not match (the caller reroutes it to the
/// `--unassigned` bin and tallies it unassigned). Pure (no writer, no metrics),
/// so it runs on the worker.
fn prepare_body(
    engine: &Engine,
    target: Option<&Target>,
    hits: &GroupHits,
    segments: &[Segment],
    fragment: &Fragment,
) -> Result<BodyOutcome> {
    let plan = engine.plan;
    let pool = engine.pool;
    let input_formats = engine.input_formats;
    let mut streams: HashMap<&str, Extracted> = HashMap::new();
    for spec in &plan.extracts {
        if let Some(mut extracted) = extract(&spec.span, segments, hits, engine.group_index)? {
            // Attach the matched tag's canonical sequence as the stream's
            // corrected form, so a `raw=false` `--tag`/`--template` emits the
            // corrected barcode while a `raw=true` one (or any stream with no
            // corrected form) emits the observed bases. Only the `@grp`
            // self-span is corrected (anchored-offset and between-anchor spans
            // are unaffected). The corrected form is always the declared
            // canonical (matching=only `both_strands=true`); orientation is a
            // use-site choice (`--tag T=~stream` / `--template ~stream`).
            if let SpanBody::AnchorMatch { group } = &spec.span.body {
                if let Some((canonical, antisense)) = canonical_correction(engine, hits, group) {
                    extracted.corrected = Some(corrected_form(&extracted, canonical, antisense));
                }
            }
            streams.insert(spec.name.as_str(), extracted);
        }
    }

    let bodies = match assemble_bodies(plan, fragment, &mut streams, engine.movable_streams) {
        Ok(bodies) => bodies,
        // A templated stream is absent (its anchoring group did not match): the
        // caller reroutes to unassigned, naming the missing stream in the QC
        // slug.
        Err(missing) => return Ok(BodyOutcome::MissingStream(missing)),
    };

    let tag_data = build_tag_data(&plan.tags, &streams);
    // An assigned read carries its sample's read group; a pass-through read (no
    // target) carries the default pool read group, so a SAM/BAM/CRAM `--out`
    // matches the default `@RG` in its header.
    let read_group = target.map(Target::label).or_else(|| Some(pool.to_string()));
    let sample = target.map(|t| t.sample.as_str());
    let sub_sample = target.and_then(|t| t.sub_sample.as_deref());

    let mut datas: Vec<Data> = bodies
        .iter()
        .map(|body| merge_data(body.carried.as_ref(), &tag_data))
        .collect();
    // With `--qc-tag`, every output body of this record carries the
    // routing-outcome slug (assigned with its sample, or pass-through). Built
    // once from the read's group hits and observed bases.
    if let Some(tag) = &plan.qc_tag {
        let base_segments: Vec<&[u8]> = segments.iter().map(|segment| segment.bases).collect();
        let slug = routed_slug(engine.groups, hits, &base_segments, sample, sub_sample);
        for data in &mut datas {
            insert_string_tag(data, tag, &slug);
        }
    }

    // The format is fixed by the path extension (or the stdout mirror),
    // independent of the ordinal.
    let probe = PathContext {
        pool,
        sample,
        sub_sample,
        ordinal: Some(1),
        source: None,
    };
    let probe_path = plan.out.as_ref().and_then(|p| resolve_pattern(p, &probe));
    let format = output_format(probe_path.as_deref(), input_formats)?;
    let is_alignment = matches!(
        format,
        OutputFormat::Sam | OutputFormat::Bam | OutputFormat::Cram
    );

    // Multi-read FASTX to a file fans out one file per `%ordinal`; SAM/BAM (a
    // pair in one file), single-read FASTX, and stdout all write every read to
    // one destination.
    let multi_fastx = (!is_alignment && bodies.len() > 1)
        .then_some(plan.out.as_ref())
        .flatten();
    let writes = if let Some(pattern) = multi_fastx {
        if !pattern.uses(Placeholder::Ordinal) {
            bail!(
                "unmux: multi-read FASTX output needs the %ordinal placeholder in --out to \
                 separate the {} template reads",
                bodies.len()
            );
        }
        bodies
            .into_iter()
            .zip(datas)
            .enumerate()
            .map(|(index, (body, data))| {
                let ctx = PathContext {
                    pool,
                    sample,
                    sub_sample,
                    ordinal: Some(index + 1),
                    source: None,
                };
                PreparedWrite {
                    path: resolve_pattern(pattern, &ctx),
                    payload: WritePayload::Records(vec![prepared_body_read(body, data, None)]),
                }
            })
            .collect()
    } else {
        let ctx = PathContext {
            pool,
            sample,
            sub_sample,
            ordinal: (bodies.len() == 1).then_some(1),
            source: None,
        };
        let path = plan.out.as_ref().and_then(|p| resolve_pattern(p, &ctx));
        let reads = bodies
            .into_iter()
            .zip(datas)
            .map(|(body, data)| {
                let read_group = if is_alignment {
                    read_group.clone()
                } else {
                    None
                };
                prepared_body_read(body, data, read_group)
            })
            .collect();
        vec![PreparedWrite {
            path,
            payload: WritePayload::Records(reads),
        }]
    };
    Ok(BodyOutcome::Assembled(writes))
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
    }
}

/// One assembled output record: its name, bases, qualities, and (for a raw
/// pass-through body) the carried input tags it should keep.
struct Body {
    name: Vec<u8>,
    bases: Vec<u8>,
    quals: Option<Vec<u8>>,
    carried: Option<Data>,
}

/// Assemble a record's output bodies: the raw input segments (each carrying its
/// tags) when there is no `--template`, else one body per `--template`. Returns
/// `Err(missing_stream_name)` when a template references a stream that was not
/// extracted (its anchoring group did not match), naming the absent stream for
/// the `--qc-tag` slug.
fn assemble_bodies(
    plan: &DemuxPlan,
    fragment: &Fragment,
    streams: &mut HashMap<&str, Extracted>,
    movable: &HashSet<&str>,
) -> std::result::Result<Vec<Body>, String> {
    if plan.templates.is_empty() {
        return Ok(fragment
            .records
            .iter()
            .map(|record| Body {
                name: record.name.clone(),
                bases: record.bases.clone(),
                quals: record.quals.clone(),
                carried: record.tags.clone(),
            })
            .collect());
    }
    let name = fragment
        .records
        .first()
        .map(|record| record.name.clone())
        .unwrap_or_default();
    let mut bodies = Vec::with_capacity(plan.templates.len());
    for template in &plan.templates {
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
        bodies.push(Body {
            name: name.clone(),
            bases,
            quals,
            carried: None,
        });
    }
    Ok(bodies)
}

/// Merge the carried input tags (if any) with the assembled `--tag` fields,
/// with `--tag` winning a key collision.
fn merge_data(carried: Option<&Data>, tag_data: &Data) -> Data {
    let mut data = carried.cloned().unwrap_or_default();
    for (tag, value) in tag_data.iter() {
        data.insert(tag, value.clone());
    }
    data
}

/// Build the corrected form of a matched `@grp` stream: always the declared
/// canonical sequence, in the orientation the user supplied it (orientation is
/// a use-site choice via `~`, not baked in here). Degenerate / IUPAC positions
/// are emitted verbatim. The observed qualities are copied and resized to the
/// corrected length; for an antisense match the observed quals are in
/// read-strand order (the reverse complement of the declared orientation), so
/// they are reversed to co-orient with the declared bases. Non-destructive: the
/// observed bases/quals on `extracted` are left intact for the record body and
/// for a `raw=true` tag/template.
fn corrected_form(extracted: &Extracted, canonical: &[u8], antisense: bool) -> Corrected {
    let bases = canonical.to_vec();
    let quals = extracted.quals.as_ref().map(|quals| {
        let mut quals = quals.clone();
        if antisense {
            quals.reverse();
        }
        fit_quals(quals, bases.len())
    });
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
    fn test_partition_threads_reserves_roles() {
        // partition_threads is only consulted for N >= 3 (drive_parallel); 1
        // and 2 use the serial drivers. No pool: every spare thread is a
        // matcher.
        for threads in [3, 4, 5, 6, 7, 8, 16] {
            let budget = partition_threads(threads, 0);
            assert_eq!(budget.reader, 1, "one read-ahead thread");
            assert_eq!(budget.main, 1, "one main/consumer thread");
            assert_eq!(budget.writer, 0, "no BGZF pool");
            assert!(budget.matchers >= 1, "at least one matcher");
            assert_eq!(
                budget.reader + budget.main + budget.writer + budget.matchers,
                threads,
                "the partition accounts for every thread"
            );
        }
        // With a BGZF pool, the writer threads come out of the matcher budget.
        for (threads, writer) in [(4, 1), (8, 2), (16, 4)] {
            let budget = partition_threads(threads, writer);
            assert_eq!(budget.writer, writer);
            assert!(budget.matchers >= 1, "at least one matcher remains");
            assert_eq!(
                budget.reader + budget.main + budget.writer + budget.matchers,
                threads
            );
        }
    }

    #[test]
    fn test_writer_pool_size_even_split() {
        // No compressed output: no BGZF pool, every spare thread is a matcher.
        assert_eq!(writer_pool_size(8, false), 0);
        // Compressed output: spare threads split evenly with the matchers. At
        // --threads 8 this is the measured optimum writer=3/matchers=3
        // (matchers never below the writer count).
        assert_eq!(writer_pool_size(8, true), 3);
        let budget = partition_threads(8, writer_pool_size(8, true));
        assert_eq!((budget.writer, budget.matchers), (3, 3));
        // Too few threads to afford a pool (reader + main + >=1 writer + >=1
        // matcher needs 4): the pool is empty and compression runs inline on
        // the consumer.
        assert_eq!(writer_pool_size(2, true), 0);
        assert_eq!(writer_pool_size(3, true), 0);
        // Small thread counts still reserve a matcher; high counts are capped
        // at MAX_WRITER_THREADS.
        assert_eq!(writer_pool_size(4, true), 1);
        assert_eq!(writer_pool_size(16, true), MAX_WRITER_THREADS);
        assert!(writer_pool_size(64, true) <= MAX_WRITER_THREADS);
    }

    #[test]
    fn test_low_thread_output_matches_high_thread() {
        // --threads 1 (serial) and 2 (pipelined) must produce the same
        // demultiplexed records as a high-thread run. Plain FASTQ is
        // byte-identical across every thread count. Compressed output (inline
        // BAM at 1/2/3 vs pooled BGZF at 16) can frame BGZF blocks differently,
        // so those compare decoded records via read_back rather than raw bytes.
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

        // Compressed BAM: compare decoded records (inline vs pooled BGZF
        // framing may differ).
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
        // itself stays orientation-neutral.
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
        // A `both_strands=true` group matches antisense, but the corrected tag
        // is always the declared canonical (both_strands is matching-only); the
        // raw CR still shows the observed read-strand bases.
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
            "g::loc=0:10:16,both_strands=true,dist=1".to_string(),
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
        // CR (raw): the observed read bases CACGAT (a 1-mismatch instance of
        // the rc, read strand).
        assert_eq!(string_tag(&record, *b"CR").as_deref(), Some("CACGAT"));
    }

    #[test]
    fn test_both_strands_corrected_quals_co_orient_with_declared() {
        // The corrected qualities must co-orient with the DECLARED canonical
        // bases. For an antisense match the observed quals are in read-strand
        // order (the rc of the declared orientation), so they are reversed.
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
            "g::loc=0:10:16,both_strands=true,dist=1".to_string(),
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
        // A `both_strands=true` group also searches the forward strand. On a
        // forward match the observed quals already co-orient with the declared
        // canonical, so they must NOT be reversed: the per-match strand governs
        // the correction, not the group's `both_strands` attribute.
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
            "g::loc=0:10:16,both_strands=true,dist=1".to_string(),
        ];
        args.extracts = vec!["bc=@g".to_string()];
        // CB defaults to the CY quality tag, carrying the corrected qualities.
        args.tags = vec!["CB=bc".to_string()];
        args.out = Some(out.display().to_string());
        run_demux(args).unwrap();

        let record = first_bam_record(&out);
        assert_eq!(string_tag(&record, *b"CB").as_deref(), Some("AACGTG"));
        // CY stays "56789:" (read order); reversing it to ":98765" would be the
        // bug of using the group's `both_strands` config instead of the match
        // strand.
        assert_eq!(string_tag(&record, *b"CY").as_deref(), Some("56789:"));
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
    fn test_pooled_fastx_gz_is_bgzf_and_round_trips() {
        // A gzipped-FASTX fan-out goes through the BGZF compressor pool: each
        // .fq.gz is BGZF-framed (valid gzip), round-trips its reads, and a
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
        // Per the spec, --compression sets the BAM BGZF level. The same data at
        // level 0 (stored) must produce a clearly larger file than at level 9
        // (max), and both must round-trip every record.
        let dir = tempfile::tempdir().unwrap();
        let bases = [b'A', b'C', b'G', b'T'];
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
    fn test_pooled_bam_unassigned_bin_round_trips() {
        // A BAM `--unassigned` bin is pooled (file-backed BAM), so it exercises
        // the pooled BAM writer's bytes directly: it must be a valid BAM that
        // reads back, proving the pooled bam writer carries no double-BGZF and
        // the EOF is sound.
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
    fn test_pooled_error_path_returns_clean_error_not_panic() {
        // A consumer write error while a BGZF pool is live must surface a clean
        // Err, not panic via the pool/PooledWriter drop order (the pool is
        // stopped only after every writer is finalized). The run errors because
        // a 3-input pass-through emits 3 reads per record, which alignment
        // output cannot hold; a poolable `.fq.gz` unassigned bin forces the
        // pool to be built.
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
    fn test_pooled_fastx_gz_deterministic_across_threads() {
        // BGZF output is byte-identical regardless of --threads: the consumer
        // feeds the pool in input order and BGZF block compression is
        // deterministic at a fixed level.
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
        // Matching fans across the rayon pool, but the consumer applies results
        // in input order, so a single-file --out is byte-identical regardless
        // of --threads. The input spans several READ_CHUNK chunks so
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
        // The encode-on-workers SAM/BAM path: records are encoded on the rayon
        // matcher workers and appended by the consumer in input order, so a
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

        // 4 is the pooled-BGZF floor: read-ahead + main + >=1 writer + >=1
        // matcher. Below it (Task 2), compression runs inline; here both runs
        // use the pool so the bytes are identical.
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
}
