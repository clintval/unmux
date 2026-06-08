//! Sample fan-out: turn the parsed selectors into a routing table that assigns
//! each record to a fan-out target (a `(sample, sub_sample)` pair), removes
//! skiplisted records, or leaves the rest unassigned.
//!
//! The grammar parser has already checked that every selector names a declared
//! group; what it could not check is anything that needs the group's loaded tag
//! set. This module closes that gap: it resolves each selector token to a
//! concrete tag (by id or by sequence), enforces the one-tag-one-target
//! collision rule, optionally requires full coverage of a group's tags, and
//! parses the table form of `--sample` from a sheet. Routing itself is a pure
//! function over a record's per-group best matches, so it is unit-testable
//! without any I/O.
//!
//! The fan-out source is normalized to a flat list of [`Sample`] targets before
//! compiling: an inline `--sample` list is used as-is, a `--sample-sheet` is
//! parsed by [`parse_sample_sheet`], and a `--sample-from-group` is expanded by
//! [`expand_from_group`]. The pipeline glue does that normalization and then
//! calls [`compile_routing`]; keeping the steps separate keeps each one
//! testable on its own.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};

use crate::extract::GroupHits;
use crate::grammar::{
    GroupSelector, OutputPattern, RemoveRule, RemoveSelector, Sample, Selector, SubSample,
};
use crate::tags::TagSet;

/// A resolved fan-out target: the sample (read-group `SM`) and an optional
/// sub_sample (read-group `LB`). The fan-out key is the whole `(sample,
/// sub_sample)` pair, so one sample may appear both with and without a
/// sub_sample.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Target {
    /// The sample id (`SM`).
    pub sample: String,
    /// The sub_sample id (`LB`), already resolved (`%pool` expanded to the pool
    /// id); `None` when the sample has no sub_sample.
    pub sub_sample: Option<String>,
}

impl Target {
    /// The display label `sample[.sub_sample]`, matching the read-group id
    /// layout.
    pub fn label(&self) -> String {
        match &self.sub_sample {
            Some(sub) => format!("{}.{}", self.sample, sub),
            None => self.sample.clone(),
        }
    }
}

/// A compiled `--remove` rule: the group whose match triggers removal (by
/// canonical index), an optional specific tag (the resolved `group::id`), and
/// the optional output pattern the removed records are written to as raw input
/// segments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemoveTarget {
    /// The canonical index of the group whose match triggers removal.
    pub group: usize,
    /// The specific tag index to remove on, or `None` to remove on any tag of
    /// the group (a bare group).
    pub tag: Option<usize>,
    /// Where to write the removed records, or `None` to drop them.
    pub pattern: Option<OutputPattern>,
    /// The user's selector echoed verbatim (`group` or `group::id`), for the
    /// `--qc-tag` removed slug.
    pub selector: String,
}

/// What the router decided for one record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Disposition<'a> {
    /// Routed to a fan-out target.
    Assigned(&'a Target),
    /// Matched a `--remove` rule; removed and tallied as `removed`, and written
    /// to the rule's pattern when it has one.
    Removed(&'a RemoveTarget),
    /// Demux is configured but the record matched no sample (the `--unassigned`
    /// bin).
    Unassigned,
    /// No demux is configured; the whole pool passes through to `--out`.
    PassThrough,
}

/// One AND term of a route: the record's best match in `group` (the group's
/// canonical index) must be one of `members`, or any tag of the group when
/// `members` is `None` (a bare-group selector).
#[derive(Debug, Clone)]
struct Term {
    group: usize,
    members: Option<Vec<usize>>,
}

/// One fan-out route: a target plus the AND-joined terms a record must satisfy
/// to reach it.
#[derive(Debug, Clone)]
struct Route {
    target: Target,
    terms: Vec<Term>,
}

/// The compiled routing table consumed by the engine: an ordered set of routes,
/// the resolved remove rules, and the distinct fan-out targets (for the writer
/// to pre-create outputs and `@RG` lines).
#[derive(Debug, Clone)]
pub struct Routing {
    /// Whether any sample fan-out is configured; `false` is a pure pass-through
    /// run.
    demux: bool,
    /// Routes in declaration order; the first a record satisfies wins.
    routes: Vec<Route>,
    /// Resolved `--remove` rules, checked before sample assignment.
    removes: Vec<RemoveTarget>,
    /// The distinct fan-out targets, in first-seen order.
    targets: Vec<Target>,
}

impl Routing {
    /// Route one record given its per-group best matches, indexed positionally
    /// by canonical group index: a `Some` slot carries that group's single
    /// best-matched tag (an absent / `None` slot did not match). `--remove`
    /// rules are checked first and take precedence over sample assignment.
    pub fn route(&self, hits: &GroupHits) -> Disposition<'_> {
        for remove in &self.removes {
            if term_satisfied(
                hits,
                remove.group,
                remove.tag.as_ref().map(std::slice::from_ref),
            ) {
                return Disposition::Removed(remove);
            }
        }
        if !self.demux {
            return Disposition::PassThrough;
        }
        for route in &self.routes {
            if route
                .terms
                .iter()
                .all(|term| term_satisfied(hits, term.group, term.members.as_deref()))
            {
                return Disposition::Assigned(&route.target);
            }
        }
        Disposition::Unassigned
    }

    /// The distinct fan-out targets in first-seen order (one `@RG` / output bin
    /// each).
    pub fn targets(&self) -> &[Target] {
        &self.targets
    }
}

/// Whether the record's best match in `group` (by canonical index) satisfies a
/// membership constraint: any matched tag when `members` is `None`, else the
/// best tag index must appear in `members`. A `None` slot (or an index past the
/// end) means the group did not match.
fn term_satisfied(hits: &GroupHits, group: usize, members: Option<&[usize]>) -> bool {
    match hits.get(group).and_then(Option::as_ref) {
        None => false,
        Some(hit) => members.is_none_or(|m| m.contains(&hit.tag_idx)),
    }
}

/// Build a routing table from the normalized fan-out targets and remove rules.
///
/// `samples` is the fully expanded target list (inline `--sample` as-is, a
/// sheet already parsed by [`parse_sample_sheet`], or `--sample-from-group`
/// already expanded by [`expand_from_group`]). `demux` is `false` only for a
/// pure pass-through run (no sample flags), where every non-removed record
/// flows to `--out`. `tag_sets` must hold the loaded set for every group
/// referenced by a selector or a remove rule. `group_index` maps each group
/// name to its canonical index (its position in the engine's group slice), so a
/// compiled route/remove records the index the per-record router indexes
/// positionally. `pool` resolves `%pool` sub_samples. With
/// `require_samples_explain_all_tags`, every tag of every group named by a
/// selector must be claimed by some target, else this fails fast (the caller
/// passes `false` for `--sample-from-group`, which is 1:1 by construction).
///
/// The compile-time collision and coverage checks stay keyed by group NAME and
/// emit human-readable name-based errors; only the per-record
/// [`Routing::route`] path is index-based.
pub fn compile_routing(
    samples: &[Sample],
    removes: &[RemoveRule],
    tag_sets: &HashMap<String, TagSet>,
    group_index: &HashMap<String, usize>,
    pool: &str,
    demux: bool,
    require_samples_explain_all_tags: bool,
) -> Result<Routing> {
    let mut routes = Vec::new();
    let mut targets = Vec::new();
    // The one-tag-one-target collision map: a unit is the AND-tuple of `(group,
    // tag)` claims; an OR pool and a bare group each expand to one unit per
    // tag. The same unit claimed by two distinct targets is the collision we
    // fail on.
    let mut claims: HashMap<Vec<(String, usize)>, Target> = HashMap::new();

    for sample in samples {
        let target = resolve_target(sample, pool);
        if !targets.contains(&target) {
            targets.push(target.clone());
        }
        let terms = compile_terms(&sample.selector, tag_sets, group_index)?;
        for unit in claim_units(&sample.selector, &terms, tag_sets) {
            match claims.get(&unit) {
                Some(existing) if *existing != target => {
                    bail!(
                        "{} is routed to more than one sample: `{}` and `{}`",
                        render_unit(&unit, tag_sets),
                        existing.label(),
                        target.label()
                    );
                }
                Some(_) => {} // the same target re-claiming a unit (an OR pool) is fine.
                None => {
                    claims.insert(unit, target.clone());
                }
            }
        }
        routes.push(Route { target, terms });
    }

    if require_samples_explain_all_tags {
        check_coverage(samples, &claims, tag_sets)?;
    }

    let removes = removes
        .iter()
        .map(|rule| compile_remove(rule, tag_sets, group_index))
        .collect::<Result<Vec<_>>>()?;

    Ok(Routing {
        demux,
        routes,
        removes,
        targets,
    })
}

/// Resolve a parsed sample's `(sample, sub_sample)` into a concrete [`Target`],
/// expanding `%pool`.
fn resolve_target(sample: &Sample, pool: &str) -> Target {
    let sub_sample = match &sample.sub_sample {
        None => None,
        Some(SubSample::Literal(sub)) => Some(sub.clone()),
        Some(SubSample::Pool) => Some(pool.to_string()),
    };
    Target {
        sample: sample.sample.clone(),
        sub_sample,
    }
}

/// Resolve a selector's terms against the loaded tag sets, turning each token
/// into a tag index and each group name into its canonical index (for the
/// per-record router). A bare-group term (no listed members) becomes a
/// whole-group term (`None`); a token that names neither a tag id nor a tag
/// sequence in its group is a fail-fast error.
fn compile_terms(
    selector: &Selector,
    tag_sets: &HashMap<String, TagSet>,
    group_index: &HashMap<String, usize>,
) -> Result<Vec<Term>> {
    let mut terms = Vec::new();
    for term in &selector.terms {
        // Every referenced group must have its tag set loaded, bare-group
        // selectors included: a missing set is a caller contract violation, and
        // resolving it here keeps a bare group from silently claiming zero
        // units (which would disable its collision and coverage checks).
        let set = require_tag_set(tag_sets, &term.group)?;
        let members = if term.members.is_empty() {
            None
        } else {
            let mut idxs = Vec::new();
            for token in &term.members {
                let idx = set.resolve(token).ok_or_else(|| {
                    anyhow!(
                        "selector references `{}::{token}`, but `{token}` is not a tag id or sequence in group `{}`",
                        term.group,
                        term.group
                    )
                })?;
                if !idxs.contains(&idx) {
                    idxs.push(idx);
                }
            }
            Some(idxs)
        };
        terms.push(Term {
            group: require_group_index(group_index, &term.group)?,
            members,
        });
    }
    Ok(terms)
}

/// Enumerate the selector units a route claims: the cartesian product over its
/// terms of each term's claimed tag indices (a bare-group term claims every tag
/// of the group). Each unit is canonicalized by sorting its `(group, tag)`
/// pairs, so units compare equal regardless of term order. Keyed by group NAME
/// (taken from the original `selector`, term-aligned with `terms`) so the
/// collision and coverage errors stay human-readable; `terms` only supplies the
/// resolved member indices.
fn claim_units(
    selector: &Selector,
    terms: &[Term],
    tag_sets: &HashMap<String, TagSet>,
) -> Vec<Vec<(String, usize)>> {
    let per_term: Vec<(String, Vec<usize>)> = selector
        .terms
        .iter()
        .zip(terms)
        .map(|(selector_term, term)| {
            let group = selector_term.group.clone();
            let idxs = match &term.members {
                Some(members) => members.clone(),
                None => (0..tag_sets.get(&group).map_or(0, |s| s.entries.len())).collect(),
            };
            (group, idxs)
        })
        .collect();

    let mut units: Vec<Vec<(String, usize)>> = vec![Vec::new()];
    for (group, idxs) in &per_term {
        let mut next = Vec::new();
        for prefix in &units {
            for &idx in idxs {
                let mut unit = prefix.clone();
                unit.push((group.clone(), idx));
                next.push(unit);
            }
        }
        units = next;
    }
    for unit in &mut units {
        unit.sort();
    }
    units
}

/// Render a selector unit for an error message, naming each claim as
/// `group::tag-id`. A single claim reads as `tag \`group::id\``; an AND-tuple
/// reads as `the combination \`a::x+b::y\``.
fn render_unit(unit: &[(String, usize)], tag_sets: &HashMap<String, TagSet>) -> String {
    let parts: Vec<String> = unit
        .iter()
        .map(|(group, idx)| {
            let id = tag_sets
                .get(group)
                .and_then(|set| set.entries.get(*idx))
                .map_or("?", |entry| entry.id.as_str());
            format!("{group}::{id}")
        })
        .collect();
    if parts.len() == 1 {
        format!("tag `{}`", parts[0])
    } else {
        format!("the combination `{}`", parts.join("+"))
    }
}

/// Require full coverage: every tag of every group named by a selector must
/// appear in some claimed unit, else fail fast and name the unclaimed tags.
fn check_coverage(
    samples: &[Sample],
    claims: &HashMap<Vec<(String, usize)>, Target>,
    tag_sets: &HashMap<String, TagSet>,
) -> Result<()> {
    let mut referenced: HashSet<&str> = HashSet::new();
    for sample in samples {
        for term in &sample.selector.terms {
            referenced.insert(term.group.as_str());
        }
    }
    let mut claimed: HashSet<(&str, usize)> = HashSet::new();
    for unit in claims.keys() {
        for (group, idx) in unit {
            claimed.insert((group.as_str(), *idx));
        }
    }
    let mut unclaimed = Vec::new();
    for group in referenced {
        if let Some(set) = tag_sets.get(group) {
            for (idx, entry) in set.entries.iter().enumerate() {
                if !claimed.contains(&(group, idx)) {
                    unclaimed.push(format!("{group}::{}", entry.id));
                }
            }
        }
    }
    if !unclaimed.is_empty() {
        bail!(
            "--require-samples-explain-all-tags: these tags are claimed by no sample: {}",
            unclaimed.join(", ")
        );
    }
    Ok(())
}

/// Resolve a `--remove` rule's selector against its group tag set: a bare group
/// removes on any of its tags; a `group::id` removes on that one tag, resolved
/// by id or sequence. The group name is resolved to its canonical index for the
/// per-record router.
fn compile_remove(
    rule: &RemoveRule,
    tag_sets: &HashMap<String, TagSet>,
    group_index: &HashMap<String, usize>,
) -> Result<RemoveTarget> {
    let (group, tag) = match &rule.selector {
        RemoveSelector::Group(group) => (group.as_str(), None),
        RemoveSelector::GroupTag { group, id } => {
            let set = require_tag_set(tag_sets, group)?;
            let idx = set.resolve(id).ok_or_else(|| {
                anyhow!(
                    "--remove references `{group}::{id}`, but `{id}` is not a tag id or sequence in group `{group}`"
                )
            })?;
            (group.as_str(), Some(idx))
        }
    };
    let selector = match &rule.selector {
        RemoveSelector::Group(group) => group.clone(),
        RemoveSelector::GroupTag { group, id } => format!("{group}::{id}"),
    };
    Ok(RemoveTarget {
        group: require_group_index(group_index, group)?,
        tag,
        pattern: rule.pattern.clone(),
        selector,
    })
}

/// Look up a group's loaded tag set, with a guard error if the caller forgot to
/// load it.
fn require_tag_set<'a>(tag_sets: &'a HashMap<String, TagSet>, group: &str) -> Result<&'a TagSet> {
    tag_sets
        .get(group)
        .with_context(|| format!("no loaded tag set for group `{group}`"))
}

/// Look up a group's canonical index, with a guard error if the caller forgot
/// to map it.
fn require_group_index(group_index: &HashMap<String, usize>, group: &str) -> Result<usize> {
    group_index
        .get(group)
        .copied()
        .with_context(|| format!("no canonical index for group `{group}`"))
}

/// Expand `--sample-from-group GROUP` into one [`Sample`] per tag: `SM` is the
/// tag id, `LB` is the tag's `sub_sample` column (if any), and the selector
/// routes that single tag to that sample. This is 1:1 by construction, so it
/// never collides and needs no coverage check.
pub fn expand_from_group(group: &str, tag_set: &TagSet) -> Vec<Sample> {
    tag_set
        .entries
        .iter()
        .map(|entry| Sample {
            sample: entry.id.clone(),
            sub_sample: entry.sub_sample.clone().map(SubSample::Literal),
            selector: Selector {
                terms: vec![GroupSelector {
                    group: group.to_string(),
                    members: vec![entry.id.clone()],
                }],
            },
        })
        .collect()
}

/// Parse a `--sample-sheet` TSV into the same [`Sample`] targets `--sample`
/// produces.
///
/// Reserved columns are `sample` (required) and `sub_sample` (optional); every
/// other non-empty column header must name a declared group, and its cells
/// select a tag in that group. Multiple group columns on a row are AND-joined;
/// rows sharing `(sample, sub_sample)` form an OR pool (so, unlike inline
/// `--sample`, a repeated target is allowed and accumulates members). The
/// tokens themselves are resolved against the group tag sets later, at
/// routing-compile time. `#` comments and blank lines are skipped and CRLF line
/// endings are tolerated.
pub fn parse_sample_sheet(text: &str, declared_groups: &HashSet<&str>) -> Result<Vec<Sample>> {
    let mut rows = text
        .lines()
        .map(|line| line.trim_end_matches('\r'))
        .filter(|line| !line.trim().is_empty() && !line.trim_start().starts_with('#'));

    let header = rows
        .next()
        .context("sample sheet is empty (a header row is required)")?;
    let columns: Vec<&str> = header.split('\t').map(str::trim).collect();

    let sample_col = columns
        .iter()
        .position(|c| *c == "sample")
        .context("sample sheet has no `sample` column")?;
    let sub_sample_col = columns.iter().position(|c| *c == "sub_sample");

    // Every remaining non-empty column must name a declared group; its cells
    // select within it.
    let mut group_cols = Vec::new();
    let mut seen_groups = HashSet::new();
    for (i, column) in columns.iter().enumerate() {
        if i == sample_col || Some(i) == sub_sample_col || column.is_empty() {
            continue;
        }
        if !declared_groups.contains(column) {
            bail!(
                "sample sheet column `{column}` is not a declared --group (every non-reserved column must name a group)"
            );
        }
        // One column per group: a repeated group column AND-joins a group with
        // itself, an unsatisfiable route (OR pools come from repeated rows, not
        // repeated columns).
        if !seen_groups.insert(*column) {
            bail!(
                "sample sheet names group `{column}` in more than one column; build an OR pool with repeated rows, not repeated columns"
            );
        }
        group_cols.push((i, *column));
    }

    let mut samples = Vec::new();
    for (line_no, row) in rows.enumerate() {
        let cells: Vec<&str> = row.split('\t').collect();
        let sample = cells.get(sample_col).map(|c| c.trim()).unwrap_or("");
        if sample.is_empty() {
            bail!("sample sheet row {} has an empty `sample`", line_no + 1);
        }
        let sub_sample = sub_sample_col
            .and_then(|c| cells.get(c))
            .map(|c| c.trim())
            .filter(|s| !s.is_empty())
            .map(|s| {
                if s == "%pool" {
                    SubSample::Pool
                } else {
                    SubSample::Literal(s.to_string())
                }
            });

        let mut terms = Vec::new();
        for (i, group) in &group_cols {
            let cell = cells.get(*i).map(|c| c.trim()).unwrap_or("");
            if cell.is_empty() {
                continue;
            }
            terms.push(GroupSelector {
                group: group.to_string(),
                members: vec![cell.to_string()],
            });
        }
        if terms.is_empty() {
            bail!(
                "sample sheet row {} selects no barcode in any group column",
                line_no + 1
            );
        }
        samples.push(Sample {
            sample: sample.to_string(),
            sub_sample,
            selector: Selector { terms },
        });
    }
    if samples.is_empty() {
        bail!("sample sheet has a header but no rows");
    }
    Ok(samples)
}

/// Load and parse a `--sample-sheet` from a file.
pub fn load_sample_sheet(path: &Path, declared_groups: &HashSet<&str>) -> Result<Vec<Sample>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read sample sheet: {}", path.display()))?;
    parse_sample_sheet(&text, declared_groups)
        .with_context(|| format!("failed to parse sample sheet: {}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::extract::{GroupHit, MatchSpan};
    use crate::grammar::{PatternSegment, Placeholder};
    use crate::tags::TagEntry;

    /// Build a tag set from `(id, seq)` pairs (no sub_samples).
    fn tagset(pairs: &[(&str, &str)]) -> TagSet {
        TagSet {
            entries: pairs
                .iter()
                .map(|(id, seq)| TagEntry {
                    id: id.to_string(),
                    seq: seq.to_string(),
                    sub_sample: None,
                })
                .collect(),
        }
    }

    /// A `tag_sets` map from `(group, TagSet)` pairs.
    fn sets(pairs: &[(&str, TagSet)]) -> HashMap<String, TagSet> {
        pairs
            .iter()
            .map(|(g, s)| (g.to_string(), s.clone()))
            .collect()
    }

    /// A deterministic canonical `name -> index` map for a set of groups,
    /// mirroring the engine's positional indexing (the production map is built
    /// from the final group order; here, sorted names give a stable order
    /// independent of `HashMap` iteration). Both [`compile`] and [`matched`]
    /// use this same derivation so route indices line up.
    fn group_index(tag_sets: &HashMap<String, TagSet>) -> HashMap<String, usize> {
        let mut names: Vec<&String> = tag_sets.keys().collect();
        names.sort();
        names
            .into_iter()
            .enumerate()
            .map(|(idx, name)| (name.clone(), idx))
            .collect()
    }

    /// Compile a routing table for a test, deriving the canonical group index
    /// from the tag sets, and return both so [`matched`] can build
    /// positionally-aligned hits.
    fn compile(
        samples: &[Sample],
        removes: &[RemoveRule],
        tag_sets: &HashMap<String, TagSet>,
        pool: &str,
        demux: bool,
        require_samples_explain_all_tags: bool,
    ) -> Result<(Routing, HashMap<String, usize>)> {
        let index = group_index(tag_sets);
        let routing = compile_routing(
            samples,
            removes,
            tag_sets,
            &index,
            pool,
            demux,
            require_samples_explain_all_tags,
        )?;
        Ok((routing, index))
    }

    /// Build a sample with an optional sub_sample and AND terms given as
    /// `(group, &[member tokens])` (an empty member slice is a bare-group /
    /// whole-group term).
    fn sample(name: &str, sub: Option<&str>, terms: &[(&str, &[&str])]) -> Sample {
        Sample {
            sample: name.to_string(),
            sub_sample: sub.map(|s| SubSample::Literal(s.to_string())),
            selector: Selector {
                terms: terms
                    .iter()
                    .map(|(group, members)| GroupSelector {
                        group: group.to_string(),
                        members: members.iter().map(|m| m.to_string()).collect(),
                    })
                    .collect(),
            },
        }
    }

    /// A record's per-group best matches as positional hits, from `(group,
    /// tag_idx)` pairs and the canonical group index. A group absent from the
    /// pairs gets a `None` slot (it did not match); matched groups carry a
    /// `span` (a normal record-window group), which a `match=` group would not,
    /// but these routing tests never anchor on the span, so a placeholder span
    /// is fine.
    fn matched(index: &HashMap<String, usize>, pairs: &[(&str, usize)]) -> Vec<Option<GroupHit>> {
        let mut hits: Vec<Option<GroupHit>> = vec![None; index.len()];
        for (group, tag_idx) in pairs {
            if let Some(&idx) = index.get(*group) {
                hits[idx] = Some(GroupHit {
                    tag_idx: *tag_idx,
                    span: Some(MatchSpan {
                        file: 0,
                        start: 0,
                        end: 0,
                    }),
                    subs: 0,
                    indels: 0,
                    revcomp: false,
                });
            }
        }
        hits
    }

    /// A `match=`-style hit for `group` (matched but anchoring nothing: `span:
    /// None`), used to assert such a group still routes while leaving its span
    /// absent.
    fn matched_no_span(
        index: &HashMap<String, usize>,
        group: &str,
        tag_idx: usize,
    ) -> Vec<Option<GroupHit>> {
        let mut hits: Vec<Option<GroupHit>> = vec![None; index.len()];
        if let Some(&idx) = index.get(group) {
            hits[idx] = Some(GroupHit {
                tag_idx,
                span: None,
                subs: 0,
                indels: 0,
                revcomp: false,
            });
        }
        hits
    }

    #[test]
    fn test_route_single_group_or_pool() {
        // dna01 claims cbt01/cbt02; dna02 claims cbt03. Routing partitions one
        // group's tags.
        let tag_sets = sets(&[(
            "grp",
            tagset(&[("cbt01", "AAAA"), ("cbt02", "CCCC"), ("cbt03", "GGGG")]),
        )]);
        let samples = vec![
            sample("dna01", None, &[("grp", &["cbt01", "cbt02"])]),
            sample("dna02", None, &[("grp", &["cbt03"])]),
        ];
        let (routing, index) = compile(&samples, &[], &tag_sets, "pool1", true, false).unwrap();

        assert_eq!(
            routing.route(&matched(&index, &[("grp", 0)])),
            Disposition::Assigned(&Target {
                sample: "dna01".to_string(),
                sub_sample: None
            })
        );
        assert_eq!(
            routing.route(&matched(&index, &[("grp", 2)])),
            Disposition::Assigned(&Target {
                sample: "dna02".to_string(),
                sub_sample: None
            })
        );
        // A group that did not match routes to no sample.
        assert_eq!(
            routing.route(&matched(&index, &[])),
            Disposition::Unassigned
        );
    }

    #[test]
    fn test_route_and_across_groups() {
        // `both = grp_x::X + grp_y::Y` requires a best match in both groups.
        let tag_sets = sets(&[
            ("grp_x", tagset(&[("X", "AAAAAA")])),
            ("grp_y", tagset(&[("Y", "TTTTTT")])),
        ]);
        let samples = vec![sample(
            "both",
            None,
            &[("grp_x", &["X"]), ("grp_y", &["Y"])],
        )];
        let (routing, index) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();

        let both = matched(&index, &[("grp_x", 0), ("grp_y", 0)]);
        assert!(matches!(routing.route(&both), Disposition::Assigned(_)));
        // Only one of the two groups matched: unassigned.
        assert_eq!(
            routing.route(&matched(&index, &[("grp_x", 0)])),
            Disposition::Unassigned
        );
        assert_eq!(
            routing.route(&matched(&index, &[("grp_y", 0)])),
            Disposition::Unassigned
        );
    }

    #[test]
    fn test_route_bare_group_matches_any_tag() {
        // A bare-group selector routes any read that matched the group,
        // regardless of which tag.
        let tag_sets = sets(&[("grp", tagset(&[("a", "AAAA"), ("b", "CCCC")]))]);
        let samples = vec![sample("s1", None, &[("grp", &[])])];
        let (routing, index) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();

        assert!(matches!(
            routing.route(&matched(&index, &[("grp", 0)])),
            Disposition::Assigned(_)
        ));
        assert!(matches!(
            routing.route(&matched(&index, &[("grp", 1)])),
            Disposition::Assigned(_)
        ));
        assert_eq!(
            routing.route(&matched(&index, &[])),
            Disposition::Unassigned
        );
    }

    #[test]
    fn test_positional_hits_absent_and_match_only() {
        // The positional regression: an absent (None) slot must NOT satisfy a
        // bare-group selector, and a `match=`-style hit (matched but span=None)
        // must still route. Two groups so an intervening None slot is exercised
        // (gx absent while gy matched).
        let tag_sets = sets(&[
            ("gx", tagset(&[("x", "AAAA")])),
            ("gy", tagset(&[("y", "CCCC")])),
        ]);
        let samples = vec![
            sample("sx", None, &[("gx", &[])]),
            sample("sy", None, &[("gy", &["y"])]),
        ];
        let (routing, index) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();

        // gx's slot is None (it did not match): the bare-group selector for sx
        // is NOT satisfied, and gy's `match=`-style hit (span=None) still
        // routes to sy.
        assert_eq!(
            routing.route(&matched_no_span(&index, "gy", 0)),
            Disposition::Assigned(&Target {
                sample: "sy".to_string(),
                sub_sample: None
            }),
            "a span=None (match=) hit still routes; the absent gx does not claim it"
        );
        // Neither group matched: nothing claims the read.
        assert_eq!(
            routing.route(&matched(&index, &[])),
            Disposition::Unassigned
        );
    }

    #[test]
    fn test_unresolved_selector_token_is_error() {
        let tag_sets = sets(&[("grp", tagset(&[("cbt01", "AAAA")]))]);
        let samples = vec![sample("s1", None, &[("grp", &["nope"])])];
        let err = compile(&samples, &[], &tag_sets, "p", true, false).unwrap_err();
        assert!(
            err.to_string().contains("not a tag id or sequence"),
            "{err}"
        );
    }

    #[test]
    fn test_selector_token_resolves_by_sequence() {
        // A sheet may list pool members by sequence; that resolves against the
        // group tag set.
        let tag_sets = sets(&[("grp", tagset(&[("cbt01", "AACTGT")]))]);
        let samples = vec![sample("s1", None, &[("grp", &["AACTGT"])])];
        let (routing, index) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();
        assert!(matches!(
            routing.route(&matched(&index, &[("grp", 0)])),
            Disposition::Assigned(_)
        ));
    }

    #[test]
    fn test_collision_one_tag_two_samples() {
        // cbt02 routed to both dna01 and dna02 is a fail-fast collision naming
        // the tag and samples.
        let tag_sets = sets(&[(
            "grp_cbt",
            tagset(&[("c1", "AAAA"), ("c2", "CCCC"), ("c3", "GGGG")]),
        )]);
        let samples = vec![
            sample("dna01", None, &[("grp_cbt", &["c1", "c2"])]),
            sample("dna02", None, &[("grp_cbt", &["c2", "c3"])]),
        ];
        let err = compile(&samples, &[], &tag_sets, "p", true, false).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("c2"), "names the colliding tag: {msg}");
        assert!(
            msg.contains("dna01") && msg.contains("dna02"),
            "names both samples: {msg}"
        );
    }

    #[test]
    fn test_no_collision_when_same_target_or_pools_tag() {
        // The same target claiming a tag across two sheet rows (an OR pool) is
        // not a collision.
        let tag_sets = sets(&[("grp", tagset(&[("a", "AAAA"), ("b", "CCCC")]))]);
        let samples = vec![
            sample("s1", Some("lib1"), &[("grp", &["a"])]),
            sample("s1", Some("lib1"), &[("grp", &["b"])]),
        ];
        let (routing, index) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();
        // One distinct target, reachable by either tag.
        assert_eq!(routing.targets().len(), 1);
        assert!(matches!(
            routing.route(&matched(&index, &[("grp", 0)])),
            Disposition::Assigned(_)
        ));
        assert!(matches!(
            routing.route(&matched(&index, &[("grp", 1)])),
            Disposition::Assigned(_)
        ));
    }

    #[test]
    fn test_and_tuple_collision_names_combination() {
        // Two samples both claim the AND-tuple grp_x::X+grp_y::Y.
        let tag_sets = sets(&[
            ("grp_x", tagset(&[("X", "AAAA")])),
            ("grp_y", tagset(&[("Y", "TTTT")])),
        ]);
        let samples = vec![
            sample("s1", None, &[("grp_x", &["X"]), ("grp_y", &["Y"])]),
            sample("s2", None, &[("grp_y", &["Y"]), ("grp_x", &["X"])]),
        ];
        let err = compile(&samples, &[], &tag_sets, "p", true, false).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("combination"), "{msg}");
        assert!(msg.contains("s1") && msg.contains("s2"), "{msg}");
    }

    #[test]
    fn test_pool_sub_sample_resolves_to_pool_id() {
        let tag_sets = sets(&[("grp", tagset(&[("a", "AAAA")]))]);
        let samples = vec![Sample {
            sample: "dna01".to_string(),
            sub_sample: Some(SubSample::Pool),
            selector: Selector {
                terms: vec![GroupSelector {
                    group: "grp".to_string(),
                    members: vec!["a".to_string()],
                }],
            },
        }];
        let (routing, _index) = compile(&samples, &[], &tag_sets, "lib01", true, false).unwrap();
        assert_eq!(
            routing.targets()[0],
            Target {
                sample: "dna01".to_string(),
                sub_sample: Some("lib01".to_string())
            }
        );
    }

    #[test]
    fn test_remove_group_tag_takes_precedence() {
        // `--remove grp::c` removes reads matching that tag, before any sample
        // assignment.
        let tag_sets = sets(&[("grp", tagset(&[("a", "AAAA"), ("c", "GGGG")]))]);
        let samples = vec![sample("s1", None, &[("grp", &["a"])])];
        let removes = vec![RemoveRule {
            selector: RemoveSelector::GroupTag {
                group: "grp".to_string(),
                id: "c".to_string(),
            },
            pattern: None,
        }];
        let (routing, index) = compile(&samples, &removes, &tag_sets, "p", true, false).unwrap();
        assert!(matches!(
            routing.route(&matched(&index, &[("grp", 1)])),
            Disposition::Removed(_)
        ));
        assert!(matches!(
            routing.route(&matched(&index, &[("grp", 0)])),
            Disposition::Assigned(_)
        ));
    }

    #[test]
    fn test_remove_bare_group_drops_any_match() {
        let tag_sets = sets(&[("grp_y", tagset(&[("y1", "TTTT"), ("y2", "CCCC")]))]);
        let removes = vec![RemoveRule {
            selector: RemoveSelector::Group("grp_y".to_string()),
            pattern: None,
        }];
        // Pass-through run (no samples): non-removed reads pass through.
        let (routing, index) = compile(&[], &removes, &tag_sets, "p", false, false).unwrap();
        assert!(matches!(
            routing.route(&matched(&index, &[("grp_y", 0)])),
            Disposition::Removed(_)
        ));
        assert!(matches!(
            routing.route(&matched(&index, &[("grp_y", 1)])),
            Disposition::Removed(_)
        ));
        assert_eq!(
            routing.route(&matched(&index, &[])),
            Disposition::PassThrough
        );
    }

    #[test]
    fn test_remove_target_echoes_selector_label() {
        // The compiled RemoveTarget carries the user's selector verbatim, for
        // the --qc-tag removed slug: `group::id` for a specific tag, the bare
        // group name for a whole-group rule.
        let tag_sets = sets(&[("grp", tagset(&[("a", "AAAA"), ("c", "GGGG")]))]);

        let removes = vec![RemoveRule {
            selector: RemoveSelector::GroupTag {
                group: "grp".to_string(),
                id: "c".to_string(),
            },
            pattern: None,
        }];
        let (routing, index) = compile(&[], &removes, &tag_sets, "p", false, false).unwrap();
        match routing.route(&matched(&index, &[("grp", 1)])) {
            Disposition::Removed(target) => assert_eq!(target.selector, "grp::c"),
            other => panic!("expected Removed, got {other:?}"),
        }

        let removes = vec![RemoveRule {
            selector: RemoveSelector::Group("grp".to_string()),
            pattern: None,
        }];
        let (routing, index) = compile(&[], &removes, &tag_sets, "p", false, false).unwrap();
        match routing.route(&matched(&index, &[("grp", 0)])) {
            Disposition::Removed(target) => assert_eq!(target.selector, "grp"),
            other => panic!("expected Removed, got {other:?}"),
        }
    }

    #[test]
    fn test_pass_through_without_samples() {
        // No sample flags: every read passes through to --out.
        let (routing, index) = compile(&[], &[], &HashMap::new(), "p", false, false).unwrap();
        assert_eq!(
            routing.route(&matched(&index, &[])),
            Disposition::PassThrough
        );
        assert_eq!(
            routing.route(&matched(&index, &[("grp", 0)])),
            Disposition::PassThrough
        );
        assert!(routing.targets().is_empty());
    }

    #[test]
    fn test_remove_unresolved_tag_is_error() {
        let tag_sets = sets(&[("grp", tagset(&[("a", "AAAA")]))]);
        let removes = vec![RemoveRule {
            selector: RemoveSelector::GroupTag {
                group: "grp".to_string(),
                id: "missing".to_string(),
            },
            pattern: None,
        }];
        let err = compile(&[], &removes, &tag_sets, "p", false, false).unwrap_err();
        assert!(
            err.to_string().contains("not a tag id or sequence"),
            "{err}"
        );
    }

    #[test]
    fn test_coverage_fails_on_unclaimed_tag() {
        // grp has 3 tags; the samples claim only 2, so the check fails and
        // names the third.
        let tag_sets = sets(&[(
            "grp",
            tagset(&[("a", "AAAA"), ("b", "CCCC"), ("c", "GGGG")]),
        )]);
        let samples = vec![
            sample("s1", None, &[("grp", &["a"])]),
            sample("s2", None, &[("grp", &["b"])]),
        ];
        let err = compile(&samples, &[], &tag_sets, "p", true, true).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("require-samples-explain-all-tags"), "{msg}");
        assert!(msg.contains("grp::c"), "names the unclaimed tag: {msg}");
    }

    #[test]
    fn test_coverage_passes_when_all_claimed() {
        let tag_sets = sets(&[("grp", tagset(&[("a", "AAAA"), ("b", "CCCC")]))]);
        let samples = vec![
            sample("s1", None, &[("grp", &["a"])]),
            sample("s2", None, &[("grp", &["b"])]),
        ];
        assert!(compile(&samples, &[], &tag_sets, "p", true, true).is_ok());
    }

    #[test]
    fn test_coverage_bare_group_covers_all() {
        // A bare-group selector claims every tag, so full coverage holds.
        let tag_sets = sets(&[("grp", tagset(&[("a", "AAAA"), ("b", "CCCC")]))]);
        let samples = vec![sample("s1", None, &[("grp", &[])])];
        assert!(compile(&samples, &[], &tag_sets, "p", true, true).is_ok());
    }

    #[test]
    fn test_targets_dedup_and_order() {
        let tag_sets = sets(&[(
            "grp",
            tagset(&[("a", "AAAA"), ("b", "CCCC"), ("c", "GGGG")]),
        )]);
        let samples = vec![
            sample("dna01", Some("sub1"), &[("grp", &["a"])]),
            sample("dna01", Some("sub1"), &[("grp", &["b"])]),
            sample("dna01", Some("sub2"), &[("grp", &["c"])]),
        ];
        let (routing, _index) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();
        let targets = routing.targets();
        assert_eq!(targets.len(), 2);
        assert_eq!(targets[0].sub_sample.as_deref(), Some("sub1"));
        assert_eq!(targets[1].sub_sample.as_deref(), Some("sub2"));
    }

    #[test]
    fn test_expand_from_group_one_sample_per_tag() {
        let set = TagSet {
            entries: vec![
                TagEntry {
                    id: "s1".to_string(),
                    seq: "AAAA".to_string(),
                    sub_sample: Some("lib1".to_string()),
                },
                TagEntry {
                    id: "s2".to_string(),
                    seq: "CCCC".to_string(),
                    sub_sample: None,
                },
            ],
        };
        let samples = expand_from_group("sample_bc", &set);
        assert_eq!(samples.len(), 2);
        assert_eq!(samples[0].sample, "s1");
        assert_eq!(
            samples[0].sub_sample,
            Some(SubSample::Literal("lib1".to_string()))
        );
        assert_eq!(samples[1].sub_sample, None);
        // The selector routes the tag to its own sample.
        assert_eq!(samples[0].selector.terms[0].group, "sample_bc");
        assert_eq!(samples[0].selector.terms[0].members, vec!["s1".to_string()]);
    }

    #[test]
    fn test_parse_sample_sheet_basic_and_or_pool() {
        let groups: HashSet<&str> = ["grp_cbt"].into_iter().collect();
        let sheet = "sample\tsub_sample\tgrp_cbt\n\
                     dna01\tsub1\tcbt01\n\
                     dna01\tsub1\tcbt02\n\
                     dna01\tsub2\tcbt03\n\
                     dna02\tsub1\tcbt07\n";
        let samples = parse_sample_sheet(sheet, &groups).unwrap();
        assert_eq!(samples.len(), 4);
        // dna01/sub1 appears twice (an OR pool of cbt01 and cbt02).
        assert_eq!(samples[0].sample, "dna01");
        assert_eq!(
            samples[0].sub_sample,
            Some(SubSample::Literal("sub1".to_string()))
        );
        assert_eq!(
            samples[0].selector.terms[0].members,
            vec!["cbt01".to_string()]
        );
        assert_eq!(
            samples[1].selector.terms[0].members,
            vec!["cbt02".to_string()]
        );
    }

    #[test]
    fn test_parse_sample_sheet_and_across_group_columns() {
        // Two group columns on one row are AND-joined.
        let groups: HashSet<&str> = ["grp_x", "grp_y"].into_iter().collect();
        let sheet = "sample\tgrp_x\tgrp_y\ns1\tX\tY\n";
        let samples = parse_sample_sheet(sheet, &groups).unwrap();
        assert_eq!(samples[0].selector.terms.len(), 2);
        assert_eq!(samples[0].sub_sample, None);
    }

    #[test]
    fn test_parse_sample_sheet_unknown_column_is_error() {
        let groups: HashSet<&str> = ["grp_x"].into_iter().collect();
        let sheet = "sample\tgrp_x\tnonsense\ns1\tX\tZ\n";
        let err = parse_sample_sheet(sheet, &groups).unwrap_err();
        assert!(err.to_string().contains("nonsense"), "{err}");
        assert!(err.to_string().contains("declared --group"), "{err}");
    }

    #[test]
    fn test_parse_sample_sheet_missing_sample_column_is_error() {
        let groups: HashSet<&str> = ["grp_x"].into_iter().collect();
        let err = parse_sample_sheet("id\tgrp_x\ns1\tX\n", &groups).unwrap_err();
        assert!(err.to_string().contains("no `sample` column"), "{err}");
    }

    #[test]
    fn test_parse_sample_sheet_skips_comments_and_blanks() {
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let sheet = "# header comment\n\nsample\tgrp\n# row comment\ns1\ta\n\n";
        let samples = parse_sample_sheet(sheet, &groups).unwrap();
        assert_eq!(samples.len(), 1);
        assert_eq!(samples[0].sample, "s1");
    }

    #[test]
    fn test_parse_sample_sheet_pool_sub_sample_cell() {
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let sheet = "sample\tsub_sample\tgrp\ndna01\t%pool\ta\n";
        let samples = parse_sample_sheet(sheet, &groups).unwrap();
        assert_eq!(samples[0].sub_sample, Some(SubSample::Pool));
    }

    #[test]
    fn test_sheet_round_trips_through_routing() {
        // A parsed sheet feeds compile_routing exactly like inline --sample.
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let sheet = "sample\tgrp\ndna01\tcbt01\ndna02\tcbt02\n";
        let samples = parse_sample_sheet(sheet, &groups).unwrap();
        let tag_sets = sets(&[("grp", tagset(&[("cbt01", "AAAA"), ("cbt02", "CCCC")]))]);
        let (routing, index) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();
        assert_eq!(
            routing.route(&matched(&index, &[("grp", 0)])),
            Disposition::Assigned(&Target {
                sample: "dna01".to_string(),
                sub_sample: None
            })
        );
    }

    #[test]
    fn test_load_sample_sheet_round_trip() {
        use std::io::Write;
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let mut file = tempfile::NamedTempFile::new().unwrap();
        write!(file, "sample\tgrp\ns1\tcbt01\n").unwrap();
        file.flush().unwrap();
        let samples = load_sample_sheet(file.path(), &groups).unwrap();
        assert_eq!(samples[0].sample, "s1");
    }

    #[test]
    fn test_remove_pattern_is_carried_through() {
        // A `--remove SEL=PATTERN` rule keeps its output pattern for the
        // writer.
        let tag_sets = sets(&[("grp", tagset(&[("a", "AAAA")]))]);
        let pattern = OutputPattern {
            segments: vec![PatternSegment::Placeholder(Placeholder::Source)],
        };
        let removes = vec![RemoveRule {
            selector: RemoveSelector::Group("grp".to_string()),
            pattern: Some(pattern.clone()),
        }];
        let (routing, index) = compile(&[], &removes, &tag_sets, "p", false, false).unwrap();
        match routing.route(&matched(&index, &[("grp", 0)])) {
            Disposition::Removed(target) => assert_eq!(target.pattern, Some(pattern)),
            other => panic!("expected Removed, got {other:?}"),
        }
    }

    #[test]
    fn test_parse_sample_sheet_duplicate_group_column_is_error() {
        // Two columns headed by the same group AND-join it with itself (a dead
        // route); reject it.
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let err = parse_sample_sheet("sample\tgrp\tgrp\ns1\ta\tb\n", &groups).unwrap_err();
        assert!(err.to_string().contains("more than one column"), "{err}");
    }

    #[test]
    fn test_bare_group_missing_tag_set_is_error() {
        // A bare-group selector over a group whose tag set was not loaded fails
        // fast (like the member path), rather than silently claiming nothing
        // and disabling its collision check.
        let samples = vec![sample("s1", None, &[("grp", &[])])];
        let err = compile(&samples, &[], &HashMap::new(), "p", true, false).unwrap_err();
        assert!(err.to_string().contains("no loaded tag set"), "{err}");
    }
}
