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

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};

use crate::extract::GroupHits;
use crate::grammar::{
    escape_column, GroupSelector, OutputPattern, RemoveRule, RemoveSelector, Sample, Selector,
    SubSample,
};
use crate::tags::{find_header, is_column_value, known_columns, TagSet};

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
    /// Routed to a fan-out target by this route, which for a sample sheet is one
    /// of its rows.
    Assigned(&'a Route),
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
#[derive(Debug, Clone, PartialEq, Eq)]
struct Term {
    group: usize,
    members: Option<Vec<usize>>,
}

/// One fan-out route: a target plus the AND-joined terms a record must satisfy
/// to reach it, and the value of each `%sample::COLUMN` stream for the records
/// it routes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Route {
    target: Target,
    target_index: usize,
    terms: Vec<Term>,
    columns: Vec<Option<String>>,
}

impl Route {
    /// The fan-out target this route assigns records to.
    pub fn target(&self) -> &Target {
        &self.target
    }

    /// The position of this route's target in [`Routing::targets`].
    pub fn target_index(&self) -> usize {
        self.target_index
    }

    /// The value of the `slot`th `%sample::COLUMN` stream, in extract order, or
    /// `None` when the cell is blank.
    pub(crate) fn column(&self, slot: usize) -> Option<&str> {
        self.columns.get(slot)?.as_deref()
    }

    /// The canonical indices of the groups this route reads, sorted.
    fn groups(&self) -> Vec<usize> {
        let mut groups: Vec<usize> = self.terms.iter().map(|term| term.group).collect();
        groups.sort_unstable();
        groups.dedup();
        groups
    }

    /// The tags this route accepts in `group`, or `None` when it accepts any.
    fn members(&self, group: usize) -> Option<&[usize]> {
        let term = self.terms.iter().find(|term| term.group == group)?;
        term.members.as_deref()
    }
}

/// A route that loses a record it matches to an earlier route of another target.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Shadowed {
    /// The route that takes such a record.
    taker: usize,
    /// The route that loses it.
    loser: usize,
    /// The record's best tag as `(group, tag)` in each group either route reads.
    tags: Vec<(usize, usize)>,
}

/// The most routes [`shadowed_warnings`] lists before counting the rest.
const LISTED_SHADOWED: usize = 10;

/// The sorted groups a set of routes fills, and those routes' indices.
type Shape = (Vec<usize>, Vec<usize>);

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
        if let Some(remove) = self.removal(hits) {
            return Disposition::Removed(remove);
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
                return Disposition::Assigned(route);
            }
        }
        Disposition::Unassigned
    }

    /// The distinct fan-out targets in first-seen order (one `@RG` / output bin
    /// each).
    pub fn targets(&self) -> &[Target] {
        &self.targets
    }

    /// The `--remove` rules, in the order they are checked.
    pub(crate) fn removes(&self) -> &[RemoveTarget] {
        &self.removes
    }

    /// The first `--remove` rule that `hits` satisfy, if any.
    pub(crate) fn removal(&self, hits: &GroupHits) -> Option<&RemoveTarget> {
        self.removes.iter().find(|remove| {
            term_satisfied(
                hits,
                remove.group,
                remove.tag.as_ref().map(std::slice::from_ref),
            )
        })
    }

    /// Whether `hits` satisfy a `--remove` rule, or one on group `next` or later still may.
    pub(crate) fn may_remove(&self, hits: &GroupHits, next: usize) -> bool {
        self.removes.iter().any(|remove| remove.group >= next) || self.removal(hits).is_some()
    }

    /// Each route's target, in route order.
    pub(crate) fn route_targets(&self) -> impl ExactSizeIterator<Item = &Target> {
        self.routes.iter().map(|route| &route.target)
    }

    /// Whether a route reads the match of the group at canonical index `group`.
    pub(crate) fn routes_on(&self, group: usize) -> bool {
        self.routes
            .iter()
            .any(|route| route.terms.iter().any(|term| term.group == group))
    }

    /// Store the `%sample::COLUMN` stream values on the routes, one list per
    /// stream in extract order, each holding one value per route.
    pub(crate) fn set_columns(&mut self, columns: Vec<Vec<Option<String>>>) -> Result<()> {
        if let Some(values) = columns.iter().find(|v| v.len() != self.routes.len()) {
            bail!(
                "a sample column has {} values for {} routes",
                values.len(),
                self.routes.len()
            );
        }
        for route in &mut self.routes {
            route.columns.clear();
        }
        for values in columns {
            for (route, value) in self.routes.iter_mut().zip(values) {
                route.columns.push(value);
            }
        }
        Ok(())
    }

    /// Fail when two sample-sheet rows can route one record to the same target
    /// but hold different values in a `%sample::COLUMN` stream, whose column
    /// headers are `names` in extract order. A sheet row selects one tag per
    /// filled group cell, so two rows can route one record when they select the
    /// same tags in the groups both fill. Rows whose overlap an earlier row of
    /// another sample always claims first are still reported.
    pub(crate) fn check_column_conflicts(&self, names: &[&str]) -> Result<()> {
        let mut shapes: HashMap<(&Target, Vec<usize>), Vec<usize>> = HashMap::new();
        for (index, route) in self.routes.iter().enumerate() {
            let mut groups: Vec<usize> = route.terms.iter().map(|term| term.group).collect();
            groups.sort_unstable();
            shapes
                .entry((&route.target, groups))
                .or_default()
                .push(index);
        }
        let mut by_target: HashMap<&Target, Vec<Shape>> = HashMap::new();
        for ((target, groups), rows) in shapes {
            by_target.entry(target).or_default().push((groups, rows));
        }
        let values = |index: usize| &self.routes[index].columns;
        let mut conflict: Option<(usize, usize)> = None;
        for shapes in by_target.values_mut() {
            shapes.sort_unstable_by_key(|(_, rows)| rows[0]);
            for (i, (first, first_rows)) in shapes.iter().enumerate() {
                for (second, second_rows) in &shapes[i..] {
                    let shared: Vec<usize> = first
                        .iter()
                        .filter(|group| second.contains(group))
                        .copied()
                        .collect();
                    let key = |index: usize| -> Vec<Option<&[usize]>> {
                        let terms = &self.routes[index].terms;
                        shared
                            .iter()
                            .map(|group| {
                                terms
                                    .iter()
                                    .find(|term| term.group == *group)
                                    .and_then(|term| term.members.as_deref())
                            })
                            .collect()
                    };
                    let mut seen: HashMap<Vec<Option<&[usize]>>, Vec<usize>> = HashMap::new();
                    for &index in second_rows {
                        let rows = seen.entry(key(index)).or_default();
                        if rows.len() < 2 && rows.iter().all(|&row| values(row) != values(index)) {
                            rows.push(index);
                        }
                    }
                    for &index in first_rows {
                        let Some(rows) = seen.get(&key(index)) else {
                            continue;
                        };
                        if let Some(&row) = rows.iter().find(|&&row| values(row) != values(index)) {
                            let pair = (row.min(index), row.max(index));
                            conflict = Some(conflict.map_or(pair, |found| found.min(pair)));
                        }
                    }
                }
            }
        }
        let Some((first, second)) = conflict else {
            return Ok(());
        };
        let (a, b) = (values(first), values(second));
        let slot = (0..a.len()).find(|&slot| a[slot] != b[slot]).unwrap_or(0);
        let quote = |value: &Option<String>| match value {
            Some(value) => format!("{value:?}"),
            None => "a blank cell".to_string(),
        };
        bail!(
            "sample sheet rows {} and {} can both route a record to `{}` but disagree on {:?} \
             ({} vs {})",
            first + 1,
            second + 1,
            self.routes[first].target.label(),
            names.get(slot).copied().unwrap_or_default(),
            quote(&a[slot]),
            quote(&b[slot])
        )
    }

    /// Each route that loses a record it matches to an earlier route of another
    /// target, with the route that takes the record. Routes over the same
    /// groups never share a record, as compiling refuses a tag combination that
    /// two targets claim, so each route is compared with earlier routes over
    /// other groups, from the first, until one takes a record they share.
    fn shadowed(&self) -> Vec<Shadowed> {
        let mut shapes: Vec<Shape> = Vec::new();
        for (index, route) in self.routes.iter().enumerate() {
            let groups = route.groups();
            match shapes.iter_mut().find(|(known, _)| *known == groups) {
                Some((_, routes)) => routes.push(index),
                None => shapes.push((groups, vec![index])),
            }
        }
        if shapes.len() < 2 {
            return Vec::new();
        }
        let mut found: Vec<Option<(usize, Shadowed)>> = vec![None; self.routes.len()];
        for (earlier_groups, earlier) in &shapes {
            for (later_groups, later) in &shapes {
                if earlier_groups == later_groups {
                    continue;
                }
                // Index the earlier routes by their tags in a group that both shapes read.
                let pivot = earlier_groups
                    .iter()
                    .copied()
                    .find(|g| later_groups.contains(g));
                let mut by_tag: HashMap<usize, Vec<usize>> = HashMap::new();
                let mut any_tag = Vec::new();
                if let Some(pivot) = pivot {
                    for &index in earlier {
                        match self.routes[index].members(pivot) {
                            Some(tags) => {
                                for &tag in tags {
                                    by_tag.entry(tag).or_default().push(index);
                                }
                            }
                            None => any_tag.push(index),
                        }
                    }
                }
                for &loser in later {
                    let tags = pivot.and_then(|pivot| self.routes[loser].members(pivot));
                    let candidates: Cow<[usize]> = match tags {
                        Some(tags) => {
                            let mut listed: Vec<usize> = tags
                                .iter()
                                .filter_map(|tag| by_tag.get(tag))
                                .flatten()
                                .chain(&any_tag)
                                .copied()
                                .collect();
                            listed.sort_unstable();
                            listed.dedup();
                            Cow::Owned(listed)
                        }
                        None => Cow::Borrowed(earlier),
                    };
                    let target = &self.routes[loser].target;
                    let lost = candidates
                        .iter()
                        .copied()
                        .take_while(|&candidate| candidate < loser)
                        .filter(|&candidate| self.routes[candidate].target != *target)
                        .find_map(|candidate| {
                            let tags = shared_record(&self.routes[candidate], &self.routes[loser])?;
                            let taker = self.first_route(&tags)?;
                            (self.routes[taker].target != *target)
                                .then_some((candidate, Shadowed { taker, loser, tags }))
                        });
                    if let Some((candidate, shadowed)) = lost {
                        if found[loser]
                            .as_ref()
                            .is_none_or(|(known, _)| candidate < *known)
                        {
                            found[loser] = Some((candidate, shadowed));
                        }
                    }
                }
            }
        }
        found
            .into_iter()
            .flatten()
            .map(|(_, shadowed)| shadowed)
            .collect()
    }

    /// The first route that a record with the best tags `tags` satisfies.
    fn first_route(&self, tags: &[(usize, usize)]) -> Option<usize> {
        self.routes.iter().position(|route| {
            route.terms.iter().all(|term| {
                tags.iter().any(|&(group, tag)| {
                    group == term.group && term.members.as_ref().is_none_or(|m| m.contains(&tag))
                })
            })
        })
    }
}

/// The best tags of a record that satisfies both routes, one per group either
/// reads, or `None` when no record does. Where both accept any tag of a group,
/// the record has its first.
fn shared_record(a: &Route, b: &Route) -> Option<Vec<(usize, usize)>> {
    let mut groups: Vec<usize> = a.terms.iter().chain(&b.terms).map(|t| t.group).collect();
    groups.sort_unstable();
    groups.dedup();
    groups
        .into_iter()
        .map(|group| {
            let mut lists = a
                .terms
                .iter()
                .chain(&b.terms)
                .filter(|term| term.group == group)
                .filter_map(|term| term.members.as_deref());
            let Some(first) = lists.next() else {
                return Some((group, 0));
            };
            let rest: Vec<&[usize]> = lists.collect();
            let tag = first
                .iter()
                .find(|tag| rest.iter().all(|list| list.contains(tag)))?;
            Some((group, *tag))
        })
        .collect()
}

/// Warnings for the routes that lose a record they match to an earlier route of
/// another target, for the rows of a sample sheet when `sheet` is set and for
/// `--sample` targets otherwise. The first ten are listed and the rest counted.
pub(crate) fn shadowed_warnings(
    routing: &Routing,
    tag_sets: &HashMap<String, TagSet>,
    group_index: &HashMap<String, usize>,
    sheet: bool,
) -> Vec<String> {
    let shadowed = routing.shadowed();
    let names: HashMap<usize, &str> = group_index
        .iter()
        .map(|(name, &index)| (index, name.as_str()))
        .collect();
    let mut warnings: Vec<String> = shadowed
        .iter()
        .take(LISTED_SHADOWED)
        .map(|shadowed| {
            let mut unit: Vec<(String, usize)> = shadowed
                .tags
                .iter()
                .map(|&(group, tag)| (names[&group].to_string(), tag))
                .collect();
            unit.sort();
            let read = render_unit(&unit, tag_sets);
            let taker = routing.routes[shadowed.taker].target.label();
            let loser = routing.routes[shadowed.loser].target.label();
            match sheet {
                true => format!(
                    "sample sheet rows {} (`{taker}`) and {} (`{loser}`) can both route one read, such as one with {read}; rows are tried in order, so such reads go to `{taker}`",
                    shadowed.taker + 1,
                    shadowed.loser + 1
                ),
                false => format!(
                    "samples `{taker}` and `{loser}` can both route one read, such as one with {read}; --sample targets are tried in order, so such reads go to `{taker}`"
                ),
            }
        })
        .collect();
    let more = shadowed.len().saturating_sub(LISTED_SHADOWED);
    if more > 0 {
        let plural = if more == 1 { "" } else { "s" };
        warnings.push(match sheet {
            true => format!(
                "{more} more sample sheet row{plural} can lose reads to an earlier row of another sample"
            ),
            false => {
                format!("{more} more sample{plural} can lose reads to an earlier --sample target")
            }
        });
    }
    warnings
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
        let target_index = match targets.iter().position(|known| *known == target) {
            Some(index) => index,
            None => {
                targets.push(target.clone());
                targets.len() - 1
            }
        };
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
        routes.push(Route {
            target,
            target_index,
            terms,
            columns: Vec::new(),
        });
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
pub(crate) fn resolve_target(sample: &Sample, pool: &str) -> Target {
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

/// What a `%sample::COLUMN` extract reads from each sample-sheet row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum SheetColumn {
    /// The routed target's sample ID.
    Sample,
    /// The routed target's sub-sample ID, when it has one.
    SubSample,
    /// Each row's trimmed cell, `None` when blank.
    Cells(Vec<Option<String>>),
}

/// The sample-sheet columns that `%sample::COLUMN` extracts read, keyed by the
/// column as the extract spells it.
pub(crate) type SheetColumns = HashMap<String, SheetColumn>;

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
///
/// Columns named in `ignored_columns` (`--sample-sheet-ignore`) are skipped
/// unread. Each must appear in the header and may not be a reserved column or a
/// declared group, so a typo still fails rather than silently skipping nothing.
pub fn parse_sample_sheet(
    text: &str,
    declared_groups: &HashSet<&str>,
    ignored_columns: &HashSet<&str>,
) -> Result<Vec<Sample>> {
    parse_sample_sheet_columns(text, declared_groups, ignored_columns, &[]).map(|(rows, _)| rows)
}

/// Parse a `--sample-sheet` TSV as [`parse_sample_sheet`] does, along with the
/// columns that `reads` pairs of `(extract, column)` read. A read column may
/// not be a group column, and its cells hold printable ASCII and spaces.
pub(crate) fn parse_sample_sheet_columns(
    text: &str,
    declared_groups: &HashSet<&str>,
    ignored_columns: &HashSet<&str>,
    reads: &[(&str, &str)],
) -> Result<(Vec<Sample>, SheetColumns)> {
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

    let mut ignored: Vec<&str> = ignored_columns.iter().copied().collect();
    ignored.sort_unstable();
    for column in ignored {
        if column == "sample" || column == "sub_sample" {
            bail!("--sample-sheet-ignore cannot skip the reserved column `{column}`");
        }
        if declared_groups.contains(column) {
            bail!(
                "--sample-sheet-ignore names `{column}`, a declared --group; a group column cannot be ignored"
            );
        }
        if !columns.contains(&column) {
            bail!(
                "--sample-sheet-ignore names `{column}`, which is not a column in the sample sheet header"
            );
        }
    }

    let mut sheet_columns = SheetColumns::new();
    let mut read_cols: Vec<(String, usize)> = Vec::new();
    for &(extract, column) in reads {
        if sheet_columns.contains_key(column) {
            continue;
        }
        let found = find_header(&columns, column).map_err(|_| {
            anyhow!(
                "extract `{extract}` reads sample-sheet column {column:?}, but more than one \
                 column matches it regardless of case"
            )
        })?;
        let read = match found {
            Some(i) if i == sample_col => SheetColumn::Sample,
            Some(i) if Some(i) == sub_sample_col => SheetColumn::SubSample,
            None if column.eq_ignore_ascii_case("sub_sample") => SheetColumn::SubSample,
            None => bail!(
                "extract `{extract}` reads sample-sheet column {column:?}; the columns are: {}",
                known_columns(&columns)
            ),
            Some(i) if declared_groups.contains(columns[i]) => bail!(
                "extract `{extract}` reads sample-sheet column {column:?}, a group column; use \
                 `@{}::id` for the matched tag",
                columns[i]
            ),
            Some(i) => {
                if ignored_columns.contains(columns[i]) {
                    log::warn!(
                        "--sample-sheet-ignore skips {:?}, which `%sample::{}` reads; the column \
                         is read",
                        columns[i],
                        escape_column(column)
                    );
                }
                read_cols.push((column.to_string(), i));
                SheetColumn::Cells(Vec::new())
            }
        };
        sheet_columns.insert(column.to_string(), read);
    }

    // Every remaining non-empty column must be a declared group or be read by
    // an extract; a group column's cells select within it.
    let mut group_cols = Vec::new();
    let mut seen_groups = HashSet::new();
    for (i, column) in columns.iter().enumerate() {
        if i == sample_col
            || Some(i) == sub_sample_col
            || column.is_empty()
            || ignored_columns.contains(column)
            || read_cols.iter().any(|(_, read)| *read == i)
        {
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
        for (column, i) in &read_cols {
            let value = cells
                .get(*i)
                .map(|cell| cell.trim_matches(' '))
                .filter(|cell| !cell.is_empty());
            if let Some(value) = value.filter(|value| !is_column_value(value)) {
                bail!(
                    "sample sheet row {}: value {value:?} in column {:?} may hold only printable \
                     ASCII and spaces",
                    line_no + 1,
                    columns[*i]
                );
            }
            if let Some(SheetColumn::Cells(values)) = sheet_columns.get_mut(column) {
                values.push(value.map(str::to_string));
            }
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
    Ok((samples, sheet_columns))
}

/// Load and parse a `--sample-sheet` from a file.
pub fn load_sample_sheet(
    path: &Path,
    declared_groups: &HashSet<&str>,
    ignored_columns: &HashSet<&str>,
) -> Result<Vec<Sample>> {
    load_sample_sheet_columns(path, declared_groups, ignored_columns, &[]).map(|(rows, _)| rows)
}

/// Load and parse a `--sample-sheet` from a file, along with the columns that
/// `reads` pairs of `(extract, column)` read.
pub(crate) fn load_sample_sheet_columns(
    path: &Path,
    declared_groups: &HashSet<&str>,
    ignored_columns: &HashSet<&str>,
    reads: &[(&str, &str)],
) -> Result<(Vec<Sample>, SheetColumns)> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read sample sheet: {}", path.display()))?;
    parse_sample_sheet_columns(&text, declared_groups, ignored_columns, reads)
        .with_context(|| format!("failed to parse sample sheet: {}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::extract::{GroupHit, MatchSpan};
    use crate::grammar::{PatternSegment, Placeholder};
    use crate::tags::TagEntry;

    /// The target a record was assigned to, if any.
    fn assigned(disposition: Disposition<'_>) -> Option<&Target> {
        match disposition {
            Disposition::Assigned(route) => Some(route.target()),
            _ => None,
        }
    }

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
            assigned(routing.route(&matched(&index, &[("grp", 0)]))),
            Some(&Target {
                sample: "dna01".to_string(),
                sub_sample: None
            })
        );
        assert_eq!(
            assigned(routing.route(&matched(&index, &[("grp", 2)]))),
            Some(&Target {
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
            assigned(routing.route(&matched_no_span(&index, "gy", 0))),
            Some(&Target {
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
        let samples = parse_sample_sheet(sheet, &groups, &HashSet::new()).unwrap();
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
        let samples = parse_sample_sheet(sheet, &groups, &HashSet::new()).unwrap();
        assert_eq!(samples[0].selector.terms.len(), 2);
        assert_eq!(samples[0].sub_sample, None);
    }

    #[test]
    fn test_parse_sample_sheet_unknown_column_is_error() {
        let groups: HashSet<&str> = ["grp_x"].into_iter().collect();
        let sheet = "sample\tgrp_x\tnonsense\ns1\tX\tZ\n";
        let err = parse_sample_sheet(sheet, &groups, &HashSet::new()).unwrap_err();
        assert!(err.to_string().contains("nonsense"), "{err}");
        assert!(err.to_string().contains("declared --group"), "{err}");
    }

    #[test]
    fn test_parse_sample_sheet_missing_sample_column_is_error() {
        let groups: HashSet<&str> = ["grp_x"].into_iter().collect();
        let err = parse_sample_sheet("id\tgrp_x\ns1\tX\n", &groups, &HashSet::new()).unwrap_err();
        assert!(err.to_string().contains("no `sample` column"), "{err}");
    }

    #[test]
    fn test_parse_sample_sheet_skips_comments_and_blanks() {
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let sheet = "# header comment\n\nsample\tgrp\n# row comment\ns1\ta\n\n";
        let samples = parse_sample_sheet(sheet, &groups, &HashSet::new()).unwrap();
        assert_eq!(samples.len(), 1);
        assert_eq!(samples[0].sample, "s1");
    }

    #[test]
    fn test_parse_sample_sheet_pool_sub_sample_cell() {
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let sheet = "sample\tsub_sample\tgrp\ndna01\t%pool\ta\n";
        let samples = parse_sample_sheet(sheet, &groups, &HashSet::new()).unwrap();
        assert_eq!(samples[0].sub_sample, Some(SubSample::Pool));
    }

    #[test]
    fn test_sheet_round_trips_through_routing() {
        // A parsed sheet feeds compile_routing exactly like inline --sample.
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let sheet = "sample\tgrp\ndna01\tcbt01\ndna02\tcbt02\n";
        let samples = parse_sample_sheet(sheet, &groups, &HashSet::new()).unwrap();
        let tag_sets = sets(&[("grp", tagset(&[("cbt01", "AAAA"), ("cbt02", "CCCC")]))]);
        let (routing, index) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();
        assert_eq!(
            assigned(routing.route(&matched(&index, &[("grp", 0)]))),
            Some(&Target {
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
        let samples = load_sample_sheet(file.path(), &groups, &HashSet::new()).unwrap();
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
        let err = parse_sample_sheet("sample\tgrp\tgrp\ns1\ta\tb\n", &groups, &HashSet::new())
            .unwrap_err();
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

    #[test]
    fn test_parse_sample_sheet_ignored_column_is_skipped() {
        // An ignored column's cells are never read, so any value is fine.
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let ignored: HashSet<&str> = ["assay"].into_iter().collect();
        let sheet = "sample\tassay\tgrp\ns1\tscATAC\ta\ns2\t%not a tag, ::+\tb\ns3\t\tc\n";
        let samples = parse_sample_sheet(sheet, &groups, &ignored).unwrap();
        assert_eq!(samples.len(), 3);
        for sample in &samples {
            assert_eq!(sample.selector.terms.len(), 1);
            assert_eq!(sample.selector.terms[0].group, "grp");
        }
        assert_eq!(samples[1].selector.terms[0].members, vec!["b".to_string()]);
    }

    #[test]
    fn test_parse_sample_sheet_ignored_column_absent_from_header_is_error() {
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let ignored: HashSet<&str> = ["asay"].into_iter().collect();
        let sheet = "sample\tassay\tgrp\ns1\tx\ta\n";
        let err = parse_sample_sheet(sheet, &groups, &ignored).unwrap_err();
        assert!(err.to_string().contains("`asay`"), "{err}");
        assert!(err.to_string().contains("not a column"), "{err}");
    }

    #[test]
    fn test_parse_sample_sheet_ignoring_reserved_column_is_error() {
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let sheet = "sample\tsub_sample\tgrp\ns1\tlib\ta\n";
        for reserved in ["sample", "sub_sample"] {
            let ignored: HashSet<&str> = [reserved].into_iter().collect();
            let err = parse_sample_sheet(sheet, &groups, &ignored).unwrap_err();
            assert!(err.to_string().contains("reserved column"), "{err}");
            assert!(err.to_string().contains(reserved), "{err}");
        }
    }

    #[test]
    fn test_parse_sample_sheet_ignoring_declared_group_is_error() {
        let groups: HashSet<&str> = ["grp", "other"].into_iter().collect();
        let ignored: HashSet<&str> = ["other"].into_iter().collect();
        let sheet = "sample\tgrp\tother\ns1\ta\tb\n";
        let err = parse_sample_sheet(sheet, &groups, &ignored).unwrap_err();
        assert!(err.to_string().contains("declared --group"), "{err}");
        assert!(err.to_string().contains("`other`"), "{err}");
    }

    #[test]
    fn test_parse_sample_sheet_ignore_leaves_other_columns_strict() {
        // Ignoring one extra column does not relax the check on another.
        let groups: HashSet<&str> = ["grp"].into_iter().collect();
        let ignored: HashSet<&str> = ["assay"].into_iter().collect();
        let sheet = "sample\tassay\tgrp\tnotes\ns1\tx\ta\ty\n";
        let err = parse_sample_sheet(sheet, &groups, &ignored).unwrap_err();
        assert!(err.to_string().contains("`notes`"), "{err}");
        assert!(err.to_string().contains("declared --group"), "{err}");
    }

    /// Parse a sheet whose extracts read `columns`, with groups `a` and `b`.
    fn sheet_columns(
        sheet: &str,
        ignored: &[&str],
        columns: &[&str],
    ) -> Result<(Vec<Sample>, SheetColumns)> {
        let groups: HashSet<&str> = ["a", "b"].into_iter().collect();
        let ignored: HashSet<&str> = ignored.iter().copied().collect();
        let reads: Vec<(&str, &str)> = columns.iter().map(|column| ("v", *column)).collect();
        parse_sample_sheet_columns(sheet, &groups, &ignored, &reads)
    }

    fn cells(values: &[Option<&str>]) -> SheetColumn {
        SheetColumn::Cells(values.iter().map(|v| v.map(str::to_string)).collect())
    }

    #[test]
    fn test_parse_sample_sheet_columns_reads_one_trimmed_value_per_row() {
        let sheet = "sample\ta\tnucleus\n\
                     s1\tx\t  007 Ab  \n\
                     # a comment\n\
                     \n\
                     s1\ty\t\n\
                     s2\tz\tx,y;z+w\n\
                     s2\tw\n";
        let (samples, columns) = sheet_columns(sheet, &[], &["nucleus"]).unwrap();
        assert_eq!(samples.len(), 4);
        assert_eq!(
            columns["nucleus"],
            cells(&[Some("007 Ab"), None, Some("x,y;z+w"), None])
        );
        assert_eq!(
            samples,
            parse_sample_sheet(
                "sample\ta\ns1\tx\ns1\ty\ns2\tz\ns2\tw\n",
                &["a", "b"].into_iter().collect(),
                &HashSet::new()
            )
            .unwrap()
        );
    }

    #[test]
    fn test_parse_sample_sheet_columns_reserved_columns() {
        let sheet = "sample\tsub_sample\ta\ns1\t%pool\tx\n";
        let (_, columns) = sheet_columns(sheet, &[], &["sample", "SAMPLE", "sub_sample"]).unwrap();
        assert_eq!(columns["sample"], SheetColumn::Sample);
        assert_eq!(columns["SAMPLE"], SheetColumn::Sample);
        assert_eq!(columns["sub_sample"], SheetColumn::SubSample);
        let (_, columns) = sheet_columns("sample\ta\ns1\tx\n", &[], &["sub_sample"]).unwrap();
        assert_eq!(columns["sub_sample"], SheetColumn::SubSample);
    }

    #[test]
    fn test_parse_sample_sheet_columns_finds_headers() {
        let read = |header: &str, column: &str| {
            sheet_columns(
                &format!("sample\ta\t{header}\ns1\tx\tv1\tv2\n"),
                &[],
                &[column],
            )
            .map(|(_, columns)| columns[column].clone())
            .map_err(|err| err.to_string())
        };
        assert_eq!(read("nucleus", "nucleus").unwrap(), cells(&[Some("v1")]));
        assert_eq!(read("Nucleus", "nucleus").unwrap(), cells(&[Some("v1")]));
        assert_eq!(read("a b", "a b").unwrap(), cells(&[Some("v1")]));
        assert_eq!(
            read("Nucleus\tnucleus", "nucleus").unwrap_err(),
            "sample sheet column `Nucleus` is not a declared --group (every non-reserved column \
             must name a group)"
        );
        assert_eq!(
            read("nucleus", "nucleu").unwrap_err(),
            "extract `v` reads sample-sheet column \"nucleu\"; the columns are: sample, a, nucleus"
        );
        assert_eq!(
            read("Nucleus\tNUCLEUS", "nucleus").unwrap_err(),
            "extract `v` reads sample-sheet column \"nucleus\", but more than one column matches \
             it regardless of case"
        );
        assert_eq!(
            read("b", "b").unwrap_err(),
            "extract `v` reads sample-sheet column \"b\", a group column; use `@b::id` for the \
             matched tag"
        );
        let (_, columns) = sheet_columns(
            "sample\ta\tNucleus\tnucleus\ns1\tx\tv1\tv2\n",
            &[],
            &["Nucleus", "nucleus"],
        )
        .unwrap();
        assert_eq!(columns["Nucleus"], cells(&[Some("v1")]));
        assert_eq!(columns["nucleus"], cells(&[Some("v2")]));
    }

    #[test]
    fn test_parse_sample_sheet_columns_read_and_ignored_is_read() {
        let (_, columns) =
            sheet_columns("sample\ta\tnucleus\ns1\tx\tv\n", &["nucleus"], &["nucleus"]).unwrap();
        assert_eq!(columns["nucleus"], cells(&[Some("v")]));
    }

    #[test]
    fn test_parse_sample_sheet_columns_rejects_a_value_that_is_not_printable_ascii() {
        let err = sheet_columns(
            "sample\ta\tnucleus\ns1\tx\tok\ns2\ty\t\u{c1}1\n",
            &[],
            &["nucleus"],
        )
        .unwrap_err()
        .to_string();
        assert_eq!(
            err,
            "sample sheet row 2: value \"\u{c1}1\" in column \"nucleus\" may hold only printable \
             ASCII and spaces"
        );
    }

    /// Compile a two-group sheet and check `values` of column `v` for conflicts.
    fn check_rows(sheet: &str, values: &[Option<&str>]) -> std::result::Result<(), String> {
        check_columns(sheet, &[("v", values)])
    }

    /// Compile a two-group sheet and check its `(column, values)` for conflicts.
    fn check_columns(
        sheet: &str,
        columns: &[(&str, &[Option<&str>])],
    ) -> std::result::Result<(), String> {
        let groups: HashSet<&str> = ["a", "b"].into_iter().collect();
        let samples = parse_sample_sheet(sheet, &groups, &HashSet::new()).unwrap();
        let tag_sets = sets(&[
            ("a", tagset(&[("a1", "AAAA"), ("a2", "CCCC")])),
            ("b", tagset(&[("b1", "GGGG"), ("b2", "TTTT")])),
        ]);
        let (mut routing, _) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();
        let names: Vec<&str> = columns.iter().map(|(name, _)| *name).collect();
        let values = columns
            .iter()
            .map(|(_, values)| values.iter().map(|v| v.map(str::to_string)).collect())
            .collect();
        routing.set_columns(values).unwrap();
        routing
            .check_column_conflicts(&names)
            .map_err(|err| err.to_string())
    }

    #[test]
    fn test_check_column_conflicts_compares_every_distinct_value_of_an_overlap() {
        let sheet = "sample\ta\tb\ns1\ta1\t\ns1\ta1\tb1\ns1\ta1\tb2\n";
        assert_eq!(
            check_rows(sheet, &[Some("x"), Some("x"), Some("y")]).unwrap_err(),
            "sample sheet rows 1 and 3 can both route a record to `s1` but disagree on \"v\" \
             (\"x\" vs \"y\")"
        );
        let reversed = "sample\ta\tb\ns1\ta1\tb2\ns1\ta1\tb1\ns1\ta1\t\n";
        assert_eq!(
            check_rows(reversed, &[Some("y"), Some("x"), Some("x")]).unwrap_err(),
            "sample sheet rows 1 and 3 can both route a record to `s1` but disagree on \"v\" \
             (\"y\" vs \"x\")"
        );
    }

    #[test]
    fn test_check_column_conflicts_checks_every_column_at_once() {
        let sheet = "sample\ta\tb\ns1\ta1\tb1\ns1\ta1\tb1\ns2\ta2\tb2\n";
        let same: &[Option<&str>] = &[Some("x"), Some("x"), Some("z")];
        let differ: &[Option<&str>] = &[Some("p"), Some("q"), Some("z")];
        assert_eq!(check_columns(sheet, &[("u", same), ("w", same)]), Ok(()));
        let err = check_columns(sheet, &[("u", same), ("w", differ)]).unwrap_err();
        assert!(err.ends_with("disagree on \"w\" (\"p\" vs \"q\")"), "{err}");
    }

    #[test]
    fn test_set_columns_needs_one_value_per_route() {
        let groups: HashSet<&str> = ["a"].into_iter().collect();
        let samples =
            parse_sample_sheet("sample\ta\ns1\ta1\ns2\ta2\n", &groups, &HashSet::new()).unwrap();
        let tag_sets = sets(&[("a", tagset(&[("a1", "AAAA"), ("a2", "CCCC")]))]);
        let (mut routing, index) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();
        assert!(routing.set_columns(vec![vec![None]]).is_err());
        let values = vec![
            vec![Some("x".to_string()), None],
            vec![Some("s1".to_string()), Some("s2".to_string())],
        ];
        routing.set_columns(values).unwrap();
        let Disposition::Assigned(route) = routing.route(&matched(&index, &[("a", 1)])) else {
            panic!("not assigned");
        };
        assert_eq!((route.column(0), route.column(1)), (None, Some("s2")));
    }

    #[test]
    fn fuzz_parse_sample_sheet_columns() {
        let declared: HashSet<&str> = ["bc", "grp", "cbt"].into_iter().collect();
        let reads = [("v", "note"), ("s", "sample"), ("u", "sub_sample")];
        bolero::check!()
            .with_type::<String>()
            .for_each(|text: &String| {
                let _ = parse_sample_sheet_columns(text, &declared, &HashSet::new(), &reads);
            });
    }

    #[test]
    fn test_check_row_values_rejects_rows_that_can_route_one_record_to_one_target() {
        let duplicate = "sample\ta\tb\ns1\ta1\tb1\ns1\ta1\tb1\n";
        assert_eq!(
            check_rows(duplicate, &[Some("x"), Some("y")]).unwrap_err(),
            "sample sheet rows 1 and 2 can both route a record to `s1` but disagree on \"v\" \
             (\"x\" vs \"y\")"
        );
        let subset = "sample\ta\tb\ns2\ta2\tb2\ns1\ta1\tb1\ns1\ta1\t\n";
        let err = check_rows(subset, &[Some("x"), Some("x"), None]).unwrap_err();
        assert!(
            err.contains("rows 2 and 3 can both route a record to `s1`")
                && err.ends_with("(\"x\" vs a blank cell)"),
            "{err}"
        );
        let crossing = "sample\ta\tb\ns1\ta1\t\ns1\t\tb1\n";
        assert!(check_rows(crossing, &[Some("x"), Some("y")]).is_err());
        let by_sequence = "sample\ta\tb\ns1\ta1\tb1\ns1\tAAAA\tGGGG\n";
        assert!(check_rows(by_sequence, &[Some("x"), Some("y")]).is_err());
        let third = "sample\ta\tb\ns1\ta1\tb1\ns1\ta1\tb2\ns1\ta1\t\n";
        let err = check_rows(third, &[Some("x"), Some("y"), Some("x")]).unwrap_err();
        assert!(err.contains("rows 2 and 3"), "{err}");
    }

    #[test]
    fn test_check_row_values_accepts_rows_that_agree_or_cannot_overlap() {
        let sheet = "sample\ta\tb\n\
                     s1\ta1\tb1\n\
                     s1\ta1\tb1\n\
                     s1\ta1\t\n\
                     s1\ta2\tb2\n\
                     s2\ta1\tb2\n";
        let values = [Some("x"), Some("x"), Some("x"), Some("y"), Some("z")];
        assert_eq!(check_rows(sheet, &values), Ok(()));
        let blanks = "sample\ta\tb\ns1\ta1\tb1\ns1\ta1\tb1\n";
        assert_eq!(check_rows(blanks, &[None, None]), Ok(()));
    }

    #[test]
    fn test_route_reports_the_route_that_matched() {
        let groups: HashSet<&str> = ["a", "b"].into_iter().collect();
        let sheet = "sample\ta\tb\ns1\ta1\tb1\ns1\ta1\tb2\ns2\ta2\t\n";
        let samples = parse_sample_sheet(sheet, &groups, &HashSet::new()).unwrap();
        let tag_sets = sets(&[
            ("a", tagset(&[("a1", "AAAA"), ("a2", "CCCC")])),
            ("b", tagset(&[("b1", "GGGG"), ("b2", "TTTT")])),
        ]);
        let (routing, index) = compile(&samples, &[], &tag_sets, "p", true, false).unwrap();
        let route = |pairs: &[(&str, usize)]| match routing.route(&matched(&index, pairs)) {
            Disposition::Assigned(route) => Some((
                route.target().sample.clone(),
                routing
                    .routes
                    .iter()
                    .position(|r| std::ptr::eq(r, route))
                    .unwrap(),
            )),
            _ => None,
        };
        assert_eq!(route(&[("a", 0), ("b", 0)]), Some(("s1".to_string(), 0)));
        assert_eq!(route(&[("a", 0), ("b", 1)]), Some(("s1".to_string(), 1)));
        assert_eq!(route(&[("a", 1), ("b", 1)]), Some(("s2".to_string(), 2)));
        assert_eq!(route(&[("a", 1)]), Some(("s2".to_string(), 2)));
        assert_eq!(route(&[("a", 0)]), None);
        let targets: Vec<&str> = routing.route_targets().map(|t| t.sample.as_str()).collect();
        assert_eq!(targets, ["s1", "s1", "s2"]);
    }

    /// A route that loses a record, as `(taker, loser, tags)`.
    type Lost = (usize, usize, Vec<(usize, usize)>);

    /// Each route of `samples` that loses a record, over groups `a`, `b`, and
    /// `c` (indices 0, 1, 2) of two tags each.
    fn shadowed_of(samples: &[Sample]) -> Vec<Lost> {
        let tag_sets = sets(&[
            ("a", tagset(&[("a0", "AAAA"), ("a1", "CCCC")])),
            ("b", tagset(&[("b0", "GGGG"), ("b1", "TTTT")])),
            ("c", tagset(&[("c0", "ACAC"), ("c1", "GTGT")])),
        ]);
        let (routing, _) = compile(samples, &[], &tag_sets, "p", true, false).unwrap();
        let shadowed = routing.shadowed().into_iter();
        shadowed.map(|s| (s.taker, s.loser, s.tags)).collect()
    }

    #[test]
    fn test_shadowed_finds_routes_over_other_groups_that_share_a_record() {
        let crossed = [
            sample("s1", None, &[("a", &["a0"])]),
            sample("s2", None, &[("b", &["b0"])]),
        ];
        assert_eq!(shadowed_of(&crossed), [(0, 1, vec![(0, 0), (1, 0)])]);
        let covered = [
            sample("s1", None, &[("a", &["a0"])]),
            sample("s2", None, &[("a", &["a0"]), ("b", &["b1"])]),
        ];
        assert_eq!(shadowed_of(&covered), [(0, 1, vec![(0, 0), (1, 1)])]);
        let narrower_first = [
            sample("s2", None, &[("a", &["a0"]), ("b", &["b1"])]),
            sample("s1", None, &[("a", &["a0"])]),
        ];
        assert_eq!(shadowed_of(&narrower_first), [(0, 1, vec![(0, 0), (1, 1)])]);
        let bare = [
            sample("s1", None, &[("a", &[])]),
            sample("s2", None, &[("a", &["a1"]), ("b", &["b1"])]),
        ];
        assert_eq!(shadowed_of(&bare), [(0, 1, vec![(0, 1), (1, 1)])]);
        let pooled = [
            sample("s1", None, &[("a", &["a0", "a1"])]),
            sample("s2", None, &[("a", &["a1"]), ("b", &["b0"])]),
        ];
        assert_eq!(shadowed_of(&pooled), [(0, 1, vec![(0, 1), (1, 0)])]);
        let sub_samples = [
            sample("s1", Some("x"), &[("a", &["a0"])]),
            sample("s1", Some("y"), &[("b", &["b0"])]),
        ];
        assert_eq!(shadowed_of(&sub_samples), [(0, 1, vec![(0, 0), (1, 0)])]);
        let several = [
            sample("s1", None, &[("a", &["a0"])]),
            sample("s2", None, &[("b", &["b0"])]),
            sample("s3", None, &[("c", &["c0"])]),
        ];
        assert_eq!(
            shadowed_of(&several),
            [(0, 1, vec![(0, 0), (1, 0)]), (0, 2, vec![(0, 0), (2, 0)])]
        );
    }

    #[test]
    fn test_shadowed_skips_routes_that_share_no_record_or_a_target() {
        let same_groups = [
            sample("s1", None, &[("a", &["a0"]), ("b", &["b0"])]),
            sample("s2", None, &[("a", &["a1"]), ("b", &["b0"])]),
        ];
        assert!(shadowed_of(&same_groups).is_empty());
        let disjoint = [
            sample("s1", None, &[("a", &["a0"])]),
            sample("s2", None, &[("a", &["a1"]), ("b", &["b0"])]),
        ];
        assert!(shadowed_of(&disjoint).is_empty());
        let one_sample = [
            sample("s1", None, &[("a", &["a0"])]),
            sample("s1", None, &[("b", &["b0"])]),
        ];
        assert!(shadowed_of(&one_sample).is_empty());
    }

    #[test]
    fn test_shadowed_routes_the_shared_record_to_find_its_taker() {
        // Route 0 of the loser's own sample takes every record that route 2 shares with route 1.
        let own_first = [
            sample("s2", None, &[("b", &["b0"])]),
            sample("s1", None, &[("a", &["a0"])]),
            sample("s2", None, &[("a", &["a0"]), ("b", &["b0"])]),
        ];
        assert_eq!(shadowed_of(&own_first), [(0, 1, vec![(0, 0), (1, 0)])]);
        // Route 3 loses to route 1 the record it shares with route 2, not the one route 0 takes.
        let earlier_taker = [
            sample("s2", None, &[("c", &["c0"])]),
            sample("s1", None, &[("a", &[]), ("c", &[])]),
            sample("s3", None, &[("c", &["c1"])]),
            sample("s2", None, &[("a", &["a0"])]),
        ];
        assert_eq!(
            shadowed_of(&earlier_taker),
            [
                (0, 1, vec![(0, 0), (2, 0)]),
                (1, 2, vec![(0, 0), (2, 1)]),
                (1, 3, vec![(0, 0), (2, 1)]),
            ]
        );
    }
}
