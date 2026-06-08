//! The optional per-record demux-provenance slug (`--qc-tag`): a compact UTF-8
//! JSON object describing how one record was demultiplexed, written to a SAM
//! tag (default `ZS`). It is off by default and is debug/QC output, not a
//! stable API. The same top-level shape is used for assigned, pass-through,
//! unassigned, and removed records (an `outcome` discriminator), so one `jq`
//! filter works across every output bin.
//!
//! Schema (v1): `{"v":1,"outcome":..,..}`. Assigned carries `sample` and
//! `sub_sample` (the value or JSON `null`); pass-through carries neither. All
//! outcomes carry a `groups` array, one entry per matched group:
//! `{"g":<group>,"tag":<canonical
//! bases>,"loc":"<file>:<start>:<end>","obs":<observed
//! bases>,"sub":<substitutions>,"ind":<indels>[,"str":"-"]}`. `loc` is the
//! matched record span (0-based, file-relative) and `obs` the observed bases
//! there; both are absent for a `match=` group (a synthetic joined value, not a
//! record span). `str` is present only for a reverse-complement match.
//! Unassigned also carries a `reason`: `find_constraint` adds the offending
//! `group`, the failed `constraint`, and the observed `found` vs the `limit`;
//! `missing_stream` adds the `stream`. Removed carries `rule` (the `--remove`
//! selector verbatim, e.g. `grp::C` or a bare `grp`), the resolved `group`, and
//! the triggering `tag` (canonical bases; omitted for a bare-group rule). JSON
//! is built directly (no serde dependency); all string values are escaped.

use crate::extract::GroupHits;
use crate::matcher::{CompiledGroup, FindFailure};

/// Why a record was not assigned to a sample, for the unassigned slug.
pub enum Unassigned {
    /// A group's find constraint (`min`/`maxFindsPerGroup`/`PerTag`) was not
    /// met; carries the group's canonical index (resolved to its name at write
    /// time, so this stays `Send` across the matcher boundary) and the specific
    /// [`FindFailure`] (which bound, the observed count, the limit).
    FindConstraint {
        /// Canonical index of the offending group (resolved to its name at
        /// write time).
        group: usize,
        /// Which bound failed, the observed count, and the limit.
        failure: FindFailure,
    },
    /// The matched groups were claimed by no `--sample` selector.
    NoSample,
    /// A `--template` stream's anchoring group did not match, so the body could
    /// not be assembled.
    MissingStream(String),
}

/// Append `value` to `out` as a quoted, escaped JSON string.
fn push_json_string(out: &mut String, value: &str) {
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

/// Append a `"key":"value"` JSON member (escaped) to `out`.
fn push_member_str(out: &mut String, key: &str, value: &str) {
    push_json_string(out, key);
    out.push(':');
    push_json_string(out, value);
}

/// The `groups` array body: one object per matched group, with its canonical
/// tag, the observed record bases at the match span, the edit breakdown, and a
/// `str` marker for a reverse-complement match.
fn push_groups(out: &mut String, groups: &[CompiledGroup], hits: &GroupHits, segments: &[&[u8]]) {
    out.push_str(",\"groups\":[");
    let mut first = true;
    for (idx, hit) in hits.iter().enumerate() {
        let Some(hit) = hit else { continue };
        let Some(group) = groups.get(idx) else {
            continue;
        };
        if !first {
            out.push(',');
        }
        first = false;
        out.push('{');
        push_member_str(out, "g", &group.name);
        out.push(',');
        let canonical = group
            .tags
            .get(hit.tag_idx)
            .map(|t| String::from_utf8_lossy(t).into_owned())
            .unwrap_or_default();
        push_member_str(out, "tag", &canonical);
        // The matched read span and the observed bases there. `loc`
        // (`file:start:end`) is per-read information that `obs` alone does not
        // convey for variable-length, sliding, or `next`-relative groups (where
        // the offset/extent varies); both are absent for a `match=` group,
        // which matched a synthetic joined value rather than a read span.
        if let Some(span) = &hit.span {
            out.push_str(",\"loc\":\"");
            out.push_str(&span.file.to_string());
            out.push(':');
            out.push_str(&span.start.to_string());
            out.push(':');
            out.push_str(&span.end.to_string());
            out.push('"');
            if let Some(seg) = segments.get(span.file) {
                if let Some(obs) = seg.get(span.start..span.end) {
                    out.push(',');
                    push_member_str(out, "obs", &String::from_utf8_lossy(obs));
                }
            }
        }
        out.push_str(",\"sub\":");
        out.push_str(&hit.subs.to_string());
        out.push_str(",\"ind\":");
        out.push_str(&hit.indels.to_string());
        if hit.revcomp {
            out.push_str(",\"str\":\"-\"");
        }
        out.push('}');
    }
    out.push(']');
}

/// The slug for a record routed to a sample (`outcome:"assigned"`) or passed
/// through with no sample fan-out (`outcome:"pass_through"`, `sample` omitted).
pub fn routed_slug(
    groups: &[CompiledGroup],
    hits: &GroupHits,
    segments: &[&[u8]],
    sample: Option<&str>,
    sub_sample: Option<&str>,
) -> String {
    let mut out = String::from("{\"v\":1,");
    match sample {
        Some(sample) => {
            push_member_str(&mut out, "outcome", "assigned");
            out.push(',');
            push_member_str(&mut out, "sample", sample);
            // `sub_sample` completes the routing destination
            // (`sample.sub_sample` = the `@RG` ID); emitted as the value or
            // JSON `null` when the sample has no sub_sample.
            out.push_str(",\"sub_sample\":");
            match sub_sample {
                Some(sub_sample) => push_json_string(&mut out, sub_sample),
                None => out.push_str("null"),
            }
        }
        None => push_member_str(&mut out, "outcome", "pass_through"),
    }
    push_groups(&mut out, groups, hits, segments);
    out.push('}');
    out
}

/// The slug for a record that matched no sample (`outcome:"unassigned"`), with
/// the reason and, where known, the offending group or stream.
pub fn unassigned_slug(
    groups: &[CompiledGroup],
    hits: &GroupHits,
    segments: &[&[u8]],
    reason: &Unassigned,
) -> String {
    let mut out = String::from("{\"v\":1,");
    push_member_str(&mut out, "outcome", "unassigned");
    out.push(',');
    match reason {
        Unassigned::FindConstraint { group, failure } => {
            push_member_str(&mut out, "reason", "find_constraint");
            if let Some(group) = groups.get(*group) {
                out.push(',');
                push_member_str(&mut out, "group", &group.name);
            }
            // The specific bound that failed, the observed count, and the
            // limit, so an over-match (found 2, max 1) is distinguishable from
            // a no-match (found 0, min 1).
            out.push(',');
            push_member_str(&mut out, "constraint", failure.kind.as_str());
            out.push_str(",\"found\":");
            out.push_str(&failure.observed.to_string());
            out.push_str(",\"limit\":");
            out.push_str(&failure.bound.to_string());
        }
        Unassigned::NoSample => push_member_str(&mut out, "reason", "no_sample"),
        Unassigned::MissingStream(stream) => {
            push_member_str(&mut out, "reason", "missing_stream");
            out.push(',');
            push_member_str(&mut out, "stream", stream);
        }
    }
    // The groups that DID match before the read failed (empty for a read whose
    // first required group found nothing): for a no_sample failure this shows
    // the barcode combination that no selector claimed, and for a
    // find_constraint it shows the upstream groups that matched.
    push_groups(&mut out, groups, hits, segments);
    out.push('}');
    out
}

/// The slug for a record dropped by a `--remove` rule (`outcome:"removed"`):
/// `rule` echoes the user's selector verbatim (`group` or `group::id`), `group`
/// is the resolved group whose match fired the rule, and `tag` is the
/// triggering tag's canonical bases (omitted for a bare-group rule, which
/// removes on any tag). The trailing `groups` array is the standard
/// matched-group breakdown.
pub fn removed_slug(
    groups: &[CompiledGroup],
    hits: &GroupHits,
    segments: &[&[u8]],
    selector: &str,
    group: usize,
    tag: Option<usize>,
) -> String {
    let mut out = String::from("{\"v\":1,");
    push_member_str(&mut out, "outcome", "removed");
    out.push(',');
    push_member_str(&mut out, "rule", selector);
    if let Some(group) = groups.get(group) {
        out.push(',');
        push_member_str(&mut out, "group", &group.name);
        // The triggering tag's canonical bases for a `group::id` rule; omitted
        // for a bare-group rule, matching the `tag` semantics inside the
        // `groups` array (bases, not the id).
        if let Some(bases) = tag.and_then(|idx| group.tags.get(idx)) {
            out.push(',');
            push_member_str(&mut out, "tag", &String::from_utf8_lossy(bases));
        }
    }
    push_groups(&mut out, groups, hits, segments);
    out.push('}');
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::extract::{GroupHit, MatchSpan};

    fn group(name: &str, tags: &[&str]) -> CompiledGroup {
        CompiledGroup::new(name, tags.iter().map(|t| t.as_bytes().to_vec()).collect())
    }

    #[test]
    fn test_routed_slug_assigned_with_observed_and_correction() {
        let groups = vec![group("grp", &["GAAGAG"])];
        let hits: Vec<Option<GroupHit>> = vec![Some(GroupHit {
            tag_idx: 0,
            span: Some(MatchSpan {
                file: 0,
                start: 0,
                end: 6,
            }),
            subs: 1,
            indels: 0,
            revcomp: false,
        })];
        let segments: [&[u8]; 1] = [b"GAAGGGTTTT"]; // observed GAAGGG, corrected to GAAGAG
        let slug = routed_slug(&groups, &hits, &segments, Some("dna01"), Some("sub01"));
        assert_eq!(
            slug,
            r#"{"v":1,"outcome":"assigned","sample":"dna01","sub_sample":"sub01","groups":[{"g":"grp","tag":"GAAGAG","loc":"0:0:6","obs":"GAAGGG","sub":1,"ind":0}]}"#
        );
    }

    #[test]
    fn test_routed_slug_pass_through_omits_sample() {
        // No sample (a pure pass-through run): outcome is pass_through and
        // `sample` is absent.
        let groups = vec![group("grp", &["ACGT"])];
        let hits: Vec<Option<GroupHit>> = vec![Some(GroupHit {
            tag_idx: 0,
            span: Some(MatchSpan {
                file: 0,
                start: 0,
                end: 4,
            }),
            subs: 0,
            indels: 0,
            revcomp: false,
        })];
        let segments: [&[u8]; 1] = [b"ACGTACGT"];
        let slug = routed_slug(&groups, &hits, &segments, None, None);
        assert_eq!(
            slug,
            r#"{"v":1,"outcome":"pass_through","groups":[{"g":"grp","tag":"ACGT","loc":"0:0:4","obs":"ACGT","sub":0,"ind":0}]}"#
        );
    }

    #[test]
    fn test_routed_slug_match_group_has_no_observed_span_and_revcomp_marker() {
        // A `match=` group hit carries no span (it matched a synthetic joined
        // value, not a read window), so `obs` is omitted; a reverse-complement
        // match adds the `str:"-"` marker.
        let groups = vec![group("idx", &["AACC"])];
        let hits: Vec<Option<GroupHit>> = vec![Some(GroupHit {
            tag_idx: 0,
            span: None,
            subs: 0,
            indels: 1,
            revcomp: true,
        })];
        let slug = routed_slug(&groups, &hits, &[], Some("s1"), None);
        assert_eq!(
            slug,
            r#"{"v":1,"outcome":"assigned","sample":"s1","sub_sample":null,"groups":[{"g":"idx","tag":"AACC","sub":0,"ind":1,"str":"-"}]}"#
        );
    }

    #[test]
    fn test_unassigned_slug_reasons() {
        use crate::matcher::FindKind;
        let groups = vec![group("grp_cb3", &["ACGT"])];
        // No groups matched before the failure, so `groups` is empty in each
        // case.
        let no_hits: [Option<GroupHit>; 0] = [];
        // find_constraint names the offending group, the failed bound, and
        // found-vs-limit.
        assert_eq!(
            unassigned_slug(
                &groups,
                &no_hits,
                &[],
                &Unassigned::FindConstraint {
                    group: 0,
                    failure: FindFailure {
                        kind: FindKind::MinPerGroup,
                        observed: 0,
                        bound: 1
                    },
                },
            ),
            r#"{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"grp_cb3","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}"#
        );
        assert_eq!(
            unassigned_slug(&groups, &no_hits, &[], &Unassigned::NoSample),
            r#"{"v":1,"outcome":"unassigned","reason":"no_sample","groups":[]}"#
        );
        assert_eq!(
            unassigned_slug(
                &groups,
                &no_hits,
                &[],
                &Unassigned::MissingStream("bc".to_string())
            ),
            r#"{"v":1,"outcome":"unassigned","reason":"missing_stream","stream":"bc","groups":[]}"#
        );
    }

    #[test]
    fn test_removed_slug_group_tag_carries_rule_group_and_tag() {
        // A `group::id` remove rule: the slug echoes the selector in `rule`,
        // names the resolved `group`, and carries the triggering tag's
        // canonical bases in `tag`, plus the groups array.
        let groups = vec![group("grp", &["CCCC"])];
        let hits: Vec<Option<GroupHit>> = vec![Some(GroupHit {
            tag_idx: 0,
            span: Some(MatchSpan {
                file: 0,
                start: 0,
                end: 4,
            }),
            subs: 0,
            indels: 0,
            revcomp: false,
        })];
        let segments: [&[u8]; 1] = [b"CCCCAAAA"];
        let slug = removed_slug(&groups, &hits, &segments, "grp::C", 0, Some(0));
        assert_eq!(
            slug,
            r#"{"v":1,"outcome":"removed","rule":"grp::C","group":"grp","tag":"CCCC","groups":[{"g":"grp","tag":"CCCC","loc":"0:0:4","obs":"CCCC","sub":0,"ind":0}]}"#
        );
    }

    #[test]
    fn test_removed_slug_bare_group_omits_tag() {
        // A bare-group rule (removes on any tag) has no specific tag, so the
        // `tag` field is omitted entirely (not `null`); `rule` and `group` both
        // read the bare group name.
        let groups = vec![group("grp", &["CCCC"])];
        let hits: Vec<Option<GroupHit>> = vec![Some(GroupHit {
            tag_idx: 0,
            span: Some(MatchSpan {
                file: 0,
                start: 0,
                end: 4,
            }),
            subs: 0,
            indels: 0,
            revcomp: false,
        })];
        let segments: [&[u8]; 1] = [b"CCCCAAAA"];
        let slug = removed_slug(&groups, &hits, &segments, "grp", 0, None);
        assert_eq!(
            slug,
            r#"{"v":1,"outcome":"removed","rule":"grp","group":"grp","groups":[{"g":"grp","tag":"CCCC","loc":"0:0:4","obs":"CCCC","sub":0,"ind":0}]}"#
        );
    }
}
