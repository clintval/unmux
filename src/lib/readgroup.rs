//! Resolve an `@RG` group source against the input header: turn the header's
//! read groups into a synthetic tag set, the 1:1 fan-out samples, and an
//! RG-id to tag-index map the matcher uses to route each record by its `RG:Z`.
use std::collections::HashMap;

use anyhow::{bail, Result};

use crate::grammar::{GroupSelector, ReadGroupKey, Sample, Selector, SubSample};
use crate::input::{InputReadGroup, ReadGroupInfo};
use crate::tags::{TagEntry, TagSet};

/// The synthetic tag set, samples, and RG-id map for one `@RG` group.
pub struct ResolvedReadGroup {
    /// One entry per distinct target key; `entry.id` is a unique routing token.
    pub tag_set: TagSet,
    /// One sample per entry (`sample`/`sub_sample` are the real SM/LB).
    pub samples: Vec<Sample>,
    /// Each RG id mapped to its entry index (its `tag_idx`).
    pub rg_to_idx: HashMap<Vec<u8>, usize>,
}

/// The (sample, sub_sample) a read group maps to under `key`. A missing subfield
/// falls back to the RG id so no record is silently dropped.
fn key_of(rg: &InputReadGroup, key: ReadGroupKey) -> (String, Option<String>) {
    let id = String::from_utf8_lossy(&rg.id).into_owned();
    let sm = rg.sample.clone().unwrap_or_else(|| id.clone());
    let lb = rg.library.clone().unwrap_or_else(|| id.clone());
    match key {
        ReadGroupKey::Id => (id, None),
        ReadGroupKey::Sm => (sm, None),
        ReadGroupKey::Lb => (lb, None),
        ReadGroupKey::SmLb => (sm, Some(lb)),
    }
}

/// The unique routing token / synthetic tag id for a target key.
fn label(tier1: &str, tier2: Option<&str>) -> String {
    match tier2 {
        Some(t2) => format!("{tier1}.{t2}"),
        None => tier1.to_string(),
    }
}

/// Resolve the read groups declared in the input header against `key`, one
/// tag-set entry / sample per distinct target, and a map from each RG id to
/// its entry index.
pub fn resolve(group: &str, info: &ReadGroupInfo, key: ReadGroupKey) -> Result<ResolvedReadGroup> {
    if info.read_groups.is_empty() {
        bail!("group `{group}` uses `@RG` but the input header declares no @RG lines (the `@RG` source needs a SAM/BAM/CRAM header)");
    }
    let mut tag_set = TagSet::default();
    let mut samples = Vec::new();
    let mut rg_to_idx = HashMap::new();
    let mut key_to_idx: HashMap<(String, Option<String>), usize> = HashMap::new();

    for rg in &info.read_groups {
        let (tier1, tier2) = key_of(rg, key);
        let idx = *key_to_idx
            .entry((tier1.clone(), tier2.clone()))
            .or_insert_with(|| {
                let idx = tag_set.entries.len();
                // Fold `idx` into the token so it stays unique even when two
                // distinct (SM, LB) pairs format to the same `label()` string
                // (e.g. `("sA", "l1.2")` and `("sA.l1", "2")` both read `sA.l1.2`).
                let token = format!("{}#{idx}", label(&tier1, tier2.as_deref()));
                tag_set.entries.push(TagEntry {
                    id: token.clone(),
                    seq: token.clone(),
                    sub_sample: tier2.clone(),
                });
                samples.push(Sample {
                    sample: tier1.clone(),
                    sub_sample: tier2.clone().map(SubSample::Literal),
                    selector: Selector {
                        terms: vec![GroupSelector {
                            group: group.to_string(),
                            members: vec![token],
                        }],
                    },
                });
                idx
            });
        rg_to_idx.insert(rg.id.clone(), idx);
    }
    Ok(ResolvedReadGroup {
        tag_set,
        samples,
        rg_to_idx,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn info(rgs: &[(&str, Option<&str>, Option<&str>)]) -> ReadGroupInfo {
        ReadGroupInfo {
            read_groups: rgs
                .iter()
                .map(|(id, sm, lb)| InputReadGroup {
                    id: id.as_bytes().to_vec(),
                    sample: sm.map(String::from),
                    library: lb.map(String::from),
                    map: Default::default(),
                })
                .collect(),
            has_sq: false,
        }
    }

    #[test]
    fn resolve_by_id_one_target_per_read_group() {
        let r = resolve(
            "rg",
            &info(&[("rg1", Some("sA"), None), ("rg2", Some("sB"), None)]),
            ReadGroupKey::Id,
        )
        .unwrap();
        assert_eq!(
            r.tag_set
                .entries
                .iter()
                .map(|e| e.id.clone())
                .collect::<Vec<_>>(),
            vec!["rg1#0", "rg2#1"]
        );
        assert_eq!(
            r.samples
                .iter()
                .map(|s| s.sample.clone())
                .collect::<Vec<_>>(),
            vec!["rg1", "rg2"]
        );
        assert_eq!(r.rg_to_idx[b"rg1".as_slice()], 0);
        assert_eq!(r.rg_to_idx[b"rg2".as_slice()], 1);
    }

    #[test]
    fn resolve_by_sm_merges_read_groups_sharing_a_sample() {
        let r = resolve(
            "sm",
            &info(&[
                ("rg1", Some("sA"), None),
                ("rg2", Some("sB"), None),
                ("rg3", Some("sA"), None),
            ]),
            ReadGroupKey::Sm,
        )
        .unwrap();
        assert_eq!(
            r.samples
                .iter()
                .map(|s| s.sample.clone())
                .collect::<Vec<_>>(),
            vec!["sA", "sB"]
        );
        assert_eq!(r.rg_to_idx[b"rg1".as_slice()], 0);
        assert_eq!(r.rg_to_idx[b"rg3".as_slice()], 0); // merged with rg1
        assert_eq!(r.rg_to_idx[b"rg2".as_slice()], 1);
    }

    #[test]
    fn resolve_two_tier_keeps_sample_and_subsample() {
        let r = resolve(
            "b",
            &info(&[
                ("rg1", Some("sA"), Some("l1")),
                ("rg2", Some("sA"), Some("l2")),
            ]),
            ReadGroupKey::SmLb,
        )
        .unwrap();
        assert_eq!(r.samples.len(), 2);
        assert_eq!(r.samples[0].sample, "sA");
        assert!(
            matches!(r.samples[0].sub_sample, Some(crate::grammar::SubSample::Literal(ref s)) if s == "l1")
        );
        assert_ne!(r.tag_set.entries[0].id, r.tag_set.entries[1].id); // unique routing ids
    }

    #[test]
    fn resolve_missing_subfield_falls_back_to_rg_id() {
        let r = resolve("sm", &info(&[("rg1", None, None)]), ReadGroupKey::Sm).unwrap();
        assert_eq!(r.samples[0].sample, "rg1"); // no SM: use the RG id
    }

    #[test]
    fn resolve_empty_header_errors() {
        assert!(resolve("rg", &info(&[]), ReadGroupKey::Id).is_err());
    }

    #[test]
    fn resolve_two_tier_tokens_unique_even_with_dotted_values() {
        let r = resolve(
            "b",
            &info(&[
                ("rg1", Some("sA"), Some("l1.2")),
                ("rg2", Some("sA.l1"), Some("2")),
            ]),
            ReadGroupKey::SmLb,
        )
        .unwrap();
        assert_eq!(r.tag_set.entries.len(), 2);
        assert_ne!(r.tag_set.entries[0].id, r.tag_set.entries[1].id);
        assert_ne!(
            r.rg_to_idx[b"rg1".as_slice()],
            r.rg_to_idx[b"rg2".as_slice()]
        );
        // each Sample still carries the real SM as its sample name
        assert_eq!(r.samples[0].sample, "sA");
        assert_eq!(r.samples[1].sample, "sA.l1");
    }
}
