//! IUPAC degeneracy validation for tag groups: a combinatorial-explosion guard
//! (warn plus a hard cap) and a mode-aware collision check for two tags whose
//! degenerate sets overlap within `dist`. Matching itself uses sassy's native
//! `Iupac` profile with no expansion; this module conceptually expands the
//! degenerate codes only to validate the tag set up front, so a barcode that
//! explodes combinatorially or two barcodes that cannot be told apart are
//! caught before a run rather than silently mis-assigning records. It also
//! provides the base-level IUPAC primitives ([`complement`],
//! [`reverse_complement`]) shared by extraction, input parsing, and tag
//! correction.

use anyhow::{bail, Result};

use crate::grammar::{Dist, MatchMode};

/// Warn when one tag expands to at least this many concrete sequences; a
/// barcode this degenerate (six fully-degenerate positions is already 4096) is
/// almost always a mistake.
const TAG_DEGEN_WARN: u64 = 4096;
/// Reject a single tag that expands beyond this: it is not a barcode.
const TAG_DEGEN_CAP: u64 = 1 << 20;
/// Warn / reject on a group's total expanded size (the sum over its tags).
const GROUP_DEGEN_WARN: u64 = 100_000;
const GROUP_DEGEN_CAP: u64 = 10_000_000;

/// The 4-bit ACGT membership mask for an IUPAC code (A=1, C=2, G=4, T=8). A
/// byte that is not a recognized code maps to 0, so it matches nothing and
/// never overlaps another tag.
pub(crate) fn mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101, // A G
        b'Y' => 0b1010, // C T
        b'S' => 0b0110, // C G
        b'W' => 0b1001, // A T
        b'K' => 0b1100, // G T
        b'M' => 0b0011, // A C
        b'B' => 0b1110, // C G T
        b'D' => 0b1101, // A G T
        b'H' => 0b1011, // A C T
        b'V' => 0b0111, // A C G
        b'N' => 0b1111,
        _ => 0,
    }
}

/// The number of concrete ACGT sequences a tag expands to: the product of its
/// per-position option counts, saturating so a long degenerate tag cannot
/// overflow.
pub fn degeneracy(seq: &[u8]) -> u64 {
    seq.iter().fold(1u64, |acc, &b| {
        acc.saturating_mul(mask(b).count_ones().max(1) as u64)
    })
}

/// Validate one group's tags. Returns soft warnings (the caller logs them); an
/// explosion past the hard cap or a disqualifying collision is a fail-fast
/// error.
///
/// Collision is mode-aware. Two equal-length tags whose IUPAC sets intersect at
/// every position are indistinguishable (a record matches both with zero
/// errors) and always fail. Otherwise, with `D` the number of positions where
/// the sets are disjoint, a record can stay within `k = dist.mismatch` of both
/// iff `D <= 2*k` (at a disjoint position a record base satisfies at most one
/// of the two, costing the other a mismatch): under `mode=all` that is a silent
/// mis-assignment and fails; under `mode=nearest` the tie is resolved or
/// dropped by `delta`, so it only warns. Cross-length (indel-mediated) overlaps
/// are not checked.
pub fn validate_group(
    name: &str,
    seqs: &[Vec<u8>],
    dist: &Dist,
    mode: MatchMode,
) -> Result<Vec<String>> {
    let mut warnings = Vec::new();

    let mut group_total: u64 = 0;
    for seq in seqs {
        let d = degeneracy(seq);
        group_total = group_total.saturating_add(d);
        let shown = String::from_utf8_lossy(seq);
        // An all-`N` tag matches any sequence of its length: it is a
        // fixed-length positional span, not a barcode. Reject it (regardless of
        // length, so a short `NNNN` that slips under the degeneracy cap is
        // still caught) and point at location-based slicing.
        if !seq.is_empty() && seq.iter().all(|&b| mask(b) == 0b1111) {
            bail!(
                "group `{name}`: tag `{shown}` is all-N (it matches any {} bases), which is a \
                 fixed-length span, not a barcode; slice it with a location window \
                 (loc=file:start:stop) or an anchored extract (@group+offset:len) instead",
                seq.len()
            );
        }
        if d > TAG_DEGEN_CAP {
            bail!("group `{name}`: tag `{shown}` expands to {d} sequences, over the {TAG_DEGEN_CAP} cap (too degenerate to be a barcode)");
        }
        if d > TAG_DEGEN_WARN {
            warnings.push(format!(
                "group `{name}`: tag `{shown}` is highly degenerate ({d} expansions)"
            ));
        }
    }
    if group_total > GROUP_DEGEN_CAP {
        bail!("group `{name}`: tags expand to {group_total} sequences total, over the {GROUP_DEGEN_CAP} cap");
    }
    if group_total > GROUP_DEGEN_WARN {
        warnings.push(format!(
            "group `{name}`: tags expand to {group_total} sequences total (consider tightening the degenerate positions)"
        ));
    }

    let k = dist.mismatch;
    for (i, a) in seqs.iter().enumerate() {
        for b in &seqs[i + 1..] {
            if a.len() != b.len() {
                continue;
            }
            let disjoint = a
                .iter()
                .zip(b)
                .filter(|(&x, &y)| mask(x) & mask(y) == 0)
                .count();
            let (a_s, b_s) = (String::from_utf8_lossy(a), String::from_utf8_lossy(b));
            if disjoint == 0 {
                bail!("group `{name}`: tags `{a_s}` and `{b_s}` are indistinguishable (a read matches both with zero errors); remove the duplicate or make them distinct");
            }
            if disjoint <= 2 * k {
                let overlap = format!(
                    "group `{name}`: tags `{a_s}` and `{b_s}` can both match one read within dist={k} (they differ at only {disjoint} position(s))"
                );
                match mode {
                    MatchMode::All => bail!("{overlap}; with mode=all this silently mis-assigns - separate the tags, lower dist, or use mode=nearest with a delta"),
                    MatchMode::Nearest => warnings.push(format!(
                        "{overlap}; mode=nearest leaves ambiguous reads unassigned unless delta separates them"
                    )),
                }
            }
        }
    }

    Ok(warnings)
}

/// IUPAC-aware complement of a single base (case-preserving). Unrecognized
/// bytes pass through.
pub(crate) fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'u' => b'a',
        b'r' => b'y',
        b'y' => b'r',
        b's' => b's',
        b'w' => b'w',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        b'n' => b'n',
        other => other,
    }
}

/// Reverse-complement a sequence in place: reverse the bytes and
/// IUPAC-[`complement`] each.
pub(crate) fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    for base in seq.iter_mut() {
        *base = complement(*base);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dist(mismatch: usize) -> Dist {
        Dist {
            mismatch,
            indel: 0,
            total: mismatch,
        }
    }

    fn seqs(items: &[&str]) -> Vec<Vec<u8>> {
        items.iter().map(|s| s.as_bytes().to_vec()).collect()
    }

    #[test]
    fn test_degeneracy_products() {
        assert_eq!(degeneracy(b"ACGT"), 1);
        assert_eq!(degeneracy(b"NNNNNNNN"), 65536); // 4^8
        assert_eq!(degeneracy(b"ACGTR"), 2); // one two-way code
        assert_eq!(degeneracy(b"ACGTB"), 3); // one three-way code
    }

    #[test]
    fn test_distinct_barcodes_pass() {
        let w = validate_group(
            "g",
            &seqs(&["AAAA", "CCCC", "GGGG"]),
            &dist(1),
            MatchMode::All,
        )
        .unwrap();
        assert!(w.is_empty(), "well-separated barcodes are clean: {w:?}");
    }

    #[test]
    fn test_indistinguishable_pair_fails_regardless_of_mode() {
        // AAAW (W = A|T) covers AAAA at every position -> a read AAAA matches
        // both at dist 0.
        for mode in [MatchMode::All, MatchMode::Nearest] {
            let err = validate_group("g", &seqs(&["AAAA", "AAAW"]), &dist(0), mode)
                .unwrap_err()
                .to_string();
            assert!(err.contains("indistinguishable"), "{err}");
        }
    }

    #[test]
    fn test_within_dist_overlap_is_mode_aware() {
        // AAAA vs AATT differ at 2 positions; at dist=1 (k=1) D=2 <= 2*k ->
        // overlap.
        let tags = seqs(&["AAAA", "AATT"]);
        // mode=all: hard error.
        let err = validate_group("g", &tags, &dist(1), MatchMode::All)
            .unwrap_err()
            .to_string();
        assert!(err.contains("mode=all"), "{err}");
        // mode=nearest: warning, not error.
        let w = validate_group("g", &tags, &dist(1), MatchMode::Nearest).unwrap();
        assert_eq!(w.len(), 1, "one collision warning: {w:?}");
        assert!(w[0].contains("within dist=1"));
    }

    #[test]
    fn test_well_separated_pair_ok_at_dist() {
        // AAAA vs TTTT differ at all 4 positions; at dist=1, D=4 > 2 -> no
        // collision.
        assert!(
            validate_group("g", &seqs(&["AAAA", "TTTT"]), &dist(1), MatchMode::All)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn test_explosion_hard_cap_errors() {
        // 10 fully-degenerate positions plus a 3-way `B` = 4^10 * 3 ~= 3.1M,
        // over the 1<<20 cap. The `B` keeps it from being all-N (which has its
        // own, earlier error), so this exercises the cap.
        let err = validate_group("g", &seqs(&["NNNNNNNNNNB"]), &dist(0), MatchMode::All)
            .unwrap_err()
            .to_string();
        assert!(err.contains("cap"), "{err}");
    }

    #[test]
    fn test_all_n_tag_is_rejected_with_slicing_hint() {
        // A short all-N tag (NNNN = 256, under the degeneracy cap) is still
        // rejected: it matches any 4 bases, so it is a positional span, not a
        // barcode. This is the gap the cap alone misses.
        let err = validate_group("bc", &seqs(&["NNNN"]), &dist(0), MatchMode::All)
            .unwrap_err()
            .to_string();
        assert!(err.contains("all-N"), "{err}");
        assert!(
            err.contains("loc=") || err.contains("@group"),
            "suggests slicing: {err}"
        );
    }

    #[test]
    fn test_lowercase_all_n_tag_is_rejected() {
        // The mask is case-insensitive, so lowercase nnnn is caught too.
        assert!(validate_group("bc", &seqs(&["nnnn"]), &dist(0), MatchMode::All).is_err());
    }

    #[test]
    fn test_partially_degenerate_tag_is_not_all_n() {
        // A tag with any non-N base is not all-N; NNNANNN (4^6 = 4096, at the
        // warn threshold, not over it) is accepted, so the all-N guard does not
        // fire on it.
        assert!(validate_group("bc", &seqs(&["NNNANNN"]), &dist(0), MatchMode::All).is_ok());
    }
}
