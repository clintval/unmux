//! IUPAC degeneracy validation for tag groups: a combinatorial-explosion guard
//! (warn plus a hard cap) and a mode-aware check for two tags that one read can
//! match together. Matching itself uses sassy's native `Iupac` profile with no
//! expansion; this module reads the degenerate codes as base sets only to
//! validate the tag set up front, so a barcode that explodes combinatorially or
//! two barcodes that one read can match are caught before a run rather than
//! silently mis-assigning records. It also provides the base-level IUPAC
//! primitives ([`complement`], [`reverse_complement`]) shared by extraction,
//! input parsing, and tag correction.

use std::borrow::Cow;
use std::cell::Cell;
use std::cmp::Ordering;
use std::collections::HashMap;

use anyhow::{bail, Result};
use smallvec::SmallVec;

use crate::grammar::{Anchor, Dist, MatchMode};

/// Warn when one tag expands to at least this many concrete sequences; a
/// barcode this degenerate (six fully-degenerate positions is already 4096) is
/// almost always a mistake.
const TAG_DEGEN_WARN: u64 = 4096;
/// Reject a single tag that expands beyond this: it is not a barcode.
const TAG_DEGEN_CAP: u64 = 1 << 20;
/// Warn / reject on a group's total expanded size (the sum over its tags).
const GROUP_DEGEN_WARN: u64 = 100_000;
const GROUP_DEGEN_CAP: u64 = 10_000_000;

/// The colliding pairs a group warns about one by one.
const LISTED_PAIRS: usize = 10;
/// The colliding pairs a group counts; past them, only the pairs that would
/// still be refused are looked for.
const COUNTED_PAIRS: usize = 1_000;
/// What `mode=nearest` does with a read that two tags match.
const UNASSIGNED: &str =
    "mode=nearest leaves ambiguous reads unassigned unless delta separates them";
/// What an explicit `mode=all` does with a read that two tags match.
const TIE_BREAK: &str = "with mode=all such a read goes to one of the tags by tie-break";
/// The mode recommended for tags that one read can match only through an indel.
const NEAREST: &str = "set mode=nearest to leave such reads unassigned (recommended)";
/// The most cells a table of one tag's alignments within another may hold; a
/// pair that needs a larger one is taken to collide through an indel.
const TABLE_CELLS: usize = 1 << 20;
/// The most cells the tables of one group may hold in all; past them, a pair
/// that needs a table is taken to collide through an indel.
const GROUP_CELLS: usize = 1 << 27;
/// Why a pair taken to collide gives no bases that match both tags.
const TOO_COSTLY: &str = "(too costly to compare exactly under this indel budget)";

thread_local! {
    /// The cells left for the tables of the group checked on this thread, kept
    /// here so that the comparison every pair goes through takes no budget.
    static CELLS_LEFT: Cell<usize> = const { Cell::new(GROUP_CELLS) };
}
/// The most pairs compared one by one for two tags that one read matches with
/// zero errors before a refusal; plain tags of one length are looked up instead.
const EXACT_PAIRS: usize = 1 << 20;

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

/// How a group searches records for its tags, which decides whether one read
/// can match two of them.
#[derive(Debug, Clone, Copy)]
pub struct Search {
    /// The error budget each tag is matched within.
    pub dist: Dist,
    /// How matches within the budget are accepted, or `None` when the group
    /// leaves it to the default, `mode=all`.
    pub mode: Option<MatchMode>,
    /// The window edge each tag is pinned to, or `None` to let tags slide.
    pub anchor: Option<Anchor>,
    /// Whether tags also match the reverse-complement strand.
    pub both_strands: bool,
    /// The widest the search window can be, or `None` when nothing bounds it.
    pub window: Option<usize>,
}

/// Validate one group's tags. Returns soft warnings (the caller logs them); an
/// explosion past the hard cap or a disqualifying collision is a fail-fast
/// error.
///
/// Two tags collide when one read can match both at one place: some bases, no
/// longer than the search window, that the first tag matches end to end within
/// `dist` and the second matches within (at their start under `anchor=5p`, at
/// their end under `anchor=3p`), either tag taken first, and the second also
/// reverse-complemented when `bothStrands` is set. Each tag keeps its own
/// substitution, indel, and total caps, and the message gives bases checked to
/// match both, without an indel when that is enough; a pair whose table of
/// alignments would hold more than `TABLE_CELLS` cells, or more than the
/// group has left of `GROUP_CELLS`, is taken to collide through an indel
/// instead.
/// A collision with zero errors always fails, before any other refusal, and
/// points at `bothStrands` when only opposite strands match.
/// One that needs no indel fails under `mode=all`, written or not, which would
/// assign such reads by tie-break; one that needs an indel fails only when no
/// mode is written, asking for one. The rest are warnings: the first pairs are
/// listed and the rest counted up to a limit, past which only collisions that
/// would fail are looked for, and under `mode=all` a last warning recommends
/// `mode=nearest`. Truncated matches (`partial5`/`partial3`) are not considered.
pub fn validate_group(name: &str, seqs: &[Vec<u8>], search: &Search) -> Result<Vec<String>> {
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

    let shape = Shape::of(search.anchor);
    let tags: Vec<Tag> = seqs.iter().map(|seq| Tag::new(seq)).collect();
    let flipped: Vec<Tag> = match search.both_strands {
        true => seqs
            .iter()
            .map(|seq| {
                let mut seq = seq.clone();
                reverse_complement(&mut seq);
                Tag::new(&seq)
            })
            .collect(),
        false => Vec::new(),
    };
    let subs = search.dist.mismatch.min(search.dist.total);
    let substitutions = Search {
        dist: Dist {
            mismatch: subs,
            indel: 0,
            total: subs,
        },
        ..*search
    };
    let exact_only = Search {
        dist: Dist::default(),
        ..*search
    };
    let nearest = search.mode == Some(MatchMode::Nearest);
    let (still_refused, consequence) = match nearest {
        true => (&exact_only, UNASSIGNED),
        false => (&substitutions, TIE_BREAK),
    };
    let budget = dist_spec(&search.dist);
    CELLS_LEFT.set(GROUP_CELLS);
    let (mut colliding, mut compared, mut example) = (0, 0, None);
    for (i, a) in tags.iter().enumerate() {
        for (j, b) in tags.iter().enumerate().skip(i + 1) {
            let strands = search.both_strands.then(|| (&flipped[i], &flipped[j]));
            let checked = match colliding < COUNTED_PAIRS {
                true => search,
                false => still_refused,
            };
            let Some(found) = collision(a, b, strands, checked, shape) else {
                continue;
            };
            let plain = match search.dist.indel {
                0 => None,
                _ => collision(a, b, strands, &substitutions, shape),
            };
            let needs_indel = search.dist.indel > 0 && plain.is_none();
            let found = plain.unwrap_or(found);
            let (a_s, b_s) = (
                String::from_utf8_lossy(&seqs[i]),
                String::from_utf8_lossy(&seqs[j]),
            );
            let bases = found.bases.as_deref().map(String::from_utf8_lossy);
            let strand = if found.flipped {
                ", one on each strand"
            } else {
                ""
            };
            if found.exact {
                let found = collision(a, b, None, &exact_only, shape).unwrap_or(found);
                bail!("{}", indistinguishable(name, &seqs[i], &seqs[j], &found));
            }
            let (verb, witness) = match bases {
                Some(bases) => ("can", format!("(the bases `{bases}` match both{strand})")),
                None => ("may", TOO_COSTLY.to_string()),
            };
            let overlap = format!(
                "group `{name}`: tags `{a_s}` and `{b_s}` {verb} both match one read within dist={budget} {witness}"
            );
            if !nearest && !needs_indel {
                let exact = exact_pair(seqs, &tags, &flipped, &exact_only, shape, EXACT_PAIRS);
                if let Some((i, j, found)) = exact {
                    bail!("{}", indistinguishable(name, &seqs[i], &seqs[j], &found));
                }
                bail!("{overlap}; with mode=all this silently mis-assigns reads: separate the tags, lower dist, or use mode=nearest, which only warns and lists the colliding pairs");
            }
            colliding += 1;
            compared += usize::from(found.bases.is_some());
            if search.mode.is_none() {
                let first_compared = found.bases.is_some() && compared == 1;
                if example.is_none() || first_compared {
                    example = Some(format!("`{a_s}` and `{b_s}` {witness}"));
                }
            } else if colliding <= LISTED_PAIRS {
                warnings.push(format!("{overlap}; {consequence}"));
            }
        }
    }
    let complete = colliding < COUNTED_PAIRS;
    let pairs = |count: usize, kind: &str| {
        let at_least = if complete { "" } else { "at least " };
        let noun = if count == 1 { "pair" } else { "pairs" };
        format!("{at_least}{count} {kind}{noun}")
    };
    let verb = if compared > 0 { "can" } else { "may" };
    if let Some(pair) = example {
        bail!(
            "group `{name}`: {} of tags {verb} both match one read within dist={budget} only through an indel, such as {pair}; {NEAREST}, or write mode=all to keep assigning each to one of the tags by tie-break",
            pairs(colliding, "")
        );
    }
    if colliding > LISTED_PAIRS {
        let stopped = match complete {
            true => String::new(),
            false => format!(" (counting stopped at {COUNTED_PAIRS})"),
        };
        warnings.push(format!(
            "group `{name}`: {}{stopped} of tags {verb} both match one read within dist={budget}; {consequence}",
            pairs(colliding - LISTED_PAIRS, "more ")
        ));
    }
    if !nearest && colliding > 0 {
        warnings.push(format!(
            "group `{name}`: {} of tags {verb} match one read only through an indel; {NEAREST}",
            pairs(colliding, "")
        ));
    }

    Ok(warnings)
}

/// Two tags, by index, that one read matches with zero errors, with their
/// overlap on one strand when they have one. Plain ACGT tags of one length that
/// fit the window match that way only when equal or, under `bothStrands`,
/// reverse complements, so they are looked up by sequence; other sets are
/// compared pair by pair, up to `limit` pairs.
fn exact_pair(
    seqs: &[Vec<u8>],
    tags: &[Tag],
    flipped: &[Tag],
    exact_only: &Search,
    shape: Shape,
    limit: usize,
) -> Option<(usize, usize, Overlap)> {
    let compare = |(i, j): (usize, usize)| {
        let strands = exact_only.both_strands.then(|| (&flipped[i], &flipped[j]));
        let found = collision(&tags[i], &tags[j], None, exact_only, shape)
            .or_else(|| collision(&tags[i], &tags[j], strands, exact_only, shape))?;
        Some((i, j, found))
    };
    if looked_up(seqs, exact_only.window) {
        return equal_pair(seqs, exact_only.both_strands).and_then(compare);
    }
    (0..tags.len())
        .flat_map(|i| (i + 1..tags.len()).map(move |j| (i, j)))
        .take(limit)
        .find_map(compare)
}

/// Whether `seqs` are all plain ACGT of one length that fits `window`.
fn looked_up(seqs: &[Vec<u8>], window: Option<usize>) -> bool {
    let Some(len) = seqs.first().map(Vec::len) else {
        return false;
    };
    let plain =
        |seq: &Vec<u8>| seq.len() == len && seq.iter().all(|&base| mask(base).count_ones() == 1);
    window.is_none_or(|width| width >= len) && seqs.iter().all(plain)
}

/// The first of `seqs` equal to an earlier one or, with `both_strands`, to its
/// reverse complement, by index after that earlier one.
fn equal_pair(seqs: &[Vec<u8>], both_strands: bool) -> Option<(usize, usize)> {
    fn canonical(seq: &[u8]) -> Cow<'_, [u8]> {
        match seq.iter().all(|base| b"ACGT".contains(base)) {
            true => Cow::Borrowed(seq),
            false => seq
                .iter()
                .map(|&base| match base.to_ascii_uppercase() {
                    b'U' => b'T',
                    base => base,
                })
                .collect(),
        }
    }
    let mut seen: HashMap<Cow<'_, [u8]>, usize> = HashMap::with_capacity(seqs.len());
    for (j, seq) in seqs.iter().enumerate() {
        let forward = canonical(seq);
        let mut back = Vec::new();
        if both_strands {
            back.extend_from_slice(&forward);
            reverse_complement(&mut back);
        }
        let earlier = seen.get(&*forward).or_else(|| match both_strands {
            true => seen.get(back.as_slice()),
            false => None,
        });
        if let Some(&i) = earlier {
            return Some((i, j));
        }
        seen.insert(forward, j);
    }
    None
}

/// The refusal for two tags that one read matches with zero errors, which
/// points at `bothStrands` when only opposite strands match.
fn indistinguishable(name: &str, a: &[u8], b: &[u8], found: &Overlap) -> String {
    let (a_s, b_s) = (String::from_utf8_lossy(a), String::from_utf8_lossy(b));
    let bases = String::from_utf8_lossy(found.bases.as_deref().unwrap_or_default());
    if !found.flipped {
        return format!("group `{name}`: tags `{a_s}` and `{b_s}` are indistinguishable (the bases `{bases}` match both with zero errors); remove one or make them distinct");
    }
    let relation = match a.len().cmp(&b.len()) {
        Ordering::Equal => {
            format!("tags `{a_s}` and `{b_s}` are reverse complements of each other")
        }
        Ordering::Less => format!("tag `{a_s}` lies within the reverse complement of `{b_s}`"),
        Ordering::Greater => format!("tag `{b_s}` lies within the reverse complement of `{a_s}`"),
    };
    format!("group `{name}`: {relation} (the bases `{bases}` match both with zero errors, one on each strand), so bothStrands=true cannot tell them apart; drop bothStrands to match each tag on the forward strand only")
}

/// The `dist` attribute as written, leaving out the fields at their defaults.
fn dist_spec(dist: &Dist) -> String {
    let Dist {
        mismatch,
        indel,
        total,
    } = *dist;
    match (indel, total) {
        (0, total) if total == mismatch => mismatch.to_string(),
        (indel, total) if total == mismatch + indel => format!("{mismatch}:{indel}"),
        _ => format!("{mismatch}:{indel}:{total}"),
    }
}

/// A tag's IUPAC base sets, with the positions meeting each set as bits in
/// words of 64 positions, and its 3-grams: the concrete ones as bits, the
/// windows that may match more than those, and the 3-grams it can hold.
struct Tag {
    masks: Vec<u8>,
    sites: Vec<[u64; 16]>,
    concrete: [u32; 4],
    admitting: [u32; 4],
    empty: u32,
    windows: usize,
    grams: u64,
    loose: usize,
    holds: u64,
}

impl Tag {
    fn new(seq: &[u8]) -> Self {
        let masks: Vec<u8> = seq.iter().map(|&base| mask(base)).collect();
        let mut sites = vec![[0u64; 16]; masks.len().div_ceil(64)];
        let (mut concrete, mut admitting) = ([0; 4], [0; 4]);
        for (at, &set) in masks.iter().enumerate() {
            for (meets, bits) in sites[at / 64].iter_mut().enumerate() {
                if usize::from(set) & meets != 0 {
                    *bits |= 1 << (at % 64);
                }
            }
            for (base, count) in admitting.iter_mut().enumerate() {
                *count += u32::from(set >> base & 1);
            }
            if set.count_ones() == 1 {
                concrete[set.trailing_zeros() as usize] += 1;
            }
        }
        let (mut grams, mut kept, mut wild) = (0u64, 0, 0);
        for window in masks.windows(3).filter(|window| !window.contains(&0)) {
            if window.iter().any(|set| set.count_ones() > 1) {
                wild += 1;
            } else {
                let code = window
                    .iter()
                    .fold(0, |code, set| code << 2 | set.trailing_zeros());
                grams |= 1 << code;
                kept += 1;
            }
        }
        let degenerate = masks.iter().any(|set| set.count_ones() > 1);
        Self {
            empty: masks.iter().map(|&set| u32::from(set == 0)).sum(),
            sites,
            windows: masks.len().saturating_sub(2),
            loose: wild + kept - grams.count_ones() as usize,
            holds: if degenerate { u64::MAX } else { grams },
            grams,
            masks,
            concrete,
            admitting,
        }
    }

    fn len(&self) -> usize {
        self.masks.len()
    }

    /// The positions meeting each set, when the tag fits one word.
    fn word(&self) -> Option<&[u64; 16]> {
        match self.sites.as_slice() {
            [sites] => Some(sites),
            _ => None,
        }
    }

    /// The 64 positions from `from` on that meet `set`, as bits, with none
    /// outside the tag.
    fn bits(&self, set: usize, from: isize) -> u64 {
        let word = |at: isize| {
            let sites = usize::try_from(at).ok().and_then(|at| self.sites.get(at));
            sites.map_or(0, |sites| sites[set])
        };
        let (at, shift) = (from.div_euclid(64), from.rem_euclid(64) as u32);
        word(at) >> shift | word(at + 1).checked_shl(64 - shift).unwrap_or(0)
    }
}

/// Where one tag may sit within the bases another tag spans: anywhere, or
/// pinned at their start or their end.
#[derive(Clone, Copy)]
struct Shape {
    free_start: bool,
    free_end: bool,
}

impl Shape {
    const PINNED: Shape = Shape {
        free_start: false,
        free_end: false,
    };

    fn of(anchor: Option<Anchor>) -> Self {
        Shape {
            free_start: anchor != Some(Anchor::FivePrime),
            free_end: anchor != Some(Anchor::ThreePrime),
        }
    }
}

/// How uneven and how large two tags' error counts are, then how long the
/// bases are; smaller is a clearer example of a read that matches both.
type Spread = (usize, usize, usize);

/// Bases one read can carry that two tags both match.
struct Overlap {
    /// The bases, or `None` when the tags were too costly to compare exactly
    /// and are taken to collide through an indel.
    bases: Option<Vec<u8>>,
    exact: bool,
    flipped: bool,
    spread: Spread,
}

impl Overlap {
    /// Two tags taken to collide, ranked behind any overlap with bases.
    fn unchecked() -> Self {
        Overlap {
            bases: None,
            exact: false,
            flipped: false,
            spread: (usize::MAX, usize::MAX, usize::MAX),
        }
    }
}

/// The overlap of two tags over the orders and strands the search allows,
/// preferring one with zero errors, then the most even one.
fn collision(
    a: &Tag,
    b: &Tag,
    flipped: Option<(&Tag, &Tag)>,
    search: &Search,
    shape: Shape,
) -> Option<Overlap> {
    let symmetric = search.dist.indel == 0 && a.len() == b.len();
    let (flipped_a, flipped_b) = flipped.unzip();
    let orders = [
        Some((a, b, false)),
        (!symmetric).then_some((b, a, false)),
        flipped_b.map(|flipped_b| (a, flipped_b, true)),
        flipped_a
            .filter(|_| !symmetric)
            .map(|flipped_a| (b, flipped_a, true)),
    ];
    orders
        .into_iter()
        .flatten()
        .filter_map(|(outer, inner, flipped)| {
            overlap(outer, inner, search, shape).map(|found| Overlap { flipped, ..found })
        })
        .min_by_key(|found| (!found.exact, found.spread))
}

/// Bases that `outer` matches end to end and `inner` matches within, checked
/// against the budget.
fn overlap(outer: &Tag, inner: &Tag, search: &Search, shape: Shape) -> Option<Overlap> {
    let dist = &search.dist;
    let reach = outer.len() + dist.indel;
    let longest = search.window.map_or(reach, |width| width.min(reach));
    if inner.len() == 0 || inner.len().max(outer.len()) > longest + dist.indel {
        return None;
    }
    let found = match dist.indel {
        0 => ungapped(outer, inner, dist, shape, longest)?,
        _ => match gapped(outer, inner, dist, shape, longest)? {
            Gapped::Found(found) => found,
            Gapped::TooLarge => return Some(Overlap::unchecked()),
        },
    };
    let bases = found.bases(outer, inner);
    let fits = aligns(&outer.masks, &bases, dist, Shape::PINNED)
        && aligns(&inner.masks, &bases, dist, shape);
    fits.then_some(Overlap {
        bases: Some(bases),
        exact: found.exact,
        flipped: false,
        spread: found.spread,
    })
}

/// One column of an alignment of an inner tag within an outer one, by outer
/// and inner position: an outer base outside the inner tag, a pair, an outer
/// base the inner tag lacks, and an inner base the outer tag lacks.
#[derive(Clone, Copy)]
enum Step {
    Flank(usize),
    Pair(usize, usize),
    Skip(usize),
    Extra(usize),
}

/// The disjoint pairs, skipped outer bases, and extra inner bases charged to
/// the outer tag; the rest are charged to the inner one.
#[derive(Clone, Copy)]
struct Split {
    subs: usize,
    skipped: usize,
    extra: usize,
}

/// An alignment of an inner tag within an outer one and how its edits split.
struct Candidate {
    steps: Vec<Step>,
    split: Split,
    spread: Spread,
    exact: bool,
}

impl Candidate {
    /// The bases where each edit counts against only the tag it is charged to.
    fn bases(&self, outer: &Tag, inner: &Tag) -> Vec<u8> {
        let Split {
            mut subs,
            mut skipped,
            mut extra,
        } = self.split;
        let mut bases = Vec::with_capacity(outer.len() + extra);
        for &step in &self.steps {
            match step {
                Step::Flank(p) => bases.push(base(outer.masks[p])),
                Step::Pair(p, j) => {
                    let (x, y) = (outer.masks[p], inner.masks[j]);
                    let set = if x & y != 0 {
                        x & y
                    } else if subs > 0 {
                        subs -= 1;
                        y
                    } else {
                        x
                    };
                    bases.push(base(set));
                }
                Step::Skip(_) if skipped > 0 => skipped -= 1,
                Step::Skip(p) => bases.push(base(outer.masks[p])),
                Step::Extra(j) if extra > 0 => {
                    extra -= 1;
                    bases.push(base(inner.masks[j]));
                }
                Step::Extra(_) => {}
            }
        }
        bases
    }
}

/// A concrete base in the IUPAC set `set`, or `A` for the empty set.
fn base(set: u8) -> u8 {
    match set.trailing_zeros() {
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'A',
    }
}

/// The placement of `inner` within `outer` with the fewest disjoint pairs and
/// no indels.
fn ungapped(
    outer: &Tag,
    inner: &Tag,
    dist: &Dist,
    shape: Shape,
    longest: usize,
) -> Option<Candidate> {
    let (n, k) = (outer.len(), inner.len());
    if n > longest || k > n {
        return None;
    }
    let last = n - k;
    let offsets = match (shape.free_start, shape.free_end) {
        (false, _) => 0..=0,
        (true, false) => last..=last,
        (true, true) => 0..=last,
    };
    let (offset, subs) = offsets
        .map(|offset| (offset, disjoint(outer, inner, offset)))
        .min_by_key(|&(_, subs)| subs)?;
    if subs > 2 * dist.mismatch {
        return None;
    }
    let (spread, split) = split_edits(subs, 0, 0, dist, n, longest)?;
    let steps = (0..offset)
        .map(Step::Flank)
        .chain((0..k).map(|j| Step::Pair(offset + j, j)))
        .chain((offset + k..n).map(Step::Flank))
        .collect();
    Some(Candidate {
        steps,
        split,
        spread,
        exact: subs == 0,
    })
}

/// How many bases of `inner` share no base with those of `outer` from `offset`.
fn disjoint(outer: &Tag, inner: &Tag, offset: usize) -> usize {
    match (outer.word(), inner.word()) {
        (Some(x), Some(y)) => {
            let met = [1, 2, 4, 8]
                .into_iter()
                .fold(0, |met, set| met | (x[set] >> offset) & y[set]);
            inner.len() - met.count_ones() as usize
        }
        _ => inner
            .masks
            .iter()
            .zip(&outer.masks[offset..])
            .filter(|&(y, x)| x & y == 0)
            .count(),
    }
}

/// What the table over two tags' alignments found, when there is one.
enum Gapped {
    /// An alignment whose edits split within the budget.
    Found(Candidate),
    /// No table was filled, as it would hold more than [`TABLE_CELLS`] cells or
    /// more than the group has left of [`GROUP_CELLS`].
    TooLarge,
}

/// The alignment of `inner` within `outer`, indels allowed, whose edits split
/// most evenly within `dist`.
fn gapped(outer: &Tag, inner: &Tag, dist: &Dist, shape: Shape, longest: usize) -> Option<Gapped> {
    let budget = 2 * dist.total;
    if unmatched(outer, inner) > budget
        || !shares_grams(outer, inner, budget)
        || !covered(outer, inner, shape, 2 * dist.indel, budget)
        || !edits_within(outer, inner, shape, budget)
    {
        return None;
    }
    tabled(outer, inner, dist, shape, longest)
}

/// The alignment that a table finds for a pair the filters let through, kept
/// out of line so that the path every pair takes stays small.
#[inline(never)]
fn tabled(outer: &Tag, inner: &Tag, dist: &Dist, shape: Shape, longest: usize) -> Option<Gapped> {
    let (size, left) = (Grid::size(outer, inner, dist, shape), CELLS_LEFT.get());
    if size > TABLE_CELLS || size > left {
        return Some(Gapped::TooLarge);
    }
    CELLS_LEFT.set(left - size);
    let grid = Grid::fill(outer, inner, dist, shape);
    let (gaps, skips, extras) = (grid.gaps, grid.skips, grid.extras);
    let ends = match shape.free_end {
        true => 0..=outer.len(),
        false => outer.len()..=outer.len(),
    };
    let (spread, split, end, exact) = ends
        .flat_map(|end| {
            (0..=skips).flat_map(move |skipped| {
                (0..=extras.min(gaps - skipped)).map(move |extra| (end, skipped, extra))
            })
        })
        .filter_map(|(end, skipped, extra)| {
            let subs = grid.get(inner.len(), end, skipped, extra)?;
            let (spread, split) = split_edits(subs, skipped, extra, dist, outer.len(), longest)?;
            Some((
                spread,
                split,
                (end, skipped, extra),
                subs + skipped + extra == 0,
            ))
        })
        .min_by_key(|&(spread, ..)| spread)?;
    Some(Gapped::Found(Candidate {
        steps: grid.trace(outer, inner, shape, end),
        split,
        spread,
        exact,
    }))
}

/// How many bases of `inner` no base of `outer` can pair with, a lower bound on
/// the edits of any alignment of the two.
fn unmatched(outer: &Tag, inner: &Tag) -> usize {
    let short = (0..4).map(|base| inner.concrete[base].saturating_sub(outer.admitting[base]));
    (inner.empty + short.sum::<u32>()) as usize
}

/// Whether `inner` can have enough 3-grams that `outer` holds for an alignment
/// with at most `budget` edits, as each edit breaks at most three of them.
fn shares_grams(outer: &Tag, inner: &Tag, budget: usize) -> bool {
    let held = (inner.grams & outer.holds).count_ones() as usize + inner.loose;
    held + 3 * budget >= inner.windows
}

/// Whether `gaps + 1` neighboring diagonals, where an alignment with at most
/// `gaps` unpaired bases stays, pair all but at most `budget` bases of `inner`
/// with outer bases they meet; assumed when that takes many diagonals.
fn covered(outer: &Tag, inner: &Tag, shape: Shape, gaps: usize, budget: usize) -> bool {
    let gaps = gaps.min(outer.len() + inner.len()) as isize;
    let last = match shape.free_start {
        true => outer.len() as isize - inner.len() as isize,
        false => 0,
    };
    let words = inner.sites.len();
    let diagonals = usize::try_from(last + 2 * gaps + 1).unwrap_or(0);
    if gaps > 8 || diagonals * words > 1024 {
        return true;
    }
    let met: SmallVec<[u64; 32]> = (-gaps..=last + gaps)
        .flat_map(|diagonal| (0..words).map(move |word| (diagonal, word)))
        .map(|(diagonal, word)| {
            let from = 64 * word as isize + diagonal;
            [1, 2, 4, 8].into_iter().fold(0, |met, set| {
                met | outer.bits(set, from) & inner.sites[word][set]
            })
        })
        .collect();
    (-gaps..=last).any(|first| {
        let start = (first + gaps) as usize * words;
        let paired = (0..words).map(|word| {
            let window = (0..=gaps as usize).map(|offset| met[start + offset * words + word]);
            window.fold(0u64, |paired, met| paired | met).count_ones() as usize
        });
        inner.len() - paired.sum::<usize>() <= budget
    })
}

/// Whether an alignment of `inner` within `outer` has at most `budget`
/// disjoint pairs and unpaired bases, by Myers' bit-vector algorithm when
/// `inner` fits one word.
fn edits_within(outer: &Tag, inner: &Tag, shape: Shape, budget: usize) -> bool {
    let Some(sites) = inner.word() else {
        return banded_within(outer, inner, shape, budget);
    };
    let last = 1u64 << (inner.len() - 1);
    let (mut pv, mut mv, mut score) = (u64::MAX, 0u64, inner.len());
    let mut fewest = score;
    for &set in &outer.masks {
        let eq = sites[usize::from(set & 15)];
        let xv = eq | mv;
        let xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq;
        let ph = mv | !(xh | pv);
        let mh = pv & xh;
        score = score + usize::from(ph & last != 0) - usize::from(mh & last != 0);
        let ph = ph << 1 | u64::from(!shape.free_start);
        pv = mh << 1 | !(xv | ph);
        mv = ph & xv;
        fewest = fewest.min(score);
    }
    match shape.free_end {
        true => fewest <= budget,
        false => score <= budget,
    }
}

/// Whether an alignment of `inner` within `outer` has at most `budget`
/// disjoint pairs and unpaired bases, by a table over the diagonals that such
/// an alignment can reach.
fn banded_within(outer: &Tag, inner: &Tag, shape: Shape, budget: usize) -> bool {
    let (n, reach) = (outer.len() as isize, budget as isize);
    let lo = -reach;
    let hi = reach
        + match shape.free_start {
            true => n - inner.len() as isize,
            false => 0,
        };
    let Ok(width) = usize::try_from(hi - lo + 1) else {
        return false;
    };
    let mut row: Vec<usize> = (lo..=hi)
        .map(|p| match (0..=n).contains(&p) {
            true if shape.free_start => 0,
            true => p as usize,
            false => usize::MAX,
        })
        .collect();
    let mut next = vec![usize::MAX; width];
    for (j, &base) in inner.masks.iter().enumerate() {
        for (slot, diagonal) in (lo..=hi).enumerate() {
            let p = j as isize + 1 + diagonal;
            next[slot] = usize::MAX;
            if !(0..=n).contains(&p) {
                continue;
            }
            if p > 0 {
                let missed = usize::from(outer.masks[p as usize - 1] & base == 0);
                next[slot] = row[slot].saturating_add(missed);
            }
            if slot + 1 < width {
                next[slot] = next[slot].min(row[slot + 1].saturating_add(1));
            }
            if slot > 0 && p > 0 {
                next[slot] = next[slot].min(next[slot - 1].saturating_add(1));
            }
        }
        if next.iter().all(|&edits| edits > budget) {
            return false;
        }
        std::mem::swap(&mut row, &mut next);
    }
    let end = n - inner.len() as isize - lo;
    match shape.free_end {
        true => row.iter().any(|&edits| edits <= budget),
        false => usize::try_from(end)
            .ok()
            .and_then(|slot| row.get(slot))
            .is_some_and(|&edits| edits <= budget),
    }
}

/// The fewest disjoint pairs over the alignments of an inner tag within an
/// outer one, by inner bases placed, outer bases passed, and the outer and
/// inner bases left unpaired, on the diagonals that such alignments can reach.
struct Grid {
    lo: isize,
    width: usize,
    span: usize,
    skips: usize,
    extras: usize,
    gaps: usize,
    cells: Vec<usize>,
}

impl Grid {
    fn index(&self, j: usize, p: usize, skipped: usize, extra: usize) -> Option<usize> {
        let offset = usize::try_from(p as isize - j as isize - self.lo).ok();
        let offset = offset.filter(|&offset| offset < self.width && p <= self.span)?;
        Some(((j * self.width + offset) * (self.skips + 1) + skipped) * (self.extras + 1) + extra)
    }

    fn get(&self, j: usize, p: usize, skipped: usize, extra: usize) -> Option<usize> {
        let subs = self.cells[self.index(j, p, skipped, extra)?];
        (subs != usize::MAX).then_some(subs)
    }

    fn paired(
        &self,
        outer: &Tag,
        inner: &Tag,
        j: usize,
        p: usize,
        skipped: usize,
        extra: usize,
    ) -> Option<usize> {
        let subs = self.get(j - 1, p - 1, skipped, extra)?;
        Some(subs + usize::from(outer.masks[p - 1] & inner.masks[j - 1] == 0))
    }

    /// The lowest and highest diagonals, and the most skipped and extra bases,
    /// that the alignments of `inner` within `outer` can reach.
    fn bounds(outer: &Tag, inner: &Tag, dist: &Dist, shape: Shape) -> (isize, isize, usize, usize) {
        let gaps = 2 * dist.indel;
        let (skips, extras) = (gaps.min(outer.len()), gaps.min(inner.len()));
        let hi = match shape.free_start {
            true => outer.len() as isize - inner.len() as isize + extras as isize,
            false => skips as isize,
        };
        (-(extras as isize), hi, skips, extras)
    }

    /// How many cells the table of `inner` within `outer` holds.
    fn size(outer: &Tag, inner: &Tag, dist: &Dist, shape: Shape) -> usize {
        let (lo, hi, skips, extras) = Grid::bounds(outer, inner, dist, shape);
        let width = usize::try_from(hi - lo + 1).unwrap_or(0);
        [inner.len() + 1, width, skips + 1, extras + 1]
            .into_iter()
            .fold(1, usize::saturating_mul)
    }

    fn fill(outer: &Tag, inner: &Tag, dist: &Dist, shape: Shape) -> Self {
        let gaps = 2 * dist.indel;
        let (lo, hi, skips, extras) = Grid::bounds(outer, inner, dist, shape);
        let width = usize::try_from(hi - lo + 1).unwrap_or(0);
        let size = (inner.len() + 1) * width * (skips + 1) * (extras + 1);
        let mut grid = Grid {
            lo,
            width,
            span: outer.len(),
            skips,
            extras,
            gaps,
            cells: vec![usize::MAX; size],
        };
        for j in 0..=inner.len() {
            let first = (j as isize + lo).max(0);
            let last = (j as isize + hi).min(outer.len() as isize);
            for p in (first..=last).map(|p| p as usize) {
                for skipped in 0..=skips {
                    for extra in 0..=extras.min(gaps - skipped) {
                        let start = j == 0 && skipped + extra == 0 && (p == 0 || shape.free_start);
                        let mut subs = if start { 0 } else { usize::MAX };
                        if j > 0 && p > 0 {
                            let paired = grid.paired(outer, inner, j, p, skipped, extra);
                            subs = subs.min(paired.unwrap_or(usize::MAX));
                        }
                        if p > 0 && skipped > 0 {
                            let from = grid.get(j, p - 1, skipped - 1, extra);
                            subs = subs.min(from.unwrap_or(usize::MAX));
                        }
                        if j > 0 && extra > 0 {
                            let from = grid.get(j - 1, p, skipped, extra - 1);
                            subs = subs.min(from.unwrap_or(usize::MAX));
                        }
                        if subs > 2 * dist.mismatch || subs + skipped + extra > 2 * dist.total {
                            subs = usize::MAX;
                        }
                        if let Some(at) = grid.index(j, p, skipped, extra) {
                            grid.cells[at] = subs;
                        }
                    }
                }
            }
        }
        grid
    }

    /// The steps of an alignment that ends at `end`, an outer position and its
    /// skipped and extra bases, with the fewest disjoint pairs.
    fn trace(
        &self,
        outer: &Tag,
        inner: &Tag,
        shape: Shape,
        end: (usize, usize, usize),
    ) -> Vec<Step> {
        let (mut p, mut skipped, mut extra) = end;
        let mut j = inner.len();
        let mut steps: Vec<Step> = (p..outer.len()).rev().map(Step::Flank).collect();
        while j > 0 || skipped + extra > 0 || !(p == 0 || shape.free_start) {
            let subs = self.get(j, p, skipped, extra);
            if j > 0 && p > 0 && self.paired(outer, inner, j, p, skipped, extra) == subs {
                steps.push(Step::Pair(p - 1, j - 1));
                (j, p) = (j - 1, p - 1);
            } else if p > 0 && skipped > 0 && self.get(j, p - 1, skipped - 1, extra) == subs {
                steps.push(Step::Skip(p - 1));
                (p, skipped) = (p - 1, skipped - 1);
            } else {
                steps.push(Step::Extra(j - 1));
                (j, extra) = (j - 1, extra - 1);
            }
        }
        steps.extend((0..p).rev().map(Step::Flank));
        steps.reverse();
        steps
    }
}

/// The most even way to charge an alignment's disjoint pairs, skipped outer
/// bases, and extra inner bases to the two tags within `dist`, with the bases
/// it leaves between one and `longest` long.
fn split_edits(
    subs: usize,
    skipped: usize,
    extra: usize,
    dist: &Dist,
    outer_len: usize,
    longest: usize,
) -> Option<(Spread, Split)> {
    let mut best: Option<(Spread, Split)> = None;
    for outer_skipped in 0..=skipped {
        for outer_extra in 0..=extra {
            let len = (outer_len + outer_extra).saturating_sub(outer_skipped);
            if len == 0 || len > longest {
                continue;
            }
            let outer_gaps = outer_skipped + outer_extra;
            let inner_gaps = skipped + extra - outer_gaps;
            if outer_gaps.max(inner_gaps) > dist.indel.min(dist.total) {
                continue;
            }
            let fewest = subs
                .saturating_sub(dist.mismatch)
                .max(subs.saturating_sub(dist.total - inner_gaps));
            let most = subs.min(dist.mismatch).min(dist.total - outer_gaps);
            for outer_subs in fewest..=most {
                let charged = outer_subs + outer_gaps;
                let rest = subs - outer_subs + inner_gaps;
                let spread = (charged.abs_diff(rest), charged.max(rest), len);
                if best.is_none_or(|(seen, _)| spread < seen) {
                    let split = Split {
                        subs: outer_subs,
                        skipped: outer_skipped,
                        extra: outer_extra,
                    };
                    best = Some((spread, split));
                }
            }
        }
    }
    best
}

/// Whether `tag`, as IUPAC sets, aligns to `bases` within `dist`, with the
/// bases before and after it free as `shape` allows.
fn aligns(tag: &[u8], bases: &[u8], dist: &Dist, shape: Shape) -> bool {
    let gaps = dist.indel.min(tag.len() + bases.len()) + 1;
    let (rows, cols) = (tag.len() + 1, bases.len() + 1);
    let index = |j: usize, i: usize, q: usize| (j * cols + i) * gaps + q;
    let mut fewest = vec![usize::MAX; rows * cols * gaps];
    for j in 0..rows {
        for i in 0..cols {
            for q in 0..gaps {
                let start = j == 0 && q == 0 && (i == 0 || shape.free_start);
                let mut subs = if start { 0 } else { usize::MAX };
                if j > 0 && i > 0 {
                    let missed = usize::from(tag[j - 1] & mask(bases[i - 1]) == 0);
                    subs = subs.min(fewest[index(j - 1, i - 1, q)].saturating_add(missed));
                }
                if q > 0 && j > 0 {
                    subs = subs.min(fewest[index(j - 1, i, q - 1)]);
                }
                if q > 0 && i > 0 {
                    subs = subs.min(fewest[index(j, i - 1, q - 1)]);
                }
                fewest[index(j, i, q)] = subs;
            }
        }
    }
    let ends = match shape.free_end {
        true => 0..cols,
        false => cols - 1..cols,
    };
    ends.flat_map(|i| (0..gaps).map(move |q| (i, q)))
        .any(|(i, q)| {
            let subs = fewest[index(rows - 1, i, q)];
            subs <= dist.mismatch && subs + q <= dist.total
        })
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

    fn search(dist: Dist, mode: MatchMode) -> Search {
        Search {
            mode: Some(mode),
            ..unset(dist)
        }
    }

    fn unset(dist: Dist) -> Search {
        Search {
            dist,
            mode: None,
            anchor: None,
            both_strands: false,
            window: None,
        }
    }

    fn check(tags: &[&str], dist: Dist, mode: MatchMode) -> Result<Vec<String>> {
        validate_group("g", &seqs(tags), &search(dist, mode))
    }

    fn budget(mismatch: usize, indel: usize) -> Dist {
        Dist {
            mismatch,
            indel,
            total: mismatch + indel,
        }
    }

    struct Rng(u64);

    impl Rng {
        fn below(&mut self, n: usize) -> usize {
            self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            ((z ^ (z >> 31)) % n as u64) as usize
        }

        fn sample(&mut self, len: usize, letters: &[u8]) -> Vec<u8> {
            (0..len)
                .map(|_| letters[self.below(letters.len())])
                .collect()
        }
    }

    /// Random barcodes, each at least `apart` substitutions from every other.
    fn spread_barcodes(count: usize, len: usize, apart: usize) -> Vec<Vec<u8>> {
        let mut rng = Rng(0x2545_f491);
        let mut tags: Vec<Vec<u8>> = Vec::new();
        while tags.len() < count {
            let tag = rng.sample(len, b"ACGT");
            let far = |other: &Vec<u8>| other.iter().zip(&tag).filter(|(a, b)| a != b).count();
            if tags.iter().all(|other| far(other) >= apart) {
                tags.push(tag);
            }
        }
        tags
    }

    #[test]
    fn test_indel_budget_overlap_is_a_warning_under_either_mode_written_out() {
        // Deleting the last base of either tag leaves ACGTA, one indel from both.
        let tags = ["ACGTAC", "ACGTAG"];
        let pair = "within dist=0:1 (the bases `ACGTA` match both)";
        let w = check(&tags, budget(0, 1), MatchMode::All).unwrap();
        assert_eq!(w.len(), 2, "a pair and a count: {w:?}");
        let tie_break = "with mode=all such a read goes to one of the tags by tie-break";
        assert!(w[0].contains(pair) && w[0].ends_with(tie_break), "{w:?}");
        let count = "group `g`: 1 pair of tags can match one read only through an indel; set mode=nearest to leave such reads unassigned (recommended)";
        assert_eq!(w[1], count);
        let w = check(&tags, budget(0, 1), MatchMode::Nearest).unwrap();
        assert_eq!(w.len(), 1, "a pair: {w:?}");
        let unassigned = "mode=nearest leaves ambiguous reads unassigned";
        assert!(w[0].contains(pair) && w[0].contains(unassigned), "{w:?}");
    }

    #[test]
    fn test_indel_budget_overlap_without_a_mode_asks_for_one() {
        let err = validate_group("g", &seqs(&["ACGTAC", "ACGTAG"]), &unset(budget(0, 1)))
            .unwrap_err()
            .to_string();
        let pair = "group `g`: 1 pair of tags can both match one read within dist=0:1 only through an indel, such as `ACGTAC` and `ACGTAG` (the bases `ACGTA` match both)";
        let hint = "set mode=nearest to leave such reads unassigned (recommended), or write mode=all to keep assigning each to one of the tags by tie-break";
        assert_eq!(err, format!("{pair}; {hint}"));
    }

    #[test]
    fn test_without_a_mode_the_refusal_counts_the_pairs_that_need_an_indel() {
        let tags = spread_barcodes(96, 8, 3);
        let err = validate_group("g", &tags, &unset(budget(1, 1)))
            .unwrap_err()
            .to_string();
        let count = "group `g`: at least 1000 pairs of tags can both match one read within dist=1:1 only through an indel, such as";
        assert!(err.starts_with(count), "{err}");
    }

    #[test]
    fn test_without_a_mode_substitutions_alone_refuse_a_pair_as_under_mode_all() {
        let err = validate_group("g", &seqs(&["AAAACCCC", "AAAAGGCC"]), &unset(budget(1, 1)))
            .unwrap_err()
            .to_string();
        let refusal =
            "(the bases `AAAAGCCC` match both); with mode=all this silently mis-assigns reads";
        assert!(err.contains(refusal), "{err}");
    }

    #[test]
    fn test_substitutions_alone_still_refuse_a_pair_under_mode_all_with_an_indel_budget() {
        // AAAAGCCC is one substitution from each tag, so the read needs no indel.
        let tags = ["AAAACCCC", "AAAAGGCC"];
        let err = check(&tags, budget(1, 1), MatchMode::All)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("within dist=1:1 (the bases `AAAAGCCC` match both)"),
            "{err}"
        );
        assert!(
            err.contains("with mode=all this silently mis-assigns reads"),
            "{err}"
        );
        let w = check(&tags, budget(1, 1), MatchMode::Nearest).unwrap();
        assert_eq!(w.len(), 1, "no count of pairs that need an indel: {w:?}");
    }

    #[test]
    fn test_pair_matched_without_an_indel_is_shown_with_bases_that_need_none() {
        // ACGTA is one deletion from both tags, and ACGTAC is one substitution from ACGTAG.
        let err = check(&["ACGTAC", "ACGTAG"], budget(1, 1), MatchMode::All)
            .unwrap_err()
            .to_string();
        let bases = "within dist=1:1 (the bases `ACGTAC` match both)";
        assert!(err.contains(bases), "{err}");
    }

    #[test]
    fn test_mode_all_lists_the_pairs_that_need_an_indel_and_counts_them() {
        // With no substitutions every pair needs an indel: AAAA is one deletion from each tag.
        let tags = ["AAAAA", "AAAAC", "AAAAG", "AAAAT", "AAACA", "AAAGA"];
        let w = check(&tags, budget(0, 1), MatchMode::All).unwrap();
        assert_eq!(w.len(), 12, "ten pairs and two counts: {w:?}");
        let tie_break = "with mode=all such a read goes to one of the tags by tie-break";
        assert!(w[..10].iter().all(|w| w.ends_with(tie_break)), "{w:?}");
        let more = "group `g`: 5 more pairs of tags can both match one read within dist=0:1;";
        assert!(w[10].starts_with(more), "{w:?}");
        let count = "group `g`: 15 pairs of tags can match one read only through an indel; set mode=nearest to leave such reads unassigned (recommended)";
        assert_eq!(w[11], count);
    }

    #[test]
    fn test_mode_all_stops_counting_but_still_refuses_pairs_without_an_indel() {
        // Barcodes three substitutions apart collide only through an indel.
        let mut tags = spread_barcodes(96, 8, 3);
        let w = validate_group("g", &tags, &search(budget(1, 1), MatchMode::All)).unwrap();
        assert_eq!(w.len(), 12, "ten pairs and two counts: {w:?}");
        let more = "at least 990 more pairs (counting stopped at 1000)";
        assert!(w[10].contains(more), "{w:?}");
        let count = "at least 1000 pairs of tags can match one read only through an indel";
        assert!(w[11].contains(count), "{w:?}");
        let other = |base: u8| if base == b'C' { b'A' } else { b'C' };
        let last = tags[95].clone();
        let mut near = last.clone();
        near[0] = other(near[0]);
        near[7] = other(near[7]);
        tags.push(near.clone());
        let err = validate_group("g", &tags, &search(budget(1, 1), MatchMode::All))
            .unwrap_err()
            .to_string();
        let pair = format!(
            "tags `{}` and `{}` can both match one read within dist=1:1",
            String::from_utf8_lossy(&last),
            String::from_utf8_lossy(&near)
        );
        assert!(err.contains(&pair), "{err}");
        assert!(
            err.contains("with mode=all this silently mis-assigns reads"),
            "{err}"
        );
    }

    #[test]
    fn test_substitution_only_budget_tells_apart_tags_differing_at_one_base() {
        assert!(check(&["ACGTAC", "ACGTAG"], dist(0), MatchMode::All)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn test_tag_inside_a_longer_tag_is_indistinguishable() {
        // Every read that starts with ACGTAC starts with ACGT too.
        for tags in [["ACGT", "ACGTAC"], ["ACGTAC", "ACGT"]] {
            for mode in [MatchMode::All, MatchMode::Nearest] {
                let err = check(&tags, dist(0), mode).unwrap_err().to_string();
                assert!(
                    err.contains(
                        "are indistinguishable (the bases `ACGTAC` match both with zero errors)"
                    ),
                    "{err}"
                );
            }
        }
    }

    #[test]
    fn test_tag_near_part_of_a_longer_tag_is_mode_aware() {
        // ACGA is one substitution from the ACGT inside TTACGTTT.
        let tags = ["ACGA", "TTACGTTT"];
        let err = check(&tags, dist(1), MatchMode::All)
            .unwrap_err()
            .to_string();
        assert!(err.contains("within dist=1"), "{err}");
        assert!(err.contains("mode=all"), "{err}");
        let w = check(&tags, dist(1), MatchMode::Nearest).unwrap();
        assert_eq!(w.len(), 1, "one collision warning: {w:?}");
        assert!(check(&tags, dist(0), MatchMode::All).unwrap().is_empty());
    }

    #[test]
    fn test_well_separated_barcodes_pass() {
        // No read is within one substitution of two barcodes three substitutions apart.
        let tags = spread_barcodes(96, 8, 3);
        let w = validate_group("g", &tags, &search(dist(1), MatchMode::All)).unwrap();
        assert!(w.is_empty(), "well-separated barcodes are clean: {w:?}");
    }

    fn shaped(anchor: Option<Anchor>, both_strands: bool, window: Option<usize>) -> Search {
        Search {
            anchor,
            both_strands,
            window,
            ..search(dist(0), MatchMode::All)
        }
    }

    #[test]
    fn test_anchored_tags_are_compared_at_the_pinned_edge_only() {
        // CGTA sits one base into ACGTAC, ACGT starts it, and GTAC ends it.
        let (five, three) = (Some(Anchor::FivePrime), Some(Anchor::ThreePrime));
        for (tags, at_five, at_three) in [
            (["CGTA", "ACGTAC"], false, false),
            (["ACGT", "ACGTAC"], true, false),
            (["GTAC", "ACGTAC"], false, true),
        ] {
            let tags = seqs(&tags);
            let refused = |search: Search| validate_group("g", &tags, &search).is_err();
            assert!(refused(shaped(None, false, None)), "{tags:?} slide");
            assert_eq!(refused(shaped(five, false, None)), at_five, "{tags:?} 5p");
            assert_eq!(
                refused(shaped(three, false, Some(8))),
                at_three,
                "{tags:?} 3p"
            );
        }
    }

    #[test]
    fn test_both_strands_compares_each_tag_with_the_other_reverse_complemented() {
        let tags = seqs(&["ACGTACGG", "CCGTACGT"]);
        assert!(validate_group("g", &tags, &shaped(None, false, None))
            .unwrap()
            .is_empty());
        let err = validate_group("g", &tags, &shaped(None, true, None))
            .unwrap_err()
            .to_string();
        assert!(err.contains("one on each strand"), "{err}");
    }

    #[test]
    fn test_window_narrower_than_a_tag_holds_no_read_matching_it() {
        // Without indels an 8 bp tag never fits a 6 bp window.
        let tags = seqs(&["ACGTACGT", "ACGTACGA"]);
        let narrow = Search {
            dist: dist(1),
            ..shaped(None, false, Some(6))
        };
        assert!(validate_group("g", &tags, &narrow).unwrap().is_empty());
        let wide = Search {
            window: Some(8),
            ..narrow
        };
        assert!(validate_group("g", &tags, &wide).is_err());
    }

    #[test]
    fn test_bases_may_run_an_indel_past_the_longer_tag() {
        // ACTGTACGT is ACGTACGT with a base inserted and ACTGTACA with its last base deleted.
        let tags = seqs(&["ACGTACGT", "ACTGTACA"]);
        let five = Search {
            dist: budget(0, 1),
            ..shaped(Some(Anchor::FivePrime), false, None)
        };
        let w = validate_group("g", &tags, &five).unwrap();
        assert!(
            w[0].contains("within dist=0:1 (the bases `ACTGTACGT` match both)"),
            "{w:?}"
        );
    }

    /// Two 150 bp tags two substitutions apart, which one read matches through indels.
    fn long_pair() -> Vec<Vec<u8>> {
        let tag = Rng(0x1093).sample(150, b"ACGT");
        let mut other = tag.clone();
        for at in [40, 110] {
            other[at] = if other[at] == b'A' { b'C' } else { b'A' };
        }
        vec![tag, other]
    }

    #[test]
    fn test_long_tags_too_costly_to_compare_are_taken_to_collide() {
        let tags = long_pair();
        let w = validate_group("g", &tags, &search(budget(0, 8), MatchMode::Nearest)).unwrap();
        let unchecked = "may both match one read within dist=0:8 (too costly to compare exactly under this indel budget)";
        assert!(w.len() == 1 && w[0].contains(unchecked), "{w:?}");
        let err = validate_group("g", &tags, &unset(budget(0, 8)))
            .unwrap_err()
            .to_string();
        let refusal = "group `g`: 1 pair of tags may both match one read within dist=0:8 only through an indel, such as";
        assert!(err.starts_with(refusal), "{err}");
    }

    #[test]
    fn test_long_tags_under_a_smaller_indel_budget_are_compared_exactly() {
        let w = validate_group("g", &long_pair(), &search(budget(0, 7), MatchMode::Nearest));
        let w = w.unwrap();
        let checked = "can both match one read within dist=0:7 (the bases `";
        assert!(w.len() == 1 && w[0].contains(checked), "{w:?}");
    }

    #[test]
    fn test_reverse_complement_tags_point_at_both_strands() {
        // ATACCAACGC and GCGTTGGTAT, two Illumina UD i7 indexes, are reverse complements.
        let tags = seqs(&["ATACCAACGC", "GCGTTGGTAT"]);
        for mode in [MatchMode::All, MatchMode::Nearest] {
            let both = Search {
                both_strands: true,
                ..search(dist(0), mode)
            };
            let err = validate_group("g", &tags, &both).unwrap_err().to_string();
            let reason = "tags `ATACCAACGC` and `GCGTTGGTAT` are reverse complements of each other (the bases `ATACCAACGC` match both with zero errors, one on each strand), so bothStrands=true cannot tell them apart; drop bothStrands to match each tag on the forward strand only, or remove one of them";
            assert!(err.ends_with(reason), "{err}");
        }
        // ACGT is its own reverse complement, but it also lies within ACGTAC on one strand.
        let both = Search {
            both_strands: true,
            ..unset(dist(0))
        };
        let err = validate_group("g", &seqs(&["ACGT", "ACGTAC"]), &both)
            .unwrap_err()
            .to_string();
        let same = err.contains("are indistinguishable") && !err.contains("bothStrands");
        assert!(same, "{err}");
    }

    #[test]
    fn test_a_zero_error_pair_is_reported_before_a_pair_that_needs_substitutions() {
        // AAAA and AATT differ by two substitutions, and every read holding GGGGCC holds GGGG.
        let tags = seqs(&["AAAA", "AATT", "GGGGCC", "GGGG"]);
        for checked in [unset(dist(1)), search(dist(1), MatchMode::All)] {
            let err = validate_group("g", &tags, &checked)
                .unwrap_err()
                .to_string();
            let reason = "tags `GGGGCC` and `GGGG` are indistinguishable";
            assert!(err.contains(reason), "{err}");
        }
        let tags = seqs(&["AAACCCGG", "CCGGGTTA", "ATACCAACGC", "GCGTTGGTAT"]);
        let both = Search {
            both_strands: true,
            ..unset(dist(1))
        };
        let err = validate_group("g", &tags, &both).unwrap_err().to_string();
        let reason = "tags `ATACCAACGC` and `GCGTTGGTAT` are reverse complements of each other";
        assert!(err.contains(reason), "{err}");
    }

    #[test]
    fn test_iupac_tags_that_meet_on_opposite_strands_are_not_called_reverse_complements() {
        // The reverse complement of CCTY is RAGG, and R holds the A of AAGG.
        let both = Search {
            both_strands: true,
            ..unset(dist(0))
        };
        let err = validate_group("g", &seqs(&["AAGG", "CCTY"]), &both)
            .unwrap_err()
            .to_string();
        let reason = "tags `AAGG` and `CCTY` each match the other's reverse complement (the bases `AAGG` match both with zero errors, one on each strand)";
        assert!(err.contains(reason), "{err}");
    }

    /// Plain 16 bp tags at least three substitutions apart, after two that one
    /// substitution each turns into the same read.
    fn after_a_close_pair(count: usize) -> Vec<Vec<u8>> {
        let mut tags = seqs(&["ACGTACGTACGTACGT", "ACGTACGTACGTACCA"]);
        tags.extend(spread_barcodes(count, 16, 3));
        tags
    }

    #[test]
    fn test_a_zero_error_pair_among_thousands_of_plain_tags_is_found() {
        let mut tags = after_a_close_pair(3_000);
        tags.push(tags[1_500].clone());
        let shown = String::from_utf8_lossy(&tags[1_500]).into_owned();
        let err = validate_group("g", &tags, &unset(dist(1)))
            .unwrap_err()
            .to_string();
        let reason = format!("tags `{shown}` and `{shown}` are indistinguishable");
        assert!(err.contains(&reason), "{err}");
        tags.pop();
        let mut flipped = tags[1_500].clone();
        reverse_complement(&mut flipped);
        let flipped_shown = String::from_utf8_lossy(&flipped).into_owned();
        tags.push(flipped);
        let both = Search {
            both_strands: true,
            ..unset(dist(1))
        };
        let err = validate_group("g", &tags, &both).unwrap_err().to_string();
        let reason = format!("tags `{shown}` and `{flipped_shown}` are reverse complements");
        assert!(err.contains(&reason), "{err}");
    }

    #[test]
    fn test_a_zero_error_pair_past_the_pair_limit_gives_the_plain_refusal() {
        // A 12 bp tag at the start of a 16 bp one is the last pair, past the first 2^20.
        let mut tags = after_a_close_pair(1_496);
        let long = Rng(0x7a11).sample(16, b"ACGT");
        tags.extend([long.clone(), long[..12].to_vec()]);
        assert!(tags.len() * (tags.len() - 1) / 2 > 1 << 20);
        let err = validate_group("g", &tags, &unset(dist(1)))
            .unwrap_err()
            .to_string();
        let plain =
            "tags `ACGTACGTACGTACGT` and `ACGTACGTACGTACCA` can both match one read within dist=1";
        assert!(err.contains(plain), "{err}");
    }

    #[test]
    fn test_plain_tags_of_one_length_are_looked_up_without_comparing_pairs() {
        let mut seqs = after_a_close_pair(3_000);
        seqs.push(seqs[1_500].to_ascii_lowercase());
        let tags: Vec<Tag> = seqs.iter().map(|seq| Tag::new(seq)).collect();
        let shape = Shape::of(None);
        let pair = |seqs: &[Vec<u8>], limit| {
            exact_pair(seqs, &tags, &[], &unset(dist(0)), shape, limit).map(|(i, j, _)| (i, j))
        };
        assert_eq!(
            pair(&seqs, 0),
            Some((1_500, 3_002)),
            "found with no pair compared"
        );
        let mut flipped = seqs[1_500].clone();
        reverse_complement(&mut flipped);
        let both = Search {
            both_strands: true,
            ..unset(dist(0))
        };
        let mut seqs = seqs[..3_002].to_vec();
        seqs.push(flipped);
        let flipped_tags: Vec<Tag> = seqs.iter().map(|seq| Tag::new(&flip(seq))).collect();
        let tags: Vec<Tag> = seqs.iter().map(|seq| Tag::new(seq)).collect();
        let found = exact_pair(&seqs, &tags, &flipped_tags, &both, shape, 0);
        assert_eq!(found.map(|(i, j, _)| (i, j)), Some((1_500, 3_002)));
    }

    #[test]
    fn test_other_sets_are_compared_pair_by_pair_up_to_the_limit() {
        // ACGN and ACGT, the first and last tags, form the second pair compared.
        let seqs = seqs(&["ACGN", "TTTT", "ACGT"]);
        let tags: Vec<Tag> = seqs.iter().map(|seq| Tag::new(seq)).collect();
        let pair = |limit| {
            exact_pair(&seqs, &tags, &[], &unset(dist(0)), Shape::of(None), limit)
                .map(|(i, j, _)| (i, j))
        };
        assert_eq!(pair(1), None);
        assert_eq!(pair(2), Some((0, 2)));
    }

    #[test]
    fn test_a_pair_that_cannot_fit_the_window_is_never_taken_to_collide() {
        // A 150 bp tag needs at least 142 bases under eight indels, more than the window holds.
        let mut tags = long_pair();
        tags[1].truncate(120);
        let narrow = Search {
            window: Some(141),
            ..search(budget(0, 8), MatchMode::Nearest)
        };
        assert!(validate_group("g", &tags, &narrow).unwrap().is_empty());
        let narrow = Search {
            window: Some(141),
            ..unset(budget(0, 8))
        };
        assert!(validate_group("g", &tags, &narrow).is_ok());
    }

    #[test]
    fn test_totals_say_may_when_any_counted_pair_was_not_compared_exactly() {
        let mut tags = long_pair();
        tags.push(b"ACGTTGCAACGT".to_vec());
        let w = validate_group("g", &tags, &search(budget(0, 8), MatchMode::All)).unwrap();
        let total = "group `g`: 3 pairs of tags may match one read only through an indel";
        assert!(w.last().is_some_and(|w| w.starts_with(total)), "{w:?}");
    }

    /// Every 4-mer; two of them collide at dist=1 when they differ at two or fewer positions.
    fn every_4mer() -> Vec<Vec<u8>> {
        let base = |code: usize, at: usize| b"ACGT"[code >> (2 * at) & 3];
        (0..256)
            .map(|code| (0..4).map(|at| base(code, at)).collect())
            .collect()
    }

    #[test]
    fn test_nearest_lists_the_first_colliding_pairs_and_counts_the_rest() {
        let tags = seqs(&["AAAAA", "AAAAC", "AAAAG", "AAAAT", "AAACA", "AAAGA"]);
        let w = validate_group("g", &tags, &search(dist(1), MatchMode::Nearest)).unwrap();
        assert_eq!(w.len(), 11, "ten pairs and a count: {w:?}");
        assert!(w[..10]
            .iter()
            .all(|w| w.contains("can both match one read")));
        assert!(
            w[10]
                .contains("group `g`: 5 more pairs of tags can both match one read within dist=1;"),
            "{}",
            w[10]
        );
    }

    #[test]
    fn test_nearest_stops_counting_but_still_refuses_zero_error_pairs() {
        let mut tags = every_4mer();
        let w = validate_group("g", &tags, &search(dist(1), MatchMode::Nearest)).unwrap();
        assert_eq!(w.len(), 11, "ten pairs and a count: {w:?}");
        assert!(w[10].contains("at least 990 more pairs"), "{}", w[10]);
        tags.push(tags[255].clone());
        let err = validate_group("g", &tags, &search(dist(1), MatchMode::Nearest))
            .unwrap_err()
            .to_string();
        assert!(err.contains("indistinguishable"), "{err}");
    }

    #[test]
    fn test_mode_all_points_at_mode_nearest_to_list_the_pairs() {
        let err = check(&["AAAA", "AATT"], dist(1), MatchMode::All)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("or use mode=nearest, which only warns and lists the colliding pairs"),
            "{err}"
        );
    }

    #[test]
    fn test_base_matching_nothing_is_charged_to_its_own_tag() {
        // X matches no read base, so both tags pay for it and no read is within one of both.
        let tags = seqs(&["AXA", "CXA"]);
        let w = validate_group("g", &tags, &search(dist(1), MatchMode::All)).unwrap();
        assert!(w.is_empty(), "{w:?}");
        assert!(validate_group("g", &tags, &search(dist(2), MatchMode::All)).is_err());
    }

    #[test]
    fn test_dist_spec_leaves_out_default_fields() {
        assert_eq!(dist_spec(&dist(2)), "2");
        assert_eq!(dist_spec(&budget(1, 1)), "1:1");
        let capped = Dist {
            total: 1,
            ..budget(1, 1)
        };
        assert_eq!(dist_spec(&capped), "1:1:1");
    }

    /// A random search over short tags, within what the grammar accepts.
    fn random_search(rng: &mut Rng) -> Search {
        let (mismatch, indel, total) = [
            (0, 0, 0),
            (1, 0, 1),
            (2, 0, 2),
            (0, 1, 1),
            (1, 1, 1),
            (1, 1, 2),
            (0, 2, 2),
            (2, 1, 3),
            (1, 2, 3),
        ][rng.below(9)];
        let anchor = [None, Some(Anchor::FivePrime), Some(Anchor::ThreePrime)][rng.below(3)]
            .filter(|&anchor| anchor != Anchor::ThreePrime || indel == 0);
        let bounded = anchor != Some(Anchor::FivePrime) && rng.below(3) == 0;
        Search {
            dist: Dist {
                mismatch,
                indel,
                total,
            },
            mode: Some(MatchMode::All),
            anchor,
            both_strands: anchor.is_none() && rng.below(2) == 0,
            window: bounded.then(|| 1 + rng.below(8)),
        }
    }

    fn flip(seq: &[u8]) -> Vec<u8> {
        let mut seq = seq.to_vec();
        reverse_complement(&mut seq);
        seq
    }

    fn found(a: &[u8], b: &[u8], search: &Search) -> Option<Overlap> {
        let flipped = (Tag::new(&flip(a)), Tag::new(&flip(b)));
        let strands = search.both_strands.then_some((&flipped.0, &flipped.1));
        let shape = Shape::of(search.anchor);
        collision(&Tag::new(a), &Tag::new(b), strands, search, shape)
    }

    #[test]
    fn test_reported_bases_match_both_tags_in_the_matcher() {
        let mut rng = Rng(0x5eed);
        let mut reported = 0;
        for _ in 0..4000 {
            let search = random_search(&mut rng);
            let letters: &[u8] = [&b"ACGT"[..], b"AC", b"ACGTNRY"][rng.below(3)];
            let (len_a, len_b) = (1 + rng.below(8), 1 + rng.below(8));
            let (a, b) = (rng.sample(len_a, letters), rng.sample(len_b, letters));
            let Some(overlap) = found(&a, &b, &search) else {
                continue;
            };
            reported += 1;
            let bases = overlap.bases.expect("short tags are compared exactly");
            let mut group = crate::matcher::CompiledGroup::new("g", vec![a.clone(), b.clone()]);
            group.dist = search.dist;
            group.revcomp = search.both_strands;
            group.anchor = search.anchor;
            let read = [bases.as_slice()];
            let mut scratch = crate::matcher::Scratch::new();
            let outcome = crate::matcher::match_group(&group, &mut scratch, &read);
            for tag_idx in 0..2 {
                assert!(
                    outcome.accepted.iter().any(|m| m.tag_idx == tag_idx),
                    "tag {tag_idx} of {:?} misses {:?} under {search:?}: {outcome:?}",
                    seqs_shown(&[&a, &b]),
                    String::from_utf8_lossy(&bases),
                );
            }
        }
        assert!(reported > 500, "the trials exercise the check: {reported}");
    }

    fn seqs_shown(tags: &[&[u8]]) -> Vec<String> {
        tags.iter()
            .map(|tag| String::from_utf8_lossy(tag).into_owned())
            .collect()
    }

    /// `seq` with up to `edits` random substitutions, insertions, and deletions.
    fn mutate(rng: &mut Rng, seq: &[u8], edits: usize, letters: &[u8]) -> Vec<u8> {
        let mut seq = seq.to_vec();
        for _ in 0..rng.below(edits + 1) {
            let (at, base) = (rng.below(seq.len() + 1), letters[rng.below(letters.len())]);
            match rng.below(3) {
                0 if at < seq.len() => seq[at] = base,
                1 => seq.insert(at, base),
                _ if at < seq.len() && seq.len() > 1 => drop(seq.remove(at)),
                _ => {}
            }
        }
        seq
    }

    /// The fewest disjoint pairs plus unpaired bases over the alignments of
    /// `inner` within `outer`, by a plain table.
    fn plain_edits(outer: &[u8], inner: &[u8], shape: Shape) -> usize {
        let mut row: Vec<usize> = (0..=outer.len())
            .map(|i| if shape.free_start { 0 } else { i })
            .collect();
        for (j, &base) in inner.iter().enumerate() {
            let mut next = vec![j + 1; outer.len() + 1];
            for (i, &other) in outer.iter().enumerate() {
                let missed = usize::from(mask(base) & mask(other) == 0);
                next[i + 1] = (row[i] + missed).min(row[i + 1] + 1).min(next[i] + 1);
            }
            row = next;
        }
        match shape.free_end {
            true => *row.iter().min().unwrap(),
            false => row[outer.len()],
        }
    }

    /// The fewest disjoint pairs plus unpaired bases over the alignments of
    /// `inner` within `outer` with at most `gaps` unpaired bases, by a plain table.
    fn plain_fewest(outer: &[u8], inner: &[u8], shape: Shape, gaps: usize) -> usize {
        let mut fewest: Vec<Vec<usize>> = (0..=outer.len())
            .map(|i| {
                let mut cell = vec![usize::MAX; gaps + 1];
                match shape.free_start {
                    true => cell[0] = 0,
                    false if i <= gaps => cell[i] = i,
                    false => {}
                }
                cell
            })
            .collect();
        for &base in inner {
            let mut next = vec![vec![usize::MAX; gaps + 1]; outer.len() + 1];
            for i in 0..=outer.len() {
                for g in 0..=gaps {
                    let mut best = usize::MAX;
                    if g > 0 {
                        best = best.min(fewest[i][g - 1].saturating_add(1));
                    }
                    if i > 0 {
                        let missed = usize::from(mask(base) & mask(outer[i - 1]) == 0);
                        best = best.min(fewest[i - 1][g].saturating_add(missed));
                        if g > 0 {
                            best = best.min(next[i - 1][g - 1].saturating_add(1));
                        }
                    }
                    next[i][g] = best;
                }
            }
            fewest = next;
        }
        let ends = match shape.free_end {
            true => 0..=outer.len(),
            false => outer.len()..=outer.len(),
        };
        ends.flat_map(|i| fewest[i].clone()).min().unwrap()
    }

    #[test]
    fn test_edit_bounds_agree_with_a_plain_table() {
        let mut rng = Rng(0xed17);
        for trial in 0..1600 {
            let letters: &[u8] = [&b"ACGT"[..], b"AC", b"ACGTNRX"][rng.below(3)];
            let len_x = match trial % 8 {
                0 => 65 + rng.below(70),
                _ => 1 + rng.below(20),
            };
            let x = rng.sample(len_x, letters);
            let len_y = [len_x, 1 + rng.below(len_x + 2)][rng.below(2)];
            let y = match rng.below(2) {
                0 => mutate(&mut rng, &x, 6, letters),
                _ => rng.sample(len_y, letters),
            };
            let shape = Shape {
                free_start: rng.below(2) == 0,
                free_end: rng.below(2) == 0,
            };
            let (outer, inner) = (Tag::new(&x), Tag::new(&y));
            let fewest = plain_edits(&x, &y, shape);
            let by_gaps: Vec<usize> = (0..5).map(|g| plain_fewest(&x, &y, shape, g)).collect();
            for budget in 0..8 {
                let shown = seqs_shown(&[&y, &x]);
                let within = edits_within(&outer, &inner, shape, budget);
                assert_eq!(within, fewest <= budget, "{shown:?} within {budget}");
                if within {
                    assert!(unmatched(&outer, &inner) <= budget, "{shown:?} counts");
                    assert!(shares_grams(&outer, &inner, budget), "{shown:?} 3-grams");
                }
                for (gaps, &capped) in by_gaps.iter().enumerate() {
                    if capped <= budget {
                        let banded = covered(&outer, &inner, shape, gaps, budget);
                        assert!(banded, "{shown:?} with {gaps} gaps within {budget}");
                    }
                }
            }
        }
    }

    #[test]
    fn test_long_tags_far_apart_are_ruled_out_before_the_table() {
        let mut rng = Rng(0x10ad);
        let sliding = Shape::of(None);
        let mut banded = 0;
        for _ in 0..50 {
            let (x, y) = (rng.sample(70, b"ACGT"), rng.sample(70, b"ACGT"));
            let (outer, inner) = (Tag::new(&x), Tag::new(&y));
            let shown = seqs_shown(&[&y, &x]);
            assert!(!edits_within(&outer, &inner, sliding, 8), "{shown:?}");
            banded += usize::from(!covered(&outer, &inner, sliding, 4, 8));
        }
        assert!(banded >= 40, "the diagonal band rules out {banded} of 50");
    }

    #[test]
    fn test_table_for_an_indel_budget_past_the_tags_stays_small() {
        let (outer, inner) = (Tag::new(b"ACGTACGT"), Tag::new(b"TGCATGCA"));
        let huge = Dist {
            mismatch: 0,
            indel: 40,
            total: 40,
        };
        let shape = Shape::of(None);
        let grid = Grid::fill(&outer, &inner, &huge, shape);
        assert!(grid.cells.len() <= 9 * 17 * 9 * 9, "{}", grid.cells.len());
    }

    #[test]
    fn test_a_table_past_the_cell_cap_is_never_filled() {
        // One base in fifty differs, so the filters pass and only the table could tell.
        let tag = Rng(0x500).sample(500, b"ACGT");
        let other: Vec<u8> = (tag.iter().enumerate())
            .map(|(at, &base)| match at % 50 {
                25 if base == b'A' => b'C',
                25 => b'A',
                _ => base,
            })
            .collect();
        let huge = budget(0, 50);
        let (outer, inner) = (Tag::new(&tag), Tag::new(&other));
        let shape = Shape::of(None);
        let size = Grid::size(&outer, &inner, &huge, shape);
        assert!(size > 500 * TABLE_CELLS, "{size}");
        let checked = search(huge, MatchMode::Nearest);
        CELLS_LEFT.set(GROUP_CELLS);
        let found = collision(&outer, &inner, None, &checked, shape).unwrap();
        assert!(found.bases.is_none());
        assert_eq!(CELLS_LEFT.get(), GROUP_CELLS, "no table was filled");
        let tags = vec![tag, other];
        let w = validate_group("g", &tags, &checked).unwrap();
        assert!(w.len() == 1 && w[0].contains(TOO_COSTLY), "{w:?}");
    }

    #[test]
    fn test_pairs_past_the_group_cell_budget_are_taken_to_collide() {
        // Deleting the last base of either tag leaves ACGTA, one indel from both.
        let (a, b) = (Tag::new(b"ACGTAC"), Tag::new(b"ACGTAG"));
        let (checked, shape) = (search(budget(0, 1), MatchMode::Nearest), Shape::of(None));
        let size = Grid::size(&a, &b, &checked.dist, shape);
        CELLS_LEFT.set(size - 1);
        let found = collision(&a, &b, None, &checked, shape).unwrap();
        assert!(
            found.bases.is_none(),
            "taken to collide once the budget runs short"
        );
        assert_eq!(CELLS_LEFT.get(), size - 1, "no table was filled");
        CELLS_LEFT.set(4 * size);
        let found = collision(&a, &b, None, &checked, shape).unwrap();
        assert_eq!(found.bases.as_deref(), Some(&b"ACGTA"[..]));
        assert_eq!(CELLS_LEFT.get(), 2 * size, "both orders drew a table");
        let tags = seqs(&["ACGTAC", "ACGTAG"]);
        assert!(validate_group("g", &tags, &checked).unwrap()[0].contains("can both match"));
        assert_eq!(
            CELLS_LEFT.get(),
            GROUP_CELLS - 2 * size,
            "each group starts afresh"
        );
    }

    /// Whether bases exist, no longer than the search window, that `outer`
    /// matches end to end and `inner` matches within as the search places it.
    fn bases_exist(outer: &[u8], inner: &[u8], search: &Search) -> bool {
        let (outer, inner) = (Tag::new(outer), Tag::new(inner));
        let reach = outer.len() + search.dist.indel;
        let longest = search.window.map_or(reach, |w| w.min(reach));
        let shape = Shape::of(search.anchor);
        (1..=longest).any(|len| {
            (0..1usize << (2 * len)).any(|code| {
                let bases: Vec<u8> = (0..len).map(|at| b"ACGT"[code >> (2 * at) & 3]).collect();
                aligns(&outer.masks, &bases, &search.dist, Shape::PINNED)
                    && aligns(&inner.masks, &bases, &search.dist, shape)
            })
        })
    }

    #[test]
    fn test_bases_matching_both_tags_are_found_whenever_they_exist() {
        let mut rng = Rng(0xb0a7);
        let mut collided = 0;
        for _ in 0..600 {
            let search = random_search(&mut rng);
            let letters: &[u8] = [&b"ACGT"[..], b"AC", b"ACGTNR"][rng.below(3)];
            let (len_a, len_b) = (1 + rng.below(4), 1 + rng.below(4));
            let (a, b) = (rng.sample(len_a, letters), rng.sample(len_b, letters));
            let mut orders = vec![(a.clone(), b.clone()), (b.clone(), a.clone())];
            if search.both_strands {
                orders.extend([(a.clone(), flip(&b)), (b.clone(), flip(&a))]);
            }
            let exists = orders
                .iter()
                .any(|(outer, inner)| bases_exist(outer, inner, &search));
            collided += usize::from(exists);
            assert_eq!(
                found(&a, &b, &search).is_some(),
                exists,
                "{:?} under {search:?}",
                seqs_shown(&[&a, &b])
            );
        }
        assert!(collided > 100, "the trials exercise collisions: {collided}");
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
            &search(dist(1), MatchMode::All),
        )
        .unwrap();
        assert!(w.is_empty(), "well-separated barcodes are clean: {w:?}");
    }

    #[test]
    fn test_indistinguishable_pair_fails_regardless_of_mode() {
        // AAAW (W = A|T) covers AAAA at every position -> a read AAAA matches
        // both at dist 0.
        for mode in [MatchMode::All, MatchMode::Nearest] {
            let err = validate_group("g", &seqs(&["AAAA", "AAAW"]), &search(dist(0), mode))
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
        let err = validate_group("g", &tags, &search(dist(1), MatchMode::All))
            .unwrap_err()
            .to_string();
        assert!(err.contains("mode=all"), "{err}");
        // mode=nearest: warning, not error.
        let w = validate_group("g", &tags, &search(dist(1), MatchMode::Nearest)).unwrap();
        assert_eq!(w.len(), 1, "one collision warning: {w:?}");
        assert!(w[0].contains("within dist=1"));
    }

    #[test]
    fn test_well_separated_pair_ok_at_dist() {
        // AAAA vs TTTT differ at all 4 positions; at dist=1, D=4 > 2 -> no
        // collision.
        assert!(validate_group(
            "g",
            &seqs(&["AAAA", "TTTT"]),
            &search(dist(1), MatchMode::All)
        )
        .unwrap()
        .is_empty());
    }

    #[test]
    fn test_explosion_hard_cap_errors() {
        // 10 fully-degenerate positions plus a 3-way `B` = 4^10 * 3 ~= 3.1M,
        // over the 1<<20 cap. The `B` keeps it from being all-N (which has its
        // own, earlier error), so this exercises the cap.
        let err = validate_group(
            "g",
            &seqs(&["NNNNNNNNNNB"]),
            &search(dist(0), MatchMode::All),
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("cap"), "{err}");
    }

    #[test]
    fn test_all_n_tag_is_rejected_with_slicing_hint() {
        // A short all-N tag (NNNN = 256, under the degeneracy cap) is still
        // rejected: it matches any 4 bases, so it is a positional span, not a
        // barcode. This is the gap the cap alone misses.
        let err = validate_group("bc", &seqs(&["NNNN"]), &search(dist(0), MatchMode::All))
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
        assert!(validate_group("bc", &seqs(&["nnnn"]), &search(dist(0), MatchMode::All)).is_err());
    }

    #[test]
    fn test_partially_degenerate_tag_is_not_all_n() {
        // A tag with any non-N base is not all-N; NNNANNN (4^6 = 4096, at the
        // warn threshold, not over it) is accepted, so the all-N guard does not
        // fire on it.
        assert!(
            validate_group("bc", &seqs(&["NNNANNN"]), &search(dist(0), MatchMode::All)).is_ok()
        );
    }

    #[test]
    fn test_mask_covers_every_iupac_code() {
        // Concrete bases occupy one bit, U aliases T, each ambiguity code is the
        // union of its member bits, and N is all four. The mask upper-cases first,
        // so case is irrelevant; any byte that is not a code is the empty set.
        let expected: &[(u8, u8)] = &[
            (b'A', 0b0001),
            (b'C', 0b0010),
            (b'G', 0b0100),
            (b'T', 0b1000),
            (b'U', 0b1000),
            (b'R', 0b0101), // A G
            (b'Y', 0b1010), // C T
            (b'S', 0b0110), // C G
            (b'W', 0b1001), // A T
            (b'K', 0b1100), // G T
            (b'M', 0b0011), // A C
            (b'B', 0b1110), // C G T
            (b'D', 0b1101), // A G T
            (b'H', 0b1011), // A C T
            (b'V', 0b0111), // A C G
            (b'N', 0b1111),
        ];
        for &(code, bits) in expected {
            assert_eq!(mask(code), bits, "mask({})", code as char);
            assert_eq!(
                mask(code.to_ascii_lowercase()),
                bits,
                "mask is case-insensitive for {}",
                code as char
            );
        }
        for &junk in b"X-. 0" {
            assert_eq!(
                mask(junk),
                0,
                "non-code {:?} maps to the empty set",
                junk as char
            );
        }
    }

    #[test]
    fn test_complement_covers_every_iupac_code() {
        // Watson-Crick on ACGT, U complements like T, and each ambiguity code maps
        // to the code of its complemented member set (e.g. R={A,G} -> {T,C}=Y).
        // Case is preserved; unrecognized bytes pass through unchanged.
        let pairs: &[(u8, u8)] = &[
            (b'A', b'T'),
            (b'C', b'G'),
            (b'G', b'C'),
            (b'T', b'A'),
            (b'U', b'A'),
            (b'R', b'Y'),
            (b'Y', b'R'),
            (b'S', b'S'),
            (b'W', b'W'),
            (b'K', b'M'),
            (b'M', b'K'),
            (b'B', b'V'),
            (b'V', b'B'),
            (b'D', b'H'),
            (b'H', b'D'),
            (b'N', b'N'),
        ];
        for &(base, comp) in pairs {
            assert_eq!(complement(base), comp, "complement({})", base as char);
            assert_eq!(
                complement(base.to_ascii_lowercase()),
                comp.to_ascii_lowercase(),
                "complement preserves case for {}",
                base as char
            );
        }
        for &junk in b"X-* " {
            assert_eq!(
                complement(junk),
                junk,
                "unknown byte {:?} passes through",
                junk as char
            );
        }
    }

    #[test]
    fn test_reverse_complement_mixed_case_and_ambiguity() {
        // Reverse the bytes, then complement each, preserving case. AACG ->
        // (reverse) GCAA -> (complement) CGTT.
        let mut seq = b"AACG".to_vec();
        reverse_complement(&mut seq);
        assert_eq!(&seq, b"CGTT");
        // Lowercase plus ambiguity codes: rySw -> (reverse) wSyr -> (complement)
        // wSry (w<->w, S<->S self-complement; y->r, r->y).
        let mut mixed = b"rySw".to_vec();
        reverse_complement(&mut mixed);
        assert_eq!(&mixed, b"wSry");
    }

    #[test]
    fn test_single_tag_degeneracy_warns_under_cap() {
        // 4^7 = 16384 expansions: past the per-tag warn threshold (4096) but well
        // under the hard cap, and not all-N (the trailing A), so it is accepted
        // with a soft warning rather than rejected.
        let w =
            validate_group("g", &seqs(&["NNNNNNNA"]), &search(dist(0), MatchMode::All)).unwrap();
        assert_eq!(w.len(), 1, "one per-tag degeneracy warning: {w:?}");
        assert!(w[0].contains("highly degenerate"), "{}", w[0]);
    }

    #[test]
    fn test_group_total_degeneracy_warns_under_cap() {
        // Two tags of 4^8 = 65536 sum to 131072: each is under the per-tag cap, but
        // the group total clears the group warn threshold (100_000) while staying
        // under the group cap. They differ only at the last position, so at dist=0
        // they do not collide.
        let w = validate_group(
            "g",
            &seqs(&["NNNNNNNNA", "NNNNNNNNC"]),
            &search(dist(0), MatchMode::All),
        )
        .unwrap();
        assert!(
            w.iter().any(|m| m.contains("sequences total")),
            "a group-total warning is present: {w:?}"
        );
    }

    #[test]
    fn test_group_total_hard_cap_errors() {
        // Ten tags of 4^10 = 1_048_576 each (exactly the per-tag cap, so none trips
        // it individually) sum to 10_485_760, over the 10_000_000 group cap.
        let big = vec![b"NNNNNNNNNNA".to_vec(); 10];
        let err = validate_group("g", &big, &search(dist(0), MatchMode::All))
            .unwrap_err()
            .to_string();
        assert!(err.contains("sequences total"), "{err}");
        assert!(err.contains("cap"), "{err}");
    }
}
