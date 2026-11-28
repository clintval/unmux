//! Approximate matching of tag groups against records, built on the `sassy`
//! SIMD matcher with the `Iupac` profile, so IUPAC ambiguity codes in a TAG
//! (pattern side) match natively at zero cost. A record U matches as T, and any
//! other non-ACGT record base (an `N` no-call or a record-side degenerate) is
//! instead masked to a never-match sentinel so it is charged against the `dist`
//! budget rather than matching any tag base for free (splitcode semantics; an
//! all-`N` record falls to unmatched, not mis-assigned).
//!
//! This is the record-window search primitive: given a [`CompiledGroup`] (its
//! tag sequences plus the `loc`/`dist`/`bothStrands`/`mode`/`delta`
//! attributes) and a record's per-file base segments, it returns every
//! accepted match within
//! the error budget and the single best match under a deterministic total order
//! (so `@grp` anchoring is reproducible). The `dist=mismatch:indel:total`
//! budget is enforced by passing `total` as sassy's `k` and then post-filtering
//! each alignment's CIGAR so the substitution and indel counts stay within
//! their own caps. The search reports one alignment per local minimum of the
//! edit cost over the end positions, so wherever that filter drops an alignment,
//! the tag is searched again for alignments within the caps.
//!
//! The sequential orchestration (`next`/`prev`), `match=` over joined streams,
//! `partial5`/`partial3` overhang, and the canonicalization (corrected-form)
//! layer on top of this primitive: the engine runs groups in dependency order
//! ([`CompiledGroup::prereq`]) and calls [`match_verdict_after`] when a
//! relative `next` link locates a downstream group past an upstream match.

use std::borrow::Cow;
use std::cell::LazyCell;
use std::cmp::{Ordering, Reverse};

use ahash::{AHashMap, AHashSet};
use pa_types::CigarOp;
use sassy::profiles::Iupac;
use sassy::{EncodedPatterns, Match as SassyMatch, Searcher, Strand};
use smallvec::SmallVec;

use crate::extract::MatchSpan;
use crate::grammar::{Anchor, Dist, Endpoint, Location, MatchMode, NextLink, Partial};

/// The smallest equal-length tag bucket that uses the `sassy` v2 batched search
/// path; smaller buckets stay on the per-tag `search()`.
const BATCHED_MIN_TAGS: usize = 9;
/// The longest pattern the batched path encodes (one pattern per SIMD lane);
/// longer tags stay per-tag.
const BATCHED_MAX_LEN: usize = 64;

/// A precomputed batch of equal-length tags for the `sassy` v2 batched search
/// path: `encode_patterns` is run once (here), and `search_encoded_patterns`
/// distributes the batch across SIMD lanes per record. `tag_indices` maps a
/// match's `pattern_idx` back to the group's tag index.
#[derive(Debug, Clone)]
pub struct EncodedBatch {
    /// The group-local tag indices in this batch, in `pattern_idx` order.
    pub tag_indices: Vec<usize>,
    /// The encoded patterns (one per tag, all the same length).
    pub encoded: EncodedPatterns<Iupac>,
}

/// A group's sequential prerequisites: every upstream group that must match
/// before this one (the group is skipped if any did not match), and an optional
/// relative window from a `next=GROUP:lo-hi` link that replaces this group's
/// `loc`. Both a `next=` link targeting this group and this group's own `prev=`
/// contribute upstreams, so neither prerequisite is ever dropped; a bare
/// `next`/`prev` link carries no window (the group searches at its own `loc`).
///
/// Upstreams are stored as canonical group indices (each group's position in
/// the engine's group slice), resolved once at compile time, so the per-record
/// matching path indexes the positional hits directly rather than re-keying a
/// map by group name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Prereq {
    /// Every upstream group that must have matched first, by canonical index
    /// (deduped).
    pub upstreams: Vec<usize>,
    /// The single relative window (the grammar rejects two `next=` links to the
    /// same group), or `None` for a bare/`prev` prerequisite.
    pub window: Option<RelativeWindow>,
}

/// The relative search window a `next=GROUP:lo-hi` link gives its downstream
/// group, measured past the upstream group's match in that match's
/// orientation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RelativeWindow {
    /// The upstream group whose match the window is measured from, by
    /// canonical index.
    pub upstream: usize,
    /// Bases past the upstream match where the window starts.
    pub lo: usize,
    /// Bases past the upstream match where the window ends.
    pub hi: usize,
}

/// The strand a tag matched on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchStrand {
    /// Forward strand.
    Forward,
    /// Reverse complement of the record.
    ReverseComplement,
}

/// One accepted match of a group's tag within a record. Coordinates are
/// absolute offsets into the matched input file's bases (the search-window
/// offset is already folded in).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TagMatch {
    /// Index of the matched tag within the group (definition order).
    pub tag_idx: usize,
    /// The 0-based input file the match lies in.
    pub file: usize,
    /// 0-based start offset within that file's bases.
    pub start: usize,
    /// 0-based exclusive end offset within that file's bases.
    pub end: usize,
    /// Total edit cost of the alignment.
    pub cost: usize,
    /// Substitutions in the alignment.
    pub subs: usize,
    /// Insertions plus deletions in the alignment.
    pub indels: usize,
    /// The strand the tag matched on.
    pub strand: MatchStrand,
}

impl TagMatch {
    /// This match as seen on the reverse complement of a text `len` bases long.
    fn mirrored(self, len: usize) -> Self {
        let strand = match self.strand {
            MatchStrand::Forward => MatchStrand::ReverseComplement,
            MatchStrand::ReverseComplement => MatchStrand::Forward,
        };
        Self {
            start: len - self.end,
            end: len - self.start,
            strand,
            ..self
        }
    }
}

/// A group compiled for matching: the tag sequences plus the attributes
/// governing a record-window search. Built from a parsed group and its loaded
/// tag set; held independently of the parser so it can be unit-tested with
/// inline tags.
#[derive(Debug, Clone)]
pub struct CompiledGroup {
    /// The group name (for diagnostics).
    pub name: String,
    /// Tag sequences in definition order (bases / IUPAC, uppercase).
    pub tags: Vec<Vec<u8>>,
    /// The tags as `sassy` searches them, with U read as T, since it complements no U.
    pub(crate) patterns: Vec<Vec<u8>>,
    /// The error budget.
    pub dist: Dist,
    /// The optional search window; `None` searches every file's whole record.
    pub loc: Option<Location>,
    /// Search the reverse-complement strand too (`bothStrands=`).
    pub revcomp: bool,
    /// Match acceptance policy.
    pub mode: MatchMode,
    /// Minimum best-vs-runner-up cost gap (only meaningful with
    /// `mode=nearest`).
    pub delta: usize,
    /// `minFindsPerTag`: the minimum occurrences of a tag that does occur.
    pub min_finds_per_tag: Option<usize>,
    /// `maxFindsPerTag`: the maximum occurrences of any single tag.
    pub max_finds_per_tag: Option<usize>,
    /// `minFindsPerGroup`: the minimum total matches across the group.
    pub min_finds_per_group: Option<usize>,
    /// `maxFindsPerGroup`: the maximum total matches across the group.
    pub max_finds_per_group: Option<usize>,
    /// The sequential prerequisite (from a `next`/`prev` link), or `None` for
    /// an independent group.
    pub prereq: Option<Prereq>,
    /// `match=`: the `--extract` stream names whose joined value this group
    /// matches against, instead of a record window. `None` is a normal
    /// record-window search.
    pub match_streams: Option<Vec<String>>,
    /// `anchor=`: when set, pins each tag's 5'/3' edge at the `loc` boundary
    /// and matches over exactly its own length by a direct IUPAC Hamming
    /// compare (no `sassy` search); `None` is the sliding window.
    pub anchor: Option<Anchor>,
    /// `partial5=`: allow a tag truncated at the record's 5' end (sassy
    /// overhang).
    pub partial5: Option<Partial>,
    /// `partial3=`: allow a tag truncated at the record's 3' end (sassy
    /// overhang).
    pub partial3: Option<Partial>,
    /// Precomputed equal-length tag batches for the `sassy` v2 batched path,
    /// populated by [`CompiledGroup::encode`]. Empty until then; tags not in
    /// any batch use per-tag `search()`.
    pub encoded: Vec<EncodedBatch>,
    /// The group-local tag indices covered by `encoded` (searched as a batch;
    /// the rest fall back to per-tag `search()`), precomputed by
    /// [`CompiledGroup::encode`] alongside `encoded` so the matcher reads it
    /// per search instead of rebuilding the set for every record.
    pub batched_tags: AHashSet<usize>,
    /// For an `@RG` group: RG id -> tag index, used to route a record by its
    /// `RG:Z` instead of a sequence search. `None` for a normal group.
    pub read_group: Option<std::collections::HashMap<Vec<u8>, usize>>,
}

/// How a group outcome fails its find bounds, for the optional `--qc-tag`
/// unassigned slug and the unassigned metrics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FindFailure {
    /// A find-count constraint was violated.
    Bound {
        /// The constraint that was violated.
        kind: FindKind,
        /// The observed count: group matches for a per-group bound, or
        /// occurrences of the offending tag for a per-tag bound.
        observed: usize,
        /// The configured bound that was crossed.
        bound: usize,
    },
    /// Under `mode=nearest` no match beat the runner-up by `delta`, which fails
    /// every find bound of the group.
    NearestTie,
}

/// The kind of find-count constraint, for [`FindFailure`]; [`FindKind::as_str`]
/// is the attribute name used in the slug.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum FindKind {
    /// `minFindsPerGroup`: too few of the group's tags matched the record.
    MinPerGroup,
    /// `maxFindsPerGroup`: too many of the group's tags matched the record.
    MaxPerGroup,
    /// `minFindsPerTag`: a single tag matched fewer times than required.
    MinPerTag,
    /// `maxFindsPerTag`: a single tag matched more times than allowed.
    MaxPerTag,
}

impl FindKind {
    /// The snake_case attribute name, matching the `--group` find-count
    /// attribute that was violated.
    pub fn as_str(self) -> &'static str {
        match self {
            FindKind::MinPerGroup => "min_finds_per_group",
            FindKind::MaxPerGroup => "max_finds_per_group",
            FindKind::MinPerTag => "min_finds_per_tag",
            FindKind::MaxPerTag => "max_finds_per_tag",
        }
    }
}

impl CompiledGroup {
    /// A bare group: exact, forward-strand, whole-record, `mode=all`, a `delta`
    /// of 1, no find bounds. Tests and callers set the attributes they care
    /// about on top of this.
    pub fn new(name: impl Into<String>, tags: Vec<Vec<u8>>) -> Self {
        Self {
            name: name.into(),
            patterns: search_patterns(&tags),
            tags,
            dist: Dist::default(),
            loc: None,
            revcomp: false,
            mode: MatchMode::All,
            delta: 1,
            min_finds_per_tag: None,
            max_finds_per_tag: None,
            min_finds_per_group: None,
            max_finds_per_group: None,
            prereq: None,
            match_streams: None,
            partial5: None,
            partial3: None,
            anchor: None,
            encoded: Vec::new(),
            batched_tags: AHashSet::new(),
            read_group: None,
        }
    }

    /// Precompute the `sassy` v2 batched search batches: group the tags by
    /// length and, for each length bucket of at least [`BATCHED_MIN_TAGS`] tags
    /// up to [`BATCHED_MAX_LEN`] bp, `encode_patterns` once for reuse across
    /// every record. Idempotent: it rebuilds the search patterns and
    /// [`CompiledGroup::encoded`] from the current tags.
    ///
    /// Forward-strand groups only. Two cases skip batching and stay on per-tag
    /// `search()`: partial (overhang) groups, because the batched path is
    /// non-overhang; and reverse-complement groups, because `sassy` 0.2.3
    /// computes the rc strand differently between its two entry points.
    /// `search()` reverse-complements the *text* and remaps coordinates (dedup
    /// in the reversed-text frame), whereas `encode_patterns(include_rc)` +
    /// `search_encoded_patterns` encodes each *pattern*'s reverse complement
    /// against the forward text (dedup in the forward frame). The two frames
    /// produce different rightmost-local-minimum sets, so an rc batch can drop,
    /// add, or shift a match relative to per-tag search and flip the demux
    /// outcome. Forward-strand batching is bit-identical to per-tag, so the
    /// speedup is preserved for the common case; rc stays per-tag until the two
    /// `sassy` paths agree. A small or over-long bucket also stays per-tag.
    pub fn encode(&mut self) {
        self.patterns = search_patterns(&self.tags);
        self.encoded.clear();
        self.batched_tags.clear();
        // TODO(upstream sassy): file an issue against `sassy` (0.2.3) that
        // `search_encoded_patterns` and `search` disagree on the
        // reverse-complement strand: the batched (encoded) path tiles each
        // pattern's reverse complement over the forward text and dedups local
        // minima in the forward frame, while `search` reverse-complements the
        // text and dedups in the reversed frame, so the two report different rc
        // match sets (drop/add/shifted matches). A 60k-trial differential fuzz
        // found 0 forward divergences and thousands of rc ones. The
        // `self.revcomp` short-circuit below is the workaround (rc groups fall
        // back to the verified per-tag path); REMOVE IT once a fixed sassy
        // makes the two paths agree on rc.
        if self.partial5.is_some() || self.partial3.is_some() || self.revcomp {
            return;
        }
        let mut by_len: AHashMap<usize, Vec<usize>> = AHashMap::new();
        for (index, tag) in self.tags.iter().enumerate() {
            if !tag.is_empty() && tag.len() <= BATCHED_MAX_LEN {
                by_len.entry(tag.len()).or_default().push(index);
            }
        }
        // A throwaway forward searcher encodes the patterns; the encoding is
        // reused with each read's per-call searcher.
        let mut searcher = Searcher::<Iupac>::new_fwd();
        for tag_indices in by_len.into_values() {
            if tag_indices.len() < BATCHED_MIN_TAGS {
                continue;
            }
            let patterns: Vec<Vec<u8>> = tag_indices
                .iter()
                .map(|&i| self.patterns[i].clone())
                .collect();
            let encoded = searcher.encode_patterns(&patterns);
            self.encoded.push(EncodedBatch {
                tag_indices,
                encoded,
            });
        }
        self.batched_tags = self
            .encoded
            .iter()
            .flat_map(|batch| batch.tag_indices.iter().copied())
            .collect();
    }

    /// Whether a group outcome satisfies the per-tag and per-group find-count
    /// bounds. An unset bound is unconstrained. Per-group bounds gate the total
    /// accepted-match count; per-tag bounds apply to each tag that occurs, by
    /// its occurrence count. Under `mode=nearest` only the kept match counts.
    pub fn satisfies_finds(&self, outcome: &GroupOutcome) -> bool {
        self.find_failure(outcome).is_none()
    }

    /// The find-count constraint this outcome violates (with the observed count
    /// and the bound), or `None` when every bound is satisfied. Per-group
    /// bounds gate the total accepted-match count; per-tag bounds apply to each
    /// tag that occurs, by its occurrence count. Under `mode=nearest` the only
    /// find is the match it keeps, and a tie fails every bound the group has.
    /// The per-tag check runs in ascending tag-index order so the reported
    /// failure is deterministic when several tags are out of bounds (HashMap
    /// iteration order is not).
    pub fn find_failure(&self, outcome: &GroupOutcome) -> Option<FindFailure> {
        let finds = match (&outcome.best, self.mode) {
            (Some(best), MatchMode::Nearest) => std::slice::from_ref(best),
            (None, MatchMode::Nearest) if !outcome.accepted.is_empty() => {
                return self.has_find_bounds().then_some(FindFailure::NearestTie);
            }
            _ => outcome.accepted.as_slice(),
        };
        let group_count = finds.len();
        if let Some(min) = self.min_finds_per_group {
            if group_count < min {
                return Some(FindFailure::Bound {
                    kind: FindKind::MinPerGroup,
                    observed: group_count,
                    bound: min,
                });
            }
        }
        if let Some(max) = self.max_finds_per_group {
            if group_count > max {
                return Some(FindFailure::Bound {
                    kind: FindKind::MaxPerGroup,
                    observed: group_count,
                    bound: max,
                });
            }
        }
        if self.min_finds_per_tag.is_some() || self.max_finds_per_tag.is_some() {
            let mut counts: AHashMap<usize, usize> = AHashMap::new();
            for m in finds {
                *counts.entry(m.tag_idx).or_default() += 1;
            }
            let mut counts: Vec<(usize, usize)> = counts.into_iter().collect();
            counts.sort_unstable();
            for (_tag_idx, count) in counts {
                if let Some(max) = self.max_finds_per_tag {
                    if count > max {
                        return Some(FindFailure::Bound {
                            kind: FindKind::MaxPerTag,
                            observed: count,
                            bound: max,
                        });
                    }
                }
                if let Some(min) = self.min_finds_per_tag {
                    if count < min {
                        return Some(FindFailure::Bound {
                            kind: FindKind::MinPerTag,
                            observed: count,
                            bound: min,
                        });
                    }
                }
            }
        }
        None
    }

    /// Whether a find bound can fail the group: a minimum of at least one, or any maximum.
    fn has_find_bounds(&self) -> bool {
        self.min_finds_per_group.is_some_and(|min| min >= 1)
            || self.max_finds_per_group.is_some()
            || self.min_finds_per_tag.is_some_and(|min| min >= 1)
            || self.max_finds_per_tag.is_some()
    }
}

/// The downstream search window implied by a `next=GROUP:lo-hi` link, given the
/// chosen upstream match span on the forward strand: `[matchEnd+lo,
/// matchEnd+hi]` in the upstream match's file. Returns `None` when the link has
/// no relative window (a bare `next=` leaves the downstream group searching at
/// its own `loc`). [`match_verdict_after`] also handles a reverse-strand match.
pub fn next_window(upstream: &MatchSpan, link: &NextLink) -> Option<Location> {
    let (lo, hi) = link.window?;
    Some(window_after(upstream, lo, hi))
}

/// The window `[upstream.end + lo, upstream.end + hi]` in the upstream match's
/// file, the search slice a relative `next=GROUP:lo-hi` link gives its
/// downstream group past a forward-strand match. Only the upstream's `file` and
/// `end` are needed, so it takes a [`MatchSpan`] (the positional anchor) rather
/// than a full [`TagMatch`].
pub fn window_after(upstream: &MatchSpan, lo: usize, hi: usize) -> Location {
    Location {
        file: upstream.file,
        start: Endpoint::FromStart(upstream.end + lo),
        end: Some(Endpoint::FromStart(upstream.end + hi)),
    }
}

/// The outcome of matching a group against one record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroupOutcome {
    /// Every match accepted within the budget, across all tags, files, and
    /// strands.
    pub accepted: Vec<TagMatch>,
    /// The single best match (by the total order), or `None` when the group is
    /// unmatched: no candidate, or `mode=nearest` and the best did not beat the
    /// runner-up by `delta`.
    pub best: Option<TagMatch>,
}

impl GroupOutcome {
    /// An unmatched outcome (no accepted matches). Used when a group is skipped
    /// because its sequential prerequisite did not match, or a `match=` stream
    /// was absent.
    pub fn unmatched() -> Self {
        Self {
            accepted: Vec::new(),
            best: None,
        }
    }
}

/// What routing and extraction read from a group's outcome on one record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GroupVerdict {
    /// The single best match, or `None` when the group is unmatched.
    pub best: Option<TagMatch>,
    /// The find-count constraint the outcome violates, if any.
    pub failure: Option<FindFailure>,
}

impl GroupVerdict {
    /// The verdict of `outcome` under the find-count bounds of `group`.
    pub fn of(group: &CompiledGroup, outcome: &GroupOutcome) -> Self {
        Self {
            best: outcome.best,
            failure: group.find_failure(outcome),
        }
    }

    /// The verdict of a group that is skipped, as if it had no match.
    pub fn unmatched(group: &CompiledGroup) -> Self {
        Self::of(group, &GroupOutcome::unmatched())
    }
}

/// Per-worker reusable matcher scratch: one `sassy` [`Searcher`] per `(revcomp,
/// partial)` variant a group can need, built on first use and reused across
/// every record the worker handles. A `Searcher` owns mutable internal buffers
/// (`hp`/`hm`/`cost_matrices`/`lanes`/`matches`) and its search methods take
/// `&mut self`, so it must be thread-local: build exactly one `Scratch` per
/// worker, never construct a `Searcher` per record, and never share one
/// across threads. The group's precomputed `encoded` patterns stay read-only
/// shared and are not part of `Scratch`.
#[derive(Default)]
pub struct Scratch {
    fwd: Option<Searcher<Iupac>>,
    rc: Option<Searcher<Iupac>>,
    fwd_overhang: Option<Searcher<Iupac>>,
    rc_overhang: Option<Searcher<Iupac>>,
    /// Per-worker verdict cache, group name -> (window key -> verdict). The
    /// same barcode window recurs across millions of records, so a hit turns
    /// the per-record `sassy` search into a hash lookup. The key holds the
    /// bytes the kernel reads, not where they are, so windows at different
    /// offsets share an entry and the cached match is kept relative to its
    /// window. Bounded per group ([`MATCH_CACHE_CAP`]).
    cache: AHashMap<String, GroupCache>,
    /// Reused scratch for building a cache lookup key, so a cache hit (the
    /// common case) does not allocate a fresh key per record.
    key_buf: WindowKey,
}

/// A window key, inline for windows of up to 20 bases.
type WindowKey = SmallVec<[u8; 24]>;

/// One group's cached verdicts and the [`key_reach`] that their keys cover.
struct GroupCache {
    reach: usize,
    verdicts: AHashMap<WindowKey, CachedVerdict>,
}

/// A cached verdict, with its best match relative to the window start.
#[derive(Clone, Copy)]
struct CachedVerdict {
    best: Option<RelativeMatch>,
    failure: Option<FindFailure>,
}

/// A best match with offsets from its window start.
#[derive(Clone, Copy)]
struct RelativeMatch {
    tag_idx: u32,
    start: u32,
    end: u32,
    cost: u32,
    subs: u32,
    indels: u32,
    strand: MatchStrand,
}

impl CachedVerdict {
    /// The verdict with its best match made relative to `win_start`, unless an
    /// offset or count does not fit.
    fn relative(verdict: &GroupVerdict, win_start: usize) -> Option<Self> {
        let narrow = |value: usize| u32::try_from(value).ok();
        let best = match &verdict.best {
            None => None,
            Some(m) => Some(RelativeMatch {
                tag_idx: narrow(m.tag_idx)?,
                start: narrow(m.start.checked_sub(win_start)?)?,
                end: narrow(m.end.checked_sub(win_start)?)?,
                cost: narrow(m.cost)?,
                subs: narrow(m.subs)?,
                indels: narrow(m.indels)?,
                strand: m.strand,
            }),
        };
        Some(Self {
            best,
            failure: verdict.failure,
        })
    }

    /// The verdict for a window that starts at `win_start` in input `file`.
    fn placed(&self, file: usize, win_start: usize) -> GroupVerdict {
        GroupVerdict {
            best: self.best.map(|m| TagMatch {
                tag_idx: m.tag_idx as usize,
                file,
                start: win_start + m.start as usize,
                end: win_start + m.end as usize,
                cost: m.cost as usize,
                subs: m.subs as usize,
                indels: m.indels as usize,
                strand: m.strand,
            }),
            failure: self.failure,
        }
    }
}

/// Per-(group, worker) cache entry cap, about 1.6 MB. Barcode windows recur
/// with low cardinality (tags plus a few error variants), so this is not
/// reached for them; past it a group searches each new window uncached.
const MATCH_CACHE_CAP: usize = 1 << 14;

impl Scratch {
    /// A fresh scratch with no searchers yet built; build one per worker
    /// thread.
    pub fn new() -> Self {
        Self::default()
    }

    /// The searcher for a group's `(revcomp, partial)` configuration, built
    /// once and reused. A `partial` group searches with sassy's overhang mode
    /// (a tag may run off the record's 5'/3' end, the truncated bases charged
    /// at [`OVERHANG_ALPHA`] each).
    fn searcher_for(&mut self, revcomp: bool, partial: bool) -> &mut Searcher<Iupac> {
        let slot = match (revcomp, partial) {
            (false, false) => &mut self.fwd,
            (true, false) => &mut self.rc,
            (false, true) => &mut self.fwd_overhang,
            (true, true) => &mut self.rc_overhang,
        };
        slot.get_or_insert_with(|| match (revcomp, partial) {
            (false, false) => Searcher::<Iupac>::new_fwd(),
            (true, false) => Searcher::<Iupac>::new_rc(),
            (false, true) => Searcher::<Iupac>::new_fwd_with_overhang(OVERHANG_ALPHA),
            (true, true) => Searcher::<Iupac>::new_rc_with_overhang(OVERHANG_ALPHA),
        })
    }
}

/// Match a compiled group against a record's per-file base segments, using its
/// own `loc` window.
pub fn match_group(
    group: &CompiledGroup,
    scratch: &mut Scratch,
    segments: &[&[u8]],
) -> GroupOutcome {
    match_group_at(group, scratch, segments, group.loc.as_ref())
}

/// The verdict of matching a group in the window `loc` (the whole record when
/// `None`), cached for a windowed or `match=` group. A whole-record search is
/// not cached, as its window is unique per record. The verdict is a pure
/// function of the group and the keyed bytes, which cover every base the
/// kernel reads, so the cache is exact.
pub fn match_verdict(
    group: &CompiledGroup,
    scratch: &mut Scratch,
    segments: &[&[u8]],
    loc: Option<&Location>,
) -> GroupVerdict {
    let windows = resolve_windows(loc, segments);
    let whole_record = loc.is_none() && group.match_streams.is_none();
    let &[(file, start, end)] = windows.as_slice() else {
        return GroupVerdict::of(group, &match_group_at(group, scratch, segments, loc));
    };
    if whole_record {
        return GroupVerdict::of(group, &match_group_at(group, scratch, segments, loc));
    }
    let cached = scratch.cache.get(group.name.as_str());
    let reach = cached.map_or_else(|| key_reach(group), |cached| cached.reach);
    scratch.key_buf.clear();
    write_window_key(&mut scratch.key_buf, segments[file], start, end, reach);
    if let Some(hit) = cached.and_then(|cached| cached.verdicts.get(&scratch.key_buf)) {
        return hit.placed(file, start);
    }
    let verdict = GroupVerdict::of(group, &match_group_at(group, scratch, segments, loc));
    let Some(relative) = CachedVerdict::relative(&verdict, start) else {
        return verdict;
    };
    let key = scratch.key_buf.clone();
    let cached = scratch
        .cache
        .entry(group.name.clone())
        .or_insert_with(|| GroupCache {
            reach,
            verdicts: AHashMap::new(),
        });
    if cached.verdicts.len() < MATCH_CACHE_CAP {
        cached.verdicts.insert(key, relative);
    }
    verdict
}

/// The verdict of a relative `next=GROUP:lo-hi` target: `group` searched
/// `lo..hi` bases past `upstream` in the orientation of that match. Past a
/// reverse-strand match the window lies before it in the read, and `group` is
/// searched on the reverse complement there, so its strands are taken relative
/// to the upstream tag and a read and its reverse complement give the same
/// verdict past the same upstream match span.
pub fn match_verdict_after(
    group: &CompiledGroup,
    scratch: &mut Scratch,
    segments: &[&[u8]],
    upstream: &MatchSpan,
    reverse: bool,
    window: &RelativeWindow,
) -> GroupVerdict {
    if !reverse {
        let loc = window_after(upstream, window.lo, window.hi);
        return match_verdict(group, scratch, segments, Some(&loc));
    }
    let before = &segments[upstream.file][..upstream.start];
    // The search reads no further than `hi + key_reach` bases past the match.
    let needed = window.hi.saturating_add(key_reach(group)).min(before.len());
    let mut flipped = before[before.len() - needed..].to_vec();
    crate::iupac::reverse_complement(&mut flipped);
    let mut oriented: SmallVec<[&[u8]; 4]> = segments.iter().copied().collect();
    oriented[upstream.file] = &flipped;
    let loc = Location {
        file: upstream.file,
        start: Endpoint::FromStart(window.lo),
        end: Some(Endpoint::FromStart(window.hi)),
    };
    let verdict = match_verdict(group, scratch, &oriented, Some(&loc));
    GroupVerdict {
        best: verdict.best.map(|best| best.mirrored(upstream.start)),
        failure: verdict.failure,
    }
}

/// The bases from a window start that the kernel reads whatever the window end.
fn key_reach(group: &CompiledGroup) -> usize {
    if group.anchor != Some(Anchor::FivePrime) {
        return 0;
    }
    let longest = group.tags.iter().map(Vec::len).max().unwrap_or(0);
    longest.saturating_add(group.dist.indel)
}

/// Build the cache key for the window `[start, end)` of `segment`: the window
/// length, then the bytes the kernel reads, which are the window bytes and any
/// of the `reach` bases from the window start that lie past its end (see
/// [`key_reach`]). Two windows hit only when the same bytes are read, wherever
/// they start.
fn write_window_key(key: &mut WindowKey, segment: &[u8], start: usize, end: usize, reach: usize) {
    let key_end = end.max(start.saturating_add(reach).min(segment.len()));
    key.extend_from_slice(&((end - start) as u32).to_le_bytes());
    key.extend_from_slice(&segment[start..key_end]);
}

/// Match a compiled group with an explicit search window `loc`, overriding the
/// group's own, as for the window a relative `next=GROUP:lo-hi` link computes
/// from the upstream match; all other attributes
/// (`dist`/`bothStrands`/`mode`/`delta`/find counts) come from `group`. The
/// reusable [`Scratch`] supplies the thread-local searcher.
pub fn match_group_at(
    group: &CompiledGroup,
    scratch: &mut Scratch,
    segments: &[&[u8]],
    loc: Option<&Location>,
) -> GroupOutcome {
    // `anchor=5p`/`3p` bypass `sassy` entirely: each tag is matched over a span
    // pinned to the window's anchored edge rather than searched for anywhere in
    // it. This is the bounded case `sassy`'s find-anywhere search is unreliable
    // on (short windows), so the direct kernel both fixes that and is cheaper.
    // The grammar validates anchored groups to be forward-strand and free of
    // partial overhang. The substitution-only cases (`anchor=3p`, or
    // `anchor=5p` with `indel=0`) are a direct IUPAC Hamming compare over the
    // tag's own length; `anchor=5p` with an indel budget uses an anchored
    // banded edit-distance DP (see `match_anchored`), so it counts indels too.
    if group.anchor.is_some() {
        let accepted = match_anchored(group, segments, loc);
        let best = select_best(&accepted, group);
        return GroupOutcome { accepted, best };
    }
    let partial = group.partial5.is_some() || group.partial3.is_some();
    let searcher = scratch.searcher_for(group.revcomp, partial);
    // Tags handled by a precomputed batched bucket; the rest fall back to
    // per-tag search.
    let batched = &group.batched_tags;

    let mut accepted = Vec::new();
    for (file, win_start, win_end) in resolve_windows(loc, segments) {
        let text = &segments[file][win_start..win_end];
        if text.is_empty() {
            continue;
        }
        // Search a masked copy of the same length; extraction and output keep the read's bases.
        let masked: Vec<u8> = text.iter().map(|&base| mask_no_call(base)).collect();
        let text = masked.as_slice();
        let window_first = accepted.len();
        let mut dropped: SmallVec<[(usize, MatchStrand); 32]> = SmallVec::new();

        // Batched path: equal-length tag batches searched together across SIMD
        // lanes. Only populated for forward-strand non-overhang groups (see
        // `encode`), so every match is a full forward match and the per-type
        // `dist` budget gates directly, matching `accept_match`'s non-overhang
        // branch that the per-tag fallback applies.
        for batch in &group.encoded {
            for m in searcher.search_encoded_patterns(&batch.encoded, text, group.dist.total) {
                let (subs, indels) = count_edits(m);
                let tag_idx = batch.tag_indices[m.pattern_idx];
                if subs > group.dist.mismatch || indels > group.dist.indel {
                    dropped.push((tag_idx, strand_of(m.strand)));
                    continue;
                }
                accepted.push(TagMatch {
                    tag_idx,
                    file,
                    start: win_start + m.text_start,
                    end: win_start + m.text_end,
                    cost: m.cost as usize,
                    subs,
                    indels,
                    strand: strand_of(m.strand),
                });
            }
        }

        // Per-tag path for tags not in a batch (small or over-long buckets, and
        // partial groups).
        for (tag_idx, pattern) in group.patterns.iter().enumerate() {
            if pattern.is_empty() || batched.contains(&tag_idx) {
                continue;
            }
            for m in searcher.search(pattern, text, search_budget(group, pattern.len())) {
                let (subs, indels) = count_edits(&m);
                if !accept_match(group, pattern.len(), &m, subs, indels) {
                    dropped.push((tag_idx, strand_of(m.strand)));
                    continue;
                }
                accepted.push(TagMatch {
                    tag_idx,
                    file,
                    start: win_start + m.text_start,
                    end: win_start + m.text_end,
                    cost: m.cost as usize,
                    subs,
                    indels,
                    strand: strand_of(m.strand),
                });
            }
        }

        if !dropped.is_empty() {
            accept_hidden_matches(
                group,
                (file, win_start),
                text,
                &mut dropped,
                window_first,
                &mut accepted,
            );
        }
    }

    let best = select_best(&accepted, group);
    GroupOutcome { accepted, best }
}

/// Accept the cheapest alignment within the caps in each edit-cost valley that has no match.
fn accept_hidden_matches(
    group: &CompiledGroup,
    origin: (usize, usize),
    text: &[u8],
    dropped: &mut SmallVec<[(usize, MatchStrand); 32]>,
    window_first: usize,
    accepted: &mut Vec<TagMatch>,
) {
    let (file, win_start) = origin;
    dropped.sort_unstable_by_key(|&(tag_idx, strand)| (strand_rank(strand), tag_idx));
    dropped.dedup();
    for of_strand in dropped.chunk_by(|a, b| a.1 == b.1) {
        let strand = of_strand[0].1;
        let mut oriented = Cow::Borrowed(text);
        if strand == MatchStrand::ReverseComplement {
            crate::iupac::reverse_complement(oriented.to_mut());
        }
        let flip = |start: usize, end: usize| match strand {
            MatchStrand::Forward => (start, end),
            MatchStrand::ReverseComplement => (text.len() - end, text.len() - start),
        };
        let sites = base_sites(&oriented);
        let mut capped = EndSearch::new(&sites, text.len(), capped_states(&group.dist));
        for &(tag_idx, _) in of_strand {
            let tag = &group.tags[tag_idx];
            let mut hidden = capped.end_costs(tag);
            if hidden.is_empty() {
                continue;
            }
            let mut found: SmallVec<[(usize, usize); 4]> = accepted[window_first..]
                .iter()
                .filter(|m| m.tag_idx == tag_idx && m.strand == strand)
                .map(|m| (flip(m.start - win_start, m.end - win_start).1, m.cost))
                .collect();
            hidden.retain(|known| !found.contains(known));
            hidden.sort_unstable_by_key(|&(end, cost)| (cost, Reverse(end)));
            let edit_costs = LazyCell::new(|| {
                let states = edit_states(search_budget(group, tag.len()));
                EndSearch::new(&sites, text.len(), states).end_costs(tag)
            });
            for (end, cost) in hidden {
                let mut cheaper = found.iter().filter(|other| other.1 <= cost);
                if cheaper.any(|other| share_valley(&edit_costs, other.0, end)) {
                    continue;
                }
                let Some(hit) = align_ending_at(tag, &oriented, end, &group.dist) else {
                    continue;
                };
                found.push((end, cost));
                let (start, end) = flip(end - hit.end, end);
                accepted.push(TagMatch {
                    tag_idx,
                    file,
                    start: win_start + start,
                    end: win_start + end,
                    cost,
                    subs: hit.subs,
                    indels: hit.indels,
                    strand,
                });
            }
        }
    }
}

/// The best alignment of the tag that ends at `end`; the hit's `end` counts the bases it spans.
fn align_ending_at(tag: &[u8], text: &[u8], end: usize, dist: &Dist) -> Option<AnchoredHit> {
    if dist.indel == 0 {
        return anchored_hamming(tag, text, end.checked_sub(tag.len())?, dist.mismatch);
    }
    align_from_edge(tag.len(), end, dist, |j, i| {
        iupac_overlap(text[end - 1 - i], tag[tag.len() - 1 - j])
    })
}

/// For each 4-bit IUPAC mask, the bit set of text positions that match it.
fn base_sites(text: &[u8]) -> Vec<u64> {
    let words = text.len() / 64 + 1;
    let mut sites = vec![0u64; 16 * words];
    for (w, chunk) in text.chunks(64).enumerate() {
        for (bit, &base) in chunk.iter().enumerate() {
            let mask = crate::iupac::mask(base) as usize;
            sites[mask * words + w] |= 1 << bit;
        }
        sites[w] = 0;
        for mask in 1..16 {
            let rest = mask & (mask - 1);
            sites[mask * words + w] = sites[rest * words + w] | sites[(mask ^ rest) * words + w];
        }
    }
    sites
}

/// One [`EndSearch`] state: its cost, and the states it extends by one more edit.
struct EditState {
    cost: usize,
    substitute: Option<usize>,
    indel: Option<usize>,
}

/// The states for every count of substitutions and indels within the caps.
fn capped_states(dist: &Dist) -> Vec<EditState> {
    let mut states = Vec::new();
    let mut fewer_subs = 0;
    for s in 0..=dist.mismatch.min(dist.total) {
        let first = states.len();
        for d in 0..=dist.indel.min(dist.total - s) {
            states.push(EditState {
                cost: s + d,
                substitute: (s > 0).then(|| fewer_subs + d),
                indel: (d > 0).then(|| first + d - 1),
            });
        }
        fewer_subs = first;
    }
    states
}

/// The states that allow up to `0..=total` edits of any type.
fn edit_states(total: usize) -> Vec<EditState> {
    let state = |cost: usize| EditState {
        cost,
        substitute: cost.checked_sub(1),
        indel: cost.checked_sub(1),
    };
    (0..=total).map(state).collect()
}

/// An exact search for where a tag's alignments end, bit-parallel over the text offsets.
struct EndSearch<'a> {
    sites: &'a [u64],
    text_len: usize,
    states: Vec<EditState>,
    rows: Vec<u64>,
}

impl<'a> EndSearch<'a> {
    fn new(sites: &'a [u64], text_len: usize, states: Vec<EditState>) -> Self {
        let rows = vec![0; (2 * states.len() + 1) * (sites.len() / 16)];
        Self {
            sites,
            text_len,
            states,
            rows,
        }
    }

    /// The `(end, cost)` of the cheapest alignment at each end that has one, in end order.
    fn end_costs(&mut self, tag: &[u8]) -> Vec<(usize, usize)> {
        let words = self.sites.len() / 16;
        let (none, rows) = self.rows.split_at_mut(words);
        let (mut row, mut next) = rows.split_at_mut(self.states.len() * words);
        row.fill(!0);
        for &base in tag {
            let mask = crate::iupac::mask(base) as usize;
            let same = &self.sites[mask * words..][..words];
            for (i, state) in self.states.iter().enumerate() {
                let (done, out) = next.split_at_mut(i * words);
                let out = &mut out[..words];
                let matched = &row[i * words..][..words];
                let substituted = match state.substitute {
                    Some(from) => &row[from * words..][..words],
                    None => &none[..words],
                };
                let (deleted, inserted) = match state.indel {
                    Some(from) => (
                        &row[from * words..][..words],
                        &done[from * words..][..words],
                    ),
                    None => (&none[..words], &none[..words]),
                };
                let mut carry = 0;
                for w in 0..words {
                    let step = (matched[w] & same[w]) | substituted[w] | inserted[w];
                    out[w] = step << 1 | carry | deleted[w];
                    carry = step >> 63;
                }
            }
            std::mem::swap(&mut row, &mut next);
        }

        let mut costs = Vec::new();
        for w in 0..words {
            let states = || self.states.iter().zip(row[w..].iter().step_by(words));
            let in_text = !0 >> (64 - (self.text_len + 1 - w * 64).min(64));
            let mut ends = states().fold(0, |ends, (_, row)| ends | row) & in_text;
            while ends != 0 {
                let bit = ends.trailing_zeros();
                ends &= ends - 1;
                let within = states().filter(|(_, row)| *row >> bit & 1 == 1);
                let cost = within.map(|(state, _)| state.cost).min();
                costs.extend(cost.map(|cost| (w * 64 + bit as usize, cost)));
            }
        }
        costs
    }
}

/// Whether the edit cost, given in end order, has no peak between the ends `a` and `b`.
fn share_valley(costs: &[(usize, usize)], a: usize, b: usize) -> bool {
    let index = |end| costs.binary_search_by_key(&end, |&(e, _)| e);
    let (Ok(i), Ok(j)) = (index(a.min(b)), index(a.max(b))) else {
        return false;
    };
    let mut steps = costs[i..=j]
        .windows(2)
        .map(|pair| pair[0].1.cmp(&pair[1].1));
    let mut rising = false;
    j - i == a.abs_diff(b)
        && steps.all(|step| {
            rising |= step == Ordering::Less;
            !(rising && step == Ordering::Greater)
        })
}

/// The best anchored alignment of one tag: its substitution and indel counts
/// and the window offset (bases consumed) where the tag's 3' end lands.
struct AnchoredHit {
    subs: usize,
    indels: usize,
    end: usize,
}

/// Whether an observed record base and a tag base share ACGT membership (a
/// match). A record no-call/`N` or record-side degenerate is masked to
/// never-match (charged as a substitution, like the `sassy` path); a tag-side
/// IUPAC code stays a free wildcard.
fn iupac_overlap(observed: u8, tag_base: u8) -> bool {
    crate::iupac::mask(mask_no_call(observed)) & crate::iupac::mask(tag_base) != 0
}

/// Match an anchored group (`anchor=5p`/`3p`): each tag's anchored edge is
/// fixed at the window boundary and matched over its own length, with **no
/// `sassy` search** (so the short-window result-selection limitation cannot
/// apply; deterministic). The `sassy` kernel is semi-global with a free
/// text-start, which would let the tag float within the window (the very
/// sliding that anchoring prevents), so this is a direct kernel: a Hamming
/// compare over the anchored span for the substitution-only case (`indel=0`),
/// or, for `anchor=5p`, an anchored banded edit-distance DP (5' pinned, 3'
/// floats by the indel budget) when indels are allowed. A tag that cannot fit
/// anchored does not match.
fn match_anchored(
    group: &CompiledGroup,
    segments: &[&[u8]],
    loc: Option<&Location>,
) -> Vec<TagMatch> {
    let mut accepted = Vec::new();
    for (file, win_start, win_end) in resolve_windows(loc, segments) {
        let segment = segments[file];
        for (tag_idx, tag) in group.tags.iter().enumerate() {
            if tag.is_empty() {
                continue;
            }
            // The anchored start: the window start for anchor=5p, or `win_end -
            // len(tag)` for anchor=3p (the tag's 3' base pinned at the window
            // end). A tag too long to fit anchored is skipped. anchor=3p is
            // substitution-only (the grammar rejects an indel budget for it).
            let (start, hit) = if group.anchor == Some(Anchor::ThreePrime) {
                let Some(start) = win_end.checked_sub(tag.len()).filter(|&s| s >= win_start) else {
                    continue;
                };
                (
                    start,
                    anchored_hamming(tag, segment, start, group.dist.mismatch),
                )
            } else if group.dist.indel == 0 {
                (
                    win_start,
                    anchored_hamming(tag, segment, win_start, group.dist.mismatch),
                )
            } else {
                let avail = segment.len().saturating_sub(win_start);
                let n = (tag.len() + group.dist.indel).min(avail);
                (
                    win_start,
                    anchored_align(tag, &segment[win_start..win_start + n], &group.dist),
                )
            };
            if let Some(hit) = hit {
                accepted.push(TagMatch {
                    tag_idx,
                    file,
                    start,
                    end: start + hit.end,
                    cost: hit.subs + hit.indels,
                    subs: hit.subs,
                    indels: hit.indels,
                    strand: MatchStrand::Forward,
                });
            }
        }
    }
    accepted
}

/// The substitution-only anchored match (`indel=0`): a direct Hamming compare
/// of `tag` against `segment[win_start .. win_start+len(tag)]`, or `None` if
/// the tag runs off the record end or exceeds the mismatch budget.
/// Allocation-free; this is the common, hot anchored path.
fn anchored_hamming(
    tag: &[u8],
    segment: &[u8],
    win_start: usize,
    max_subs: usize,
) -> Option<AnchoredHit> {
    if win_start + tag.len() > segment.len() {
        return None;
    }
    let mut subs = 0;
    for (offset, &tag_base) in tag.iter().enumerate() {
        if !iupac_overlap(segment[win_start + offset], tag_base) {
            subs += 1;
            if subs > max_subs {
                return None;
            }
        }
    }
    Some(AnchoredHit {
        subs,
        indels: 0,
        end: tag.len(),
    })
}

/// The best anchored edit-distance alignment of `tag` against `window` (the
/// record bases from the anchor onward), with the tag's 5' base pinned at
/// `window[0]` (no free text-start) and the 3' end floating by the indel
/// budget. Returns the lowest-cost alignment within the
/// `dist=mismatch:indel:total` caps (then fewest indels, then the 3' end
/// closest to the nominal tag length, then the shortest span), or `None` if
/// none fits the caps.
fn anchored_align(tag: &[u8], window: &[u8], dist: &Dist) -> Option<AnchoredHit> {
    align_from_edge(tag.len(), window.len(), dist, |j, i| {
        iupac_overlap(window[i], tag[j])
    })
}

/// The best alignment from a pinned edge; `overlap(j, i)` compares tag base `j` and text base `i`.
fn align_from_edge(
    tag_len: usize,
    text_len: usize,
    dist: &Dist,
    overlap: impl Fn(usize, usize) -> bool,
) -> Option<AnchoredHit> {
    const UNREACHED: usize = usize::MAX;
    let mismatch = dist.mismatch.min(tag_len);
    // Deleting the whole tag costs `tag_len` indels, so the best alignment never has more.
    let indel = dist.indel.min(tag_len);
    let within = |subs: usize, indels: usize| subs <= mismatch && subs + indels <= dist.total;

    // Fewest substitutions by count of deleted tag bases (row) and inserted text bases (column).
    let width = indel + 1;
    let mut cells = SmallVec::<[usize; 16]>::from_elem(UNREACHED, width * width);
    let fewest = cells.as_mut_slice();
    fewest[0] = 0;
    for j in 0..tag_len {
        let mut alive = false;
        for del in (0..=indel).rev() {
            for ins in 0..=indel - del {
                let at = del * width + ins;
                let subs = std::mem::replace(&mut fewest[at], UNREACHED);
                if !within(subs, del + ins) {
                    continue;
                }
                alive = true;
                let i = j + ins - del;
                if i < text_len {
                    fewest[at] = subs + usize::from(!overlap(j, i));
                }
                if del + ins < indel {
                    fewest[at + width] = fewest[at + width].min(subs);
                    if i < text_len {
                        fewest[at + 1] = fewest[at + 1].min(subs);
                    }
                }
            }
        }
        if !alive {
            return None;
        }
    }

    (0..=indel)
        .flat_map(|del| (0..=indel - del).map(move |ins| (del, ins)))
        .map(|(del, ins)| (fewest[del * width + ins], del, ins))
        .filter(|&(subs, del, ins)| within(subs, del + ins))
        .map(|(subs, del, ins)| AnchoredHit {
            subs,
            indels: del + ins,
            end: tag_len + ins - del,
        })
        .min_by_key(|hit| {
            let shift = hit.end.abs_diff(tag_len);
            (hit.subs + hit.indels, hit.indels, shift, hit.end)
        })
}

/// Map a `sassy` strand to our [`MatchStrand`].
fn strand_of(strand: Strand) -> MatchStrand {
    match strand {
        Strand::Fwd => MatchStrand::Forward,
        Strand::Rc => MatchStrand::ReverseComplement,
    }
}

/// Resolve the search windows as `(file, start, end)` triples. A `loc` pins one
/// file and a slice; no `loc` searches every file's whole record.
fn resolve_windows(
    loc: Option<&Location>,
    segments: &[&[u8]],
) -> SmallVec<[(usize, usize, usize); 4]> {
    let mut windows = SmallVec::new();
    match loc {
        Some(loc) => {
            // A `loc` pins exactly one file (when it is in range); no inner Vec
            // is needed.
            if loc.file < segments.len() {
                let len = segments[loc.file].len();
                let start = loc.start.resolve(len);
                let end = loc.end.map(|e| e.resolve(len)).unwrap_or(len);
                if start <= end {
                    windows.push((loc.file, start, end));
                }
            }
        }
        None => {
            for (file, segment) in segments.iter().enumerate() {
                windows.push((file, 0, segment.len()));
            }
        }
    }
    windows
}

/// Map an observed record base to a never-match sentinel (`X`, which the
/// `Iupac` profile encodes as matching nothing) unless it is a concrete
/// A/C/G/T, reading an RNA U as T. An `N` no-call or a record-side degenerate
/// code is thus a mismatch against a concrete tag base, matching splitcode (a
/// tag-side IUPAC code remains a free wildcard, since only the record window is
/// masked).
fn mask_no_call(base: u8) -> u8 {
    match base {
        b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' => base,
        b'U' => b'T',
        b'u' => b't',
        _ => b'X',
    }
}

/// The [`CompiledGroup::patterns`] for `tags`.
pub(crate) fn search_patterns(tags: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let read_u_as_t = |&base: &u8| match base {
        b'U' => b'T',
        b'u' => b't',
        base => base,
    };
    tags.iter()
        .map(|tag| tag.iter().map(read_u_as_t).collect())
        .collect()
}

/// The per-overhang-base cost for partial matching. At `1.0` a truncated tag
/// base costs the same as a mismatch, so a full match always outranks a
/// truncated one and the reported `cost` stays integral (the
/// `min_match`/`max_mismatch_freq` post-filter does the real gating).
const OVERHANG_ALPHA: f32 = 1.0;

/// The search budget (sassy `k`) for one tag: the `dist` total, plus the
/// overhang the group's partial5/3 policies permit (so an edge-truncated tag is
/// found before [`accept_match`] applies the min-match / frequency
/// post-filter). With no partial policy the budget is just `dist.total`.
fn search_budget(group: &CompiledGroup, tag_len: usize) -> usize {
    match [group.partial5, group.partial3]
        .into_iter()
        .flatten()
        .map(|partial| partial.min_match)
        .min()
    {
        Some(min_match) => group.dist.total + tag_len.saturating_sub(min_match),
        None => group.dist.total,
    }
}

/// Whether a sassy match passes the group's acceptance rules. A full
/// (non-overhang) match is gated by the per-type `dist` budget. A match
/// truncated at the record's 5'/3' end (sassy overhang) is allowed only when
/// each overhanging end has a `partial5`/`partial3` policy it satisfies: at
/// least `min_match` matched bases and (no indels) a substitution count within
/// the policy's `max_mismatch_freq` per matched base, or, when that frequency
/// is unset, within the group's `dist` substitution budget.
fn accept_match(
    group: &CompiledGroup,
    tag_len: usize,
    m: &SassyMatch,
    subs: usize,
    indels: usize,
) -> bool {
    let left_overhang = m.pattern_start;
    let right_overhang = tag_len.saturating_sub(m.pattern_end);
    if left_overhang == 0 && right_overhang == 0 {
        return subs <= group.dist.mismatch && indels <= group.dist.indel;
    }
    if (left_overhang > 0 && group.partial5.is_none())
        || (right_overhang > 0 && group.partial3.is_none())
    {
        return false;
    }
    if indels > 0 {
        return false;
    }
    let matched = m.text_end - m.text_start;
    for (overhang, partial) in [
        (left_overhang, group.partial5),
        (right_overhang, group.partial3),
    ] {
        if overhang > 0 {
            let partial = partial.expect("an overhanging end has a partial policy");
            let subs_ok = match partial.max_mismatch_freq {
                Some(freq) => subs as f64 <= freq * matched as f64,
                None => subs <= group.dist.mismatch,
            };
            if matched < partial.min_match || !subs_ok {
                return false;
            }
        }
    }
    matched > 0
}

/// Count substitutions and indels (insertions plus deletions) in a match's
/// CIGAR.
fn count_edits(m: &SassyMatch) -> (usize, usize) {
    let mut subs = 0;
    let mut indels = 0;
    for element in &m.cigar.ops {
        match element.op {
            CigarOp::Sub => subs += element.cnt as usize,
            CigarOp::Ins | CigarOp::Del => indels += element.cnt as usize,
            _ => {}
        }
    }
    (subs, indels)
}

/// Pick the single best match by the deterministic total order, honoring
/// `mode=nearest`/`delta`.
fn select_best(accepted: &[TagMatch], group: &CompiledGroup) -> Option<TagMatch> {
    let best_index =
        (0..accepted.len()).min_by(|&i, &j| total_order(&accepted[i], &accepted[j]))?;
    let best = &accepted[best_index];

    // In nearest mode the best must beat the next-lowest-cost match by `delta`,
    // else the read is ambiguous and left unmatched. `delta=0` imposes no gap.
    if group.mode == MatchMode::Nearest && group.delta > 0 {
        // The runner-up is the best match of a DIFFERENT tag: a second
        // occurrence of the chosen tag is the same barcode, not a competing
        // candidate, so it must not trigger the delta rejection.
        let best_tag = best.tag_idx;
        let runner_up = accepted
            .iter()
            .filter(|m| m.tag_idx != best_tag)
            .min_by(|a, b| total_order(a, b));
        if let Some(runner_up) = runner_up {
            if runner_up.cost < best.cost + group.delta {
                return None;
            }
        }
    }
    Some(*best)
}

/// The deterministic total order for picking the single best match: lowest
/// cost, then lowest file index, then leftmost start, then forward before
/// reverse complement, then tag definition order, then shortest span.
fn total_order(a: &TagMatch, b: &TagMatch) -> Ordering {
    a.cost
        .cmp(&b.cost)
        .then(a.file.cmp(&b.file))
        .then(a.start.cmp(&b.start))
        .then(strand_rank(a.strand).cmp(&strand_rank(b.strand)))
        .then(a.tag_idx.cmp(&b.tag_idx))
        .then((a.end - a.start).cmp(&(b.end - b.start)))
}

/// Forward sorts before reverse complement in the total order.
fn strand_rank(strand: MatchStrand) -> u8 {
    match strand {
        MatchStrand::Forward => 0,
        MatchStrand::ReverseComplement => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A group with a single tag and the given attributes left at defaults
    /// unless overridden.
    fn group(tags: &[&str]) -> CompiledGroup {
        CompiledGroup::new("g", tags.iter().map(|t| t.as_bytes().to_vec()).collect())
    }

    /// Match a single-segment record.
    fn outcome(group: &CompiledGroup, read: &str) -> GroupOutcome {
        match_group(group, &mut Scratch::new(), &[read.as_bytes()])
    }

    /// A budget of `n` substitutions and no indels.
    fn substitutions(n: usize) -> Dist {
        Dist {
            mismatch: n,
            indel: 0,
            total: n,
        }
    }

    /// The `[start, end)` window of the first file.
    fn window(start: usize, end: usize) -> Option<Location> {
        Some(Location {
            file: 0,
            start: Endpoint::FromStart(start),
            end: Some(Endpoint::FromStart(end)),
        })
    }

    /// A deterministic SplitMix64 generator.
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

        fn bases(&mut self, len: usize) -> Vec<u8> {
            self.sample(len, b"ACGT")
        }
    }

    #[test]
    fn test_match_group_exact_match_coordinates() {
        let g = group(&["ACGTACGT"]);
        let result = outcome(&g, "TTTACGTACGTCC");
        let best = result.best.expect("a match");
        assert_eq!(best.start, 3);
        assert_eq!(best.end, 11);
        assert_eq!(best.cost, 0);
        assert_eq!(best.strand, MatchStrand::Forward);
        assert_eq!(best.tag_idx, 0);
    }

    #[test]
    fn test_anchored_length_strict_no_interior_match() {
        // The split-seq-like grp_cbt failure: variable-length 5'-anchored tags.
        // Read AGTACTCT, true tag AGTACTC (7nt). Under an unanchored window a
        // shorter tag TACTCA (6nt) slides to offset 2 (TACTCT vs TACTCA, 1
        // mismatch) -> a spurious 2nd hit. anchor=5p pins each tag at the start
        // over its own length, so TACTCA is only tested at [0,6)=AGTACT (5
        // mismatches) -> no spurious match.
        let mut g = group(&["AGTACTC", "TACTCA"]);
        g.dist = Dist {
            mismatch: 1,
            indel: 0,
            total: 1,
        };
        g.loc = Some(Location {
            file: 0,
            start: Endpoint::FromStart(0),
            end: Some(Endpoint::FromStart(8)),
        });

        // an unanchored window (today's default): both the true tag and the
        // spurious offset-2 tag match.
        let window = outcome(&g, "AGTACTCT");
        assert_eq!(
            window.accepted.len(),
            2,
            "window slides -> spurious 2nd hit"
        );

        // anchor=5p: only the true tag matches at the anchored start.
        g.anchor = Some(Anchor::FivePrime);
        let anchored = outcome(&g, "AGTACTCT");
        assert_eq!(anchored.accepted.len(), 1, "anchored -> single hit");
        let best = anchored.best.expect("a match");
        assert_eq!(best.tag_idx, 0, "the true tag AGTACTC");
        assert_eq!((best.start, best.end), (0, 7));
        assert_eq!(best.cost, 0);
        assert_eq!(best.indels, 0);
    }

    #[test]
    fn test_anchored_uniform_length_is_window_equivalent() {
        // For equal-length tags fully inside the window, anchored and window
        // agree.
        let make = |place| {
            let mut g = group(&["ACGT", "TTTT"]);
            g.dist = Dist {
                mismatch: 1,
                indel: 0,
                total: 1,
            };
            g.loc = Some(Location {
                file: 0,
                start: Endpoint::FromStart(0),
                end: Some(Endpoint::FromStart(4)),
            });
            g.anchor = place;
            outcome(&g, "ACGA").best.map(|m| (m.tag_idx, m.cost))
        };
        assert_eq!(make(None), Some((0, 1)));
        assert_eq!(make(Some(Anchor::FivePrime)), Some((0, 1)));
    }

    #[test]
    fn test_anchored_read_n_is_a_mismatch() {
        // A read no-call (N) is charged as a substitution against a concrete
        // tag base, like the sassy path; an anchored tag that does not fit the
        // read does not match.
        let mut g = group(&["ACGT"]);
        g.anchor = Some(Anchor::FivePrime);
        g.loc = Some(Location {
            file: 0,
            start: Endpoint::FromStart(0),
            end: Some(Endpoint::FromStart(4)),
        });
        // dist exact (0): the N at position 1 is a mismatch, so no match.
        assert!(outcome(&g, "ANGT").best.is_none());
        // The tag runs past the 3-base read: no match.
        assert!(outcome(&g, "ACG").best.is_none());
    }

    #[test]
    fn test_anchor5p_anchors_at_given_window_start() {
        // anchor5p pins the tag at whatever window start match_group_at is
        // handed (e.g. a relative next=GROUP:lo-hi window), not at read
        // position 0. ACGT sits at [3,7); anchored at 3 it matches, anchored at
        // 0 (where the read is TTTA) it does not.
        let mut g = group(&["ACGT"]);
        g.anchor = Some(Anchor::FivePrime);
        let read: &[u8] = b"TTTACGTGG";
        let at_three = Location {
            file: 0,
            start: Endpoint::FromStart(3),
            end: None,
        };
        let best = match_group_at(&g, &mut Scratch::new(), &[read], Some(&at_three))
            .best
            .expect("anchored at the given window start");
        assert_eq!((best.start, best.end), (3, 7));
        let at_zero = Location {
            file: 0,
            start: Endpoint::FromStart(0),
            end: None,
        };
        assert!(
            match_group_at(&g, &mut Scratch::new(), &[read], Some(&at_zero))
                .best
                .is_none(),
            "anchored at 0 does not match (read[0..4] = TTTA)"
        );
    }

    #[test]
    fn test_anchor5p_indels() {
        // anchor=5p with an indel budget: the 5' base stays pinned at the
        // window start and the 3' end floats by the indel count, via the
        // anchored edit-distance kernel.
        let mk = |tags: &[&str], mm, ind, tot| {
            let mut g = group(tags);
            g.anchor = Some(Anchor::FivePrime);
            g.dist = Dist {
                mismatch: mm,
                indel: ind,
                total: tot,
            };
            g.loc = Some(Location {
                file: 0,
                start: Endpoint::FromStart(0),
                end: None,
            });
            g
        };
        // 1 deletion: read is ACGTACGT with the index-3 T missing (ACGACGT) ->
        // 1 indel, 3' end pulls in.
        let del = outcome(&mk(&["ACGTACGT"], 0, 1, 1), "ACGACGTGG")
            .best
            .expect("deletion match");
        assert_eq!((del.subs, del.indels), (0, 1));
        assert_eq!(
            (del.start, del.end),
            (0, 7),
            "deletion pulls the 3' end in by one"
        );
        // 1 insertion: an extra (no-call) base at index 4 -> 1 indel, 3' end
        // pushes out.
        let ins = outcome(&mk(&["ACGTACGT"], 0, 1, 1), "ACGTNACGT")
            .best
            .expect("insertion match");
        assert_eq!((ins.subs, ins.indels), (0, 1));
        assert_eq!(
            (ins.start, ins.end),
            (0, 9),
            "insertion pushes the 3' end out by one"
        );
        // A pure 1-substitution read under dist=0:1 does NOT match: the sub cap
        // is 0, and an indel detour would cost two indels (over budget), so
        // there is no valid alignment.
        assert!(
            outcome(&mk(&["ACGT"], 0, 1, 1), "AGGTAA").best.is_none(),
            "1 sub exceeds the sub cap with no cheaper indel alignment"
        );
        // sub + indel together within dist=1:1: one substitution and one
        // deletion.
        let both = outcome(&mk(&["ACGTACGT"], 1, 1, 2), "ACCACGTGG")
            .best
            .expect("sub+del match");
        assert_eq!((both.subs, both.indels), (1, 1));
    }

    /// An `anchor=5p` group over `tags` with a `loc` of `[0, end)` in file 0.
    fn anchored5p(tags: &[&str], end: usize, dist: Dist) -> CompiledGroup {
        let mut g = group(tags);
        g.anchor = Some(Anchor::FivePrime);
        g.dist = dist;
        g.loc = Some(Location {
            file: 0,
            start: Endpoint::FromStart(0),
            end: Some(Endpoint::FromStart(end)),
        });
        g
    }

    /// A best match as `(tag, start, end)`, or `None` for an unmatched read.
    type Best = Option<(usize, usize, usize)>;

    /// Assert each read's best match from one warm `Scratch`, in both orders.
    fn assert_best_in_both_orders(group: &CompiledGroup, expected: &[(&str, Best)]) {
        let mut order = expected.to_vec();
        for _ in 0..2 {
            let mut warm = Scratch::new();
            let found: Vec<(&str, Best)> = order
                .iter()
                .map(|&(read, _)| {
                    let segments = [read.as_bytes()];
                    let loc = group.loc.as_ref();
                    let best = match_verdict(group, &mut warm, &segments, loc).best;
                    (read, best.map(|m| (m.tag_idx, m.start, m.end)))
                })
                .collect();
            assert_eq!(found, order);
            order.reverse();
        }
    }

    #[test]
    fn test_anchor5p_indel_budget_outcome_ignores_earlier_reads() {
        let indel = Dist {
            mismatch: 0,
            indel: 1,
            total: 1,
        };
        // The reads differ at base 8, which the inserted G pulls into the tag.
        assert_best_in_both_orders(
            &anchored5p(&["ACGTACGT"], 8, indel),
            &[("ACGGTACGTAAAA", Some((0, 0, 9))), ("ACGGTACGAAAAA", None)],
        );
    }

    #[test]
    fn test_anchor5p_tag_past_window_end_outcome_ignores_earlier_reads() {
        assert_best_in_both_orders(
            &anchored5p(&["AAAAAACC", "AAAAAAGG"], 6, Dist::default()),
            &[
                ("AAAAAACCTTTT", Some((0, 0, 8))),
                ("AAAAAAGGTTTT", Some((1, 0, 8))),
            ],
        );
    }

    #[test]
    fn test_anchor5p_short_read_outcome_ignores_earlier_reads() {
        assert_best_in_both_orders(
            &anchored5p(&["AAAAAACC"], 6, Dist::default()),
            &[
                ("AAAAAACC", Some((0, 0, 8))),
                ("AAAAAAC", None),
                ("AAAAAA", None),
            ],
        );
        let indel = Dist {
            mismatch: 0,
            indel: 1,
            total: 1,
        };
        assert_best_in_both_orders(
            &anchored5p(&["ACGTACGT"], 8, indel),
            &[("ACGGTACGT", Some((0, 0, 9))), ("ACGGTACG", None)],
        );
    }

    #[test]
    fn test_key_reach_is_anchor5p_longest_tag_plus_indel_budget() {
        let mut g = group(&["ACGT", "ACGTAC"]);
        g.dist = Dist {
            mismatch: 1,
            indel: 2,
            total: 3,
        };
        assert_eq!(key_reach(&g), 0);
        g.anchor = Some(Anchor::ThreePrime);
        assert_eq!(key_reach(&g), 0);
        g.anchor = Some(Anchor::FivePrime);
        assert_eq!(key_reach(&g), 8);
    }

    #[test]
    fn test_write_window_key_covers_the_reach_within_the_read() {
        let key = |read: &str, start: usize, end: usize, reach: usize| {
            let mut key = WindowKey::new();
            let end = end.min(read.len());
            write_window_key(&mut key, read.as_bytes(), start.min(end), end, reach);
            key
        };
        let window = |read: &str, reach: usize| key(read, 2, 4, reach);
        assert_eq!(window("AAGTAAAA", 0), window("AAGTCCCC", 0));
        assert_eq!(window("AAGTAAAA", 2), window("AAGTCCCC", 2));
        assert_eq!(window("AAGTACGA", 5), window("AAGTACGC", 5));
        assert_ne!(window("AAGTACGA", 5), window("AAGTACTA", 5));
        assert_ne!(window("AAGTAC", 5), window("AAGTACG", 5));
        assert_ne!(window("A", 5), window("AAG", 5));
        assert_eq!(key("TTTTGTAC", 4, 6, 0), window("AAGTCCCC", 0));
        assert_ne!(key("AAGTCCCC", 2, 5, 0), window("AAGTCCCC", 0));
    }

    #[test]
    fn test_windows_at_other_offsets_share_a_cache_entry() {
        let mut g = group(&["ACGT", "TTGA"]);
        g.dist = substitutions(1);
        let mut warm = Scratch::new();
        let read: &[u8] = b"GGACGTCCTGACGTCC";
        let mut verdicts = Vec::new();
        for start in [1, 9, 1] {
            let loc = window(start, start + 6);
            verdicts.push(match_verdict(&g, &mut warm, &[read], loc.as_ref()));
            let fresh = match_group_at(&g, &mut Scratch::new(), &[read], loc.as_ref());
            assert_eq!(*verdicts.last().unwrap(), GroupVerdict::of(&g, &fresh));
        }
        assert_eq!(warm.cache["g"].verdicts.len(), 1);
        let starts: Vec<usize> = verdicts.iter().map(|v| v.best.unwrap().start).collect();
        assert_eq!(starts, [2, 10, 2]);
    }

    #[test]
    fn test_random_windows_warm_scratch_matches_fresh_scratch() {
        let mut rng = Rng(11);
        for anchor in [None, Some(Anchor::FivePrime), Some(Anchor::ThreePrime)] {
            for (mismatch, indel) in [(0, 0), (1, 0), (1, 1), (0, 2)] {
                if anchor == Some(Anchor::ThreePrime) && indel > 0 {
                    continue;
                }
                let mut g = group(&["ACGT", "AGGTA", "CCA"]);
                g.anchor = anchor;
                g.dist = Dist {
                    mismatch,
                    indel,
                    total: mismatch + indel,
                };
                g.mode = MatchMode::Nearest;
                g.max_finds_per_tag = Some(1);
                let mut warm = Scratch::new();
                for _ in 0..400 {
                    let len = rng.below(14);
                    let read = rng.sample(len, b"ACGTN");
                    let start = rng.below(len + 2);
                    let loc = window(start, start + rng.below(9));
                    let segments = [read.as_slice()];
                    let fresh = match_group_at(&g, &mut Scratch::new(), &segments, loc.as_ref());
                    assert_eq!(
                        match_verdict(&g, &mut warm, &segments, loc.as_ref()),
                        GroupVerdict::of(&g, &fresh),
                        "{read:?} in {loc:?} with {anchor:?} {:?}",
                        g.dist
                    );
                }
            }
        }
    }

    #[test]
    fn test_anchor5p_cache_splits_entries_within_the_reach_only() {
        let entries = |g: &CompiledGroup, reads: &[&str]| {
            let mut warm = Scratch::new();
            for read in reads {
                match_verdict(g, &mut warm, &[read.as_bytes()], g.loc.as_ref());
            }
            warm.cache["g"].verdicts.len()
        };
        let exact = anchored5p(&["ACGTACGT"], 8, Dist::default());
        assert_eq!(entries(&exact, &["ACGTACGTAAAA", "ACGTACGTCCCC"]), 1);
        let indel = Dist {
            mismatch: 0,
            indel: 1,
            total: 1,
        };
        let reads = ["ACGGTACGTAAAA", "ACGGTACGTCCCC", "ACGGTACGAAAAA"];
        assert_eq!(entries(&anchored5p(&["ACGTACGT"], 8, indel), &reads), 2);
    }

    /// Every read over `alphabet` of up to `max_len` bases, shortest first.
    fn reads_up_to(alphabet: &str, max_len: usize) -> Vec<String> {
        let mut reads = vec![String::new()];
        let mut longest = 0..1;
        for _ in 0..max_len {
            for index in longest.clone() {
                for base in alphabet.chars() {
                    reads.push(format!("{}{base}", reads[index]));
                }
            }
            longest = longest.end..reads.len();
        }
        reads
    }

    /// Assert a warm `Scratch` gives each read its fresh-`Scratch` verdict.
    fn assert_warm_scratch_matches_fresh(group: &CompiledGroup, reads: &[String]) {
        let fresh: Vec<GroupVerdict> = reads
            .iter()
            .map(|read| GroupVerdict::of(group, &outcome(group, read)))
            .collect();
        let mut order: Vec<usize> = (0..reads.len()).collect();
        for _ in 0..2 {
            let mut warm = Scratch::new();
            for &index in &order {
                let segments = [reads[index].as_bytes()];
                assert_eq!(
                    match_verdict(group, &mut warm, &segments, group.loc.as_ref()),
                    fresh[index],
                    "read {:?} with {:?}, {:?}",
                    reads[index],
                    group.loc,
                    group.dist
                );
            }
            order.reverse();
        }
    }

    /// Sweep `loc` windows, `(mismatch, indel)` budgets, and every short read.
    fn assert_sweep_matches_fresh(anchor: Option<Anchor>, budgets: &[(usize, usize)]) {
        let reads = reads_up_to("AC", 8);
        let windows = [
            (Endpoint::FromStart(0), Some(Endpoint::FromStart(2))),
            (Endpoint::FromStart(0), Some(Endpoint::FromStart(4))),
            (Endpoint::FromStart(1), Some(Endpoint::FromStart(3))),
            (Endpoint::FromStart(2), None),
            (Endpoint::FromEnd(4), Some(Endpoint::ReadEnd)),
            (Endpoint::FromStart(0), Some(Endpoint::FromEnd(2))),
        ];
        for (start, end) in windows {
            for &(mismatch, indel) in budgets {
                let mut g = group(&["ACA", "CCAC", "AACCA"]);
                g.anchor = anchor;
                g.dist = Dist {
                    mismatch,
                    indel,
                    total: mismatch + indel,
                };
                g.loc = Some(Location {
                    file: 0,
                    start,
                    end,
                });
                assert_warm_scratch_matches_fresh(&g, &reads);
            }
        }
    }

    #[test]
    fn test_anchor5p_warm_scratch_matches_fresh_scratch() {
        assert_sweep_matches_fresh(
            Some(Anchor::FivePrime),
            &[(0, 0), (1, 0), (0, 1), (1, 1), (0, 2)],
        );
    }

    #[test]
    fn test_anchor3p_warm_scratch_matches_fresh_scratch() {
        assert_sweep_matches_fresh(Some(Anchor::ThreePrime), &[(0, 0), (1, 0)]);
    }

    #[test]
    fn test_sliding_window_warm_scratch_matches_fresh_scratch() {
        assert_sweep_matches_fresh(None, &[(0, 0), (1, 1)]);
    }

    #[test]
    fn test_joined_stream_warm_scratch_matches_fresh_scratch() {
        let mut g = group(&["ACA", "CCAC", "AACCA"]);
        g.match_streams = Some(vec!["i7".to_string(), "i5".to_string()]);
        g.dist = Dist {
            mismatch: 1,
            indel: 1,
            total: 2,
        };
        assert_warm_scratch_matches_fresh(&g, &reads_up_to("AC", 8));
    }

    /// The `(subs, indels, end)` of the best alignment, by visiting every reachable state.
    fn exhaustive_align(tag: &[u8], window: &[u8], dist: &Dist) -> Option<(usize, usize, usize)> {
        let mut seen = AHashSet::new();
        let mut pending = vec![(0, 0, 0, 0)];
        while let Some(state) = pending.pop() {
            let (j, i, subs, indels) = state;
            let over_caps =
                subs > dist.mismatch || indels > dist.indel || subs + indels > dist.total;
            if over_caps || j > tag.len() || i > window.len() || !seen.insert(state) {
                continue;
            }
            pending.push((j + 1, i, subs, indels + 1));
            pending.push((j, i + 1, subs, indels + 1));
            if j < tag.len() && i < window.len() {
                let sub = usize::from(!iupac_overlap(window[i], tag[j]));
                pending.push((j + 1, i + 1, subs + sub, indels));
            }
        }
        seen.into_iter()
            .filter(|state| state.0 == tag.len())
            .map(|(_, end, subs, indels)| (subs, indels, end))
            .min_by_key(|&(subs, indels, end)| {
                (subs + indels, indels, end.abs_diff(tag.len()), end)
            })
    }

    #[test]
    fn test_anchored_alignments_equal_an_exhaustive_search() {
        let reversed = |seq: &[u8]| seq.iter().rev().copied().collect::<Vec<u8>>();
        let counts = |hit: AnchoredHit| (hit.subs, hit.indels, hit.end);
        let mut rng = Rng(0xa11e);
        for trial in 0..1200 {
            let letters: &[u8] = [&b"ACGT"[..], b"AC", b"ACGTNRY"][rng.below(3)];
            let (longest, most) = [(16, 1), (16, 2), (12, 4), (6, 20)][trial % 4];
            let len = rng.below(longest + 1);
            let tag = rng.sample(len, letters);
            let (mismatch, indel) = (rng.below(most + 1), rng.below(most + 1));
            let dist = Dist {
                mismatch,
                indel,
                total: mismatch + rng.below(indel + 2),
            };

            let (left, right) = (1 + rng.below(3), 3 + rng.below(4));
            let left = rng.sample(left, letters);
            let mut text = [left, tag.clone(), rng.sample(right, letters)].concat();
            for _ in 0..rng.below(4) {
                let at = rng.below(text.len());
                match rng.below(3) {
                    0 => text[at] = letters[rng.below(letters.len())],
                    1 => text.insert(at, letters[rng.below(letters.len())]),
                    _ => drop(text.remove(at)),
                }
            }

            for at in 0..=text.len() {
                let (before, after) = text.split_at(at);
                let context = format!(
                    "tag {} at {at} of {} {dist:?}",
                    String::from_utf8_lossy(&tag),
                    String::from_utf8_lossy(&text)
                );
                assert_eq!(
                    anchored_align(&tag, after, &dist).map(counts),
                    exhaustive_align(&tag, after, &dist),
                    "starting: {context}"
                );
                assert_eq!(
                    align_ending_at(&tag, &text, at, &dist).map(counts),
                    exhaustive_align(&reversed(&tag), &reversed(before), &dist),
                    "ending: {context}"
                );
            }
        }
    }

    #[test]
    fn test_match_group_exact_required_rejects_mismatch() {
        // dist defaults to exact, so a one-base substitution does not match.
        let g = group(&["ACGTACGT"]);
        let result = outcome(&g, "TTTACGTCCGTCC");
        assert!(result.best.is_none());
        assert!(result.accepted.is_empty());
    }

    #[test]
    fn test_match_group_one_substitution_within_dist() {
        let mut g = group(&["ACGTACGT"]);
        g.dist = Dist {
            mismatch: 1,
            indel: 0,
            total: 1,
        };
        // One substitution at the 5th base (A -> C).
        let result = outcome(&g, "TTTACGTCCGTCC");
        let best = result.best.expect("a match");
        assert_eq!(best.cost, 1);
        assert_eq!(best.subs, 1);
        assert_eq!(best.indels, 0);
    }

    #[test]
    fn test_substitution_is_matched_when_an_indel_alignment_hides_it() {
        // The last base is substituted; an equal-cost indel alignment ends one base later.
        let mut g = group(&["ACGTACGT"]);
        g.dist = substitutions(1);
        for loc in [None, window(0, 12), window(0, 8)] {
            g.loc = loc;
            let best = outcome(&g, "ACGTACGATGGCCAA").best.expect("a match");
            assert_eq!((best.start, best.end), (0, 8), "loc={:?}", g.loc);
            assert_eq!((best.cost, best.subs, best.indels), (1, 1, 0));
        }
    }

    #[test]
    fn test_substitutions_are_matched_when_a_cheaper_indel_alignment_hides_them() {
        // Hamming distance 2, where an alignment that deletes one tag base costs only 1.
        let mut g = group(&["AACCGGTT"]);
        g.dist = substitutions(2);
        g.loc = window(0, 8);
        let best = outcome(&g, "TAACGGTTCCCCAA").best.expect("a match");
        assert_eq!((best.start, best.end), (0, 8));
        assert_eq!((best.cost, best.subs, best.indels), (2, 2, 0));
    }

    #[test]
    fn test_planted_substitutions_within_budget_are_always_matched() {
        let mut rng = Rng(0x5eed);
        for trial in 0..600 {
            let len = 6 + rng.below(15);
            let mismatch = 1 + rng.below(2);
            let n_tags = [1, BATCHED_MIN_TAGS + 2][trial % 2];
            let mut g = CompiledGroup::new("g", (0..n_tags).map(|_| rng.bases(len)).collect());
            g.dist = substitutions(mismatch);
            g.revcomp = trial % 3 == 0;
            g.encode();

            let target = rng.below(n_tags);
            let mut planted = g.tags[target].clone();
            let mut planted_subs = 0;
            for _ in 0..rng.below(mismatch + 1) {
                let at = rng.below(len);
                planted[at] = b"ACGTN"[rng.below(5)];
            }
            for (a, b) in planted.iter().zip(&g.tags[target]) {
                planted_subs += usize::from(a != b);
            }
            if g.revcomp && rng.below(2) == 0 {
                crate::iupac::reverse_complement(&mut planted);
            }
            let left = rng.below(90);
            let right = rng.below(90);
            let read = [rng.bases(left), planted, rng.bases(right)].concat();
            let read = std::str::from_utf8(&read).unwrap();
            let end = left + len;

            for loc in [
                None,
                window(left, end),
                window(left.saturating_sub(4), end + 4),
                window(left, end + 4),
                window(left.saturating_sub(4), end),
                Some(Location {
                    file: 0,
                    start: Endpoint::FromEnd(right + len),
                    end: None,
                }),
            ] {
                g.loc = loc;
                let result = outcome(&g, read);
                assert!(
                    result.accepted.iter().any(|m| m.tag_idx == target
                        && m.cost <= planted_subs
                        && m.start < end
                        && left < m.end),
                    "tag {} with {planted_subs} substitutions at {left} in {read} \
                     dist={mismatch} loc={:?}: {result:?}",
                    String::from_utf8_lossy(&g.tags[target]),
                    g.loc
                );
            }
        }
    }

    #[test]
    fn test_each_tag_cost_equals_its_cheapest_alignment_within_caps() {
        let mut rng = Rng(0xcafe);
        for trial in 0..400 {
            let (mismatch, indel, total) = [
                (1, 0, 1),
                (2, 0, 2),
                (0, 1, 1),
                (1, 1, 1),
                (1, 1, 2),
                (2, 1, 2),
                (1, 2, 2),
                (2, 1, 3),
            ][trial % 8];
            let len = 6 + rng.below(8);
            // N, R, and Y are wildcards in a tag and mismatches in a read.
            let letters: &[u8] = [&b"ACGT"[..], b"AC", b"ACGTNRY"][rng.below(3)];
            let n_tags = [1, BATCHED_MIN_TAGS][rng.below(2)];
            let mut g =
                CompiledGroup::new("g", (0..n_tags).map(|_| rng.sample(len, letters)).collect());
            g.dist = Dist {
                mismatch,
                indel,
                total,
            };
            g.revcomp = trial % 3 == 0;
            g.encode();

            let mut planted = g.tags[0].clone();
            for _ in 0..rng.below(4) {
                let at = rng.below(planted.len());
                match rng.below(3) {
                    0 => planted[at] = letters[rng.below(letters.len())],
                    1 => planted.insert(at, letters[rng.below(letters.len())]),
                    _ => drop(planted.remove(at)),
                }
            }
            let (left, right) = (rng.below(20), rng.below(20));
            let left = rng.sample(left, letters);
            let mut read = [left, planted, rng.sample(right, letters)].concat();
            let mut other_strand = read.clone();
            crate::iupac::reverse_complement(&mut other_strand);
            if rng.below(2) == 0 {
                std::mem::swap(&mut read, &mut other_strand);
            }

            let result = outcome(&g, std::str::from_utf8(&read).unwrap());
            let strands = if g.revcomp {
                vec![&read, &other_strand]
            } else {
                vec![&read]
            };
            for (tag_idx, tag) in g.tags.iter().enumerate() {
                let cheapest = strands
                    .iter()
                    .flat_map(|text| {
                        (0..text.len()).filter_map(|start| {
                            let end = (start + len + indel).min(text.len());
                            anchored_align(tag, &text[start..end], &g.dist)
                        })
                    })
                    .map(|hit| hit.subs + hit.indels)
                    .min();
                let reported = result.accepted.iter().filter(|m| m.tag_idx == tag_idx);
                assert_eq!(
                    reported.map(|m| m.cost).min(),
                    cheapest,
                    "tag {} in {} dist={mismatch}:{indel}:{total}",
                    String::from_utf8_lossy(tag),
                    String::from_utf8_lossy(&read)
                );
            }
        }
    }

    #[test]
    fn test_hidden_matches_add_at_most_one_find_per_edit_cost_minimum() {
        let mut rng = Rng(0xf1d5);
        let mut searcher = Searcher::<Iupac>::new_fwd();
        for trial in 0..400 {
            let mut g = CompiledGroup::new("g", vec![rng.sample(8, b"AC")]);
            g.dist = [
                substitutions(1),
                substitutions(2),
                Dist {
                    mismatch: 1,
                    indel: 1,
                    total: 2,
                },
            ][trial % 3];
            let read = rng.sample(40, b"AC");
            let result = outcome(&g, std::str::from_utf8(&read).unwrap());
            let reported = searcher.search(&g.tags[0], &read, g.dist.total);
            assert!(result.accepted.len() <= reported.len(), "{g:?} {result:?}");
            for m in &reported {
                let (subs, indels) = count_edits(m);
                if subs <= g.dist.mismatch && indels <= g.dist.indel {
                    assert!(result
                        .accepted
                        .iter()
                        .any(|a| (a.start, a.end) == (m.text_start, m.text_end)));
                }
            }
        }
    }

    #[test]
    fn test_hidden_match_is_found_in_a_partial_group() {
        let mut g = group(&["ACGTACGT"]);
        g.dist = substitutions(1);
        g.partial5 = Some(Partial {
            min_match: 7,
            max_mismatch_freq: None,
        });
        let best = outcome(&g, "GGACGTACGATGGCCAA").best.expect("a match");
        assert_eq!((best.start, best.end, best.cost), (2, 10, 1));
    }

    #[test]
    fn test_match_group_indel_budget_rejects_substitution() {
        // dist=0:1:1 allows one indel but zero substitutions, so a sub-only
        // match is rejected.
        let mut g = group(&["ACGTACGT"]);
        g.dist = Dist {
            mismatch: 0,
            indel: 1,
            total: 1,
        };
        let result = outcome(&g, "TTTACGTCCGTCC");
        assert!(
            result.best.is_none(),
            "a substitution must not satisfy an indel-only budget: {result:?}"
        );
    }

    #[test]
    fn test_match_group_deletion_within_indel_budget() {
        // The read barcode region drops one base of the tag; allow one indel.
        let mut g = group(&["ACGTACGT"]);
        g.dist = Dist {
            mismatch: 0,
            indel: 1,
            total: 1,
        };
        let result = outcome(&g, "TTACGACGTTT");
        let best = result.best.expect("a deletion match");
        assert_eq!(best.cost, 1);
        assert_eq!(best.subs, 0);
        assert_eq!(best.indels, 1);
    }

    #[test]
    fn test_match_group_loc_window_restricts_search() {
        // The tag occurs only outside the window, so the windowed search finds
        // nothing.
        let mut g = group(&["ACGTACGT"]);
        g.loc = Some(Location {
            file: 0,
            start: Endpoint::FromStart(0),
            end: Some(Endpoint::FromStart(8)),
        });
        let result = outcome(&g, "GGGGGGGGGGACGTACGT");
        assert!(
            result.best.is_none(),
            "match lies past the window: {result:?}"
        );
    }

    #[test]
    fn test_match_group_loc_window_finds_in_window() {
        let mut g = group(&["ACGTACGT"]);
        g.loc = Some(Location {
            file: 0,
            start: Endpoint::FromStart(0),
            end: Some(Endpoint::FromStart(8)),
        });
        let result = outcome(&g, "ACGTACGTGGGGGG");
        let best = result.best.expect("a match in the window");
        assert_eq!(best.start, 0);
        assert_eq!(best.end, 8);
    }

    #[test]
    fn test_match_group_revcomp_strand() {
        // Tag AAAACCCC has reverse complement GGGGTTTT; place that in the read.
        let mut g = group(&["AAAACCCC"]);
        g.revcomp = true;
        let result = outcome(&g, "TTGGGGTTTTAA");
        let best = result.best.expect("a reverse-complement match");
        assert_eq!(best.strand, MatchStrand::ReverseComplement);
        assert_eq!(best.cost, 0);
    }

    #[test]
    fn test_match_group_forward_only_misses_revcomp() {
        let g = group(&["AAAACCCC"]);
        let result = outcome(&g, "TTGGGGTTTTAA");
        assert!(
            result.best.is_none(),
            "forward search must miss the rc occurrence"
        );
    }

    #[test]
    fn test_match_group_iupac_n_is_free() {
        // An expected N in the tag (pattern side) matches any base at zero cost
        // (barcode-N rule). The read here is all concrete; read-side N is
        // tested separately below.
        let g = group(&["ACGTNCGT"]);
        let result = outcome(&g, "TTACGTACGTCC");
        let best = result.best.expect("a match through the N");
        assert_eq!(best.cost, 0);
    }

    #[test]
    fn test_match_group_read_n_is_no_call() {
        // An observed N in the read is a no-call (a mismatch), NOT a free match
        // of a concrete tag.
        let g = group(&["ACGTACGT"]);
        let exact = outcome(&g, "TTACGTNCGTCC"); // N at the 5th barcode base
        assert!(
            exact.best.is_none(),
            "a read N must not match a concrete tag for free at dist=0: {exact:?}"
        );
        // Within a one-mismatch budget the no-call is charged as that mismatch.
        let mut tolerant = group(&["ACGTACGT"]);
        tolerant.dist = Dist {
            mismatch: 1,
            indel: 0,
            total: 1,
        };
        let result = outcome(&tolerant, "TTACGTNCGTCC");
        let best = result
            .best
            .expect("a read N matches within a one-mismatch budget");
        assert_eq!(best.cost, 1);
        assert_eq!(best.subs, 1);
    }

    #[test]
    fn test_match_group_read_iupac_is_no_call() {
        // A read-side degenerate code (R) is a no-call too, not a free
        // degenerate match.
        let g = group(&["ACGT"]);
        let result = outcome(&g, "ARGT"); // R at base 1
        assert!(
            result.best.is_none(),
            "a read-side R must be a no-call at dist=0: {result:?}"
        );
    }

    #[test]
    fn test_match_group_read_u_matches_a_tag_t() {
        let g = group(&["AATTC"]);
        for read in ["CCAAUUCCC", "ccaauuccc"] {
            let best = outcome(&g, read).best.expect("a read U matches a tag T");
            assert_eq!((best.start, best.end, best.cost), (2, 7, 0), "{read}");
        }
    }

    #[test]
    fn test_match_group_tag_u_matches_on_both_strands() {
        for tag in ["ACGTACGG", "ACGUACGG", "acguacgg"] {
            let mut g = group(&[tag]);
            g.revcomp = true;
            for (read, strand, span) in [
                ("ACGTACGGAAAA", MatchStrand::Forward, (0, 8)),
                ("TTTTCCGTACGT", MatchStrand::ReverseComplement, (4, 12)),
            ] {
                let best = outcome(&g, read)
                    .best
                    .unwrap_or_else(|| panic!("{tag} finds no match in {read}"));
                assert_eq!(
                    (best.strand, (best.start, best.end), best.cost),
                    (strand, span, 0),
                    "{tag} in {read}"
                );
            }
        }
    }

    /// The read with each T and t written as U and u.
    fn as_rna(read: &[u8]) -> Vec<u8> {
        read.iter()
            .map(|&base| match base {
                b'T' => b'U',
                b't' => b'u',
                base => base,
            })
            .collect()
    }

    #[test]
    fn test_u_bases_in_reads_and_tags_match_as_t_bases_do() {
        let mut rng = Rng(0x75);
        let codes = b"ACGTACGTACGTRYKMSWBDHVN";
        for trial in 0..2000 {
            let batched = trial % 4 == 0;
            let length = 3 + rng.below(5);
            let count = if batched {
                BATCHED_MIN_TAGS + rng.below(3)
            } else {
                1 + rng.below(4)
            };
            let tags = (0..count)
                .map(|_| {
                    let len = if batched { length } else { 3 + rng.below(5) };
                    let mut tag = rng.sample(len, codes);
                    if rng.below(4) == 0 {
                        tag.make_ascii_lowercase();
                    }
                    tag
                })
                .collect();
            let mut g = CompiledGroup::new("g", tags);
            let (mismatch, indel) = (rng.below(3), rng.below(3));
            g.dist = Dist {
                mismatch,
                indel,
                total: mismatch + indel,
            };
            g.mode = [MatchMode::All, MatchMode::Nearest][rng.below(2)];
            let partial = Partial {
                min_match: 2 + rng.below(3),
                max_mismatch_freq: None,
            };
            match rng.below(6) {
                _ if batched => {}
                0 => g.anchor = Some(Anchor::FivePrime),
                1 if indel == 0 => g.anchor = Some(Anchor::ThreePrime),
                2 => g.partial5 = Some(partial),
                3 => g.partial3 = Some(partial),
                _ => {}
            }
            g.revcomp = !batched && g.anchor.is_none() && rng.below(2) == 1;
            if rng.below(2) == 1 {
                let start = rng.below(4);
                g.loc = window(start, start + 4 + rng.below(16));
            } else if rng.below(2) == 1 {
                g.match_streams = Some(vec!["i7".to_string(), "i5".to_string()]);
            }
            g.encode();
            let mut rna_tags = g.clone();
            rna_tags.tags = g.tags.iter().map(|tag| as_rna(tag)).collect();
            rna_tags.encode();
            let len = rng.below(24);
            let read = rng.sample(len, b"ACGTacgtN");
            let rna = as_rna(&read);
            let dna = match_group(&g, &mut Scratch::new(), &[&read]);
            let tags: Vec<_> = rna_tags
                .tags
                .iter()
                .map(|t| String::from_utf8_lossy(t))
                .collect();
            let context = format!(
                "{} as {} for {tags:?} {:?} {:?} {:?} rc={} {:?} {:?}",
                String::from_utf8_lossy(&read),
                String::from_utf8_lossy(&rna),
                g.dist,
                g.loc,
                g.anchor,
                g.revcomp,
                g.partial5,
                g.partial3,
            );
            for (pair, group, read) in [
                ("T tags, U read", &g, &rna),
                ("U tags, T read", &rna_tags, &read),
                ("U tags, U read", &rna_tags, &rna),
            ] {
                assert_eq!(
                    match_group(group, &mut Scratch::new(), &[read]),
                    dna,
                    "{pair}: {context}"
                );
                assert_eq!(
                    match_verdict(group, &mut Scratch::new(), &[read], group.loc.as_ref()),
                    GroupVerdict::of(group, &dna),
                    "{pair}: {context}"
                );
            }
        }
    }

    #[test]
    fn test_match_group_best_is_lowest_cost() {
        // Two tags both occur; the exact one wins over the one-mismatch one.
        let mut g = group(&["ACGTACGT", "TTTTGGGG"]);
        g.dist = Dist {
            mismatch: 1,
            indel: 0,
            total: 1,
        };
        // ACGTACGT exact at 0..8; TTTTGGGG with one sub later.
        let result = outcome(&g, "ACGTACGTTTTTGGGC");
        let best = result.best.expect("a match");
        assert_eq!(best.tag_idx, 0);
        assert_eq!(best.cost, 0);
        assert!(result.accepted.len() >= 2, "both tags accepted: {result:?}");
    }

    #[test]
    fn test_match_group_variable_length_tags() {
        // A group may hold barcodes of differing lengths; each is searched
        // independently and the best across them is chosen, with its own span
        // length (splitcode supports this, so must we).
        let g = group(&["AAAAAA", "CCCCCCC", "GGGGGGGG"]); // 6, 7, 8 nt
        let seven = outcome(&g, "TTCCCCCCCTT");
        let best = seven.best.expect("the 7 nt tag matches");
        assert_eq!(best.tag_idx, 1);
        assert_eq!(best.cost, 0);
        assert_eq!(best.end - best.start, 7);
        let eight = outcome(&g, "TTGGGGGGGGTT");
        let best = eight.best.expect("the 8 nt tag matches");
        assert_eq!(best.tag_idx, 2);
        assert_eq!(best.end - best.start, 8);
    }

    #[test]
    fn test_match_group_nearest_rejects_ambiguous_ties() {
        // Two tags match the same region equally well; nearest mode with delta
        // cannot disambiguate.
        let mut g = group(&["ACGTACGT", "ACGTACGA"]);
        g.dist = Dist {
            mismatch: 1,
            indel: 0,
            total: 1,
        };
        g.mode = MatchMode::Nearest;
        g.delta = 1;
        // ACGTACGT exact (cost 0); ACGTACGA has one sub vs the same region
        // (cost 1). Gap is 1 -> ok.
        let clear = outcome(&g, "ACGTACGTCC");
        assert!(clear.best.is_some(), "a clear winner by delta 1: {clear:?}");

        // Now both at cost 0 against different equally-good regions ->
        // ambiguous.
        let g2 = {
            let mut g2 = group(&["AAAACCCC", "GGGGTTTT"]);
            g2.mode = MatchMode::Nearest;
            g2.delta = 1;
            g2
        };
        let tie = outcome(&g2, "AAAACCCCGGGGTTTT");
        assert!(
            tie.best.is_none(),
            "two cost-0 matches are ambiguous: {tie:?}"
        );
    }

    #[test]
    fn test_match_group_nearest_bare_group_rejects_ties() {
        let mut g = group(&["AAAACCCC", "GGGGTTTT"]);
        g.mode = MatchMode::Nearest;
        let tie = outcome(&g, "AAAACCCCGGGGTTTT");
        assert_eq!(tie.accepted.len(), 2, "both tags match: {tie:?}");
        assert!(
            tie.best.is_none(),
            "two cost-0 matches are ambiguous: {tie:?}"
        );
    }

    #[test]
    fn test_match_group_total_order_leftmost_breaks_cost_tie() {
        // The same tag occurs twice at cost 0; the leftmost wins.
        let g = group(&["ACGTACGT"]);
        let result = outcome(&g, "ACGTACGTACGTACGT");
        let best = result.best.expect("a match");
        assert_eq!(best.start, 0, "leftmost occurrence wins: {result:?}");
    }

    #[test]
    fn test_match_group_nearest_same_tag_twice_not_ambiguous() {
        // The same tag occurs twice at cost 0. In nearest mode the runner-up
        // considered for the delta gap must be the best of a DIFFERENT tag; a
        // second occurrence of the chosen tag is the same barcode, not a
        // competitor, so it must not trigger the delta rejection. (Regression:
        // the runner-up was once selected by excluding only the best
        // occurrence's index, so an equal-cost duplicate of the winning tag
        // killed an otherwise unambiguous match.)
        let mut g = group(&["ACGTACGT"]);
        g.mode = MatchMode::Nearest;
        g.delta = 1;
        let result = outcome(&g, "ACGTACGTACGTACGT");
        let best = result
            .best
            .expect("a repeated single tag is unambiguous: {result:?}");
        assert_eq!(best.tag_idx, 0);
        assert_eq!(best.cost, 0);
        assert_eq!(best.start, 0, "leftmost occurrence wins: {result:?}");
    }

    #[test]
    fn test_satisfies_finds_group_bounds() {
        // The tag occurs twice; min/maxFindsPerGroup gate the total.
        let mut g = group(&["ACGT"]);
        g.min_finds_per_group = Some(1);
        g.max_finds_per_group = Some(1);
        let result = outcome(&g, "ACGTACGT");
        assert!(result.accepted.len() >= 2, "two occurrences: {result:?}");
        assert!(
            !g.satisfies_finds(&result),
            "two matches exceed maxFindsPerGroup=1"
        );

        let mut single = group(&["ACGT"]);
        single.min_finds_per_group = Some(1);
        single.max_finds_per_group = Some(1);
        let one = outcome(&single, "ACGTTTTT");
        assert!(single.satisfies_finds(&one), "exactly one match satisfies");
    }

    #[test]
    fn test_satisfies_finds_min_group_zero_matches() {
        let mut g = group(&["ACGT"]);
        g.min_finds_per_group = Some(1);
        let none = outcome(&g, "TTTTTTTT");
        assert!(
            !g.satisfies_finds(&none),
            "no match fails minFindsPerGroup=1"
        );
    }

    #[test]
    fn test_satisfies_finds_max_per_tag() {
        // maxFindsPerTag=1: a tag appearing twice violates.
        let mut g = group(&["ACGT"]);
        g.max_finds_per_tag = Some(1);
        let result = outcome(&g, "ACGTACGT");
        assert!(
            !g.satisfies_finds(&result),
            "a tag found twice exceeds maxFindsPerTag=1"
        );
    }

    /// A find bound by its attribute, and how it sets a group.
    type FindBound = (&'static str, fn(&mut CompiledGroup));

    /// Each find bound a group can have.
    fn find_bounds() -> [FindBound; 5] {
        [
            ("minFindsPerGroup=1", |g| g.min_finds_per_group = Some(1)),
            ("maxFindsPerGroup=1", |g| g.max_finds_per_group = Some(1)),
            ("minFindsPerTag=1", |g| g.min_finds_per_tag = Some(1)),
            ("maxFindsPerTag=1", |g| g.max_finds_per_tag = Some(1)),
            ("findOne", |g| {
                g.max_finds_per_tag = Some(1);
                g.min_finds_per_group = Some(1);
                g.max_finds_per_group = Some(1);
            }),
        ]
    }

    /// `group` with the find bound `bound` set.
    fn bounded(group: &CompiledGroup, bound: fn(&mut CompiledGroup)) -> CompiledGroup {
        let mut bounded = group.clone();
        bound(&mut bounded);
        bounded
    }

    /// A `mode=nearest` group and a read that tags a and b both match exactly, a twice.
    fn nearest_tie() -> (CompiledGroup, GroupOutcome) {
        let mut g = group(&["AAAACCCC", "GGGGTTTT"]);
        g.mode = MatchMode::Nearest;
        let tie = outcome(&g, "AAAACCCCGGGGTTTTAAAACCCC");
        assert_eq!((tie.accepted.len(), tie.best), (3, None), "{tie:?}");
        (g, tie)
    }

    #[test]
    fn test_find_failure_of_a_nearest_tie_is_the_same_under_each_bound() {
        let (g, tie) = nearest_tie();
        let failures = find_bounds().map(|(_, bound)| bounded(&g, bound).find_failure(&tie));
        assert!(failures[0].is_some(), "a tie fails minFindsPerGroup=1");
        for ((name, _), failure) in find_bounds().iter().zip(&failures) {
            assert_eq!(failure, &failures[0], "{name}");
        }
    }

    #[test]
    fn test_find_failure_of_a_nearest_tie_without_find_bounds_is_none() {
        let (g, tie) = nearest_tie();
        assert_eq!(g.find_failure(&tie), None);
        let unbounded = bounded(&g, |g| g.min_finds_per_group = Some(0));
        assert_eq!(unbounded.find_failure(&tie), None);
    }

    #[test]
    fn test_find_failure_counts_a_nearest_winner_as_one_find() {
        // Tag a matches exactly at [1,7) and [12,18), tag b one substitution away at [0,6) and [11,17).
        let mut g = group(&["AACCGG", "GAACCC"]);
        g.dist = substitutions(1);
        g.mode = MatchMode::Nearest;
        let win = outcome(&g, "GAACCGGTTTTGAACCGGTTTT");
        assert_eq!(win.best.map(|best| best.tag_idx), Some(0), "{win:?}");
        assert_eq!(win.accepted.len(), 4, "{win:?}");
        for (name, bound) in find_bounds() {
            assert_eq!(bounded(&g, bound).find_failure(&win), None, "{name}");
        }
        assert!(!bounded(&g, |g| g.min_finds_per_group = Some(2)).satisfies_finds(&win));
        assert!(!bounded(&g, |g| g.min_finds_per_tag = Some(2)).satisfies_finds(&win));
    }

    #[test]
    fn test_find_failure_counts_every_match_under_mode_all() {
        let g = group(&["AAAACCCC", "GGGGTTTT"]);
        let tie = outcome(&g, "AAAACCCCGGGGTTTTAAAACCCC");
        assert_eq!(tie.best.map(|best| best.tag_idx), Some(0), "{tie:?}");
        let passes = find_bounds().map(|(_, bound)| bounded(&g, bound).satisfies_finds(&tie));
        assert_eq!(passes, [true, false, true, false, false]);
        assert!(bounded(&g, |g| g.min_finds_per_group = Some(3)).satisfies_finds(&tie));
    }

    #[test]
    fn test_partial5_matches_tag_truncated_at_read_start() {
        // partial5=4:0.1: a tag whose 5' bases run off the read start still
        // matches; the extracted span is the matched region, not padded to the
        // full tag length.
        let mut g = group(&["ATCGATCG"]);
        g.partial5 = Some(Partial {
            min_match: 4,
            max_mismatch_freq: Some(0.1),
        });
        // The read begins with CGATCG = the tag's last 6 bases (its first 2,
        // AT, are truncated).
        let result = outcome(&g, "CGATCGTTTT");
        let best = result.best.expect("a 5'-truncated match");
        assert_eq!(best.start, 0);
        assert_eq!(best.end, 6, "matched 6 nt, not padded to 8");
    }

    #[test]
    fn test_partial3_matches_tag_truncated_at_read_end() {
        // partial3=4:0.1: a tag whose 3' bases run off the read end matches the
        // truncated 3' side.
        let mut g = group(&["ACGTACGT"]);
        g.partial3 = Some(Partial {
            min_match: 4,
            max_mismatch_freq: Some(0.1),
        });
        // The read ends with ACGTA = the tag's first 5 bases (its last 3, CGT,
        // are truncated).
        let result = outcome(&g, "TTTTTACGTA");
        let best = result.best.expect("a 3'-truncated match");
        assert_eq!(best.start, 5);
        assert_eq!(best.end, 10, "matched 5 nt at the read end");
    }

    #[test]
    fn test_partial3_only_rejects_5prime_truncation() {
        // With only partial3 set, a 5'-truncated occurrence (CGATCG at the read
        // start) is not a permitted truncation, and the tag's full form is
        // absent, so there is no match.
        let mut g = group(&["ATCGATCG"]);
        g.partial3 = Some(Partial {
            min_match: 4,
            max_mismatch_freq: Some(0.1),
        });
        let result = outcome(&g, "CGATCGTTTT");
        assert!(
            result.best.is_none(),
            "a 5' truncation is rejected when only partial3 is set: {result:?}"
        );
    }

    #[test]
    fn test_partial5_only_rejects_3prime_truncation() {
        // The reciprocal: with only partial5 set, a 3'-truncated occurrence
        // (ACGTA at the read end) is not permitted, so it does not match.
        let mut g = group(&["ACGTACGT"]);
        g.partial5 = Some(Partial {
            min_match: 4,
            max_mismatch_freq: Some(0.1),
        });
        let result = outcome(&g, "TTTTTACGTA");
        assert!(
            result.best.is_none(),
            "a 3' truncation is rejected when only partial5 is set: {result:?}"
        );
    }

    #[test]
    fn test_partial5_rejects_over_mismatched_truncation() {
        // partial5=4:0.1 tolerates at most 10% mismatches over the matched
        // region. A 6 bp truncation with one mismatch is ~16.7% > 10%, so it is
        // rejected.
        let mut g = group(&["ATCGATCG"]);
        g.partial5 = Some(Partial {
            min_match: 4,
            max_mismatch_freq: Some(0.1),
        });
        // The read starts with CTATCG: the tag's 5' two bases truncated, then
        // one mismatch (G->T) in the 6 matched bases.
        let result = outcome(&g, "CTATCGTTTT");
        assert!(
            result.best.is_none(),
            "a truncation exceeding the mismatch frequency is rejected: {result:?}"
        );
    }

    #[test]
    fn test_partial5_min_only_reuses_dist_mismatch_budget() {
        // partial5=4 with no `:freq` reuses the group's `dist` substitution
        // count. The read is 6 bp, shorter than the 8 bp tag, so only a
        // 5'-truncated placement fits: CTATCG aligns to the tag's last 6 bases
        // CGATCG with one mismatch. With dist=1 that mismatch is within budget.
        let mut g = group(&["ATCGATCG"]);
        g.dist = Dist {
            mismatch: 1,
            indel: 0,
            total: 1,
        };
        g.partial5 = Some(Partial {
            min_match: 4,
            max_mismatch_freq: None,
        });
        let result = outcome(&g, "CTATCG");
        let best = result
            .best
            .expect("a 5'-truncated match within the dist budget");
        assert_eq!(best.start, 0);
        assert_eq!(
            best.end, 6,
            "matched all 6 nt, one sub within the dist budget"
        );
    }

    #[test]
    fn test_partial5_min_only_requires_exact_overlap_without_dist() {
        // The same 6 bp read, but with no `:freq` and no `dist` (mismatch
        // budget 0): a bare partial5 then demands an exact overlap, so the
        // one-mismatch truncation is rejected.
        let mut g = group(&["ATCGATCG"]);
        g.partial5 = Some(Partial {
            min_match: 4,
            max_mismatch_freq: None,
        });
        let result = outcome(&g, "CTATCG");
        assert!(
            result.best.is_none(),
            "a mismatched truncation is rejected when no dist budget is set: {result:?}"
        );
    }

    #[test]
    fn test_match_group_at_overrides_loc() {
        // With no own `loc`, an explicit window (as a relative `next` would
        // compute) restricts the search: the tag is found only when the window
        // covers it.
        let g = group(&["ACGT"]);
        let read: &[u8] = b"GGGGACGTGGGG";
        let window = |start, end| Location {
            file: 0,
            start: Endpoint::FromStart(start),
            end: Some(Endpoint::FromStart(end)),
        };
        let inside = match_group_at(&g, &mut Scratch::new(), &[read], Some(&window(4, 8)));
        assert_eq!(inside.best.expect("match in window").start, 4);
        let outside = match_group_at(&g, &mut Scratch::new(), &[read], Some(&window(0, 4)));
        assert!(
            outside.best.is_none(),
            "window excludes the tag: {outside:?}"
        );
    }

    #[test]
    fn test_next_window_from_upstream_match() {
        let upstream = MatchSpan {
            file: 2,
            start: 0,
            end: 10,
        };
        let link = NextLink {
            group: "g2".to_string(),
            window: Some((3, 5)),
        };
        let window = next_window(&upstream, &link).expect("a relative window");
        assert_eq!(window.file, 2);
        assert_eq!(window.start, Endpoint::FromStart(13));
        assert_eq!(window.end, Some(Endpoint::FromStart(15)));
    }

    #[test]
    fn test_next_window_bare_link_is_none() {
        let upstream = MatchSpan {
            file: 0,
            start: 0,
            end: 8,
        };
        let link = NextLink {
            group: "g2".to_string(),
            window: None,
        };
        assert!(next_window(&upstream, &link).is_none());
    }

    #[test]
    fn test_a_window_past_a_reverse_strand_match_mirrors_the_reverse_complement() {
        let mut rng = Rng(0xa11e);
        let mut found = 0;
        for trial in 0..2000 {
            let len = 3 + rng.below(5);
            let tags = (0..1 + rng.below(3)).map(|_| rng.bases(len)).collect();
            let mut g = CompiledGroup::new("g", tags);
            g.dist = [
                substitutions(0),
                substitutions(1),
                Dist {
                    mismatch: 0,
                    indel: 1,
                    total: 1,
                },
                Dist {
                    mismatch: 1,
                    indel: 1,
                    total: 2,
                },
            ][trial % 4];
            match rng.below(4) {
                0 => g.anchor = Some(Anchor::FivePrime),
                1 if g.dist.indel == 0 => g.anchor = Some(Anchor::ThreePrime),
                2 => g.revcomp = true,
                _ => {}
            }
            let mut planted = g.tags[rng.below(g.tags.len())].clone();
            if rng.below(2) == 0 {
                crate::iupac::reverse_complement(&mut planted);
            }
            let (left, right) = (rng.below(25), rng.below(25));
            let read = [
                rng.sample(left, b"ACGTN"),
                planted,
                rng.sample(right, b"ACGTN"),
            ]
            .concat();
            let mut flipped = read.clone();
            crate::iupac::reverse_complement(&mut flipped);
            let start = rng.below(read.len() + 1);
            let end = (start + rng.below(6)).min(read.len());
            let lo = rng.below(12);
            let window = RelativeWindow {
                upstream: 0,
                lo,
                hi: lo + rng.below(16),
            };
            let upstream = MatchSpan {
                file: 0,
                start,
                end,
            };
            let mirrored = MatchSpan {
                file: 0,
                start: read.len() - end,
                end: read.len() - start,
            };
            let past_reverse =
                match_verdict_after(&g, &mut Scratch::new(), &[&read], &upstream, true, &window);
            let past_forward = match_verdict_after(
                &g,
                &mut Scratch::new(),
                &[&flipped],
                &mirrored,
                false,
                &window,
            );
            assert_eq!(past_reverse.failure, past_forward.failure, "{g:?}");
            assert_eq!(
                past_reverse.best,
                past_forward.best.map(|best| best.mirrored(read.len())),
                "{g:?} {} {upstream:?} {window:?}",
                String::from_utf8_lossy(&read)
            );
            found += usize::from(past_reverse.best.is_some());
        }
        assert!(found > 200, "{found} trials matched");
    }

    #[test]
    fn test_match_group_multi_file_loc_picks_right_file() {
        let mut g = group(&["ACGTACGT"]);
        g.loc = Some(Location {
            file: 1,
            start: Endpoint::FromStart(0),
            end: None,
        });
        let result = match_group(&g, &mut Scratch::new(), &[b"ACGTACGT", b"GGACGTACGTGG"]);
        let best = result.best.expect("a match in file 1");
        assert_eq!(best.file, 1);
        assert_eq!(best.start, 2);
    }

    /// Assert two outcomes are observably equal: the same `best` and the same
    /// set of accepted matches. The accepted-Vec order is not part of the
    /// contract (the batched path may interleave strands differently); only
    /// `best` and the per-tag/group counts are consumed downstream.
    fn assert_same_outcome(a: &GroupOutcome, b: &GroupOutcome) {
        assert_eq!(a.best, b.best, "best is identical");
        let key = |m: &TagMatch| {
            (
                m.tag_idx,
                m.file,
                m.start,
                m.end,
                m.cost,
                m.subs,
                m.indels,
                strand_rank(m.strand),
            )
        };
        let mut left = a.accepted.clone();
        let mut right = b.accepted.clone();
        left.sort_by_key(key);
        right.sort_by_key(key);
        assert_eq!(left, right, "same accepted set");
    }

    /// Ten distinct equal-length tags (more than the batched threshold).
    fn ten_tags() -> CompiledGroup {
        group(&[
            "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT", "ACACAC", "GTGTGT", "TGTGTG", "CACACA",
            "ACGTAC", "TGCATG",
        ])
    }

    #[test]
    fn test_encode_builds_batch_above_threshold_only() {
        // A >8 equal-length bucket is batched; an 8-tag group is not (per-tag).
        let mut big = ten_tags();
        big.encode();
        assert_eq!(big.encoded.len(), 1, "one length bucket, batched");
        assert_eq!(big.encoded[0].tag_indices.len(), 10);

        let mut small = group(&[
            "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT", "ACACAC", "GTGTGT", "TGTGTG", "CACACA",
        ]);
        small.encode();
        assert!(small.encoded.is_empty(), "8 tags stay per-tag");
    }

    #[test]
    fn test_scratch_reuse_matches_fresh_searcher() {
        // One Scratch reused across many sequential calls (and across searcher
        // variants) must give the same outcomes as a fresh Searcher each time:
        // the reused internal buffers are scratch, not state (warm-reuse
        // correctness).
        let mut fwd = group(&["ACGTACGT"]);
        fwd.dist = Dist {
            mismatch: 1,
            indel: 0,
            total: 1,
        };
        let mut rc = fwd.clone();
        rc.revcomp = true;
        let reads = ["TTTACGTACGTCC", "GGACGTACGTAA", "ACGTACGTACGT", "NNNNNNNN"];

        let mut shared = Scratch::new();
        for read in reads {
            let seg = [read.as_bytes()];
            // Forward then reverse-complement on the SAME scratch exercises two
            // searcher slots plus buffer reuse across reads.
            let warm_fwd = match_group(&fwd, &mut shared, &seg);
            let warm_rc = match_group(&rc, &mut shared, &seg);
            assert_same_outcome(&warm_fwd, &match_group(&fwd, &mut Scratch::new(), &seg));
            assert_same_outcome(&warm_rc, &match_group(&rc, &mut Scratch::new(), &seg));
        }
    }

    #[test]
    fn test_batched_path_matches_per_tag_results() {
        // The batched path (after encode) must produce identical results to
        // per-tag search.
        let mut g = ten_tags();
        g.dist = Dist {
            mismatch: 1,
            indel: 0,
            total: 1,
        };
        let read = "TTACGTACGG"; // ACGTAC (tag 8) sits at [2, 8)
        let per_tag = outcome(&g, read);
        assert!(g.encoded.is_empty(), "not yet encoded -> per-tag");

        g.encode();
        let batched = outcome(&g, read);
        assert_same_outcome(&per_tag, &batched);
        let best = batched.best.expect("a match");
        assert_eq!(
            g.tags[best.tag_idx], b"ACGTAC",
            "pattern_idx maps to the right tag"
        );
    }

    #[test]
    fn test_batched_path_parity_with_mixed_lengths() {
        // A forward group with a >8 batched 6-mer bucket plus three per-tag
        // 8-mers: the encoded and unencoded forms agree across several reads
        // (the batched bucket and the per-tag fallback bucket coexist).
        let mut base = group(&[
            "AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT", "ACACAC", "GTGTGT", "TGTGTG", "CACACA",
            "ACGTAC", "TGCATGCA", "GGGGTTTT", "AACCGGTT", // three 8-mers (per-tag bucket)
        ]);
        base.dist = Dist {
            mismatch: 1,
            indel: 0,
            total: 1,
        };
        let mut encoded = base.clone();
        encoded.encode();
        assert_eq!(
            encoded.encoded.len(),
            1,
            "only the >8 6-mer bucket is batched"
        );

        for read in [
            "TTACGTACGG",       // a 6-mer forward
            "TTGTACGTAA",       // rc of ACGTAC region
            "AAAAGGGGTTTTCCCC", // contains an 8-mer (per-tag bucket)
            "NNNNNNNNNN",       // no-calls only
        ] {
            assert_same_outcome(&outcome(&base, read), &outcome(&encoded, read));
        }
    }

    #[test]
    fn test_revcomp_group_is_not_batched() {
        // Reverse-complement groups never take the batched path: sassy 0.2.3
        // computes the rc strand differently between `search` (per-tag) and
        // `search_encoded_patterns` (batched), so an rc batch can
        // drop/add/shift a match and flip the demux outcome. `encode` must
        // leave rc groups empty so they fall back to the verified per-tag
        // `search()`.
        //
        // This is the regression for a found parity break: nine length-6 IUPAC
        // tags with an indel budget over an N-containing read. Before the fix
        // the batched path emitted an extra Forward+ReverseComplement duplicate
        // of one tag, so with max_finds_per_tag=1 / max_finds_per_group=3 the
        // per-tag outcome satisfied the bounds (read assigned) while the
        // batched outcome violated them (read mis-routed to unassigned). With
        // rc batching disabled the two paths are the same code, so the demux
        // decision is identical.
        let mut base = group(&[
            "RNKTNA", "RMCYRK", "KWRNTA", "NCGWCW", "MAGTYN", "TAAYAR", "TYNRCA", "YNKYNA",
            "TAWKKR",
        ]);
        base.dist = Dist {
            mismatch: 0,
            indel: 1,
            total: 1,
        };
        base.revcomp = true;
        base.max_finds_per_tag = Some(1);
        base.max_finds_per_group = Some(3);

        let mut encoded = base.clone();
        encoded.encode();
        assert!(
            encoded.encoded.is_empty(),
            "an rc group stays entirely on per-tag search"
        );

        let read = "TTNTTAAYARNCCNTCNGTTCGCNGANGNANCTCCANTT";
        let per_tag = outcome(&base, read);
        let batched = outcome(&encoded, read);
        assert_same_outcome(&per_tag, &batched);
        assert_eq!(
            base.satisfies_finds(&per_tag),
            encoded.satisfies_finds(&batched),
            "the demux decision is identical with and without encode()"
        );
    }
}
