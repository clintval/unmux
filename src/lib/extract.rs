//! Extraction: resolve an `--extract` span against a record's per-file
//! segments, carrying base qualities alongside the bases. Because we slice the
//! original per-base qualities by the matched coordinates (rather than
//! re-emitting bases through a separate writer), the qualities are never
//! mangled, so no separate quality-repair step is needed downstream.
//!
//! A span is one of: a `file:start:end` window (the trailing number is an end
//! position; negatives count from the record end; `end` is the record length),
//! the matched span of a group (`@grp`), a fixed length `offset` past a group's
//! matched end (`@grp+offset:len`), or the region between two matched anchors
//! (`@grpA..@grpB`). A span anchored on a reverse-strand match is measured in
//! the matched tag's own orientation and reverse-complemented into it, with its
//! qualities reversed, so a read and its reverse complement give the same
//! stream wherever the tag matches the same bases in both (the match itself is
//! chosen in read coordinates); a `~` at the point of use (`--tag T=~stream` /
//! `--template ~stream`) reverse-complements a stream again. Anchored spans
//! that reference a group which did not match resolve to no stream (`None`),
//! so the caller can route the record.

use std::collections::HashMap;

use anyhow::Result;

use crate::grammar::{AnchorLen, Span, SpanBody};

/// Resolve an anchor group name to its match span and whether that match is on
/// the reverse strand, looking the group up by name in the canonical `name ->
/// index` map and reading that slot of the positional hits. Returns `None` when
/// the group is unknown, did not match, or is a `match=` group (which records a
/// tag but anchors no span).
fn anchor_span(
    group: &str,
    hits: &GroupHits,
    index: &HashMap<String, usize>,
) -> Option<(MatchSpan, bool)> {
    let hit = hits.get(*index.get(group)?)?.as_ref()?;
    Some((hit.span?, hit.revcomp))
}

/// A record's bases and (optional) per-base qualities for one input file. FASTA
/// inputs have no qualities, so `quals` is `None`.
#[derive(Debug, Clone, Copy)]
pub struct Segment<'a> {
    /// The bases.
    pub bases: &'a [u8],
    /// The per-base qualities, when the input format carries them.
    pub quals: Option<&'a [u8]>,
}

/// The resolved location of a group's chosen match, used to anchor extraction
/// spans. Coordinates are forward-strand offsets into the named input file's
/// bases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatchSpan {
    /// The 0-based input file the match lies in.
    pub file: usize,
    /// 0-based start offset.
    pub start: usize,
    /// 0-based exclusive end offset.
    pub end: usize,
}

/// What matching resolved for one group on one record, stored positionally by
/// the group's canonical index (its position in the engine's group slice). The
/// presence of a `GroupHit` in a slot IS the match predicate: an absent slot
/// (`None`) means the group did not match.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroupHit {
    /// Index of the chosen tag within the group (definition order).
    pub tag_idx: usize,
    /// The match span, used to anchor `@grp` extraction. `Some` only for a
    /// normal record-window group; a `match=` group hits but anchors nothing,
    /// so its span is `None`.
    pub span: Option<MatchSpan>,
    /// Substitutions in the chosen match (carried for the `--qc-tag` provenance
    /// slug).
    pub subs: usize,
    /// Insertions plus deletions in the chosen match.
    pub indels: usize,
    /// Whether the chosen match was on the reverse-complement strand.
    pub revcomp: bool,
}

/// A record's per-group matches, indexed by canonical group position. A slot is
/// `Some` only for a group that matched; the anchor for `@grp` extraction is
/// its `span`.
pub type GroupHits = [Option<GroupHit>];

/// The error-corrected form of a matched `@grp` stream: the canonical tag bases
/// with the observed qualities length-fitted. Carried alongside the observed
/// bases so a `--tag`/`--template` can emit the corrected sequence (the
/// default) or, with `raw=true`, the observed one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Corrected {
    /// The canonical (corrected) bases.
    pub bases: Vec<u8>,
    /// The observed qualities, fitted to the corrected length.
    pub quals: Option<Vec<u8>>,
}

/// Where an extracted stream is cut from: bases `[start, end)` of input file `file`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StreamSource {
    /// The 0-based input file.
    pub file: usize,
    /// The 0-based start offset in the file's read.
    pub start: usize,
    /// The 0-based exclusive end offset in the file's read.
    pub end: usize,
    /// Whether the stream's bases were reverse-complemented from the read's orientation.
    pub reversed: bool,
}

/// An extracted stream: the observed bases and (optional) qualities, both
/// already trimmed and, for a span measured from a reverse-strand match,
/// reverse-complemented into the matched tag's orientation, plus an optional
/// error-corrected form and where in the input it is cut from. `bases`/`quals`
/// always hold the observed, uncorrected values; `corrected` is populated only
/// for a matched keeplist `@grp` self-span, and a `--tag`/`--template` consults
/// it unless `raw=true` selects the observed bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Extracted {
    /// The observed (uncorrected) bases, in the matched tag's orientation for a
    /// span measured from a reverse-strand match.
    pub bases: Vec<u8>,
    /// The observed qualities, in the same order as `bases`, when the segment
    /// carried them.
    pub quals: Option<Vec<u8>>,
    /// The error-corrected form, populated only for a matched keeplist `@grp`
    /// self-span; `None` otherwise.
    pub corrected: Option<Corrected>,
    /// Where the stream is cut from; `None` for text or a stream with no region.
    pub source: Option<StreamSource>,
}

impl Extracted {
    /// A column stream's text, which has no qualities and no corrected form.
    pub fn from_text(text: &str) -> Self {
        Self {
            bases: text.as_bytes().to_vec(),
            quals: None,
            corrected: None,
            source: None,
        }
    }

    /// The bases a SAM tag should emit: the corrected form when present and not
    /// raw-overridden, else the observed bases.
    pub fn tag_bases(&self, raw: bool) -> &[u8] {
        match &self.corrected {
            Some(corrected) if !raw => &corrected.bases,
            _ => &self.bases,
        }
    }

    /// The qualities a SAM tag should emit, mirroring
    /// [`Extracted::tag_bases`]'s corrected/observed choice.
    pub fn tag_quals(&self, raw: bool) -> Option<&[u8]> {
        match &self.corrected {
            Some(corrected) if !raw => corrected.quals.as_deref(),
            _ => self.quals.as_deref(),
        }
    }
}

/// Resolve a span into an extracted stream.
///
/// Returns `Ok(None)` when an anchored span references a group that did not
/// match (the stream is absent); `Ok(Some(_))` otherwise, including a
/// possibly-empty stream for a degenerate window or an empty between-anchors
/// region. A column span has no source range, so it also resolves to `Ok(None)`
/// here; the engine looks its text up by the matched tag.
pub fn extract(
    span: &Span,
    segments: &[Segment],
    hits: &GroupHits,
    index: &HashMap<String, usize>,
) -> Result<Option<Extracted>> {
    let resolved = match resolve_range(&span.body, segments, hits, index)? {
        SpanRange::Bounded {
            file,
            start,
            end,
            reverse,
        } => Some(slice(
            segment_at(segments, file)?,
            file,
            start,
            end,
            reverse,
        )),
        // No region between the anchors: the stream is empty but present.
        SpanRange::EmptyPresent => Some(Extracted {
            bases: Vec::new(),
            quals: segments.first().and_then(|s| s.quals.map(|_| Vec::new())),
            corrected: None,
            source: None,
        }),
        // An `@grp` whose group did not match: the stream is absent.
        SpanRange::Absent => None,
    };

    Ok(resolved)
}

/// The source byte range `(file, start, end)` a span draws from, or `None` when
/// it resolves to no concrete range (an unmatched `@grp` anchor, or an empty
/// between-anchors region). Used to compute the `frac_bases_unextracted`
/// coverage.
pub fn span_range(
    span: &Span,
    segments: &[Segment],
    hits: &GroupHits,
    index: &HashMap<String, usize>,
) -> Result<Option<(usize, usize, usize)>> {
    Ok(match resolve_range(&span.body, segments, hits, index)? {
        SpanRange::Bounded {
            file, start, end, ..
        } => Some((file, start, end)),
        SpanRange::EmptyPresent | SpanRange::Absent => None,
    })
}

/// Where a span resolves: a concrete byte range, an empty-but-present stream (a
/// between-anchors span that bounds no region), or absent (an `@grp` whose
/// group did not match).
enum SpanRange {
    /// A concrete `[start, end)` byte range in input file `file`.
    Bounded {
        /// Input file index.
        file: usize,
        /// 0-based start offset.
        start: usize,
        /// 0-based exclusive end offset.
        end: usize,
        /// Whether the range is read reverse-complemented, in the orientation
        /// of the reverse-strand match it is measured from.
        reverse: bool,
    },
    /// The stream exists but is empty: a between-anchors span with an unmatched
    /// anchor, or whose anchors overlap, are out of order, or lie in different
    /// input files.
    EmptyPresent,
    /// The stream is absent (an `@grp` whose group did not match).
    Absent,
}

/// Resolve a span body to its source byte range, without slicing. Shared by
/// [`extract`] (which slices) and [`span_range`] (which reports coverage for
/// metrics).
fn resolve_range(
    body: &SpanBody,
    segments: &[Segment],
    hits: &GroupHits,
    index: &HashMap<String, usize>,
) -> Result<SpanRange> {
    Ok(match body {
        SpanBody::File { file, start, end } => {
            let segment = segment_at(segments, *file)?;
            let len = segment.bases.len();
            let start = start.resolve(len);
            let end = end.resolve(len).max(start);
            SpanRange::Bounded {
                file: *file,
                start,
                end,
                reverse: false,
            }
        }
        SpanBody::AnchorMatch { group } => match anchor_span(group, hits, index) {
            Some((anchor, reverse)) => SpanRange::Bounded {
                file: anchor.file,
                start: anchor.start,
                end: anchor.end,
                reverse,
            },
            None => SpanRange::Absent,
        },
        SpanBody::Column { .. } | SpanBody::SampleColumn { .. } => SpanRange::Absent,
        SpanBody::AnchorOffset {
            group,
            before,
            offset,
            len,
        } => match anchor_span(group, hits, index) {
            Some((anchor, reverse)) => {
                let total = segment_at(segments, anchor.file)?.bases.len();
                // In a reverse-strand match's orientation, past it is toward the read start.
                let (start, end) = if *before != reverse {
                    // Leftward: the span ends `offset` bases before the matched
                    // start and extends `len` bases back, saturating at the
                    // read start (0); `:end` reaches the start.
                    let end = anchor.start.saturating_sub(*offset);
                    let start = match len {
                        AnchorLen::Bases(n) => end.saturating_sub(*n),
                        AnchorLen::ToEnd => 0,
                    };
                    (start, end)
                } else {
                    // Rightward: `len` bases starting `offset` past the matched
                    // end, clamped to the read end.
                    let start = (anchor.end + offset).min(total);
                    let end = match len {
                        AnchorLen::Bases(n) => (start + n).min(total),
                        AnchorLen::ToEnd => total,
                    };
                    (start, end)
                };
                SpanRange::Bounded {
                    file: anchor.file,
                    start,
                    end,
                    reverse,
                }
            }
            None => SpanRange::Absent,
        },
        SpanBody::Between { from, to } => {
            match (anchor_span(from, hits, index), anchor_span(to, hits, index)) {
                (Some((a, false)), Some((b, _))) if a.file == b.file && a.end <= b.start => {
                    SpanRange::Bounded {
                        file: a.file,
                        start: a.end,
                        end: b.start,
                        reverse: false,
                    }
                }
                // In a reverse-strand `from`'s orientation, `to` follows toward the read start.
                (Some((a, true)), Some((b, _))) if a.file == b.file && b.end <= a.start => {
                    SpanRange::Bounded {
                        file: a.file,
                        start: b.end,
                        end: a.start,
                        reverse: true,
                    }
                }
                _ => SpanRange::EmptyPresent,
            }
        }
    })
}

/// Look up a segment by file index.
fn segment_at<'a>(segments: &'a [Segment<'a>], file: usize) -> Result<Segment<'a>> {
    segments.get(file).copied().ok_or_else(|| {
        anyhow::anyhow!("extraction references input file {file}, which has no read")
    })
}

/// Slice `[start, end)` of input file `file`'s segment, carrying qualities when present, and
/// reverse-complement it (reversing the qualities) when `reverse` is set.
fn slice(segment: Segment, file: usize, start: usize, end: usize, reverse: bool) -> Extracted {
    let mut bases = segment.bases[start..end].to_vec();
    let mut quals = segment.quals.map(|q| q[start..end].to_vec());
    if reverse {
        crate::iupac::reverse_complement(&mut bases);
        if let Some(quals) = quals.as_mut() {
            quals.reverse();
        }
    }
    Extracted {
        bases,
        quals,
        corrected: None,
        source: Some(StreamSource {
            file,
            start,
            end,
            reversed: false,
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::{parse_span, Endpoint};

    /// Parse an extraction span for a test, panicking on a malformed span.
    fn sp(value: &str) -> Span {
        parse_span(value).unwrap()
    }

    /// One single-file record with bases and matching ASCII qualities.
    fn segs<'a>(bases: &'a [u8], quals: &'a [u8]) -> Vec<Segment<'a>> {
        vec![Segment {
            bases,
            quals: Some(quals),
        }]
    }

    /// Build the positional hits plus the `name -> index` map for a set of
    /// forward-strand `(group, span)` anchors, assigning each named group a
    /// canonical index in first-seen order.
    fn anchors(pairs: &[(&str, MatchSpan)]) -> (Vec<Option<GroupHit>>, HashMap<String, usize>) {
        let pairs: Vec<_> = pairs
            .iter()
            .map(|&(name, span)| (name, span, false))
            .collect();
        stranded(&pairs)
    }

    /// The same as [`anchors`] for `(group, span, revcomp)` anchors on either strand.
    fn stranded(
        pairs: &[(&str, MatchSpan, bool)],
    ) -> (Vec<Option<GroupHit>>, HashMap<String, usize>) {
        let mut index = HashMap::new();
        let mut hits: Vec<Option<GroupHit>> = Vec::new();
        for &(name, span, revcomp) in pairs {
            let idx = *index.entry(name.to_string()).or_insert_with(|| {
                hits.push(None);
                hits.len() - 1
            });
            hits[idx] = Some(GroupHit {
                tag_idx: 0,
                span: Some(span),
                subs: 0,
                indels: 0,
                revcomp,
            });
        }
        (hits, index)
    }

    /// A span in input file 0.
    fn at(start: usize, end: usize) -> MatchSpan {
        MatchSpan {
            file: 0,
            start,
            end,
        }
    }

    #[test]
    fn test_extract_file_window_carries_qualities() {
        // bases[2..6] with matching qualities.
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = Span {
            body: SpanBody::File {
                file: 0,
                start: Endpoint::FromStart(2),
                end: Endpoint::FromStart(6),
            },
        };
        let got = extract(&span, &segments, &[], &HashMap::new())
            .unwrap()
            .unwrap();
        assert_eq!(got.bases, b"CCGG");
        assert_eq!(got.quals.unwrap(), b"2345");
    }

    #[test]
    fn test_extract_file_to_end() {
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("0:8:end");
        let got = extract(&span, &segments, &[], &HashMap::new())
            .unwrap()
            .unwrap();
        assert_eq!(got.bases, b"AC");
        assert_eq!(got.quals.unwrap(), b"89");
    }

    #[test]
    fn test_extract_file_negative_end() {
        // 0:5:-2 keeps bases [5, len-2) of a length-10 read -> [5, 8).
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("0:5:-2");
        let got = extract(&span, &segments, &[], &HashMap::new())
            .unwrap()
            .unwrap();
        assert_eq!(got.bases, b"GTT");
        assert_eq!(got.quals.unwrap(), b"567");
    }

    #[test]
    fn test_extract_anchor_match_span() {
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("@grp");
        let (hits, index) = anchors(&[(
            "grp",
            MatchSpan {
                file: 0,
                start: 2,
                end: 6,
            },
        )]);
        let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
        assert_eq!(got.bases, b"CCGG");
    }

    #[test]
    fn test_extract_anchor_offset_fixed_length() {
        // @grp+2:3 -> 3 bases starting 2 past the match end (end=6) -> [8, 11)
        // clamped to [8, 10).
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("@grp+2:3");
        let (hits, index) = anchors(&[(
            "grp",
            MatchSpan {
                file: 0,
                start: 0,
                end: 6,
            },
        )]);
        let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
        assert_eq!(got.bases, b"AC");
    }

    #[test]
    fn test_extract_anchor_offset_to_end() {
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("@grp+0:end");
        let (hits, index) = anchors(&[(
            "grp",
            MatchSpan {
                file: 0,
                start: 0,
                end: 8,
            },
        )]);
        let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
        assert_eq!(got.bases, b"AC");
    }

    #[test]
    fn test_extract_anchor_unmatched_group_is_none() {
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("@grp");
        let got = extract(&span, &segments, &[], &HashMap::new()).unwrap();
        assert!(got.is_none());
    }

    #[test]
    fn test_extract_between_anchors() {
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("@a..@b");
        let (hits, index) = anchors(&[
            (
                "a",
                MatchSpan {
                    file: 0,
                    start: 0,
                    end: 2,
                },
            ),
            (
                "b",
                MatchSpan {
                    file: 0,
                    start: 6,
                    end: 8,
                },
            ),
        ]);
        let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
        // region [a.end=2, b.start=6) -> "CCGG".
        assert_eq!(got.bases, b"CCGG");
    }

    #[test]
    fn test_extract_between_adjacent_anchors_is_empty() {
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("@a..@b");
        let (hits, index) = anchors(&[
            (
                "a",
                MatchSpan {
                    file: 0,
                    start: 0,
                    end: 2,
                },
            ),
            (
                "b",
                MatchSpan {
                    file: 0,
                    start: 2,
                    end: 4,
                },
            ),
        ]);
        let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
        assert!(got.bases.is_empty());
        assert_eq!(got.quals.unwrap(), b"");
    }

    #[test]
    fn test_extract_between_missing_anchor_is_empty() {
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("@a..@b");
        let (hits, index) = anchors(&[(
            "a",
            MatchSpan {
                file: 0,
                start: 0,
                end: 2,
            },
        )]);
        let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
        assert!(got.bases.is_empty());
    }

    #[test]
    fn test_extract_between_out_of_order_is_empty() {
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("@a..@b");
        // Anchor `b` nested in `a`, wholly before it, and overlapping its end.
        for (a, b) in [((0, 8), (2, 4)), ((6, 8), (0, 2)), ((2, 6), (4, 8))] {
            let (hits, index) = anchors(&[
                (
                    "a",
                    MatchSpan {
                        file: 0,
                        start: a.0,
                        end: a.1,
                    },
                ),
                (
                    "b",
                    MatchSpan {
                        file: 0,
                        start: b.0,
                        end: b.1,
                    },
                ),
            ]);
            let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
            assert!(got.bases.is_empty(), "a={a:?} b={b:?}");
            assert_eq!(got.quals.unwrap(), b"", "a={a:?} b={b:?}");
            assert_eq!(span_range(&span, &segments, &hits, &index).unwrap(), None);
        }
    }

    #[test]
    fn test_extract_between_different_files_is_empty() {
        let segments = [
            Segment {
                bases: b"AACCGGTTAC",
                quals: Some(b"0123456789"),
            },
            Segment {
                bases: b"TTGGCCAATG",
                quals: Some(b"abcdefghij"),
            },
        ];
        let span = sp("@a..@b");
        let (hits, index) = anchors(&[
            (
                "a",
                MatchSpan {
                    file: 0,
                    start: 0,
                    end: 2,
                },
            ),
            (
                "b",
                MatchSpan {
                    file: 1,
                    start: 6,
                    end: 8,
                },
            ),
        ]);
        let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
        assert!(got.bases.is_empty());
        assert_eq!(got.quals.unwrap(), b"");
        assert_eq!(span_range(&span, &segments, &hits, &index).unwrap(), None);
    }

    #[test]
    fn test_extract_fasta_has_no_qualities() {
        let segments = vec![Segment {
            bases: b"AACCGGTTAC",
            quals: None,
        }];
        let span = sp("0:0:4");
        let got = extract(&span, &segments, &[], &HashMap::new())
            .unwrap()
            .unwrap();
        assert_eq!(got.bases, b"AACC");
        assert!(got.quals.is_none());
    }

    #[test]
    fn test_extract_anchor_offset_uses_match_len_not_end() {
        // Confirm the trailing number is a LENGTH for anchored spans: @grp+0:2
        // -> 2 bases from end.
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("@grp+0:2");
        let (hits, index) = anchors(&[(
            "grp",
            MatchSpan {
                file: 0,
                start: 0,
                end: 4,
            },
        )]);
        let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
        assert_eq!(got.bases, b"GG");
        assert!(matches!(
            span.body,
            SpanBody::AnchorOffset {
                len: AnchorLen::Bases(2),
                ..
            }
        ));
    }

    #[test]
    fn test_extract_anchor_offset_before_clamps_at_read_start() {
        // The match sits at [4, 8); `@grp-offset:len` reads LEFTWARD, ending
        // `offset` before the start.
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let (hits, index) = anchors(&[(
            "grp",
            MatchSpan {
                file: 0,
                start: 4,
                end: 8,
            },
        )]);
        // @grp-0:3 -> the 3 bases ending at the match start (4): read[1, 4) =
        // "ACC".
        let near = extract(&sp("@grp-0:3"), &segments, &hits, &index)
            .unwrap()
            .unwrap();
        assert_eq!(near.bases, b"ACC");
        // @grp-0:9 underflows the read start, so it clamps to read[0, 4) =
        // "AACC".
        let clamped = extract(&sp("@grp-0:9"), &segments, &hits, &index)
            .unwrap()
            .unwrap();
        assert_eq!(clamped.bases, b"AACC");
    }

    #[test]
    fn test_column_span_has_no_source_range() {
        let (hits, index) = anchors(&[(
            "grp",
            MatchSpan {
                file: 0,
                start: 0,
                end: 4,
            },
        )]);
        let segments = segs(b"ACGTACGT", b"IIIIIIII");
        let span = sp("@grp::cell");
        assert_eq!(extract(&span, &segments, &hits, &index).unwrap(), None);
        assert_eq!(span_range(&span, &segments, &hits, &index).unwrap(), None);
    }

    #[test]
    fn test_extract_reverse_strand_match_is_in_the_tag_orientation() {
        // The reverse-strand match over ACCG is CGGT in the tag's orientation.
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let (hits, index) = stranded(&[("grp", at(1, 5), true)]);
        let got = extract(&sp("@grp"), &segments, &hits, &index)
            .unwrap()
            .unwrap();
        assert_eq!(got.bases, b"CGGT");
        assert_eq!(got.quals.unwrap(), b"4321");
    }

    #[test]
    fn test_extract_reverse_strand_offsets_step_in_the_tag_orientation() {
        // From a reverse-strand match, `+` steps toward the read start and `-` toward its end.
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let (hits, index) = stranded(&[("grp", at(4, 8), true)]);
        for (value, bases, quals, range) in [
            ("@grp+0:3", "GGT", "321", (1, 4)),
            ("@grp+1:2", "GT", "21", (1, 3)),
            ("@grp+2:9", "TT", "10", (0, 2)),
            ("@grp+0:end", "GGTT", "3210", (0, 4)),
            ("@grp-0:2", "GT", "98", (8, 10)),
            ("@grp-1:5", "G", "9", (9, 10)),
            ("@grp-0:end", "GT", "98", (8, 10)),
            ("@grp-3:1", "", "", (10, 10)),
        ] {
            let span = sp(value);
            let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
            assert_eq!(got.bases, bases.as_bytes(), "{value}");
            assert_eq!(got.quals.unwrap(), quals.as_bytes(), "{value}");
            assert_eq!(
                span_range(&span, &segments, &hits, &index).unwrap(),
                Some((0, range.0, range.1)),
                "{value}"
            );
        }
    }

    #[test]
    fn test_extract_between_is_in_the_orientation_of_the_first_anchor() {
        // [2,7) is CCGGT forward and ACCGG in the orientation of a reverse-strand `a`.
        let segments = segs(b"AACCGGTTAC", b"0123456789");
        let span = sp("@a..@b");
        for (a, a_reverse, b, b_reverse, bases, quals) in [
            (at(0, 2), false, at(7, 9), false, "CCGGT", "23456"),
            (at(0, 2), false, at(7, 9), true, "CCGGT", "23456"),
            (at(7, 9), true, at(0, 2), true, "ACCGG", "65432"),
            (at(7, 9), true, at(0, 2), false, "ACCGG", "65432"),
            (at(0, 2), true, at(7, 9), true, "", ""),
            (at(7, 9), false, at(0, 2), false, "", ""),
        ] {
            let (hits, index) = stranded(&[("a", a, a_reverse), ("b", b, b_reverse)]);
            let got = extract(&span, &segments, &hits, &index).unwrap().unwrap();
            let case = format!("a={a:?} {a_reverse} b={b:?} {b_reverse}");
            assert_eq!(got.bases, bases.as_bytes(), "{case}");
            assert_eq!(got.quals.unwrap(), quals.as_bytes(), "{case}");
        }
        let (hits, index) = stranded(&[("a", at(7, 9), true), ("b", at(0, 2), true)]);
        assert_eq!(
            span_range(&span, &segments, &hits, &index).unwrap(),
            Some((0, 2, 7))
        );
    }

    #[test]
    fn test_reverse_strand_extracts_equal_forward_extracts_of_the_reverse_complement() {
        let bases = b"AACGTTGCANCG";
        let quals = b"0123456789:;";
        let mut flipped_bases = bases.to_vec();
        crate::iupac::reverse_complement(&mut flipped_bases);
        let flipped_quals: Vec<u8> = quals.iter().rev().copied().collect();
        let read = segs(bases, quals);
        let flipped = segs(&flipped_bases, &flipped_quals);
        let len = bases.len();
        let check = |value: &str, pairs: &[(&str, MatchSpan, bool)]| {
            let (hits, index) = stranded(pairs);
            let mirrored: Vec<_> = pairs
                .iter()
                .map(|&(name, span, reverse)| {
                    (name, at(len - span.end, len - span.start), !reverse)
                })
                .collect();
            let (flipped_hits, flipped_index) = stranded(&mirrored);
            let span = sp(value);
            let unsourced = |extracted: Option<Extracted>| {
                extracted.map(|extracted| Extracted {
                    source: None,
                    ..extracted
                })
            };
            let mirror = |source: StreamSource| StreamSource {
                start: len - source.end,
                end: len - source.start,
                reversed: !source.reversed,
                ..source
            };
            let forward = extract(&span, &read, &hits, &index).unwrap();
            let backward = extract(&span, &flipped, &flipped_hits, &flipped_index).unwrap();
            if let (Some(forward), Some(backward)) = (&forward, &backward) {
                if !forward.bases.is_empty() {
                    assert_eq!(
                        forward.source.map(mirror),
                        backward.source,
                        "{value} {pairs:?}"
                    );
                }
            }
            assert_eq!(unsourced(forward), unsourced(backward), "{value} {pairs:?}");
        };
        let spans: Vec<MatchSpan> = (0..len)
            .flat_map(|start| (start + 1..=len).map(move |end| at(start, end)))
            .collect();
        for &a in &spans {
            for a_reverse in [false, true] {
                for value in [
                    "@a", "@a+0:3", "@a+2:4", "@a+1:end", "@a-0:3", "@a-2:4", "@a-1:end",
                ] {
                    check(value, &[("a", a, a_reverse)]);
                }
                for &b in &spans {
                    for b_reverse in [false, true] {
                        check("@a..@b", &[("a", a, a_reverse), ("b", b, b_reverse)]);
                    }
                }
            }
        }
    }
}
