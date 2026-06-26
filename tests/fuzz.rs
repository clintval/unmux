//! Fuzz / property targets for unmux's text parsers.
//!
//! These run under plain `cargo test` on stable Rust: `bolero::check!` replays
//! the on-disk corpus and exercises a bounded set of random inputs (raise the
//! count with `BOLERO_RANDOM_ITERATIONS`). For coverage-guided fuzzing, install
//! `cargo-bolero` and run on a nightly toolchain, e.g. `cargo bolero test
//! fuzz_parse_span_structured --engine libfuzzer`.
//!
//! The contract under test is "no panics": these parsers take untrusted CLI /
//! sheet text and must return `Ok`/`Err`, never crash. They are the
//! highest-value surface to fuzz (the binary decoding is noodles' job, fuzzed
//! upstream); the input reader's own robustness is covered by the
//! malformed-FASTQ acceptance scenario in tests/features/formats.feature.

use std::collections::HashSet;
use std::fmt::Write as _;

use bolero::generator::TypeGenerator;

/// Raw-bytes span fuzzing: arbitrary strings, almost all rejected at the first
/// character. Cheap insurance that total garbage never panics. The structured
/// target below is what actually drives the deep grammar paths.
#[test]
fn fuzz_parse_span() {
    bolero::check!()
        .with_type::<String>()
        .for_each(|value: &String| {
            let _ = unmux::grammar::parse_span(value);
        });
}

/// The sample-sheet TSV parser (header detection, group columns, OR pools
/// across rows). Multi-line, column-indexed text from an untrusted file.
#[test]
fn fuzz_parse_sample_sheet() {
    let declared: HashSet<&str> = ["bc", "grp", "cbt"].into_iter().collect();
    bolero::check!()
        .with_type::<String>()
        .for_each(|text: &String| {
            let _ = unmux::parse_sample_sheet(text, &declared);
        });
}

// Structure-aware span generator.
//
// Random strings hit `parse_span`'s first character and bounce; without
// coverage guidance (which needs a nightly libfuzzer run) they rarely reach the
// parser's interior. This generator emits grammar-SHAPED spans instead, so
// every input is syntactically a span and the run explores the semantic edge
// space: huge / overflowing offsets, `end` in the wrong slot, leftward anchors,
// between-anchor spans, and identifiers that are sometimes valid and sometimes
// not. It mirrors the grammar in `src/lib/grammar.rs::parse_span`.

/// A grammar-shaped span: `~`? then one of the four body forms.
#[derive(Debug, TypeGenerator)]
struct SpanInput {
    revcomp: bool,
    body: SpanBodyShape,
}

#[derive(Debug, TypeGenerator)]
enum SpanBodyShape {
    /// `file:start:end`
    File {
        file: Num,
        start: Endpoint,
        end: Endpoint,
    },
    /// `@group`
    AnchorMatch { group: Ident },
    /// `@group+offset:len` / `@group-offset:len`
    AnchorOffset {
        group: Ident,
        before: bool,
        offset: Num,
        len: Len,
    },
    /// `@from..@to`
    Between { from: Ident, to: Ident },
}

#[derive(Debug, TypeGenerator)]
enum Endpoint {
    End,
    FromStart(Num),
    FromEnd(Num),
}

#[derive(Debug, TypeGenerator)]
enum Len {
    End,
    Bases(Num),
}

/// Numbers span the small range, the full `u64` range, and a digit string past
/// `u64`/`usize` so the parser's integer-overflow error branch is reachable.
#[derive(Debug, TypeGenerator)]
enum Num {
    Value(u64),
    Overflow,
}

/// Identifiers are mostly drawn from a small valid pool (so the parser reaches
/// its deep accept paths), but sometimes arbitrary, to exercise
/// `validate_ident` and separator-confusion error branches.
#[derive(Debug, TypeGenerator)]
enum Ident {
    Pool(u8),
    Raw(String),
}

const IDENT_POOL: [&str; 6] = ["g", "grp", "bc", "cbt", "umi_1", "x"];

fn num_str(n: &Num) -> String {
    match n {
        Num::Value(v) => v.to_string(),
        Num::Overflow => "999999999999999999999999999999".to_string(),
    }
}

fn ident_str(i: &Ident) -> String {
    match i {
        Ident::Pool(n) => IDENT_POOL[(*n as usize) % IDENT_POOL.len()].to_string(),
        Ident::Raw(s) => s.clone(),
    }
}

fn endpoint_str(e: &Endpoint) -> String {
    match e {
        Endpoint::End => "end".to_string(),
        Endpoint::FromStart(n) => num_str(n),
        Endpoint::FromEnd(n) => format!("-{}", num_str(n)),
    }
}

fn len_str(l: &Len) -> String {
    match l {
        Len::End => "end".to_string(),
        Len::Bases(n) => num_str(n),
    }
}

/// Render a generated AST to its span string.
fn render(span: &SpanInput) -> String {
    let mut s = String::new();
    if span.revcomp {
        s.push('~');
    }
    match &span.body {
        SpanBodyShape::File { file, start, end } => {
            let _ = write!(
                s,
                "{}:{}:{}",
                num_str(file),
                endpoint_str(start),
                endpoint_str(end)
            );
        }
        SpanBodyShape::AnchorMatch { group } => {
            let _ = write!(s, "@{}", ident_str(group));
        }
        SpanBodyShape::AnchorOffset {
            group,
            before,
            offset,
            len,
        } => {
            let sep = if *before { '-' } else { '+' };
            let _ = write!(
                s,
                "@{}{}{}:{}",
                ident_str(group),
                sep,
                num_str(offset),
                len_str(len)
            );
        }
        SpanBodyShape::Between { from, to } => {
            let _ = write!(s, "@{}..@{}", ident_str(from), ident_str(to));
        }
    }
    s
}

#[test]
fn fuzz_parse_span_structured() {
    bolero::check!()
        .with_type::<SpanInput>()
        .for_each(|span: &SpanInput| {
            let _ = unmux::grammar::parse_span(&render(span));
        });
}
