//! Flexible record parsing and demultiplexing to FASTX/SAM/BAM/CRAM.
use std::path::PathBuf;
use std::process;

use anyhow::{Error, Result};
use clap::builder::styling::{AnsiColor, Effects, Style, Styles};
use clap::{CommandFactory, FromArgMatches, Parser};
use env_logger::Env;
use log::*;
use mimalloc::MiMalloc;

use unmux::DemuxArgs;

/// A fast general-purpose allocator for the whole binary; demux is
/// allocation-heavy in its hot loop (per-record segment buffers, tag joins,
/// output records), so the global allocator matters. Measured ~15-18% faster
/// than the system allocator on a dual-index demux at equal RSS.
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

pub(crate) const HEADER: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const USAGE: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const LITERAL: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const PLACEHOLDER: Style = AnsiColor::Cyan.on_default();
pub(crate) const ERROR: Style = AnsiColor::Red.on_default().effects(Effects::BOLD);
pub(crate) const VALID: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const INVALID: Style = AnsiColor::Yellow.on_default().effects(Effects::BOLD);

/// The tool's one-line description at the top of the help: primary color
/// (green), bold.
pub(crate) const TITLE: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
/// Indented code examples and backtick-wrapped terms in help text: tertiary
/// color (yellow). Like the title color, this is a named ANSI slot rather than
/// a hard-coded RGB value, so it follows the user's terminal palette.
pub(crate) const CODE: Style = AnsiColor::Yellow.on_default();
/// The license/attribution footer at the bottom of the help: magenta, a named
/// ANSI slot distinct from the other help colors, so it follows the terminal.
pub(crate) const FOOTER: Style = AnsiColor::Magenta.on_default();
/// The right-hand description column of an indented two-column table: a faded
/// gray (bright-black), so the left-hand code term (in [`CODE`]) stands out.
pub(crate) const FADED: Style = AnsiColor::BrightBlack.on_default();

/// Cargo's color style.
/// [source](https://github.com/crate-ci/clap-cargo/blob/master/src/style.rs)
pub(crate) const CARGO_STYLING: Styles = Styles::styled()
    .header(HEADER)
    .usage(USAGE)
    .literal(LITERAL)
    .placeholder(PLACEHOLDER)
    .error(ERROR)
    .valid(VALID)
    .invalid(INVALID);

/// Flexible record parsing and demultiplexing to FASTX/SAM/BAM/CRAM.
///
/// This tool reads multiple FASTX/SAM/BAM/CRAM inputs, identifies and extracts
/// technical sequences (barcodes, UMIs, adapters) with error tolerance using an
/// approximate matcher, and writes FASTX/SAM/BAM/CRAM data with preserved
/// per-record segment qualities, fanning a record pool out into per-sample, and
/// optionally per-sample and per-library, files in a single pass.
///
/// Mental model:
///
///  1. scan the inputs and match each tag `--group`(s) against each record
///  2. pull those matched bases into named `--extract` streams
///  3. route the record: assign to a `--sample`, else unassigned, or `--remove` it
///  4. set streams as the primary record sequences with `--template`
///  5. also set streams into SAM tags with `--tag`(s)
///  6. write records by fanning them to files by sample, sub-sample, ordinal
///
/// A "pool" is all input records for one run of unmux.
///
/// Quick start:
///
///   unmux in.fq --out out.bam   # simply converts FASTQ to uBAM
///   unmux in.fq --extract myslice=0:0:9 --template myslice > out.fq
///   unmux r1.fq r2.fq --group bc=bc.tsv --sample s=bc::t01 --out %sample.bam
///
/// Notation (also see param docs below with expressive examples):
///
///   file:start:end   0-based, half-open; `end` is record LENGTH; neg counts
///                    from the record END. The FIRST number is the input file
///                    index (0:0:8 = file 0, bases [0,8); 1:0:8 = file 1).
///   @grp             a group's matched span. `@grp+off:len` & `@grp-off:len`
///                    step past/before it (trailing number is a LENGTH!).
///   @grpA..@grpB     the region between two matched spans.
///   +                concat streams (`cb+umi`) or 'AND' samples (`gA::a+gB::b`).
///   ~                reverse-complement the stream (`~cb`, `BC=~bc`).
///   ,                list separator ('OR' a tag pool & attribute lists).
///   %XX              percent-escape a byte in a tag sep/qual-sep value or an
///                    output path (`%20`=space, `%09`=tab, `%2C`=comma).
///   %pool            the pool ID (see `--pool`).
///   %sample          the sample ID (`--out` only)
///   %sub_sample      the sub_sample ID (`@RG LB`) (`--out` only).
///   %ordinal         1-based read ordinal (R%ordinal → R1, R2) (`--out` only).
///   %source          0-based input file idx (`--unassigned` or `--remove` only).
#[derive(Debug, Parser)]
#[command(author, version, color = clap::ColorChoice::Always, verbatim_doc_comment, override_usage = "unmux [READS]... [OPTIONS]")]
#[clap(styles = CARGO_STYLING)]
struct Cli {
    #[command(flatten)]
    demux: DemuxCmd,
}

/// Arguments for the `demux` subcommand.
#[derive(Debug, Parser)]
#[command(about, rename_all = "kebab-case")]
struct DemuxCmd {
    /// Input files, positional (FASTX/SAM/BAM/CRAM; auto-detected).
    ///
    /// 0-based by order, the first file is index 0 (splitcode-style). With no
    /// inputs at all, file 0 defaults to stdin, so bare unmux reads stdin.
    /// Mutually exclusive with `--in`. Inputs must agree on yes/no qualities
    /// (FASTA cannot blend with a quality-containing FASTQ/SAM/BAM/CRAM).
    ///
    ///   unmux r1.fq r2.fq i1.fq   files 0, 1, 2
    ///   unmux < reads.fq          file 0 from stdin (no args = stdin)
    #[arg(value_name = "READS", num_args = 0.., verbatim_doc_comment)]
    inputs_positional: Vec<PathBuf>,

    /// Identifier for the input pool; fills the placeholder %pool.
    ///
    /// Optional, and defaults to the common stem of the input filenames.
    ///
    ///   --pool lib01   %pool placeholder is now set to 'lib01'
    #[arg(long = "pool", value_name = "ID", verbatim_doc_comment)]
    pool: Option<String>,

    /// Input files set with `N=PATH` for explicit 0-based file indices.
    ///
    /// Repeatable; indices must be unique and contiguous from 0 (flag order
    /// is free; a gap or non-zero start is an error). `PATH` may be '`-`' for
    /// stdin (at most once). Mutually exclusive with positional inputs.
    ///
    ///   --in 0=r1.fq.gz --in 1=r2.fq.gz   file 0 = r1, file 1 = r2
    ///   --in 0=- --in 1=r2.fq.gz          file 0 from stdin, file 1 = r2
    #[arg(
        long = "in",
        value_name = "N=PATH",
        conflicts_with = "inputs_positional",
        verbatim_doc_comment
    )]
    inputs: Vec<String>,

    /// Output path for demuxed records.
    ///
    /// Format set by extension (FASTX/SAM/BAM/CRAM). '`-`' or `/dev/stdout`
    /// writes standard output in the input format. Read groups and SAM tags
    /// are SAM/BAM/CRAM-only; FASTX puts `--tag` values in the record-name
    /// comment. Missing parent dirs are created. Placeholders fan-out the
    /// pool into multiple files including `%pool`, `%sample`, `%sub_sample`,
    /// and `%ordinal`.
    ///
    ///   --out out.bam                   one file, all assigned records
    ///   --out %sample.bam               one file per sample
    ///   --out %sample.%sub_sample.bam   per sample and sub-sample
    ///   --out %pool.R%ordinal.fq.gz     per pool, per template record
    ///
    /// [default: /dev/stdout]
    #[arg(long = "out", value_name = "PATTERN", verbatim_doc_comment)]
    out: Option<String>,

    /// Tag-group definition or attributes (repeatable; accumulates).
    ///
    /// `NAME=SOURCE` is the UMI/barcode/adapter set; `NAME::ATTRS` constrains
    /// matching. Tags may use IUPAC codes. A matched `@grp` span is
    /// error-corrected by default.
    ///
    /// Sources:
    ///
    ///   bc=tags.tsv        from a TSV file with `id` and `seq` columns
    ///   bc={AAC,ACG,TTG}   inline set, auto tag IDs
    ///   bc={a=AAC,b=ACG}   inline set, explicit tag IDs
    ///
    /// Attributes (`NAME::key=val`):
    ///
    ///   bc::loc=1:0:8            file 1 bases [0,8) (default: whole record)
    ///   bc::dist=1               allow 1 substitution
    ///   bc::dist=1:1:2           1 sub + 1 indel, total <= 2
    ///   bc::mode=nearest         keep best only if it beats runner-up...
    ///   bc::delta=2              ...by >= 2 (needs mode=nearest)
    ///   bc::next=bc2:0-4         bc2 follows, 0-4 bp past this match
    ///   bc::prev=bc1             require bc1 to have matched earlier
    ///   bc::minFindsPerGroup=1   group matches >= once (keeplist)
    ///   bc::maxFindsPerGroup=1   group matches <= once
    ///   bc::minFindsPerTag=1     per-tag bounds (also use maxFindsPerTag)
    ///   bc::findOne              exactly one match (unambiguous single tag)
    ///   bc::both_strands=true    match forward and reverse-complement
    ///   bc::partial5=3:0.1       5' truncation ok: >=3 bp, <=10% mismatches
    ///   bc::partial3=3:0.1       same, at the 3' end
    ///   bc::anchor=5p            anchor tags' 5-prime base at `loc.start`
    ///   bc::anchor=3p            anchor tags' 3-prime base at `loc.end`
    ///   bc::match=i7+i5          match tags on joined `--extract` streams
    #[arg(long = "group", value_name = "SPEC", verbatim_doc_comment)]
    groups: Vec<String>,

    /// Make a named record segment into a stream `NAME=[SPEC]` (repeatable).
    ///
    /// `NAME` becomes a stream for `--template`, `--tag`, and group `match=`. An
    /// extracted stream carries both bases AND qualities.
    ///
    /// `file:start:end` is 0-based half-open; the trailing number is an `END`
    /// (full record length); negatives count from the record end. The
    /// anchored `@grp` forms take a `LENGTH` as the trailing number instead.
    ///
    ///   r=0:9:end          file 0, base 9 to the end
    ///   bc=1:0:8           file 1, bases [0,8)
    ///   tail=0:-10:end     last 10 bases of file 0
    ///   mid=0:5:-2         file 0, base 5 to (length - 2)
    ///   cb=@grp            the group's own matched span
    ///   umi=@grp+19:9      9 bases, 19 past grp's match end
    ///   up=@grp-0:9        9 bases just left of grp's match
    ///   ins=@grpA..@grpB   region between two anchors
    #[arg(long = "extract", value_name = "SPEC", verbatim_doc_comment)]
    extracts: Vec<String>,

    /// Set which streams become the primary record sequences (repeatable).
    ///
    /// Each name is an `--extract` stream (an input record is not a stream
    /// until extracted). Concatenate with '`+`'; one value per output record
    /// (R1, R2, ...). Optional; with none the full input is the raw output.
    ///
    /// SAM/BAM/CRAM allows at most two ordinals; multi-FASTX may have more.
    /// `::raw=true` emits observed bases for a stream from a corrected `@grp`
    /// (default: corrected; no effect on a never-corrected stream).
    ///
    ///   --template cdna               one record = the whole cdna stream
    ///   --template cb+umi             concatenate two streams into one
    ///   --template r1 r2              two output records (a pair)
    ///   --template r1 --template r2   same, repeated-flag form
    ///   --template cb::raw=true       observed bases of a corrected stream
    ///   --template ~r1                reverse-complement the stream
    #[arg(long = "template", value_name = "SPEC", num_args = 1.., verbatim_doc_comment)]
    templates: Vec<String>,

    /// SAM tag binding or attributes (repeatable; accumulates).
    ///
    /// `TAG=STREAM[+STREAM]` binds record bases (join with '`+`'); `TAG::ATTRS`
    /// sets qual/sep/raw. A multi-stream tag joins sequences with `sep`
    /// (default: '`-`') and qualities with `qual-sep` (default: a space).
    /// Default qual tags pre-exist for: CB/CY CR/CY RX/QX BC/QT OX/BZ.
    ///
    ///   --tag RX=umi             UMI tag (auto quality tag `QX`)
    ///   --tag CB=bc1+bc2+bc3     join three barcode streams
    ///   --tag CB::sep=_          join sequences with '`_`' not '`-`'
    ///   --tag CB::qual=CY        name the paired quality tag
    ///   --tag CB::qual=none      emit no quality tag
    ///   --tag CB::qual-sep=%20   join qualities with a space
    ///   --tag CR=bc::raw=true    set the observed (uncorrected) bases
    ///   --tag BC=~bc             reverse-complement the stream
    #[arg(long = "tag", value_name = "SPEC", verbatim_doc_comment)]
    tags: Vec<String>,

    /// Shared `@RG` header fields for every output read group (repeatable).
    ///
    ///   --rg-tag PL=ILLUMINA CN=Acme   sequencing platform and center
    ///   --rg-tag PU=run1.lane1         platform unit
    #[arg(long = "rg-tag", value_name = "K=V", num_args = 1.., verbatim_doc_comment)]
    rg_tags: Vec<String>,

    /// Per-record QC in a UTF-8 JSON tag value for how the record was routed.
    ///
    /// Bare `--qc-tag` uses tag `ZS`; `--qc-tag=TAG` names another 2-char tag which
    /// must start with X/Y/Z or contain a lowercase char. QC is available for
    /// assigned, unassigned, and removed records.
    ///
    ///   ...           no `--qc-tag` means no QC is written
    ///   --qc-tag      JSON in tag `ZS` (the default)
    ///   --qc-tag=ZQ   JSON in tag `ZQ`
    #[arg(long = "qc-tag", value_name = "TAG", num_args = 0..=1, default_missing_value = "ZS", verbatim_doc_comment)]
    qc_tag: Option<String>,

    /// Sample fan-out target `SAMPLE[::SUB_SAMPLE]=SELECTOR` (repeatable).
    ///
    /// SELECTOR is `group::id-or-seq[,...]` (comma is OR pool), a bare `group`
    /// (any of its tags), or several joined with '`+`' (AND across groups).
    /// SUB_SAMPLE → `@RG LB`.
    ///
    /// Exclusive with `--sample-sheet` and `--sample-from-group`.
    ///
    ///   --sample s1=bc::t01          route tag t01 to sample s1
    ///   --sample s1=bc::t01,t02      OR pool: any listed tag
    ///   --sample s1=i7::a+i5::b      AND: needs both indices
    ///   --sample s1::lib9=bc::t01    sub_sample lib9 (→ `@RG LB`)
    ///   --sample s1::%pool=bc::t01   sub_sample from pool ID
    #[arg(long = "sample", value_name = "SPEC", verbatim_doc_comment)]
    samples: Vec<String>,

    /// Input sample sheet in TSV format (the table form of `--sample`).
    ///
    /// Columns for `sample` (→ `@RG SM`, required), optional `sub_sample` (→
    /// `@RG LB`), and one column per group (cell = a tag ID or sequence). Rows
    /// sharing a key OR; multiple group columns AND.
    ///
    /// Exclusive with `--sample` and `--sample-from-group`.
    ///
    ///   --sample-sheet samples.tsv
    #[arg(
        long = "sample-sheet",
        value_name = "FILE",
        conflicts_with_all = ["samples", "sample_from_group"],
        verbatim_doc_comment
    )]
    sample_sheet: Option<PathBuf>,

    /// Make every tag in GROUP its own sample, 1-to-1 (a shortcut).
    ///
    /// The "just split by barcode" mode: each tag in `GROUP` becomes a sample
    /// with no `--sample` lines and no sheet to maintain. `@RG SM` is the tag ID;
    /// an optional `sub_sample` column in the group's tag file sets `@RG LB`.
    /// Records whose `GROUP` tag matches nothing are unassigned, as with
    /// `--sample`. Pair with `--out %sample.bam` to write one file per tag.
    ///
    /// Exclusive with `--sample` and `--sample-sheet`.
    ///
    ///   --sample-from-group bc   one sample per tag in group bc
    #[arg(long = "sample-from-group", value_name = "GROUP", conflicts_with_all = ["samples", "sample_sheet"], verbatim_doc_comment)]
    sample_from_group: Option<String>,

    /// Output path pattern for records matching no sample.
    ///
    /// The only placeholders are `%pool` and `%source` (input file index);
    /// `%source` fans these to one file per input record. Without this flag,
    /// unassigned records are dropped. Unmux warns when unassigned % reaches
    /// >=20% of the pool.
    ///
    ///   --unassigned unmatched.%source.fq.gz   one file per input file
    ///   --unassigned unmatched.fa              all segments in one FASTA
    #[arg(long = "unassigned", value_name = "PATTERN", verbatim_doc_comment)]
    unassigned: Option<String>,

    /// Records to remove via `SEL[=PATH_PATTERN]` (repeatable).
    ///
    /// A record matching `SEL` (a `group` or `group::id`) is removed and tallied as
    /// removed (distinct from `--unassigned`). With `=PATH_PATTERN` the removed
    /// records are written to the output path otherwise they are simply
    /// ignored. `SEL` is required. The only placeholders allowed in
    /// `PATH_PATTERN` are `%pool` and `%source` (input file index).
    ///
    ///   --remove phiX                      drop records matching phiX group
    ///   --remove bc::t99                   drop one specific tag's records
    ///   --remove phiX=phiX.%source.fq.gz   ...and write them out
    #[arg(
        long = "remove",
        value_name = "SEL[=PATH_PATTERN]",
        verbatim_doc_comment
    )]
    remove: Vec<String>,

    /// Per-sample metrics TSV with one data row per fan-out target.
    ///
    /// Only `%pool` is a valid path placeholder.
    ///
    ///   --metrics-per-sample %pool.unmux.per_sample.tsv
    #[arg(
        long = "metrics-per-sample",
        value_name = "PATTERN",
        verbatim_doc_comment
    )]
    metrics_per_sample: Option<PathBuf>,

    /// Pool-level summary metrics TSV.
    ///
    /// Only `%pool` is a valid path placeholder.
    ///
    ///   --metrics-summary %pool.unmux.summary.tsv
    #[arg(long = "metrics-summary", value_name = "PATTERN", verbatim_doc_comment)]
    metrics_summary: Option<PathBuf>,

    /// Fail fast unless every tag of every group named by a `--sample` and
    /// `--sample-sheet` selector is claimed by some sample (no tag in
    /// a sampled group left unclaimed). Off by default; applies only when
    /// samples are declared.
    #[arg(long = "require-samples-explain-all-tags", verbatim_doc_comment)]
    require_samples_explain_all_tags: bool,

    /// Disable auto pair-detection. Treat each record as single-end, even
    /// when one input looks interleaved (mate-name pairing). Interleaving is
    /// otherwise auto-detected.
    #[arg(long, verbatim_doc_comment)]
    per_record: bool,

    /// Compression level 0-9 for BGZF (BAM) and gzip (FASTX.gz). CRAM uses
    /// its own codecs and ignores this.
    #[arg(long, value_name = "LEVEL", default_value_t = 5, value_parser = clap::value_parser!(u8).range(0..=9), verbatim_doc_comment)]
    compression: u8,

    /// Worker thread count.
    ///
    /// More threads speed up large inputs; 1 runs fully serially.
    #[arg(long, default_value_t = 8, value_parser = clap::value_parser!(u16).range(1..), verbatim_doc_comment)]
    threads: u16,
}

/// The ANSI escape that starts `style`, or an empty string when `color` is off
/// (honoring `NO_COLOR`). Its matching reset comes from [`reset`].
fn esc(style: Style, color: bool) -> String {
    if color {
        style.render().to_string()
    } else {
        String::new()
    }
}

/// The reset escape for `style` (a full SGR reset), or empty when `color` is off.
fn esc_reset(style: Style, color: bool) -> String {
    if color {
        style.render_reset().to_string()
    } else {
        String::new()
    }
}

/// Paint the first line of `text` in the primary [`TITLE`] color (the tool's
/// one-line description), leaving the rest of the text untouched. With `color`
/// off, the text is returned unchanged.
fn green_first_line(text: &str, color: bool) -> String {
    let title = esc(TITLE, color);
    let reset = esc_reset(TITLE, color);
    match text.split_once('\n') {
        Some((first, rest)) => format!("{title}{first}{reset}\n{rest}"),
        None => format!("{title}{text}{reset}"),
    }
}

/// Drop the backticks around terms (`` `like this` ``), painting each term with
/// the `code` escape. `after` is the escape that restores the surrounding text's
/// color once a term ends: a reset for default prose, or the code/faded color
/// again inside a code or description column. With color off, `code` and `after`
/// are empty, so this just strips the backticks.
fn paint_backtick_terms(text: &str, code: &str, after: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(open) = rest.find('`') {
        out.push_str(&rest[..open]);
        let tail = &rest[open + 1..];
        match tail.find('`') {
            Some(close) => {
                out.push_str(code);
                out.push_str(&tail[..close]);
                out.push_str(after);
                rest = &tail[close + 1..];
            }
            // An unmatched backtick has no closing partner; keep it verbatim.
            None => {
                out.push('`');
                rest = tail;
                break;
            }
        }
    }
    out.push_str(rest);
    out
}

/// Split an indented two-column row into (left, gap, right) at the first run of
/// three or more spaces that follows the left-hand term. The 3-space minimum is
/// the table delimiter; single-column example lines (no such gap) return `None`.
fn split_two_column(content: &str) -> Option<(&str, &str, &str)> {
    let bytes = content.as_bytes();
    let mut seen_term = false;
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b' ' {
            let start = i;
            while i < bytes.len() && bytes[i] == b' ' {
                i += 1;
            }
            if seen_term && i - start >= 3 {
                return Some((&content[..start], &content[start..i], &content[i..]));
            }
        } else {
            seen_term = true;
            i += 1;
        }
    }
    None
}

/// Style one help string, dropping backticks throughout:
///
/// - indented two-column rows: the left code term in [`CODE`] (yellow), the
///   right description column in [`FADED`] (gray);
/// - indented description continuations (deeper indent, no term): all [`FADED`];
/// - indented single-column examples (a bare command): all [`CODE`];
/// - prose: default color, with backtick terms painted in [`CODE`].
///
/// Backtick terms stay [`CODE`] everywhere, so code stands out even inside a
/// faded description.
fn style_help_text(text: &str, color: bool) -> String {
    let code = esc(CODE, color);
    let faded = esc(FADED, color);
    let reset = esc_reset(CODE, color);
    text.split_inclusive('\n')
        .map(|line| {
            let (content, newline) = match line.strip_suffix('\n') {
                Some(content) => (content, "\n"),
                None => (line, ""),
            };
            if !content.starts_with("  ") || content.trim().is_empty() {
                // Prose (or a blank line): only the backtick terms are painted.
                return format!("{}{newline}", paint_backtick_terms(content, &code, &reset));
            }
            if let Some((left, gap, right)) = split_two_column(content) {
                let left = paint_backtick_terms(left, &code, &code);
                let right = paint_backtick_terms(right, &code, &faded);
                format!("{code}{left}{reset}{gap}{faded}{right}{reset}{newline}")
            } else if content.len() - content.trim_start().len() > 2 {
                // A wrapped description continuation aligned past the left column.
                let painted = paint_backtick_terms(content, &code, &faded);
                format!("{faded}{painted}{reset}{newline}")
            } else {
                // A single-column code example (a bare command, no description).
                let painted = paint_backtick_terms(content, &code, &code);
                format!("{code}{painted}{reset}{newline}")
            }
        })
        .collect()
}

/// Style the command's help: paint the first line of the about text in the
/// primary [`TITLE`] color and apply [`style_help_text`] (code-color the
/// examples, paint backtick terms, drop the backticks). Applied to the short
/// about (drives `-h`), the long about (drives `--help`), and each option's
/// short and long help.
fn decorate_help(cmd: clap::Command, color: bool) -> clap::Command {
    let about = cmd.get_about().map(ToString::to_string);
    let long_about = cmd
        .get_long_about()
        .map(ToString::to_string)
        .or_else(|| about.clone());

    let footer = format!(
        "{}MIT License 2026 · Clint Valentine{}",
        esc(FOOTER, color),
        esc_reset(FOOTER, color),
    );

    let choice = if color {
        clap::ColorChoice::Always
    } else {
        clap::ColorChoice::Never
    };
    let mut cmd = cmd.color(choice);
    if let Some(about) = about {
        cmd = cmd.about(green_first_line(&style_help_text(&about, color), color));
    }
    if let Some(long_about) = long_about {
        cmd = cmd.long_about(green_first_line(
            &style_help_text(&long_about, color),
            color,
        ));
    }
    cmd = cmd
        .next_line_help(true)
        .after_help(footer.clone())
        .after_long_help(footer);
    cmd.mut_args(|mut arg| {
        if let Some(help) = arg.get_help().map(ToString::to_string) {
            arg = arg.help(style_help_text(&help, color));
        }
        if let Some(long_help) = arg.get_long_help().map(ToString::to_string) {
            arg = arg.long_help(style_help_text(&long_help, color));
        }
        arg
    })
}

/// Main binary entrypoint.
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), Error> {
    // NO_COLOR (https://no-color.org): if the variable is present at all, even
    // empty, suppress every ANSI color, both ours and clap's and the logger's.
    let color = std::env::var_os("NO_COLOR").is_none();

    let env = Env::default().default_filter_or("info");
    let write_style = if color {
        env_logger::WriteStyle::Auto
    } else {
        env_logger::WriteStyle::Never
    };
    env_logger::Builder::from_env(env)
        .write_style(write_style)
        .init();

    // Our doc comments are hand-wrapped to fit 80 columns (that is what
    // `verbatim_doc_comment` is for), and `decorate_help` injects ANSI color
    // into them. clap's own help wrapping would count those escape bytes as
    // visible width and mis-wrap the styled lines, so we disable it and let the
    // authored line breaks stand as the wrapping.
    let cmd = decorate_help(Cli::command().term_width(usize::MAX), color);
    let matches = cmd.get_matches();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());

    // `unmux` is the demux command: the flattened top-level args are the demux
    // invocation.
    let result = run_demux_cmd(cli.demux);

    match result {
        Ok(()) => process::exit(0),
        Err(e) => {
            error!("{e:#}");
            process::exit(1);
        }
    }
}

/// Resolve a parsed `demux` invocation into [`DemuxArgs`] and run it.
/// Positional inputs are rewritten to the `N=PATH` form so the engine sees one
/// input list regardless of entry style. With no inputs given at all, file 0
/// defaults to stdin (`0=-`), so a bare `unmux` is a stdin→stdout filter
/// (paired with --out's stdout default).
fn run_demux_cmd(cmd: DemuxCmd) -> Result<()> {
    let inputs = if !cmd.inputs_positional.is_empty() {
        cmd.inputs_positional
            .iter()
            .enumerate()
            .map(|(i, p)| format!("{i}={}", p.display()))
            .collect()
    } else if !cmd.inputs.is_empty() {
        cmd.inputs
    } else {
        vec!["0=-".to_string()]
    };
    unmux::run_demux(DemuxArgs {
        pool: cmd.pool,
        inputs,
        groups: cmd.groups,
        extracts: cmd.extracts,
        templates: cmd.templates,
        tags: cmd.tags,
        rg_tags: cmd.rg_tags,
        samples: cmd.samples,
        sample_sheet: cmd.sample_sheet,
        sample_from_group: cmd.sample_from_group,
        require_samples_explain_all_tags: cmd.require_samples_explain_all_tags,
        remove: cmd.remove,
        out: cmd.out,
        unassigned: cmd.unassigned,
        metrics_per_sample: cmd.metrics_per_sample,
        metrics_summary: cmd.metrics_summary,
        qc_tag: cmd.qc_tag,
        per_record: cmd.per_record,
        compression: cmd.compression,
        threads: cmd.threads as usize,
        command_line: Some(std::env::args().collect::<Vec<_>>().join(" ")),
    })
}
