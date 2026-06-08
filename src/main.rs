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
/// unmux reads multiple FASTX/SAM/BAM/CRAM inputs, identifies and extracts
/// technical sequences (barcodes, UMIs, adapters) with error tolerance using
/// an approximate matcher, and writes FASTX/SAM/BAM/CRAM files with preserved
/// per-record segment qualities, fanning a record pool out into per-sample, and
/// optionally per-sample and per-library, files in a single pass.
///
/// Mental model: scan the input and match each input record against tag
/// --group(s), then pull those bases into named --extract streams, then route
/// the record somewhere such as assigning it to a --sample, else unassigned,
/// or removed by --remove. Set streams as the primary record sequences with
/// --template and set --tag(s) to write streams to SAM tags. A "pool" is all
/// input records for one run of unmux.
///
/// Quick start:
///
///   unmux in.fq --out out.bam
///   unmux in.fq --extract r=0:0:9 --template r > out.fq
///   unmux r1.fq r2.fq --group bc=bc.tsv --sample s=bc::t01 --out %sample.bam
///
/// Notation (also see param docs below with expressive examples`):
///
///   file:start:end  0-based, half-open; `end`=record length; neg counts from
///                   the record end. The FIRST number is the input file index
///                   (0:0:8 = file 0, bases [0,8); 1:0:8 = file 1).
///   @grp            a group's matched span. @grp+off:len / @grp-off:len step
///                   past / before it (trailing number is a LENGTH, not end).
///   @grpA..@grpB    the region between two matched spans.
///   +               concat streams (cb+umi) or 'AND' samples (gA::a+gB::b).
///   ~               reverse-complement the stream (~cb, BC=~bc).
///   ,               list separator ('OR' a tag pool & attribute lists).
///   %XX             percent-escape a byte in a tag sep/qual-sep value or an
///                   output path (%20 space, %09 tab, %2C comma).
///   %pool           the pool ID (see --pool).
///   %sample         the sample ID (--out only)
///   %sub_sample     the sub_sample ID (@RG LB) (--out only).
///   %ordinal        1-based read ordinal (R%ordinal → R1, R2) (--out only).
///   %source         0-based input file idx (--unassigned or --remove only).
#[derive(Debug, Parser)]
#[command(author, version, color = clap::ColorChoice::Always, term_width = 80, verbatim_doc_comment, override_usage = "unmux [READS]... [OPTIONS]")]
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
    /// Mutually exclusive with --in. Inputs must agree on yes/no qualities
    /// (FASTA cannot blend with a quality-containing FASTQ/SAM/BAM/CRAM).
    ///
    ///   unmux r1.fq r2.fq i1.fq    files 0, 1, 2
    ///   unmux < reads.fq           file 0 from stdin (no args = stdin)
    #[arg(value_name = "READS", num_args = 0.., verbatim_doc_comment)]
    inputs_positional: Vec<PathBuf>,

    /// Identifier for the input pool; fills the symbol %pool.
    ///
    /// Optional, and defaults to the common stem of the input filenames.
    ///
    ///   --pool lib01   %pool → lib01
    #[arg(long = "pool", value_name = "ID", verbatim_doc_comment)]
    pool: Option<String>,

    /// Input files set with `N=PATH` for explicit 0-based file indices.
    ///
    /// Repeatable; indices must be unique and contiguous from 0 (flag
    /// order is free; a gap or non-zero start is an error). `PATH` may be
    /// `-` for stdin (at most once). Mutually exclusive with positional
    /// inputs.
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
    /// Attributes (NAME::key=val):
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
    ///   bc::anchor=5p            anchor tags' 5-prime base at loc.start
    ///   bc::anchor=3p            anchor tags' 3-prime base at loc.end
    ///   bc::match=i7+i5          match tags on joined --extract streams
    #[arg(long = "group", value_name = "SPEC", verbatim_doc_comment)]
    groups: Vec<String>,

    /// Make a named record segment into a stream `NAME=[SPEC]` (repeatable).
    ///
    /// NAME becomes a stream for --template, --tag, and group `match=`.
    /// An extracted stream carries both bases AND qualities.
    ///
    /// `file:start:end` is 0-based half-open; the trailing number is an
    /// END (= record length); negatives count from the record end. The
    /// anchored `@grp` forms take a LENGTH as the trailing number instead.
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
    /// Each name is an --extract stream (an input record is not a stream
    /// until extracted). Concatenate with `+`; one value per output record
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
    /// `TAG=STREAM[+STREAM]` binds record bases (join with `+`); `TAG::ATTRS`
    /// sets qual/sep/raw. A multi-stream tag joins sequences with `sep`
    /// (default: `-`) and qualities with `qual-sep` (default: a space).
    /// Defaults pre-exist for common tags: CB/CY CR/CY RX/QX BC/QT OX/BZ.
    ///
    ///   --tag RX=umi             UMI tag (auto quality tag QX)
    ///   --tag CB=bc1+bc2+bc3     join three barcode streams
    ///   --tag CB::sep=_          join sequences with _ not -
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

    /// Sample fan-out target `SAMPLE[::SUB_SAMPLE]=SELECTOR` (repeatable).
    ///
    /// SELECTOR is `group::id-or-seq[,...]` (comma = OR pool), a bare
    /// `group` (any of its tags), or several joined with `+` (AND across
    /// groups). SUB_SAMPLE → `LB`.
    ///
    /// Exclusive with --sample-sheet and --sample-from-group.
    ///
    ///   --sample s1=bc::t01           route tag t01 to sample s1
    ///   --sample s1=bc::t01,t02       OR pool: any listed tag
    ///   --sample s1=i7::a+i5::b       AND: needs both indices
    ///   --sample s1::lib9=bc::t01     sub_sample lib9 (→ LB)
    ///   --sample s1::%pool=bc::t01    sub_sample from pool ID
    #[arg(long = "sample", value_name = "SPEC", verbatim_doc_comment)]
    samples: Vec<String>,

    /// Input sample sheet in TSV format (the table form of --sample).
    ///
    /// Columns for `sample` (→ `SM`, required), optional `sub_sample`
    /// (→ `LB`), and one column per group (cell = a tag ID or sequence).
    /// Rows sharing a key OR; multiple group columns AND.
    ///
    /// Exclusive with --sample / --sample-from-group.
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
    /// The "just split by barcode" mode: each tag in GROUP becomes a
    /// sample with no --sample lines and no sheet to maintain. `SM`
    /// is the tag ID; an optional `sub_sample` column in the group's
    /// tag file sets `LB`. Records whose GROUP tag matches nothing are
    /// unassigned, as with --sample. Pair with --out %sample.bam
    /// to write one file per tag.
    ///
    /// Exclusive with --sample / --sample-sheet.
    ///
    ///   --sample-from-group bc   one sample per tag in group bc
    #[arg(long = "sample-from-group", value_name = "GROUP", conflicts_with_all = ["samples", "sample_sheet"], verbatim_doc_comment)]
    sample_from_group: Option<String>,

    /// Records to remove via `SEL[=PATH_PATTERN]` (repeatable).
    ///
    /// A record matching SEL (a `group` or `group::id`) is removed and
    /// tallied as `removed` (distinct from --unassigned). With
    /// `=PATH_PATTERN` the removed records are written to the output path
    /// otherwise they are simply ignored. SEL is required. The only
    /// placeholders allowed in PATH_PATTERN are `%pool` and `%source`
    /// (input file index).
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

    /// Output path for demuxed records.
    ///
    /// Format set by extension (FASTX/SAM/BAM/CRAM). `-` or `/dev/stdout`
    /// writes standard output in the input format. Read groups and SAM tags
    /// are SAM/BAM/CRAM-only; FASTX puts --tag values in the record-name
    /// comment. Missing parent dirs are created. Placeholders fan-out the
    /// pool into multiple files including %pool, %sample, %sub_sample, and
    /// %ordinal.
    ///
    ///   --out out.bam                   one file, all assigned records
    ///   --out %sample.bam               one file per sample
    ///   --out %sample.%sub_sample.bam   per sample and sub-sample
    ///   --out %pool.R%ordinal.fq.gz     per pool, per template record
    ///
    /// [default: /dev/stdout]
    #[arg(long = "out", value_name = "PATTERN", verbatim_doc_comment)]
    out: Option<String>,

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

    /// Per-record QC in a UTF-8 JSON tag value for how the record was routed.
    ///
    /// Bare --qc-tag uses tag `ZS`; --qc-tag=TAG names another 2-char tag
    /// which must start with X/Y/Z or contain a lowercase char. QC is
    /// available for assigned, unassigned, and removed records.
    ///
    ///   ...           no --qc-tag means no QC is written
    ///   --qc-tag      JSON in tag ZS (the default)
    ///   --qc-tag=ZQ   JSON in tag ZQ
    #[arg(long = "qc-tag", value_name = "TAG", num_args = 0..=1, default_missing_value = "ZS", verbatim_doc_comment)]
    qc_tag: Option<String>,

    /// Fail fast unless every tag of every group named by a
    /// --sample/--sample-sheet selector is claimed by some sample (no
    /// tag in a sampled group left unclaimed). Off by default; applies only
    /// when samples are declared.
    #[arg(long = "require-samples-explain-all-tags", verbatim_doc_comment)]
    require_samples_explain_all_tags: bool,

    /// Disable auto pair-detection. Treat each record as single-end, even
    /// when one input looks interleaved (mate-name pairing). Interleaving
    /// is otherwise auto-detected.
    #[arg(long, verbatim_doc_comment)]
    per_record: bool,

    /// Compression level 0-9 for BGZF (BAM) and gzip (FASTX.gz).
    /// CRAM uses its own codecs and ignores this.
    #[arg(long, value_name = "LEVEL", default_value_t = 5, value_parser = clap::value_parser!(u8).range(0..=9), verbatim_doc_comment)]
    compression: u8,

    /// Worker thread count.
    ///
    /// More threads speed up large inputs; 1 runs fully serially.
    #[arg(long, default_value_t = 8, value_parser = clap::value_parser!(u16).range(1..), verbatim_doc_comment)]
    threads: u16,
}

/// Main binary entrypoint.
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), Error> {
    let env = Env::default().default_filter_or("info");
    env_logger::Builder::from_env(env).init();

    let matches = Cli::command().term_width(80).get_matches();
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
