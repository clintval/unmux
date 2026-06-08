//! Demux metrics: per-fan-out-target and pool-level record counts plus
//! `frac_bases_unextracted`, written as headered, concatenable TSVs
//! (`--metrics-per-sample` / `--metrics-summary`) and always logged to stderr.
//!
//! `frac_bases_unextracted` is the fraction of bases in retained records not
//! routed into a `--template` or `--tag` stream, pooled as `sum(unextracted) /
//! sum(denominator)` over the target's records (0 for a pure pass-through run,
//! where the whole record is the output).

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};
use csv::{Terminator, WriterBuilder};
use serde::Serialize;
use thousands::Separable;

use crate::fanout::Target;

/// A running per-target tally: records routed to it, and the pooled
/// `frac_bases_unextracted` numerator and denominator.
#[derive(Default, Clone)]
struct TargetTally {
    reads: u64,
    denominator_bases: u64,
    unextracted_bases: u64,
}

/// Accumulated demux metrics for one run.
pub struct Metrics {
    pool: String,
    targets: Vec<Target>,
    per_target: HashMap<Target, TargetTally>,
    total: u64,
    pass_through: u64,
    unassigned: u64,
    removed: u64,
    retained_denominator_bases: u64,
    retained_unextracted_bases: u64,
    /// Per-input-file (segment-index) record-length sum and count, for the mean
    /// record length QC column.
    input_len_sums: Vec<u64>,
    input_len_counts: Vec<u64>,
}

impl Metrics {
    /// A fresh accumulator for `pool` with one tally slot per fan-out `target`
    /// (in first-seen order).
    pub fn new(pool: impl Into<String>, targets: &[Target]) -> Self {
        let per_target = targets
            .iter()
            .cloned()
            .map(|target| (target, TargetTally::default()))
            .collect();
        Self {
            pool: pool.into(),
            targets: targets.to_vec(),
            per_target,
            total: 0,
            pass_through: 0,
            unassigned: 0,
            removed: 0,
            retained_denominator_bases: 0,
            retained_unextracted_bases: 0,
            input_len_sums: Vec::new(),
            input_len_counts: Vec::new(),
        }
    }

    /// Count one processed fragment (every record, before routing).
    pub fn record_processed(&mut self) {
        self.total += 1;
    }

    /// Accumulate one record's per-input-file (segment) base lengths, for the
    /// mean-record-length column. Called for every processed record (assigned,
    /// pass-through, unassigned, or removed), indexed by input order.
    pub fn record_read_lengths(&mut self, lengths: impl Iterator<Item = usize>) {
        for (input, len) in lengths.enumerate() {
            if input >= self.input_len_sums.len() {
                self.input_len_sums.resize(input + 1, 0);
                self.input_len_counts.resize(input + 1, 0);
            }
            self.input_len_sums[input] += len as u64;
            self.input_len_counts[input] += 1;
        }
    }

    /// The mean record length per input file as a comma-separated list (rounded
    /// to the nearest base), e.g. `117,145,61`; empty when no records were
    /// processed.
    fn mean_read_len_by_input(&self) -> String {
        self.input_len_sums
            .iter()
            .zip(&self.input_len_counts)
            // Rounded mean read length; checked_div yields None (-> 0) for an
            // unused input slot.
            .map(|(&sum, &count)| (sum + count / 2).checked_div(count).unwrap_or(0))
            .map(|mean| mean.to_string())
            .collect::<Vec<_>>()
            .join(",")
    }

    /// Count a record routed to `target`, with its `frac_bases_unextracted`
    /// contribution.
    pub fn record_assigned(&mut self, target: &Target, denominator: u64, unextracted: u64) {
        let tally = self.per_target.entry(target.clone()).or_default();
        tally.reads += 1;
        tally.denominator_bases += denominator;
        tally.unextracted_bases += unextracted;
        self.retain_bases(denominator, unextracted);
    }

    /// Count a passed-through record (no sample fan-out configured), with its
    /// bases contribution.
    pub fn record_pass_through(&mut self, denominator: u64, unextracted: u64) {
        self.pass_through += 1;
        self.retain_bases(denominator, unextracted);
    }

    /// Count a record that matched no sample.
    pub fn record_unassigned(&mut self) {
        self.unassigned += 1;
    }

    /// Count a record removed by a `--remove` rule.
    pub fn record_removed(&mut self) {
        self.removed += 1;
    }

    fn retain_bases(&mut self, denominator: u64, unextracted: u64) {
        self.retained_denominator_bases += denominator;
        self.retained_unextracted_bases += unextracted;
    }

    /// Total records routed to a fan-out target.
    fn assigned(&self) -> u64 {
        self.targets
            .iter()
            .map(|target| self.tally(target).reads)
            .sum()
    }

    /// The fraction of processed records that matched no sample (the `log`
    /// warning fires at `>= 0.20`).
    pub fn unassigned_rate(&self) -> f64 {
        frac(self.unassigned, self.total)
    }

    fn tally(&self, target: &Target) -> TargetTally {
        self.per_target.get(target).cloned().unwrap_or_default()
    }

    /// Write the per-sample TSV: leading `pool`/`sample`/`sub_sample`, then the
    /// wide metric columns, one row per fan-out target.
    pub fn write_per_sample_tsv(&self, path: &Path) -> Result<()> {
        let mut writer = tsv_writer(path)?;
        for target in &self.targets {
            let tally = self.tally(target);
            writer
                .serialize(PerSampleRow {
                    pool: &self.pool,
                    sample: &target.sample,
                    sub_sample: target.sub_sample.as_deref().unwrap_or(""),
                    reads: tally.reads,
                    frac_of_pool: frac(tally.reads, self.total),
                    frac_bases_unextracted: frac(tally.unextracted_bases, tally.denominator_bases),
                })
                .context("failed to write a per-sample metrics row")?;
        }
        writer.flush().context("failed to flush metrics")
    }

    /// Write the pool-level summary TSV: leading `pool`, then the wide metric
    /// columns, one row.
    pub fn write_summary_tsv(&self, path: &Path) -> Result<()> {
        let assigned = self.assigned();
        let mut writer = tsv_writer(path)?;
        writer
            .serialize(SummaryRow {
                pool: &self.pool,
                total_reads: self.total,
                reads_assigned: assigned,
                reads_pass_through: self.pass_through,
                reads_unassigned: self.unassigned,
                reads_removed: self.removed,
                frac_assigned: frac(assigned, self.total),
                frac_bases_unextracted: frac(
                    self.retained_unextracted_bases,
                    self.retained_denominator_bases,
                ),
                mean_read_len_by_input: self.mean_read_len_by_input(),
            })
            .context("failed to write the summary metrics row")?;
        writer.flush().context("failed to flush metrics")
    }

    /// Log the metrics to stderr as a long-form table (always emitted, with or
    /// without the TSV flags).
    pub fn log(&self) {
        let assigned = self.assigned();
        log::info!(
            "metrics pool={}: {} reads ({} assigned, {} pass-through, {} unassigned, {} removed) frac_assigned={:.4} frac_bases_unextracted={:.4} mean_read_len_by_input={}",
            self.pool,
            self.total.separate_with_commas(),
            assigned.separate_with_commas(),
            self.pass_through.separate_with_commas(),
            self.unassigned.separate_with_commas(),
            self.removed.separate_with_commas(),
            frac(assigned, self.total),
            frac(self.retained_unextracted_bases, self.retained_denominator_bases),
            self.mean_read_len_by_input(),
        );
        for target in &self.targets {
            let tally = self.tally(target);
            log::info!(
                "metrics pool={} sample={} sub_sample={}: reads={} frac_of_pool={:.4} frac_bases_unextracted={:.4}",
                self.pool,
                target.sample,
                target.sub_sample.as_deref().unwrap_or(""),
                tally.reads.separate_with_commas(),
                frac(tally.reads, self.total),
                frac(tally.unextracted_bases, tally.denominator_bases),
            );
        }
        let unassigned_rate = self.unassigned_rate();
        if unassigned_rate >= 0.20 {
            log::warn!(
                "metrics pool={}: {:.1}% of reads were unassigned (>= 20%)",
                self.pool,
                unassigned_rate * 100.0,
            );
        }
    }
}

/// `numerator / denominator` as a fraction, or `0.0` when the denominator is
/// zero. Every caller passes `numerator <= denominator` (a count over its own
/// total), so the result is in `[0, 1]`; the assert enforces that invariant in
/// debug/test builds rather than silently emitting a fraction > 1.
fn frac(numerator: u64, denominator: u64) -> f64 {
    debug_assert!(
        numerator <= denominator,
        "frac numerator {numerator} exceeds denominator {denominator}"
    );
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

/// Serialize an `f64` as a fixed 6-decimal string, the metrics TSV's fraction
/// format.
fn six_dp<S: serde::Serializer>(
    value: &f64,
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    serializer.serialize_str(&format!("{value:.6}"))
}

/// One per-sample metrics row; the field names are the TSV header, in column
/// order.
#[derive(Serialize)]
struct PerSampleRow<'a> {
    pool: &'a str,
    sample: &'a str,
    sub_sample: &'a str,
    reads: u64,
    #[serde(serialize_with = "six_dp")]
    frac_of_pool: f64,
    #[serde(serialize_with = "six_dp")]
    frac_bases_unextracted: f64,
}

/// The pool-level summary metrics row; the field names are the TSV header, in
/// column order.
#[derive(Serialize)]
struct SummaryRow<'a> {
    pool: &'a str,
    total_reads: u64,
    reads_assigned: u64,
    reads_pass_through: u64,
    reads_unassigned: u64,
    reads_removed: u64,
    #[serde(serialize_with = "six_dp")]
    frac_assigned: f64,
    #[serde(serialize_with = "six_dp")]
    frac_bases_unextracted: f64,
    mean_read_len_by_input: String,
}

/// A tab-delimited, `\n`-terminated `csv::Writer` at `path`. The header row is
/// the serialized struct's field names (written before the first record).
/// Missing parent directories are created.
fn tsv_writer(path: &Path) -> Result<csv::Writer<BufWriter<File>>> {
    Ok(WriterBuilder::new()
        .delimiter(b'\t')
        .terminator(Terminator::Any(b'\n'))
        .from_writer(create(path)?))
}

/// Create a buffered TSV writer at `path`, creating missing parent directories.
fn create(path: &Path) -> Result<BufWriter<File>> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).with_context(|| {
                format!("failed to create metrics directory: {}", parent.display())
            })?;
        }
    }
    let file = File::create(path)
        .with_context(|| format!("failed to create metrics file: {}", path.display()))?;
    Ok(BufWriter::new(file))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn target(sample: &str, sub: Option<&str>) -> Target {
        Target {
            sample: sample.to_string(),
            sub_sample: sub.map(str::to_string),
        }
    }

    /// Read a written TSV back as lines.
    fn lines(path: &Path) -> Vec<String> {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .map(str::to_string)
            .collect()
    }

    #[test]
    fn test_per_sample_tsv_rows_and_fracs() {
        let dir = tempfile::tempdir().unwrap();
        let targets = vec![target("dna01", Some("lib1")), target("dna02", None)];
        let mut metrics = Metrics::new("pool1", &targets);
        for _ in 0..10 {
            metrics.record_processed();
        }
        // dna01: 6 reads, 9/100 bases unextracted each. dna02: 2 reads, fully
        // extracted.
        for _ in 0..6 {
            metrics.record_assigned(&targets[0], 100, 9);
        }
        for _ in 0..2 {
            metrics.record_assigned(&targets[1], 50, 0);
        }
        metrics.record_unassigned();
        metrics.record_unassigned();

        let path = dir.path().join("per_sample.tsv");
        metrics.write_per_sample_tsv(&path).unwrap();
        let rows = lines(&path);
        assert_eq!(
            rows[0],
            "pool\tsample\tsub_sample\treads\tfrac_of_pool\tfrac_bases_unextracted"
        );
        // dna01: 6 reads, frac_of_pool 0.6, frac_bases_unextracted 54/600 =
        // 0.09.
        assert_eq!(rows[1], "pool1\tdna01\tlib1\t6\t0.600000\t0.090000");
        // dna02: 2 reads, 0.2, fully extracted 0.0.
        assert_eq!(rows[2], "pool1\tdna02\t\t2\t0.200000\t0.000000");
    }

    #[test]
    fn test_summary_tsv_totals() {
        let dir = tempfile::tempdir().unwrap();
        let targets = vec![target("s1", None)];
        let mut metrics = Metrics::new("p", &targets);
        // Two inputs; lengths vary so the mean rounds: input 0 ->
        // mean(10,10,11,11,12), input 1 fixed 20.
        for lens in [[10usize, 20], [10, 20], [11, 20], [11, 20], [12, 20]] {
            metrics.record_processed();
            metrics.record_read_lengths(lens.into_iter());
        }
        metrics.record_assigned(&targets[0], 100, 10);
        metrics.record_assigned(&targets[0], 100, 10);
        metrics.record_unassigned();
        metrics.record_removed();
        metrics.record_removed();

        let path = dir.path().join("summary.tsv");
        metrics.write_summary_tsv(&path).unwrap();
        let rows = lines(&path);
        assert!(rows[0].ends_with("\tmean_read_len_by_input"));
        // total 5, assigned 2, pass_through 0, unassigned 1, removed 2,
        // frac_assigned 0.4, frac_bases_unextracted 20/200 = 0.1; input 0 mean
        // (10+10+11+11+12)/5 = 10.8 -> 11, input 1 = 20.
        assert_eq!(rows[1], "p\t5\t2\t0\t1\t2\t0.400000\t0.100000\t11,20");
    }

    #[test]
    fn test_unassigned_rate_threshold() {
        // The stderr warning fires at >= 20% unassigned; the rate predicate
        // drives it.
        let mut low = Metrics::new("p", &[]);
        for _ in 0..10 {
            low.record_processed();
        }
        low.record_unassigned();
        assert!(low.unassigned_rate() < 0.20, "1/10 is below the threshold");

        let mut high = Metrics::new("p", &[]);
        for _ in 0..10 {
            high.record_processed();
        }
        for _ in 0..3 {
            high.record_unassigned();
        }
        assert!(high.unassigned_rate() >= 0.20, "3/10 meets the threshold");
    }

    #[test]
    fn test_pass_through_is_fully_extracted() {
        let dir = tempfile::tempdir().unwrap();
        let mut metrics = Metrics::new("p", &[]);
        for _ in 0..3 {
            metrics.record_processed();
            metrics.record_pass_through(80, 0);
        }
        let path = dir.path().join("summary.tsv");
        metrics.write_summary_tsv(&path).unwrap();
        let rows = lines(&path);
        // all pass-through, nothing assigned, frac_bases_unextracted 0; no read
        // lengths recorded here, so the mean-read-len column is empty.
        assert_eq!(rows[1], "p\t3\t0\t3\t0\t0\t0.000000\t0.000000\t");
    }
}
