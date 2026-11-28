"""Read the metrics files: the three files of the PhiX-bin run add up to the pool and match the files written."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def body(ctx, w, runs, templates):
    bins = ctx.work("bins")
    s = cb.read_tsv(w / "metrics.summary.tsv")[0]
    rows = cb.read_tsv(w / "metrics.per_sample.tsv")
    reasons = cb.read_tsv(w / "metrics.unassigned.tsv")
    assigned = sum(int(r["reads"]) for r in rows)
    un = sum(int(r["reads"]) for r in reasons)
    total = int(s["total_reads"])
    cb.check(assigned == int(s["reads_assigned"]) and un == int(s["reads_unassigned"]), "per-sample and reason rows add up")
    cb.check(assigned + un + int(s["reads_removed"]) + int(s["reads_pass_through"]) == total, "the four outcomes add up")
    pool = sum(1 for _ in cb.read_fastq(bins / "R1.fastq.gz"))
    cb.check(total == pool, "total_reads is the number of pairs in the pool")
    for row in rows:
        written = sum(1 for _ in cb.read_fastq(bins / "fastq" / f"{row['sample']}.R1.fastq.gz"))
        cb.check(int(row["reads"]) == written, f"{row['sample']}: per-sample reads equal the pairs in its files")
    removed = sum(1 for _ in cb.read_fastq(bins / "removed" / "PhiX.0.fastq.gz"))
    undetermined = sum(1 for _ in cb.read_fastq(bins / "fastq" / "Undetermined.0.fastq.gz"))
    cb.check(int(s["reads_removed"]) == removed and int(s["reads_unassigned"]) == undetermined, "bins match the summary")
    values = {
        "summary_cmd": cb.fence(runs["summary"].command, "bash"),
        "summary": cb.fence(runs["summary"].out, "tsv"),
        "per_cmd": cb.fence(runs["per"].command, "bash"),
        "per": cb.fence(runs["per"].out, "tsv"),
        "unassigned_cmd": cb.fence(runs["unassigned"].command, "bash"),
        "unassigned": cb.fence(runs["unassigned"].out, "tsv"),
        "assigned": cb.fmt(assigned),
        "un": cb.fmt(un),
        "removed": cb.fmt(int(s["reads_removed"])),
        "total": cb.fmt(total),
    }
    return values, []


cb.main("metrics", body)
