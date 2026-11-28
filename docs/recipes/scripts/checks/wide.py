"""Search windows wider than their tags: the README windows against windows 1 bp wider under each mode."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb
import splitseq


def body(ctx, w, runs, templates):
    cells = splitseq.called_cells(ctx, w)
    base = ctx.work("splitseq")
    reference = {f[0]: t.get("CB") for f, t in cb.sam_records(ctx, base, "splitseq.unmapped.bam")}
    rows = [
        ("README windows, `dist=1`", splitseq.splitseq_scores(ctx, base, "splitseq.unmapped.bam", cells)),
        ("1 bp wider, `dist=1:1:2,mode=all`", splitseq.splitseq_scores(ctx, w, "splitseq.unmapped.bam", cells, reference)),
        ("1 bp wider, `dist=1:1:2,mode=nearest`",
         splitseq.splitseq_scores(ctx, w / "nearest", "splitseq.unmapped.bam", cells, reference)),
    ]
    table = "\n".join([
        "| Command | Written | With a CB | CB of a called cell | Other CB | No CB | CB changed from the README windows |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ] + [f"| {label} | {cb.fmt(n)} | {cb.fmt(found)} | {cb.fmt(inc)} | {cb.fmt(found - inc)} | {cb.fmt(n - found)} | "
         + ("" if label.startswith("README") else cb.fmt(chg)) + " |" for label, (n, found, inc, chg) in rows])
    (_, readme), (_, wide_all), (_, wide_nearest) = rows
    cb.check(readme[0] == readme[1] and wide_all[0] == wide_all[1], "every written read has a CB except under mode=nearest")
    missing = wide_nearest[0] - wide_nearest[1]
    cb.check(missing > 0, "mode=nearest writes reads without a CB")
    counted = runs["count"]
    cb.check(counted.out.split() == [str(wide_nearest[0]), str(wide_nearest[1])], "samtools counts match")
    tied = [json.loads(t["ZS"]) for _, t in cb.sam_records(ctx, w / "nearest", "qc.unmapped.bam") if "CB" not in t]
    cb.check(len(tied) == missing and all(len(z["groups"]) < 3 for z in tied), "QC of a read without CB lists fewer than three rounds")
    run, run2 = runs["wide_all"], runs["wide_nearest"]
    warned = [line for line in run.err.split("\n") if "wider than" in line]
    pairs = len([line for line in run.err.split("\n") if "can both match" in line])
    quiet = "wider than" not in run2.err
    oracle = splitseq.splitseq_oracle(w)
    values = {
        "missing_round": cb.fmt(sum(1 for v in oracle.values() if not v)),
        "wide_all": cb.fence(run.command, "bash"),
        "warned": cb.fence("\n".join(warned), "text"),
        "pairs": cb.fmt(pairs),
        "wide_nearest": cb.details("The command", cb.fence(run2.command, "bash")),
        "quiet": ("It prints no window warning: under `mode=nearest`, a window wider than the tags by twice "
                  "the indel budget is as wide as unmux advises."
                  if quiet else cb.fence(cb.log_lines(run2.err, ["wider than"]), "text")),
        "table": table,
        "more_cells": cb.fmt(wide_all[2] - readme[2]),
        "changed_all": cb.fmt(wide_all[3]),
        "extra_other": cb.fmt((wide_all[1] - wide_all[2]) - (readme[1] - readme[2])),
        "changed_nearest": cb.fmt(wide_nearest[3]),
        "missing": cb.fmt(missing),
        "count_cmd": cb.fence(counted.command, "bash"),
        "counted": cb.fence(counted.out, "text"),
    }
    return values, []


cb.main("wide", body)
