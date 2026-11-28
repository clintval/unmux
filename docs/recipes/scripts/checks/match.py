"""Keep the reads of called cells with a joined match: the kept reads and their names against a lookup of the windows."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb
import splitseq


def body(ctx, w, runs, templates):
    run = runs["match"]
    cells = splitseq.called_cells(ctx, w)
    oracle = splitseq.splitseq_oracle(w)
    kept = cb.sam_records(ctx, w, "cells.unmapped.bam")
    other = {f[0] for f, _ in cb.sam_records(ctx, w, "other.unmapped.bam")}
    want = {n: cells[c[0].replace("-", "")] for n, c in oracle.items() if c and c[0].replace("-", "") in cells}
    cb.check({f[0]: t.get("CB") for f, t in kept} == want, "kept reads are the oracle's called-cell reads, each with its cell")
    cb.check(len(kept) + len(other) == len(oracle) and not other & set(want), "every other read goes to other.unmapped.bam")
    r2 = {n: s for n, s, _ in cb.read_fastq(w / "R2.fastq.gz")}
    observed = {n: "".join(s[a:b] for a, b in ((78, 86), (48, 56), (10, 18))) for n, s in r2.items()}
    corrected = sum(1 for f, _ in kept if oracle[f[0]][0].replace("-", "") != observed[f[0]])
    names = len(set(want.values()))
    warn = cb.log_lines(run.err, ["WARN"])
    share = cb.pct(len(other), len(oracle))
    cb.check("no `--sample`" in warn and "`CY` is left out" in warn and f"{share} of reads were unassigned" in warn,
             "unmux warns about the missing sample target, the qualities of CB and the unassigned share")
    values = {
        "called": cb.fmt(len(set(cells.values()))),
        "geo_link": cb.link("GSM5169184", ctx.sources[("geo", splitseq.GEO_CELLS)].url),
        "cells_cmd": cb.fence(runs["cells"].command, "bash"),
        "cells": cb.fence(runs["cells"].out, "tsv"),
        "match": cb.fence(run.command, "bash"),
        "view_cmd": cb.fence(runs["view"].command, "bash"),
        "view": cb.fence(runs["view"].out, "text"),
        "kept": cb.fmt(len(kept)),
        "names": cb.fmt(names),
        "other": cb.fmt(len(other)),
        "corrected": cb.fmt(corrected),
        "share": share,
        "warnings": cb.fence(warn, "text"),
    }
    return values, []


cb.main("match", body)
