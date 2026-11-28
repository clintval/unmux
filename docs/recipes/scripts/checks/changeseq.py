"""Join read pairs end to end for CHANGE-seq: every merged read against the pipeline's mergeReads()."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def body(ctx, w, runs, templates):
    counts, mine, theirs = {}, {}, {}
    for sample, merged in (("TEST", "TEST_merged"), ("TEST_control", "control_TEST_merged")):
        r1 = list(cb.read_fastq(w / f"{sample}.r1.fastq.gz"))
        r2 = list(cb.read_fastq(w / f"{sample}.r2.fastq.gz"))
        out = list(cb.read_fastq(w / f"{merged}.fastq.gz"))
        cb.check(len(out) == len(r1) == len(r2), f"{merged}: one read per pair")
        for (n1, s1, q1), (_, s2, q2), (name, seq, qual) in zip(r1, r2, out):
            cb.check(name == n1.split()[0], f"{merged}: the name of {n1}")
            cb.check(seq == cb.revcomp(s1) + s2 and qual == q1[::-1] + q2, f"{merged}: {name} equals mergeReads()")
            mine[(sample, name)] = (seq, qual)
            theirs[(sample, name)] = (cb.revcomp(s1) + s2, q1[::-1] + q2)
        counts[sample] = len(out)
    files = [("changeseq", f"{s}.{r}.fastq.gz") for s in ("TEST", "TEST_control") for r in ("r1", "r2")]
    values = {
        "treated": cb.fmt(counts["TEST"]),
        "control": cb.fmt(counts["TEST_control"]),
        "links": cb.source_links(ctx, *files),
        "inputs_cmd": cb.fence(runs["inputs"].command, "bash"),
        "inputs": cb.fence(runs["inputs"].out, "text"),
        "merge": cb.fence(runs["merge"].command, "bash"),
        "view_cmd": cb.fence(runs["view"].command, "bash"),
        "view": cb.fence(runs["view"].out, "text"),
        "merge_link": cb.link("mergeReads.py", ctx.sources[("changeseq", "mergeReads.py")].url),
    }
    return values, [cb.concord("mergeReads() of CHANGE-seq", mine, theirs, "")]


cb.main("changeseq", body)
