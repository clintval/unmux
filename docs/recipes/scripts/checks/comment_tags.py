"""Lift SAM tags from FASTQ comments: every record of the round trip against the original BAM's, tags included."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def body(ctx, w, runs, templates):
    original = cb.sam_records(ctx, w, "modbam.unmapped.bam")
    back = cb.sam_records(ctx, w, "modbam.roundtrip.bam")
    cb.check(len(original) == len(back), "the round trip keeps every record")
    for (fa, ta), (fb, tb) in zip(original, back):
        cb.check(fa[0] == fb[0] and fa[9] == fb[9] and fa[10] == fb[10],
                 f"{fa[0]}: the name, bases and qualities survive the round trip")
        cb.check({k: v for k, v in ta.items() if k != "RG"} == {k: v for k, v in tb.items() if k != "RG"},
                 f"{fa[0]}: every tag but RG survives the round trip")
        cb.check(tb["RG"] == "modbam" != ta["RG"], f"{fa[0]}: RG is the new pool's read group")
    cb.check(all("NM" not in t and "MD" not in t for _, t in cb.sam_records(ctx, w, "direct.bam")),
             "read from BAM directly, the records drop NM and MD")
    values = {
        "links": cb.source_links(ctx, ("modkit", "bc_anchored_10_reads.unmapped.bam")),
        "prep": cb.fence(runs["prep"].command, "bash"),
        "header": cb.fence(runs["prep"].out, "text"),
        "convert": cb.fence(runs["convert"].command, "bash"),
        "view_cmd": cb.fence(runs["view"].command, "bash"),
        "view": cb.fence(runs["view"].out, "text"),
        "reads": cb.fmt(len(back)),
    }
    return values, []


cb.main("comment_tags", body)
