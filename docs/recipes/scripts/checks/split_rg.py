"""Split an aligned BAM by read group: every read group's records against samtools split."""

import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def body(ctx, w, runs, templates):
    same = runs["compare"]
    cb.check(same.out.split() == ["grp1", "identical", "grp2", "identical", "unknown", "identical"],
             f"each read group holds the records samtools split writes ({same.out.split()})")
    total = int(runs["total"].out)

    def placed(files):
        calls = collections.defaultdict(list)
        for rg, path in files:
            for name, rows in cb.by_name(ctx, w, path).items():
                calls[name] += [(rg, r) for r in rows]
        return {name: tuple(sorted(rows)) for name, rows in calls.items()}

    groups = ("grp1", "grp2", "unknown")
    row = cb.concord("samtools split", placed((g, f"by_rg/{g}.bam") for g in groups),
                     placed((g, f"{g}.samtools.bam") for g in groups), "headers aside")
    cb.check(sum(len(v) for v in placed((g, f"by_rg/{g}.bam") for g in groups).values()) == total, "every record is in one output")
    values = {
        "total": cb.fmt(total),
        "links": cb.source_links(ctx, ("samtools", "split.sam")),
        "split": cb.fence(runs["split"].command, "bash"),
        "counts_cmd": cb.fence(runs["counts"].command, "bash"),
        "counts": cb.fence(runs["counts"].out, "text"),
        "header_cmd": cb.fence(runs["header"].command, "bash"),
        "header": cb.fence(runs["header"].out, "text"),
        "compare_cmd": cb.fence(same.command, "bash"),
        "compare": cb.fence(same.out, "text"),
        "hdiff_cmd": cb.fence(runs["hdiff"].command, "bash"),
        "hdiff": cb.fence(runs["hdiff"].out, "diff"),
    }
    return values, [row]


cb.main("split_rg", body)
