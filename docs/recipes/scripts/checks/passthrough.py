"""Plain SAM or BAM pass-through: the records against samtools reset, with MM, ML and the reverse strand checked."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def body(ctx, w, runs, templates):
    same = runs["compare"]
    cb.check(same.out.strip() == "identical", f"the records equal what samtools reset writes ({same.out.strip()!r})")
    aligned = cb.sam_records(ctx, w, "modbam.aligned.bam")
    unmapped = cb.sam_records(ctx, w, "modbam.unmapped.bam")
    cb.check(all(a[1]["MM"] == u[1]["MM"] and a[1]["ML"] == u[1]["ML"] for a, u in zip(aligned, unmapped)),
             "MM and ML are unchanged on both records")
    cb.check(all(u[0][9] == (cb.revcomp(a[0][9]) if int(a[0][1]) & 16 else a[0][9]) for a, u in zip(aligned, unmapped)),
             "the reverse-strand read is the reverse complement of its aligned SEQ")
    row = cb.concord("samtools reset", cb.by_name(ctx, w, "modbam.unmapped.bam", cb.NO_MAPQ),
                     cb.by_name(ctx, w, "reset.bam", cb.NO_MAPQ), "MAPQ and the tags left out")
    dropped = sorted(set(aligned[0][1]) - set(unmapped[0][1]))
    values = {
        "links": cb.source_links(ctx, ("modkit", "fwd_rev_modbase_records.bam")),
        "convert": cb.fence(runs["convert"].command, "bash"),
        "sq": str(int(runs["sq"].out)),
        "before_cmd": cb.fence(runs["before"].command, "bash"),
        "before": cb.fence(runs["before"].out, "text"),
        "after_cmd": cb.fence(runs["after"].command, "bash"),
        "after": cb.fence(runs["after"].out, "text"),
        "dropped": ", ".join(dropped),
        "compare_cmd": cb.fence(same.command, "bash"),
        "compare": cb.fence(same.out, "text"),
    }
    return values, [row]


cb.main("passthrough", body)
