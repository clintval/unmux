"""Paired FASTQ to an unmapped BAM: every record against its FASTQ read, and the records against samtools import."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def body(ctx, w, runs, templates):
    run, same = runs["convert"], runs["compare"]
    mapq = runs["mapq"].out.split()
    records = cb.sam_records(ctx, w, "SRR390728.unmapped.bam")
    mates = [list(cb.read_fastq(w / f"SRR390728_{r}.fastq.gz")) for r in (1, 2)]
    cb.check(len(records) == 2 * len(mates[0]), "one record per read")
    for k, (fields, tags) in enumerate(records):
        name, seq, qual = mates[k % 2][k // 2]
        cb.check((fields[0], fields[9], fields[10]) == (name, seq, qual), f"record {k} equals its FASTQ read")
        cb.check(int(fields[1]) == (77 if k % 2 == 0 else 141) and tags == {"RG": "SRR390728"}, f"record {k} flags and tags")
    cb.check(same.out.strip() == "identical" and mapq == ["0"], f"samtools import comparison: {same.out!r} {mapq}")
    row = cb.concord("samtools import", cb.by_name(ctx, w, "SRR390728.unmapped.bam", cb.NO_MAPQ),
                     cb.by_name(ctx, w, "samtools.unmapped.bam", cb.NO_MAPQ), "MAPQ (255 against 0) and the tags left out")
    values = {
        "pairs": cb.fmt(len(mates[0])),
        "run_link": cb.link("SRR390728", ctx.sources[("sra-SRR390728", "SRR390728_1.fastq")].url),
        "fetch": cb.fence(cb.fetch_command(ctx, "sra-SRR390728", "SRR390728_1.fastq"), "bash"),
        "convert": cb.fence(run.command, "bash"),
        "log": cb.fence(cb.log_lines(run.err, ["input", "output", "Demultiplexed", "metrics"]), "text"),
        "view_cmd": cb.fence(runs["view"].command, "bash"),
        "view": cb.fence(runs["view"].out, "text"),
        "records": cb.fmt(len(records)),
        "compare_cmd": cb.fence(same.command, "bash"),
        "compare": cb.fence(same.out, "text"),
    }
    return values, [row]


cb.main("fastq_to_ubam", body)
