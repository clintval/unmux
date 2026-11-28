"""Unmapped BAM back to paired FASTQ: the FASTQs against samtools fastq and against the reads of the first recipe."""

import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def body(ctx, w, runs, templates):
    same = runs["compare"]
    cb.check(same.out.split("\n")[:2] == ["R1 identical", "R2 identical"],
             f"R1 and R2 equal what samtools fastq writes ({same.out.strip()!r})")
    back = [(s, q) for _, s, q in cb.read_fastq(w / "SRR390728.R1.fastq.gz")]
    orig = [(s, q) for _, s, q in cb.read_fastq(ctx.work("fastq_to_ubam") / "SRR390728_1.fastq.gz")]
    cb.check(back == orig, "the FASTQ written back holds the bases and qualities of the input R1")

    def pairs(prefix):
        reads = collections.defaultdict(list)
        for r in (1, 2):
            for name, seq, qual in cb.read_fastq(w / f"{prefix}.R{r}.fastq.gz"):
                reads[name.split()[0]].append((name.replace("\t", " "), seq, qual))
        return {name: tuple(mates) for name, mates in reads.items()}

    row = cb.concord("samtools fastq -T RG", pairs("SRR390728"), pairs("samtools"), "samtools' tab before the tags read as a space")
    values = {
        "convert": cb.fence(runs["convert"].command, "bash"),
        "listing": cb.fence(runs["listing"].out, "text"),
        "view_cmd": cb.fence(runs["view"].command, "bash"),
        "view": cb.fence(runs["view"].out, "text"),
        "compare_cmd": cb.fence(same.command, "bash"),
        "compare": cb.fence(same.out, "text"),
        "pairs": cb.fmt(len(back)),
    }
    return values, [row]


cb.main("ubam_to_fastq", body)
