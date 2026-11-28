"""Direct RNA reads with U bases to BAM: T for U as samtools import stores it, and a T-only tag that finds the U read."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def body(ctx, w, runs, templates):
    seq_in = next(cb.read_fastq(w / "direct_rna.fastq"))[1]
    rec = cb.sam_records(ctx, w, "direct_rna.bam")[0]
    cb.check(rec[0][9] == seq_in.replace("U", "T"), "the BAM record holds the read with T for U")
    cb.check(cb.sam_records(ctx, w, "samtools.bam")[0][0][9] == rec[0][9], "BAM bases equal samtools import's")
    probe = cb.revcomp(seq_in[100:124].replace("U", "T"))
    ctx.sh(w, templates["probe"].replace("{{probe}}", probe))
    cb.check(len(list(cb.read_fastq(w / "probe.dna.fastq"))) == 1, "a T-only tag from the other strand finds the U read")
    values = {
        "links": cb.source_links(ctx, ("dorado", "fastq_with_us_and_tags.fq")),
        "head_cmd": cb.fence(runs["head"].command, "bash"),
        "head": cb.fence(runs["head"].out, "text"),
        "convert": cb.fence(runs["convert"].command, "bash"),
        "view_cmd": cb.fence(runs["view"].command, "bash"),
        "view": cb.fence(runs["view"].out, "text"),
        "u": cb.fmt(seq_in.count("U")),
        "text_cmd": cb.fence(runs["text"].command, "bash"),
        "text": cb.fence(runs["text"].out, "text"),
    }
    return values, []


cb.main("direct_rna", body)
