"""Inline barcodes into unmapped BAMs: every read of 32 pooled RAD-seq runs against its run, and against fqtk."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb
import illumina


def body(ctx, w, runs, templates):
    origin, records = illumina.inline_pool(ctx)
    reads = {name: (seq, qual) for name, seq, qual in records}
    barcodes = {sample: barcode for _, sample, barcode in illumina.CARABUS}
    correct = wrong = 0
    calls = {}
    mismatched = None
    for path in sorted((w / "bam").glob("*.unmapped.bam")):
        sample = path.name.split(".")[0]
        for fields, tags in cb.sam_records(ctx, w, f"bam/{path.name}"):
            seq, qual = reads[fields[0]]
            correct += origin[fields[0]] == sample
            wrong += origin[fields[0]] != sample
            calls[fields[0]] = sample
            cb.check((tags["BC"], tags["QT"]) == (seq[:5], qual[:5]), f"{fields[0]}: BC and QT are the read's first 5 bases")
            cb.check((fields[9], fields[10]) == (seq[5:], qual[5:]), f"{fields[0]}: the record is the read from base 5 on")
            if mismatched is None and tags["BC"] != barcodes[sample]:
                mismatched = (fields[0], sample)
    cb.check(wrong == 0 and correct == len(records), f"every read goes to its run's sample ({correct}, {wrong})")
    view = ctx.sh(w, templates["view"].replace("{{sample}}", mismatched[1]).replace("{{read}}", mismatched[0]))
    exact_n = sum(len(cb.sam_records(ctx, w, f"exact/{p.name}")) for p in (w / "exact").glob("*.bam"))
    site = sum(1 for s, _ in reads.values() if s[5:10] == "TGCAG")
    mine = {n: calls.get(n) for n in origin}
    truth_row = cb.concord("the run each read came from", mine, origin, "")
    version = cb.tool_version(runs["fqtk_version"], "fqtk")
    fqtk_row = cb.concord(version, mine, illumina.fqtk_calls(w / "fqtk"), "with `--max-mismatches 1 --min-mismatch-delta 1`")
    cb.check(site == len(records), "every read has the remnant, so the site group changes no call")
    values = {
        "study_link": cb.link("DRP001067", ctx.sources[("sra-carabus", "DRR021959.fastq")].url),
        "fetch": cb.fence(cb.fetch_command(ctx, "sra-carabus", "DRR021959.fastq"), "bash"),
        "head_cmd": cb.fence(runs["head"].command, "bash"),
        "head": cb.fence(runs["head"].out, "text"),
        "demux": cb.fence(runs["demux"].command, "bash"),
        "view_cmd": cb.fence(view.command, "bash"),
        "view": cb.fence(view.out, "text"),
        "reads": cb.fmt(len(records)),
        "site": cb.fmt(site),
        "exact": cb.fmt(exact_n),
        "inexact": cb.fmt(len(records) - exact_n),
        "fqtk_version": version,
        "fqtk": cb.fence(runs["fqtk"].command, "bash"),
        "fqtk_table": cb.agreement_table("fqtk", [("fqtk", fqtk_row)]),
        "fqtk_reads": ("all " if fqtk_row["same"] == fqtk_row["reads"] else f"{cb.fmt(fqtk_row['same'])} of ")
                      + cb.fmt(fqtk_row["reads"]),
    }
    return values, [truth_row, fqtk_row]


cb.main("inline", body)
