"""Nanopore native barcodes at both ends: every read's call against its run, the trims, and a trim's per-base tags."""

import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb
import nanopore


def body(ctx, w, runs, templates):
    origin, records = nanopore.pooled(ctx)
    cb.check(len(origin) == len(records), "read names are unique across the 24 runs")
    run, mods = runs["demux"], runs["mods"]
    reads = {name: seq for name, seq, _ in records}
    correct = wrong = 0
    calls = {}
    cut = []
    for path in sorted((w / "demux").glob("NB*.fastq.gz")):
        sample = path.name.split(".")[0]
        for name, seq, _ in cb.read_fastq(path):
            correct += origin[name] == sample
            wrong += origin[name] != sample
            calls[name] = sample
            start = reads[name].find(seq)
            cb.check(start >= 0, f"{name}: the trimmed read is part of the read")
            cut.append((start, len(reads[name]) - start - len(seq)))
    unclassified = sum(1 for _ in cb.read_fastq(w / "demux" / "unclassified.fastq.gz"))
    cb.check(wrong == 0, f"no read goes to another barcode ({wrong})")
    cb.check(correct + unclassified == len(records), "every read is assigned or unclassified")
    cb.check(all(a <= 150 and b <= 150 for a, b in cut), "the cuts fall inside the search windows")
    front_cut = sorted(a for a, _ in cut)
    rear_cut = sorted(b for _, b in cut)
    tight = collections.Counter()
    for path in sorted((w / "demux8").glob("NB*.fastq.gz")):
        for name, _, _ in cb.read_fastq(path):
            tight["right" if origin[name] == path.name.split(".")[0] else "wrong"] += 1
    tight["unclassified"] = sum(1 for _ in cb.read_fastq(w / "demux8" / "unclassified.fastq.gz"))
    cb.check(sum(tight.values()) == len(records), "with 8 edits every read is assigned or unclassified")
    cb.check(tight["wrong"] == 0 and tight["right"] < correct, "8 edits assign fewer reads and none to a wrong barcode")
    kept_tags, trimmed_tags = (line.split() for line in mods.out.strip().split("\n"))
    per_base = [tag for tag in kept_tags if tag not in trimmed_tags]
    cb.check(set(per_base) == {"MM", "ML", "mv"} and trimmed_tags == [tag for tag in kept_tags if tag not in per_base],
             "a trimming template keeps the read's tags except MM, ML and mv")
    left_out = cb.log_lines(mods.err, ["left out per-base tags"])
    cb.check(bool(left_out) and all(tag in left_out for tag in per_base), "unmux warns which per-base tags it left out")
    warn, listed = cb.pair_warnings(run.err)
    cb.check(listed > 0 and "window" not in warn, "unmux warns about tag pairs and not about the windows")
    row = cb.concord("the run each read came from", {n: calls.get(n) for n in origin}, origin,
                     f"reads left unclassified ({cb.fmt(unclassified)})")
    values = {
        "project_link": cb.link("PRJNA730334", ctx.sources[("sra-verd6", "SRR14560737.fastq")].url),
        "reads": cb.fmt(len(records)),
        "fetch": cb.fence(cb.fetch_command(ctx, "sra-verd6", "SRR14560737.fastq"), "bash"),
        "kit_links": cb.source_links(ctx, ("dorado", "barcode_kits.cpp")),
        "tags_cmd": cb.fence(runs["tags"].command, "bash"),
        "tags": cb.fence(runs["tags"].out, "tsv"),
        "demux": cb.fence(run.command, "bash"),
        "counts_cmd": cb.fence(runs["counts"].command, "bash"),
        "counts": cb.fence(runs["counts"].out, "text"),
        "listed": str(listed),
        "warn": cb.fence(warn, "text"),
        "correct": cb.fmt(correct),
        "share": cb.pct(correct, len(records)),
        "wrong": cb.fmt(wrong) if wrong else "none",
        "unclassified": cb.fmt(unclassified),
        "front_cut": str(front_cut[len(front_cut) // 2]),
        "rear_cut": str(rear_cut[len(rear_cut) // 2]),
        "tight_right": cb.fmt(tight["right"]),
        "tight_wrong": cb.fmt(tight["wrong"]) if tight["wrong"] else "none",
        "tight_unclassified": cb.fmt(tight["unclassified"]),
        "mods_cmd": cb.fence(mods.command, "bash"),
        "mods": cb.fence(mods.out, "text"),
        "left_out": cb.fence(left_out, "text"),
    }
    return values, [row]


cb.main("nanopore", body)
