"""Dual-index demultiplexing with a sample sheet: every call against UDiTaS's demultiplex() rules and against fqtk."""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb
import illumina


def fqtk_comparison(w, runs, mine):
    version = cb.tool_version(runs["fqtk_version"], "fqtk")
    close = cb.concord(version, mine, illumina.fqtk_calls(w / "fqtk"), "with `--max-mismatches 2 --min-mismatch-delta 1`")
    loose = illumina.fqtk_calls(w / "fqtk-defaults")
    defaults = cb.concord(version, mine, loose, "")
    i7 = {r["id"]: r["seq"] for r in cb.read_tsv(w / "i7.tsv")}
    i5 = {r["id"]: r["seq"] for r in cb.read_tsv(w / "i5.tsv")}
    sheet_rows = {r["sample"]: (i7[r["i7"]], i5[r["i5"]]) for r in cb.read_tsv(w / "samples.tsv")}
    seen = {n.split()[0]: s for n, s, _ in cb.read_fastq(w / "Undetermined_S0_L001_I1_001.fastq.gz")}
    seen5 = {n.split()[0]: s[:8] for n, s, _ in cb.read_fastq(w / "Undetermined_S0_L001_I2_001.fastq.gz")}

    def mismatches(read, sample):
        return cb.hamming(seen[read], sheet_rows[sample][0]), cb.hamming(seen5[read], sheet_rows[sample][1])

    lost = [n for n in mine if mine[n] and not loose.get(n)]
    cb.check(close["different"] == close["unmux"] == close["other"] == 0, "with two mismatches and a lead of one, fqtk makes unmux's calls")
    cb.check(defaults["different"] == defaults["other"] == 0 and len(lost) == defaults["unmux"], "fqtk's defaults only lose reads")
    cb.check(all(mismatches(n, mine[n]) == (1, 1) for n in lost), "every read fqtk's defaults lose has one mismatch in each index")
    twice = sum(1 for n in mine if mine[n] and 2 in mismatches(n, mine[n]))
    values = {
        "fqtk_version": version,
        "fqtk_sheet": cb.fence(runs["fqtk_sheet"].command + "\n" + runs["fqtk_head"].command, "bash"),
        "fqtk_head": cb.fence(runs["fqtk_head"].out, "tsv"),
        "fqtk": cb.fence(runs["fqtk"].command, "bash"),
        "fqtk_table": cb.agreement_table("fqtk", [("these settings", close), ("fqtk's defaults", defaults)]),
        "fqtk_reads": cb.fmt(close["reads"]),
        "twice": f"; {cb.fmt(twice)} reads here have that." if twice else "; no read here has that.",
        "lost": cb.fmt(len(lost)),
    }
    return values, close


def body(ctx, w, runs, templates):
    with open(w / "sample_info.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    i1 = list(cb.read_fastq(w / "Undetermined_S0_L001_I1_001.fastq.gz"))
    i2 = list(cb.read_fastq(w / "Undetermined_S0_L001_I2_001.fastq.gz"))
    by_name = {name.split()[0]: (s, q) for name, s, q in i2}
    i1 = [(n.split()[0], s, q) for n, s, q in i1]
    i2 = [(n.split()[0], s, q) for n, s, q in i2]
    calls = {}
    for path in sorted((w / "bam").glob("*.bam")):
        sample = path.name.split(".")[0]
        for fields, tags in cb.sam_records(ctx, w, f"bam/{path.name}"):
            seq, qual = by_name[fields[0]]
            cb.check(tags["RX"] == seq[8:] and tags["QX"] == qual[8:], f"{fields[0]}: RX and QX equal I2 bases 8 to 18")
            if int(fields[1]) & 64:
                cb.check(fields[0] not in calls, f"{fields[0]} is in one sample only")
                calls[fields[0]] = sample
    for name, _, _ in i1:
        calls.setdefault(name, "unassigned")
    plain = illumina.uditas_calls(i1, i2, rows, masked=False)
    masked = illumina.uditas_calls(i1, i2, rows)

    def same(a, b):
        return a == b or (a == "unassigned" and b in ("mismatched_adapters", "not_in_exp_list"))

    cb.check(all(same(calls[n], plain[n]) for n in calls), "unmux equals the UDiTaS rules without masking")
    agree = sum(1 for n in calls if same(calls[n], masked[n]))
    differ = [n for n in calls if not same(calls[n], masked[n])]
    cb.check(all(masked[n] == "mismatched_adapters" and calls[n] != "unassigned" for n in differ),
             "every difference is a read UDiTaS leaves unmatched and unmux assigns")
    quals = {n: (q1, q2[:8]) for (n, _, q1), (_, _, q2) in zip(i1, i2)}
    low = [sum(1 for n in differ if any(ord(c) - 33 < 12 for c in quals[n][k])) for k in (0, 1)]
    assigned = sum(1 for c in calls.values() if c != "unassigned")
    mine = {n: None if c == "unassigned" else c for n, c in calls.items()}
    rules = cb.concord("UDiTaS's demultiplex() rules", mine,
                       {n: None if c in ("mismatched_adapters", "not_in_exp_list") else c for n, c in masked.items()},
                       f"quality masking ({cb.fmt(len(differ))} reads)" if differ else "")
    fqtk_values, fqtk_row = fqtk_comparison(w, runs, mine)
    masking = (
        f"unmux makes the same call on all {cb.fmt(agree)} reads. UDiTaS also masks bases below "
        "quality 12 to N before matching, which unmux does not do; no read here has such a base "
        "that pushes an index past one mismatch, and on a run that has them, UDiTaS would leave "
        "those reads unmatched where unmux assigns them."
        if not differ else
        f"unmux makes the same call on {cb.fmt(agree)} reads. UDiTaS also masks bases below quality "
        f"12 to N before matching, which unmux does not do, and the other {cb.fmt(len(differ))} reads "
        "are all ones that masking pushes over one mismatch, so UDiTaS leaves them unmatched while "
        f"unmux assigns them ({cb.fmt(low[0])} have a masked base in I1, {cb.fmt(low[1])} in the i5 part "
        "of I2). With the masking left out of the reimplementation, every read agrees."
    )
    files = [("uditas", f"Undetermined_S0_L001_{r}_001.fastq") for r in ("R1", "R2", "I1", "I2")]
    values = {
        "reads": cb.fmt(len(i1)),
        "links": cb.source_links(ctx, *files, ("uditas", "sample_info.csv")),
        "samples": str(len(rows)),
        "index_cmd": cb.fence(runs["index"].command, "bash"),
        "index": cb.fence(runs["index"].out, "text"),
        "sheet": cb.fence(runs["sheet"].command, "bash"),
        "sheets_cmd": cb.fence(runs["sheets"].command, "bash"),
        "sheets": cb.fence(runs["sheets"].out, "tsv"),
        "demux": cb.fence(runs["demux"].command, "bash"),
        "counts_cmd": cb.fence(runs["counts"].command, "bash"),
        "counts": cb.fence(runs["counts"].out, "text"),
        "view_cmd": cb.fence(runs["view"].command, "bash"),
        "view": cb.fence(runs["view"].out, "text"),
        "assigned": cb.fmt(assigned),
        "helpers": cb.link("uditas_helpers.py", ctx.sources[("uditas", "uditas_helpers.py")].url),
        "masking": masking,
        **fqtk_values,
    }
    return values, [rules, fqtk_row]


cb.main("uditas", body)
