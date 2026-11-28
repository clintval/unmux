"""PacBio barcode pairs: every read's sample against lima's, the trims, and why the two leave out the reads they do."""

import collections
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def body(ctx, w, runs, templates):
    run = runs["demux"]
    tag_ids = {row["seq"]: row["id"] for row in cb.read_tsv(w / "barcodes.tsv")}
    tag_len = {len(seq) for seq in tag_ids}
    cb.check(len(tag_len) == 1, "every tag has one length")
    tag_len = tag_len.pop()
    lima_version = runs["version"].out.split("\n")[0].strip()
    cb.check(re.fullmatch(r"lima \d+\.\d+\.\d+", lima_version) is not None, f"lima reports its version ({lima_version!r})")
    with open(w / "biosamples.csv", newline="") as fh:
        pairs = {}
        for row in csv.DictReader(fh):
            a, b = row["Barcode"].split("--")
            pairs[(a, b)] = pairs[(b, a)] = row["Bio Sample Name"]
    samples = len(set(pairs.values()))
    report = {row["ZMW"].split("/")[1]: row for row in cb.read_tsv(w / "lima.lima.report")}
    lima = {z: pairs.get((row["IdxFirstNamed"], row["IdxCombinedNamed"])) if row["PassedFilters"] == "1" else None
            for z, row in report.items()}
    inputs = cb.sam_records(ctx, w, "hifi_reads.bam", "--no-PG")
    reads = {f[0].split("/")[1]: f[9] for f, _ in inputs}
    carried = {f[0].split("/")[1]: {k: v for k, v in tags.items() if k != "RG"} for f, tags in inputs}
    cb.check(len(reads) == len(inputs), "one record per ZMW")
    unmux_calls, cuts = {}, []
    for path in sorted(p for p in (w / "unmux").glob("*.bam") if p.name != "unassigned.bam"):
        for fields, tags in cb.sam_records(ctx, w, f"unmux/{path.name}"):
            z = fields[0].split("/")[1]
            unmux_calls[z] = path.name.split(".")[0]
            start = reads[z].find(fields[9])
            cb.check(start >= 0, f"{fields[0]}: the trimmed read is part of the read")
            cuts.append((start, len(reads[z]) - start - len(fields[9])))
            cb.check({k: v for k, v in tags.items() if k != "RG"} == carried[z] and tags["RG"] == unmux_calls[z],
                     f"{fields[0]}: the record keeps its read's tags, with the sample's read group")
    cb.check(all(a <= 60 and b <= 60 for a, b in cuts), "the cuts fall inside the search windows")
    table = collections.Counter()
    for z in reads:
        a, b = lima.get(z), unmux_calls.get(z)
        table[("same" if a == b else "different") if a and b else "lima only" if a else "unmux only" if b else "neither"] += 1
    cb.check(table["different"] == 0, f"lima and unmux never call different samples ({table['different']})")
    lima_n = sum(1 for z in reads if lima.get(z))
    unmux_only = [z for z in reads if unmux_calls.get(z) and not lima.get(z)]
    lima_only = [z for z in reads if lima.get(z) and not unmux_calls.get(z)]
    neither = [z for z in reads if not lima.get(z) and not unmux_calls.get(z)]
    cb.check(all(z in report and report[z]["PassedFilters"] != "1" for z in unmux_only),
             "lima removes every read only unmux assigns")
    scores = sorted(int(report[z]["ScoreCombined"]) for z in unmux_only)
    cb.check(scores[-1] <= 80, "every read only unmux assigns has a lima combined score of 80 or less")
    qc = {f[0].split("/")[1]: json.loads(t["ZS"]) for f, t in cb.sam_records(ctx, w, "qc/unassigned.bam")}
    cb.check(set(qc) == set(lima_only) | set(neither), "unmux leaves out exactly the reads it does not assign")

    def ends(z):
        found = {g["g"] for g in qc[z]["groups"]}
        return "both" if {"start", "end"} <= found else "one" if found else "none"

    sides = collections.Counter(ends(z) for z in lima_only)
    left = collections.Counter(ends(z) for z in neither)
    lima_short = sum(1 for z in neither if z in report and int(report[z]["NumBarcodedRegionsPassed"]) < 2)
    twice = unlisted = 0
    for z in (z for z in neither if ends(z) == "both"):
        found = {g["g"]: tag_ids[g["tag"]] for g in qc[z]["groups"]}
        row = report[z]
        cb.check({found["start"], found["end"]} == {row["IdxFirstNamed"], row["IdxCombinedNamed"]},
                 f"{z}: lima reports the barcodes unmux found")
        cb.check((found["start"], found["end"]) not in pairs, f"{z}: the barcodes are not a sample's pair")
        twice += found["start"] == found["end"]
        unlisted += found["start"] != found["end"]
    missing = [z for z in neither if z not in report]
    cb.check(all(z in report for z in reads if z not in missing), "lima's report has every other read")
    warn, listed = cb.pair_warnings(run.err)
    cb.check(listed > 0 and "unassigned" in warn and "per-base" not in run.err,
             "unmux warns about tag pairs and the unassigned share, and leaves out no per-base tag")
    lima_out = sorted(w.glob("lima.bc*.bam"))
    lima_tags = cb.sam_records(ctx, w, lima_out[0].name)[0][1]
    pacbio_tags = sorted(t for t in inputs[0][1] if t != "RG")
    cb.check(all(t in lima_tags for t in pacbio_tags), "lima keeps the input tags")
    front = sorted(a for a, _ in cuts)[len(cuts) // 2]
    rear = sorted(b for _, b in cuts)[len(cuts) // 2]
    missing_text = ("" if not missing else
                    f"One read, of {len(reads[missing[0]])} bases, is missing from lima's report." if len(missing) == 1 else
                    f"{cb.fmt(len(missing))} reads are missing from lima's report.")
    row = cb.concord(lima_version, {z: unmux_calls.get(z) for z in reads}, {z: lima.get(z) for z in reads},
                     f"lima's score filter ({cb.fmt(len(unmux_only))}); one end beyond unmux's budget ({cb.fmt(len(lima_only))})")
    values = {
        "reads": cb.fmt(len(reads)),
        "samples": str(samples),
        "links": cb.source_links(ctx, ("pacbio-cyp2d6", "prefix.bam"), ("pacbio-cyp2d6", "biosamples.csv"),
                                 ("pacbio-cyp2d6", "Sequel_96_barcodes_v1.fasta")),
        "fetch": cb.fence(cb.fetch_command(ctx, "pacbio-cyp2d6", "prefix.bam"), "bash"),
        "prep": cb.fence(runs["prep"].command, "bash"),
        "heads_cmd": cb.fence(runs["heads"].command, "bash"),
        "heads": cb.fence(runs["heads"].out, "tsv"),
        "lima_version": lima_version,
        "lima": cb.fence(runs["lima"].command, "bash"),
        "demux": cb.fence(run.command, "bash"),
        "listed": str(listed),
        "warn": cb.fence(warn, "text"),
        "lima_n": cb.fmt(lima_n),
        "same": cb.fmt(table["same"]),
        "lima_only": cb.fmt(table["lima only"]),
        "unmux_only": cb.fmt(table["unmux only"]),
        "neither": cb.fmt(table["neither"]),
        "unmux_only_reads": cb.fmt(len(unmux_only)),
        "score_low": str(scores[0]),
        "score_high": str(scores[-1]),
        "lima_only_reads": cb.fmt(len(lima_only)),
        "lima_only_ends": cb.tally([(sides["one"], "at one end only"), (sides["both"], "at both ends"),
                                    (sides["none"], "at neither end")], len(lima_only)),
        "neither_reads": cb.fmt(len(neither)),
        "neither_ends": cb.tally([(left["one"], "at one end only"), (left["none"], "at neither end")]),
        "lima_short": cb.fmt(lima_short),
        "both_ends": cb.fmt(left["both"]),
        "both_kinds": cb.tally([(twice, "the same barcode twice"), (unlisted, "a pair no sample has")]),
        "missing": missing_text,
        "front": str(front),
        "rear": str(rear),
        "tag_len": str(tag_len),
        "pacbio_tags": ", ".join(pacbio_tags),
    }
    return values, [row]


cb.main("pacbio", body)
