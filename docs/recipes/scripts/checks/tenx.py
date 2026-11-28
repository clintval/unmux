"""10x Chromium cell barcodes and UMIs: CB against a lookup in the whole keeplist, and CR, CY, UR, UB and UY against R1."""

import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb

FILES = ("S10_L001_R1_001.fastq.gz", "S10_L001_R2_001.fastq.gz", "10xv2_whitelist.txt")


def body(ctx, w, runs, templates):
    run = runs["tag"]
    cells = [row["seq"] for row in cb.read_tsv(w / "cells.tsv")]
    keep = {line.strip() for line in (w / "10xv2_whitelist.txt").read_text().splitlines() if line.strip()}
    records = cb.sam_records(ctx, w, "S10.unmapped.bam")
    r1 = {name.split()[0]: (seq, qual) for name, seq, qual in cb.read_fastq(w / "S10_L001_R1_001.fastq.gz")}
    r2 = {name.split()[0]: (seq, qual) for name, seq, qual in cb.read_fastq(w / "S10_L001_R2_001.fastq.gz")}
    cb.check(len(records) == len(r1), "one record per read pair")

    def neighbors(observed):
        return {observed[:p] + b + observed[p + 1:] for p in range(16) for b in cb.BASES if b != observed[p]}

    seen = {seq[:16] for seq, _ in r1.values()}
    reach = {s for o in seen for s in neighbors(o) | {o} if s in keep}
    cb.check(len(cells) == len(set(cells)) and set(cells) == reach,
             "cells.tsv holds each keeplist barcode within one substitution of a barcode in R1, once")

    def oracle(observed):
        if observed in keep:
            return observed, "exact"
        hits = neighbors(observed) & keep
        return (hits.pop(), "corrected") if len(hits) == 1 else (None, "tied" if hits else "far")

    agree = exact = corrected = 0
    kinds = collections.Counter()
    mine, theirs = {}, {}
    pick_fixed = pick_none = None
    for fields, tags in records:
        seq, qual = r1[fields[0]]
        want, kind = oracle(seq[:16])
        kinds[kind] += 1
        got = tags.get("CB")
        agree += got == want
        mine[fields[0]], theirs[fields[0]] = got, want
        cb.check(tags["CR"] == seq[:16] and tags["CY"] == qual[:16], f"{fields[0]}: CR and CY are R1 bases 0 to 16")
        cb.check(tags["UR"] == tags["UB"] == seq[16:26] and tags["UY"] == qual[16:26], f"{fields[0]}: UR, UB and UY are R1 bases 16 to 26")
        cb.check((fields[9], fields[10]) == r2[fields[0]], f"{fields[0]}: the record is R2")
        if got is None:
            pick_none = pick_none or (fields[0] if kind == "far" else None)
        elif got == seq[:16]:
            exact += 1
        else:
            corrected += 1
            pick_fixed = pick_fixed or fields[0]
    cb.check(agree == len(records), f"CB agrees with the whole keeplist on every read ({agree} of {len(records)})")
    none = kinds["far"] + kinds["tied"]
    view_cmd = templates["view"].replace("{{fixed}}", str(pick_fixed)).replace("{{none}}", str(pick_none))
    view = ctx.sh(w, view_cmd)
    shown = [line for line in run.err.split("\n") if "WARN" in line]
    excerpt = "\n".join(shown[:3] + (["…"] + shown[-1:] if len(shown) > 4 else shown[3:]))
    values = {
        "reads": cb.fmt(len(records)),
        "links": cb.source_links(ctx, *(("nfcore", name) for name in FILES)),
        "keeplist": cb.fmt(len(keep)),
        "cells": cb.fmt(len(cells)),
        "seen": cb.fmt(len(seen)),
        "cells_cmd": cb.fence(runs["cells"].command, "bash"),
        "tag": cb.fence(run.command, "bash"),
        "view_cmd": cb.fence(view_cmd, "bash"),
        "view": cb.fence(view.out, "text"),
        "exact": cb.fmt(exact),
        "corrected": cb.fmt(corrected),
        "none": cb.fmt(none),
        "far": cb.fmt(kinds["far"]),
        "tied": f" and {cb.fmt(kinds['tied'])} with two or more." if kinds["tied"] else ".",
        "warnings": cb.fence(excerpt, "text"),
    }
    return values, [cb.concord("the whole 10x keeplist", mine, theirs, "")]


cb.main("tenx", body)
