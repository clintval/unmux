"""CITE-seq antibody barcodes: fx and fb against a lookup of the panel, and CR, CY, UR and UY against R1."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb

NAMES = (
    "subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz",
    "subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz",
    "5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv",
)


def body(ctx, w, runs, templates):
    panel = {row["id"]: row["seq"] for row in cb.read_tsv(w / "antibodies.tsv")}
    gap = min(cb.hamming(a, b) for a in panel.values() for b in panel.values() if a != b)
    r1 = {n.split()[0]: (s, q) for n, s, q in cb.read_fastq(w / NAMES[0])}
    r2 = {n.split()[0]: (s, q) for n, s, q in cb.read_fastq(w / NAMES[1])}

    def lookup(window):
        hits = [i for i, s in panel.items() if cb.hamming(window, s) <= 1 and len(window) == len(s)]
        return hits[0] if len(hits) == 1 else None

    want = {n: lookup(s[10:25]) for n, (s, _) in r2.items()}
    records = cb.sam_records(ctx, w, "adt.unmapped.bam")
    for fields, tags in records:
        (s1, q1), (s2, q2) = r1[fields[0]], r2[fields[0]]
        cb.check(tags["fx"] == want[fields[0]] and tags["fb"] == panel[tags["fx"]], f"{fields[0]}: fx and fb match the panel lookup")
        cb.check((tags["CR"], tags["CY"], tags["UR"], tags["UY"]) == (s1[:16], q1[:16], s1[16:], q1[16:]), f"{fields[0]}: CR, CY, UR, UY are R1")
        cb.check((fields[9], fields[10]) == (s2, q2), f"{fields[0]}: the record is R2")
    cb.check(len(records) == sum(1 for v in want.values() if v), "every read whose window holds a panel barcode is written")
    exact = sum(1 for n, (s, _) in r2.items() if s[10:25] in panel.values())
    values = {
        "reads": cb.fmt(len(r2)),
        "antibodies": str(len(panel)),
        "links": cb.source_links(ctx, *(("nfcore-modules", n) for n in NAMES)),
        "input_cmd": cb.fence(runs["input"].command, "bash"),
        "input": cb.fence(runs["input"].out, "text"),
        "prep_cmd": cb.fence(runs["prep"].command, "bash"),
        "prep": cb.fence(runs["prep"].out, "tsv"),
        "tag": cb.fence(runs["tag"].command, "bash"),
        "view_cmd": cb.fence(runs["view"].command, "bash"),
        "view": cb.fence(runs["view"].out, "text"),
        "count_cmd": cb.fence(runs["count"].command, "bash"),
        "count": cb.fence(runs["count"].out, "text"),
        "written": cb.fmt(len(records)),
        "exact": cb.fmt(exact),
        "gap": str(gap),
    }
    return values, []


cb.main("cite", body)
