"""Shared helpers of the SPLiT-seq recipes: their staged inputs, the read layout, and the lookups that score unmux."""

import csv
import gzip
import re

import cookbook as cb

GEO_CELLS = "GSM5169184_C2C12_short_1k_cell_metadata.csv.gz"
INPUTS = [
    ("sra-SRR13948565", "SRR13948565_1.fastq", "R1.fastq.gz"),
    ("sra-SRR13948565", "SRR13948565_2.fastq", "R2.fastq.gz"),
    ("splitseq", "bc_8nt_v1.csv", "bc_8nt_v1.csv"),
    ("splitseq", "bc_8nt_v2.csv", "bc_8nt_v2.csv"),
    ("geo", GEO_CELLS, GEO_CELLS),
]
STAGED = [path for _, _, path in INPUTS] + ["round1.tags.txt", "round2_3.tags.txt", "cells.tsv"]

SPLITSEQ_TAGS = """awk -F , 'BEGIN { print "id\\tseq" } { print "bc" $1 "\\t" $2 }' bc_8nt_v2.csv > round1.tags.txt
awk -F , 'BEGIN { print "id\\tseq" } { print "bc" $1 "\\t" $2 }' bc_8nt_v1.csv > round2_3.tags.txt"""

SPLITSEQ_CELLS = """gzip -dc GSM5169184_C2C12_short_1k_cell_metadata.csv.gz | awk -F , 'NR == FNR { rt[$1] = $2; next }
  FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; print "id\\tseq\\tcell"; next }
  { cb = $c["cell_barcode"]; w = $c["rnd1_well"]; r3 = substr(cb, 1, 8); r2 = substr(cb, 9, 8)
    print cb "_dT\\t" rt[w] r2 r3 "\\t" cb; print cb "_hex\\t" rt[w + 48] r2 r3 "\\t" cb }' bc_8nt_v2.csv - > cells.tsv"""

SPLITCODE_BARCODES = {"tutorial": ("r1_R.txt", "r1_T.txt", "r2_r3.txt"), "library": ("r2_r3.txt",)}
SPLITCODE_STAGED = [
    f"splitcode-{run}/{name}"
    for run, barcodes in SPLITCODE_BARCODES.items()
    for name in ("R1.fastq.gz", "R2.fastq.gz", "bc_8nt_v2.csv", "config.txt", *barcodes)
]


def splitseq_inputs(ctx, w):
    """Make the round tag files and cells.tsv in w, which holds the files of INPUTS."""
    ctx.sh(w, SPLITSEQ_TAGS)
    ctx.sh(w, SPLITSEQ_CELLS)


def splitcode_inputs(ctx, w):
    """Give each splitcode run a directory in w; the run with this library's lists writes its own round-1 files."""
    config = splitcode_config(ctx.fetch("splitcode", "tutorials_splitseq.rst").read_text())
    for run, barcodes in SPLITCODE_BARCODES.items():
        d = w / f"splitcode-{run}"
        for name in ("R1.fastq.gz", "R2.fastq.gz", "bc_8nt_v2.csv"):
            cb.place(w / name, d / name)
        for name in barcodes:
            cb.place(ctx.fetch("splitcode", name), d / name)
        cb.write_text(d / "config.txt", config)


def splitcode_config(rst):
    """The config.txt of the tutorial's barcode reformatting example, from its reStructuredText source."""
    lines = rst.split("\n")
    start = lines.index("  :caption: config.txt")
    body = []
    for line in lines[start + 1:]:
        if not line.strip():
            if body:
                break
            continue
        body.append(line[3:])
    cb.check(body[0].startswith("@extract") and len(body) == 6, "the tutorial's config.txt has its @extract line, a header and four tags")
    return "\n".join(body) + "\n"


def splitseq_layout(ctx):
    tools = ctx.fetch("splitseq", "tools.py").read_text()
    block = tools.split("elif chemistry=='v2':", 1)[1]
    amp = re.search(r"amp_seq = '([NIACGT]+)'", block).group(1)
    parts = [p for p in re.split(r"(N+|I+)", amp) if p]
    linkers = [p for p in parts if p[0] not in "NI"]
    return linkers[0], linkers[1], amp


def splitseq_oracle(w, windows=((78, 86), (48, 56), (10, 18))):
    """Each read's CB, CY, RX and QX under the README windows: a unique tag within one substitution per round."""
    round1 = [r["seq"] for r in cb.read_tsv(w / "round1.tags.txt")]
    round23 = [r["seq"] for r in cb.read_tsv(w / "round2_3.tags.txt")]

    def fix(observed, tags):
        hits = [t for t in tags if len(observed) == len(t) and cb.hamming(observed, t) <= 1]
        return hits[0] if len(hits) == 1 else None

    calls = {}
    for name, seq, qual in cb.read_fastq(w / "R2.fastq.gz"):
        found = [fix(seq[a:b], tags) for (a, b), tags in zip(windows, (round1, round23, round23))]
        calls[name] = None if None in found else (
            "-".join(found), " ".join(qual[a:b] for a, b in windows), seq[:10], qual[:10])
    return calls


def called_cells(ctx, w):
    """The cells split-pipe called, by joined round 1, 2 and 3 barcode, checked against the published wells."""
    rows = cb.read_tsv(w / "cells.tsv")
    lines = (w / "bc_8nt_v1.csv").read_text().splitlines()
    v1 = {int(line.split(",")[0]): line.strip().split(",")[1] for line in lines if line.strip()}
    with gzip.open(w / GEO_CELLS, "rt") as fh:
        for row in csv.DictReader(fh):
            barcode = row["cell_barcode"]
            cb.check(barcode[:8] == v1[int(row["rnd3_well"]) + 1] and barcode[8:16] == v1[int(row["rnd2_well"]) + 1],
                     f"{barcode}: rounds 3 and 2 are the barcodes of its wells")
    return {r["seq"]: r["cell"] for r in rows}


def splitseq_scores(ctx, w, path, cells, reference=None):
    """Written reads, with a CB, with a called cell's CB, and with a CB other than the reference run's."""
    written = with_cb = in_cells = changed = 0
    for fields, tags in cb.sam_records(ctx, w, path):
        written += 1
        barcode = tags.get("CB")
        if barcode is None:
            continue
        with_cb += 1
        in_cells += barcode.replace("-", "") in cells
        if reference is not None and reference.get(fields[0]) not in (None, barcode):
            changed += 1
    return written, with_cb, in_cells, changed
