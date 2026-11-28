"""Explain one read's assignment with --qc-tag: four reads of the PhiX-bin lane, each with the outcome of its file."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def body(ctx, w, runs, templates):
    barcodes = {row["id"]: row["seq"] for row in cb.read_tsv(w / "i7.tsv")}
    fixed = sample = None
    i1 = {n.split()[0]: s for n, s, _ in cb.read_fastq(w / "I1.fastq.gz")}
    for path in sorted((w / "fastq").glob("N*.R1.fastq.gz")):
        i7 = barcodes[path.name.split("_")[0]]
        hit = next((n.split()[0] for n, _, _ in cb.read_fastq(path) if cb.hamming(i1[n.split()[0]], i7) == 1), None)
        if hit:
            fixed, sample = hit, path.name.split(".")[0]
            break
    other = next(n.split()[0] for n, _, _ in cb.read_fastq(w / "fastq" / "Undetermined.0.fastq.gz"))
    removed = [n.split()[0] for n, _, _ in cb.read_fastq(w / "removed" / "PhiX.0.fastq.gz")]
    phix = next(n for n in removed if n.startswith("A00488"))
    qc = {h.split()[0]: json.loads(h.split("ZS:Z:", 1)[1]) for h, _, _ in cb.read_fastq(w / "removed" / "PhiX.0.fastq.gz")}
    stray = next(n for n in removed if n.startswith("M03958") and all(g["sub"] > 0 for g in qc[n]["groups"] if g["g"] == "phix"))
    show = (templates["show"].replace("{{sample}}", sample).replace("{{fixed}}", fixed).replace("{{other}}", other)
            .replace("{{phix}}", phix).replace("{{stray}}", stray))
    out = ctx.sh(w, show)
    lines = out.out.splitlines()
    cb.check(len(lines) == 4, "one header line per read shown")
    outcomes = [json.loads(line.split("ZS:Z:", 1)[1])["outcome"] for line in lines]
    cb.check(outcomes == ["assigned", "unassigned", "removed", "removed"], f"QC outcomes {outcomes}")
    cb.check(json.loads(lines[0].split("ZS:Z:", 1)[1])["sample"] == sample, "the assigned read's QC sample")
    matched = [g["g"] for g in json.loads(lines[1].split("ZS:Z:", 1)[1])["groups"]]
    why = {(): "matched neither index", ("i7",): "matched an i7 of the sheet but not its i5"}[tuple(matched)]
    values = {
        "command": cb.fence(runs["qc"].command, "bash"),
        "show_cmd": cb.fence(show, "bash"),
        "show": cb.fence(out.out, "text"),
        "why": why,
    }
    return values, []


cb.main("qc_tag", body)
