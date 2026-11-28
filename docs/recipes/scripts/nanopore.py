"""Inputs of the Nanopore recipe: the 24 VERD6 runs pooled into one FASTQ, and dorado's native barcodes as tags."""

import re

import cookbook as cb

VERD6 = [
    ("NB01", "SRR14560737"), ("NB02", "SRR14560738"), ("NB03", "SRR14560739"), ("NB04", "SRR14560652"),
    ("NB05", "SRR14560583"), ("NB06", "SRR14560584"), ("NB07", "SRR14560585"), ("NB08", "SRR14560586"),
    ("NB09", "SRR14560587"), ("NB10", "SRR14560588"), ("NB11", "SRR14560763"), ("NB12", "SRR14560764"),
    ("NB13", "SRR14560765"), ("NB14", "SRR14560523"), ("NB15", "SRR14560525"), ("NB16", "SRR14560631"),
    ("NB17", "SRR14560632"), ("NB18", "SRR14560633"), ("NB19", "SRR14560634"), ("NB20", "SRR14560635"),
    ("NB21", "SRR14560769"), ("NB22", "SRR14560772"), ("NB23", "SRR14560526"), ("NB24", "SRR14560531"),
]


def native_barcodes(ctx):
    text = ctx.fetch("dorado", "barcode_kits.cpp").read_text()
    front = re.search(r'NB_1st_FRONT = "([ACGT]+)"', text).group(1)
    rear = re.search(r'NB_1st_REAR = "([ACGT]+)"', text).group(1)
    codes = dict(re.findall(r'\{"(NB\d\d)", "([ACGT]+)"\}', text))
    return front, rear, [(f"NB{k:02d}", codes[f"NB{k:02d}"]) for k in range(1, 25)]


def pooled(ctx):
    """The reads of the 24 runs one run after another, and the barcode of the run each read came from."""
    origin, records = {}, []
    for barcode, run in VERD6:
        for name, seq, qual in cb.read_fastq(ctx.fetch("sra-verd6", f"{run}.fastq")):
            records.append((name, seq, qual))
            origin[name] = barcode
    return origin, records


def stage(ctx, w):
    """Write a tag of each native barcode between its flanks, a sheet of each barcode at both ends, and the pooled reads."""
    front, rear, barcodes = native_barcodes(ctx)
    cb.write_tsv(w / "native_barcodes.tsv", ["id", "seq"], [(i, front + s + rear) for i, s in barcodes])
    cb.write_tsv(w / "barcode_pairs.tsv", ["sample", "front", "rear"], [(i, i, i) for i, _ in barcodes])
    cb.write_fastq(w / "reads.fastq.gz", pooled(ctx)[1])
