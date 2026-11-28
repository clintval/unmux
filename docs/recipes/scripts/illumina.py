"""Shared helpers of the Illumina recipes: the UDiTaS rules, fqtk's calls, and the pools of reads of known origin."""

import cookbook as cb

CARABUS_BARCODES = [
    "AAAAA", "AACCC", "AAGGG", "AATTT", "ACACG", "ACCAT", "ACGTA", "ACTGC",
    "CAACT", "CACAG", "CAGTC", "CATGA", "CCAAC", "CCCCA", "CCGGT", "CCTTG",
    "GAAGC", "GACTA", "GAGAT", "GATCG", "GCATT", "GCCGG", "GCGCC", "GCTAA",
    "TAATG", "TACGT", "TAGCA", "TATAC", "TCAGA", "TCCTC", "TCGAG", "TCTCT",
]
CARABUS = [(f"DRR0219{59 + k}", f"CR{k + 1:02d}", barcode) for k, barcode in enumerate(CARABUS_BARCODES)]


def uditas_calls(i1, i2, rows, masked=True):
    """Each read's call under the rules of UDiTaS's demultiplex(): sample, or the bin it goes to."""
    bc1 = [r["barcode_I1"] for r in rows]
    bc2 = [r["barcode_I2"] for r in rows]
    name1 = dict(zip(bc1, (r["index_I1"] for r in rows)))
    name2 = dict(zip(bc2, (r["index_I2"] for r in rows)))
    pairs = {(r["index_I1"], r["index_I2"]) for r in rows}

    def mask(seq, qual):
        return "".join(b if ord(q) - 33 >= 12 else "N" for b, q in zip(seq, qual)) if masked else seq

    def closest(seq, barcodes):
        near = [(cb.hamming(seq, b), b) for b in barcodes if cb.hamming(seq, b) <= 1]
        best = min((d for d, _ in near), default=None)
        return [b for d, b in near if d == best]

    calls = {}
    for (name, s1, q1), (_, s2, q2) in zip(i1, i2):
        s1, s2 = mask(s1, q1), mask(s2[:8], q2[:8])
        if not (s1 in bc1 and s2 in bc2):
            m1, m2 = closest(s1, bc1), closest(s2, bc2)
            if not (m1 and m2):
                calls[name] = "mismatched_adapters"
                continue
            s1, s2 = m1[0], m2[0]
        pair = (name1[s1], name2[s2])
        calls[name] = f"{pair[0]}_{pair[1]}" if pair in pairs else "not_in_exp_list"
    return calls


def fqtk_calls(directory):
    """Each read's sample from fqtk's per-sample FASTQs, or None for its unmatched reads."""
    calls = {}
    for path in sorted(directory.glob("*.R1.fq.gz")):
        sample = path.name[: -len(".R1.fq.gz")]
        for name, _, _ in cb.read_fastq(path):
            calls[name.split()[0]] = None if sample == "unmatched" else sample
    return calls


def bins_pool(ctx):
    """The lane: UDiTaS fig2b and fig2c reads, then PhiX pairs with index reads built from their spot group."""
    origin = {}
    files = {r: [] for r in ("R1", "R2", "I1", "I2")}
    for key in ("uditas", "uditas-fig2c"):
        for r, records in files.items():
            for name, seq, qual in cb.read_fastq(ctx.fetch(key, f"Undetermined_S0_L001_{r}_001.fastq")):
                records.append((name, seq, qual))
                origin[name.split()[0]] = key
    mates = [list(cb.read_fastq(ctx.fetch("sra-SRR27862880", f"SRR27862880_{m}.fastq"))) for m in (1, 2)]
    for (n1, s1, q1), (n2, s2, q2) in zip(*mates):
        name, comment = n1.split(" ", 1)
        i7, i5 = comment.split(":")[3].split("+")
        files["R1"].append((n1, s1, q1))
        files["R2"].append((n2, s2, q2))
        files["I1"].append((name, i7, "F" * len(i7)))
        files["I2"].append((name, i5, "F" * len(i5)))
        origin[name] = "phix-run"
    cb.check(len(origin) == len(files["R1"]), "read names are unique across the pooled runs")
    return origin, files


def stage_bins(ctx, w):
    _, files = bins_pool(ctx)
    for r, records in files.items():
        cb.write_fastq(w / f"{r}.fastq.gz", records)


def inline_pool(ctx):
    """The 32 RAD-seq runs one after another, and the sample of the run each read came from."""
    origin, records = {}, []
    for run, sample, _ in CARABUS:
        for name, seq, qual in cb.read_fastq(ctx.fetch("sra-carabus", f"{run}.fastq")):
            records.append((name, seq, qual))
            origin[name] = sample
    cb.check(len(origin) == len(records), "read names are unique across the 32 runs")
    return origin, records


def stage_inline(ctx, w):
    _, records = inline_pool(ctx)
    cb.write_fastq(w / "pool.fastq.gz", records)
    cb.write_tsv(w / "barcodes.tsv", ["id", "seq"], [(sample, barcode) for _, sample, barcode in CARABUS])
