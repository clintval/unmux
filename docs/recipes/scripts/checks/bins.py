"""Undetermined reads and a PhiX bin: every pair of a lane pooled from three runs against the run it came from."""

import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb
import illumina


def phix_reads(ctx, names_and_reads):
    """Names of read pairs with an exact PhiX 24-mer, either strand, in R1."""
    phix = cb.read_fasta(ctx.fetch("ncbi", "NC_001422.1.fasta"))
    kmers = {s[i:i + 24] for s in (phix, cb.revcomp(phix)) for i in range(len(s) - 23)}
    return {n for n, seq in names_and_reads if any(seq[i:i + 24] in kmers for i in range(len(seq) - 23))}


def body(ctx, w, runs, templates):
    origin, files = illumina.bins_pool(ctx)
    removed = {name.split()[0] for name, _, _ in cb.read_fastq(w / "removed" / "PhiX.0.fastq.gz")}
    undetermined = {name.split()[0] for name, _, _ in cb.read_fastq(w / "fastq" / "Undetermined.0.fastq.gz")}
    assigned = {}
    for path in sorted((w / "fastq").glob("N*.R1.fastq.gz")):
        for name, _, _ in cb.read_fastq(path):
            assigned[name.split()[0]] = path.name.split(".")[0]
    phix = phix_reads(ctx, [(n.split()[0], s) for n, s, _ in files["R1"]])
    with open(w / "sample_info.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    i1 = [(n.split()[0], s, q) for n, s, q in cb.read_fastq(ctx.fetch("uditas", "Undetermined_S0_L001_I1_001.fastq"))]
    i2 = [(n.split()[0], s, q) for n, s, q in cb.read_fastq(ctx.fetch("uditas", "Undetermined_S0_L001_I2_001.fastq"))]
    expected = illumina.uditas_calls(i1, i2, rows, masked=False)
    by_origin = collections.Counter(origin.values())
    removed_by = collections.Counter(origin[n] for n in removed)
    cb.check(len(assigned) + len(removed) + len(undetermined) == len(origin), "each pair is in exactly one output")
    cb.check(removed_by["phix-run"] == by_origin["phix-run"], "every pair of the PhiX run is in the PhiX bin")
    cb.check(all(origin[n] == "uditas" and expected[n] == s for n, s in assigned.items()), "fig2b reads go to their sheet sample")
    cb.check(all(origin[n] == "uditas-fig2c" for n in undetermined), "Undetermined holds only fig2c pairs")
    cb.check(len(undetermined) == by_origin["uditas-fig2c"] - removed_by["uditas-fig2c"], "every fig2c pair not removed is undetermined")
    cb.check(len(assigned) == by_origin["uditas"] - removed_by["uditas"], "every fig2b pair not removed is assigned")
    kept_out = sorted(n for n in removed if origin[n] != "phix-run")
    exact = sum(1 for n in kept_out if n in phix)
    exact_removed = collections.Counter(origin[n.split()[0]] for n, _, _ in cb.read_fastq(w / "exact" / "removed.PhiX.0.fastq.gz"))
    strict_removed = w / "strict" / "removed.PhiX.0.fastq.gz"
    removed_strict = {n.split()[0] for n, _, _ in cb.read_fastq(strict_removed)} if strict_removed.exists() else set()
    undetermined_strict = {n.split()[0] for n, _, _ in cb.read_fastq(w / "strict" / "Undetermined.0.fastq.gz")}
    run_pairs = {n for n, o in origin.items() if o == "phix-run"}
    cb.check(run_pairs <= undetermined_strict and not run_pairs & removed_strict,
             "with minFindsPerGroup on the index groups, the PhiX run's pairs go to Undetermined, not the PhiX bin")
    truth = {n: expected[n] if o == "uditas" and expected[n] not in ("mismatched_adapters", "not_in_exp_list") else None
             for n, o in origin.items()}
    truth_row = cb.concord("the run of each pair, and the UDiTaS rules within fig2b", {n: assigned.get(n) for n in origin}, truth,
                           f"fig2b pairs that match the PhiX k-mers ({cb.fmt(len(kept_out))})")
    values = {
        "fig2b": cb.fmt(by_origin["uditas"]),
        "fig2c": cb.fmt(by_origin["uditas-fig2c"]),
        "fig2c_link": cb.source_links(ctx, ("uditas-fig2c", "Undetermined_S0_L001_R1_001.fastq")),
        "phix_run": cb.fmt(by_origin["phix-run"]),
        "run_link": cb.link("SRR27862880", ctx.sources[("sra-SRR27862880", "SRR27862880_1.fastq")].url),
        "fetch": cb.fence(cb.fetch_command(ctx, "sra-SRR27862880", "SRR27862880_1.fastq"), "bash"),
        "phix_link": cb.link("NC_001422.1", ctx.sources[("ncbi", "NC_001422.1.fasta")].url),
        "kmers": cb.fence(runs["kmers"].command, "bash"),
        "demux": cb.fence(runs["demux"].command, "bash"),
        "listing_cmd": cb.fence(runs["listing"].command, "bash"),
        "listing": cb.fence(runs["listing"].out, "text"),
        "unassigned_cmd": cb.fence(runs["unassigned"].command, "bash"),
        "unassigned": cb.fence(runs["unassigned"].out, "tsv"),
        "kept_out": cb.fmt(len(kept_out)),
        "exact": cb.fmt(exact),
        "near": cb.fmt(len(kept_out) - exact),
        "undetermined": cb.fmt(len(undetermined)),
        "assigned": cb.fmt(len(assigned)),
        "exact_phix": cb.fmt(exact_removed["phix-run"]),
        "exact_others": cb.fmt(sum(exact_removed.values()) - exact_removed["phix-run"]),
    }
    return values, [truth_row]


cb.main("bins", body)
