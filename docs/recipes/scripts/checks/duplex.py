"""Random duplex UMIs into RX: every record's RX, QX and bases against its read pair, and the swapped UMI pairs counted."""

import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb
import fastquorum


def swapped(umi_pair):
    a, b = umi_pair.split("-")
    return f"{b}-{a}"


def body(ctx, w, runs, templates):
    structure = fastquorum.read_structure(ctx.fetch("fastquorum", "samplesheet.tiny.csv"), "SRR6109255_1.fastq.gz")
    cb.check(structure == "10M1S+T 10M1S+T", f"the test sheet gives the read structure 10M1S+T 10M1S+T ({structure})")
    r1 = list(cb.read_fastq(w / "SRR6109255_1.fastq.gz"))
    r2 = list(cb.read_fastq(w / "SRR6109255_2.fastq.gz"))
    records = cb.sam_records(ctx, w, "SRR6109255.unmapped.bam")
    cb.check(len(records) == 2 * len(r1) == 2 * len(r2), "one record for each read of each pair")
    pairs = collections.Counter()
    for k, (a, b) in enumerate(zip(r1, r2)):
        (f1, t1), (f2, t2) = records[2 * k], records[2 * k + 1]
        rx, qx = f"{a[1][:10]}-{b[1][:10]}", f"{a[2][:10]} {b[2][:10]}"
        cb.check(f1[0] == f2[0] == a[0].split()[0] and int(f1[1]) & 64 and int(f2[1]) & 128,
                 f"{a[0].split()[0]}: the records are the pair, first then second")
        cb.check(t1["RX"] == t2["RX"] == rx and t1["QX"] == t2["QX"] == qx,
                 f"{f1[0]}: RX and QX hold the UMI of R1, then of R2 ({rx}, {qx})")
        cb.check((f1[9], f1[10], f2[9], f2[10]) == (a[1][11:], a[2][11:], b[1][11:], b[2][11:]),
                 f"{f1[0]}: each record holds its read after the UMI and skipped base")
        pairs[rx] += 1
    both = {p for p in pairs if swapped(p) in pairs and swapped(p) != p}
    count = runs["swap"]
    cb.check(int(count.out) == len(both), f"the shell pipeline counts the swapped UMI pairs ({count.out.strip()} and {len(both)})")
    spacers = collections.Counter(s[10] for _, s, _ in r1 + r2)
    files = [("fastquorum", n) for n in ("SRR6109255_1.fastq.gz", "SRR6109255_2.fastq.gz", "samplesheet.tiny.csv")]
    values = {
        "pairs": cb.fmt(len(r1)),
        "links": cb.source_links(ctx, *files),
        "structure": structure,
        "spacer_t": cb.pct(spacers["T"], len(r1) + len(r2)),
        "inputs_cmd": cb.fence(runs["inputs"].command, "bash"),
        "inputs": cb.fence(runs["inputs"].out, "text"),
        "duplex": cb.fence(runs["duplex"].command, "bash"),
        "view_cmd": cb.fence(runs["view"].command, "bash"),
        "view": cb.fence(runs["view"].out, "text"),
        "swapped": cb.fmt(len(both)),
        "umi_pairs": cb.fmt(len(pairs)),
        "molecules": cb.fmt(len(both) // 2),
        "swapped_pairs": cb.fmt(sum(pairs[p] for p in both)),
        "swap_cmd": cb.fence(count.command, "bash"),
        "swap": cb.fence(count.out, "text"),
    }
    return values, []


cb.main("duplex", body)
