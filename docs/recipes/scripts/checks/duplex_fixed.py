"""Fixed duplex UMIs corrected to a list: RX, OX and BZ against the nearest list UMI, then dist=0, dist=2 and --unassigned."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb
import fastquorum


def body(ctx, w, runs, templates):
    sheet = ctx.fetch("fastquorum", "samplesheet.nonrandom.csv")
    structure = fastquorum.read_structure(sheet, "synthetic_nonrandom_umi_1.fastq.gz")
    cb.check(structure == "8M+T 8M+T", f"the test sheet gives the read structure 8M+T 8M+T ({structure})")
    umis = [row["seq"] for row in cb.read_tsv(w / "umis.tsv")]
    r1 = list(cb.read_fastq(w / "synthetic_nonrandom_umi_1.fastq.gz"))
    r2 = list(cb.read_fastq(w / "synthetic_nonrandom_umi_2.fastq.gz"))
    for made, end in ((r1, 1), (r2, 2)):
        real = list(cb.read_fastq(ctx.fetch("fastquorum", f"SRR6109255_{end}.fastq.gz")))
        cb.check([m[0] for m in made] == [r[0] for r in real], f"the R{end} reads are those of SRR6109255, in order")
        cb.check(all(m[1][8:] == r[1][11:] and m[2][8:] == r[2][11:] and m[2][:8] == "I" * 8 for m, r in zip(made, real)),
                 f"each R{end} read is 8 bases of quality I and the SRR6109255 read after its first 11 bases")
    records = cb.sam_records(ctx, w, "umis.unmapped.bam")
    cb.check(len(records) == 2 * len(r1), "one record for each read of each pair")

    def nearest(seen):
        best = sorted((cb.hamming(seen, u), u) for u in umis)
        cb.check(best[0][0] <= 1 and best[1][0] > best[0][0],
                 f"the UMI {seen} is within one substitution of one list UMI only ({best[:2]})")
        return best[0][1]

    fixed = [0, 0]
    pick = None
    for k, (a, b) in enumerate(zip(r1, r2)):
        (f1, t1), (f2, t2) = records[2 * k], records[2 * k + 1]
        n1, n2 = nearest(a[1][:8]), nearest(b[1][:8])
        fixed[0] += n1 != a[1][:8]
        fixed[1] += n2 != b[1][:8]
        for t in (t1, t2):
            cb.check(t["RX"] == f"{n1}-{n2}" and t["OX"] == f"{a[1][:8]}-{b[1][:8]}" and "QX" not in t,
                     f"{f1[0]}: RX holds the corrected UMIs, OX the UMIs as read, and no QX is written")
            cb.check(t["BZ"] == f"{a[2][:8]} {b[2][:8]}", f"{f1[0]}: BZ holds the qualities of the UMIs as read")
        cb.check((f1[9], f2[9]) == (a[1][8:], b[1][8:]), f"{f1[0]}: each record holds its read after the UMI")
        if pick is None and n1 != a[1][:8]:
            pick = f1[0]
    view = ctx.sh(w, templates["view"].replace("{{pick}}", pick))
    counted = runs["counted"]
    either = sum(1 for a, b in zip(r1, r2) if a[1][:8] not in umis or b[1][:8] not in umis)
    exact_records = cb.sam_records(ctx, w, "umis.exact.bam")
    cb.check(all("OX" in t and "BZ" in t for _, t in exact_records), "with dist=0 every record still has OX and BZ")
    cb.check(counted.out.split() == [str(len(r1)), str(len(r1) - either)],
             f"with dist=0 every pair is written and only those with both UMIs on the list have RX ({counted.out.split()})")
    cb.check("can both match" in runs["wide"].err, "dist=2 is refused as a collision")
    kept = cb.sam_records(ctx, w, "umis.kept.bam", "-f 64")
    set_aside = cb.sam_records(ctx, w, "umis.aside.bam", "-f 64")
    cb.check(len(set_aside) == either and len(kept) == len(r1) - either and all("RX" in t for _, t in kept),
             "with minFindsPerGroup=1 the pairs without a corrected UMI go to --unassigned")
    names = ("synthetic_nonrandom_umi_1.fastq.gz", "synthetic_nonrandom_umi_2.fastq.gz", "idt-xgen-cfdna-umis.txt",
             "samplesheet.nonrandom.csv")
    values = {
        "links": cb.source_links(ctx, *(("fastquorum", n) for n in names)),
        "pairs": cb.fmt(len(r1)),
        "umis": str(len(umis)),
        "corrected": cb.fmt(fixed[0] + fixed[1]),
        "reads": cb.fmt(2 * len(r1)),
        "structure": structure,
        "inputs_cmd": cb.fence(runs["inputs"].command, "bash"),
        "inputs": cb.fence(runs["inputs"].out, "text"),
        "prep_cmd": cb.fence(runs["prep"].command, "bash"),
        "prep": cb.fence(runs["prep"].out, "tsv"),
        "fixed": cb.fence(runs["fixed"].command, "bash"),
        "view_cmd": cb.fence(view.command, "bash"),
        "view": cb.fence(view.out, "text"),
        "corrected_r1": cb.fmt(fixed[0]),
        "corrected_r2": cb.fmt(fixed[1]),
        "exact": cb.fence(runs["exact"].command, "bash"),
        "counted_cmd": cb.fence(counted.command, "bash"),
        "counted": cb.fence(counted.out, "text"),
        "either": cb.fmt(either),
    }
    return values, []


cb.main("duplex_fixed", body)
