"""SPLiT-seq barcodes in three rounds: every record against a lookup of its three windows, and the calls against splitcode."""

import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb
import splitseq


def splitcode_calls(d, names):
    """Each read's CB from splitcode's stitched barcodes (rounds 3, 2 and 1), or None for a read it drops."""
    windows = {n.split()[0]: s for n, s, _ in cb.read_fastq(d / "R2.fastq.gz")}
    calls = dict.fromkeys(names)
    for name, seq, _ in cb.read_fastq(d / "output_R2.fastq.gz"):
        n = name.split()[0].split("::")[0]
        r3, r2, r1 = seq[:8], seq[8:16], seq[16:24]
        cb.check(len(seq) == 24 and all(cb.hamming(a, b) <= 1 for a, b in ((r3, windows[n][10:18]), (r2, windows[n][48:56]), (r1, windows[n][78:86]))),
                 f"{n}: splitcode's barcode is rounds 3, 2 and 1, each within one substitution of its window")
        calls[n] = f"{r1}-{r2}-{r3}"
    return calls


def splitcode_comparison(ctx, w, runs, mine):
    version = cb.tool_version(runs["splitcode_version"], "splitcode")
    as_written = runs["as_written"]
    cb.check("Incorrect number of output files" in as_written.err + as_written.out,
             "splitcode refuses the tutorial's command as written")
    calls = {}
    for label in ("tutorial", "library"):
        d = w / f"splitcode-{label}"
        cb.check(not list(cb.read_fastq(d / "unused_R2.fastq.gz")), "the second -o file stays empty")
        calls[label] = splitcode_calls(d, mine)
    tutorial1 = [s for n in ("r1_R.txt", "r1_T.txt") for s in ctx.fetch("splitcode", n).read_text().split()]
    library1 = {r["seq"] for r in cb.read_tsv(w / "round1.tags.txt")}
    cb.check(set(ctx.fetch("splitcode", "r2_r3.txt").read_text().split()) == {r["seq"] for r in cb.read_tsv(w / "round2_3.tags.txt")},
             "the tutorial's round-2 and round-3 list is this library's")
    shared = len(set(tutorial1) & library1)
    theirs = calls["tutorial"]
    window = {n.split()[0]: s[78:86] for n, s, _ in cb.read_fastq(w / "R2.fastq.gz")}
    unmux_only = [n for n in mine if mine[n] and not theirs[n]]
    outside = sum(1 for n in unmux_only if mine[n].split("-")[0] not in tutorial1)
    two = sum(1 for n in unmux_only if mine[n].split("-")[0] in tutorial1 and sum(cb.hamming(window[n], b) <= 1 for b in tutorial1) > 1)
    cb.check(outside + two == len(unmux_only),
             f"each read only unmux calls has a round-1 barcode outside the tutorial's lists or two within reach ({outside}, {two}, {len(unmux_only)})")
    other_only = [n for n in mine if theirs[n] and not mine[n]]
    cb.check(all(theirs[n].split("-")[0] not in library1 and min(cb.hamming(window[n], b) for b in library1) > 1 for n in other_only),
             "each read only splitcode calls holds a tutorial round-1 barcode more than one substitution from this library's")
    exact = sum(1 for n in other_only if window[n] == theirs[n].split("-")[0])
    differ = [n for n in mine if mine[n] and theirs[n] and mine[n] != theirs[n]]
    cb.check(all(mine[n].split("-")[1:] == theirs[n].split("-")[1:] and theirs[n].split("-")[0] not in library1 for n in differ),
             "the calls that differ differ in round 1, where splitcode takes a tutorial barcode")
    differ_exact = sum(1 for n in differ if window[n] == theirs[n].split("-")[0])
    tutorial = cb.concord(f"{version}, tutorial files", mine, theirs,
                          f"{shared} of the tutorial's {len(tutorial1)} round-1 barcodes are this library's")
    library = cb.concord(f"{version}, this library's lists", mine, calls["library"], "")
    cb.check(library["different"] == library["unmux"] == library["other"] == 0, "with this library's lists, splitcode makes unmux's calls")
    values = {
        "splitcode": version,
        "rst_link": cb.link("tutorials_splitseq.rst", ctx.sources[("splitcode", "tutorials_splitseq.rst")].url),
        "splitcode_links": cb.source_links(ctx, ("splitcode", "r1_R.txt"), ("splitcode", "r1_T.txt"), ("splitcode", "r2_r3.txt")),
        "config": cb.fence((w / "splitcode-tutorial" / "config.txt").read_text(), "text"),
        "splitcode_cmd": cb.fence(runs["tutorial"].command, "bash"),
        "tutorial1": cb.fmt(len(tutorial1)),
        "shared": str(shared),
        "library1": str(len(library1)),
        "lists_cmd": cb.fence(runs["lists"].command, "bash"),
        "agreement": cb.agreement_table("splitcode", [("tutorial files", tutorial), ("this library's lists", library)]),
        "library_same": cb.fmt(library["same"]),
        "library_neither": cb.fmt(library["neither"]),
        "unmux_only": cb.fmt(len(unmux_only)),
        "outside": cb.fmt(outside),
        "two": f" and {cb.fmt(two)} have round-1 bases within one substitution of two of its barcodes" if two else "",
        "other_only": cb.fmt(len(other_only)),
        "other_exact": cb.fmt(exact),
        "differ": cb.fmt(len(differ)),
        "differ_round1": f"is a tutorial barcode exactly in {cb.fmt(differ_exact)}" if differ_exact else "within one of a tutorial barcode",
    }
    return values, [tutorial, library]


def body(ctx, w, runs, templates):
    _, _, amp = splitseq.splitseq_layout(ctx)
    run = runs["demux"]
    oracle = splitseq.splitseq_oracle(w)
    r1 = {n: (s, q) for n, s, q in cb.read_fastq(w / "R1.fastq.gz")}
    records = cb.sam_records(ctx, w, "splitseq.unmapped.bam")
    for fields, tags in records:
        want = oracle[fields[0]]
        cb.check(want is not None, f"{fields[0]}: written only when every round matches")
        cb.check((tags["CB"], tags["CY"], tags["RX"], tags["QX"]) == want, f"{fields[0]}: CB, CY, RX and QX match the oracle")
        cb.check((fields[9], fields[10]) == r1[fields[0]], f"{fields[0]}: the record is R1")
    cb.check(len(records) == sum(1 for v in oracle.values() if v), "every read the oracle calls is written")
    cells = splitseq.called_cells(ctx, w)
    kept_cells = collections.Counter(cells[t["CB"].replace("-", "")] for _, t in records if t["CB"].replace("-", "") in cells)
    by_primer = collections.defaultdict(set)
    for _, t in records:
        joined = t["CB"].replace("-", "")
        if joined in cells:
            by_primer[cells[joined]].add(joined)
    both = sum(1 for seqs in by_primer.values() if len(seqs) == 2)
    r2 = next(cb.read_fastq(w / "R2.fastq.gz"))
    reads = len(oracle)
    mine = dict.fromkeys(oracle)
    for fields, tags in records:
        mine[fields[0]] = tags["CB"]
    splitcode_values, splitcode_rows = splitcode_comparison(ctx, w, runs, mine)
    values = {
        "reads": cb.fmt(reads),
        "run_link": cb.link("SRR13948565", ctx.sources[("sra-SRR13948565", "SRR13948565_1.fastq")].url),
        "amp": amp,
        "links": cb.source_links(ctx, ("splitseq", "bc_8nt_v2.csv"), ("splitseq", "bc_8nt_v1.csv"), ("splitseq", "tools.py")),
        "fetch": cb.fence(cb.fetch_command(ctx, "sra-SRR13948565", "SRR13948565_1.fastq"), "bash"),
        "r2": cb.fence(f"@{r2[0]}\n{r2[1]}\n+\n{r2[2]}", "text"),
        "tags_cmd": cb.fence(runs["tags"].command, "bash"),
        "tags": cb.fence(runs["tags"].out, "tsv"),
        "demux": cb.fence(run.command, "bash"),
        "view_cmd": cb.fence(runs["view"].command, "bash"),
        "view": cb.fence(runs["view"].out, "text"),
        "written": cb.fmt(len(records)),
        "in_cells": cb.fmt(sum(kept_cells.values())),
        "kept_cells": cb.fmt(len(kept_cells)),
        "called": cb.fmt(len(set(cells.values()))),
        "geo_link": cb.link("GSM5169184", ctx.sources[("geo", splitseq.GEO_CELLS)].url),
        "warnings": cb.fence(cb.log_lines(run.err, ["WARN"]), "text"),
        "both": cb.fmt(both),
        **splitcode_values,
    }
    return values, splitcode_rows


cb.main("splitseq", body)
