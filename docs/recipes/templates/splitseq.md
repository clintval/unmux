Use this for split-pool data such as SPLiT-seq, where the cell barcode is three 8 bp barcodes, one per round of barcoding, separated by fixed linkers in R2. This is the README's Quick Start command.

**Data:** {{reads}} read pairs of SRA run {{run_link}}, the C2C12_short_1k library of LR-Split-seq (Rebboah et al. 2021) made with the Parse Biosciences v1 kit, whose R2 has the chemistry v2 layout of the split-seq-pipeline: a 10 bp UMI, round 3, linker 1, round 2, linker 2 and round 1 (`{{amp}}`, where N is the UMI and I a barcode). The barcodes are the pipeline's published lists ({{links}}).

{{fetch}}

{{r2}}

Make the two tag files from the published lists:

{{tags_cmd}}

{{tags}}

{{demux}}

{{view_cmd}}

{{view}}

CB joins the three corrected barcodes with `-`, and CY their qualities with spaces; RX and QX hold the UMI. A read that misses any round fails `minFindsPerGroup=1` and is dropped.

{{written}} of {{reads}} reads were written. The generator looks up each read's three windows itself, a published barcode within one substitution in each, and every written read has the CB, CY, RX and QX that lookup gives, while every read it drops misses a round; each record holds its R1. {{in_cells}} of the written reads belong to {{kept_cells}} of the {{called}} cells that split-pipe called in this library ({{geo_link}}).

{{warnings}}

The warning is expected: this command tags reads and does not split them by sample.

> [!NOTE]
> The Parse kit primes each round-1 well with an oligo-dT and a random-hexamer barcode, so one cell's reads carry either, and this command gives the cell two CBs; {{both}} called cells appear under both here. {{ref:match}} joins them through a cell list.

> [!NOTE]
> Without `--pool`, the pool is the common stem of the input names, here `R` from R1 and R2, so every record gets `RG:Z:R`; add `--pool` to set it.

**Compared with {{splitcode}}.** splitcode's SPLiT-seq tutorial ({{rst_link}}) uses the same layout as this command: rounds 1, 2 and 3 at R2 bases 78 to 86, 48 to 56 and 10 to 18, one substitution, and the UMI in R2's first 10 bases. Its config.txt for barcode reformatting, with its barcode files ({{splitcode_links}}):

{{config}}

{{splitcode_cmd}}

As the tutorial gives the command, with one `-o` file, splitcode 0.31.6 stops with `Error: Incorrect number of output files`, since it counts the `-o` files against `--nFastqs` even with `--select=0`. A second name, whose file stays empty, is the only change. The file output_R2.fastq.gz then holds each kept read's corrected rounds 3, 2 and 1.

The tutorial's round-1 files hold {{tutorial1}} barcodes, of which {{shared}} are among the {{library1}} of this library's list, while its round-2 and round-3 file is this library's list. A second run keeps the configuration and puts this library's oligo-dT and random-hexamer barcodes in r1_T.txt and r1_R.txt:

{{lists_cmd}}

Per read, against unmux:

{{agreement}}

With this library's lists, splitcode makes unmux's call on every read: {{library_same}} with a CB and {{library_neither}} with none. With the tutorial's files, of the {{unmux_only}} reads only unmux calls, {{outside}} have a round-1 barcode outside the tutorial's lists{{two}}; the {{other_only}} only splitcode calls hold one of the tutorial's round-1 barcodes, {{other_exact}} of them exactly, more than one substitution from any of this library's; and the {{differ}} calls that differ agree on rounds 2 and 3, while round 1 is one substitution from this library's barcode and {{differ_round1}}, which splitcode takes.

How unmux expresses splitcode's settings:

| splitcode | unmux |
| --- | --- |
| tag files such as `r1_R.txt$`, one tag per line, in group `round1` | `--group "round1=round1.tags.txt"`, one tag file for both primers |
| `distances` 1 | `dist=1` |
| `locations` 1,78,86 | `loc=1:78:86` |
| `minFindsG` 1 for round 1, and 2 for rounds 2 and 3 together | `minFindsPerGroup=1` on each round |
| `@extract 1:0<umi[10]>` | `--extract "umi=1:0:10"` and `--tag "RX=umi"` |
| `@extract <output_R2{*}>`, the corrected tags stitched together | `--tag "CB=bc1+bc2+bc3"` |
| `--select=0`, writing R1 only | `--template "cdna"` |
| `--min-delta` | `mode=nearest`, with `delta` |
| `--keep-grp`, files of tag-group arrangements to keep | `--sample-sheet` |
| `--subs`, such as to turn R barcodes into T barcodes | a `match=` group with a column, as in {{ref:match}} |
| `--revcomp` | `bothStrands=true` |
| `next` and `previous`, with the bases between | `next=` |
| `--partial5` and `--partial3` | `partial5` and `partial3` |
| `--qtrim` and `--filter-len` | no equivalent |
| `--assign` and `--mapping`, an ID for each combination of tags | no equivalent |

What splitcode does that unmux does not: it trims bases by quality and filters reads by length, and it can give each combination of tags a short ID with a mapping file. It also nests configurations and remultiplexes reads.
