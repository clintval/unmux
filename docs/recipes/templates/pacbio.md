Use this to demultiplex PacBio HiFi reads that still carry their barcodes, one at each end, as lima does, here with a different barcode at each end of a sample's amplicons.

**Data:** the first {{reads}} HiFi reads of PacBio's CYP2D6 amplicon set, {{samples}} samples each with a pair of barcodes, before lima ran on them, cut from the first 10 MB of the BAM by an HTTP range request ({{links}}). Each read starts with a pad of a few bases, a barcode and the M13 forward or reverse primer, and ends with the same for the other barcode, reverse-complemented.

{{fetch}}

Keep the complete records of the cut, and build a tag of each barcode with its M13 primer and a sheet with each pair in both orders:

{{prep}}

{{heads_cmd}}

{{heads}}

{{lima_version}}, with its HiFi preset for a different barcode at each end:

{{lima}}

and unmux:

{{demux}}

unmux warns about pairs of tags that one read could match within this budget, listing {{listed}} for each group, and about the share of reads it leaves unassigned:

{{warn}}

Per read, against lima's report ({{lima_n}} reads with a sample):

| Call | Reads |
| --- | ---: |
| Same sample from lima and unmux | {{same}} |
| A sample from lima only | {{lima_only}} |
| A sample from unmux only | {{unmux_only}} |
| No sample from either | {{neither}} |

The two never call different samples for one read. lima scores each barcode by alignment, and its ASYMMETRIC preset removes a read whose combined score is 80 or less: the {{unmux_only_reads}} reads only unmux assigns score from {{score_low}} to {{score_high}}, though both of their ends fit one sheet row within unmux's budget. unmux needs a tag within its budget at both ends, and in the {{lima_only_reads}} reads only lima assigns, its `--qc-tag` found a tag {{lima_only_ends}}.

Of the {{neither_reads}} reads neither assigns, unmux found a tag {{neither_ends}}, and lima passed fewer than two barcode regions in {{lima_short}}. The other {{both_ends}} carry tags at both ends but not a sample's pair: {{both_kinds}}, and lima reports the same barcodes. {{missing}}

Every trimmed read is part of its read, cut within the search windows: a median of {{front}} bases off the start and {{rear}} off the end, each a {{tag_len}}-base tag and its pad.

How unmux expresses lima's settings:

| Setting | lima | unmux |
| --- | --- | --- |
| The same barcode at both ends | `--hifi-preset SYMMETRIC`, which sets `--same` | a sheet row with that barcode in both groups, as in {{ref:nanopore}} |
| A different barcode at each end | `--hifi-preset ASYMMETRIC`, which sets `--different` | a sheet row for the pair in each order, since a read can come from either strand |
| The reverse-complemented barcode at the far end | searched by default | `bothStrands=true` on the end group |
| Where to search | `-W` in bases, or `-w` as a multiple of the barcode length | `loc=0:0:N` for the start and `loc=0:-N:end` for the end |
| How well a barcode must fit | `--min-score` and `--min-end-score`, alignment scores out of 100 | `dist`, a budget of substitutions, indels and edits in all |
| How far the best must lead | `--min-score-lead` | `mode=nearest`, with `delta` for a larger lead |
| Sample names and files | `--biosample-csv`, and `--split-named` for a file per barcode pair | `--sample-sheet`, and `--out "%sample.bam"` for a file per sample |
| Barcode removal | clipped by default, kept with `--no-clip` | cut with `--template` and `@start..@end`, kept with no template |
| Read filters | `--min-ref-span`, `--min-scoring-regions` and `--min-length` | no equivalent |
| Finding which barcodes are in use | `--peek-guess` | no equivalent: the sheet lists them |
| Per-read reports | `bc` and `bq` tags, and the `.lima.report` file | `--qc-tag` JSON, and the metrics files |

> [!NOTE]
> Both tools keep the tags PacBio wrote ({{pacbio_tags}}) on the reads they clip, and lima adds its own. unmux sets RG to the sample's read group.
