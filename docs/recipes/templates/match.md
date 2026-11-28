Use this to keep only reads from cells you have called and write each cell's name: a `match=` group looks up the three corrected round barcodes, joined, in a list of cells, and a column of the list gives the name.

**Data:** the SPLiT-seq reads of {{ref:splitseq}} and the {{called}} cells split-pipe called in that library ({{geo_link}}). Each cell's `cell_barcode` is its round-3 and round-2 barcodes, then its round-1 well, and the well's oligo-dT and random-hexamer barcodes are rows w and w + 48 of the round-1 list, so the list gets a row for each, with the same `cell`:

{{cells_cmd}}

{{cells}}

{{match}}

{{view_cmd}}

{{view}}

{{kept}} reads were kept, of {{names}} cells, and {{other}} went to other.unmapped.bam, both records of each pair as they were read. Checked against the generator's own lookup of the three windows: the kept reads are exactly those whose joined barcode is a called cell's, each with that cell's name, both primers of a cell under one name. {{corrected}} of them needed a correction in some round and are kept because `match=` joins each round as corrected, the way `--tag` writes it.

unmux warns, as in {{ref:splitseq}}, that nothing is split by sample; that CB, a column of cells.tsv, has no base qualities to write to CY, which `CB::qual=none` silences; and that {{share}} of the reads were unassigned, the {{other}} in other.unmapped.bam:

{{warnings}}

> [!TIP]
> Without `minFindsPerGroup=1` on the `match=` group, reads of other barcodes are written too, without a CB.
