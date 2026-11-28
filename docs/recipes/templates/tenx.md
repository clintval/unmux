Use this to tag 10x Chromium reads with their cell barcode and UMI in an unmapped BAM, correcting each barcode to a list of known cells: CB holds the corrected barcode, CR and CY the barcode as read with its qualities, and UR, UY and UB the UMI.

**Data:** the {{reads}} 10x v2 read pairs of the nf-core/scrnaseq test set (R1: 16 bp barcode and 10 bp UMI; R2: cDNA) and the 10x v2 barcode keeplist ({{links}}).

The keeplist holds {{keeplist}} barcodes, and this recipe uses a subset of it: the {{cells}} within one substitution of one of the {{seen}} barcodes seen in R1. A read can be corrected only to such a barcode, so the subset gives every read the call the whole list would. unmux reads tag files with `id` and `seq` columns, so the barcode is both:

{{cells_cmd}}

{{tag}}

A read whose barcode was corrected and one with no keeplist barcode within one substitution, in file order:

{{view_cmd}}

{{view}}

Of {{reads}} reads, {{exact}} matched a keeplist barcode exactly, {{corrected}} were corrected, and {{none}} have no CB but keep CR, UR and UB: {{far}} with no keeplist barcode within one substitution{{tied}} On every read, CB agrees with a lookup of the barcode and its one-substitution neighbors in the whole keeplist, CR, CY, UR, UB and UY equal the R1 bases and qualities they come from, and the record holds R2.

unmux warns about tags one read can match equally well, which `mode=nearest` leaves uncorrected, and that no reads are split by sample, which is expected here:

{{warnings}}

> [!NOTE]
> unmux does not correct UMIs, so UB holds the same bases as UR; tools such as velocyto read UB.

> [!WARNING]
> Without `--pool S10` the pool, and with it the read group, would be `S10_L001_R`: the common stem of the two file names.
