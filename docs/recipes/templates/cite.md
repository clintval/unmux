Use this for antibody-derived tags (ADT) in CITE-seq or 10x feature barcoding: find the antibody barcode in R2, correct it, and write the antibody's name next to the cell barcode and UMI.

**Data:** the {{reads}} antibody read pairs of the 10x 5k anti-CMV T cell set in nf-core/test-datasets, with their feature reference of {{antibodies}} TotalSeq-C antibodies ({{links}}). R1 is a 16 bp cell barcode and the UMI; in R2 the 15 bp antibody barcode follows 10 other bases, as the reference's pattern `^NNNNNNNNNN(BC)NNNNNNNNN` gives:

{{input_cmd}}

{{input}}

Give the feature reference the `id` and `seq` columns unmux reads:

{{prep_cmd}}

{{prep}}

{{tag}}

{{view_cmd}}

{{view}}

The `@ab::id` stream is the matched antibody's `id` column, written as text to fx; fb holds the corrected barcode.

{{count_cmd}}

{{count}}

{{written}} of {{reads}} reads were written: every read whose R2 bases 10 to 25 hold a panel barcode within one substitution ({{exact}} of them exactly), checked against the generator's own lookup. On every record fx and fb match that lookup, CR, CY, UR and UY equal the R1 bases and qualities, and the record holds R2. The panel's barcodes are at least {{gap}} substitutions apart, so one substitution cannot reach a wrong antibody.

> [!NOTE]
> The tags fb and fx follow Cell Ranger's names for a feature barcode and its feature. Any two-letter tag works; lowercase ones are free for local use.
