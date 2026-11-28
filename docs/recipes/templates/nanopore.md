Use this to demultiplex Nanopore reads from a native barcoding kit, where each read starts with the barcode and ends with its reverse complement, trimming the barcodes off.

**Data:** 200 reads of each of 24 SRA runs from one MinION flow cell of SARS-CoV-2 ARTIC amplicons ({{project_link}}, batch VERD6, kits SQK-LSK109 with EXP-NBD104 and EXP-NBD114). Guppy split the reads into one run per native barcode, requiring the barcode at both ends and leaving it in the reads, so the barcode of every read is known. The {{reads}} reads are pooled one run after another:

{{fetch}}

Each tag is a native barcode between its flanks, from dorado ({{kit_links}}), and the sheet asks for the same barcode at both ends:

{{tags_cmd}}

{{tags}}

{{demux}}

The `front` group searches the first 150 bases and the `rear` group the last 150 on both strands, so it finds the reverse-complemented copy. A budget of 12 edits in a 46 bp tag suits the error rate of these R9.4.1 reads, and `mode=nearest` leaves a read that two barcodes fit equally well unclassified. The span `@front..@rear` keeps the bases between the two tags.

{{counts_cmd}}

{{counts}}

With a budget this large, unmux warns about pairs of tags that one read could match, listing {{listed}} for each group and counting the rest:

{{warn}}

The windows are much wider than the tags, as long reads need, and unmux prints no window warning: it warns about one only where it estimates that over 0.3% of reads with an error could match another tag at a shifted position.

{{correct}} of {{reads}} reads ({{share}}) went to the barcode of their run and {{wrong}} to another; {{unclassified}} are unclassified. Every trimmed read is part of its read, cut within the search windows: a median of {{front_cut}} bases off the front and {{rear_cut}} off the rear.

> [!TIP]
> With `dist=8:8:8`, {{tight_right}} reads go to the barcode of their run, {{tight_wrong}} to another, and {{tight_unclassified}} are unclassified: a tighter budget costs reads here without saving any from a wrong barcode.

What a trim does to the basecaller's tags, on the ten modkit reads of {{ref:comment_tags}}, before and after a 20-base trim:

{{mods_cmd}}

{{mods}}

{{left_out}}

> [!WARNING]
> A `--template` record keeps the tags of its read, but per-base tags such as MM, ML and the move table mv describe the read's bases one by one, so unmux leaves them out of a record that does not copy its read exactly, with one warning for the run. Demultiplex without a template to keep the modification calls; the barcodes then stay in the reads.
