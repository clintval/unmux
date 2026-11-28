Use this to separate an aligned SAM or BAM that holds several read groups, such as lanes or samples, into one file per read group with the alignments kept whole.

**Data:** {{total}} aligned records in two read groups, some with no read group or one missing from the header ({{links}}).

{{split}}

{{counts_cmd}}

{{counts}}

{{header_cmd}}

{{header}}

Each file has only its own `@RG` line, and the records are those `samtools split` writes:

{{compare_cmd}}

{{compare}}

The headers differ only in where `@SQ` sits and the added unmux `@PG` line:

{{hdiff_cmd}}

{{hdiff}}

> [!TIP]
> With `--group "rg=@RG::SM"` the files are per sample, merging the read groups of one sample, and `@RG::SM::LB` adds a library level for `--out "%sample.%sub_sample.bam"`.
