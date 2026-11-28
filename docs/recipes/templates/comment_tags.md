Use this when reads went through FASTQ with their SAM tags in the name comment, as `samtools fastq -T` and basecallers write them, and the tags, MM and ML among them, should be SAM tags again.

**Data:** ten unmapped Nanopore reads with modified-base calls ({{links}}), written to FASTQ with every tag in the comment:

{{prep}}

{{header}}

{{convert}}

{{view_cmd}}

{{view}}

On all {{reads}} reads the bases, qualities and every tag, MM and ML included, equal the original BAM's, except RG:Z, which becomes the pool as for any input. Even NM and MD come back: tags lifted from a FASTQ comment are kept as they are, while the same tags on a SAM or BAM input record are dropped as alignment tags.

> [!NOTE]
> A comment that is not a SAM tag, such as Illumina's `1:N:0:ACGT`, is not kept, and a malformed tag stops the run.
