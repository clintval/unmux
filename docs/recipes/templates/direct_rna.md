Use this to store direct RNA reads, which basecallers can write with U for uracil, in an unmapped BAM.

**Data:** a dorado test read written with U, with SAM tags in its FASTQ comment ({{links}}).

{{head_cmd}}

{{head}}

{{convert}}

{{view_cmd}}

{{view}}

BAM has no code for U, so all {{u}} U bases are stored as T, as `samtools import` stores them. The tags in the comment become SAM tags (more in {{ref:comment_tags}}), and RG:Z is replaced by the pool, as for any input.

SAM and FASTQ output keep the U as read:

{{text_cmd}}

{{text}}

> [!TIP]
> When matching, U in a read or a tag counts as T, so DNA barcodes find RNA reads written with U, on either strand. The generator checks this with a 24-mer of this read's reverse complement written with T.
