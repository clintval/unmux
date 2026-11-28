Use this to get an unmapped copy of aligned reads, to realign them or share them without the alignment, as `samtools reset` does.

**Data:** two aligned Nanopore reads with modified-base calls (MM/ML), one on each strand ({{links}}).

{{convert}}

The input's header without its {{sq}} `@SQ` lines, the first six columns of its records, and the tags of the first record:

{{before_cmd}}

{{before}}

The same for the output:

{{after_cmd}}

{{after}}

What changes:

- Each record becomes unmapped: flag 4, no position or CIGAR, MAPQ 255. The read on the reverse strand (flag 16) gets back its sequenced orientation, bases and qualities.
- Tags that describe the alignment are dropped, here {{dropped}}; per-read tags such as MM, ML and the basecaller's tags are kept.
- The header keeps `@HD VN` and gets the pool's read group and an unmux `@PG` line. The input `@RG`, `@PG`, `@CO` and `@SQ` lines are not kept, and each record's RG:Z becomes the pool.

Apart from MAPQ and RG:Z, the records equal those of `samtools reset`:

{{compare_cmd}}

{{compare}}

The generator also checks that MM and ML are unchanged on both records and that the reverse-strand read is the reverse complement of its aligned SEQ.

> [!TIP]
> To keep the alignments and split them into files instead, use a read-group split, as in {{ref:split_rg}}.
