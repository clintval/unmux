Use this to store raw read pairs as one unmapped BAM (uBAM), the input many alignment workflows such as GATK's expect, with a read group on every record.

**Data:** the first {{pairs}} read pairs of SRA run {{run_link}}, fetched with fg-sra:

{{fetch}}

{{convert}}

The log reports what unmux read and wrote:

{{log}}

{{view_cmd}}

{{view}}

Every record carries `RG:Z:SRR390728`. The pool, which becomes the read group's ID, SM and LB, defaults to the common stem of the input names; set it with `--pool`.

Each of the {{records}} records has the name, bases and qualities of its FASTQ read, flags 77 and 141 for the two reads of a pair, and no tag but RG.

Compared with `samtools import`, the records differ only in MAPQ and the read group:

{{compare_cmd}}

{{compare}}

> [!NOTE]
> MAPQ is left out of the comparison: unmux writes 255 (unavailable) on unmapped records, where `samtools import` writes 0. All {{records}} records match otherwise.
