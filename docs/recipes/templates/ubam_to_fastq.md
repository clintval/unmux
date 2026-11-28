Use this to hand reads stored in an unmapped BAM to a tool that wants FASTQ. The `%ordinal` placeholder writes one file per read of the pair.

**Data:** the unmapped BAM from {{ref:fastq_to_ubam}}.

{{convert}}

{{listing}}

{{view_cmd}}

{{view}}

The tags of each record, here the read group, go into the FASTQ comment. The files match `samtools fastq -T RG` once its tab before the tags is read as a space:

{{compare_cmd}}

{{compare}}

Bases and qualities equal the original FASTQ for all {{pairs}} pairs.
