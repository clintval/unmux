Use this to see why a read went where it did. The `--qc-tag` option writes a JSON summary of the match into a tag, ZS unless another is given, in the FASTQ comment or on the SAM record.

**Data:** the pooled lane of {{ref:bins}}.

{{command}}

A fig2b read whose i7 had one mismatch, a fig2c read, a PhiX read, and the fig2b read the PhiX rule took with one substitution:

{{show_cmd}}

{{show}}

Each matched group lists the tag it matched, where (`loc` is file, start and end), the observed bases, and the substitutions and indels it took. The fig2c read {{why}}, so no sample claims it and its reason is `no_sample`; a removed read gives the `--remove` rule that took it, and `str` marks a match on the reverse strand. The last record is how to find out why a read was removed by mistake. The generator checks each record's outcome against the file it is in.
