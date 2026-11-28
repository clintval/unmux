Use this to keep the reads that match no sample, the way `bcl-convert` writes Undetermined files, and to move spike-in reads such as PhiX into a bin of their own.

**Data:** a lane pooled from three real runs, so the origin of every read is known: the {{fig2b}} fig2b reads of {{ref:uditas}}, the first {{fig2c}} reads of UDiTaS's fig2c run, a library with an i5 that fig2b's sheet does not list ({{fig2c_link}} and its other reads), and {{phix_run}} read pairs of the PhiX control run {{run_link}}, fetched with fg-sra:

{{fetch}}

SRA kept the PhiX run's index only in its spot group, `GGGGGGGG+AGATCTCG` (the i7 reads as G with no index on a two-color instrument), so the pool's I1 and I2 for those pairs are built from it, with every quality F. That is the only constructed part.

The PhiX tags are 24-mers taken every 12 bases of the PhiX174 genome, {{phix_link}}:

{{kmers}}

{{demux}}

The `phix` group looks for any of the 24-mers on either strand of R1. The `--remove` rule sends the reads it matches to their own files, and `--unassigned` keeps the rest that match no sample, one file per input since `%source` is the input's number from 0 (R1, R2, I1, I2).

{{listing_cmd}}

{{listing}}

{{unassigned_cmd}}

{{unassigned}}

Checked against the origin of every pair: the PhiX bin holds all {{phix_run}} pairs of the PhiX run and {{kept_out}} UDiTaS pairs, {{exact}} whose R1 holds an exact PhiX 24-mer and {{near}} with one a substitution away from a tag. Undetermined holds the other {{undetermined}} fig2c pairs and nothing else, and the samples got the other {{assigned}} fig2b pairs, each in the sample {{ref:uditas}} gives it.

> [!TIP]
> With `dist=0` in the `phix` group, the bin takes {{exact_phix}} PhiX pairs and {{exact_others}} others: an exact match spares those UDiTaS reads but misses PhiX reads in which every tag that fits has an error.

> [!WARNING]
> Do not add `minFindsPerGroup` to the index groups here. A read that fails a find constraint is unassigned before a `--remove` rule sees it, so with `minFindsPerGroup=1` on i7 and i5 all the PhiX run's pairs go to Undetermined and none to the PhiX bin.

> [!TIP]
> Pass the `--out` pattern itself to `--unassigned` to keep unassigned reads in the `--out` file instead of a file of their own.
