Widening a window lets a barcode that an insertion or deletion shifted still be found. But in a window wider than its tags, a read with an error in its barcode can also match another tag, shifted, as closely as its own. unmux estimates how often that happens and warns.

**Data:** the SPLiT-seq reads of {{ref:splitseq}}, of which {{missing_round}} miss a round in the README windows.

The README command with each window 1 bp wider on both sides, one indel allowed, and `mode=all` as the collision check of {{ref:collision}} asks:

{{wide_all}}

{{warned}}

The run also prints {{pairs}} warnings that list tag pairs one indel apart, as `mode=all` does. Following the advice, the same command with `mode=nearest` in each round:

{{wide_nearest}}

{{quiet}}

Real reads have no known truth, so the table scores each run two ways: against the 568 cells split-pipe called, whose barcodes a correct CB most likely belongs to, and against the README windows, which take no shifted match:

{{table}}

The wider windows find {{more_cells}} more reads of called cells, but under `mode=all` they also change the CB of {{changed_all}} reads that the README windows assign, and add {{extra_other}} reads whose CB is no called cell's. Under `mode=nearest` fewer CBs change ({{changed_nearest}}), but {{missing}} reads that tie in some round are written with no CB at all:

{{count_cmd}}

{{counted}}

> [!WARNING]
> Under `mode=nearest`, a read that ties between two tags of a round is still written, without the CB built from that round, though every round sets `minFindsPerGroup=1`; its `--qc-tag` lists only the rounds that matched. Keep only reads with a CB, for example with `samtools view -d CB`.
