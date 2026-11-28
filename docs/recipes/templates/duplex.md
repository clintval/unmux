Use this for duplex sequencing libraries with a random UMI at the start of each read: one RX tag gets both UMIs, R1's first, the UMIs and spacers come off the reads, and the unmapped BAM is ready for fgbio's duplex workflow.

**Data:** the {{pairs}} read pairs of SRR6109255 in the nf-core/fastquorum test set ({{links}}). The sample sheet gives the read structure `{{structure}}`: in each read a 10 bp UMI and a 1 bp spacer, then the template. The spacer is a T in {{spacer_t}} of reads.

{{inputs_cmd}}

{{inputs}}

{{duplex}}

{{view_cmd}}

{{view}}

On all {{pairs}} pairs, RX is R1's first 10 bases, a hyphen, then R2's, QX holds their qualities joined by a space, and both reads start after their spacer, at base 11 counting from 0. The pool, and so the read group, is the common stem of the input names, `SRR6109255`.

In this sample {{swapped}} of the {{umi_pairs}} UMI pairs also appear swapped, as B-A for an A-B: {{molecules}} molecules read from both strands, in {{swapped_pairs}} read pairs. These are the pairs the paired strategy joins:

{{swap_cmd}}

{{swap}}

> [!IMPORTANT]
> Keep R1's UMI first. The read pairs from one strand of a molecule carry A-B and those from the other strand B-A, since their R1 starts from the other end. fgbio's paired strategy joins the two by that swap and tells the strands apart by it; sorting the two UMIs, or putting R2's first, would lose the strand.

> [!TIP]
> Next, align the reads, bring the tags back onto the alignments with fgbio ZipperBams, and run fgbio GroupReadsByUmi with `--strategy paired`, then CallDuplexConsensusReads.
