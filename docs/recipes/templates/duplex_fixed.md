Use this for duplex kits with a fixed set of UMIs, such as IDT xGen cfDNA: each read's UMI is corrected to the list, RX gets the corrected pair, and OX keeps the bases as read with their qualities in BZ.

**Data:** the nf-core/fastquorum test set for fixed UMIs and the UMI list ({{links}}). Its {{pairs}} read pairs are those of {{ref:duplex}} with the first 11 bases of each read, the random UMI and the base after it, replaced by 8 bases of quality I: one of the {{umis}} IDT xGen cfDNA UMIs or, in {{corrected}} of the {{reads}} reads, one substitution away from one. The rest of each read is unchanged. The read structure is `{{structure}}`: an 8 bp UMI, then the template.

{{inputs_cmd}}

{{inputs}}

The list holds one UMI per line, so give it the `id` and `seq` columns unmux reads:

{{prep_cmd}}

{{prep}}

{{fixed}}

A pair whose R1 UMI had one substitution:

{{view_cmd}}

{{view}}

On all {{pairs}} pairs RX holds the two list UMIs and OX the bases as read, with BZ their qualities. {{corrected_r1}} R1 UMIs and {{corrected_r2}} R2 UMIs had one substitution and were corrected. Every UMI in this data is within one substitution of the list, so no pair is left uncorrected here.

To see what happens to a UMI that cannot be corrected, allow none (`dist=0`):

<details>
<summary>The command</summary>

{{exact}}

</details>

{{counted_cmd}}

{{counted}}

The {{either}} pairs with a substitution in either UMI are still written, with OX and BZ but no RX: a tag is left out when any stream in it is missing.

> [!TIP]
> To set such pairs aside instead, add `minFindsPerGroup=1` to both groups and give `--unassigned`.

> [!NOTE]
> The 32 UMIs are at least 3 substitutions apart, so `dist=1` corrects without ties. With `dist=2` unmux refuses the set, since one read could then match two UMIs equally well.
