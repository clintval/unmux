Use the three metrics files to check a run: how many reads each sample got, how many matched no sample and why, and the totals.

**Data:** the pooled lane of {{ref:bins}}, whose command wrote all three.

{{summary_cmd}}

{{summary}}

A read is a template: one pair here. Every read is assigned, passed through (when no sample is declared), unassigned or removed, and `frac_assigned` is the assigned share. The `mean_read_len_by_input` column gives the mean length of each input in order, which catches an index read shorter than its tags.

{{per_cmd}}

{{per}}

The `frac_bases_unextracted` column is the share of a sample's input bases left out of its output, here the index reads.

{{unassigned_cmd}}

{{unassigned}}

The per-sample reads add up to `reads_assigned` ({{assigned}}), the reasons to `reads_unassigned` ({{un}}), and with the {{removed}} removed they make up all {{total}} pairs of the pool; each sample's count, the removed count and the unassigned count also equal the pairs in its files. A reason of `no_sample` means no sample claimed what the groups matched; a group with `minFindsPerGroup` adds rows that give the group and the bound it missed.
