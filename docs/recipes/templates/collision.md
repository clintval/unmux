Before it reads a record, unmux checks the tags of each group against each other. A set in which one read could match two tags equally well is refused, with the ways out.

**Data:** the 10x and SPLiT-seq inputs from {{ref:tenx}} and {{ref:splitseq}}.

The 10x command with `dist=1` and no mode. Cell barcodes are not designed to be far apart, so some differ in only two places, and a read with one error between them matches both:

{{command}}

{{error}}

Under the default `mode=all` the first match would silently win. With `mode=nearest`, as in the 10x recipe, such reads are left uncorrected and the check only warns, listing up to ten pairs.

The SPLiT-seq command with one indel allowed (`dist=1:1:2`). Its barcodes are 4 substitutions apart, but an insertion or deletion can bring a read within reach of two:

{{indel}}

{{error2}}

Here unmux asks for a mode. With `mode=nearest` such reads are left unassigned, and with `mode=all` each is still assigned by a tie-break, with warnings that list the pairs; {{ref:wide}} compares the two. For both commands the generator checks the exit status, the message, and that no file was written.

> [!TIP]
> For a set designed to be far apart, such as Illumina indexes, a refusal usually means `dist` is too large for the set.
