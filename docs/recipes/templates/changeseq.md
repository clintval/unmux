Use this to turn read pairs into the single reads that CIRCLE-seq and CHANGE-seq align: the pipeline's `mergeReads()` writes the reverse complement of R1 followed by R2, with R1's qualities reversed, so both reads of a cleaved site line up on one strand.

**Data:** the treated ({{treated}} pairs) and control ({{control}} pairs) test libraries of the CHANGE-seq pipeline ({{links}}).

{{inputs_cmd}}

{{inputs}}

{{merge}}

{{view_cmd}}

{{view}}

The `~` reverse-complements the R1 stream, qualities reversed with it, and `+` joins the two streams into one read. Every merged read of both libraries equals what {{merge_link}} writes for its pair, bases and qualities, and keeps the pair's name without its comment.
