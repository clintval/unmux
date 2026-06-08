Feature: Read routing, removal, and fan-out primitives
  Routing decides which bin each read lands in: the --remove skiplist (a bare group, or a narrower
  group::id combination, optionally rerouted to a destination with =PATTERN), --sample fan-out with the
  no-sample reads collected in --unassigned, the AND-across-groups (+) selector, and the fail-fast error
  when one tag is routed to more than one sample. Each scenario is built so the routing feature is
  load-bearing: removing or broadening it changes which bin a distinguishable read lands in, rather than
  passing on a default that would route the same way anyway. (The unassigned and removed bins hold raw
  input segments and never carry SAM tags, so non-routing is shown by the raw read landing - or not -
  in that bin.)

  Scenario: --remove of a group::id combination drops only that tag, not its siblings
    Given a file "reads.fq" containing:
      """
      @a1
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      @b1
      TTTTTTTTGGCCATGC
      +
      IIIIIIIIIIIIIIII
      @c1
      GGGGGGGGTTAACGAT
      +
      IIIIIIIIIIIIIIII
      """
    # The three reads carry sibling tags A, B, C of one group. --remove grp::C is narrow: it drops only
    # the C read while its siblings A and B pass through. Broadening it to the bare group `--remove grp`
    # would drop A and B as well, so asserting A still passes is load-bearing on the ::C qualifier.
    When I run `unmux reads.fq --group grp={A=AAAAAAAA,B=TTTTTTTT,C=GGGGGGGG} --remove grp::C`
    Then the exit code is 0
    And stdout contains "AAAAAAAACCCCGGGG"
    And stdout contains "TTTTTTTTGGCCATGC"
    And stdout does not contain "GGGGGGGGTTAACGAT"

  Scenario: --remove of a bare group drops every read matching any of its tags
    Given a file "reads.fq" containing:
      """
      @x1
      AAAAAATAGCATGCA
      +
      IIIIIIIIIIIIIII
      @x2
      GGGGGGTAGCATGCT
      +
      IIIIIIIIIIIIIII
      @y1
      TTTTTTGATCGATCG
      +
      IIIIIIIIIIIIIII
      @y2
      CCCCCCGATCGATCC
      +
      IIIIIIIIIIIIIII
      """
    # grp_y has two tags Y1 and Y2; --remove grp_y (bare) drops both members, while grp_x reads pass.
    # Narrowing the selector to grp_y::Y1 would let the Y2 read reappear, so asserting the Y2 read is
    # absent is load-bearing on the bare-group form claiming the whole group.
    When I run `unmux reads.fq --group grp_x={X1=AAAAAA,X2=GGGGGG} --group grp_y={Y1=TTTTTT,Y2=CCCCCC} --remove grp_y`
    Then the exit code is 0
    And stdout contains "AAAAAATAGCATGCA"
    And stdout contains "GGGGGGTAGCATGCT"
    And stdout does not contain "TTTTTTGATCGATCG"
    And stdout does not contain "CCCCCCGATCGATCC"

  Scenario: --remove SEL=PATTERN reroutes the removed read to a destination file
    Given a file "reads.fq" containing:
      """
      @a1
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      @c1
      GGGGGGGGTTAACGAT
      +
      IIIIIIIIIIIIIIII
      """
    # With =removed.fq the removed C read is written there as its raw segment (the A read still passes
    # to stdout). Dropping the =removed.fq destination (a bare --remove grp::C) drops C silently and
    # never creates the file, so asserting removed.fq holds the C bases is load-bearing on =PATTERN.
    When I run `unmux reads.fq --group grp={A=AAAAAAAA,C=GGGGGGGG} --remove grp::C=removed.fq`
    Then the exit code is 0
    And stdout contains "AAAAAAAACCCCGGGG"
    And stdout does not contain "GGGGGGGGTTAACGAT"
    And a file "removed.fq" exists
    And the file "removed.fq" contains "GGGGGGGGTTAACGAT"

  Scenario: the unassigned bin collects the raw read that matches no sample
    Given a file "reads.fq" containing:
      """
      @d1a
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      @d2a
      TTTTTTTTACGTACGT
      +
      IIIIIIIIIIIIIIII
      """
    # grp has tags D1 and D2 but only D1 is routed to a sample, so the D2 read matches no sample and its
    # raw segment lands in un.fq (un.fq carries no tags, so non-assignment is shown by the raw read
    # landing there, not by a missing tag). Widening the selector to grp::D1,D2 would route D2 to S1 and
    # remove it from un.fq, so asserting the D2 read in un.fq is load-bearing on the selector membership.
    When I run `unmux reads.fq --group grp={D1=AAAAAAAA,D2=TTTTTTTT} --sample S1=grp::D1 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "S1.fq" contains "AAAAAAAACCCCGGGG"
    And the file "un.fq" contains "TTTTTTTTACGTACGT"
    And the file "S1.fq" does not contain "TTTTTTTTACGTACGT"

  Scenario: an AND-across-groups (+) selector routes only the read matching both groups
    Given a file "reads.fq" containing:
      """
      @both1
      AAAAAATTTTTTCAGTGCAC
      +
      IIIIIIIIIIIIIIIIIIII
      @xonly
      AAAAAAGGGGGGCAGTGCAC
      +
      IIIIIIIIIIIIIIIIIIII
      @yonly
      CAGTGCTTTTTTAAGGCCTT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    # The + selector requires BOTH grp_x::X (AAAAAA) and grp_y::Y (TTTTTT). Only both1 carries both, so
    # it alone routes to both.fq; the X-only and Y-only reads match a single group each and fall to the
    # unassigned bin. Dropping the +grp_y::Y term (selecting just grp_x::X) would route the X-only read
    # into both.fq, so asserting the X-only raw read in un.fq is load-bearing on the AND term.
    When I run `unmux reads.fq --group grp_x={X=AAAAAA} --group grp_y={Y=TTTTTT} --sample both=grp_x::X+grp_y::Y --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "both.fq" contains "AAAAAATTTTTTCAGTGCAC"
    And the file "un.fq" contains "AAAAAAGGGGGGCAGTGCAC"
    And the file "un.fq" contains "CAGTGCTTTTTTAAGGCCTT"
    And the file "both.fq" does not contain "AAAAAAGGGGGGCAGTGCAC"

  Scenario: a tag routed to more than one sample is a fail-fast error
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAAAACAGT
      +
      IIIIIIIIIIII
      """
    # c2 is listed in both dna01 and dna02, so the same tag would route to two fan-out targets. The check
    # runs at parse time (before any read) and names the colliding tag and the two samples. Making the
    # selectors disjoint (dna01=c1, dna02=c2,c3) removes the overlap and the run succeeds, so the error
    # is load-bearing on the shared c2 membership.
    When I run `unmux reads.fq --group grp={c1=AAAAAAAA,c2=TTTTTTTT,c3=GGGGGGGG} --sample dna01=grp::c1,c2 --sample dna02=grp::c2,c3`
    Then the exit code is 1
    And stderr contains "more than one sample"
    And stderr contains "c2"
    And stderr contains "dna01"
    And stderr contains "dna02"

  Scenario: minFindsPerGroup plus no --unassigned expresses keeplisting (there is no --keep flag)
    Given a file "reads.fq" containing:
      """
      @m
      AAAAAACGCGCGCG
      +
      IIIIIIIIIIIIII
      @n
      TTTTTTCGCGCGCG
      +
      IIIIIIIIIIIIII
      """
    # There is no --keep flag: keeplisting is minFindsPerGroup=1 (the group must match at least once)
    # plus the absence of --unassigned (so the failing read is dropped, not collected). @m carries the
    # X1 tag and survives to stdout; @n has no grp_x tag, fails the group, and is dropped. Without
    # minFindsPerGroup the body template would still pass @n through, so the constraint is load-bearing.
    When I run `unmux reads.fq --group grp_x={X1=AAAAAA} --group grp_x::loc=0:0:6,minFindsPerGroup=1 --extract body=0:0:end --template body`
    Then the exit code is 0
    And stdout contains "AAAAAACGCGCGCG"
    And stdout does not contain "TTTTTTCGCGCGCG"

  Scenario: a combined barcode is assembled with a multi-stream --tag (there is no --assign flag)
    Given a file "reads.fq" containing:
      """
      @r1
      AAAATTTTGGGG
      +
      IIIIIIIIIIII
      @r2
      AAAACCCCGGGG
      +
      IIIIIIIIIIII
      """
    # There is no --assign flag: a combined barcode is just a multi-stream --tag over per-group extracts.
    # @r1 matches both groups so CB joins them as AAAA-TTTT; @r2 matches grp_a only, so its CB degrades to
    # the single AAAA stream. The sep=- is load-bearing on the join, and the @r2 result shows the missing
    # grp_b stream simply drops out rather than erroring.
    When I run `unmux reads.fq --group grp_a={A=AAAA} --group grp_a::loc=0:0:4 --group grp_b={B=TTTT} --group grp_b::loc=0:4:8 --extract bc_a=@grp_a --extract bc_b=@grp_b --extract body=0:0:end --template body --tag CB=bc_a+bc_b --tag CB::sep=- --tag CB::qual=none`
    Then the exit code is 0
    And stdout contains "@r1 CB:Z:AAAA-TTTT"
    And stdout contains "@r2 CB:Z:AAAA"
    And stdout does not contain "@r2 CB:Z:AAAA-TTTT"

  Scenario: a literal sub_sample and the %pool placeholder set different read-group LB values
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAAAACGCG
      +
      IIIIIIIIIIII
      """
    # The sub_sample is the @RG LB and the %sub_sample path component. A literal s1::sub_a writes
    # s1.sub_a.bam with LB:sub_a; s1::%pool resolves %sub_sample to the pool id, writing s1.lib01.bam with
    # LB:lib01. Same selector, same read - only the sub_sample binding changes both the file name and the
    # LB, so the two runs are distinguishable on each.
    When I run `unmux reads.fq --group grp={A=AAAAAAAA} --group grp::loc=0:0:8 --extract body=0:0:end --template body --sample s1::sub_a=grp::A --out litsub/%sample.%sub_sample.bam`
    Then the exit code is 0
    And a file "litsub/s1.sub_a.bam" exists
    And the BAM header of "litsub/s1.sub_a.bam" contains "LB:sub_a"
    And the BAM header of "litsub/s1.sub_a.bam" does not contain "LB:lib01"
    When I run `unmux reads.fq --pool lib01 --group grp={A=AAAAAAAA} --group grp::loc=0:0:8 --extract body=0:0:end --template body --sample s1::%pool=grp::A --out poolsub/%sample.%sub_sample.bam`
    Then the exit code is 0
    And a file "poolsub/s1.lib01.bam" exists
    And the BAM header of "poolsub/s1.lib01.bam" contains "LB:lib01"
    And the BAM header of "poolsub/s1.lib01.bam" does not contain "LB:sub_a"
