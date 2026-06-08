Feature: Tag groups and sequential matching across multiple groups
  Demultiplexing often needs more than one tag group, with downstream groups depending on upstream
  ones. These scenarios exercise the real binary on small fixtures and are designed so the group
  feature under test is load-bearing: removing the link (next/prev), the count constraint
  (minFindsPerGroup/maxFindsPerGroup), or the second selector flips the outcome (which read assigns,
  which sample wins, or which span the anchor resolves), rather than passing on a whole-read or
  leftmost default. The unassigned bin holds raw input segments and never carries SAM tags, so
  non-assignment is shown by the raw read landing in the unassigned file, and routing into an
  eagerly-created per-sample file is confirmed with a --metrics-per-sample frac_of_pool assertion.

  Scenario: next=GROUP:lo-hi constrains the downstream search to a relative spacer window
    Given a file "reads.fq" containing:
      """
      @r1
      AACGTAGGGCCCAAGGGCCCTTTTGGGG
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # g1 (AACGTA) matches at [0,6). The g2 tag GGGCCC occurs twice: a decoy at [6,12) (followed by
    # AAGG) and the true copy at [14,20) (followed by TTTT). next=g2:8-20 searches g2 only in
    # [g1.end+8, g1.end+20) = [14,26), excluding the decoy, so @g2+0:4 reads TTTT. A whole-read g2
    # search takes the leftmost decoy at [6,12) and the +0:4 payload is AAGG instead.
    When I run `unmux reads.fq --group g1={s1=AACGTA} --group g1::loc=0:0:6,next=g2:8-20 --group g2={s2=GGGCCC} --extract pay=@g2+0:4 --extract body=0:0:end --template body --tag PY=pay`
    Then the exit code is 0
    And stdout contains "PY:Z:TTTT"
    And stdout does not contain "PY:Z:AAGG"

  Scenario: prev=GROUP makes an upstream match a prerequisite for searching the downstream group
    Given a file "reads.fq" containing:
      """
      @rhit
      AACGTAGGGCCCAAAA
      +
      IIIIIIIIIIIIIIII
      @rmiss
      TTTTTTGGGCCCAAAA
      +
      IIIIIIIIIIIIIIII
      """
    # Both reads carry the g2 tag GGGCCC at [6,12). Only rhit carries the g1 tag AACGTA at [0,6).
    # g2::prev=g1 skips g2 for rmiss (g1 absent), so rmiss is unassigned and its raw segment lands in
    # un.fq even though GGGCCC is present. Without prev, rmiss would match g2 and route to S instead.
    When I run `unmux reads.fq --group g1={s1=AACGTA} --group g1::loc=0:0:6 --group g2={s2=GGGCCC} --group g2::loc=0:6:12,prev=g1 --extract body=0:0:end --template body --sample S=g2::s2 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "S.fq" contains "AACGTAGGGCCCAAAA"
    And the file "un.fq" contains "TTTTTTGGGCCCAAAA"
    And the file "S.fq" does not contain "TTTTTTGGGCCCAAAA"

  Scenario: minFindsPerGroup requires a minimum number of matches before the group passes
    Given a file "reads.fq" containing:
      """
      @two
      AAAAAACCCCCCTTTTTTTT
      +
      IIIIIIIIIIIIIIIIIIII
      @one
      AAAAAATTTTTTTTTTTTTT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    # The read "two" contains both AAAAAA and CCCCCC (2 group matches); "one" contains only AAAAAA
    # (1 match). minFindsPerGroup=2 fails the group for "one", so its raw segment lands in un.fq while
    # "two" routes to S. Without the constraint (default min 0) a single match suffices and "one"
    # would route to S too.
    When I run `unmux reads.fq --group g={A=AAAAAA,B=CCCCCC,C=GGGGGG} --group g::minFindsPerGroup=2 --extract body=0:0:end --template body --sample S=g --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "S.fq" contains "AAAAAACCCCCCTTTTTTTT"
    And the file "un.fq" contains "AAAAAATTTTTTTTTTTTTT"
    And the file "S.fq" does not contain "AAAAAATTTTTTTTTTTTTT"

  Scenario: maxFindsPerGroup=1 fails a group whose match count exceeds the maximum
    Given a file "reads.fq" containing:
      """
      @two
      AAAAAACCCCCCTTTTTTTT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    # The read carries two group matches (AAAAAA and CCCCCC); maxFindsPerGroup=1 rejects a group with
    # more than one match, so the group fails and the raw read lands in un.fq. Removing the cap lets
    # both matches count and the read routes to S. (minFindsPerGroup is exercised separately above.)
    When I run `unmux reads.fq --group g={A=AAAAAA,B=CCCCCC} --group g::maxFindsPerGroup=1 --extract body=0:0:end --template body --sample S=g --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "un.fq" contains "AAAAAACCCCCCTTTTTTTT"
    And the file "S.fq" does not contain "AAAAAACCCCCCTTTTTTTT"

  Scenario: An AND selector across two groups routes by the combination, not the leftmost group alone
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAATTTTTTTTCCCC
      +
      IIIIIIIIIIIIIIIIII
      """
    # g1 matches t1 (AAAAAA) at [0,6); g2 matches p2 (TTTTTTTT) at [6,14). Two samples share the g1
    # selector and differ only in their g2 selector: SA needs g1::t1+g2::p1, SB needs g1::t1+g2::p2.
    # The g2 match (p2) is decisive, so the read routes to SB. Were g2 ignored both samples would tie
    # on g1::t1; the read assigns to SB because the AND tuple resolves to p2.
    When I run `unmux reads.fq --group g1={t1=AAAAAA,t2=CCCCCC} --group g1::loc=0:0:6 --group g2={p1=GGGGGGGG,p2=TTTTTTTT} --group g2::loc=0:6:14 --extract body=0:0:end --template body --sample SA=g1::t1+g2::p1 --sample SB=g1::t1+g2::p2 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "SB.fq" contains "AAAAAATTTTTTTTCCCC"
    And the file "SA.fq" does not contain "AAAAAATTTTTTTTCCCC"

  Scenario: three independent groups each pinned to one match route by the AND of all three
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAACCCCCCGGGGGG
      +
      IIIIIIIIIIIIIIIIII
      @r2
      AAAAAACCCCCCTTTTTT
      +
      IIIIIIIIIIIIIIIIII
      """
    # Three groups at adjacent windows, each with minFindsPerGroup=1,maxFindsPerGroup=1 (exactly one
    # match required). @r1 matches all three (grp_z = GGGGGG) and routes to S via the AND selector; @r2
    # carries TTTTTT where grp_z expects GGGGGG, so grp_z fails and the whole AND fails - @r2 lands in
    # un.fq. The grp_z term is decisive: dropping it would route @r2 to S as well.
    When I run `unmux reads.fq --group grp_x={x1=AAAAAA} --group grp_x::loc=0:0:6,minFindsPerGroup=1,maxFindsPerGroup=1 --group grp_y={y1=CCCCCC} --group grp_y::loc=0:6:12,minFindsPerGroup=1,maxFindsPerGroup=1 --group grp_z={z1=GGGGGG} --group grp_z::loc=0:12:18,minFindsPerGroup=1,maxFindsPerGroup=1 --extract body=0:0:end --template body --sample S=grp_x::x1+grp_y::y1+grp_z::z1 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "S.fq" contains "AAAAAACCCCCCGGGGGG"
    And the file "un.fq" contains "AAAAAACCCCCCTTTTTT"
    And the file "S.fq" does not contain "AAAAAACCCCCCTTTTTT"

  Scenario: a bare next=GROUP keeps the downstream group searching at its own loc
    Given a file "reads.fq" containing:
      """
      @r1
      AACGTAGGGCCCAAGGGCCCTTTTGGGG
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # g1 (AACGTA) matches at [0,6) and links to g2 with a BARE next=g2 (no :lo-hi window). A bare next is
    # only an ordering link, so g2 keeps its own loc=0:6:12 and matches the GGGCCC at [6,12); @g2+0:4 then
    # reads AAGG at [12,16). The relative-window form next=g2:8-20 (exercised in the first scenario of this
    # file) would instead search [14,26) and read TTTT, so the bare-vs-windowed distinction is load-bearing.
    When I run `unmux reads.fq --group g1={s1=AACGTA} --group g1::loc=0:0:6,next=g2 --group g2={s2=GGGCCC} --group g2::loc=0:6:12 --extract pay=@g2+0:4 --extract body=0:0:end --template body --tag PY=pay`
    Then the exit code is 0
    And stdout contains "PY:Z:AAGG"
    And stdout does not contain "PY:Z:TTTT"

  Scenario: maxFindsPerTag bounds a single tag's own match count (distinct from the per-group total)
    Given a file "reads.fq" containing:
      """
      @r
      ACGTACGTGGGGGGGGGGGGGGGGGGGGGGGGACGTACGT
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # The one tag ACGTACGT occurs twice (at [0,8) and [32,40)). maxFindsPerTag=1 caps each individual tag
    # at one match, so this read trips the PER-TAG bound and is unassigned; the QC slug names the
    # max_finds_per_tag constraint (not max_finds_per_group), which is what distinguishes the per-tag bound
    # from the per-group total exercised elsewhere in this file.
    When I run `unmux reads.fq --group g={t=ACGTACGT} --group g::dist=0,maxFindsPerTag=1 --sample S=g::t --out %sample.fq --unassigned un.fq --qc-tag`
    Then the exit code is 0
    And the file "un.fq" contains "max_finds_per_tag"
    And the file "un.fq" contains ""found":2"
    And the file "S.fq" does not contain "ACGTACGT"

  Scenario: maxFindsPerGroup counts the total across tags and fails the group when exceeded
    Given a file "reads.fq" containing:
      """
      @r
      AAAAAACCCCCCGGGGGG
      +
      IIIIIIIIIIIIIIIIII
      """
    # The read matches three different tags (a, b, c) once each, for three total group matches.
    # maxFindsPerGroup=1 counts the total across all tags and fails the group outright - it does NOT keep
    # the lowest-cost match and discard the rest. The QC slug records found=3 against the limit of 1, so
    # the read is unassigned. (Selecting a single best match across candidates is a mode=nearest behavior,
    # covered in matching.feature, not something maxFindsPerGroup does.)
    When I run `unmux reads.fq --group g={a=AAAAAA,b=CCCCCC,c=GGGGGG} --group g::dist=0,maxFindsPerGroup=1 --sample S=g --out %sample.fq --unassigned un.fq --qc-tag`
    Then the exit code is 0
    And the file "un.fq" contains "max_finds_per_group"
    And the file "un.fq" contains ""found":3"
    And the file "S.fq" does not contain "AAAAAACCCCCCGGGGGG"

  Scenario: a three-group chain assembles its CB tag by concatenating the per-group extracts
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAACCCCCCGGGGGG
      +
      IIIIIIIIIIIIIIIIII
      @r2
      AAAAAACCCCCCTTTTTT
      +
      IIIIIIIIIIIIIIIIII
      """
    # Three required groups, one per six-base window. @r1 matches all three, so each region is extracted
    # and CB joins them in order as AAAAAA-CCCCCC-GGGGGG with sep=-. @r2 fails the third group (TTTTTT is
    # not GGGGGG) and lands in un.fq, showing that any failed group sends the read to --unassigned rather
    # than emitting a partial CB.
    When I run `unmux reads.fq --group cbt={x=AAAAAA} --group cbt::loc=0:0:6,minFindsPerGroup=1,maxFindsPerGroup=1 --group cb2={y=CCCCCC} --group cb2::loc=0:6:12,minFindsPerGroup=1,maxFindsPerGroup=1 --group cb3={z=GGGGGG} --group cb3::loc=0:12:18,minFindsPerGroup=1,maxFindsPerGroup=1 --extract cbt=@cbt --extract cb2=@cb2 --extract cb3=@cb3 --extract body=0:0:end --template body --tag CB=cbt+cb2+cb3 --tag CB::sep=- --tag CB::qual=none --sample S=cbt::x+cb2::y+cb3::z --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "S.fq" contains "CB:Z:AAAAAA-CCCCCC-GGGGGG"
    And the file "un.fq" contains "AAAAAACCCCCCTTTTTT"
    And the file "S.fq" does not contain "AAAAAACCCCCCTTTTTT"

  Scenario: a variable-length first group resolves an anchored offset from the match end
    Given a file "reads.fq" containing:
      """
      @r6
      AGGGGGGCACACACACACACACACACACGTACGTA
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      @r7
      ATTTTTTTCACACACACACACACACACACGTACGTA
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # The group is 5'-anchored at loc 0:1:end with a 6 nt tag (GGGGGG) and a 7 nt tag (TTTTTTT). @r6 ends
    # its match at position 7 and @r7 at position 8, but @vlen+19:9 anchors the UMI 19 bp past each match
    # END, so both resolve to the same ACGTACGTA span. The shared UMI proves the offset tracks the matched
    # tag's length rather than a fixed coordinate.
    When I run `unmux reads.fq --group vlen={s1=GGGGGG,s2=TTTTTTT} --group vlen::loc=0:1:end,dist=0 --extract umi=@vlen+19:9 --extract body=0:0:6 --template body --tag RX=umi --tag RX::qual=none`
    Then the exit code is 0
    And stdout contains "@r6 RX:Z:ACGTACGTA"
    And stdout contains "@r7 RX:Z:ACGTACGTA"

  Scenario: a per-group maxFindsPerGroup fires ahead of a more generous per-tag maxFindsPerTag
    Given a file "reads.fq" containing:
      """
      @r
      AAAAGGAAAAGGCCCCGGCCCCGGCCCCGG
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # Tag a (AAAA) matches twice and tag b (CCCC) three times: no single tag exceeds maxFindsPerTag=3, but
    # the group total of five exceeds maxFindsPerGroup=1. The per-group cap is enforced on the total, so the
    # read fails on max_finds_per_group while the per-tag bound is never reached - the QC slug naming
    # max_finds_per_group (not max_finds_per_tag) is load-bearing on which cap dominates.
    When I run `unmux reads.fq --group g={a=AAAA,b=CCCC} --group g::dist=0,maxFindsPerTag=3,maxFindsPerGroup=1 --sample S=g --out %sample.fq --unassigned un.fq --qc-tag`
    Then the exit code is 0
    And the file "un.fq" contains "max_finds_per_group"
    And the file "S.fq" does not contain "AAAAGGAAAAGG"
