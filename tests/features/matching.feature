Feature: Tag matching with distance tolerance and error handling
  Tag matching behavior for `unmux`, exercised over the real binary on small fixtures.
  Covers exact match; Hamming distance tolerance; substitution/indel distinctions; the
  three-part distance format; nearest-neighbor disambiguation via delta; IUPAC tag wildcards;
  search-window locations; variable-length anchor resolution; and error correction.

  Most scenarios match a single-tag group, extract the matched span with `@group`, and surface
  it as a SAM tag in the FASTQ read-name comment (`TAG:Z:VALUE`), so the matched/extracted bases
  are directly observable on stdout. A read that fails to match produces no output record.

  Scenario: Exact tag match at location with zero tolerance
    Given a file "reads.fq" containing:
      """
      @r1
      TGGGTGTTAACCGGTT
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group mybc={s1=TGGGTGTT} --group mybc::loc=0:0:8,dist=0 --extract bc=@mybc --tag CB=bc --template bc`
    Then the exit code is 0
    And stdout contains "CB:Z:TGGGTGTT"

  Scenario: Hamming mismatch is rejected when dist=0 and the read routes to unassigned
    Given a file "reads.fq" containing:
      """
      @hit
      TGGGTGTTAACCGGTT
      +
      IIIIIIIIIIIIIIII
      @miss
      TGGGTGATAACCGGTT
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group mybc={s1=TGGGTGTT} --group mybc::loc=0:0:8,dist=0 --extract rest=0:0:end --template rest --sample S1=mybc::s1 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And a file "S1.fq" exists
    And the file "S1.fq" contains "TGGGTGTTAACCGGTT"
    And the file "un.fq" contains "TGGGTGATAACCGGTT"
    And the file "S1.fq" does not contain "TGGGTGATAACCGGTT"

  Scenario: Single Hamming mismatch accepted with dist=1
    Given a file "reads.fq" containing:
      """
      @exact
      ATGAGACCCCDDDD
      +
      IIIIIIIIIIIIII
      @onesub
      TTGAGACCCCDDDD
      +
      IIIIIIIIIIIIII
      @far
      GGGGGGCCCCDDDD
      +
      IIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group bc3x={s1=ATGAGA} --group bc3x::loc=0:0:6,dist=1 --extract bc=@bc3x --tag BC=bc --tag BC::raw=true --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:ATGAGA"
    And stdout contains "BC:Z:TTGAGA"
    And stdout does not contain "BC:Z:GGGGGG"

  Scenario: Two substitutions accepted at dist=2 but rejected at dist=1
    Given a file "reads.fq" containing:
      """
      @twosub
      GGTAGGCCCCDDDD
      +
      IIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group bc6={s1=GGGGGG} --group bc6::loc=0:0:6,dist=2 --extract bc=@bc6 --tag BC=bc --tag BC::raw=true --tag CB=bc --tag CB::qual=none --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:GGTAGG"
    And stdout contains "CB:Z:GGGGGG"

  Scenario: Two substitutions exceed dist=1 and the read does not match
    Given a file "reads.fq" containing:
      """
      @twosub
      GGTAGGCCCCDDDD
      +
      IIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group bc6={s1=GGGGGG} --group bc6::loc=0:0:6,dist=1 --extract bc=@bc6 --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout does not contain "BC:Z:"

  Scenario: A pure-indel read matches once a 1-indel budget is set (dist=0:1)
    Given a file "reads.fq" containing:
      """
      @exact
      ATGACACCCCDDDD
      +
      IIIIIIIIIIIIII
      @del
      ATACACCCCDDDD
      +
      IIIIIIIIIIIII
      @ins
      ATGGACACCCCDDDD
      +
      IIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group b={s1=ATGACA} --group b::loc=0:0:8,dist=0:1 --extract bc=@b --tag BC=bc --tag BC::raw=true --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:ATGACA"
    And stdout contains "BC:Z:ATACA"
    And stdout contains "BC:Z:ATGGACA"

  Scenario: A bare dist=2 counts substitutions only and rejects a single deletion
    Given a file "reads.fq" containing:
      """
      @del
      ATACACCCCDDDD
      +
      IIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group b={s1=ATGACA} --group b::loc=0:0:8,dist=2 --extract bc=@b --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout does not contain "BC:Z:"

  Scenario: The per-type substitution cap gates independently of the total budget
    Given a file "reads.fq" containing:
      """
      @twosub
      GGTAGGCCCCDDDD
      +
      IIIIIIIIIIIIII
      """
    # The 2-substitution read GGTAGG matches the tag GGGGGG at dist=2 (sub cap 2; see the scenario
    # above). Here dist=1:1 sets subs=1, indel=1, and total defaults to their sum (2) - the same total
    # budget as dist=2 - but the per-type sub cap of 1 rejects the same read: the cap is a separate
    # knob from the total.
    When I run `unmux reads.fq --group bc6={s1=GGGGGG} --group bc6::loc=0:0:6,dist=1:1 --extract bc=@bc6 --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout does not contain "BC:Z:"

  Scenario: Nearest mode returns the single best match when one tag clearly wins
    Given a file "reads.fq" containing:
      """
      @r
      TGAGCCCCCDDDD
      +
      IIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group g={s1=TGAGC,s2=TGAGA} --group g::loc=0:0:5,dist=1,mode=nearest,delta=1 --extract bc=@g --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:TGAGC"

  Scenario: Nearest mode leaves an equidistant tie unassigned
    Given a file "reads.fq" containing:
      """
      @r
      AATTCCCCDDDD
      +
      IIIIIIIIIIII
      """
    # AATT is 2 mismatches from both AAAA and TTTT; mode=nearest cannot break the tie within delta, so
    # the read is unassigned. (The IUPAC guard warns that this barcode pair overlaps at dist=2, but
    # under mode=nearest it is only a warning, not a fail-fast.)
    When I run `unmux reads.fq --group g={s1=AAAA,s2=TTTT} --group g::loc=0:0:4,dist=2,mode=nearest,delta=1 --extract bc=@g --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout does not contain "BC:Z:AATT"

  Scenario: An ambiguous barcode pair under mode=all is rejected by the collision guard
    Given a file "reads.fq" containing:
      """
      @r
      AATTCCCCDDDD
      +
      IIIIIIIIIIII
      """
    # The same pair under mode=all would silently assign the equidistant read AATT to one of them, so
    # the IUPAC collision guard rejects the config up front rather than mis-assigning (the mode=nearest
    # scenario above instead warns and drops the tie).
    When I run `unmux reads.fq --group g={s1=AAAA,s2=TTTT} --group g::loc=0:0:4,dist=2,mode=all --extract bc=@g --tag BC=bc --template bc`
    Then the exit code is 1
    And stderr contains "within dist=2"
    And stderr contains "mode=all"

  Scenario: A larger delta requirement leaves a one-gap winner unassigned
    Given a file "reads.fq" containing:
      """
      @r
      TGAGCCCCCDDDD
      +
      IIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group g={s1=TGAGC,s2=TGAGA} --group g::loc=0:0:5,dist=1,mode=nearest,delta=2 --extract bc=@g --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout does not contain "BC:Z:"

  Scenario: IUPAC ambiguity in the tag definition matches any base at the N positions
    Given a file "reads.fq" containing:
      """
      @aa
      TGGGTGAACCCC
      +
      IIIIIIIIIIII
      @cc
      TGGGTGCCCCCC
      +
      IIIIIIIIIIII
      @ta
      TGGGTGTACCCC
      +
      IIIIIIIIIIII
      """
    When I run `unmux reads.fq --group iu={s1=TGGGTGNN} --group iu::loc=0:0:8,dist=0 --extract bc=@iu --tag BC=bc --tag BC::raw=true --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:TGGGTGAA"
    And stdout contains "BC:Z:TGGGTGCC"
    And stdout contains "BC:Z:TGGGTGTA"

  Scenario: A from-end window selects the 3-prime occurrence so the anchored extract changes
    Given a file "reads.fq" containing:
      """
      @r
      AACGAACAGAGTTTTCCCCCCCCCCAACGAACAGAGGGGG
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # The tag AACGAACAGAG occurs twice: at the read start (followed by TTTT) and in the last 15 bases
    # (followed by GGGG). A whole-read search would take the leftmost (start) occurrence and extract
    # TTTT; the from-end window loc=0:-15:end excludes the start, so the match anchors on the 3-prime
    # occurrence and the +0:4 extract is GGGG instead.
    When I run `unmux reads.fq --group cbt={s1=AACGAACAGAG} --group cbt::loc=0:-15:end,dist=0 --extract pay=@cbt+0:4 --tag PY=pay --template pay`
    Then the exit code is 0
    And stdout contains "PY:Z:GGGG"
    And stdout does not contain "PY:Z:TTTT"

  Scenario: Variable-length barcodes resolve the anchor end so an offset extraction lines up
    Given a file "reads.fq" containing:
      """
      @r6
      AGGGGGGCCCCCCCCCCCCCCCCCCCACGTACGTATTTTTTTTTT
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      @r7
      ATTTTTTTCCCCCCCCCCCCCCCCCCCACGTACGTAGGGGGGGGGG
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group vlen={s1=GGGGGG,s2=TTTTTTT} --group vlen::loc=0:1:end,dist=0 --extract pay=@vlen+19:9 --tag PY=pay --template pay`
    Then the exit code is 0
    And stdout contains "PY:Z:ACGTACGTA"

  Scenario: anchor=5p matches a variable-length 5'-anchored group at each tag's own length
    Given a file "reads.fq" containing:
      """
      @r
      AGTACTCTGGGG
      +
      IIIIIIIIIIII
      """
    # Two 5'-anchored tags of different lengths. Under an unanchored window the 6nt s2 (TACTCA) would slide to
    # offset 2 (TACTCT, 1 mismatch) for a spurious 2nd hit; anchor=5p pins each tag at loc.start
    # over its own length, so only the 7nt s1 (AGTACTC) matches and just one tag is assigned.
    When I run `unmux reads.fq --group grp={s1=AGTACTC,s2=TACTCA} --group grp::loc=0:0:8,dist=1,maxFindsPerGroup=1,anchor=5p --extract bc=@grp --tag BC=bc --tag BC::raw=true --template bc::raw=true`
    Then the exit code is 0
    And stdout contains "BC:Z:AGTACTC"
    And stdout does not contain "BC:Z:TACTCA"

  Scenario: A leftward anchor offset extracts a UMI upstream of the matched barcode
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTGGGGGGTTTT
      +
      IIIIIIIIIIIIIIIIII
      """
    # The barcode GGGGGG matches at [8,14). @bc-0:8 reads the 8 bp ending at the match START (the UMI
    # to its left), while @bc+0:end takes the tail to its right. The minus direction is load-bearing:
    # +0:8 would read rightward into the tail instead.
    When I run `unmux reads.fq --group bc={s1=GGGGGG} --group bc::loc=0:8:14 --extract umi=@bc-0:8 --extract body=@bc+0:end --tag RX=umi --template body`
    Then the exit code is 0
    And stdout contains "RX:Z:ACGTACGT"
    And stdout contains "TTTT"

  Scenario: error correction is the default and stores the canonical tag sequence
    Given a file "reads.fq" containing:
      """
      @mm
      TTGACACCCCDDDD
      +
      IIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group ec={s1=ATGACA} --group ec::loc=0:0:8,dist=1 --extract bc=@ec --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:ATGACA"
    And stdout does not contain "BC:Z:TTGACA"

  Scenario: raw=true keeps the observed bases for a one-mismatch read
    Given a file "reads.fq" containing:
      """
      @mm
      TTGACACCCCDDDD
      +
      IIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group wc={s1=ATGACA} --group wc::loc=0:0:8,dist=1 --extract bc=@wc --tag BC=bc --tag BC::raw=true --template bc::raw=true`
    Then the exit code is 0
    And stdout contains "BC:Z:TTGACA"
    And stdout does not contain "BC:Z:ATGACA"

  Scenario: The split sub/indel budget dist=1:0 accepts a substitution but rejects an indel
    Given a file "reads.fq" containing:
      """
      @exact
      TGAGCAAAA
      +
      IIIIIIIII
      @sub
      TGCGCAAAA
      +
      IIIIIIIII
      @del
      TGACAAAA
      +
      IIIIIIII
      """
    # dist=1:0 is subs=1, indels=0. TGCGC is the tag TGAGC with one substitution (accepted); TGAC is
    # the tag with one deletion, which can only align via an indel (forbidden), so it is rejected and
    # produces no record. Under dist=1:1 (one indel allowed) the deletion read would also match, so the
    # indels=0 budget is load-bearing on dropping @del.
    When I run `unmux reads.fq --group bc={s1=TGAGC} --group bc::loc=0:0:6,dist=1:0 --extract x=@bc --tag BC=x --tag BC::raw=true --template x`
    Then the exit code is 0
    And stdout contains "BC:Z:TGAGC"
    And stdout contains "BC:Z:TGCGC"
    And stdout does not contain "BC:Z:TGAC"
    And stdout does not contain "@del"

  Scenario: A three-way equidistant tie is left unassigned under mode=nearest
    Given a file "reads.fq" containing:
      """
      @r
      ACGTTGGGG
      +
      IIIIIIIII
      """
    # ACGTT is one substitution from each of ACGTA, ACGTC, and ACGTG, so all three tags tie at cost 1.
    # mode=nearest needs the best to beat the runner-up by delta=1, but every candidate is equal-best, so
    # no delta gap exists and the read is ambiguous and dropped. (mode=nearest warns about the overlap
    # rather than failing fast, as it does under mode=all.)
    When I run `unmux reads.fq --group g={s1=ACGTA,s2=ACGTC,s3=ACGTG} --group g::loc=0:0:5,dist=1,mode=nearest,delta=1 --extract bc=@g --tag BC=bc --tag BC::raw=true --template bc`
    Then the exit code is 0
    And stdout does not contain "BC:Z:"

  Scenario: A minimal delta=1 admits a one-gap best-cost winner under mode=nearest
    Given a file "reads.fq" containing:
      """
      @r
      ATGAGAGGGG
      +
      IIIIIIIIII
      """
    # ATGAGA matches the tag ATGAGA at cost 0 and the tag ATGAGG at cost 1, a one-cost gap. delta=1 is
    # exactly met, so the best-cost tag wins and the read assigns to ATGAGA. The companion scenario above
    # with delta=2 leaves the same one-gap winner unassigned, so the delta threshold is load-bearing.
    When I run `unmux reads.fq --group g={s1=ATGAGA,s2=ATGAGG} --group g::loc=0:0:6,dist=2,mode=nearest,delta=1 --extract bc=@g --tag BC=bc --tag BC::raw=true --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:ATGAGA"

  Scenario: a lowercase barcode still matches an uppercase read (case-insensitive matching)
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTTTTT
      +
      IIIIIIIIIIII
      """
    # A barcode entered in lowercase (a common slip) still matches the uppercase read at dist=0:
    # sassy's Iupac profile case-folds, so no manual uppercasing of the tag is needed. The extracted
    # bases are the read's own (uppercase) bytes.
    When I run `unmux reads.fq --group bc={s1=acgtacgt} --group bc::loc=0:0:8,dist=0 --extract bc=@bc --extract body=0:8:end --tag BC=bc --tag BC::raw=true --template body`
    Then the exit code is 0
    And stdout contains "BC:Z:ACGTACGT"

  Scenario: anchor=3p pins each tag's 3' base at the loc window end
    Given a file "reads.fq" containing:
      """
      @r1
      TTTTTACGTA
      +
      IIIIIIIIII
      """
    # The barcode ACGTA sits flush with the 3' end of the 0:0:10 window. anchor=3p pins each
    # tag's 3' base at loc.end, so it matches the trailing ACGTA at [5,10) and routes the read to
    # dna01; a anchor=5p group would instead test [0,5)=TTTTT and find nothing. anchor3p is the
    # 3' mirror of anchor5p.
    When I run `unmux reads.fq --group bc={s1=ACGTA} --group bc::loc=0:0:10,anchor=3p,dist=0 --extract bc=@bc --extract body=0:0:end --tag BC=bc --tag BC::raw=true --sample dna01=bc::s1 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "dna01.fq" contains "BC:Z:ACGTA"
