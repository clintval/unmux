Feature: Extraction spans and preserved qualities
  Extraction carves named streams (bases and their original per-base qualities) out of the input
  reads for `unmux`, exercised over the real binary on small fixtures. Covers fixed start-anchored
  spans, end-relative and negative-end spans, reverse-complement, between-anchor regions, anchored
  offsets, and multi-stream tag assembly with separators. sassy never touches qualities: an extracted
  stream slices the original quals by the matched coordinates, so no separate quality-repair step is needed.

  Each scenario surfaces the extracted bases (and, where relevant, the sliced qualities) as SAM tags
  in the FASTQ read-name comment (`TAG:Z:VALUE` plus a paired quality tag), so they are directly
  observable on stdout. Every fixture uses non-repeating bases and per-position-distinct quality
  characters so that a wrong slice position, a missing reverse, or a leaked separator would change the
  asserted value.

  Scenario: A fixed start-anchored span carries the original qualities of those positions
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTAGTGCATGCAT
      +
      012345678abcdefghij
      """
    When I run `unmux reads.fq --extract umi=0:0:9 --extract body=0:9:end --tag RX=umi --tag RX::qual=QX --template body`
    Then the exit code is 0
    And stdout contains "RX:Z:ACGTACGTA"
    And stdout contains "QX:Z:012345678"
    # Pin the span END: a one-base overrun (umi=0:0:10) would emit these grown values, which the
    # prefix `contains` checks above would still accept.
    And stdout does not contain "RX:Z:ACGTACGTAG"
    And stdout does not contain "QX:Z:012345678a"

  Scenario: An end-relative span resolves the tail and slices its qualities
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTAGTGCATGCATG
      +
      IIIIIIIIIIabcdefghij
      """
    # 0:-10:end takes the last 10 bases; their qualities are abcdefghij. A whole-read span (0:0:end)
    # would instead carry the full IIIIIIIIIIabcdefghij, so the from-end coordinate is load-bearing.
    When I run `unmux reads.fq --extract tail=0:-10:end --extract body=0:0:10 --tag TL=tail --tag TL::qual=TQ --template body`
    Then the exit code is 0
    And stdout contains "TL:Z:TGCATGCATG"
    And stdout contains "TQ:Z:abcdefghij"
    And stdout does not contain "TQ:Z:IIIIIIIIIIabcdefghij"

  Scenario: A negative end gives a fixed 5' / variable middle / fixed 3' layout
    Given a file "reads.fq" containing:
      """
      @short
      AACCGGATCGATCGAT
      +
      IIIIIIIIIIIIIIII
      @long
      AACCGTTTAGCATTGCAGTCAT
      +
      IIIIIIIIIIIIIIIIIIIIII
      """
    # umi=0:0:5 is the fixed 5'; insert=0:5:-2 is the variable middle (stops 2 bp before each read's
    # own end); adapter3=0:-2:end is the fixed 3'. The middle length differs per read (9 vs 15 bp), so
    # the -2 end is load-bearing: with a plain 0:5:end the body would swallow the AT adapter (body
    # GATCGATCGAT) instead of stopping before it. This reproduces an enhanced 5M...2S read structure.
    When I run `unmux reads.fq --extract umi=0:0:5 --extract insert=0:5:-2 --extract adapter3=0:-2:end --template insert --tag A3=adapter3`
    Then the exit code is 0
    And stdout contains "GATCGATCG"
    And stdout does not contain "GATCGATCGAT"
    And stdout contains "TTTAGCATTGCAGTC"
    And stdout contains "A3:Z:AT"

  Scenario: A ~ on a tag stream reverse-complements its bases and reverses its qualities
    Given a file "reads.fq" containing:
      """
      @r1
      AAACCCGTGTGCATGCATGC
      +
      01234567abcdefghijkl
      """
    # The extract pulls the forward AAACCCGT / 01234567; the `~` on the tag stream (RX=~umi)
    # reverse-complements that contribution to ACGGGTTT and reverses its qualities to 76543210.
    # Dropping the ~ would emit the forward AAACCCGT / 01234567, so the ~ decides both the bases and
    # the quality order. The sequence is non-palindromic, so the reverse complement is distinguishable.
    When I run `unmux reads.fq --extract umi=0:0:8 --extract body=0:8:end --tag RX=~umi --tag RX::qual=QX --template body`
    Then the exit code is 0
    And stdout contains "RX:Z:ACGGGTTT"
    And stdout contains "QX:Z:76543210"
    And stdout does not contain "RX:Z:AAACCCGT"
    And stdout does not contain "QX:Z:01234567"

  Scenario: A between-anchor span captures only the sequence intervening two matched groups
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAACGCGCGTACGTTTTTTGGGG
      +
      IIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # gA (AAAAAA) ends at position 6, gB (TTTTTT) starts at position 16, so @gA..@gB is bases [6,16) =
    # CGCGCGTACG. The second anchor is load-bearing: a single-anchor @gA+0:end would over-capture
    # CGCGCGTACGTTTTTTGGGG; the `..` form stops exactly at gB's matched start.
    When I run `unmux reads.fq --group gA={s1=AAAAAA} --group gB={s1=TTTTTT} --extract spacer=@gA..@gB --extract body=0:0:6 --tag SP=spacer --template body`
    Then the exit code is 0
    And stdout contains "SP:Z:CGCGCGTACG"
    And stdout does not contain "SP:Z:CGCGCGTACGTTTTTT"

  Scenario: A ~ tag stream reverse-complements an anchored-offset span
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGGGGGGTTTTTTTTTTATCGATCGATCGCC
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # The barcode GGGGGG ends at position 12; @bc+10:12 takes the 12 bases starting 10 past the match
    # end ([22,34) = ATCGATCGATCG). The `~` on the tag stream (PY=~pay) reverse-complements that to
    # CGATCGATCGAT. Dropping the ~ would leave the forward ATCGATCGATCG, so the ~ is decisive over an
    # anchored, non-palindromic span.
    When I run `unmux reads.fq --group bc={s1=GGGGGG} --extract pay=@bc+10:12 --extract body=0:0:6 --tag PY=~pay --template body`
    Then the exit code is 0
    And stdout contains "PY:Z:CGATCGATCGAT"
    And stdout does not contain "PY:Z:ATCGATCGATCG"

  Scenario: A multi-stream tag joins its components with the configured separator
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAAAACCCCCCCCGGGGGGGGTTTT
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # CB binds three 8 bp streams in order, joined with sep=-; the separator is load-bearing (sep=_
    # would render CB:Z:AAAAAAAA_CCCCCCCC_GGGGGGGG instead).
    When I run `unmux reads.fq --extract s1=0:0:8 --extract s2=0:8:16 --extract s3=0:16:24 --extract body=0:24:end --tag CB=s1+s2+s3 --tag CB::sep=- --tag CB::qual=none --template body`
    Then the exit code is 0
    And stdout contains "CB:Z:AAAAAAAA-CCCCCCCC-GGGGGGGG"
    And stdout does not contain "CB:Z:AAAAAAAA_CCCCCCCC_GGGGGGGG"

  Scenario: A multi-stream quality tag joins component qualities with its own quality separator
    Given a file "reads.fq" containing:
      """
      @r1
      AAAACCCCGGGGGGGG
      +
      0123abcdIIIIIIII
      """
    # RX binds two streams: the sequence joins with sep=- and the qualities join with qual-sep=_, so
    # QX:Z:0123_abcd. The qual-sep is load-bearing and a separate knob from sep: qual-sep=%20 would
    # render QX:Z:0123 abcd (a space), and a multi-stream quality tag with no qual-sep is an error.
    When I run `unmux reads.fq --extract u1=0:0:4 --extract u2=0:4:8 --extract body=0:8:end --tag RX=u1+u2 --tag RX::sep=- --tag RX::qual=QX --tag RX::qual-sep=_ --template body`
    Then the exit code is 0
    And stdout contains "RX:Z:AAAA-CCCC"
    And stdout contains "QX:Z:0123_abcd"
    And stdout does not contain "QX:Z:0123 abcd"

  Scenario: Independent spans keep their own qualities with no cross-contamination
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAAAAAACCCCCCGGGGGGGGGG
      +
      0123456789ABCDEFqrstuvwxyz
      """
    # Three disjoint regions feed three channels: umi [0:10) -> UY, barcode [10:16) -> QT, template
    # [16:end) -> QUAL. Each quality tag slices only its own region (0123456789 / ABCDEF / qrstuvwxyz);
    # the barcode span [10:16) is load-bearing for QT, since a [0:6) span would make QT read 012345.
    When I run `unmux reads.fq --extract umi=0:0:10 --extract barcode=0:10:16 --extract template=0:16:end --tag UZ=umi --tag UZ::qual=UY --tag BC=barcode --tag BC::qual=QT --template template`
    Then the exit code is 0
    And stdout contains "UY:Z:0123456789"
    And stdout contains "QT:Z:ABCDEF"
    And stdout does not contain "QT:Z:012345"

  Scenario: qual=none suppresses the quality tag the default map would emit
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAAAACCCCCCCCGGGG
      +
      abcdefghIIIIIIIIIIII
      """
    # A single-stream CB tag defaults to emitting its paired quality tag CY (the default seq->qual
    # map). qual=none is load-bearing: without it the record would carry CY:Z:abcdefgh.
    When I run `unmux reads.fq --extract s1=0:0:8 --extract body=0:8:end --tag CB=s1 --tag CB::qual=none --template body`
    Then the exit code is 0
    And stdout contains "CB:Z:AAAAAAAA"
    And stdout does not contain "CY:Z:"

  Scenario: An anchored offset on a second-segment barcode resolves across the I2 read
    Given a file "r1.fq" containing:
      """
      @q
      GGGGGGGGGGGGGGGG
      +
      IIIIIIIIIIIIIIII
      """
    And a file "i2.fq" containing:
      """
      @q
      AACGTACGCACACACACACACACACACTTTTGGGGC
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    And a file "cbt.tags.tsv" containing:
      """
      id	seq
      cbt7	ACGTACG
      """
    # The cbt group is sourced from a TSV and pinned to segment 1 (the I2 read) at loc 1:1:end. The 7 nt
    # ACGTACG matches I2 at [1,8); @grp_cbt captures it (CB:Z:ACGTACG) and @grp_cbt+19:9 reads the UMI 19 bp
    # past the match end at I2[27,36) = TTTTGGGGC (RX). The extract coordinates are resolved within segment
    # 1 even though the matched read is the second input, so the cross-segment anchor is load-bearing.
    When I run `unmux --in 0=r1.fq --in 1=i2.fq --group grp_cbt=cbt.tags.tsv --group grp_cbt::loc=1:1:end,dist=0 --extract cbt=@grp_cbt --extract umi=@grp_cbt+19:9 --extract body=0:0:end --template body --tag CB=cbt --tag CB::qual=none --tag RX=umi --tag RX::qual=none`
    Then the exit code is 0
    And stdout contains "CB:Z:ACGTACG"
    And stdout contains "RX:Z:TTTTGGGGC"
