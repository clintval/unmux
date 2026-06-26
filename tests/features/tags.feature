Feature: SAM tag emission with quality tag pairing and separators
  Tag streams are assembled into SAM tags with paired quality tags per the default map or an
  override. A tag binding two or more streams joins them with `sep` (default `-`); when it also
  emits a quality tag the qualities join with `qual-sep` (default a single space). A single-stream
  tag needs neither. Streams concatenate with `+` (`CB=cbt+cb2+cb3`); a comma is an error. Default
  seq->qual map: CB/CY, RX/QX, BC/QT, OX/BZ. With no --out, tags render into the FASTQ read-name
  comment as `TAG:Z:VALUE`, so these scenarios assert on stdout.

  Scenario: CB tag from multiple streams with no quality tag
    Given a file "test.fq" containing:
      """
      @r1
      AAAATTTTCCCCGGGGACGTACGTTTTTAAAA
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux test.fq --extract cbt=0:0:8 --extract cb2=0:8:16 --extract cb3=0:16:24 --extract body=0:24:end --tag CB=cbt+cb2+cb3 --tag CB::sep=- --tag CB::qual=none --template body`
    Then the exit code is 0
    And stdout contains "CB:Z:AAAATTTT-CCCCGGGG-ACGTACGT"
    And stdout does not contain "CY:Z:"
    And stdout contains "TTTTAAAA"

  Scenario: RX tag with custom quality tag override
    Given a file "r1.fq" containing:
      """
      @r1
      ACGTACGTAGGGGGGGGGGG
      +
      IIIIIIIIIBBBBBBBBBBB
      """
    And a file "r2.fq" containing:
      """
      @r1
      TTTTCCCCAGGGGGGGGGGG
      +
      JJJJJJJJJCCCCCCCCCCC
      """
    When I run `unmux r1.fq r2.fq --extract umir1=0:0:9 --extract umir2=1:0:9 --extract body=0:9:end --tag RX=umir1+umir2 --tag RX::sep=- --tag RX::qual=BZ --tag RX::qual-sep=%20 --template body`
    Then the exit code is 0
    And stdout contains "RX:Z:ACGTACGTA-TTTTCCCCA"
    And stdout contains "BZ:Z:IIIIIIIII JJJJJJJJJ"
    And stdout does not contain "QX:Z:"

  Scenario: default quality tag pairing for single-stream tags
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGGTTTTACGTACGTAAAACCCC
      +
      0123456789ABCDEFGHIJKLMNOPQRSTUV
      """
    When I run `unmux test.fq --extract bc=0:0:8 --extract rx=0:8:16 --extract body=0:16:end --tag BC=bc --tag RX=rx --template body`
    Then the exit code is 0
    And stdout contains "BC:Z:AAAACCCC"
    And stdout contains "QT:Z:01234567"
    And stdout contains "RX:Z:GGGGTTTT"
    And stdout contains "QX:Z:89ABCDEF"

  Scenario: a multi-stream tag defaults sep to "-" and qual-sep to a space
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGGTTTTACGTACGTAAAACCCC
      +
      0123456789ABCDEFGHIJKLMNOPQRSTUV
      """
    When I run `unmux test.fq --extract bc1=0:0:8 --extract bc2=0:8:16 --extract body=0:16:end --tag BC=bc1+bc2 --template body`
    Then the exit code is 0
    And stdout contains "BC:Z:AAAACCCC-GGGGTTTT"
    And stdout contains "QT:Z:01234567 89ABCDEF"

  Scenario: a custom sequence separator joins multiple streams
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGGTTTTACGTACGTAAAACCCC
      +
      0123456789ABCDEFGHIJKLMNOPQRSTUV
      """
    When I run `unmux test.fq --extract bc1=0:0:8 --extract bc2=0:8:16 --extract body=0:16:end --tag BC=bc1+bc2 --tag BC::sep=: --tag BC::qual=none --template body`
    Then the exit code is 0
    And stdout contains "BC:Z:AAAACCCC:GGGGTTTT"
    And stdout does not contain "QT:Z:"

  Scenario: custom quality separator distinct from sequence separator
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGGTTTTACGTACGTAAAACCCC
      +
      0123456789ABCDEFGHIJKLMNOPQRSTUV
      """
    When I run `unmux test.fq --extract umi1=0:0:9 --extract umi2=0:9:18 --extract body=0:18:end --tag RX=umi1+umi2 --tag RX::sep=_ --tag RX::qual-sep=: --template body`
    Then the exit code is 0
    And stdout contains "RX:Z:AAAACCCCG_GGGTTTTAC"
    And stdout contains "QX:Z:012345678:9ABCDEFGH"

  Scenario: multi-segment tag layout with mixed quality assignments
    Given a file "r1.fq" containing:
      """
      @r1
      ACGTACGTAGGGGGGGGGGG
      +
      IIIIIIIIIBBBBBBBBBBB
      """
    And a file "r2.fq" containing:
      """
      @r1
      TTTTCCCCATGCAACGTATGC
      +
      JJJJJJJJJ01234567CCCC
      """
    # cbt (TGCAACGT) occurs once in r2 at [9,17) with distinct quals 01234567, so CB/BC pin that exact
    # match and CY/QT pin its exact quality slice (a repeated tag or uniform quals would pass even if
    # the coordinates were wrong).
    When I run `unmux r1.fq r2.fq --group grp_cbt={cbt01=TGCAACGT} --extract umir1=0:0:9 --extract umir2=1:0:9 --extract cbt=@grp_cbt --extract body=0:9:end --tag RX=umir1+umir2 --tag RX::sep=- --tag RX::qual=BZ --tag RX::qual-sep=%20 --tag CB=cbt --tag BC=cbt --template body`
    Then the exit code is 0
    And stdout contains "RX:Z:ACGTACGTA-TTTTCCCCA"
    And stdout contains "BZ:Z:IIIIIIIII JJJJJJJJJ"
    And stdout contains "CB:Z:TGCAACGT"
    And stdout contains "CY:Z:01234567"
    And stdout contains "BC:Z:TGCAACGT"
    And stdout contains "QT:Z:01234567"

  Scenario: --qc-tag writes a per-record demux-provenance slug for an assigned read
    Given a file "test.fq" containing:
      """
      @r1
      GAAGGGAAAA
      +
      IIIIIIIIII
      """
    # Observed GAAGGG is 1 mismatch from canonical GAAGAG; the slug records the sample, the matched
    # group, the matched read span (loc), the observed-vs-corrected bases, and the edit breakdown. Off
    # unless --qc-tag is given.
    When I run `unmux test.fq --group grp={dna01=GAAGAG} --group grp::loc=0:0,dist=1,minFindsPerGroup=1 --sample dna01=grp::dna01 --qc-tag`
    Then the exit code is 0
    And stdout contains "ZS:Z:{"v":1,"outcome":"assigned","sample":"dna01","sub_sample":null,"groups":[{"g":"grp","tag":"GAAGAG","loc":"0:0:6","obs":"GAAGGG","sub":1,"ind":0}]}"

  Scenario: --qc-tag is off by default
    Given a file "test.fq" containing:
      """
      @r1
      GAAGGGAAAA
      +
      IIIIIIIIII
      """
    When I run `unmux test.fq --group grp={dna01=GAAGAG} --group grp::loc=0:0,dist=1`
    Then the exit code is 0
    And stdout does not contain "ZS:Z:"

  Scenario: --qc-tag records the reason an unassigned read failed
    Given a file "test.fq" containing:
      """
      @r2
      TTTTTTAAAA
      +
      IIIIIIIIII
      """
    # The required group does not match, so the read is unassigned; its slug names the reason and the
    # offending group, written into the --unassigned bin.
    When I run `unmux test.fq --group grp={dna01=GAAGAG} --group grp::loc=0:0,dist=1,minFindsPerGroup=1 --sample dna01=grp::dna01 --unassigned un.fq --qc-tag`
    Then the exit code is 0
    And the file "un.fq" contains "ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"grp","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}"

  Scenario: --qc-tag rejects a SAM-reserved (two-uppercase) tag
    Given a file "test.fq" containing:
      """
      @r1
      ACGT
      +
      IIII
      """
    When I run `unmux test.fq --qc-tag=SC`
    Then the exit code is 1
    And stderr contains "local-use"
