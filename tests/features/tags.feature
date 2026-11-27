Feature: SAM tag emission with quality tag pairing and separators
  Tag streams are assembled into SAM tags with paired quality tags per the default map or an
  override. A tag binding two or more streams joins them with `sep` (default `-`); when it also
  emits a quality tag the qualities join with `qual-sep` (default a single space). A single-stream
  tag needs neither. Streams concatenate with `+` (`CB=cbt+cb2+cb3`); a comma is an error. Default
  seq->qual map: CB/CY, CR/CY, RX/QX, BC/QT, OX/BZ. A tag that two writers would set is an error,
  except that CB and CR may both default to CY, which then holds the CR qualities on a record where
  CR writes them and the CB qualities elsewhere. With no --out, tags render into the FASTQ
  read-name comment as `TAG:Z:VALUE`, so most of these scenarios assert on stdout.

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

  Scenario Outline: a tag that two writers would set is an error
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGG
      +
      IIII5555????
      """
    When I run `unmux test.fq --extract u=0:0:4 --extract v=0:4:8 <flags> --out out.sam`
    Then the exit code is 1
    And stderr contains "<message>"
    And a file "out.sam" does not exist

    Examples:
      | flags                                   | message                                                                                                               |
      | --tag ZS=u --qc-tag                     | tag `ZS` would be written twice, by `--tag ZS=u` and by `--qc-tag=ZS`; choose another tag for one of them              |
      | --tag XA=u --tag XA::qual=XA            | tag `XA` would be written twice, by `--tag XA=u` and by `--tag XA::qual=XA`; choose another tag for one of them        |
      | --tag XA=u --tag XB=v --tag XB::qual=XA | tag `XA` would be written twice, by `--tag XA=u` and by `--tag XB::qual=XA`; choose another tag for one of them        |
      | --tag RX=u --tag RX::qual=ZS --qc-tag   | tag `ZS` would be written twice, by `--tag RX::qual=ZS` and by `--qc-tag=ZS`; choose another tag for one of them       |
      | --tag RG=u                              | tag `RG` is reserved for the read group, which unmux sets itself on SAM/BAM/CRAM output; choose another tag for `--tag RG=u` |
      | --tag XA=u --tag XA::qual=RG            | tag `RG` is reserved for the read group, which unmux sets itself on SAM/BAM/CRAM output; choose another tag for `--tag XA::qual=RG` |
      | --tag CY=u --tag CB=v                   | tag `CY` would be written twice, by `--tag CY=u` and as the default quality tag of `--tag CB=v`; set `CB::qual=none`   |
      | --tag CB=u --tag CB::qual=CY --tag CR=v | tag `CY` would be written twice, by `--tag CB::qual=CY` and as the default quality tag of `--tag CR=v`; set `CR::qual=none` |
      | --tag CB=u --tag CR=v --tag CY=u        | tag `CY` would be written twice, as the default quality tag of `--tag CB=u` and by `--tag CY=u`; set `CB::qual=none`   |

  Scenario: a clash error gives a tag's streams with their fallbacks
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGG
      +
      IIII5555????
      """
    When I run `unmux test.fq --extract u=0:0:4 --extract v=0:4:8 --tag XA=v|u --tag XB=u --tag XB::qual=XA --out out.sam`
    Then the exit code is 1
    And stderr contains "tag `XA` would be written twice, by `--tag XA=v|u` and by `--tag XB::qual=XA`; choose another tag for one of them"
    And a file "out.sam" does not exist

  Scenario: RG is reserved even when every output is FASTX
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGG
      +
      IIII5555????
      """
    When I run `unmux test.fq --extract u=0:0:4 --tag RG=u --out out.fq`
    Then the exit code is 1
    And stderr contains "tag `RG` is reserved for the read group, which unmux sets itself on SAM/BAM/CRAM output; choose another tag for `--tag RG=u`"
    And a file "out.fq" does not exist

  Scenario Outline: CY holds the CR qualities when CB and CR both default to it, in either order
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGG
      +
      IIII5555????
      """
    When I run `unmux test.fq --pool p --extract u=0:0:4 --extract v=0:4:8 <flags> --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	AAAACCCCGGGG	IIII5555????	<tags>	RG:Z:p
      """

    Examples:
      | flags                                     | tags                            |
      | --tag CB=u --tag CR=v                     | CB:Z:AAAA	CY:Z:5555	CR:Z:CCCC |
      | --tag CR=v --tag CB=u                     | CR:Z:CCCC	CY:Z:5555	CB:Z:AAAA |
      | --tag CB=u                                | CB:Z:AAAA	CY:Z:IIII           |
      | --tag CB=u --tag CR=v --tag CR::qual=none | CB:Z:AAAA	CY:Z:IIII	CR:Z:CCCC |
      | --tag CB=u --tag CB::qual=none --tag CR=v | CB:Z:AAAA	CR:Z:CCCC	CY:Z:5555 |

  Scenario Outline: a corrected CB and a raw CR of one barcode keep the observed qualities in CY
    Given a file "test.fq" containing:
      """
      @r1
      GGGGGGGGGGCACGATAAAA
      +
      IIIIIIIIII56789:IIII
      """
    # CACGAT is one substitution from the reverse complement of AACGTG, so the corrected qualities
    # of CB run in reverse (:98765) while CR and its CY keep the observed order.
    When I run `unmux test.fq --pool p --group g={AACGTG} --group g::loc=0:10:16,bothStrands=true,dist=1 --extract bc=@g <flags> --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	GGGGGGGGGGCACGATAAAA	IIIIIIIIII56789:IIII	<tags>	RG:Z:p
      """

    Examples:
      | flags                                      | tags                                  |
      | --tag CB=bc --tag CR=bc --tag CR::raw=true | CB:Z:AACGTG	CY:Z:56789:	CR:Z:CACGAT |
      | --tag CR=bc --tag CR::raw=true --tag CB=bc | CR:Z:CACGAT	CY:Z:56789:	CB:Z:AACGTG |
      | --tag CB=bc                                | CB:Z:AACGTG	CY:Z::98765             |

  Scenario Outline: CY holds the CB qualities on a record where CR is left out
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGG
      +
      IIII5555????
      @r2
      AAAATTTTGGGG
      +
      IIII5555????
      """
    When I run `unmux test.fq --pool p --group g={CCCC} --group g::loc=0:4:8 --extract u=0:8:12 --extract v=@g <before> --tag CR=v <after> --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	AAAACCCCGGGG	IIII5555????	<r1>	RG:Z:p
      r2	4	*	0	255	*	*	0	0	AAAATTTTGGGG	IIII5555????	CB:Z:GGGG	CY:Z:????	RG:Z:p
      """

    Examples:
      | before     | after      | r1                              |
      | --tag CB=u |            | CB:Z:GGGG	CY:Z:5555	CR:Z:CCCC |
      |            | --tag CB=u | CR:Z:CCCC	CY:Z:5555	CB:Z:GGGG |

  Scenario Outline: CY holds the CB qualities on a record where CR falls back to text
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGG
      +
      IIII5555????
      @r2
      AAAATTTTGGGG
      +
      IIII5555????
      """
    # On r2 the stream v is absent, so CR takes the a1 tag ID, which has no qualities.
    When I run `unmux test.fq --pool p --group a={a1=AAAA} --group a::loc=0:0:4 --group g={CCCC} --group g::loc=0:4:8 --extract u=0:8:12 --extract v=@g --extract n=@a::id <before> --tag CR=v|n <after> --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	AAAACCCCGGGG	IIII5555????	<r1>	RG:Z:p
      r2	4	*	0	255	*	*	0	0	AAAATTTTGGGG	IIII5555????	<r2>	RG:Z:p
      """

    Examples:
      | before     | after      | r1                              | r2                            |
      | --tag CB=u |            | CB:Z:GGGG	CY:Z:5555	CR:Z:CCCC | CB:Z:GGGG	CY:Z:????	CR:Z:a1 |
      |            | --tag CB=u | CR:Z:CCCC	CY:Z:5555	CB:Z:GGGG | CR:Z:a1	CB:Z:GGGG	CY:Z:???? |

  Scenario: a record where CR is left out gets the CB qualities in CY, not a carried CY
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AACGTGCAGGGG	ABCDEFGHIJKL	CB:Z:OLDCB	CR:Z:OLDCR	CY:Z:OLDCY
      r2	4	*	0	0	*	*	0	0	TTTTTTTTGGGG	ABCDEFGHIJKL	CB:Z:OLDCB	CR:Z:OLDCR	CY:Z:OLDCY
      """
    When I run `unmux in.sam --pool p --group bc={x=AACGTGCA} --group bc::loc=0:0:8 --extract b=@bc --extract u=0:0:4 --tag CB=u --tag CR=b --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	AACGTGCAGGGG	ABCDEFGHIJKL	CB:Z:AACG	CR:Z:AACGTGCA	CY:Z:ABCDEFGH	RG:Z:p
      r2	4	*	0	255	*	*	0	0	TTTTTTTTGGGG	ABCDEFGHIJKL	CB:Z:TTTT	CY:Z:ABCD	RG:Z:p
      """

  Scenario: a carried CY goes when CB is rewritten without one and CR is left out
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AACGTGCAGGGG	ABCDEFGHIJKL	CB:Z:OLDCB	CR:Z:OLDCR	CY:Z:OLDCY
      r2	4	*	0	0	*	*	0	0	TTTTTTTTGGGG	ABCDEFGHIJKL	CB:Z:OLDCB	CR:Z:OLDCR	CY:Z:OLDCY
      """
    When I run `unmux in.sam --pool p --group bc={x=AACGTGCA} --group bc::loc=0:0:8 --extract b=@bc --extract u=0:0:4 --tag CB=u --tag CB::qual=none --tag CR=b --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	AACGTGCAGGGG	ABCDEFGHIJKL	CB:Z:AACG	CR:Z:AACGTGCA	CY:Z:ABCDEFGH	RG:Z:p
      r2	4	*	0	255	*	*	0	0	TTTTTTTTGGGG	ABCDEFGHIJKL	CB:Z:TTTT	RG:Z:p
      """

  Scenario: a CB qual-sep is warned about where CR writes the CY it would join
    Given a file "test.fq" containing:
      """
      @r1
      AAAACCCCGGGG
      +
      IIII5555????
      """
    When I run `unmux test.fq --pool p --extract u=0:0:4 --extract v=0:4:8 --tag CB=u+v --tag CB::qual-sep=_ --tag CR=u+v --out out.sam`
    Then the exit code is 0
    And stderr contains "tag `CB` sets `qual-sep`, which applies only to records where `CR` writes no `CY`, since `CY` holds the `CR` qualities on the others"
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	AAAACCCCGGGG	IIII5555????	CB:Z:AAAA-CCCC	CY:Z:IIII 5555	CR:Z:AAAA-CCCC	RG:Z:p
      """
