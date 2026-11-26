Feature: Input tags on --template records
  A record built by --template keeps the input tags of the read it is cut from, as a record passed
  through whole does, so UMI and barcode tags, read groups, FASTX comment tags, and per-read vendor
  tags survive trimming. Per-base tags describe the input bases one by one (MM/ML base
  modifications, OQ original qualities, PacBio kinetics, the ONT move table), so a record keeps
  them only when it is its read unchanged, with the same bases in the same orientation. Otherwise
  they are left out, and one warning per run lists them. A record that joins two reads takes the
  tags of the read that supplies its first segment, without the per-base ones. A --tag binding
  still wins over a carried tag.

  Scenario Outline: a template that copies the whole read keeps every input tag
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGTCG	IIIIIIIIII	MM:Z:C+m?,0,0;	ML:B:C,200,10	qs:f:20	ch:i:7	np:i:12	rq:f:0.999	zm:i:42
      """
    When I run `unmux in.sam --out tagged.<ext>`
    Then the exit code is 0
    When I run `unmux tagged.<ext> --pool p --extract all=0:0:end --template all --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGTCG	IIIIIIIIII	MM:Z:C+m?,0,0;	ML:B:C,200,10	qs:f:20	ch:i:7	np:i:12	rq:f:0.999	zm:i:42	RG:Z:p
      """
    And stderr does not contain "per-base tags"

    Examples:
      | ext  |
      | sam  |
      | bam  |
      | cram |

  Scenario: a trimming template keeps read-level tags and leaves out per-base ones, with one warning
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGTCG	IIIIIIIIII	qs:f:20	ch:i:7	st:Z:2023-06-22T07:17:48.308+00:00	MN:i:10	MM:Z:C+m?,0,0;	ML:B:C,200,10	mv:B:c,5,1,0,1,1,0,1,1,0,1,1,0,1,1
      r2	4	*	0	0	*	*	0	0	TTGCATGCAA	IIIIIIIIII	np:i:12	rq:f:0.999	zm:i:42	fi:B:C,1,2,3,4,5,6,7,8,9,10	ip:B:C,1,2,3,4,5,6,7,8,9,10
      r3	4	*	0	0	*	*	0	0	GGATCCATGG	IIIIIIIIII	OQ:Z:ABCDEFGHIJ	RX:Z:GGGG	QX:Z:FFFF
      """
    When I run `unmux in.sam --pool p --extract body=0:2:end --template body --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	GTACGTCG	IIIIIIII	qs:f:20	ch:i:7	st:Z:2023-06-22T07:17:48.308+00:00	RG:Z:p
      r2	4	*	0	255	*	*	0	0	GCATGCAA	IIIIIIII	np:i:12	rq:f:0.999	zm:i:42	RG:Z:p
      r3	4	*	0	255	*	*	0	0	ATCCATGG	IIIIIIII	RX:Z:GGGG	QX:Z:FFFF	RG:Z:p
      """
    And stderr contains "--template changed the bases of some records, so their per-base tags were left out: MM, ML, MN, OQ, ip, fi, mv"
    And stderr contains "per-base tags were left out" exactly once

  Scenario Outline: a record keeps per-base tags only when it is its read unchanged
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGT	IIII	ch:i:7	MM:Z:C+m?,0;	ML:B:C,200
      """
    When I run `unmux in.sam --pool p <template> --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	<bases>	IIII	<tags>	RG:Z:p
      """

    # ACGT is its own reverse complement, so only its orientation tells ~r apart from r.
    Examples:
      | template                                             | bases | tags                             |
      | --extract a=0:0:2 --extract b=0:2:end --template a+b | ACGT  | ch:i:7	MM:Z:C+m?,0;	ML:B:C,200 |
      | --extract r=0:0:end --template ~r                    | ACGT  | ch:i:7                           |
      | --extract a=0:0:2 --extract b=0:2:end --template b+a | GTAC  | ch:i:7                           |

  Scenario: each template record of a pair keeps the tags of its own read
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      q1	77	*	0	0	*	*	0	0	AAAACCCC	IIIIIIII	RX:Z:GG-TT	ZA:Z:one	OQ:Z:ABCDEFGH
      q1	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ	RX:Z:GG-TT	ZA:Z:two	OQ:Z:HGFEDCBA
      """
    When I run `unmux in.sam --pool p --extract a=0:0:end --extract b=1:2:end --template a b --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      q1	77	*	0	255	*	*	0	0	AAAACCCC	IIIIIIII	RX:Z:GG-TT	ZA:Z:one	OQ:Z:ABCDEFGH	RG:Z:p
      q1	141	*	0	255	*	*	0	0	GGTTTT	JJJJJJ	RX:Z:GG-TT	ZA:Z:two	RG:Z:p
      """
    And stderr contains "per-base tags were left out: OQ"

  Scenario: a record that joins two reads takes the tags of the read that supplies its first segment
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      q1	77	*	0	0	*	*	0	0	AAAACCCC	IIIIIIII	RX:Z:GG-TT	ZA:Z:one	OQ:Z:ABCDEFGH
      q1	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ	RX:Z:GG-TT	ZA:Z:two	OQ:Z:HGFEDCBA
      """
    When I run `unmux in.sam --pool p --extract a=0:0:4 --extract b=1:0:4 --template b+a --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      q1	4	*	0	255	*	*	0	0	GGGGAAAA	JJJJIIII	RX:Z:GG-TT	ZA:Z:two	RG:Z:p
      """
    And stderr contains "per-base tags were left out: OQ"

  Scenario: FASTQ comment tags stay on a trimmed record
    Given a file "in.fq" containing:
      """
      @r1	BC:Z:ACGT	RX:Z:TGATAA
      AAAACCCCGGGG
      +
      IIIIIIIIIIII
      """
    When I run `unmux in.fq --extract body=0:4:end --template body --out out.fq`
    Then the exit code is 0
    And the file "out.fq" is:
      """
      @r1 BC:Z:ACGT	RX:Z:TGATAA
      CCCCGGGG
      +
      IIIIIIII
      """

  Scenario: a --tag binding wins over a carried tag on a template record, as without a template
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AAAACCCCGGGG	IIIIIIIIIIII	CB:Z:OLD	CY:Z:FFF	ch:i:7
      r2	4	*	0	0	*	*	0	0	TTTTCCCCGGGG	IIIIIIIIIIII	CB:Z:OLD	CY:Z:FFF	ch:i:8
      """
    # @r1 gets CB from the tag ID and loses the carried CY, which described the old CB; @r2 misses
    # bc, so its CB binding is left out along with the carried CB and CY.
    When I run `unmux in.sam --pool p --group bc={x=AAAA} --group bc::loc=0:0:4 --extract cell=@bc::id --extract body=0:4:end --tag CB=cell --template body --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	CCCCGGGG	IIIIIIII	CB:Z:x	ch:i:7	RG:Z:p
      r2	4	*	0	255	*	*	0	0	CCCCGGGG	IIIIIIII	ch:i:8	RG:Z:p
      """
