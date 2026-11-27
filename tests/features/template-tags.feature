Feature: Input tags on --template records
  A record built by --template keeps the input tags of the read it is cut from, as a record passed
  through whole does, so UMI and barcode tags, read groups, FASTX comment tags, and per-read vendor
  tags survive trimming. Simple per-base tags (OQ, U2, E2, PacBio ip and pw, fgbio consensus and
  methylation counts) follow the record's bases: each segment takes its slice of the tag, reversed
  where the segment is flipped, and E2 and consensus bases are also complemented. Other per-base
  tags (MM/ML base modifications, the ONT move table, PacBio strand kinetics and pulse tracks) stay
  only on a record that is its read unchanged; otherwise they are left out, and one warning per run
  lists them. A record that joins two reads, or holds a corrected barcode, loses every per-base
  tag. A --tag binding still wins over a carried tag.

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

  Scenario: a trimming template slices the per-base tags it can follow and leaves out the rest, with one warning
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGTCG	IIIIIIIIII	qs:f:20	ch:i:7	st:Z:2023-06-22T07:17:48.308+00:00	sm:f:96.5	pt:i:23	MN:i:10	MM:Z:C+m?,0,0;	ML:B:C,200,10	mv:B:c,5,1,0,1,1,0,1,1,0,1,1,0,1,1
      r2	4	*	0	0	*	*	0	0	TTGCATGCAA	IIIIIIIIII	np:i:12	rq:f:0.999	zm:i:42	fi:B:C,1,2,3,4,5,6,7,8,9,10	ip:B:C,1,2,3,4,5,6,7,8,9,10
      r3	4	*	0	0	*	*	0	0	GGATCCATGG	IIIIIIIIII	OQ:Z:ABCDEFGHIJ	RX:Z:GGGG	QX:Z:FFFF
      """
    When I run `unmux in.sam --pool p --extract body=0:2:end --template body --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	GTACGTCG	IIIIIIII	qs:f:20	ch:i:7	st:Z:2023-06-22T07:17:48.308+00:00	sm:f:96.5	pt:i:23	RG:Z:p
      r2	4	*	0	255	*	*	0	0	GCATGCAA	IIIIIIII	np:i:12	rq:f:0.999	zm:i:42	ip:B:C,3,4,5,6,7,8,9,10	RG:Z:p
      r3	4	*	0	255	*	*	0	0	ATCCATGG	IIIIIIII	OQ:Z:CDEFGHIJ	RX:Z:GGGG	QX:Z:FFFF	RG:Z:p
      """
    And stderr contains "--template changed the bases of some records, so their per-base tags were left out: MM, ML, MN, fi, mv"
    And stderr contains "per-base tags were left out" exactly once

  Scenario: a trimming template leaves out PacBio pulse and pileup tags
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGTCG	IIIIIIIIII	zm:i:42	sa:B:I,10,3	sx:B:C,0,0,0,1,0,0,0,0,0,0	pc:Z:ACgGTACGTCG	pq:Z:IIIIIIIIIII	pg:Z:IIIIIIIIIII	pm:B:S,3,3,4,3,3,3,4,3,3,3,4	pd:B:S,8,5,5,8,8,5,5,8,8,5,5	px:B:S,2,2,4,5,2,2,4,5,2,2,4	pe:B:C,0,0,2,0,0,0,0,0,0,0,0	sf:B:I,10,15,23,28,30,35,40,45,50,55,60	bs:B:f,1.5,2.5,3.5,4.5	pb:B:I,11
      """
    When I run `unmux in.sam --pool p --extract body=0:2:end --template body --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	GTACGTCG	IIIIIIII	zm:i:42	RG:Z:p
      """
    And stderr contains "per-base tags were left out: sa, sx, pc, pq, pg, pm, pd, px, pe, sf, bs, pb"

  Scenario: a trimmed record keeps the matching slices of OQ and ip
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGTCG	IIIIIIIIII	ch:i:7	OQ:Z:ABCDEFGHIJ	ip:B:C,10,20,30,40,50,60,70,80,90,100
      """
    When I run `unmux in.sam --pool p --extract r=0:2:-2 --template r --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	GTACGT	IIIIII	ch:i:7	OQ:Z:CDEFGH	ip:B:C,30,40,50,60,70,80	RG:Z:p
      """
    And stderr does not contain "per-base tags"

  Scenario: a flipped read reverses OQ, reverse-complements E2, and still leaves out MM and ML
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AAACCCGT	ABCDEFGH	ch:i:7	OQ:Z:abcdefgh	E2:Z:CCGTTAAG	cd:B:s,1,2,3,4,5,6,7,8	MM:Z:C+m?,0;	ML:B:C,200
      """
    When I run `unmux in.sam --pool p --extract r=0:0:end --template ~r --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGGGTTT	HGFEDCBA	ch:i:7	OQ:Z:hgfedcba	E2:Z:CTTAACGG	cd:B:s,8,7,6,5,4,3,2,1	RG:Z:p
      """
    And stderr contains "--template changed the bases of some records, so their per-base tags were left out: MM, ML"

  Scenario: a per-base tag whose length is not its read's is left out of a changed record
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII	ch:i:7	OQ:Z:ABCD	ip:B:C,1,2,3
      """
    When I run `unmux in.sam --pool p --extract r=0:2:end --template r --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	GTACGT	IIIIII	ch:i:7	RG:Z:p
      """
    And stderr contains "--template changed the bases of some records, so their per-base tags were left out: OQ, ip; OQ, ip did not match the read's length"

  Scenario: fgbio consensus tags follow a trimmed, flipped read, while PacBio's same-named per-read tags stay
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGTCG	IIIIIIIIII	zm:i:42	ac:B:i,1,0,1,0	bc:B:S,3,3	bq:i:100	bt:Z:ACGTACGT
      r2	4	*	0	0	*	*	0	0	AAACGG	IIIIII	ac:Z:AAACGG	bc:Z:AAACGA	aq:Z:ABCDEF	bq:Z:abcdef	ad:B:s,1,2,3,4,5,6	cu:B:s,0,1,0,2,0,3	bt:B:s,6,5,4,3,2,1
      """
    When I run `unmux in.sam --pool p --extract r=0:1:end --template ~r --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	CGACGTACG	IIIIIIIII	zm:i:42	ac:B:i,1,0,1,0	bc:B:S,3,3	bq:i:100	bt:Z:ACGTACGT	RG:Z:p
      r2	4	*	0	255	*	*	0	0	CCGTT	IIIII	ac:Z:CCGTT	bc:Z:TCGTT	aq:Z:FEDCB	bq:Z:fedcb	ad:B:s,6,5,4,3,2	cu:B:s,3,0,2,0,1	bt:B:s,1,2,3,4,5	RG:Z:p
      """
    And stderr does not contain "per-base tags"

  Scenario Outline: a stream from a reverse-strand match turns per-base tags into its orientation
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AACGTGTCCCAAAA	ABCDEFGHIJKLMN	ch:i:7	OQ:Z:abcdefghijklmn	MM:Z:C+m?,0;	ML:B:C,200
      """
    When I run `unmux in.sam --out tagged.bam`
    Then the exit code is 0
    When I run `unmux tagged.bam --pool p --group bc={x=GGGAC} --group bc::bothStrands=true --extract ins=@bc+0:end --template <template> --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	<bases>	<quals>	ch:i:7	OQ:Z:<oq>	RG:Z:p
      """
    And stderr contains "per-base tags were left out: MM, ML"

    # GGGAC matches the reverse strand at bases 5 to 10, so ins is bases 0 to 5 in its orientation.
    Examples:
      | template | bases | quals | oq    |
      | ins      | ACGTT | EDCBA | edcba |
      | ~ins     | AACGT | ABCDE | abcde |

  Scenario: a reverse-complemented palindrome that spells its read is not an unchanged copy
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTTCCCGGGAACGT	IIIIIIIIIIIIIIII	ch:i:7	OQ:Z:abcdefghijklmnop	MM:Z:C+m?,0;	ML:B:C,200
      """
    # GGGA matches the reverse strand at bases 4 to 8, so head is bases 0 to 4, ACGT, reverse-complemented.
    When I run `unmux in.sam --pool p --group bc={x=GGGA} --group bc::bothStrands=true --extract head=@bc+0:end --extract rest=0:4:end --template head+rest --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTTCCCGGGAACGT	IIIIIIIIIIIIIIII	ch:i:7	OQ:Z:dcbaefghijklmnop	RG:Z:p
      """
    And stderr contains "per-base tags were left out: MM, ML"

  Scenario Outline: a record keeps MM and ML only when it is its read unchanged
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

  Scenario: a rearranged record carries the rearranged slices of its per-base tags
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGACGACGACG	IIIIIIIIIIII	ch:i:7	OQ:Z:ABCDEFGHIJKL	MM:Z:C+m?,0;	ML:B:C,200
      """
    When I run `unmux in.sam --pool p --extract a=0:0:3 --extract b=0:3:end --template b+a --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGACGACGACG	IIIIIIIIIIII	ch:i:7	OQ:Z:DEFGHIJKLABC	RG:Z:p
      """
    And stderr contains "per-base tags were left out: MM, ML"

  Scenario Outline: per-base tags stay with the observed bases of a whole-read barcode, not its corrected form
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII	ch:i:7	OQ:Z:ABCDEFGH
      """
    When I run `unmux in.sam --pool p --group bc={x=ACGTACGT} --extract bc=@bc --template <template> --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	<tags>	RG:Z:p
      """

    Examples:
      | template     | tags                  |
      | bc           | ch:i:7                |
      | bc::raw=true | ch:i:7	OQ:Z:ABCDEFGH |

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
      q1	141	*	0	255	*	*	0	0	GGTTTT	JJJJJJ	RX:Z:GG-TT	ZA:Z:two	OQ:Z:FEDCBA	RG:Z:p
      """
    And stderr does not contain "per-base tags"

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

  Scenario Outline: an empty stream supplies no segment of a template record
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      q1	77	*	0	0	*	*	0	0	AAAACCCC	IIIIIIII	ZA:Z:one	OQ:Z:ABCDEFGH
      q1	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ	ZA:Z:two	OQ:Z:HGFEDCBA
      """
    When I run `unmux in.sam --pool p <extracts> --template a+b --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      q1	4	*	0	255	*	*	0	0	<bases>	<tags>	RG:Z:p
      """
    And stderr does not contain "per-base tags"

    Examples:
      | extracts                              | bases             | tags                    |
      | --extract a=0:0:0 --extract b=1:0:end | GGGGTTTT	JJJJJJJJ | ZA:Z:two	OQ:Z:HGFEDCBA |
      | --extract a=0:0:end --extract b=1:0:0 | AAAACCCC	IIIIIIII | ZA:Z:one	OQ:Z:ABCDEFGH |

  Scenario: a template record with no bases takes the tags of the read its stream is cut from
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      q1	77	*	0	0	*	*	0	0	AAAACCCC	IIIIIIII	ZA:Z:one
      q1	141	*	0	0	*	*	0	0	*	*	ZA:Z:two
      """
    When I run `unmux in.sam --pool p --extract a=0:0:end --extract b=1:0:end --template a b --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      q1	77	*	0	255	*	*	0	0	AAAACCCC	IIIIIIII	ZA:Z:one	RG:Z:p
      q1	141	*	0	255	*	*	0	0	*	*	ZA:Z:two	RG:Z:p
      """

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
    # @r1 takes CB from the tag ID and drops its stale CY; @r2 misses bc, so its CB and CY go.
    When I run `unmux in.sam --pool p --group bc={x=AAAA} --group bc::loc=0:0:4 --extract cell=@bc::id --extract body=0:4:end --tag CB=cell --template body --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	CCCCGGGG	IIIIIIII	CB:Z:x	ch:i:7	RG:Z:p
      r2	4	*	0	255	*	*	0	0	CCCCGGGG	IIIIIIII	ch:i:8	RG:Z:p
      """
