Feature: Edge cases of an aligned read-group split
  Regression cases for the alignment-preserving `@RG` split: odd `@PG` chains,
  FASTX outputs beside SAM/BAM ones, `--tag` bindings that would overwrite kept
  tags, and records CRAM output cannot hold.

  Scenario: a dangling @PG PP link does not crash the split
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      @PG	ID:a	PN:a	PP:missing
      r1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 0
    And the file "out/rg1.sam" contains "PN:unmux"
    And the file "out/rg1.sam" contains "PP:missing"

  Scenario: a @PG cycle leaves the unmux @PG unchained
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      @PG	ID:a	PN:a	PP:b
      @PG	ID:b	PN:b	PP:a
      r1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 0
    And the file "out/rg1.sam" contains "@PG	ID:unmux	PN:unmux"

  Scenario: colliding @PG IDs get a counter-suffixed unmux ID
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      @PG	ID:unmux	PN:unmux
      @PG	ID:x	PN:x	PP:unmux
      @PG	ID:unmux-x	PN:other	PP:x
      r1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 0
    And the file "out/rg1.sam" contains "@PG	ID:unmux.1	PN:unmux"
    And the file "out/rg1.sam" contains "PP:unmux-x"

  Scenario Outline: FASTX beside SAM/BAM is refused before any file is written: <outputs>
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	GO:query
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      q1	99	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	RG:Z:rg1
      q1	147	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg <outputs>`
    Then the exit code is 1
    And stderr contains "this aligned @RG split cannot write both SAM/BAM, which keeps each record whole, and FASTX, which pairs mates: it writes SAM/BAM to <sam_bam> and FASTX to <fastx>; make every output SAM/BAM, or every output FASTX"
    And the directory "." holds exactly "in.sam"

    Examples:
      | outputs                                                                  | sam_bam                                            | fastx                          |
      | --out out/%sample.R%ordinal.fq --unassigned out/un.bam                   | "out/un.bam" (--unassigned)                        | "out/rg1.R1.fq" (--out)        |
      | --out out/%sample.sam --unassigned out/un.fq                             | "out/rg1.sam" (--out)                              | "out/un.fq" (--unassigned)     |
      | --out out/%sample.bam --remove rg=out/removed.fq                         | "out/rg1.bam" (--out)                              | "out/removed.fq" (--remove rg) |
      | --out out/%sample.sam --unassigned out/un.bam --remove rg=out/removed.fa | "out/rg1.sam" (--out), "out/un.bam" (--unassigned) | "out/removed.fa" (--remove rg) |
      | --unassigned out/un.fq                                                   | stdout                                             | "out/un.fq" (--unassigned)     |
      | --out - --unassigned out/un.fq                                           | stdout                                             | "out/un.fq" (--unassigned)     |
      | --out out/%sample.R%ordinal.fq --unassigned -                            | stdout                                             | "out/rg1.R1.fq" (--out)        |

  Scenario: FASTQ outputs and a FASTQ bin pair the mates and keep primary records only
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	GO:query
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      q1	99	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	RG:Z:rg1
      q1	355	chr1	300	0	8M	=	150	0	*	*	RG:Z:rg1
      q1	2145	chr1	500	9	4H4M	=	150	0	ACGT	ABCD	RG:Z:rg1
      q1	147	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	RG:Z:rg1
      q2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      q2	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.R%ordinal.fq --unassigned out/un.fq`
    Then the exit code is 0
    And the file "out/rg1.R1.fq" is:
      """
      @q1 RG:Z:rg1
      ACGTACGT
      +
      ABCDEFGH
      """
    And the file "out/rg1.R2.fq" is:
      """
      @q1 RG:Z:rg1
      AAAACCCC
      +
      ABCDEFGH
      """
    And the file "out/un.fq" is:
      """
      @q2
      CCCCAAAA
      +
      IIIIIIII
      @q2
      TTTTGGGG
      +
      JJJJJJJJ
      """

  Scenario: A bin sent to /dev/null beside FASTQ outputs is not SAM/BAM, so the mates stay paired
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	GO:query
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      q1	99	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	RG:Z:rg1
      q1	147	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	RG:Z:rg1
      q2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      q2	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.R%ordinal.fq --unassigned /dev/null`
    Then the exit code is 0
    And the file "out/rg1.R2.fq" is:
      """
      @q1 RG:Z:rg1
      AAAACCCC
      +
      ABCDEFGH
      """

  Scenario: An --out sent to /dev/null beside a FASTQ bin is not SAM/BAM, so the bin holds pairs
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	GO:query
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      q1	99	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	RG:Z:rg1
      q1	147	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	RG:Z:rg1
      q2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      q2	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out /dev/null --unassigned out/un.fq`
    Then the exit code is 0
    And the file "out/un.fq" is:
      """
      @q2
      CCCCAAAA
      +
      IIIIIIII
      @q2
      TTTTGGGG
      +
      JJJJJJJJ
      """

  Scenario: SAM and BAM outputs together keep every record whole
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	GO:query
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      q1	99	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	RG:Z:rg1
      q1	355	chr1	300	0	8M	=	150	0	*	*	RG:Z:rg1
      q1	2145	chr1	500	9	4H4M	=	150	0	ACGT	ABCD	RG:Z:rg1
      q1	147	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	RG:Z:rg1
      q2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      q2	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.sam --unassigned out/un.bam`
    Then the exit code is 0
    And the records of "out/rg1.sam" are:
      """
      q1	99	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	RG:Z:rg1
      q1	355	chr1	300	0	8M	=	150	0	*	*	RG:Z:rg1
      q1	2145	chr1	500	9	4H4M	=	150	0	ACGT	ABCD	RG:Z:rg1
      q1	147	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	RG:Z:rg1
      """
    And the records of "out/un.bam" are:
      """
      q2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      q2	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      """

  Scenario: mates next to each other still pass through one record at a time, in input order
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	GO:query
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      q1	163	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	MC:Z:8M	RG:Z:rg1
      q1	2131	chr1	400	9	4H4M	=	100	0	TTTT	IIII	RG:Z:rg1
      q1	83	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	MC:Z:8M	RG:Z:rg1
      """
    When I run `unmux in.sam --pool p1 --group rg=@RG --sample-from-group rg --out out/%sample.sam --metrics-per-sample p.tsv`
    Then the exit code is 0
    And the records of "out/rg1.sam" are:
      """
      q1	163	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	MC:Z:8M	RG:Z:rg1
      q1	2131	chr1	400	9	4H4M	=	100	0	TTTT	IIII	RG:Z:rg1
      q1	83	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	MC:Z:8M	RG:Z:rg1
      """
    And the file "p.tsv" is:
      """
      pool	sample	sub_sample	reads	frac_of_pool	frac_bases_unextracted
      p1	rg1		2	1.000000	0.000000
      """
    And stderr contains "1 secondary/supplementary records passed through"

  Scenario: a --tag that would replace a kept RG:Z is refused
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      r1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --extract head=0:0:4 --tag RG=head --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "`RG` would overwrite"
    And a file "out/rg1.sam" does not exist

  Scenario: a --tag that would replace an alignment tag is refused
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      r1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	IIIIIIII	NM:i:0	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --extract head=0:0:4 --tag NM=head --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "`NM` would overwrite"

  Scenario: a sequence-matching group is refused, so a secondary cannot be routed apart from its primary
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sA
      p1	0	chr1	100	60	12M	*	0	0	AACCGGTTACGT	IIIIIIIIIIII	RG:Z:rg1
      p1	256	chr1	300	0	12M	*	0	0	*	*	RG:Z:rg1
      p1	2048	chr1	500	9	8H4M	*	0	0	ACGT	IIII	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG::SM --group ad={X=AACCGGTT} --group ad::loc=0:0:8 --sample-from-group rg --remove ad=out/removed.sam --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "group `ad` matches record sequences"
    And stderr contains "samtools reset"
    And a file "out/sA.sam" does not exist
    And a file "out/removed.sam" does not exist

  Scenario: a sequence-matching group nothing consumes is allowed, with its no-effect warning
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sA
      p1	0	chr1	100	60	12M	*	0	0	AACCGGTTACGT	IIIIIIIIIIII	RG:Z:rg1
      p1	256	chr1	300	0	12M	*	0	0	*	*	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --group bc={x=ACGT} --group bc::loc=0:0:4 --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 0
    And stderr contains "group `bc` is defined but nothing consumes its matches"
    And the records of "out/rg1.sam" are:
      """
      p1	0	chr1	100	60	12M	*	0	0	AACCGGTTACGT	IIIIIIIIIIII	RG:Z:rg1
      p1	256	chr1	300	0	12M	*	0	0	*	*	RG:Z:rg1
      """

  Scenario: a sequence-matching group reported by --qc-tag is refused
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sA
      p1	0	chr1	100	60	12M	*	0	0	AACCGGTTACGT	IIIIIIIIIIII	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --group bc={x=ACGT} --group bc::loc=0:0:4 --sample-from-group rg --qc-tag XQ --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "group `bc` matches record sequences"

  Scenario: a --tag drawn from record sequences is refused
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      r1	0	chr1	100	60	8M	*	0	0	ACGTACGT	IIIIIIII	RG:Z:rg1
      r1	256	chr1	300	0	8M	*	0	0	*	*	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --extract head=0:0:4 --tag XH=head --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "--tag `XH` is drawn from record sequences"
    And a file "out/rg1.sam" does not exist

  Scenario: a --tag drawn from record sequences still works into FASTQ, which keeps primary records only
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      r1	0	chr1	100	60	8M	*	0	0	ACGTACGT	IIIIIIII	RG:Z:rg1
      r1	256	chr1	300	0	8M	*	0	0	*	*	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --extract head=0:0:4 --tag XH=head --tag XH::qual=none --out out/%sample.fq`
    Then the exit code is 0
    And the file "out/rg1.fq" is:
      """
      @r1 RG:Z:rg1	XH:Z:ACGT
      ACGTACGT
      +
      IIIIIIII
      """

  Scenario: mapped records with no SEQ or QUAL are refused for CRAM, leaving no files
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:40
      @RG	ID:rg1	SM:sampleA
      r1	0	chr1	1	60	8M	*	0	0	ACGTACGT	*	RG:Z:rg1
      r1	256	chr1	5	0	8M	*	0	0	*	*	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.cram --unassigned out/unassigned.cram`
    Then the exit code is 1
    And stderr contains "writing aligned records to CRAM is not supported"
    And a file "out/rg1.cram" does not exist
    And a file "out/unassigned.cram" does not exist
