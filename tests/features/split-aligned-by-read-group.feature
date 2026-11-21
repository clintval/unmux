Feature: Split an aligned input by read group, keeping alignments
  An `@RG` split of an aligned input (one whose header carries `@SQ`) into
  SAM/BAM/CRAM passes every record through unchanged: flags, position, CIGAR,
  MAPQ, mate fields, and every tag. Each output keeps the input header (`@HD`
  with its sort order, `@SQ`, `@PG`) with only its own `@RG` lines and the unmux
  `@PG` chained on, and records keep their input order. A `--template` would
  rewrite aligned records, so it is refused.

  Background:
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @SQ	SN:chr2	LN:800
      @RG	ID:rg1	SM:sampleA	LB:lib1
      @RG	ID:rg2	SM:sampleB	LB:lib2
      @PG	ID:bwa	PN:bwa	VN:0.7.17
      q1	99	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	MC:Z:8M	NM:i:0	RG:Z:rg1
      q2	163	chr1	120	30	4M1I3M	=	200	88	CCCCAAAA	IIIIHHHH	MC:Z:8M	RG:Z:rg2
      q1	147	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	MC:Z:8M	NM:i:1	RG:Z:rg1
      q2	83	chr1	200	30	8M	=	120	-88	TTTTGGGG	IIIIIIII	MC:Z:4M1I3M	RG:Z:rg2
      q3	129	chr1	500	12	8M	chr2	40	0	AACCGGTT	IIIIIIII	MC:Z:2S6M	RG:Z:rg1
      q4	256	chr1	600	0	8M	*	0	0	*	*	RG:Z:rg2
      q3	65	chr2	40	12	2S6M	chr1	500	0	AACCGGTT	IIIIIIII	MC:Z:8M	RG:Z:rg1
      q5	2048	chr2	90	7	3H5M	*	0	0	TTTTT	IIIII	SA:Z:chr1,10,+,5M3S,60,0;	RG:Z:rg2
      q6	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:rg2
      """

  Scenario: a paired aligned split by @RG::SM keeps every alignment field
    When I run `unmux in.sam --group sm=@RG::SM --sample-from-group sm --out out/%sample.bam`
    Then the exit code is 0
    And the records of "out/sampleA.bam" are:
      """
      q1	99	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	MC:Z:8M	NM:i:0	RG:Z:rg1
      q1	147	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	MC:Z:8M	NM:i:1	RG:Z:rg1
      q3	129	chr1	500	12	8M	chr2	40	0	AACCGGTT	IIIIIIII	MC:Z:2S6M	RG:Z:rg1
      q3	65	chr2	40	12	2S6M	chr1	500	0	AACCGGTT	IIIIIIII	MC:Z:8M	RG:Z:rg1
      """
    And the BAM header of "out/sampleA.bam" contains "SO:coordinate"
    And the BAM header of "out/sampleA.bam" contains "@SQ"
    And the BAM header of "out/sampleA.bam" contains "SN:chr2"
    And the BAM header of "out/sampleA.bam" contains "ID:rg1"
    And the BAM header of "out/sampleA.bam" does not contain "ID:rg2"
    And the BAM header of "out/sampleA.bam" contains "PN:bwa"
    And the BAM header of "out/sampleA.bam" contains "PP:bwa"
    And the BAM header of "out/sampleA.bam" has exactly 2 @PG lines
    And the BAM header of "out/sampleA.bam" contains "PN:unmux"

  Scenario: records keep coordinate order, with secondary, supplementary, and unmapped records
    When I run `unmux in.sam --group sm=@RG::SM --sample-from-group sm --out out/%sample.sam`
    Then the exit code is 0
    And the records of "out/sampleB.sam" are:
      """
      q2	163	chr1	120	30	4M1I3M	=	200	88	CCCCAAAA	IIIIHHHH	MC:Z:8M	RG:Z:rg2
      q2	83	chr1	200	30	8M	=	120	-88	TTTTGGGG	IIIIIIII	MC:Z:4M1I3M	RG:Z:rg2
      q4	256	chr1	600	0	8M	*	0	0	*	*	RG:Z:rg2
      q5	2048	chr2	90	7	3H5M	*	0	0	TTTTT	IIIII	SA:Z:chr1,10,+,5M3S,60,0;	RG:Z:rg2
      q6	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:rg2
      """
    And the file "out/sampleB.sam" contains "SO:coordinate"

  Scenario: an aligned split reads SAM from stdin
    When I run `unmux - --group rg=@RG --sample-from-group rg --out out/%sample.bam` with "in.sam" on stdin
    Then the exit code is 0
    And the records of "out/rg1.bam" are:
      """
      q1	99	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	MC:Z:8M	NM:i:0	RG:Z:rg1
      q1	147	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	MC:Z:8M	NM:i:1	RG:Z:rg1
      q3	129	chr1	500	12	8M	chr2	40	0	AACCGGTT	IIIIIIII	MC:Z:2S6M	RG:Z:rg1
      q3	65	chr2	40	12	2S6M	chr1	500	0	AACCGGTT	IIIIIIII	MC:Z:8M	RG:Z:rg1
      """
    And the BAM header of "out/rg1.bam" contains "SN:chr1"

  Scenario: a record with no known read group passes through to the unassigned bin
    Given a file "orphan.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      r1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      r2	16	chr1	300	50	8M	*	0	0	CCCCCCCC	IIIIIIII
      """
    When I run `unmux orphan.sam --group rg=@RG --sample-from-group rg --out out/%sample.sam --unassigned out/unassigned.sam`
    Then the exit code is 0
    And the records of "out/unassigned.sam" are:
      """
      r2	16	chr1	300	50	8M	*	0	0	CCCCCCCC	IIIIIIII
      """
    And the file "out/unassigned.sam" contains "SN:chr1"

  Scenario: metrics count every passed-through record
    When I run `unmux in.sam --pool p1 --group sm=@RG::SM --sample-from-group sm --out out/%sample.bam --metrics-per-sample p.tsv --metrics-summary s.tsv`
    Then the exit code is 0
    And the file "p.tsv" is:
      """
      pool	sample	sub_sample	reads	frac_of_pool	frac_bases_unextracted
      p1	sampleA		4	0.444444	0.000000
      p1	sampleB		5	0.555556	0.000000
      """
    And the file "s.tsv" contains "reads_assigned"

  Scenario: a --template on an aligned split into BAM is refused
    When I run `unmux in.sam --group sm=@RG::SM --sample-from-group sm --extract head=0:0:4 --template head --out out/%sample.bam`
    Then the exit code is 1
    And stderr contains "--template would rewrite"
    And a file "out/sampleA.bam" does not exist

  Scenario: a --template on an aligned split into FASTQ still works
    When I run `unmux in.sam --group sm=@RG::SM --sample-from-group sm --extract head=0:0:4 --template head --out out/%sample.fq`
    Then the exit code is 0
    And the file "out/sampleA.fq" contains "ACGT"

  Scenario: an aligned split into CRAM needs --reference
    When I run `unmux in.sam --group sm=@RG::SM --sample-from-group sm --out out/%sample.cram`
    Then the exit code is 1
    And stderr contains "--reference"

  Scenario: an aligned split into CRAM round-trips with --reference
    Given a file "ref.fa" containing:
      """
      >chr1
      ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
      """
    And a file "ref.fa.fai" containing:
      """
      chr1	40	6	40	41
      """
    And a file "small.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:40
      @RG	ID:rg1	SM:sampleA
      @RG	ID:rg2	SM:sampleB
      r1	99	chr1	1	60	8M	=	11	19	ACGTACGA	ABCDEFGH	MC:Z:8M	RG:Z:rg1
      r2	0	chr1	5	40	2S6M	*	0	0	TTACGTAC	IIIIIIII	RG:Z:rg2
      r1	147	chr1	11	60	3M1D5M	=	1	-19	GTAGTACG	HGFEDCBA	MC:Z:8M	RG:Z:rg1
      """
    When I run `unmux small.sam --group rg=@RG --sample-from-group rg --out out/%sample.cram --reference ref.fa`
    Then the exit code is 0
    When I run `unmux out/rg1.cram --group rg=@RG --sample-from-group rg --out back/%sample.sam --reference ref.fa`
    Then the exit code is 0
    And the records of "back/rg1.sam" are:
      """
      r1	99	chr1	1	60	8M	=	11	19	ACGTACGA	ABCDEFGH	MC:Z:8M	RG:Z:rg1
      r1	147	chr1	11	60	3M1D5M	=	1	-19	GTAGTACG	HGFEDCBA	MC:Z:8M	RG:Z:rg1
      """
    And the file "back/rg1.sam" contains "SN:chr1"
