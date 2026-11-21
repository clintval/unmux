Feature: Edge cases of an aligned read-group split
  Regression cases for the alignment-preserving `@RG` split: odd `@PG` chains,
  secondary and supplementary records next to FASTX outputs, `--tag` bindings
  that would overwrite kept tags, and records CRAM output cannot hold.

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

  Scenario: FASTQ outputs stay primary-only when only a bin is SAM/BAM/CRAM
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      p1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      p1	256	chr1	300	0	8M	*	0	0	*	*	RG:Z:rg1
      p1	2048	chr1	500	9	4H4M	*	0	0	AAAA	IIII	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.fq --unassigned out/unassigned.bam`
    Then the exit code is 0
    And the file "out/rg1.fq" is:
      """
      @p1 RG:Z:rg1
      AAAAAAAA
      +
      IIIIIIII
      """

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
