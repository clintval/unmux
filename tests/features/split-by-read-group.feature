Feature: Split an existing multi-read-group input by its @RG header
  A `--group NAME=@RG[::SM[::LB]]` source reads the input's read groups and,
  paired with --sample-from-group, fans records out into one file per read
  group, sample, or library, keyed on the record's own RG:Z and the
  header's @RG lines. Output is unmapped, and aligned inputs are warned
  about but processed.

  Scenario: @RG splits by read-group ID, one file per read group
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      @RG	ID:rg1	SM:sampleA
      @RG	ID:rg2	SM:sampleB
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      r2	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII	RG:Z:rg2
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.bam`
    Then the exit code is 0
    And a file "out/rg1.bam" exists
    And a file "out/rg2.bam" exists
    And the BAM header of "out/rg1.bam" contains "ID:rg1"
    And the BAM header of "out/rg1.bam" contains "SM:sampleA"
    And the BAM header of "out/rg1.bam" does not contain "ID:rg2"

  Scenario: @RG::SM merges read groups that share a sample
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      @RG	ID:rg1	SM:sampleA	LB:lib1
      @RG	ID:rg2	SM:sampleB	LB:lib2
      @RG	ID:rg3	SM:sampleA	LB:lib3
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      r2	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII	RG:Z:rg2
      r3	4	*	0	0	*	*	0	0	GGGGGGGG	IIIIIIII	RG:Z:rg3
      """
    When I run `unmux in.sam --group sm=@RG::SM --sample-from-group sm --out out/%sample.bam`
    Then the exit code is 0
    And a file "out/sampleA.bam" exists
    And a file "out/sampleB.bam" exists
    And a file "out/rg1.bam" does not exist
    And the BAM header of "out/sampleA.bam" contains "ID:rg1"
    And the BAM header of "out/sampleA.bam" contains "ID:rg3"
    And the BAM header of "out/sampleA.bam" does not contain "ID:rg2"

  Scenario: @RG::SM::LB fans out two tiers into %sample.%sub_sample
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      @RG	ID:rg1	SM:sampleA	LB:lib1
      @RG	ID:rg2	SM:sampleA	LB:lib2
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      r2	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII	RG:Z:rg2
      """
    When I run `unmux in.sam --group both=@RG::SM::LB --sample-from-group both --out out/%sample.%sub_sample.bam`
    Then the exit code is 0
    And a file "out/sampleA.lib1.bam" exists
    And a file "out/sampleA.lib2.bam" exists

  Scenario: dropping a tier merges targets into one file, with a warning
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      @RG	ID:rg1	SM:sampleA	LB:lib1
      @RG	ID:rg2	SM:sampleA	LB:lib2
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      r2	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII	RG:Z:rg2
      """
    When I run `unmux in.sam --group both=@RG::SM::LB --sample-from-group both --out out/%sample.bam`
    Then the exit code is 0
    And a file "out/sampleA.bam" exists
    And stderr contains "merge"
    And the BAM header of "out/sampleA.bam" contains "ID:rg1"
    And the BAM header of "out/sampleA.bam" contains "ID:rg2"

  Scenario: an aligned input is split but warns that output is unmapped
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      r1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.bam`
    Then the exit code is 0
    And a file "out/rg1.bam" exists
    And stderr contains "UNMAPPED"
    And the BAM header of "out/rg1.bam" does not contain "@SQ"

  Scenario: a record with no RG:Z lands in the unassigned bin
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      @RG	ID:rg1	SM:sampleA
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      r2	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.bam --unassigned out/un.bam`
    Then the exit code is 0
    And a file "out/un.bam" exists

  Scenario: @RG on a FASTX input fails fast
    Given a file "reads.fq" containing:
      """
      @q1
      AAAAAAAA
      +
      IIIIIIII
      """
    When I run `unmux reads.fq --group rg=@RG --sample-from-group rg --out out/%sample.bam`
    Then the exit code is 1
    And stderr contains "@RG"

  Scenario: an input header with no @RG lines fails fast
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII
      """
    When I run `unmux in.sam --group rg=@RG --sample-from-group rg --out out/%sample.bam`
    Then the exit code is 1
    And stderr contains "@RG"

  Scenario: a non-@RG run over inputs with a shared RG id is unaffected
    Given a file "a.sam" containing:
      """
      @HD	VN:1.6
      @RG	ID:rg1	SM:sampleA
      q1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      """
    And a file "b.sam" containing:
      """
      @HD	VN:1.6
      @RG	ID:rg1	SM:sampleB
      q1	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII	RG:Z:rg1
      """
    When I run `unmux --in 0=a.sam --in 1=b.sam --out out.bam`
    Then the exit code is 0
    And a file "out.bam" exists
