Feature: Split an input by a per-record SAM auxiliary tag
  A `--group NAME=@tag::XX` source (a literal, case-sensitive prefix) routes
  each record by the value of its own two-character SAM aux tag XX and,
  paired with --sample-from-group, writes one file per distinct value. When
  the input's header proves it is already grouped by that tag, unmux streams
  one output file at a time, closing each before the next opens; otherwise it
  keeps every value's file open at once. Output is unmapped, exactly like an
  RG-based split.

  Scenario: a grouped input streams one file per cell barcode, keeping values apart
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	SS:unsorted:CB:coordinate
      @RG	ID:rg1	SM:sampleA
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1	CB:Z:cell1
      r2	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII	RG:Z:rg1	CB:Z:cell2
      """
    When I run `unmux in.sam --group cb=@tag::CB --sample-from-group cb --out out/%sample.sam`
    Then the exit code is 0
    And a file "out/cell1.sam" exists
    And a file "out/cell2.sam" exists
    And the file "out/cell1.sam" contains "RG:Z:rg1"
    And the file "out/cell1.sam" contains "@CO"
    And the file "out/cell1.sam" does not contain "CB:Z:cell2"
    And the file "out/cell2.sam" does not contain "CB:Z:cell1"

  Scenario: a cell barcode that reappears after its group closed is a fatal error
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	SS:unsorted:CB:coordinate
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	CB:Z:cell1
      r2	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII	CB:Z:cell2
      r3	4	*	0	0	*	*	0	0	GGGGGGGG	IIIIIIII	CB:Z:cell1
      """
    When I run `unmux in.sam --group cb=@tag::CB --sample-from-group cb --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "reappeared"

  Scenario: an ungrouped input still produces one correct file per value, all open at once
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	CB:Z:cell1
      r2	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII	CB:Z:cell2
      r3	4	*	0	0	*	*	0	0	GGGGGGGG	IIIIIIII	CB:Z:cell1
      """
    When I run `unmux in.sam --group cb=@tag::CB --sample-from-group cb --out out/%sample.sam`
    Then the exit code is 0
    And a file "out/cell1.sam" exists
    And a file "out/cell2.sam" exists
    And the file "out/cell1.sam" contains "CB:Z:cell1"
    And the file "out/cell1.sam" does not contain "CB:Z:cell2"
    And the file "out/cell2.sam" contains "CB:Z:cell2"
    And the file "out/cell2.sam" does not contain "CB:Z:cell1"

  Scenario: a record with no cell barcode lands in the unassigned bin
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      @RG	ID:rg1	SM:sampleA
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1	CB:Z:cell1
      r2	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII	RG:Z:rg1
      """
    When I run `unmux in.sam --group cb=@tag::CB --sample-from-group cb --out out/%sample.sam --unassigned out/un.sam`
    Then the exit code is 0
    And a file "out/cell1.sam" exists
    And a file "out/un.sam" exists

  Scenario: an aligned input is split but warns that output is unmapped
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      r1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	IIIIIIII	CB:Z:cell1
      """
    When I run `unmux in.sam --group cb=@tag::CB --sample-from-group cb --out out/%sample.bam`
    Then the exit code is 0
    And a file "out/cell1.bam" exists
    And stderr contains "UNMAPPED"
    And the BAM header of "out/cell1.bam" does not contain "@SQ"

  Scenario: a FASTQ whose read-name comments carry the tag can still be split
    Given a file "reads.fq" containing:
      """
      @r1 CB:Z:cell1
      AAAAAAAA
      +
      IIIIIIII
      @r2 CB:Z:cell2
      CCCCCCCC
      +
      IIIIIIII
      """
    When I run `unmux reads.fq --group cb=@tag::CB --sample-from-group cb --out out/%sample.sam`
    Then the exit code is 0
    And a file "out/cell1.sam" exists
    And a file "out/cell2.sam" exists
    And the file "out/cell1.sam" contains "CB:Z:cell1"
    And the file "out/cell2.sam" contains "CB:Z:cell2"

  Scenario: a split record keeps its own original read group, not another value's
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	SS:unsorted:CB:coordinate
      @RG	ID:rg1	SM:sampleA
      @RG	ID:rg2	SM:sampleB
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1	CB:Z:cell1
      r2	4	*	0	0	*	*	0	0	CCCCCCCC	IIIIIIII	RG:Z:rg2	CB:Z:cell2
      """
    When I run `unmux in.sam --group cb=@tag::CB --sample-from-group cb --out out/%sample.sam`
    Then the exit code is 0
    And the file "out/cell1.sam" contains "RG:Z:rg1"
    And the file "out/cell1.sam" contains "@CO"
    And the file "out/cell1.sam" does not contain "RG:Z:rg2"
    And the file "out/cell2.sam" contains "RG:Z:rg2"

  Scenario: a file path that merely starts with the tag-source spelling is not a tag source
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AAAAAAAA	IIIIIIII	CB:Z:cell1
      """
    When I run `unmux in.sam --group g=./@tag::CB --sample-from-group g --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "./@tag::CB"
