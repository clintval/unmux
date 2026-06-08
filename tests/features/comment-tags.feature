Feature: Lifting SAM tags from FASTX read-name comments
  A FASTX read-name comment may carry SAM auxiliary tags (the samtools `fastq -T`
  convention that unmux itself writes on FASTX output). On input those tags are
  lifted back into real SAM tags, the inverse of that output, so a tag-bearing
  FASTQ converts to a uBAM/SAM with real tags and no extraction configured.
  Lifting is lenient and per-field, matching `samtools import -T`: a field that
  is not a well-formed `XX:T:VALUE` tag (a bare UMI, an Illumina CASAVA string,
  free text) is skipped while the valid tags in the same comment still lift.

  Scenario: a single SAM tag in a FASTQ comment becomes a real SAM tag
    Given a file "in.fq" containing:
      """
      @read1 RX:Z:ACGT
      AAAACCCC
      +
      IIIIIIII
      """
    When I run `unmux in.fq --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "RX:Z:ACGT"

  Scenario: a non-tag field is skipped while the valid tags in the same comment lift
    Given a file "in.fq" containing:
      """
      @read1	notatag	RX:Z:ACGT	BC:Z:GGGG
      AAAACCCC
      +
      IIIIIIII
      """
    When I run `unmux in.fq --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "RX:Z:ACGT"
    And the file "out.sam" contains "BC:Z:GGGG"
    And the file "out.sam" does not contain "notatag"

  Scenario: a bare UMI in the comment is left alone, not lifted into a tag
    Given a file "in.fq" containing:
      """
      @read1 GATTACAUMI
      AAAACCCC
      +
      IIIIIIII
      """
    When I run `unmux in.fq --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "AAAACCCC"
    And the file "out.sam" does not contain "GATTACAUMI"
