Feature: End-to-end demux flows over the real binary
  Concrete, executable acceptance scenarios that run the unmux binary on small fixtures and
  assert on what it emits, complementing the inline unit tests in src/lib. Every scenario across
  tests/features is executable; BAM/CRAM header internals are read back with `noodles` rather than
  parked.

  Scenario: Extract a UMI and carry it as an RX tag in the FASTQ comment
    Given a file "test.fq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux test.fq --extract umi=0:0:8 --extract body=0:8:end --tag RX=umi --template body`
    Then the exit code is 0
    And stdout contains "RX:Z:ACGTACGT"
    And stdout contains "AACCGGTTACGT"

  Scenario: Demux a pool into per-sample FASTQs with an unassigned catch-all
    Given a file "reads.fq" containing:
      """
      @a
      AAAACCCCDDDD
      +
      IIIIIIIIIIII
      @b
      CCCCGGGGTTTT
      +
      IIIIIIIIIIII
      @n
      GGGGTTTTAAAA
      +
      IIIIIIIIIIII
      """
    When I run `unmux reads.fq --group bc={s1=AAAA,s2=CCCC} --group bc::loc=0:0:4 --extract rest=0:4:end --template rest --sample S1=bc::s1 --sample S2=bc::s2 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And a file "S1.fq" exists
    And a file "S2.fq" exists
    And a file "un.fq" exists
    And the file "S1.fq" contains "CCCCDDDD"
    And the file "un.fq" contains "GGGGTTTTAAAA"

  Scenario: Write an unmapped BAM by extension (BGZF container)
    Given a file "test.fq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux test.fq --extract body=0:0:end --template body --out out.bam`
    Then the exit code is 0
    And a file "out.bam" exists
    And the file "out.bam" is gzip-compressed

  Scenario: --require-samples-explain-all-tags rejects a sample selector that leaves a barcode unclaimed
    Given a file "reads.fq" containing:
      """
      @a
      AAAACCCC
      +
      IIIIIIII
      """
    When I run `unmux reads.fq --group grp={s1=AAAA,s2=TTTT} --group grp::loc=0:0:4 --extract bc=0:0:4 --template bc --sample a=grp::s1 --require-samples-explain-all-tags --out %sample.fq`
    Then the exit code is 1
    And stderr contains "require-samples-explain-all-tags"
    And stderr contains "grp::s2"

  Scenario: Write gzip-compressed FASTQ by extension
    Given a file "test.fq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux test.fq --extract body=0:0:end --template body --out out.fq.gz`
    Then the exit code is 0
    And the file "out.fq.gz" is gzip-compressed
