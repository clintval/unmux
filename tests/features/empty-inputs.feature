Feature: Empty and zero-read inputs never fail
  A run whose inputs contain zero reads (an empty FASTA/FASTQ, an empty-gzip FASTX, or a
  header-only SAM/BAM/CRAM) must succeed and emit zero records, never a format-detection error.
  This is the no-data-from-sequencer case: the tool must not fail, and directed outputs must
  still be created. Binary inputs are produced by unmux itself (from a header-only SAM or a
  0-byte FASTQ), exercising a real writer-to-reader round-trip.

  Scenario: An empty FASTQ file is read as zero records
    Given an empty file "empty.fastq"
    When I run `unmux empty.fastq --out out.fq`
    Then the exit code is 0
    And a file "out.fq" exists

  Scenario: An empty FASTA file is read as zero records
    Given an empty file "empty.fasta"
    When I run `unmux empty.fasta --out out.fa`
    Then the exit code is 0
    And a file "out.fa" exists

  Scenario: An empty-gzip FASTQ (the empty.R1.fastq.gz case) is read as zero records
    Given an empty file "seed.fastq"
    # unmux's own FASTX.gz writer makes a valid empty .fastq.gz, which is then read back.
    When I run `unmux seed.fastq --out empty.fastq.gz`
    Then the exit code is 0
    And the file "empty.fastq.gz" is gzip-compressed
    When I run `unmux empty.fastq.gz --out roundtrip.fq`
    Then the exit code is 0
    And a file "roundtrip.fq" exists

  Scenario: A header-only SAM (zero records) is read as zero records
    Given a file "headeronly.sam" containing:
      """
      @HD	VN:1.6
      """
    When I run `unmux headeronly.sam --out out.sam`
    Then the exit code is 0
    And a file "out.sam" exists

  Scenario: A header-only BAM (zero records) is read as zero records
    Given a file "headeronly.sam" containing:
      """
      @HD	VN:1.6
      """
    # unmux writes a valid empty BAM (header + EOF), which is then read back as zero records.
    When I run `unmux headeronly.sam --out empty.bam`
    Then the exit code is 0
    When I run `unmux empty.bam --out out.bam`
    Then the exit code is 0
    And a file "out.bam" exists

  Scenario: A header-only CRAM (zero records) is read as zero records
    Given a file "headeronly.sam" containing:
      """
      @HD	VN:1.6
      """
    When I run `unmux headeronly.sam --out empty.cram`
    Then the exit code is 0
    When I run `unmux empty.cram --out out.cram`
    Then the exit code is 0
    And a file "out.cram" exists

  Scenario: A zero-read run still creates the fan-out, removed, and unassigned files
    Given an empty file "empty.fq"
    When I run `unmux empty.fq --group grp={A=AAAAAAAA,C=GGGGGGGG} --remove grp::C=removed.fq --unassigned un.fq --out out.fq`
    Then the exit code is 0
    And a file "out.fq" exists
    And a file "removed.fq" exists
    And a file "un.fq" exists
