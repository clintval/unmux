Feature: Input and output format handling
  Input auto-detection by byte-sniffing (FASTA, FASTQ, gzip, SAM, BAM, CRAM); stdin via '-';
  interleaving auto-detected from read-name pairing (--per-record opts out, no --interleaved flag);
  output format inferred from extension (FASTX/SAM/BAM/CRAM); --compression 0-9 for BGZF/gzip;
  stdout via --out - (mirrors the input format). With no --extract/--template the record body is the
  raw input read, so these pass-through scenarios isolate format detection and output framing.

  Scenario: A .bam extension produces a BGZF-framed binary file (not plain text)
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    # A .fq/.fastq output would be plain text; the .bam extension selects BGZF framing, so the file is
    # gzip-magic'd and the read bases are no longer present as plaintext (they are 4-bit packed). The
    # gzip + no-plaintext checks alone would also pass for a gzipped FASTX (unmux writes .fq.gz as BGZF
    # too), so we additionally decode the header via noodles: only a real BAM parses, carrying the @PG
    # provenance and (for a pass-through run) the default pool @RG. Full record decoding is in the
    # writer unit tests.
    When I run `unmux reads.fq --out output.bam`
    Then the exit code is 0
    And a file "output.bam" exists
    And the file "output.bam" is gzip-compressed
    And the file "output.bam" does not contain "ACGTACGTAACCGGTTACGT"
    And the BAM header of "output.bam" has exactly 1 @PG line
    And the BAM header of "output.bam" contains "@RG"

  Scenario: Write gzip-compressed FASTQ by extension
    Given a file "reads.fastq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fastq --out merged.fastq.gz`
    Then the exit code is 0
    And a file "merged.fastq.gz" exists
    And the file "merged.fastq.gz" is gzip-compressed

  Scenario: FASTA input is detected by content even with a neutral file extension
    Given a file "reads.txt" containing:
      """
      >r1
      ACGTACGTAACCGGTTACGT
      """
    # The .txt extension carries no format hint, so a successful FASTA mirror proves byte-sniffing.
    When I run `unmux reads.txt`
    Then the exit code is 0
    And stdout contains ">r1"
    And stdout contains "ACGTACGTAACCGGTTACGT"

  Scenario: A FASTQ record whose quality length differs from its sequence is a clean error
    Given a file "bad.fq" containing:
      """
      @r1
      ACGTACGTACGT
      +
      IIIIIIII
      """
    # The sequence is 12 bases but the quality string is 8 characters - a malformed FASTQ record. unmux
    # rejects it at read time with a clear message naming the offending record, rather than carrying the
    # length mismatch downstream where a quality slice would panic with an out-of-bounds index.
    When I run `unmux bad.fq`
    Then the exit code is 1
    And stderr contains "FASTQ"
    And stderr contains "quality"

  Scenario: SAM input is auto-detected and stdout mirrors the input format
    Given a file "unmapped.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGTAACCGGTTACGT	IIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux unmapped.sam`
    Then the exit code is 0
    And stdout contains "r1"
    And stdout contains "ACGTACGTAACCGGTTACGT"
    And stdout contains "ID:unmux"

  Scenario: Stream to stdout via --out -
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --out -`
    Then the exit code is 0
    And stdout contains "ACGTACGTAACCGGTTACGT"

  Scenario: Interleaved paired-end input is auto-detected from read names
    Given a file "interleaved.fq" containing:
      """
      @read1/1
      ACGTACGT
      +
      IIIIIIII
      @read1/2
      TTTTGGGG
      +
      IIIIIIII
      """
    When I run `unmux interleaved.fq --out out.R%ordinal.fq`
    Then the exit code is 0
    And a file "out.R1.fq" exists
    And a file "out.R2.fq" exists
    And the file "out.R1.fq" contains "ACGTACGT"
    And the file "out.R2.fq" contains "TTTTGGGG"

  Scenario: --per-record disables auto pair-detection so no second mate file is written
    Given a file "reads.fq" containing:
      """
      @read1/1
      ACGTACGT
      +
      IIIIIIII
      @read1/2
      TTTTGGGG
      +
      IIIIIIII
      """
    # The same paired-looking input auto-pairs into R1+R2 above; with --per-record each record is its
    # own single-end read, so only R1 is written (R2 would exist if pairing were NOT disabled) and
    # both reads land in R1.
    When I run `unmux reads.fq --per-record --out out.R%ordinal.fq`
    Then the exit code is 0
    And a file "out.R1.fq" exists
    And a file "out.R2.fq" does not exist
    And the file "out.R1.fq" contains "ACGTACGT"
    And the file "out.R1.fq" contains "TTTTGGGG"

  Scenario: CRAM output format by extension and SAM input
    Given a file "unmapped.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGTAACCGGTTACGT	IIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux unmapped.sam --out output.cram`
    Then the exit code is 0
    And a file "output.cram" exists
    And the file "output.cram" does not contain "ACGTACGTAACCGGTTACGT"

  Scenario: --threads 1 runs fully serially and still demultiplexes
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --threads 1 --out out.fq`
    Then the exit code is 0
    And a file "out.fq" exists
    And the file "out.fq" contains "ACGTACGTAACCGGTTACGT"

  Scenario: --threads 2 writes compressed output via inline compression
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --threads 2 --out out.bam`
    Then the exit code is 0
    And a file "out.bam" exists
    And the file "out.bam" is gzip-compressed
    And the file "out.bam" does not contain "ACGTACGTAACCGGTTACGT"

  Scenario: --compression 0 stores the gzip output so it is larger than the --compression 9 variant
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # The .gz extension always selects gzip/BGZF framing, so --compression only sets the deflate level,
    # not whether there is a gzip header. Level 0 stores the (highly repetitive) bases, so its file is far
    # larger than the level-9 deflated variant; both still carry the gzip magic and decode to the same
    # records. This is the only level difference the text harness can observe (it cannot decode the block).
    When I run `unmux reads.fq --compression 0 --out stored.fastq.gz`
    Then the exit code is 0
    And the file "stored.fastq.gz" is gzip-compressed
    When I run `unmux reads.fq --compression 9 --out packed.fastq.gz`
    Then the exit code is 0
    And the file "packed.fastq.gz" is gzip-compressed
    And the file "stored.fastq.gz" is larger than the file "packed.fastq.gz"

  Scenario: An extensionless --out mirrors the input format instead of erroring
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    # The --out path carries no suffix, so the output format follows the input (FASTQ here), the same
    # as stdout does, rather than failing to infer a format. An explicit suffix (e.g. --out x.bam)
    # would still win over the input format.
    When I run `unmux reads.fq --out results`
    Then the exit code is 0
    And a file "results" exists
    And the file "results" contains "@r1"
    And the file "results" contains "ACGTACGTAACCGGTTACGT"

  Scenario: FASTQ converts to FASTA, dropping the qualities
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    # A .fa output is FASTA, which has no quality line: the > record and bases carry over but the
    # per-base qualities are dropped, so there is no `+` separator in the output.
    When I run `unmux reads.fq --out out.fa`
    Then the exit code is 0
    And the file "out.fa" contains ">r1"
    And the file "out.fa" contains "ACGTACGTAACCGGTTACGT"
    And the file "out.fa" does not contain "+"

  Scenario: FASTA to FASTQ is rejected because FASTA carries no qualities
    Given a file "reads.fa" containing:
      """
      >r1
      ACGTACGTAACCGGTTACGT
      """
    # A .fq output needs per-base qualities, but a FASTA input has none. Rather than invent quality
    # scores, unmux fails fast with a clear message.
    When I run `unmux reads.fa --out out.fq`
    Then the exit code is 1
    And stderr contains "no qualities"

  Scenario: SAM input converts to a FASTQ output (alignment to FASTX)
    Given a file "unmapped.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGTAACCGGTTACGT	IIIIIIIIIIIIIIIIIIII
      """
    # An unmapped alignment record carries SEQ and QUAL, so it converts to a full FASTQ record (name,
    # bases, and qualities) when --out is .fq.
    When I run `unmux unmapped.sam --out out.fq`
    Then the exit code is 0
    And the file "out.fq" contains "@r1"
    And the file "out.fq" contains "ACGTACGTAACCGGTTACGT"
    And the file "out.fq" contains "IIIIIIIIIIIIIIIIIIII"

  Scenario: FASTQ converts to CRAM by extension
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    # The .cram extension selects CRAM container framing from a FASTX input; the bases are no longer
    # present as plaintext.
    When I run `unmux reads.fq --out out.cram`
    Then the exit code is 0
    And a file "out.cram" exists
    And the file "out.cram" does not contain "ACGTACGTAACCGGTTACGT"
