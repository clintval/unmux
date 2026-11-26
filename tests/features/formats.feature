Feature: Input and output format handling
  Input auto-detection by byte-sniffing (FASTA, FASTQ, gzip, SAM, BAM, CRAM); stdin via '-';
  interleaving auto-detected from read-name pairing and SAM/BAM/CRAM pairs from the record flags
  (--per-record opts out, no --interleaved flag); output format inferred from extension
  (FASTX/SAM/BAM/CRAM); --compression 0-9 for BGZF/gzip; stdout via --out - (mirrors the input
  format). With no --extract/--template the record body is the raw input read, so these
  pass-through scenarios isolate format detection and output framing.

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

  Scenario: Paired SAM input is auto-detected from the record flags
    Given a file "paired.sam" containing:
      """
      @HD	VN:1.6
      read1	77	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      read1	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      read2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      read2	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ
      """
    # Flags 77 and 141 mark the two mates of an unmapped pair, and each pair is one read.
    When I run `unmux paired.sam --pool p --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      read1	77	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:p
      read1	141	*	0	255	*	*	0	0	TTTTGGGG	JJJJJJJJ	RG:Z:p
      read2	77	*	0	255	*	*	0	0	CCCCAAAA	IIIIIIII	RG:Z:p
      read2	141	*	0	255	*	*	0	0	GGGGTTTT	JJJJJJJJ	RG:Z:p
      """
    And stderr contains "2 reads (0 assigned, 2 pass-through"

  Scenario: Paired SAM input writes its mates to R1 and R2 FASTQ files
    Given a file "paired.sam" containing:
      """
      @HD	VN:1.6
      read1	77	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      read1	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      read2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      read2	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ
      """
    When I run `unmux paired.sam --out out.R%ordinal.fq`
    Then the exit code is 0
    And the file "out.R1.fq" is:
      """
      @read1
      ACGTACGT
      +
      IIIIIIII
      @read2
      CCCCAAAA
      +
      IIIIIIII
      """
    And the file "out.R2.fq" is:
      """
      @read1
      TTTTGGGG
      +
      JJJJJJJJ
      @read2
      GGGGTTTT
      +
      JJJJJJJJ
      """

  Scenario Outline: A pair written to <format> reads back as a pair on <threads> threads
    Given a file "r1.fq" containing:
      """
      @read1
      ACGTACGT
      +
      IIIIIIII
      @read2
      CCCCAAAA
      +
      IIIIIIII
      """
    And a file "r2.fq" containing:
      """
      @read1
      TTTTGGGG
      +
      JJJJJJJJ
      @read2
      GGGGTTTT
      +
      JJJJJJJJ
      """
    When I run `unmux r1.fq r2.fq --out pairs.<format>`
    Then the exit code is 0
    When I run `unmux pairs.<format> --pool p --threads <threads> --out back.sam`
    Then the exit code is 0
    And the records of "back.sam" are:
      """
      read1	77	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:p
      read1	141	*	0	255	*	*	0	0	TTTTGGGG	JJJJJJJJ	RG:Z:p
      read2	77	*	0	255	*	*	0	0	CCCCAAAA	IIIIIIII	RG:Z:p
      read2	141	*	0	255	*	*	0	0	GGGGTTTT	JJJJJJJJ	RG:Z:p
      """

    Examples:
      | format | threads |
      | sam    | 1       |
      | sam    | 2       |
      | sam    | 8       |
      | bam    | 1       |
      | bam    | 2       |
      | bam    | 8       |
      | bam    | 16      |
      | cram   | 1       |
      | cram   | 2       |
      | cram   | 8       |

  Scenario: A SAM pair whose last mate comes first is read first mate first
    Given a file "swapped.sam" containing:
      """
      @HD	VN:1.6
      read1	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      read1	77	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      read2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      read2	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ
      """
    When I run `unmux swapped.sam --pool p --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      read1	77	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:p
      read1	141	*	0	255	*	*	0	0	TTTTGGGG	JJJJJJJJ	RG:Z:p
      read2	77	*	0	255	*	*	0	0	CCCCAAAA	IIIIIIII	RG:Z:p
      read2	141	*	0	255	*	*	0	0	GGGGTTTT	JJJJJJJJ	RG:Z:p
      """
    When I run `unmux swapped.sam --out out.R%ordinal.fq`
    Then the exit code is 0
    And the file "out.R1.fq" is:
      """
      @read1
      ACGTACGT
      +
      IIIIIIII
      @read2
      CCCCAAAA
      +
      IIIIIIII
      """
    And the file "out.R2.fq" is:
      """
      @read1
      TTTTGGGG
      +
      JJJJJJJJ
      @read2
      GGGGTTTT
      +
      JJJJJJJJ
      """

  Scenario: A mapped SAM pair whose last mate is leftmost is read first mate first
    Given a file "mapped.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	GO:query
      @SQ	SN:chr1	LN:1000
      q1	163	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH
      q1	83	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA
      """
    When I run `unmux mapped.sam --pool p --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      q1	77	*	0	255	*	*	0	0	AAAACCCC	ABCDEFGH	RG:Z:p
      q1	141	*	0	255	*	*	0	0	ACGTACGT	ABCDEFGH	RG:Z:p
      """

  Scenario: A mapped SAM pair outside an aligned read-group split becomes an unmapped pair
    Given a file "mapped.sam" containing:
      """
      @HD	VN:1.6	SO:unsorted	GO:query
      @SQ	SN:chr1	LN:1000
      q1	99	chr1	100	60	8M	=	150	58	ACGTACGT	ABCDEFGH	MC:Z:8M	NM:i:0
      q1	2147	chr1	400	9	4H4M	=	150	0	ACGT	ABCD	NM:i:0
      q1	147	chr1	150	60	8M	=	100	-58	GGGGTTTT	HGFEDCBA	MC:Z:8M	NM:i:1
      """
    # The supplementary record is skipped and the reverse-strand mate is written as sequenced.
    When I run `unmux mapped.sam --pool p --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      q1	77	*	0	255	*	*	0	0	ACGTACGT	ABCDEFGH	RG:Z:p
      q1	141	*	0	255	*	*	0	0	AAAACCCC	ABCDEFGH	RG:Z:p
      """

  Scenario: File index 1 is the second mate of a paired SAM input
    Given a file "paired.sam" containing:
      """
      @HD	VN:1.6
      read1	77	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      read1	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      read2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      read2	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ
      """
    When I run `unmux paired.sam --pool p --extract r1=0:0:end --extract umi=1:0:4 --template r1 --tag RX=umi --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      read1	4	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RX:Z:TTTT	QX:Z:JJJJ	RG:Z:p
      read2	4	*	0	255	*	*	0	0	CCCCAAAA	IIIIIIII	RX:Z:GGGG	QX:Z:JJJJ	RG:Z:p
      """

  Scenario Outline: The QC-fail flag of a <format> pair is kept on <threads> threads
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      read1	589	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      read1	653	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      read2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      read2	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ
      """
    When I run `unmux in.sam --pool p --out in.<format>`
    Then the exit code is 0
    When I run `unmux in.<format> --pool p --threads <threads> --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      read1	589	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:p
      read1	653	*	0	255	*	*	0	0	TTTTGGGG	JJJJJJJJ	RG:Z:p
      read2	77	*	0	255	*	*	0	0	CCCCAAAA	IIIIIIII	RG:Z:p
      read2	141	*	0	255	*	*	0	0	GGGGTTTT	JJJJJJJJ	RG:Z:p
      """

    Examples:
      | format | threads |
      | sam    | 1       |
      | sam    | 8       |
      | bam    | 1       |
      | bam    | 8       |
      | cram   | 1       |
      | cram   | 8       |

  Scenario: The QC-fail flag is kept in a sample, an unassigned bin, and a removed bin
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      a	516	*	0	0	*	*	0	0	AAAACCCCGGGG	IIIIIIIIIIII
      b	516	*	0	0	*	*	0	0	CCCCGGGGTTTT	IIIIIIIIIIII
      c	516	*	0	0	*	*	0	0	GGGGTTTTAAAA	IIIIIIIIIIII
      d	4	*	0	0	*	*	0	0	AAAACCCCGGGG	IIIIIIIIIIII
      """
    When I run `unmux in.sam --pool p --group bc={s1=AAAA,phix=GGGG} --group bc::loc=0:0:4 --sample s1=bc::s1 --remove bc::phix=rm.sam --unassigned un.sam --out %sample.sam`
    Then the exit code is 0
    And the records of "s1.sam" are:
      """
      a	516	*	0	255	*	*	0	0	AAAACCCCGGGG	IIIIIIIIIIII	RG:Z:s1
      d	4	*	0	255	*	*	0	0	AAAACCCCGGGG	IIIIIIIIIIII	RG:Z:s1
      """
    And the records of "un.sam" are:
      """
      b	516	*	0	255	*	*	0	0	CCCCGGGGTTTT	IIIIIIIIIIII	RG:Z:p
      """
    And the records of "rm.sam" are:
      """
      c	516	*	0	255	*	*	0	0	GGGGTTTTAAAA	IIIIIIIIIIII	RG:Z:p
      """

  Scenario: A templated read fails QC when any record of its fragment does
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      both	589	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      both	653	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      last	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      last	653	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ
      none	77	*	0	0	*	*	0	0	GGGGAAAA	IIIIIIII
      none	141	*	0	0	*	*	0	0	TTTTCCCC	JJJJJJJJ
      """
    When I run `unmux in.sam --pool p --extract r1=0:0:end --extract umi=1:0:4 --template r1 --tag RX=umi --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      both	516	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RX:Z:TTTT	QX:Z:JJJJ	RG:Z:p
      last	516	*	0	255	*	*	0	0	CCCCAAAA	IIIIIIII	RX:Z:GGGG	QX:Z:JJJJ	RG:Z:p
      none	4	*	0	255	*	*	0	0	GGGGAAAA	IIIIIIII	RX:Z:TTTT	QX:Z:JJJJ	RG:Z:p
      """

  Scenario: --per-record reads each mate of a paired SAM input on its own
    Given a file "paired.sam" containing:
      """
      @HD	VN:1.6
      read1	77	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      read1	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      """
    When I run `unmux paired.sam --out paired.R%ordinal.fq`
    Then the exit code is 0
    And stderr contains "1 reads (0 assigned, 1 pass-through"
    And the file "paired.R2.fq" contains "TTTTGGGG"
    When I run `unmux paired.sam --pool p --per-record --out alone.sam`
    Then the exit code is 0
    And stderr contains "2 reads (0 assigned, 2 pass-through"
    And the records of "alone.sam" are:
      """
      read1	4	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:p
      read1	4	*	0	255	*	*	0	0	TTTTGGGG	JJJJJJJJ	RG:Z:p
      """
    When I run `unmux paired.sam --per-record --out alone.R%ordinal.fq`
    Then the exit code is 0
    And a file "alone.R2.fq" does not exist

  Scenario: Paired SAM records with their mates apart are refused unless --per-record reads each alone
    Given a file "apart.sam" containing:
      """
      @HD	VN:1.6
      read1	77	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      read2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      read1	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      read2	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ
      """
    When I run `unmux apart.sam --pool p --out out.sam`
    Then the exit code is 1
    And stderr contains "paired record `read1` is followed by `read2`, not its mate"
    And stderr contains "samtools sort -n"
    And stderr contains "--per-record"
    And a file "out.sam" does not exist
    When I run `unmux apart.sam --pool p --per-record --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      read1	4	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:p
      read2	4	*	0	255	*	*	0	0	CCCCAAAA	IIIIIIII	RG:Z:p
      read1	4	*	0	255	*	*	0	0	TTTTGGGG	JJJJJJJJ	RG:Z:p
      read2	4	*	0	255	*	*	0	0	GGGGTTTT	JJJJJJJJ	RG:Z:p
      """

  Scenario: A SAM input that mixes paired and unpaired records is refused
    Given a file "mixed.sam" containing:
      """
      @HD	VN:1.6
      read1	77	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      read1	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      read2	4	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      """
    When I run `unmux mixed.sam --out out.sam`
    Then the exit code is 1
    And stderr contains "input mixes paired and unpaired records (`read2` is unpaired, unlike the first record)"
    And stderr contains "--per-record"
    And a file "out.sam" does not exist

  Scenario Outline: An empty mate of a <format> pair keeps the other mate's qualities in a template
    Given a file "trimmed.sam" containing:
      """
      @HD	VN:1.6
      a	77	*	0	0	*	*	0	0	ACGTACGT	ABCDEFGH
      a	141	*	0	0	*	*	0	0	*	*
      """
    # The second mate was trimmed to nothing, so its SEQ and QUAL are *. Like an empty FASTQ read, it
    # has an empty quality string, so the template keeps the first mate's qualities and RX gets an
    # empty QX.
    When I run `unmux trimmed.sam --pool p --out in.<format>`
    Then the exit code is 0
    When I run `unmux in.<format> --pool p --extract umi=1:0:2 --extract r=0:0:8 --template r+umi --tag RX=umi --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      a	4	*	0	255	*	*	0	0	ACGTACGT	ABCDEFGH	RX:Z:	QX:Z:	RG:Z:p
      """

    Examples:
      | format |
      | sam    |
      | bam    |
      | cram   |

  Scenario: A mate with bases and QUAL * still has no qualities
    Given a file "unscored.sam" containing:
      """
      @HD	VN:1.6
      a	77	*	0	0	*	*	0	0	ACGTACGT	ABCDEFGH
      a	141	*	0	0	*	*	0	0	TTGG	*
      """
    # Only an empty read has an empty quality string. A read with bases and QUAL * has no qualities,
    # so a template that joins it has none either, and it cannot be written as FASTQ.
    When I run `unmux unscored.sam --pool p --extract umi=1:0:2 --extract r=0:0:8 --template r+umi --tag RX=umi --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      a	4	*	0	255	*	*	0	0	ACGTACGTTT	*	RX:Z:TT	RG:Z:p
      """
    When I run `unmux unscored.sam --pool p --out out.R%ordinal.fq`
    Then the exit code is 1
    And stderr contains "cannot write FASTQ output for a record with no qualities (a FASTA input, or QUAL `*` in SAM/BAM/CRAM)"

  Scenario: Two SAM inputs that each hold one mate of every pair are read as pairs
    Given a file "r1.sam" containing:
      """
      @HD	VN:1.6
      q1	77	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      q2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      """
    And a file "r2.sam" containing:
      """
      @HD	VN:1.6
      q1	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      q2	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ
      """
    When I run `unmux r1.sam r2.sam --pool p --out pairs.sam`
    Then the exit code is 0
    And the records of "pairs.sam" are:
      """
      q1	77	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:p
      q1	141	*	0	255	*	*	0	0	TTTTGGGG	JJJJJJJJ	RG:Z:p
      q2	77	*	0	255	*	*	0	0	CCCCAAAA	IIIIIIII	RG:Z:p
      q2	141	*	0	255	*	*	0	0	GGGGTTTT	JJJJJJJJ	RG:Z:p
      """

  Scenario: Two paired SAM inputs are refused unless --per-record reads one record from each
    Given a file "reads.sam" containing:
      """
      @HD	VN:1.6
      q1	77	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      q1	141	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      q2	77	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      q2	141	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ
      """
    And a file "idx.sam" containing:
      """
      @HD	VN:1.6
      q1	77	*	0	0	*	*	0	0	AAAA	IIII
      q1	141	*	0	0	*	*	0	0	CCCC	JJJJ
      q2	77	*	0	0	*	*	0	0	GGGG	IIII
      q2	141	*	0	0	*	*	0	0	TTTT	JJJJ
      """
    # A read takes one record from each input, so the first mates of both files would make one read
    # and the last mates another, splitting the mates of every template across two reads.
    When I run `unmux reads.sam idx.sam --pool p --out two.sam`
    Then the exit code is 1
    And stderr contains "input file 0 (reads.sam) holds both first (0x40) and last (0x80) mates of paired records, as of `q1`"
    And stderr contains "would split its mates across reads"
    And stderr contains "samtools fastq -1/-2"
    And stderr contains "--per-record"
    And a file "two.sam" does not exist
    When I run `unmux reads.sam idx.sam --pool p --per-record --out two.sam`
    Then the exit code is 0
    And the records of "two.sam" are:
      """
      q1	77	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:p
      q1	141	*	0	255	*	*	0	0	AAAA	IIII	RG:Z:p
      q1	77	*	0	255	*	*	0	0	TTTTGGGG	JJJJJJJJ	RG:Z:p
      q1	141	*	0	255	*	*	0	0	CCCC	JJJJ	RG:Z:p
      q2	77	*	0	255	*	*	0	0	CCCCAAAA	IIIIIIII	RG:Z:p
      q2	141	*	0	255	*	*	0	0	GGGG	IIII	RG:Z:p
      q2	77	*	0	255	*	*	0	0	GGGGTTTT	JJJJJJJJ	RG:Z:p
      q2	141	*	0	255	*	*	0	0	TTTT	JJJJ	RG:Z:p
      """

  Scenario: Two interleaved FASTQ inputs are refused unless --per-record reads one record from each
    Given a file "reads.fq" containing:
      """
      @q1/1
      ACGTACGT
      +
      IIIIIIII
      @q1/2
      TTTTGGGG
      +
      JJJJJJJJ
      """
    And a file "idx.fq" containing:
      """
      @q1/1
      AAAA
      +
      IIII
      @q1/2
      CCCC
      +
      JJJJ
      """
    When I run `unmux reads.fq idx.fq --out out.R%ordinal.fq`
    Then the exit code is 1
    And stderr contains "input file 0 (reads.fq) is interleaved (`q1/1` is followed by its mate `q1/2`)"
    And stderr contains "would split its mates across reads"
    And stderr contains "--per-record"
    And a file "out.R1.fq" does not exist
    And a file "out.R2.fq" does not exist
    When I run `unmux reads.fq idx.fq --per-record --out out.R%ordinal.fq`
    Then the exit code is 0
    And the file "out.R1.fq" is:
      """
      @q1/1
      ACGTACGT
      +
      IIIIIIII
      @q1/2
      TTTTGGGG
      +
      JJJJJJJJ
      """
    And the file "out.R2.fq" is:
      """
      @q1/1
      AAAA
      +
      IIII
      @q1/2
      CCCC
      +
      JJJJ
      """

  Scenario Outline: Mates that follow unpaired records in one of several inputs are refused on <threads> threads
    Given a file "reads.sam" containing:
      """
      @HD	VN:1.6
      q1	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      q2	77	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      q2	141	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      """
    And a file "idx.sam" containing:
      """
      @HD	VN:1.6
      q1	4	*	0	0	*	*	0	0	AAAA	IIII
      q2	77	*	0	0	*	*	0	0	CCCC	JJJJ
      q2	141	*	0	0	*	*	0	0	GGGG	IIII
      """
    And a file "reads.fq" containing:
      """
      @q1
      ACGTACGT
      +
      IIIIIIII
      @q2/1
      TTTTGGGG
      +
      JJJJJJJJ
      @q2/2
      CCCCAAAA
      +
      IIIIIIII
      """
    And a file "idx.fq" containing:
      """
      @q1
      AAAA
      +
      IIII
      @q2/1
      CCCC
      +
      JJJJ
      @q2/2
      GGGG
      +
      IIII
      """
    When I run `unmux reads.sam idx.sam --pool p --threads <threads> --out two.sam`
    Then the exit code is 1
    And stderr contains "input file 0 (reads.sam) holds both first (0x40) and last (0x80) mates of paired records, as of `q2`"
    And a file "two.sam" does not exist
    When I run `unmux reads.fq idx.fq --threads <threads> --out out.R%ordinal.fq`
    Then the exit code is 1
    And stderr contains "input file 0 (reads.fq) is interleaved (`q2/1` is followed by its mate `q2/2`)"
    And a file "out.R1.fq" does not exist

    Examples:
      | threads |
      | 1       |
      | 2       |
      | 8       |

  Scenario Outline: Single-end FASTQ inputs are still read one record from each on <threads> threads
    Given a file "r1.fq" containing:
      """
      @read1/1
      ACGTACGT
      +
      IIIIIIII
      @read2/1
      CCCCAAAA
      +
      IIIIIIII
      """
    And a file "i1.fq" containing:
      """
      @read1
      GGGG
      +
      IIII
      @read2
      TTTT
      +
      IIII
      """
    And a file "r2.fq" containing:
      """
      @read1/2
      TTTTGGGG
      +
      JJJJJJJJ
      @read2/2
      GGGGTTTT
      +
      JJJJJJJJ
      """
    And a file "other.fq" containing:
      """
      @read1
      GGGG
      +
      IIII
      @read3
      TTTT
      +
      IIII
      """
    And a file "short.fq" containing:
      """
      @read1
      GGGG
      +
      IIII
      """
    # A /1 or /2 suffix alone does not make one of several inputs interleaved; only adjacent mates do.
    When I run `unmux r1.fq i1.fq r2.fq --threads <threads> --out out.R%ordinal.fq`
    Then the exit code is 0
    And the file "out.R1.fq" is:
      """
      @read1/1
      ACGTACGT
      +
      IIIIIIII
      @read2/1
      CCCCAAAA
      +
      IIIIIIII
      """
    And the file "out.R2.fq" is:
      """
      @read1
      GGGG
      +
      IIII
      @read2
      TTTT
      +
      IIII
      """
    And the file "out.R3.fq" is:
      """
      @read1/2
      TTTTGGGG
      +
      JJJJJJJJ
      @read2/2
      GGGGTTTT
      +
      JJJJJJJJ
      """
    When I run `unmux r1.fq other.fq --threads <threads> --out other.R%ordinal.fq`
    Then the exit code is 1
    And stderr contains "input segments have mismatched read names (`read2/1` vs `read3`)"
    And a file "other.R1.fq" does not exist
    When I run `unmux r1.fq short.fq --threads <threads> --out short.R%ordinal.fq`
    Then the exit code is 1
    And stderr contains "input files have unequal record counts"
    And a file "short.R1.fq" does not exist

    Examples:
      | threads |
      | 1       |
      | 2       |
      | 8       |

  Scenario Outline: Single-end <format> inputs are still read one record from each on <threads> threads
    Given a file "r1.sam" containing:
      """
      @HD	VN:1.6
      read1	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      read2	4	*	0	0	*	*	0	0	CCCCAAAA	IIIIIIII
      """
    And a file "r2.sam" containing:
      """
      @HD	VN:1.6
      read1	4	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ
      read2	4	*	0	0	*	*	0	0	GGGGTTTT	JJJJJJJJ
      """
    When I run `unmux r1.sam --out in1.<format>`
    Then the exit code is 0
    When I run `unmux r2.sam --out in2.<format>`
    Then the exit code is 0
    When I run `unmux in1.<format> in2.<format> --pool p --threads <threads> --out pairs.sam`
    Then the exit code is 0
    And the records of "pairs.sam" are:
      """
      read1	77	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:p
      read1	141	*	0	255	*	*	0	0	TTTTGGGG	JJJJJJJJ	RG:Z:p
      read2	77	*	0	255	*	*	0	0	CCCCAAAA	IIIIIIII	RG:Z:p
      read2	141	*	0	255	*	*	0	0	GGGGTTTT	JJJJJJJJ	RG:Z:p
      """

    Examples:
      | format | threads |
      | sam    | 1       |
      | sam    | 8       |
      | bam    | 2       |
      | bam    | 16      |
      | cram   | 1       |
      | cram   | 8       |

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

  Scenario Outline: FASTA converts to CRAM, which reads back with no qualities
    Given a file "reads.fa" containing:
      """
      >r1
      ACGTACGTAACCGGTTACGT
      >r2 XY:Z:hello
      TTTTGGGG
      """
    When I run `unmux reads.fa --pool p --threads <threads> --out out.cram`
    Then the exit code is 0
    When I run `unmux out.cram --pool p --out back.sam`
    Then the exit code is 0
    And the records of "back.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGTAACCGGTTACGT	*	RG:Z:p
      r2	4	*	0	255	*	*	0	0	TTTTGGGG	*	XY:Z:hello	RG:Z:p
      """
    When I run `unmux out.cram --out back.fa`
    Then the exit code is 0
    And the file "back.fa" is:
      """
      >r1 RG:Z:p
      ACGTACGTAACCGGTTACGT
      >r2 XY:Z:hello	RG:Z:p
      TTTTGGGG
      """

    Examples:
      | threads |
      | 1       |
      | 2       |
      | 8       |

  Scenario Outline: SAM and BAM records with no qualities convert to CRAM and back
    Given a file "unmapped.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGTAACCGGTTACGT	*	XY:Z:hello
      r2	4	*	0	0	*	*	0	0	TTTTGGGG	*
      """
    When I run `unmux unmapped.sam --out in.<ext>`
    Then the exit code is 0
    When I run `unmux in.<ext> --out out.cram`
    Then the exit code is 0
    When I run `unmux out.cram --pool p --out back.sam`
    Then the exit code is 0
    And the records of "back.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGTAACCGGTTACGT	*	XY:Z:hello	RG:Z:p
      r2	4	*	0	255	*	*	0	0	TTTTGGGG	*	RG:Z:p
      """

    Examples:
      | ext |
      | sam |
      | bam |

  Scenario: A CRAM output holds records with and without qualities
    Given a file "mixed.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGT	ABCDEFGH	XY:Z:hello
      r2	4	*	0	0	*	*	0	0	GGGGCCCC	*
      r3	4	*	0	0	*	*	0	0	TTTTAAAA	IIIIJJJJ
      r4	4	*	0	0	*	*	0	0	CC	*	XZ:i:7
      """
    When I run `unmux mixed.sam --out out.cram`
    Then the exit code is 0
    When I run `unmux out.cram --pool p --out back.sam`
    Then the exit code is 0
    And the records of "back.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGT	ABCDEFGH	XY:Z:hello	RG:Z:p
      r2	4	*	0	255	*	*	0	0	GGGGCCCC	*	RG:Z:p
      r3	4	*	0	255	*	*	0	0	TTTTAAAA	IIIIJJJJ	RG:Z:p
      r4	4	*	0	255	*	*	0	0	CC	*	XZ:i:7	RG:Z:p
      """

  Scenario Outline: A SAM record without a name keeps its place in <format> output
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      n1	4	*	0	0	*	*	0	0	AAAA	IIII
      *	4	*	0	0	*	*	0	0	CCCC	IIII
      n3	4	*	0	0	*	*	0	0	GGGG	IIII
      n4	4	*	0	0	*	*	0	0	TTTT	IIII
      """
    When I run `unmux in.sam --pool p --out out.<format>`
    Then the exit code is 0
    When I run `unmux out.<format> --pool p --out back.sam`
    Then the exit code is 0
    And the records of "back.sam" are:
      """
      n1	4	*	0	255	*	*	0	0	AAAA	IIII	RG:Z:p
      *	4	*	0	255	*	*	0	0	CCCC	IIII	RG:Z:p
      n3	4	*	0	255	*	*	0	0	GGGG	IIII	RG:Z:p
      n4	4	*	0	255	*	*	0	0	TTTT	IIII	RG:Z:p
      """

    Examples:
      | format |
      | sam    |
      | bam    |
      | cram   |

  Scenario: A FASTQ record with an empty name keeps its place in CRAM output
    Given a file "in.fq" containing:
      """
      @n1
      AAAA
      +
      IIII
      @
      CCCC
      +
      IIII
      @n3
      GGGG
      +
      IIII
      @n4
      TTTT
      +
      IIII
      """
    When I run `unmux in.fq --pool p --out out.cram`
    Then the exit code is 0
    When I run `unmux out.cram --out back.fq`
    Then the exit code is 0
    And the file "back.fq" is:
      """
      @n1 RG:Z:p
      AAAA
      +
      IIII
      @ RG:Z:p
      CCCC
      +
      IIII
      @n3 RG:Z:p
      GGGG
      +
      IIII
      @n4 RG:Z:p
      TTTT
      +
      IIII
      """

  Scenario: A record without a name in an htslib CRAM reads back without one
    Given the test file "nameless.htslib-3.1.cram"
    # samtools 1.25 wrote this CRAM 3.1 from nameless.sam, storing the missing name as * and a NUL.
    When I run `unmux nameless.htslib-3.1.cram --pool p --out back.sam`
    Then the exit code is 0
    And the records of "back.sam" are:
      """
      n1	4	*	0	255	*	*	0	0	AAAA	IIII	RG:Z:p
      *	4	*	0	255	*	*	0	0	CCCC	IIII	RG:Z:p
      n3	4	*	0	255	*	*	0	0	GGGG	IIII	RG:Z:p
      n4	4	*	0	255	*	*	0	0	TTTT	IIII	RG:Z:p
      """
    When I run `unmux nameless.htslib-3.1.cram --out back.fq`
    Then the exit code is 0
    And the file "back.fq" is:
      """
      @n1
      AAAA
      +
      IIII
      @
      CCCC
      +
      IIII
      @n3
      GGGG
      +
      IIII
      @n4
      TTTT
      +
      IIII
      """

  Scenario Outline: A read name too long for BAM fails <format> output with the same error
    Given a file "in.fq" containing:
      """
      @rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
      ACGT
      +
      IIII
      """
    When I run `unmux in.fq --out out.<format>`
    Then the exit code is 1
    And stderr contains "failed to write SAM/BAM/CRAM record: invalid name: invalid length: expected 1 <= n <= 254, got 256"
    And a file "out.<format>" does not exist

    Examples:
      | format |
      | bam    |
      | cram   |

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

  Scenario Outline: An empty read in <format> converts to FASTQ as an empty record and back
    Given a file "trimmed.sam" containing:
      """
      @HD	VN:1.6
      trimmed	4	*	0	0	*	*	0	0	*	*
      kept	4	*	0	0	*	*	0	0	ACGT	ABCD
      """
    # A read trimmed to nothing is written as samtools fastq writes it: its name, an empty sequence, a
    # +, and an empty quality line. Read back, it is empty again, with no qualities filled in.
    When I run `unmux trimmed.sam --pool p --out in.<format>`
    Then the exit code is 0
    When I run `unmux in.<format> --pool p --out out.fq`
    Then the exit code is 0
    And the file "out.fq" is:
      """
      @trimmed RG:Z:p

      +

      @kept RG:Z:p
      ACGT
      +
      ABCD
      """
    When I run `unmux out.fq --pool p --out back.<format>`
    Then the exit code is 0
    When I run `unmux back.<format> --pool p --out final.sam`
    Then the exit code is 0
    And the records of "final.sam" are:
      """
      trimmed	4	*	0	255	*	*	0	0	*	*	RG:Z:p
      kept	4	*	0	255	*	*	0	0	ACGT	ABCD	RG:Z:p
      """

    Examples:
      | format |
      | sam    |
      | bam    |
      | cram   |

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

  Scenario: BAM output stores U and u bases as T, while SAM and FASTQ output keep them
    Given a file "rna.fq" containing:
      """
      @r1
      ACGUacgu
      +
      IIIIIIII
      """
    # BAM stores each base as a 4-bit code that has no U.
    When I run `unmux rna.fq --pool p --out out.bam`
    Then the exit code is 0
    And the records of "out.bam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGT	IIIIIIII	RG:Z:p
      """
    When I run `unmux rna.fq --pool p --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGUacgu	IIIIIIII	RG:Z:p
      """
    When I run `unmux rna.fq --out out.fq`
    Then the exit code is 0
    And the file "out.fq" is:
      """
      @r1
      ACGUacgu
      +
      IIIIIIII
      """

  Scenario: Sample, unassigned, and removed BAM outputs store U bases as T
    Given a file "rna.fq" containing:
      """
      @a
      AAAAGUCC
      +
      IIIIIIII
      @b
      CCCCGUCC
      +
      IIIIIIII
      @c
      GGGGGUCC
      +
      IIIIIIII
      """
    When I run `unmux rna.fq --group bc={x=AAAA,y=GGGG} --group bc::loc=0:0:4 --sample s=bc::x --remove bc::y=rm.bam --extract body=0:4:end --template body --out %sample.bam --unassigned un.bam`
    Then the exit code is 0
    And the records of "s.bam" are:
      """
      a	4	*	0	255	*	*	0	0	GTCC	IIII	RG:Z:s
      """
    And the records of "un.bam" are:
      """
      b	4	*	0	255	*	*	0	0	CCCCGTCC	IIIIIIII	RG:Z:rna
      """
    And the records of "rm.bam" are:
      """
      c	4	*	0	255	*	*	0	0	GGGGGTCC	IIIIIIII	RG:Z:rna
      """

  Scenario: U and u bases of SAM records, aligned or not, are stored as T in BAM output
    Given a file "rna.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:100
      q1	0	chr1	10	60	8M	*	0	0	ACGUacgu	ABCDEFGH
      q2	4	*	0	0	*	*	0	0	UUGGCCAA	ABCDEFGH
      """
    When I run `unmux rna.sam --pool p --out out.bam`
    Then the exit code is 0
    And the records of "out.bam" are:
      """
      q1	4	*	0	255	*	*	0	0	ACGTACGT	ABCDEFGH	RG:Z:p
      q2	4	*	0	255	*	*	0	0	TTGGCCAA	ABCDEFGH	RG:Z:p
      """
