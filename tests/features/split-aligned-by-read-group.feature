Feature: Split an aligned input by read group, keeping alignments
  An `@RG` split of an aligned input (one whose header carries `@SQ`) into
  SAM/BAM keeps every record whole, secondary and supplementary included:
  flags, position, CIGAR, MAPQ, mate fields, and every tag. Each output keeps
  the input header (`@HD` with its sort order, `@SQ`, `@PG`, `@CO`) with only
  its own `@RG` lines and an unmux `@PG`, and records keep their input order.
  A `--template`, a `--tag`, a sequence-matching group, aligned CRAM output,
  FASTX outputs beside SAM/BAM ones, and several inputs are refused. For
  unmapped records, run `samtools reset` on the outputs.

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

  Scenario: metrics count primary records only
    When I run `unmux in.sam --pool p1 --group sm=@RG::SM --sample-from-group sm --out out/%sample.bam --metrics-per-sample p.tsv --metrics-summary s.tsv`
    Then the exit code is 0
    And the file "p.tsv" is:
      """
      pool	sample	sub_sample	reads	frac_of_pool	frac_bases_unextracted
      p1	sampleA		4	0.571429	0.000000
      p1	sampleB		3	0.428571	0.000000
      """
    And the file "s.tsv" contains "reads_assigned"
    And stderr contains "2 secondary/supplementary records passed through"

  Scenario: a --template on an aligned split into BAM is refused
    When I run `unmux in.sam --group sm=@RG::SM --sample-from-group sm --extract head=0:0:4 --template head --out out/%sample.bam`
    Then the exit code is 1
    And stderr contains "--template would rewrite"
    And stderr contains "samtools reset"
    And a file "out/sampleA.bam" does not exist

  Scenario: a --template on an aligned split into FASTQ still works
    # The mates of this coordinate-sorted input are apart, so each record is read on its own.
    When I run `unmux in.sam --per-record --group sm=@RG::SM --sample-from-group sm --extract head=0:0:4 --template head --out out/%sample.fq`
    Then the exit code is 0
    And the file "out/sampleA.fq" contains "ACGT"

  Scenario: an aligned split into CRAM is refused before any file is written
    When I run `unmux in.sam --group sm=@RG::SM --sample-from-group sm --out out/%sample.cram`
    Then the exit code is 1
    And stderr contains "writing aligned records to CRAM is not supported"
    And stderr contains "samtools view -C"
    And a file "out/sampleA.cram" does not exist
    And a file "out/sampleB.cram" does not exist

  Scenario: several aligned inputs are refused, pointing at samtools merge
    When I run `unmux --in 0=in.sam --in 1=in.sam --group sm=@RG::SM --sample-from-group sm --out out/%sample.bam`
    Then the exit code is 1
    And stderr contains "single input"
    And stderr contains "samtools merge"

  Scenario Outline: an htslib CRAM <version> input with several containers is split from a file
    Given the test file "aligned.htslib-<version>.cram"
    And the test file "ref.fa"
    And the test file "ref.fa.fai"
    When I run `unmux aligned.htslib-<version>.cram --reference ref.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 0
    And the records of "out/rg1.sam" are the "RG:Z:rg1" records of the test file "aligned.htslib.sam"
    And the records of "out/rg2.sam" are the "RG:Z:rg2" records of the test file "aligned.htslib.sam"
    And the file "out/rg1.sam" contains "SN:chr1"

    Examples:
      | version |
      | 3.0     |
      | 3.1     |

  Scenario Outline: an htslib CRAM <version> input is split from stdin
    Given the test file "aligned.htslib-<version>.cram"
    And the test file "ref.fa"
    And the test file "ref.fa.fai"
    When I run `unmux - --reference ref.fa --group rg=@RG --sample-from-group rg --out out/%sample.bam` with "aligned.htslib-<version>.cram" on stdin
    Then the exit code is 0
    And the records of "out/rg2.bam" are the "RG:Z:rg2" records of the test file "aligned.htslib.sam"

    Examples:
      | version |
      | 3.0     |
      | 3.1     |

  Scenario Outline: a CRAM split matches htslib's decode, MD and NM tags included, from <file>
    Given the test file "<file>"
    And the test file "contigs.fa"
    And the test file "contigs.fa.fai"
    When I run `unmux <file> --reference contigs.fa --group sm=@RG::SM --sample-from-group sm --out out/%sample.sam --unassigned out/unassigned.sam`
    Then the exit code is 0
    And the records of "out/sA.sam" are the "RG:Z:rg1" records of the test file "contigs.htslib.sam"
    And the records of "out/sB.sam" are the "RG:Z:rg2" records of the test file "contigs.htslib.sam"
    And the file "out/unassigned.sam" contains "RG:Z:ghost	MD:Z:20	NM:i:0"

    Examples:
      | file                   |
      | contigs.cram           |
      | contigs.multi-ref.cram |

  Scenario Outline: a CRAM record running past the end of its reference sequence decodes as htslib does on <threads> threads
    Given the test file "past-end.cram"
    And the test file "contigs.fa"
    And the test file "contigs.fa.fai"
    When I run `unmux past-end.cram --reference contigs.fa --group sm=@RG::SM --sample-from-group sm --out out/%sample.sam --threads <threads>`
    Then the exit code is 0
    And the records of "out/sA.sam" are:
      """
      pastEnd	0	ctgB	196	60	10M	*	0	0	TGCGCACGTA	IIIIIIIIII	MD:Z:5	NM:i:0	RG:Z:rg1
      insBeforeEnd	0	ctgB	193	60	3M2I8M	*	0	0	TTGGGTGCGCCAG	IIIIIIIIIIIII	MD:Z:8	NM:i:2	RG:Z:rg1
      mismatchThenEnd	0	ctgB	195	60	8M	*	0	0	GAGCGCGG	IIIIIIII	MD:Z:1T4	NM:i:1	RG:Z:rg1
      skipAcrossEnd	0	ctgB	197	60	2M3N4M	*	0	0	GCACGT	IIIIII	MD:Z:2	NM:i:0	RG:Z:rg1
      """
    And the records of "out/sB.sam" are:
      """
      clipAcrossEnd	0	ctgB	199	60	2S6M	*	0	0	AAGCTTCA	IIIIIIII	MD:Z:2	NM:i:0	RG:Z:rg2
      """

    Examples:
      | threads |
      | 1       |
      | 8       |
      | 32      |

  Scenario Outline: a CRAM record running past the end of an embedded reference decodes as htslib does on <threads> threads
    Given the test file "past-end.embedded-ref.cram"
    When I run `unmux past-end.embedded-ref.cram --group sm=@RG::SM --sample-from-group sm --out out/%sample.sam --threads <threads>`
    Then the exit code is 0
    And the records of "out/sA.sam" are:
      """
      pastEnd	0	ctgB	196	60	10M	*	0	0	TGCGCACGTA	IIIIIIIIII	MD:Z:5	NM:i:0	RG:Z:rg1
      insBeforeEnd	0	ctgB	193	60	3M2I8M	*	0	0	TTGGGTGCGCCAG	IIIIIIIIIIIII	MD:Z:8	NM:i:2	RG:Z:rg1
      mismatchThenEnd	0	ctgB	195	60	8M	*	0	0	GAGCGCGG	IIIIIIII	MD:Z:1T4	NM:i:1	RG:Z:rg1
      skipAcrossEnd	0	ctgB	197	60	2M3N4M	*	0	0	GCACGT	IIIIII	MD:Z:2	NM:i:0	RG:Z:rg1
      """
    And the records of "out/sB.sam" are:
      """
      clipAcrossEnd	0	ctgB	199	60	2S6M	*	0	0	AAGCTTCA	IIIIIIII	MD:Z:2	NM:i:0	RG:Z:rg2
      """
    And stderr does not contain "panicked"

    Examples:
      | threads |
      | 1       |
      | 5       |
      | 8       |
      | 32      |

  Scenario Outline: MD and NM are regenerated only for slices that use the external reference on <threads> threads
    Given the test file "contigs.mixed-ref.cram"
    And the test file "contigs.fa"
    And the test file "contigs.fa.fai"
    When I run `unmux contigs.mixed-ref.cram --reference contigs.fa --group sm=@RG::SM --sample-from-group sm --out out/%sample.sam --unassigned out/unassigned.sam --threads <threads>`
    Then the exit code is 0
    And the records of "out/sA.sam" are the "RG:Z:rg1" records of the test file "contigs.mixed-ref.htslib.sam"
    And the records of "out/sB.sam" are the "RG:Z:rg2" records of the test file "contigs.mixed-ref.htslib.sam"

    Examples:
      | threads |
      | 1       |
      | 5       |
      | 8       |
      | 16      |
      | 32      |

  Scenario: a --reference sequence shorter than its @SQ length fails cleanly without the check
    Given the test file "contigs.cram"
    And the test file "contigs.short.fa"
    And the test file "contigs.short.fa.fai"
    When I run `unmux contigs.cram --reference contigs.short.fa --skip-reference-check --group sm=@RG::SM --sample-from-group sm --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "reference sequence `ctgA` is 300 bases in the --reference FASTA but 400 in the @SQ header of CRAM input contigs.cram"
    And stderr does not contain "panicked"
    And a file "out/sA.sam" does not exist

  Scenario: a reference-compressed CRAM input needs --reference
    Given the test file "aligned.htslib-3.1.cram"
    When I run `unmux aligned.htslib-3.1.cram --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "CRAM input aligned.htslib-3.1.cram needs reference sequence `chr1` to decode"
    And stderr contains "--reference"
    And a file "out/rg1.sam" does not exist
    And a file "out/rg2.sam" does not exist

  Scenario: a failed run removes the outputs it created and leaves files that were already there unchanged
    Given the test file "aligned.htslib-3.0.cram"
    And a file "out/rg1.sam" containing:
      """
      existing
      """
    When I run `unmux aligned.htslib-3.0.cram --group rg=@RG --sample-from-group rg --out out/%sample.sam --unassigned out/unassigned.sam --metrics-summary out/summary.tsv`
    Then the exit code is 1
    And the file "out/rg1.sam" is:
      """
      existing
      """
    And the directory "out" holds exactly "rg1.sam"

  Scenario: a CRAM with an embedded reference splits as htslib decodes it, MD and NM included, without --reference
    Given the test file "contigs.embedded-ref.cram"
    When I run `unmux contigs.embedded-ref.cram --group sm=@RG::SM --sample-from-group sm --out out/%sample.sam --unassigned out/unassigned.sam`
    Then the exit code is 0
    And the records of "out/sA.sam" are the "RG:Z:rg1" records of the test file "contigs.htslib.sam"
    And the records of "out/sB.sam" are the "RG:Z:rg2" records of the test file "contigs.htslib.sam"
    And the file "out/unassigned.sam" contains "RG:Z:ghost	MD:Z:20	NM:i:0"

  Scenario: a CRAM written without a reference splits without --reference
    Given the test file "contigs.no-ref.cram"
    When I run `unmux contigs.no-ref.cram --group sm=@RG::SM --sample-from-group sm --out out/%sample.sam`
    Then the exit code is 0
    And the records of "out/sA.sam" are the "RG:Z:rg1" records of the test file "contigs.sam"
    And the records of "out/sB.sam" are the "RG:Z:rg2" records of the test file "contigs.sam"

  Scenario Outline: a CRAM written without a reference in multi-reference slices splits <how>
    Given the test file "contigs.no-ref.multi-ref.cram"
    And the test file "contigs.fa"
    And the test file "contigs.fa.fai"
    When I run `unmux contigs.no-ref.multi-ref.cram <args> --group sm=@RG::SM --sample-from-group sm --out out/%sample.sam`
    Then the exit code is 0
    And the records of "out/sA.sam" are the "RG:Z:rg1" records of the test file "contigs.sam"
    And the records of "out/sB.sam" are the "RG:Z:rg2" records of the test file "contigs.sam"

    Examples:
      | how                  | args                   |
      | without --reference  | --threads 1            |
      | with --reference     | --reference contigs.fa |

  Scenario: --reference without a CRAM input is unused and says so
    Given the test file "ref.fa"
    And the test file "ref.fa.fai"
    When I run `unmux in.sam --reference ref.fa --group sm=@RG::SM --sample-from-group sm --out out/%sample.bam`
    Then the exit code is 0
    And stderr contains "--reference is unused: no input is CRAM"

  Scenario: no drop warning when every primary record has a known read group
    When I run `unmux in.sam --group sm=@RG::SM --sample-from-group sm --out out/%sample.bam`
    Then the exit code is 0
    And stderr does not contain "were dropped"

  Scenario: a dropped record with an unknown read group is counted in a warning
    Given a file "orphan.sam" containing:
      """
      @HD	VN:1.6	SO:coordinate
      @SQ	SN:chr1	LN:1000
      @RG	ID:rg1	SM:sampleA
      r1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	IIIIIIII	RG:Z:rg1
      r2	0	chr1	300	50	8M	*	0	0	CCCCCCCC	IIIIIIII	RG:Z:ghost
      """
    When I run `unmux orphan.sam --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 0
    And stderr contains "1 records carried no RG:Z, or one the header does not declare, and were dropped"
