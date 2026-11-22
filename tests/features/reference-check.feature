Feature: --reference is checked against the CRAM @SQ lines
  Before any output is written, every @SQ contig of a CRAM input must be in
  the --reference FASTA with the same LN and, where @SQ has one, the same M5,
  so a wrong FASTA cannot decode to wrong bases.

  Background:
    Given the test file "ref.fa"
    And the test file "ref.fa.fai"

  Scenario Outline: a matching reference passes the check for CRAM <version>
    Given the test file "aligned.htslib-<version>.cram"
    When I run `unmux aligned.htslib-<version>.cram --reference ref.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 0
    And a file "out/rg1.sam" exists
    And stderr does not contain "have no M5"

    Examples:
      | version |
      | 3.0     |
      | 3.1     |

  Scenario: a contig with the wrong length fails before any output exists
    Given the test file "aligned.htslib-3.1.cram"
    And a copy of the test file "ref.fa" named "short.fa"
    And a copy of the test file "ref.fa.fai" named "short.fa.fai" with "3000" replaced by "2990"
    When I run `unmux aligned.htslib-3.1.cram --reference short.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "contig `chr1` is 2990 bp in the --reference FASTA short.fa but @SQ LN:3000 in aligned.htslib-3.1.cram"
    And a file "out" does not exist

  Scenario: a contig with the wrong MD5 fails before any output exists
    Given the test file "aligned.htslib-3.1.cram"
    And a copy of the test file "ref.fa" named "snv.fa" with "TTTCCTCATG" replaced by "GTTCCTCATG"
    And a copy of the test file "ref.fa.fai" named "snv.fa.fai"
    When I run `unmux aligned.htslib-3.1.cram --reference snv.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "contig `chr1` has MD5"
    And stderr contains "but @SQ M5:e1a82f39b96c07c57c14cfb6ec9e0caa in aligned.htslib-3.1.cram"
    And a file "out" does not exist

  Scenario: a contig missing from the reference fails before any output exists
    Given the test file "aligned.htslib-3.1.cram"
    And a copy of the test file "ref.fa" named "renamed.fa" with ">chr1" replaced by ">chrX"
    And a copy of the test file "ref.fa.fai" named "renamed.fa.fai" with "chr1" replaced by "chrX"
    When I run `unmux aligned.htslib-3.1.cram --reference renamed.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "contig `chr1` (@SQ in aligned.htslib-3.1.cram) is missing from the --reference FASTA renamed.fa"
    And stderr contains "--skip-reference-check"
    And a file "out" does not exist

  Scenario: a FASTA whose header line changed after indexing fails as a stale .fai
    Given the test file "aligned.htslib-3.1.cram"
    And a copy of the test file "ref.fa" named "edited.fa" with ">chr1" replaced by ">chr1 assembled"
    And a copy of the test file "ref.fa.fai" named "edited.fa.fai"
    When I run `unmux aligned.htslib-3.1.cram --reference edited.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "the .fai index edited.fa.fai is stale: contig `chr1` does not start at its indexed offset 6 in edited.fa"
    And stderr does not contain "has MD5"
    And a file "out" does not exist

  Scenario: a contig truncated after indexing fails as a stale .fai naming that contig
    Given the test file "contigs.cram"
    And a copy of the test file "contigs.fa" named "edited.fa" with "TTTCCTCATGCAATTCAAAACCATG" replaced by ""
    And a copy of the test file "contigs.fa.fai" named "edited.fa.fai"
    When I run `unmux contigs.cram --reference edited.fa --threads 4 --group sm=@RG::SM --sample-from-group sm --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "the .fai index edited.fa.fai is stale: contig `ctgA` has 375 bp, not its indexed 400"
    And a file "out" does not exist

  Scenario: a .fai that lists a contig twice is rejected before any output exists
    Given the test file "aligned.htslib-3.1.cram"
    And a file "dup.fa" containing:
      """
      >chr1
      ACGT
      >chr1
      TTTT
      """
    And a file "dup.fa.fai" containing:
      """
      chr1	4	6	4	5
      chr1	4	17	4	5
      """
    When I run `unmux aligned.htslib-3.1.cram --reference dup.fa --skip-reference-check --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "contig `chr1` is listed more than once in the .fai index dup.fa.fai"
    And a file "out" does not exist

  Scenario: a CRAM on stdin is named stdin in a failed check
    Given the test file "aligned.htslib-3.1.cram"
    And a copy of the test file "ref.fa" named "renamed.fa" with ">chr1" replaced by ">chrX"
    And a copy of the test file "ref.fa.fai" named "renamed.fa.fai" with "chr1" replaced by "chrX"
    When I run `unmux - --reference renamed.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam` with "aligned.htslib-3.1.cram" on stdin
    Then the exit code is 1
    And stderr contains "contig `chr1` (@SQ in stdin) is missing from the --reference FASTA renamed.fa"

  Scenario: a --reference FASTA without its .fai names the missing index
    Given the test file "aligned.htslib-3.1.cram"
    And a copy of the test file "ref.fa" named "unindexed.fa"
    When I run `unmux aligned.htslib-3.1.cram --reference unindexed.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "failed to read the .fai index unindexed.fa.fai of the --reference FASTA unindexed.fa; create it with samtools faidx"
    And a file "out" does not exist

  Scenario: a skipped check still fails on a record whose contig is missing
    Given the test file "aligned.htslib-3.1.cram"
    And a copy of the test file "ref.fa" named "renamed.fa" with ">chr1" replaced by ">chrX"
    And a copy of the test file "ref.fa.fai" named "renamed.fa.fai" with "chr1" replaced by "chrX"
    When I run `unmux aligned.htslib-3.1.cram --reference renamed.fa --skip-reference-check --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "invalid reference sequence name: chr1"
    And stderr does not contain "missing from the --reference FASTA"
    And a file "out" does not exist

  Scenario: a CRAM without M5 is checked by length only, with a warning
    Given the test file "aligned.htslib-3.0.no-m5.cram"
    When I run `unmux aligned.htslib-3.0.no-m5.cram --reference ref.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 0
    And stderr contains "1 of 1 @SQ contig(s) in aligned.htslib-3.0.no-m5.cram have no M5, so only their lengths were checked"
    And the records of "out/rg1.sam" are the "RG:Z:rg1" records of the test file "aligned.htslib.sam"

  Scenario: a CRAM without M5 still fails on a contig with the wrong length
    Given the test file "aligned.htslib-3.0.no-m5.cram"
    And a copy of the test file "ref.fa" named "short.fa"
    And a copy of the test file "ref.fa.fai" named "short.fa.fai" with "3000" replaced by "2990"
    When I run `unmux aligned.htslib-3.0.no-m5.cram --reference short.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "@SQ LN:3000"
    And a file "out" does not exist

  Scenario: a CRAM without M5 still fails on a stale .fai before decoding
    Given the test file "aligned.htslib-3.0.no-m5.cram"
    And a copy of the test file "ref.fa" named "edited.fa" with ">chr1" replaced by ">chr1 assembled"
    And a copy of the test file "ref.fa.fai" named "edited.fa.fai"
    When I run `unmux aligned.htslib-3.0.no-m5.cram --reference edited.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 1
    And stderr contains "the .fai index edited.fa.fai is stale: contig `chr1` does not start at its indexed offset 6 in edited.fa"
    And stderr does not contain "checksum mismatch"
    And a file "out" does not exist

  Scenario: a SAM input ignores --reference
    Given the test file "aligned.sam"
    And a copy of the test file "ref.fa" named "renamed.fa" with ">chr1" replaced by ">chrX"
    And a copy of the test file "ref.fa.fai" named "renamed.fa.fai" with "chr1" replaced by "chrX"
    When I run `unmux aligned.sam --reference renamed.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 0
    And a file "out/rg1.sam" exists
    And stderr contains "--reference is unused: no input is CRAM"

  Scenario: --skip-reference-check needs --reference
    Given the test file "aligned.htslib-3.1.cram"
    When I run `unmux aligned.htslib-3.1.cram --skip-reference-check --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 2
    And stderr contains "--reference"

  Scenario: an M5 that the htslib CRAM writer computed over non-letter bytes passes the check
    Given the test file "punctuation.cram"
    And a file "punctuation.fa" containing:
      """
      >a
      ACGT*-.ACG
      ACGT*-.ACG
      """
    And a file "punctuation.fa.fai" containing:
      """
      a	20	3	10	11
      """
    When I run `unmux punctuation.cram --reference punctuation.fa --group rg=@RG --sample-from-group rg --out out/%sample.sam`
    Then the exit code is 0
    And the records of "out/rg1.sam" are:
      """
      r1	0	a	1	60	1M	*	0	0	A	I	MD:Z:1	NM:i:0	RG:Z:rg1
      """
