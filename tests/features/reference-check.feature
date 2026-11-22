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
