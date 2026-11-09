Feature: CLI contract
  `unmux` is itself the demux command, with its options flattened at the top level. These
  scenarios hold at the CLI contract layer.

  Scenario: reports its version
    When I run `unmux --version`
    Then the exit code is 0
    And stdout contains "unmux"

  Scenario: unmux is the demux command (its options are flattened at the top level)
    When I run `unmux --help`
    Then the exit code is 0
    And stdout contains "--extract"
    And stdout contains "--compression"
    And stdout contains "--pool"

  Scenario: a bare invocation runs demux without any subcommand word
    When I run `unmux reads.fq --extract r1=0:9:end`
    Then the exit code is 1
    And stderr contains "failed to open input file"

  Scenario: mutually exclusive flags are a usage error
    When I run `unmux --sample dna01=grp::a --sample-sheet sheet.tsv`
    Then the exit code is 2

  Scenario: --sample-sheet-ignore without --sample-sheet is a usage error
    When I run `unmux reads.fq --sample-sheet-ignore assay`
    Then the exit code is 2
    And stderr contains "--sample-sheet <FILE>"

  Scenario: --sample-sheet-ignore with an inline --sample is a usage error
    When I run `unmux reads.fq --sample dna01=grp::a --sample-sheet-ignore assay`
    Then the exit code is 2
    And stderr contains "cannot be used with"

  Scenario: --sample-sheet-ignore skips a metadata column in the sample sheet
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAACGT
      +
      IIIIIIII
      @r2
      CCCCACGT
      +
      IIIIIIII
      """
    And a file "sheet.tsv" containing:
      """
      sample	assay	bc
      dna01	scATAC	A
      dna02	not, a::tag	C
      """
    When I run `unmux reads.fq --group bc={A=AAAA,C=CCCC} --group bc::loc=0:0:4 --sample-sheet sheet.tsv --sample-sheet-ignore assay --out %sample.fq`
    Then the exit code is 0
    And the file "dna01.fq" contains "@r1"
    And the file "dna02.fq" contains "@r2"

  Scenario: an extra sample sheet column is still an error unless ignored
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAACGT
      +
      IIIIIIII
      """
    And a file "sheet.tsv" containing:
      """
      sample	assay	bc
      dna01	scATAC	A
      """
    When I run `unmux reads.fq --group bc={A=AAAA} --group bc::loc=0:0:4 --sample-sheet sheet.tsv --out %sample.fq`
    Then the exit code is 1
    And stderr contains "sample sheet column `assay` is not a declared --group"

  Scenario: setting SM or LB via --rg-tag is rejected
    When I run `unmux reads.fq --rg-tag SM=foo`
    Then the exit code is 1
    And stderr contains "may not set"
