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

  Scenario: setting SM or LB via --rg-tag is rejected
    When I run `unmux reads.fq --rg-tag SM=foo`
    Then the exit code is 1
    And stderr contains "may not set"
