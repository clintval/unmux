Feature: grammar parsing and structural validation
  The quote-free command grammar is parsed and structurally validated before the engine runs, so
  user mistakes surface as fail-fast errors (exit 1) with a clear message. These scenarios exercise
  the parser through the real binary and hold at the grammar layer (no engine required).

  Scenario: non-contiguous input indices are an error
    When I run `unmux --in 0=a.fq --in 2=c.fq`
    Then the exit code is 1
    And stderr contains "contiguous"

  Scenario: redefining a group attribute is an error
    When I run `unmux a.fq --group g=f.tsv --group g::dist=1,dist=2`
    Then the exit code is 1
    And stderr contains "set more than once"

  Scenario: a comma in a tag binding is an error (streams join with +)
    When I run `unmux a.fq --extract a=0:0:4 --extract b=0:4:8 --tag CB=a,b`
    Then the exit code is 1
    And stderr contains "joins streams with"

  Scenario: an extract anchored on an undefined group is an error
    When I run `unmux a.fq --extract x=@missing`
    Then the exit code is 1
    And stderr contains "undefined group"

  Scenario: a group setting both match and loc is an error
    When I run `unmux a.fq --extract i7=0:0:8 --group sb=m.tsv --group sb::match=i7,loc=0:0`
    Then the exit code is 1
    And stderr contains "both"

  Scenario: a boolean attribute rejects numeric values
    When I run `unmux a.fq --group g=f.tsv --group g::both_strands=1`
    Then the exit code is 1
    And stderr contains "true"

  Scenario: a valid pass-through command parses and reaches the engine
    When I run `unmux a.fq --out out.bam`
    Then the exit code is 1
    And stderr contains "failed to open input file"
