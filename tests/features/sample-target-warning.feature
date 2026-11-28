Feature: The warning that no reads are split by sample
  A run that defines groups but gives no --sample, --sample-sheet or --sample-from-group writes every
  read to --out without splitting it by sample. That is intended when each group feeds something: a
  --tag or --template through an @group extract, a match= group, a next= or prev= link, a --remove
  rule, or a find constraint. So unmux warns only when a group feeds nothing, and lists those groups.

  Scenario Outline: a group that <use> raises no warning
    Given a file "reads.fq" containing:
      """
      @r1
      AAAACCCCGGGG
      +
      IIIIIIIIIIII
      """
    And a file "cells.tsv" containing:
      """
      id	seq	cell
      c1	AAAA	first
      """
    When I run `unmux reads.fq <args> --out out.sam`
    Then the exit code is 0
    And stderr does not contain "WARN"

    Examples:
      | use                                  | args                                                                                                   |
      | feeds a tag                          | --group cb=cells.tsv --group cb::loc=0:0:4 --extract cell=@cb --tag CB=cell                            |
      | feeds a tag from a tag-file column   | --group cb=cells.tsv --group cb::loc=0:0:4 --extract label=@cb::cell --tag xn=label                    |
      | feeds a template                     | --group cb=cells.tsv --group cb::loc=0:0:4 --extract cell=@cb --template cell                          |
      | feeds a match= group                 | --group cb=cells.tsv --group cb::loc=0:0:4 --group m={m1=AAAA} --group m::match=cell --extract cell=@cb |
      | leads another group with next=       | --group a={a1=AAAA} --group a::loc=0:0:4,next=b:0-4 --group b={b1=CCCC}                                |
      | follows another group with prev=     | --group a={a1=AAAA} --group a::loc=0:0:4 --group b={b1=CCCC} --group b::loc=0:4:8,prev=a              |
      | selects reads for --remove           | --group ad={x=GGGG} --group ad::loc=0:8:12 --remove ad=removed.fq                                      |
      | gates reads with a find constraint   | --group site={x=CCCC} --group site::loc=0:4:8,minFindsPerGroup=1                                       |

  Scenario: a group that feeds nothing still raises the warning, which lists it
    Given a file "reads.fq" containing:
      """
      @r1
      AAAACCCCGGGG
      +
      IIIIIIIIIIII
      """
    When I run `unmux reads.fq --group cb={c1=AAAA} --group cb::loc=0:0:4 --group g={x=CCCC} --group h={y=TTTT} --extract cell=@cb --tag CB=cell --out out.sam`
    Then the exit code is 0
    And stderr contains "groups are defined but no `--sample`, `--sample-sheet`, or `--sample-from-group` was given; reads will not be demultiplexed (they pass through to `--out` undemultiplexed). Add a sample target to name a single sample or split by sample. Unused groups: `g`, `h`."
    And stderr contains "group `g` is defined but nothing consumes its matches"
    And stderr does not contain "group `cb` is defined but nothing consumes its matches"
