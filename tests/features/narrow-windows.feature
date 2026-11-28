Feature: Search windows too narrow for every tag are refused
  A group's tags slide within its loc window, or within the window that a next=GROUP:lo-hi link on
  another group gives it. A tag matches no fewer bases than its length less the indel budget of
  dist, or than a partial5/partial3 minimum, which lets it run off the window edge. When the window
  is narrower than that for every tag, the group could never match, so unmux refuses the run at
  startup and gives the width the window needs. An anchor=5p tag runs past its window end, so its
  window is not checked.

  Background:
    Given a file "reads.fq" containing:
      """
      @r1
      TTTTAACCGGTTTTTT
      +
      IIIIIIIIIIIIIIII
      """

  Scenario: a next= window narrower than every tag is refused
    # Tag u matches [0,4), so next=bc:0-5 has group bc search [4,9), one base short of AACCGG.
    When I run `unmux reads.fq --group up={u=TTTT} --group up::loc=0:0:4,next=bc:0-5 --group bc={a=AACCGG} --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "group `bc`: the next=bc:0-5 window that group `up` gives it is at most 5 bp wide, and no tag can match in fewer than 6 bp, so `bc` could never match; widen it to next=bc:0-6"
    And a file "a.fq" does not exist

  Scenario: a next= window as wide as the tags is accepted
    When I run `unmux reads.fq --group up={u=TTTT} --group up::loc=0:0:4,next=bc:0-6 --group bc={a=AACCGG} --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr does not contain "could never match"
    And the file "a.fq" contains "@r1"

  Scenario: a next= window past its anchor is widened from its own start
    When I run `unmux reads.fq --group up={u=TTTT} --group up::loc=0:0:4,next=bc:2-6 --group bc={a=AACCGG} --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "is at most 4 bp wide, and no tag can match in fewer than 6 bp, so `bc` could never match; widen it to next=bc:2-8"

  Scenario: an anchor=5p tag runs past a narrow next= window
    When I run `unmux reads.fq --group up={u=TTTT} --group up::loc=0:0:4,next=bc:0-2 --group bc={a=AACCGG} --group bc::anchor=5p --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr does not contain "could never match"
    And the file "a.fq" contains "@r1"

  Scenario Outline: a loc window narrower than every tag is refused
    When I run `unmux reads.fq --group bc=<set> --group bc::<attrs> --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "group `bc`: the loc window is at most <width> bp wide, and no tag can match in fewer than <needed> bp, so `bc` could never match; widen loc to <needed> bp"

    Examples:
      | set                    | attrs                | width | needed |
      | {a=AACCGG}             | loc=0:4:9            | 5     | 6      |
      | {a=AACCGG,b=GTTGTTGT}  | loc=0:4:9            | 5     | 6      |
      | {a=AACCGG}             | loc=0:-5:end         | 5     | 6      |
      | {a=AACCGG}             | loc=0:-12:-7         | 5     | 6      |
      | {a=AACCGG}             | loc=0:4:8,dist=1:1:1 | 4     | 5      |
      | {a=AACCGG}             | loc=0:4:6,partial3=3 | 2     | 3      |
      | {a=AACCGG}             | loc=0:4:9,anchor=3p  | 5     | 6      |

  Scenario Outline: a window that some tag can match in is accepted
    When I run `unmux reads.fq --group bc=<set> --group bc::<attrs> --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr does not contain "could never match"
    And the file "a.fq" contains "@r1"

    Examples:
      | set                   | attrs                 |
      | {a=AACCGG}            | loc=0:4:10            |
      | {a=AACCGG,b=GTTGTTGT} | loc=0:4:10            |
      | {a=AACCGG}            | loc=0:4:9,dist=1:1:1  |
      | {a=AACCGG}            | loc=0:4:7,partial3=3  |
      | {a=AACCGG}            | loc=0:4:6,anchor=5p   |
      | {a=AACCGG}            | loc=0:4:10,anchor=3p  |
      | {a=AACCGG}            | loc=0:4               |
