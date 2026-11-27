Feature: Search windows wider than their tags are warned about
  In a window wider than a group's longest tag, a read whose barcode has an error can match another
  tag at a shifted position, on bases of the barcode and beside it, as closely as its true tag. Under
  mode=all the leftmost of those matches wins, and under mode=nearest the read is left unassigned.
  So when a group with two or more tags slides them within a dist that allows an edit, and its loc,
  or the next= window another group gives it, is wider than its longest tag, unmux warns once for
  the group at startup. With substitutions only, the warning advises a window as wide as the tags or
  mode=nearest. With an indel budget it advises mode=nearest alone, as a narrower window cuts off
  barcodes with indels, so a group that already uses mode=nearest with one is not warned about.
  Neither is a group that only whole-group --remove rules read, as any of its tags removes a read.

  Background:
    Given a file "reads.fq" containing:
      """
      @r1
      TTTTGAACCGTTTTTT
      +
      IIIIIIIIIIIIIIII
      """

  Scenario: a window one base wider than the tags lets a shifted tag win, with a warning
    # r1 holds tag a (AACCGG) at [5,11) with its last base miscalled as T. In the window [4,11), tag
    # b (GAACCC) also matches with one substitution, at [4,10), and wins as the leftmost match.
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1 --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 1 bp wider than the longest tag, so a read can match another tag at a shifted position as closely as its true tag within dist and the leftmost match wins; narrow loc to 6 bp or use mode=nearest"
    And the file "b.fq" contains "@r1"

  Scenario: a window as wide as the tags keeps the true tag, without a warning
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:5:11,dist=1 --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr does not contain "wider than the longest tag"
    And the file "a.fq" contains "@r1"

  Scenario: under mode=nearest the tie leaves the read unassigned, with a warning
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1,mode=nearest --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 1 bp wider than the longest tag, so a read can match another tag at a shifted position as closely as its true tag within dist and be left unassigned; narrow loc to 6 bp"
    And the file "un.fq" contains "@r1"

  Scenario: with an indel budget the warning advises mode=nearest alone
    When I run `unmux reads.fq --group bc={a=AACCGG,c=GTTGTT} --group bc::loc=0:4:12,dist=1:1:2 --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 2 bp wider than the longest tag, so a read can match another tag at a shifted position as closely as its true tag within dist and the leftmost match wins; use mode=nearest"
    And stderr does not contain "narrow loc"

  Scenario: a relative next= window wider than the tags is warned about
    # Tag u matches [0,4), so next=bc:0-7 has group bc search [4,11), as loc=0:4:11 does above.
    When I run `unmux reads.fq --group up={u=TTTT} --group up::loc=0:0:4,next=bc:0-7 --group bc={a=AACCGG,b=GAACCC} --group bc::dist=1 --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the next=bc:0-7 window is 1 bp wider than the longest tag"
    And stderr contains "narrow next=bc:0-7 to 6 bp or use mode=nearest"
    And stderr does not contain "group `up`: the"
    And the file "b.fq" contains "@r1"

  Scenario Outline: a window of fixed width wider than the longest tag is warned about whenever dist allows an edit
    When I run `unmux reads.fq --group <group> --group bc::<attrs> --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is <wider> bp wider than the longest tag"

    Examples:
      | group                  | attrs                              | wider |
      | bc={a=AACCGG,c=GTTGTT} | loc=0:4:12,dist=0:1                | 2     |
      | bc={a=AACCGG,c=GTTGTT} | loc=0:-12:-5,dist=1:1:2            | 1     |
      | bc={a=AACCGG,b=GAACCC} | loc=0:-10,dist=1                   | 4     |
      | bc={a=AACCGG,b=GAACCC} | loc=0:4:11,dist=1,bothStrands=true | 1     |
      | bc={a=AACCGG,b=GAACCC} | loc=0:4:11,dist=1,partial5=4       | 1     |

  Scenario Outline: no warning without an edit, a fixed window wider than the longest tag, a second tag, or a sliding tag
    When I run `unmux reads.fq --group <group> --group bc::<attrs> --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr does not contain "wider than the longest tag"

    Examples:
      | group                   | attrs                       |
      | bc={a=AACCGG,b=GAACCC}  | loc=0:4:11                  |
      | bc={a=AACCGG,b=GAACCC}  | dist=1                      |
      | bc={a=AACCGG,b=GAACCC}  | loc=0:4,dist=1              |
      | bc={a=AACCGG,b=GAACCC}  | loc=0:5:12,dist=1,anchor=5p |
      | bc={a=AACCGG,b=GAACCC}  | loc=0:4:11,dist=1,anchor=3p |
      | bc={a=AACCGG}           | loc=0:4:11,dist=1           |
      | bc={a=AACCGG,c=TGTGTGT} | loc=0:4:11,dist=1           |

  Scenario: a group that already uses mode=nearest with an indel budget is not warned about
    When I run `unmux reads.fq --group bc={a=AACCGG,c=GTTGTT} --group bc::loc=0:4:12,dist=1:1:2,mode=nearest --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr does not contain "wider than the longest tag"

  Scenario: a group that only a whole-group --remove reads is not warned about
    # Any tag of bc removes r1, so which one wins changes nothing.
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1 --remove bc=removed.fq --out out.fq`
    Then the exit code is 0
    And stderr does not contain "wider than the longest tag"
    And the file "removed.fq" contains "@r1"

  Scenario Outline: a group that a --remove rule reads by tag, or that also feeds samples or templates, is warned about
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1 <uses>`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 1 bp wider than the longest tag"

    Examples:
      | uses                                                    |
      | --remove bc::a --out out.fq                             |
      | --remove bc --sample-from-group bc --out %sample.fq     |
      | --remove bc --extract cb=@bc --template cb --out out.fq |
