Feature: Search windows wider than their tags are warned about
  In a window wider than a tag, a read whose barcode has an error can match another tag at a shifted
  position, on bases of the barcode and beside it, as closely as its true tag. Under mode=all the
  leftmost of those matches wins, and under mode=nearest the group matches no tag on that read. For a
  group whose dist allows an edit, unmux estimates the share of such reads from the slack each tag
  length has in the window (its loc, the next= window another group gives it, or its match= streams
  of fixed width), the number of tags, and the sequences within one edit of a tag, and warns once for
  the group at startup when that share passes 0.3%, giving it. Under mode=all it advises mode=nearest,
  naming the cost of a narrower window. Under mode=nearest it advises a window as wide as the tags,
  or with an indel budget one wider by twice the budget, and is quiet once the window is that wide.
  A group that only whole-group --remove rules read is not warned about under mode=all, where any of
  its tags removes the read.

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
    And stderr contains "group `bc`: the loc window is 1 bp wider than the tags, so for an estimated 0.5% of reads with an error in the barcode another tag matches at a shifted position as closely as the true tag, and the leftmost match wins; use mode=nearest, or narrow loc to 6 bp, losing reads whose barcode an indel upstream shifts"
    And the file "b.fq" contains "@r1"

  Scenario: a window as wide as the tags keeps the true tag, without a warning
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:5:11,dist=1 --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr does not contain "wider than the"
    And the file "a.fq" contains "@r1"

  Scenario: under mode=nearest the tie leaves the read unassigned, with a warning
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1,mode=nearest --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 1 bp wider than the tags, so for an estimated 0.5% of reads with an error in the barcode another tag matches at a shifted position as closely as the true tag, and mode=nearest leaves those reads unmatched; narrow loc to 6 bp if no indel upstream shifts the barcode"
    And the file "un.fq" contains "@r1"

  Scenario: each group with a wide window is warned about once
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1 --group bc2={a=AACCGG,b=GAACCC} --group bc2::loc=0:4:11,dist=1 --extract cb=@bc2 --tag CB=cb --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window" exactly 1 time
    And stderr contains "group `bc2`: the loc window" exactly 1 time

  Scenario: with an indel budget the warning advises mode=nearest alone
    When I run `unmux reads.fq --group bc={a=AACCGG,c=GTTGTT} --group bc::loc=0:4:12,dist=1:1:2 --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 2 bp wider than the tags, so for an estimated 2.2% of reads with an error in the barcode another tag matches at a shifted position as closely as the true tag, and the leftmost match wins; use mode=nearest"
    And stderr does not contain "narrow loc"

  Scenario: with an indel budget and more than twice its slack the warning advises mode=nearest in a narrower window
    When I run `unmux reads.fq --group bc={a=AACCGG,c=GTTGTT} --group bc::loc=0:2:12,dist=1:1:2 --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 4 bp wider than the tags, so for an estimated 4.3% of reads with an error in the barcode another tag matches at a shifted position as closely as the true tag, and the leftmost match wins; use mode=nearest with loc narrowed to 8 bp, which measured best with an indel budget"

  Scenario: under mode=nearest an indel budget with more than twice its slack is warned about
    When I run `unmux reads.fq --group bc={a=AACCGG,c=GTTGTT} --group bc::loc=0:2:12,dist=1:1:2,mode=nearest --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 4 bp wider than the tags, so for an estimated 4.3% of reads with an error in the barcode another tag matches at a shifted position as closely as the true tag, and mode=nearest leaves those reads unmatched; narrow loc to 8 bp, which measured best with an indel budget"

  Scenario: under mode=nearest an indel budget with twice its slack is not warned about
    When I run `unmux reads.fq --group bc={a=AACCGG,c=GTTGTT} --group bc::loc=0:4:12,dist=1:1:2,mode=nearest --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr does not contain "wider than the"

  Scenario: a relative next= window wider than the tags is warned about
    # Tag u matches [0,4), so next=bc:0-7 has group bc search [4,11), as loc=0:4:11 does above.
    When I run `unmux reads.fq --group up={u=TTTT} --group up::loc=0:0:4,next=bc:0-7 --group bc={a=AACCGG,b=GAACCC} --group bc::dist=1 --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the next=bc:0-7 window is 1 bp wider than the tags"
    And stderr contains "or narrow next=bc:0-7 to 6 bp, losing reads whose barcode an indel upstream shifts"
    And stderr does not contain "group `up`: the"
    And the file "b.fq" contains "@r1"

  Scenario: a match= group over a stream wider than the tags is warned about
    # The stream s holds [4,11) of r1, so tag b wins in it as in the loc window above.
    When I run `unmux reads.fq --extract s=0:4:11 --group bc={a=AACCGG,b=GAACCC} --group bc::match=s,dist=1 --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the match=s window is 1 bp wider than the tags"
    And stderr contains "or narrow match=s to 6 bp, losing reads whose barcode an indel upstream shifts"
    And the file "b.fq" contains "@r1"

  Scenario: a window as wide as the longest tag is warned about when shorter tags slide in it
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC,c=TGTGTGTG} --group bc::loc=0:4:12,dist=1 --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 2 bp wider than the shortest tag, so for an estimated 1.0% of reads with an error in the barcode another tag matches at a shifted position as closely as the true tag, and the leftmost match wins; use mode=nearest, or pin the tags with anchor=5p, losing reads whose barcode an indel upstream shifts"
    And the file "b.fq" contains "@r1"

  Scenario Outline: a window of fixed width that lets tags slide is warned about whenever dist allows an edit
    When I run `unmux reads.fq --group <group> --group bc::<attrs> --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is <wider> bp wider than the tags, so for an estimated <share> of reads"

    Examples:
      | group                  | attrs                              | wider | share |
      | bc={a=AACCGG,c=GTTGTT} | loc=0:4:12,dist=0:1                | 2     | 1.3%  |
      | bc={a=AACCGG,c=GTTGTT} | loc=0:-12:-5,dist=1:1:2            | 1     | 1.1%  |
      | bc={a=AACCGG,b=GAACCC} | loc=0:-10,dist=1                   | 4     | 1.8%  |
      | bc={a=AACCGG,b=GAACCC} | loc=0:4:11,dist=1,bothStrands=true | 1     | 0.9%  |
      | bc={a=AACCGG,b=GAACCC} | loc=0:4:11,dist=1,partial5=4       | 1     | 0.5%  |

  Scenario Outline: no warning without an edit, a window of fixed width, a second tag, a sliding tag, or an estimated share past 0.3%
    When I run `unmux reads.fq --group <group> --group bc::<attrs> --sample-from-group bc --out %sample.fq`
    Then the exit code is 0
    And stderr does not contain "wider than the"

    Examples:
      | group                                | attrs                       |
      | bc={a=AACCGG,b=GAACCC}               | loc=0:4:11                  |
      | bc={a=AACCGG,b=GAACCC}               | dist=1                      |
      | bc={a=AACCGG,b=GAACCC}               | loc=0:4,dist=1              |
      | bc={a=AACCGG,b=GAACCC}               | loc=0:5:12,dist=1,anchor=5p |
      | bc={a=AACCGG,b=GAACCC}               | loc=0:4:11,dist=1,anchor=3p |
      | bc={a=AACCGG}                        | loc=0:4:11,dist=1           |
      | bc={a=AACCGGTTAACC,b=GAACCCTTGGAA}   | loc=0:0:16,dist=1           |

  Scenario: a group that only a whole-group --remove reads is not warned about under mode=all
    # Any tag of bc removes r1, so which one wins changes nothing.
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1 --remove bc=removed.fq --out out.fq`
    Then the exit code is 0
    And stderr does not contain "wider than the"
    And the file "removed.fq" contains "@r1"

  Scenario: a group that only a whole-group --remove reads is warned about under mode=nearest, as a tie keeps the read
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1,mode=nearest --remove bc=removed.fq --out out.fq`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 1 bp wider than the tags"
    And stderr contains "and mode=nearest leaves those reads unmatched"
    And the file "out.fq" contains "@r1"

  Scenario Outline: a group that a --remove rule reads by tag, or that also feeds samples or templates, is warned about
    When I run `unmux reads.fq --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1 <uses>`
    Then the exit code is 0
    And stderr contains "group `bc`: the loc window is 1 bp wider than the tags"

    Examples:
      | uses                                                    |
      | --remove bc::a --out out.fq                             |
      | --remove bc --sample-from-group bc --out %sample.fq     |
      | --remove bc --extract cb=@bc --template cb --out out.fq |
