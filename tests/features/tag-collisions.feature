Feature: Tags that one read can match together are refused or warned about
  When one read can match two tags of a group, the tag it is assigned to rests on a tie-break such
  as the order the tags were declared in (mode=all), or it is left unassigned (mode=nearest).
  Before a run, unmux compares every pair of tags under the group's search: the substitution,
  indel, and total caps of its dist, a shorter tag within a longer one (at its start under
  anchor=5p and at its end under anchor=3p), and, with bothStrands=true, each tag against the
  reverse complement of the other. Such a pair is refused under mode=all and warned about under
  mode=nearest, and a pair that one read matches with zero errors is refused under either mode.
  Each message gives bases that match both tags.

  Background:
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTATGGCCCCTTTT
      +
      IIIIIIIIIIIIIIII
      """

  Scenario Outline: an indel budget that lets one read match two tags is refused under mode=all
    # Deleting the last base of either tag leaves ACGTA, so r1 is one deletion from both tags.
    When I run `unmux reads.fq --group bc=<set> --group bc::dist=0:1,loc=0:0:8 --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 1
    And stderr contains "can both match one read within dist=0:1 (the bases `ACGTA` match both)"
    And stderr contains "with mode=all this silently mis-assigns reads"

    Examples:
      | set                 |
      | {a=ACGTAC,b=ACGTAG} |
      | {b=ACGTAG,a=ACGTAC} |

  Scenario: the same indel budget under mode=nearest is a warning and the read is left unassigned
    When I run `unmux reads.fq --group bc={a=ACGTAC,b=ACGTAG} --group bc::dist=0:1,loc=0:0:8,mode=nearest --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "can both match one read within dist=0:1 (the bases `ACGTA` match both)"
    And stderr contains "mode=nearest leaves ambiguous reads unassigned"
    And the file "un.fq" contains "@r1"

  Scenario: without an indel budget the same two tags are told apart by their last base
    Given a file "pair.fq" containing:
      """
      @ra
      ACGTACGG
      +
      IIIIIIII
      @rb
      ACGTAGGG
      +
      IIIIIIII
      """
    When I run `unmux pair.fq --group bc={a=ACGTAC,b=ACGTAG} --group bc::loc=0:0:8 --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "a.fq" contains "@ra"
    And the file "b.fq" contains "@rb"

  Scenario Outline: a tag inside a longer tag is refused under either mode and declaration order
    # Every read that starts with ACGTAC starts with ACGT too, so both tags match it with zero errors.
    When I run `unmux reads.fq --group bc=<set> --group bc::mode=<mode> --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "are indistinguishable (the bases `ACGTAC` match both with zero errors)"

    Examples:
      | set               | mode    |
      | {a=ACGT,b=ACGTAC} | all     |
      | {b=ACGTAC,a=ACGT} | all     |
      | {a=ACGT,b=ACGTAC} | nearest |
      | {b=ACGTAC,a=ACGT} | nearest |

  Scenario Outline: under anchor=5p a tag at the start of a longer tag is refused
    When I run `unmux reads.fq --group bc=<set> --group bc::loc=0:0,anchor=5p --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "are indistinguishable (the bases `ACGTAC` match both with zero errors)"

    Examples:
      | set               |
      | {a=ACGT,b=ACGTAC} |
      | {b=ACGTAC,a=ACGT} |

  Scenario: under anchor=5p a tag one base into a longer tag is told apart
    # CGTA starts one base into ACGTAC, so two tags pinned at loc.start never both match one read.
    Given a file "pinned.fq" containing:
      """
      @rb
      ACGTACTT
      +
      IIIIIIII
      @ra
      CGTATTTT
      +
      IIIIIIII
      """
    When I run `unmux pinned.fq --group bc={a=CGTA,b=ACGTAC} --group bc::loc=0:0,anchor=5p --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "a.fq" contains "@ra"
    And the file "b.fq" contains "@rb"

  Scenario: the same two tags are refused when they slide along the read
    When I run `unmux reads.fq --group bc={a=CGTA,b=ACGTAC} --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "are indistinguishable (the bases `ACGTAC` match both with zero errors)"

  Scenario: under anchor=3p a tag at the end of a longer tag is refused
    When I run `unmux reads.fq --group bc={a=GTAC,b=ACGTAC} --group bc::loc=0:0:8,anchor=3p --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "are indistinguishable (the bases `ACGTAC` match both with zero errors)"

  Scenario: under anchor=3p a tag at the start of a longer tag is told apart
    Given a file "pinned.fq" containing:
      """
      @rb
      TTACGTAC
      +
      IIIIIIII
      @ra
      TTTTACGT
      +
      IIIIIIII
      """
    When I run `unmux pinned.fq --group bc={a=ACGT,b=ACGTAC} --group bc::loc=0:0:8,anchor=3p --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "a.fq" contains "@ra"
    And the file "b.fq" contains "@rb"

  Scenario Outline: tags that are reverse complements are refused with bothStrands=true
    When I run `unmux reads.fq --group bc=<set> --group bc::bothStrands=true,mode=<mode> --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "are indistinguishable (the bases `<bases>` match both with zero errors, one on each strand)"

    Examples:
      | set                     | mode    | bases    |
      | {a=ACGTACGG,b=CCGTACGT} | all     | ACGTACGG |
      | {a=AAAACCCC,b=GGGGTTTT} | all     | AAAACCCC |
      | {a=ACGTACGG,b=CCGTACGT} | nearest | ACGTACGG |

  Scenario: without bothStrands the same two tags are told apart
    Given a file "strands.fq" containing:
      """
      @ra
      TTACGTACGGTT
      +
      IIIIIIIIIIII
      @rb
      TTCCGTACGTTT
      +
      IIIIIIIIIIII
      """
    When I run `unmux strands.fq --group bc={a=ACGTACGG,b=CCGTACGT} --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "a.fq" contains "@ra"
    And the file "b.fq" contains "@rb"

  Scenario: a tag within dist of the reverse complement of another is refused under mode=all
    # The reverse complement of CCGGGTTA is TAACCCGG, one substitution from AAACCCGG.
    When I run `unmux reads.fq --group bc={a=AAACCCGG,b=CCGGGTTA} --group bc::bothStrands=true,dist=1 --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "can both match one read within dist=1 (the bases `AAACCCGG` match both, one on each strand)"
    And stderr contains "with mode=all this silently mis-assigns reads"

  Scenario: the same pair under mode=nearest is a warning
    When I run `unmux reads.fq --group bc={a=AAACCCGG,b=CCGGGTTA} --group bc::bothStrands=true,dist=1,mode=nearest --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "can both match one read within dist=1 (the bases `AAACCCGG` match both, one on each strand)"
    And stderr contains "mode=nearest leaves ambiguous reads unassigned"

  Scenario: tags far apart on both strands pass under an indel budget
    Given a file "apart.fq" containing:
      """
      @ra
      GTGTAAAAAAAAGTGT
      +
      IIIIIIIIIIIIIIII
      @rc
      GTGTCCCCCCCCGTGT
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux apart.fq --group bc={a=AAAAAAAA,c=CCCCCCCC} --group bc::dist=1:1,bothStrands=true --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr does not contain "can both match one read"
    And the file "a.fq" contains "@ra"
    And the file "c.fq" contains "@rc"
