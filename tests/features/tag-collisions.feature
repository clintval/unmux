Feature: Tags that one read can match together are refused or warned about
  When one read can match two tags of a group, the tag it is assigned to rests on a tie-break such
  as the order the tags were declared in (mode=all), or it is left unassigned (mode=nearest).
  Before a run, unmux compares every pair of tags under the group's search: the substitution,
  indel, and total caps of its dist, a shorter tag within a longer one (at its start under
  anchor=5p and at its end under anchor=3p), and, with bothStrands=true, each tag against the
  reverse complement of the other. A pair that one read matches with zero errors is refused under
  either mode, and is reported before any other pair. A pair that one read can match without an
  indel is refused under mode=all and warned about under mode=nearest. A pair that one read can
  match only through an indel is refused when the group gives no mode, asking for one, and warned
  about under mode=all or mode=nearest written out, where mode=all also counts such pairs and
  recommends mode=nearest. Each message gives bases that match both tags, unless the tags are too
  costly to compare exactly under the indel budget. Tags with no anchor slide along a whole record
  or a window wider than they are, and two of them can then also both match one read at shifted
  positions, which the check does not look for.

  Background:
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTATGGCCCCTTTT
      +
      IIIIIIIIIIIIIIII
      """

  Scenario Outline: a pair that one read can match only through an indel is refused until the group gives a mode
    # Deleting the last base of either tag leaves ACGTA, so r1 is one deletion from both tags.
    When I run `unmux reads.fq --group bc=<set> --group bc::dist=0:1,loc=0:0:8 --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 1
    And stderr contains "1 pair of tags can both match one read within dist=0:1 only through an indel, such as"
    And stderr contains "(the bases `ACGTA` match both)"
    And stderr contains "set mode=nearest to leave such reads unassigned (recommended), or write mode=all to keep assigning each to one of the tags by tie-break"

    Examples:
      | set                 |
      | {a=ACGTAC,b=ACGTAG} |
      | {b=ACGTAG,a=ACGTAC} |

  Scenario: with mode=all written out such a pair runs with a warning and the read goes to a tag by tie-break
    When I run `unmux reads.fq --group bc={a=ACGTAC,b=ACGTAG} --group bc::dist=0:1,loc=0:0:8,mode=all --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "tags `ACGTAC` and `ACGTAG` can both match one read within dist=0:1 (the bases `ACGTA` match both); with mode=all such a read goes to one of the tags by tie-break"
    And stderr contains "1 pair of tags can match one read only through an indel; set mode=nearest to leave such reads unassigned (recommended)"
    And the file "a.fq" contains "@r1"

  Scenario Outline: a pair that one read can match without an indel is refused under mode=all despite an indel budget
    # AAAAGCCC is one substitution from each tag, so the read needs no indel.
    When I run `unmux reads.fq --group bc={a=AAAACCCC,b=AAAAGGCC} --group bc::<attrs> --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "tags `AAAACCCC` and `AAAAGGCC` can both match one read within dist=1:1 (the bases `AAAAGCCC` match both)"
    And stderr contains "with mode=all this silently mis-assigns reads"

    Examples:
      | attrs             |
      | dist=1:1          |
      | dist=1:1,mode=all |

  Scenario Outline: bases one indel longer than either tag can match both, so such a pair is caught
    # ACTGTACGT is ACGTACGT with one base inserted, and ACTGTACA with its last base deleted.
    Given a file "shift.fq" containing:
      """
      @r1
      ACTGTACGTGGGGGGGG
      +
      IIIIIIIIIIIIIIIII
      """
    When I run `unmux shift.fq --group bc=<set> --group bc::dist=0:1,loc=0:0:12,anchor=5p --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "1 pair of tags can both match one read within dist=0:1 only through an indel"
    And stderr contains "(the bases `ACTGTACGT` match both)"

    Examples:
      | set                     |
      | {a=ACGTACGT,b=ACTGTACA} |
      | {b=ACTGTACA,a=ACGTACGT} |

  Scenario: under mode=nearest the first colliding pairs are listed and the rest are counted
    # The six tags differ only in their last two bases, so all fifteen pairs collide within dist=1.
    When I run `unmux reads.fq --group bc={a=AAAAA,b=AAAAC,c=AAAAG,d=AAAAT,e=AAACA,f=AAAGA} --group bc::dist=1,mode=nearest --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "tags `AAAAA` and `AAAAC` can both match one read within dist=1"
    And stderr contains "5 more pairs of tags can both match one read within dist=1"

  Scenario: under mode=all written out the pairs that need an indel are listed and counted
    # Without substitutions each pair needs an indel, as AAAA is one deletion from every tag.
    When I run `unmux reads.fq --group bc={a=AAAAA,b=AAAAC,c=AAAAG,d=AAAAT,e=AAACA,f=AAAGA} --group bc::dist=0:1,mode=all --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "tags `AAAAA` and `AAAAC` can both match one read within dist=0:1"
    And stderr contains "5 more pairs of tags can both match one read within dist=0:1; with mode=all such a read goes to one of the tags by tie-break"
    And stderr contains "15 pairs of tags can match one read only through an indel; set mode=nearest to leave such reads unassigned (recommended)"

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

  Scenario Outline: tags that are reverse complements are refused with bothStrands=true, which cannot tell them apart
    # The last two pairs come from published kits: Illumina UD i7 indexes and 10x SI-GA oligos.
    When I run `unmux reads.fq --group bc=<set> --group bc::bothStrands=true,mode=<mode> --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "are reverse complements of each other (the bases `<bases>` match both with zero errors, one on each strand), so bothStrands=true cannot tell them apart; drop bothStrands to match each tag on the forward strand only, or remove one of them"

    Examples:
      | set                         | mode    | bases      |
      | {a=ACGTACGG,b=CCGTACGT}     | all     | ACGTACGG   |
      | {a=AAAACCCC,b=GGGGTTTT}     | all     | AAAACCCC   |
      | {a=ACGTACGG,b=CCGTACGT}     | nearest | ACGTACGG   |
      | {a=ATACCAACGC,b=GCGTTGGTAT} | all     | ATACCAACGC |
      | {a=CTCCTAGT,b=ACTAGGAG}     | nearest | CTCCTAGT   |

  Scenario Outline: a tag within the reverse complement of a longer one is refused with bothStrands=true
    # The reverse complement of ACGTTT is AAACGT, which starts with AAAC.
    When I run `unmux reads.fq --group bc=<set> --group bc::bothStrands=true --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "tag `AAAC` lies within the reverse complement of `ACGTTT` (the bases `ACGTTT` match both with zero errors, one on each strand)"

    Examples:
      | set               |
      | {a=ACGTTT,b=AAAC} |
      | {a=AAAC,b=ACGTTT} |

  Scenario: a pair that one read matches with zero errors is reported before one that needs substitutions
    # AAAA and AATT differ by two substitutions, and every read holding GGGGCC holds GGGG.
    When I run `unmux reads.fq --group bc={a=AAAA,b=AATT,c=GGGGCC,d=GGGG} --group bc::dist=1 --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "tags `GGGGCC` and `GGGG` are indistinguishable (the bases `GGGGCC` match both with zero errors)"
    And stderr does not contain "use mode=nearest"

  Scenario: reverse complements are reported before a pair one substitution apart on opposite strands
    When I run `unmux reads.fq --group bc={a=AAACCCGG,b=CCGGGTTA,c=ATACCAACGC,d=GCGTTGGTAT} --group bc::dist=1,bothStrands=true --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "tags `ATACCAACGC` and `GCGTTGGTAT` are reverse complements of each other"
    And stderr does not contain "use mode=nearest"

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

  Scenario Outline: tags that slide can both match one read at shifted positions without being flagged
    # ACGTAC starts one base into TACGTACGT and CGTACG two bases in, with or without an 8 bp window.
    Given a file "shifted.fq" containing:
      """
      @r1
      TACGTACGT
      +
      IIIIIIIII
      """
    When I run `unmux shifted.fq --group bc={a=ACGTAC,b=CGTACG} --group bc::<attrs> --sample-from-group bc --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr does not contain "can both match one read"
    And stderr does not contain "may both match one read"

    Examples:
      | attrs     |
      | dist=0    |
      | loc=0:0:8 |

  Scenario: a whitelist-sized group of 100,000 tags is checked within a minute
    # Of the five billion pairs of these random 16 bp tags, at least a thousand differ at two bases or fewer.
    Given a tag file "cells.tsv" of 100000 distinct random 16 bp tags
    When I run `unmux reads.fq --group cb=cells.tsv --group cb::loc=0:0:16,dist=1,mode=nearest --extract cell=@cb --tag CB=cell --out out.sam` allowing 60 seconds
    Then the run finished in the time allowed
    And the exit code is 0
    And stderr contains "tags `TTCTGAGTTGCAGTGA` and `TTCGGAGTTCCAGTGA` can both match one read within dist=1 (the bases `TTCGGAGTTGCAGTGA` match both)"
    And stderr contains "at least 990 more pairs (counting stopped at 1000) of tags can both match one read within dist=1"
