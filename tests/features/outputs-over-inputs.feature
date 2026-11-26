Feature: An output on an input's path is refused
  Every output is written beside its destination and moved into place when the run ends, so an
  output that is one of the input files would replace that input, and the reads it held would be
  lost with a zero exit. unmux refuses such runs before any output is created, whatever the output
  option and however its path is spelled. The same holds for every other file the run reads: the
  tag files, the sample sheet, and the reference with its index. Standard input from a pipe,
  standard output and /dev/null are not files on a path, so they are exempt, while a metrics path
  is always a file, even one spelled -.

  Background:
    Given a file "in.fq" containing:
      """
      @r1
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      @r2
      TTTTTTTTGGCCATGC
      +
      IIIIIIIIIIIIIIII
      """

  Scenario: a paired demux writing each ordinal over its own input is refused
    Given a file "r1.fq" containing:
      """
      @p1
      CCCCGGGGTTTT
      +
      IIIIIIIIIIII
      @p2
      GGGGTTTTCCCC
      +
      IIIIIIIIIIII
      """
    And a file "r2.fq" containing:
      """
      @p1
      AAAAAAAACCCC
      +
      IIIIIIIIIIII
      @p2
      TTTTTTTTGGGG
      +
      IIIIIIIIIIII
      """
    When I run `unmux r1.fq r2.fq --group bc={a=AAAAAAAA} --group bc::loc=1:0:8 --sample s1=bc::a --out r%ordinal.fq`
    Then the exit code is 1
    And stderr contains "output "r1.fq" (sample s1) is input 0 (r1.fq); give it a distinct path"
    And the file "r1.fq" is:
      """
      @p1
      CCCCGGGGTTTT
      +
      IIIIIIIIIIII
      @p2
      GGGGTTTTCCCC
      +
      IIIIIIIIIIII
      """
    And the file "r2.fq" is:
      """
      @p1
      AAAAAAAACCCC
      +
      IIIIIIIIIIII
      @p2
      TTTTTTTTGGGG
      +
      IIIIIIIIIIII
      """

  Scenario: an input after the first is reported by its 0-based index
    Given a file "r1.fq" containing:
      """
      @p1
      CCCCGGGGTTTT
      +
      IIIIIIIIIIII
      """
    And a file "r2.fq" containing:
      """
      @p1
      AAAAAAAACCCC
      +
      IIIIIIIIIIII
      """
    When I run `unmux r1.fq r2.fq --out out.r%ordinal.fq --metrics-summary r2.fq`
    Then the exit code is 1
    And stderr contains "output "r2.fq" (--metrics-summary) is input 1 (r2.fq); give it a distinct path"
    And the file "r2.fq" is:
      """
      @p1
      AAAAAAAACCCC
      +
      IIIIIIIIIIII
      """
    And a file "out.r1.fq" does not exist

  Scenario Outline: a <option> file on an input's path is refused
    When I run `unmux in.fq --out x.fq <option> in.fq`
    Then the exit code is 1
    And stderr contains "output "in.fq" (<option>) is input 0 (in.fq); give it a distinct path"
    And the file "in.fq" is:
      """
      @r1
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      @r2
      TTTTTTTTGGCCATGC
      +
      IIIIIIIIIIIIIIII
      """
    And a file "x.fq" does not exist

    Examples:
      | option               |
      | --metrics-per-sample |
      | --metrics-summary    |
      | --metrics-unassigned |

  Scenario: a metrics file spelled - is a file, so it is refused on an input file called -
    Given a file "-" containing:
      """
      @r1
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux ./- --out x.fq --metrics-summary -`
    Then the exit code is 1
    And stderr contains "output "-" (--metrics-summary) is input 0 (./-); give it a distinct path"
    And the file "-" contains "@r1"

  Scenario: an --unassigned bin on an input's path is refused
    When I run `unmux in.fq --group bc={a=AAAAAAAA} --sample s1=bc::a --out s1.fq --unassigned in.fq`
    Then the exit code is 1
    And stderr contains "output "in.fq" (--unassigned) is input 0 (in.fq); give it a distinct path"
    And the file "in.fq" contains "@r1"
    And a file "s1.fq" does not exist

  Scenario: a --remove bin on an input's path is refused
    When I run `unmux in.fq --group bc={a=AAAAAAAA} --remove bc=in.fq --out x.fq`
    Then the exit code is 1
    And stderr contains "output "in.fq" (--remove bc) is input 0 (in.fq); give it a distinct path"
    And the file "in.fq" contains "@r2"
    And a file "x.fq" does not exist

  Scenario: a sample whose file is the input is refused
    When I run `unmux in.fq --group bc={in=AAAAAAAA} --sample-from-group bc --out %sample.fq`
    Then the exit code is 1
    And stderr contains "output "in.fq" (sample in) is input 0 (in.fq); give it a distinct path"
    And the file "in.fq" contains "@r2"

  Scenario Outline: an --out of <out> is the input and is refused
    When I run `unmux in.fq --extract a=0:0:4 --template a --out <out>`
    Then the exit code is 1
    And stderr contains "output "<shown>" (--out) is input 0 (in.fq); give it a distinct path"
    And the file "in.fq" is:
      """
      @r1
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      @r2
      TTTTTTTTGGCCATGC
      +
      IIIIIIIIIIIIIIII
      """
    And a file "o" does not exist

    Examples:
      | out        | shown      |
      | in.fq      | in.fq      |
      | ./in.fq    | ./in.fq    |
      | o/../in.fq | o/../in.fq |
      | %pool.fq   | in.fq      |

  Scenario: a pass-through onto its own input is refused even though the bytes would not change
    When I run `unmux in.fq --out in.fq`
    Then the exit code is 1
    And stderr contains "output "in.fq" (--out) is input 0 (in.fq); give it a distinct path"

  Scenario: a SAM pass-through onto its own input is refused
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      a	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      """
    When I run `unmux in.sam --out in.sam`
    Then the exit code is 1
    And stderr contains "output "in.sam" (--out) is input 0 (in.sam); give it a distinct path"
    And the file "in.sam" is:
      """
      @HD	VN:1.6
      a	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII
      """

  Scenario Outline: an output from <option> on a tag file is refused
    Given a file "bc.tsv" containing:
      """
      id	seq
      a	AAAAAAAA
      """
    When I run `unmux in.fq --group bc=bc.tsv --sample-from-group bc --out x.fq <option> bc.tsv`
    Then the exit code is 1
    And stderr contains "output "bc.tsv" (<option>) is the tag file of group bc (bc.tsv); give it a distinct path"
    And the file "bc.tsv" is:
      """
      id	seq
      a	AAAAAAAA
      """
    And a file "x.fq" does not exist

    Examples:
      | option            |
      | --metrics-summary |
      | --unassigned      |

  Scenario: a metrics file on the sample sheet is refused
    Given a file "sheet.tsv" containing:
      """
      sample	bc
      s1	a
      """
    When I run `unmux in.fq --group bc={a=AAAAAAAA} --sample-sheet sheet.tsv --out x.fq --metrics-per-sample sheet.tsv`
    Then the exit code is 1
    And stderr contains "output "sheet.tsv" (--metrics-per-sample) is the sample sheet (sheet.tsv); give it a distinct path"
    And the file "sheet.tsv" is:
      """
      sample	bc
      s1	a
      """

  Scenario Outline: a metrics file on the <what> is refused
    Given the test file "contigs.fa"
    And the test file "contigs.fa.fai"
    When I run `unmux in.fq --reference contigs.fa --out x.fq --metrics-summary <file>`
    Then the exit code is 1
    And stderr contains "output "<file>" (--metrics-summary) is the <what> (<file>); give it a distinct path"
    And the file "<file>" does not contain "reads_assigned"

    Examples:
      | file           | what            |
      | contigs.fa     | reference       |
      | contigs.fa.fai | reference index |

  Scenario: standard input from a pipe is exempt, whatever file fed the pipe
    When I run `unmux - --extract a=0:0:4 --template a --out in.fq` with "in.fq" on stdin
    Then the exit code is 0

  Scenario: an output beside its input is written as before
    When I run `unmux in.fq --extract a=0:0:4 --template a --out out.fq --metrics-summary in.tsv`
    Then the exit code is 0
    And the file "out.fq" contains "AAAA"
    And the file "in.tsv" contains "reads_assigned"
    And the file "in.fq" contains "AAAAAAAACCCCGGGG"

  Scenario: /dev/null as both input and outputs is exempt
    When I run `unmux /dev/null --out /dev/null --metrics-summary /dev/null`
    Then the exit code is 0

  Scenario Outline: standard input to standard output is exempt with <input> and <output>
    When I run `unmux <input> --out <output>` with "in.fq" on stdin
    Then the exit code is 0
    And stdout contains "AAAAAAAACCCCGGGG"

    Examples:
      | input      | output      |
      | -          | -           |
      | /dev/stdin | /dev/stdout |
