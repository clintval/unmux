Feature: Outputs that would be one file are refused
  Two rendered output paths that differ only in case are one file on a case-insensitive file system
  (the macOS and Windows defaults), so their writers would overwrite each other and records would be
  lost with a zero exit. unmux refuses such runs on every operating system, before any output is
  created, so a run that works on Linux cannot lose records on a laptop. The same holds for two
  spellings of one path, for a metrics file sharing a path with another output, and for an
  `--unassigned` or `--remove` bin whose path equals an `--out` path only once placeholders expand.
  Outputs written alike still share one file: samples on one `--out` pattern, bins with each other,
  and a bin given the `--out` pattern itself.

  Background:
    Given a file "reads.fq" containing:
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

  Scenario: inline --sample ids differing only in case are refused
    When I run `unmux reads.fq --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A --sample S1=bc::T --out o/%sample.fq`
    Then the exit code is 1
    And stderr contains "outputs "o/s1.fq" (sample s1) and "o/S1.fq" (sample S1) differ only in case, so they would be one file on a case-insensitive file system; rename one sample"
    And a file "o" does not exist

  Scenario: --sample-sheet ids differing only in case are refused
    Given a file "sheet.tsv" containing:
      """
      sample	bc
      dna01	A
      DNA01	T
      """
    When I run `unmux reads.fq --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample-sheet sheet.tsv --out o/%sample.fq`
    Then the exit code is 1
    And stderr contains "outputs "o/dna01.fq" (sample dna01) and "o/DNA01.fq" (sample DNA01) differ only in case"

  Scenario: --sample-from-group tag ids differing only in case are refused
    When I run `unmux reads.fq --group bc={s1=AAAAAAAA,S1=TTTTTTTT} --sample-from-group bc --out o/%sample.fq`
    Then the exit code is 1
    And stderr contains "outputs "o/s1.fq" (sample s1) and "o/S1.fq" (sample S1) differ only in case"

  Scenario: sub_samples differing only in case are refused
    When I run `unmux reads.fq --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1::lib=bc::A --sample s1::LIB=bc::T --out o/%sample.%sub_sample.fq`
    Then the exit code is 1
    And stderr contains "outputs "o/s1.lib.fq" (sample s1, sub_sample lib) and "o/s1.LIB.fq" (sample s1, sub_sample LIB) differ only in case, so they would be one file on a case-insensitive file system; rename one sub_sample"

  Scenario: a sample output and the --unassigned bin differing only in case are refused
    When I run `unmux reads.fq --group bc={A=AAAAAAAA} --sample un=bc::A --out o/%sample.fq --unassigned o/UN.fq`
    Then the exit code is 1
    And stderr contains "outputs "o/un.fq" (sample un) and "o/UN.fq" (--unassigned) differ only in case"

  Scenario: two spellings of one output path are refused
    When I run `unmux reads.fq --group bc={A=AAAAAAAA} --sample a=bc::A --out o/%sample.fq --unassigned ./o/a.fq`
    Then the exit code is 1
    And stderr contains "outputs "o/a.fq" (sample a) and "./o/a.fq" (--unassigned) are the same file; give them distinct paths"

  Scenario Outline: a bin that shares a sample's own file only once %sample expands is refused
    When I run `unmux reads.fq --pool lib --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample a=bc::A --out o/%sample.fq <bin>`
    Then the exit code is 1
    And stderr contains "outputs "o/a.fq" (sample a) and "o/a.fq" (<option>) are the same file; give them distinct paths"
    And stderr does not contain "the --out pattern itself"
    And a file "o" does not exist

    Examples:
      | bin                   | option         |
      | --unassigned o/a.fq   | --unassigned   |
      | --remove bc::T=o/a.fq | --remove bc::T |

  Scenario: paired outputs per %ordinal and a bin per %source are refused where the numbers meet
    Given a file "r1.fq" containing:
      """
      @r1
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      """
    And a file "r2.fq" containing:
      """
      @r1
      GGGGCCCCAAAATTTT
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux r1.fq r2.fq --group bc={A=AAAAAAAA,T=TTTTTTTT} --group bc::loc=0:0:8 --sample a=bc::A --out s.R%ordinal.fq --unassigned s.R%source.fq`
    Then the exit code is 1
    And stderr contains "outputs "s.R1.fq" (sample a) and "s.R1.fq" (--unassigned) are the same file; give them distinct paths"
    And stderr does not contain "the --out pattern itself"
    And a file "s.R0.fq" does not exist

  Scenario: a bin on a sample's named pipe only once %sample expands is refused
    Given a named pipe "o/a.fq"
    When I run `unmux reads.fq --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample a=bc::A --out o/%sample.fq --unassigned o/a.fq`
    Then the exit code is 1
    And stderr contains "outputs "o/a.fq" (sample a) and "o/a.fq" (--unassigned) are the same file; give them distinct paths"
    When I run `unmux reads.fq --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample a=bc::A --out o/a.fq --unassigned o/a.fq`
    Then the exit code is 0

  Scenario Outline: a bin on another spelling of a sample's named pipe is refused (<bin>)
    Given a named pipe "o/a.fq"
    When I run `unmux reads.fq --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample a=bc::A --out o/%sample.fq --unassigned <bin>`
    Then the exit code is 1
    And stderr contains "outputs "o/a.fq" (sample a) and "<bin>" (--unassigned) are the same file; give them distinct paths"

    Examples:
      | bin         |
      | ./o/a.fq    |
      | o//a.fq     |
      | o/x/../a.fq |

  Scenario Outline: a bin given another path to the one --out file is refused with the way to share it
    When I run `unmux reads.fq --pool lib --group bc={A=AAAAAAAA,T=TTTTTTTT} <outputs>`
    Then the exit code is 1
    And stderr contains "outputs "<out_path>" (<out>) and "<bin_path>" (<option>) are the same file; give them distinct paths, or give <option> the --out pattern itself to share the file"
    And a file "o" does not exist

    Examples:
      | outputs                                                 | out_path | out      | bin_path   | option         |
      | --sample a=bc::A --out o/%pool.fq --unassigned o/lib.fq | o/lib.fq | sample a | o/lib.fq   | --unassigned   |
      | --sample a=bc::A --out o/lib.fq --unassigned o/%pool.fq | o/lib.fq | sample a | o/lib.fq   | --unassigned   |
      | --out o/%pool.fq --remove bc::T=o/lib.fq                | o/lib.fq | --out    | o/lib.fq   | --remove bc::T |
      | --sample a=bc::A --out o/all.fq --unassigned o//all.fq  | o/all.fq | sample a | o//all.fq  | --unassigned   |
      | --sample a=bc::A --out o/all.fq --unassigned o/./all.fq | o/all.fq | sample a | o/./all.fq | --unassigned   |
      | --sample a=bc::A --out o/all.fq --unassigned o/all.fq/  | o/all.fq | sample a | o/all.fq/  | --unassigned   |

  Scenario: an output over an input is refused before any way to share it is suggested
    When I run `unmux reads.fq --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample a=bc::A --out %pool.fq --unassigned reads.fq`
    Then the exit code is 1
    And stderr contains "output "reads.fq" (sample a) is input 0 (reads.fq); give it a distinct path"
    And stderr does not contain "the --out pattern itself"

  Scenario Outline: a bin given the --out pattern itself shares the --out file
    When I run `unmux reads.fq --pool lib --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample a=bc::A <outputs>`
    Then the exit code is 0
    And the file "<path>" is:
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

    Examples:
      | outputs                                   | path     |
      | --out o/all.fq --unassigned o/all.fq      | o/all.fq |
      | --out o/%pool.fq --unassigned o/%pool.fq  | o/lib.fq |
      | --out o/all.fq --remove bc::T=o/all.fq    | o/all.fq |

  Scenario: bins share a file with each other, also once placeholders expand
    When I run `unmux reads.fq --pool lib --group bc={A=AAAAAAAA,T=TTTTTTTT,G=GGGGGGGG} --sample g=bc::G --out o/%sample.fq --unassigned o/lib.fq --remove bc::T=o/%pool.fq`
    Then the exit code is 0
    And the file "o/lib.fq" is:
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

  Scenario: a metrics file sharing a path with a sample output is refused
    When I run `unmux reads.fq --group bc={A=AAAAAAAA} --sample a=bc::A --out o/%sample.fq --metrics-summary o/a.fq`
    Then the exit code is 1
    And stderr contains "outputs "o/a.fq" (sample a) and "o/a.fq" (--metrics-summary) are the same file"

  Scenario: distinct sample ids each keep their own file
    When I run `unmux reads.fq --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A --sample s2=bc::T --out o/%sample.fq --metrics-summary o/summary.tsv`
    Then the exit code is 0
    And the file "o/s1.fq" contains "@r1"
    And the file "o/s1.fq" does not contain "@r2"
    And the file "o/s2.fq" contains "@r2"
    And the file "o/summary.tsv" contains "reads_assigned"

  Scenario: samples sharing one --out file still merge into it
    When I run `unmux reads.fq --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A --sample S1=bc::T --out o/all.fq`
    Then the exit code is 0
    And the file "o/all.fq" contains "@r1"
    And the file "o/all.fq" contains "@r2"
