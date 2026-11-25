Feature: Outputs that would be one file are refused
  Two rendered output paths that differ only in case are one file on a case-insensitive file system
  (the macOS and Windows defaults), so their writers would overwrite each other and records would be
  lost with a zero exit. unmux refuses such runs on every operating system, before any output is
  created, so a run that works on Linux cannot lose records on a laptop. The same holds for two
  spellings of one path, and for a metrics file sharing a path with another output.

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
