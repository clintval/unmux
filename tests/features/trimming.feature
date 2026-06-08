Feature: Trimming by omission and the unextracted-base metric
  unmux has no trim flags: you "trim" by not extracting the bases you do not want, since anything
  not routed into --template or a --tag is dropped. The per-sample metrics report
  frac_bases_unextracted so that accidentally leaving bases behind is visible.

  Scenario: dropping the 5 prime end by not extracting it
    Given a file "test.fq" containing:
      """
      @r1
      ACGTACGTAGGCCTTAACCGGTTACGT
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # Only bases [9, end) are templated; the first 9 bases are not extracted anywhere, so they are
    # genuinely dropped (not moved into a tag) and never appear in the output.
    When I run `unmux test.fq --extract body=0:9:end --template body`
    Then the exit code is 0
    And stdout contains "GGCCTTAACCGGTTACGT"
    And stdout does not contain "ACGTACGTA"
    And stdout does not contain "ACGTACGTAGGCCTTAACCGGTTACGT"

  Scenario: unextracted bases are reported in the per-sample metrics
    Given a file "reads.fq" containing:
      """
      @r1
      AAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # 100 bp read, only [0, 90) templated -> 10 of 100 bases unextracted -> 0.10.
    When I run `unmux reads.fq --group bc={s1=AAAA} --extract body=0:0:90 --template body --sample S1=bc::s1 --metrics-per-sample m.tsv`
    Then the exit code is 0
    And a file "m.tsv" exists
    And the file "m.tsv" contains "frac_bases_unextracted"
    And the file "m.tsv" contains "0.100000"

  Scenario: the unextracted fraction is zero when the extracts cover the whole read
    Given a file "reads.fq" containing:
      """
      @r1
      AAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group bc={s1=AAAA} --extract body=0:0:end --template body --sample S1=bc::s1 --metrics-per-sample m.tsv`
    Then the exit code is 0
    And the file "m.tsv" contains "0.000000"
    # The sibling scenario above (body=0:0:90) reports 0.100000; full coverage reports none of that,
    # so the zero here is genuinely frac_bases_unextracted=0, not a stray match.
    And the file "m.tsv" does not contain "0.100000"
