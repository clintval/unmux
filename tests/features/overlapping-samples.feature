Feature: Samples that one read can match together are warned about
  A record goes to the first sample whose selector it satisfies: the first --sample flag, or the
  first row of a --sample-sheet. Two samples that claim the same tag combination are refused, but
  selectors over different groups, as sheet rows that fill different group columns are, can both
  match one read. Before a run, unmux warns about each such selector that loses a read to an
  earlier one of another sample, giving the tags of one such read and the sample it goes to. It
  lists ten and counts the rest.

  Background:
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAAAAGGGGGGGG
      +
      IIIIIIIIIIIIIIII
      @r2
      AAAAAAAATTTTTTTT
      +
      IIIIIIIIIIIIIIII
      @r3
      CCCCCCCCGGGGGGGG
      +
      IIIIIIIIIIIIIIII
      """

  Scenario: rows that fill different group columns and both match one read are warned about
    # r1 holds i7 a and i5 b, so it matches both rows and goes to the first.
    Given a file "sheet.tsv" containing:
      """
      sample	i7	i5
      s1	a	
      s2		b
      """
    When I run `unmux reads.fq --group i7={a=AAAAAAAA,c=CCCCCCCC} --group i7::loc=0:0:8 --group i5={b=GGGGGGGG,d=TTTTTTTT} --group i5::loc=0:8:16 --sample-sheet sheet.tsv --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "sample sheet rows 1 (`s1`) and 2 (`s2`) can both route one read, such as one with the combination `i5::b+i7::a`; rows are tried in order, so such reads go to `s1`"
    And the file "s1.fq" contains "@r1"
    And the file "s2.fq" contains "@r3"
    And the file "s2.fq" does not contain "@r1"

  Scenario: a row that an earlier row of another sample covers is warned about
    # Every read that row 2 matches also matches row 1, so s2 gets no reads.
    Given a file "sheet.tsv" containing:
      """
      sample	i7	i5
      s1	a	
      s2	a	b
      """
    When I run `unmux reads.fq --group i7={a=AAAAAAAA,c=CCCCCCCC} --group i7::loc=0:0:8 --group i5={b=GGGGGGGG,d=TTTTTTTT} --group i5::loc=0:8:16 --sample-sheet sheet.tsv --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "sample sheet rows 1 (`s1`) and 2 (`s2`) can both route one read, such as one with the combination `i5::b+i7::a`; rows are tried in order, so such reads go to `s1`"
    And the file "s1.fq" contains "@r1"
    And the file "s2.fq" does not contain "@r1"

  Scenario: a row listed before a row it narrows is warned about too
    Given a file "sheet.tsv" containing:
      """
      sample	i7	i5
      s2	a	b
      s1	a	
      """
    When I run `unmux reads.fq --group i7={a=AAAAAAAA,c=CCCCCCCC} --group i7::loc=0:0:8 --group i5={b=GGGGGGGG,d=TTTTTTTT} --group i5::loc=0:8:16 --sample-sheet sheet.tsv --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "sample sheet rows 1 (`s2`) and 2 (`s1`) can both route one read, such as one with the combination `i5::b+i7::a`; rows are tried in order, so such reads go to `s2`"
    And the file "s2.fq" contains "@r1"
    And the file "s1.fq" contains "@r2"

  Scenario: --sample selectors over different groups that both match one read are warned about
    When I run `unmux reads.fq --group i7={a=AAAAAAAA,c=CCCCCCCC} --group i7::loc=0:0:8 --group i5={b=GGGGGGGG,d=TTTTTTTT} --group i5::loc=0:8:16 --sample s1=i7::a --sample s2=i5::b --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "samples `s1` and `s2` can both route one read, such as one with the combination `i5::b+i7::a`; --sample targets are tried in order, so such reads go to `s1`"
    And the file "s1.fq" contains "@r1"

  Scenario: a sub-sample is a sample of its own
    When I run `unmux reads.fq --group i7={a=AAAAAAAA,c=CCCCCCCC} --group i7::loc=0:0:8 --group i5={b=GGGGGGGG,d=TTTTTTTT} --group i5::loc=0:8:16 --sample s1::x=i7::a --sample s1::y=i5::b --out %sample.%sub_sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "samples `s1.x` and `s1.y` can both route one read, such as one with the combination `i5::b+i7::a`; --sample targets are tried in order, so such reads go to `s1.x`"

  Scenario: rows that fill the same group columns are not warned about
    Given a file "sheet.tsv" containing:
      """
      sample	i7	i5
      s1	a	b
      s2	c	d
      """
    When I run `unmux reads.fq --group i7={a=AAAAAAAA,c=CCCCCCCC} --group i7::loc=0:0:8 --group i5={b=GGGGGGGG,d=TTTTTTTT} --group i5::loc=0:8:16 --sample-sheet sheet.tsv --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr does not contain "can both route one read"

  Scenario: rows that differ in a group they both fill are not warned about
    Given a file "sheet.tsv" containing:
      """
      sample	i7	i5
      s1	a	
      s2	c	b
      """
    When I run `unmux reads.fq --group i7={a=AAAAAAAA,c=CCCCCCCC} --group i7::loc=0:0:8 --group i5={b=GGGGGGGG,d=TTTTTTTT} --group i5::loc=0:8:16 --sample-sheet sheet.tsv --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr does not contain "can both route one read"
    And the file "s1.fq" contains "@r1"
    And the file "s2.fq" contains "@r3"

  Scenario: rows of one sample are not warned about
    Given a file "sheet.tsv" containing:
      """
      sample	i7	i5
      s1	a	
      s1		b
      """
    When I run `unmux reads.fq --group i7={a=AAAAAAAA,c=CCCCCCCC} --group i7::loc=0:0:8 --group i5={b=GGGGGGGG,d=TTTTTTTT} --group i5::loc=0:8:16 --sample-sheet sheet.tsv --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr does not contain "can both route one read"

  Scenario: a read that an earlier row of the losing sample takes first is not warned about
    # Row 3 needs i7 a and i5 b, but row 1 of the same sample takes every read with i5 b first.
    Given a file "sheet.tsv" containing:
      """
      sample	i7	i5
      s2		b
      s1	a	
      s2	a	b
      """
    When I run `unmux reads.fq --group i7={a=AAAAAAAA,c=CCCCCCCC} --group i7::loc=0:0:8 --group i5={b=GGGGGGGG,d=TTTTTTTT} --group i5::loc=0:8:16 --sample-sheet sheet.tsv --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "sample sheet rows 1 (`s2`) and 2 (`s1`) can both route one read, such as one with the combination `i5::b+i7::a`; rows are tried in order, so such reads go to `s2`"
    And stderr contains "can both route one read" exactly 1 time

  Scenario: --sample-from-group is not warned about
    When I run `unmux reads.fq --group i7={a=AAAAAAAA,c=CCCCCCCC} --group i7::loc=0:0:8 --sample-from-group i7 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr does not contain "can both route one read"

  Scenario: past ten rows that lose reads the rest are counted
    Given a file "sheet.tsv" containing:
      """
      sample	i7	i5
      x	a	
      y1		b1
      y2		b2
      y3		b3
      y4		b4
      y5		b5
      y6		b6
      y7		b7
      y8		b8
      y9		b9
      y10		b10
      y11		b11
      """
    When I run `unmux reads.fq --group i7={a=AAAAAAAA} --group i7::loc=0:0:8 --group i5={b1=GGGGGGGG,b2=TTTTTTTT,b3=ACACACAC,b4=CACACACA,b5=AGAGAGAG,b6=GAGAGAGA,b7=ATATATAT,b8=TATATATA,b9=CGCGCGCG,b10=GCGCGCGC,b11=CTCTCTCT} --group i5::loc=0:8:16 --sample-sheet sheet.tsv --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And stderr contains "can both route one read" exactly 10 times
    And stderr contains "sample sheet rows 1 (`x`) and 11 (`y10`) can both route one read"
    And stderr does not contain "(`y11`) can both route one read"
    And stderr contains "1 more sample sheet row can lose reads to an earlier row of another sample"
