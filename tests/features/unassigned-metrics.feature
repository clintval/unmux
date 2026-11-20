Feature: unassigned records counted by reason
  The --metrics-unassigned TSV counts the records of a pool that matched no sample, one row per
  reason seen. A reason is the find-count bound a group failed (the group is given), a --template
  stream whose anchoring group did not match (the stream is given), or a barcode combination that
  no sample claims. The rows add up to the summary's reads_unassigned, and the file is not written
  unless asked for.

  Scenario: each unassigned record is counted under the reason it failed
    Given a file "reads.fq" containing:
      """
      @r1
      AAAATTTTGGGGACAC
      +
      IIIIIIIIIIIIIIII
      @r2
      CCCCTTTTGGGGACAC
      +
      IIIIIIIIIIIIIIII
      @r3
      TTTTTTTTGGGGACAC
      +
      IIIIIIIIIIIIIIII
      @r4
      AAAACCCCGGGGACAC
      +
      IIIIIIIIIIIIIIII
      @r5
      AAAATTTTTTTTACAC
      +
      IIIIIIIIIIIIIIII
      @r6
      AAAATTTTGGGGTTTT
      +
      IIIIIIIIIIIIIIII
      """
    # r1 is assigned. r2 matches bc1::s2, which no sample claims. r3 has no bc1 tag and r4 has both,
    # failing minFindsPerGroup and maxFindsPerGroup. r5 has no bc2 tag. r6 has no bc3 tag, which no
    # bound requires, but the template needs its umi stream.
    When I run `unmux reads.fq --pool lib01 --group bc1={s1=AAAA,s2=CCCC} --group bc1::loc=0:0:8,dist=0,minFindsPerGroup=1,maxFindsPerGroup=1 --group bc2={x=GGGG} --group bc2::loc=0:8:12,dist=0,minFindsPerGroup=1 --group bc3={u=ACAC} --group bc3::loc=0:12:16,dist=0 --extract body=0:0:8 --extract umi=@bc3 --template body+umi --sample dna01=bc1::s1 --out out/%sample.fq --metrics-summary lib01.summary.tsv --metrics-unassigned lib01.unassigned.tsv`
    Then the exit code is 0
    And the file "lib01.summary.tsv" contains "lib01	6	1	0	5	0"
    And the file "lib01.unassigned.tsv" is:
      """
      pool	reason	group	stream	reads	frac_of_pool
      lib01	min_finds_per_group	bc1		1	0.166667
      lib01	max_finds_per_group	bc1		1	0.166667
      lib01	min_finds_per_group	bc2		1	0.166667
      lib01	missing_stream		umi	1	0.166667
      lib01	no_sample			1	0.166667
      """

  Scenario: a pool with nothing unassigned writes only the header
    Given a file "reads.fq" containing:
      """
      @r1
      AAAATTTT
      +
      IIIIIIII
      """
    When I run `unmux reads.fq --pool lib01 --group bc={s1=AAAA} --group bc::loc=0:0:4,dist=0,minFindsPerGroup=1 --sample dna01=bc::s1 --out out/%sample.fq --metrics-unassigned %pool.unassigned.tsv`
    Then the exit code is 0
    And the file "lib01.unassigned.tsv" is:
      """
      pool	reason	group	stream	reads	frac_of_pool
      """

  Scenario: the unassigned metrics are written only when asked for
    Given a file "reads.fq" containing:
      """
      @r1
      TTTTTTTT
      +
      IIIIIIII
      """
    When I run `unmux reads.fq --pool lib01 --group bc={s1=AAAA} --group bc::loc=0:0:4,dist=0,minFindsPerGroup=1 --sample dna01=bc::s1 --out out/%sample.fq --metrics-summary lib01.summary.tsv`
    Then the exit code is 0
    And a file "lib01.unassigned.tsv" does not exist

  Scenario: the unassigned metrics path allows only the pool placeholder
    Given a file "reads.fq" containing:
      """
      @r1
      AAAATTTT
      +
      IIIIIIII
      """
    When I run `unmux reads.fq --group bc={s1=AAAA} --sample dna01=bc::s1 --metrics-unassigned %sample.unassigned.tsv`
    Then the exit code is 1
    And stderr contains "--metrics-unassigned: only `%pool` is valid in a metrics path"
