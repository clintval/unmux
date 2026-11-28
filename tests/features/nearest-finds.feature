Feature: Under mode=nearest only the settled match is a find
  Under mode=nearest the find bounds count the group's finds after nearest settles it. A clear
  winner is one find, for the group and for its tag, and a runner-up within dist is no find, so
  findOne and maxFindsPerGroup=1 keep a read that nearest resolved, and a minFindsPerGroup or
  minFindsPerTag above 1, which no read could meet, is refused before any read. A tie, where no
  match beats the runner-up by delta, fails every find bound the group has: the read is
  unassigned as nearest_tie in its --qc-tag slug and its --metrics-unassigned row, and the run's
  log counts such reads by group. A group with no find bound leaves a tied read unmatched, and
  under mode=all every match counts as before.

  Background:
    Given a file "r1.fq" containing:
      """
      @win
      TTTTGAACCGGTTTTT
      +
      IIIIIIIIIIIIIIII
      @tie
      TTTTGAACCGTTTTTT
      +
      IIIIIIIIIIIIIIII
      @none
      TTTTTTTTTTTTTTTT
      +
      IIIIIIIIIIIIIIII
      """
    And a file "r2.fq" containing:
      """
      @win
      CCCCAAAA
      +
      JJJJJJJJ
      @tie
      GGGGAAAA
      +
      JJJJJJJJ
      @none
      AAAACCCC
      +
      JJJJJJJJ
      """

  Scenario Outline: findOne keeps a won pair and unassigns a tied one from paired FASTQ on <threads> threads
    # In [4,11), @win holds tag a exactly with tag b one substitution away, and @tie is one substitution from both.
    When I run `unmux r1.fq r2.fq --pool p --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1,mode=nearest,findOne --sample-from-group bc --out %sample.sam --unassigned un.sam --qc-tag --threads <threads> --metrics-per-sample ps.tsv --metrics-summary sum.tsv --metrics-unassigned un.tsv`
    Then the exit code is 0
    And the records of "a.sam" are:
      """
      win	77	*	0	0	*	*	0	0	TTTTGAACCGGTTTTT	IIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"assigned","sample":"a","sub_sample":null,"groups":[{"g":"bc","tag":"AACCGG","loc":"0:5:11","obs":"AACCGG","sub":0,"ind":0}]}	RG:Z:a
      win	141	*	0	0	*	*	0	0	CCCCAAAA	JJJJJJJJ	ZS:Z:{"v":1,"outcome":"assigned","sample":"a","sub_sample":null,"groups":[{"g":"bc","tag":"AACCGG","loc":"0:5:11","obs":"AACCGG","sub":0,"ind":0}]}	RG:Z:a
      """
    And the records of "un.sam" are:
      """
      tie	77	*	0	0	*	*	0	0	TTTTGAACCGTTTTTT	IIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"unassigned","reason":"nearest_tie","group":"bc","groups":[]}	RG:Z:p
      tie	141	*	0	0	*	*	0	0	GGGGAAAA	JJJJJJJJ	ZS:Z:{"v":1,"outcome":"unassigned","reason":"nearest_tie","group":"bc","groups":[]}	RG:Z:p
      none	77	*	0	0	*	*	0	0	TTTTTTTTTTTTTTTT	IIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"bc","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}	RG:Z:p
      none	141	*	0	0	*	*	0	0	AAAACCCC	JJJJJJJJ	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"bc","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}	RG:Z:p
      """
    And the file "ps.tsv" contains "p	a		1	0.333333	0.000000"
    And the file "sum.tsv" contains "p	3	1	0	2	0	0.333333	0.000000	16,8"
    And the file "un.tsv" is:
      """
      pool	reason	group	stream	reads	frac_of_pool
      p	min_finds_per_group	bc		1	0.333333
      p	nearest_tie	bc		1	0.333333
      """
    And stderr contains "metrics pool=p: 1 reads unassigned as nearest_tie (bc=1)"
    And stderr does not contain "min_finds_per_group"

    Examples:
      | threads |
      | 1       |
      | 8       |

  Scenario Outline: findOne keeps a won pair and unassigns a tied one from a paired BAM on <threads> threads
    When I run `unmux r1.fq r2.fq --pool p --out pairs.bam`
    Then the exit code is 0
    When I run `unmux pairs.bam --pool p --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1,mode=nearest,findOne --sample-from-group bc --out %sample.sam --unassigned un.sam --qc-tag --threads <threads> --metrics-per-sample ps.tsv --metrics-summary sum.tsv --metrics-unassigned un.tsv`
    Then the exit code is 0
    And the records of "a.sam" are:
      """
      win	77	*	0	0	*	*	0	0	TTTTGAACCGGTTTTT	IIIIIIIIIIIIIIII	RG:Z:a	ZS:Z:{"v":1,"outcome":"assigned","sample":"a","sub_sample":null,"groups":[{"g":"bc","tag":"AACCGG","loc":"0:5:11","obs":"AACCGG","sub":0,"ind":0}]}
      win	141	*	0	0	*	*	0	0	CCCCAAAA	JJJJJJJJ	RG:Z:a	ZS:Z:{"v":1,"outcome":"assigned","sample":"a","sub_sample":null,"groups":[{"g":"bc","tag":"AACCGG","loc":"0:5:11","obs":"AACCGG","sub":0,"ind":0}]}
      """
    And the records of "un.sam" are:
      """
      tie	77	*	0	0	*	*	0	0	TTTTGAACCGTTTTTT	IIIIIIIIIIIIIIII	RG:Z:p	ZS:Z:{"v":1,"outcome":"unassigned","reason":"nearest_tie","group":"bc","groups":[]}
      tie	141	*	0	0	*	*	0	0	GGGGAAAA	JJJJJJJJ	RG:Z:p	ZS:Z:{"v":1,"outcome":"unassigned","reason":"nearest_tie","group":"bc","groups":[]}
      none	77	*	0	0	*	*	0	0	TTTTTTTTTTTTTTTT	IIIIIIIIIIIIIIII	RG:Z:p	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"bc","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}
      none	141	*	0	0	*	*	0	0	AAAACCCC	JJJJJJJJ	RG:Z:p	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"bc","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}
      """
    And the file "ps.tsv" contains "p	a		1	0.333333	0.000000"
    And the file "sum.tsv" contains "p	3	1	0	2	0	0.333333	0.000000	16,8"
    And the file "un.tsv" is:
      """
      pool	reason	group	stream	reads	frac_of_pool
      p	min_finds_per_group	bc		1	0.333333
      p	nearest_tie	bc		1	0.333333
      """
    And stderr contains "metrics pool=p: 1 reads unassigned as nearest_tie (bc=1)"
    And stderr does not contain "min_finds_per_group"

    Examples:
      | threads |
      | 1       |
      | 8       |

  Scenario Outline: a tie fails <bound> and is unassigned as nearest_tie
    Given a file "ties.fq" containing:
      """
      @tied
      AAAACCCCGGGGTTTT
      +
      IIIIIIIIIIIIIIII
      @twice
      AAAACCCCGGGGTTTTAAAACCCC
      +
      IIIIIIIIIIIIIIIIIIIIIIII
      """
    # Tags a and b both match each read exactly, and @twice holds tag a twice.
    When I run `unmux ties.fq --pool p --group bc={a=AAAACCCC,b=GGGGTTTT} --group bc::mode=nearest,<bound> --out out.sam --unassigned out.sam --qc-tag --metrics-unassigned un.tsv`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      tied	4	*	0	0	*	*	0	0	AAAACCCCGGGGTTTT	IIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"unassigned","reason":"nearest_tie","group":"bc","groups":[]}	RG:Z:p
      twice	4	*	0	0	*	*	0	0	AAAACCCCGGGGTTTTAAAACCCC	IIIIIIIIIIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"unassigned","reason":"nearest_tie","group":"bc","groups":[]}	RG:Z:p
      """
    And the file "un.tsv" is:
      """
      pool	reason	group	stream	reads	frac_of_pool
      p	nearest_tie	bc		2	1.000000
      """
    And stderr contains "metrics pool=p: 2 reads unassigned as nearest_tie (bc=2)"

    Examples:
      | bound              |
      | minFindsPerGroup=1 |
      | maxFindsPerGroup=1 |
      | minFindsPerTag=1   |
      | maxFindsPerTag=1   |
      | findOne            |

  Scenario: a tie in a group with no find bound leaves the group unmatched
    Given a file "ties.fq" containing:
      """
      @tied
      AAAACCCCGGGGTTTT
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux ties.fq --pool p --group bc={a=AAAACCCC,b=GGGGTTTT} --group bc::mode=nearest --extract cb=@bc --tag CB=cb --out out.sam --unassigned out.sam --qc-tag`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      tied	4	*	0	0	*	*	0	0	AAAACCCCGGGGTTTT	IIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"pass_through","groups":[]}	RG:Z:p
      """

  Scenario Outline: a clear winner is one find, so <bound> keeps it
    Given a file "winners.fq" containing:
      """
      @runner
      GAACCGGTTTTGAACCATTTT
      +
      IIIIIIIIIIIIIIIIIIIII
      @repeat
      AACCGGTTTTAACCGGTTTT
      +
      IIIIIIIIIIIIIIIIIIII
      """
    # Tag a wins in each read: @runner also holds tag b one substitution away twice, and @repeat holds tag a twice.
    When I run `unmux winners.fq --pool p --group bc={a=AACCGG,b=GAACCC} --group bc::dist=1,mode=nearest,<bound> --sample-from-group bc --out all.sam --unassigned all.sam --qc-tag`
    Then the exit code is 0
    And the records of "all.sam" are:
      """
      runner	4	*	0	0	*	*	0	0	GAACCGGTTTTGAACCATTTT	IIIIIIIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"assigned","sample":"a","sub_sample":null,"groups":[{"g":"bc","tag":"AACCGG","loc":"0:1:7","obs":"AACCGG","sub":0,"ind":0}]}	RG:Z:a
      repeat	4	*	0	0	*	*	0	0	AACCGGTTTTAACCGGTTTT	IIIIIIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"assigned","sample":"a","sub_sample":null,"groups":[{"g":"bc","tag":"AACCGG","loc":"0:0:6","obs":"AACCGG","sub":0,"ind":0}]}	RG:Z:a
      """

    Examples:
      | bound              |
      | findOne            |
      | maxFindsPerGroup=1 |
      | maxFindsPerTag=1   |

  Scenario Outline: <bound> under mode=nearest, which no read could meet, is refused before any read
    When I run `unmux r1.fq r2.fq --pool p --group bc={a=AACCGG,b=GAACCC} --group bc::loc=0:4:11,dist=1,mode=nearest,<bound> --sample-from-group bc --out %sample.sam --unassigned un.sam`
    Then the exit code is 1
    And stderr contains "group `bc`: <bound> can never be met under mode=nearest, which settles on one match per group"
    And a file "a.sam" does not exist
    And a file "un.sam" does not exist

    Examples:
      | bound              |
      | minFindsPerGroup=2 |
      | minFindsPerTag=3   |

  Scenario Outline: under mode=all every match counts toward <bound>, as before
    Given a file "ties.fq" containing:
      """
      @tied
      AAAACCCCGGGGTTTT
      +
      IIIIIIIIIIIIIIII
      @twice
      AAAACCCCGGGGTTTTAAAACCCC
      +
      IIIIIIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux ties.fq --pool p --group bc={a=AAAACCCC,b=GGGGTTTT} --group bc::mode=all,<bound> --out out.sam --unassigned out.sam --qc-tag --metrics-summary sum.tsv --metrics-unassigned un.tsv`
    Then the exit code is 0
    And the file "sum.tsv" contains "p	2	0	<passed>	<unassigned>	0	"
    And the file "un.tsv" is:
      """
      pool	reason	group	stream	reads	frac_of_pool
      p	<reason>	bc		<unassigned>	<frac>
      """

    Examples:
      | bound              | passed | unassigned | reason              | frac     |
      | minFindsPerGroup=3 | 1      | 1          | min_finds_per_group | 0.500000 |
      | maxFindsPerGroup=1 | 0      | 2          | max_finds_per_group | 1.000000 |
      | minFindsPerTag=2   | 0      | 2          | min_finds_per_tag   | 1.000000 |
      | maxFindsPerTag=1   | 1      | 1          | max_finds_per_tag   | 0.500000 |
      | findOne            | 0      | 2          | max_finds_per_group | 1.000000 |
