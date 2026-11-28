Feature: A --remove rule takes precedence over another group's find bounds
  A record that a --remove rule matches is removed, and counted as removed, even when a group the
  rule does not read fails a find bound such as minFindsPerGroup, whichever group is declared
  first. A group's own bounds still decide whether it matched, so a rule whose group fails its own
  bound, or follows a group that fails one, does not remove the record, which is unassigned for the
  first group to fail. The --qc-tag slug of a removed record lists the groups that matched, and that
  of an unassigned record the groups that matched before the one that failed.

  Background:
    Given a file "r1.fq" containing:
      """
      @s1
      TTTTTTTTACGTACGT
      +
      IIIIIIIIIIIIIIII
      @px
      GATCGATCACGTACGT
      +
      IIIIIIIIIIIIIIII
      @both
      GATCGATCACGTACGT
      +
      IIIIIIIIIIIIIIII
      @none
      TTTTTTTTACGTACGT
      +
      IIIIIIIIIIIIIIII
      """
    And a file "r2.fq" containing:
      """
      @s1
      AAAAAAAACCCCGGGG
      +
      JJJJJJJJJJJJJJJJ
      @px
      TTTTTTTTCCCCGGGG
      +
      JJJJJJJJJJJJJJJJ
      @both
      AAAAAAAACCCCGGGG
      +
      JJJJJJJJJJJJJJJJ
      @none
      TTTTTTTTCCCCGGGG
      +
      JJJJJJJJJJJJJJJJ
      """

  Scenario Outline: a PhiX pair from paired FASTQ is removed though the index group finds nothing, on <threads> threads
    # @px holds the PhiX tag in its first read and no index in its second, which i7 requires.
    When I run `unmux r1.fq r2.fq --pool p --group i7={s1=AAAAAAAA,s2=CCCCCCCC} --group i7::loc=1:0:8,minFindsPerGroup=1 --group phix={p=GATCGATC} --group phix::loc=0:0:8 --sample S1=i7::s1 --out %sample.sam --remove phix=removed.sam --unassigned un.sam --qc-tag --threads <threads> --metrics-per-sample ps.tsv --metrics-summary sum.tsv --metrics-unassigned un.tsv`
    Then the exit code is 0
    And the records of "S1.sam" are:
      """
      s1	77	*	0	0	*	*	0	0	TTTTTTTTACGTACGT	IIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"assigned","sample":"S1","sub_sample":null,"groups":[{"g":"i7","tag":"AAAAAAAA","loc":"1:0:8","obs":"AAAAAAAA","sub":0,"ind":0}]}	RG:Z:S1
      s1	141	*	0	0	*	*	0	0	AAAAAAAACCCCGGGG	JJJJJJJJJJJJJJJJ	ZS:Z:{"v":1,"outcome":"assigned","sample":"S1","sub_sample":null,"groups":[{"g":"i7","tag":"AAAAAAAA","loc":"1:0:8","obs":"AAAAAAAA","sub":0,"ind":0}]}	RG:Z:S1
      """
    And the records of "removed.sam" are:
      """
      px	77	*	0	0	*	*	0	0	GATCGATCACGTACGT	IIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"phix","tag":"GATCGATC","loc":"0:0:8","obs":"GATCGATC","sub":0,"ind":0}]}	RG:Z:p
      px	141	*	0	0	*	*	0	0	TTTTTTTTCCCCGGGG	JJJJJJJJJJJJJJJJ	ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"phix","tag":"GATCGATC","loc":"0:0:8","obs":"GATCGATC","sub":0,"ind":0}]}	RG:Z:p
      both	77	*	0	0	*	*	0	0	GATCGATCACGTACGT	IIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"i7","tag":"AAAAAAAA","loc":"1:0:8","obs":"AAAAAAAA","sub":0,"ind":0},{"g":"phix","tag":"GATCGATC","loc":"0:0:8","obs":"GATCGATC","sub":0,"ind":0}]}	RG:Z:p
      both	141	*	0	0	*	*	0	0	AAAAAAAACCCCGGGG	JJJJJJJJJJJJJJJJ	ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"i7","tag":"AAAAAAAA","loc":"1:0:8","obs":"AAAAAAAA","sub":0,"ind":0},{"g":"phix","tag":"GATCGATC","loc":"0:0:8","obs":"GATCGATC","sub":0,"ind":0}]}	RG:Z:p
      """
    And the records of "un.sam" are:
      """
      none	77	*	0	0	*	*	0	0	TTTTTTTTACGTACGT	IIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"i7","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}	RG:Z:p
      none	141	*	0	0	*	*	0	0	TTTTTTTTCCCCGGGG	JJJJJJJJJJJJJJJJ	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"i7","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}	RG:Z:p
      """
    And the file "ps.tsv" contains "p	S1		1	0.250000	0.000000"
    And the file "sum.tsv" contains "p	4	1	0	1	2	0.250000	0.000000	16,16"
    And the file "un.tsv" is:
      """
      pool	reason	group	stream	reads	frac_of_pool
      p	min_finds_per_group	i7		1	0.250000
      """

    Examples:
      | threads |
      | 1       |
      | 8       |

  Scenario Outline: a PhiX pair from a paired BAM is removed though the index group finds nothing, on <threads> threads
    When I run `unmux r1.fq r2.fq --pool p --out pairs.bam`
    Then the exit code is 0
    When I run `unmux pairs.bam --pool p --group i7={s1=AAAAAAAA,s2=CCCCCCCC} --group i7::loc=1:0:8,minFindsPerGroup=1 --group phix={p=GATCGATC} --group phix::loc=0:0:8 --sample S1=i7::s1 --out %sample.sam --remove phix=removed.sam --unassigned un.sam --qc-tag --threads <threads> --metrics-per-sample ps.tsv --metrics-summary sum.tsv --metrics-unassigned un.tsv`
    Then the exit code is 0
    And the records of "S1.sam" are:
      """
      s1	77	*	0	0	*	*	0	0	TTTTTTTTACGTACGT	IIIIIIIIIIIIIIII	RG:Z:S1	ZS:Z:{"v":1,"outcome":"assigned","sample":"S1","sub_sample":null,"groups":[{"g":"i7","tag":"AAAAAAAA","loc":"1:0:8","obs":"AAAAAAAA","sub":0,"ind":0}]}
      s1	141	*	0	0	*	*	0	0	AAAAAAAACCCCGGGG	JJJJJJJJJJJJJJJJ	RG:Z:S1	ZS:Z:{"v":1,"outcome":"assigned","sample":"S1","sub_sample":null,"groups":[{"g":"i7","tag":"AAAAAAAA","loc":"1:0:8","obs":"AAAAAAAA","sub":0,"ind":0}]}
      """
    And the records of "removed.sam" are:
      """
      px	77	*	0	0	*	*	0	0	GATCGATCACGTACGT	IIIIIIIIIIIIIIII	RG:Z:p	ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"phix","tag":"GATCGATC","loc":"0:0:8","obs":"GATCGATC","sub":0,"ind":0}]}
      px	141	*	0	0	*	*	0	0	TTTTTTTTCCCCGGGG	JJJJJJJJJJJJJJJJ	RG:Z:p	ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"phix","tag":"GATCGATC","loc":"0:0:8","obs":"GATCGATC","sub":0,"ind":0}]}
      both	77	*	0	0	*	*	0	0	GATCGATCACGTACGT	IIIIIIIIIIIIIIII	RG:Z:p	ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"i7","tag":"AAAAAAAA","loc":"1:0:8","obs":"AAAAAAAA","sub":0,"ind":0},{"g":"phix","tag":"GATCGATC","loc":"0:0:8","obs":"GATCGATC","sub":0,"ind":0}]}
      both	141	*	0	0	*	*	0	0	AAAAAAAACCCCGGGG	JJJJJJJJJJJJJJJJ	RG:Z:p	ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"i7","tag":"AAAAAAAA","loc":"1:0:8","obs":"AAAAAAAA","sub":0,"ind":0},{"g":"phix","tag":"GATCGATC","loc":"0:0:8","obs":"GATCGATC","sub":0,"ind":0}]}
      """
    And the records of "un.sam" are:
      """
      none	77	*	0	0	*	*	0	0	TTTTTTTTACGTACGT	IIIIIIIIIIIIIIII	RG:Z:p	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"i7","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}
      none	141	*	0	0	*	*	0	0	TTTTTTTTCCCCGGGG	JJJJJJJJJJJJJJJJ	RG:Z:p	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"i7","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}
      """
    And the file "ps.tsv" contains "p	S1		1	0.250000	0.000000"
    And the file "sum.tsv" contains "p	4	1	0	1	2	0.250000	0.000000	16,16"
    And the file "un.tsv" is:
      """
      pool	reason	group	stream	reads	frac_of_pool
      p	min_finds_per_group	i7		1	0.250000
      """

    Examples:
      | threads |
      | 1       |
      | 8       |

  Scenario Outline: a PhiX read is removed whichever group is declared first: <first> then <second>
    Given a file "lane.fq" containing:
      """
      @px
      GATCGATCTTTTTTTTCCCCGGGG
      +
      IIIIIIIIIIIIIIIIIIIIIIII
      @none
      TTTTTTTTTTTTTTTTCCCCGGGG
      +
      IIIIIIIIIIIIIIIIIIIIIIII
      """
    # Group lib matches both reads, and a slug lists the groups in the order they are declared.
    When I run `unmux lane.fq --pool p --group lib={l=CCCCGGGG} --group lib::loc=0:16:24 --group <first> --group <second> --group i7::loc=0:8:16,minFindsPerGroup=1 --group phix::loc=0:0:8 --sample S1=i7::s1 --out all.sam --remove phix=all.sam --unassigned all.sam --qc-tag`
    Then the exit code is 0
    And the records of "all.sam" are:
      """
      px	4	*	0	0	*	*	0	0	GATCGATCTTTTTTTTCCCCGGGG	IIIIIIIIIIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"lib","tag":"CCCCGGGG","loc":"0:16:24","obs":"CCCCGGGG","sub":0,"ind":0},{"g":"phix","tag":"GATCGATC","loc":"0:0:8","obs":"GATCGATC","sub":0,"ind":0}]}	RG:Z:p
      none	4	*	0	0	*	*	0	0	TTTTTTTTTTTTTTTTCCCCGGGG	IIIIIIIIIIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"i7","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[{"g":"lib","tag":"CCCCGGGG","loc":"0:16:24","obs":"CCCCGGGG","sub":0,"ind":0}]}	RG:Z:p
      """

    Examples:
      | first             | second            |
      | i7={s1=AAAAAAAA}  | phix={p=GATCGATC} |
      | phix={p=GATCGATC} | i7={s1=AAAAAAAA}  |

  Scenario Outline: a PhiX read is removed though the index group fails <bound>
    Given a file "reads.fq" containing:
      """
      @px
      GATCGATCAAAACCCCAAAACCCC
      +
      IIIIIIIIIIIIIIIIIIIIIIII
      @idx
      TTTTTTTTAAAACCCCAAAACCCC
      +
      IIIIIIIIIIIIIIIIIIIIIIII
      """
    # The i7 window [8,24) of each read holds tag s1 twice.
    When I run `unmux reads.fq --pool p --group i7={s1=AAAACCCC,s2=GGGGTTTT} --group i7::loc=0:8:24,<bound> --group phix={p=GATCGATC} --group phix::loc=0:0:8 --out out.sam --remove phix=out.sam --unassigned out.sam --qc-tag`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      px	4	*	0	0	*	*	0	0	GATCGATCAAAACCCCAAAACCCC	IIIIIIIIIIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"phix","tag":"GATCGATC","loc":"0:0:8","obs":"GATCGATC","sub":0,"ind":0}]}	RG:Z:p
      idx	4	*	0	0	*	*	0	0	TTTTTTTTAAAACCCCAAAACCCC	IIIIIIIIIIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"i7","constraint":"<constraint>","found":2,"limit":<limit>,"groups":[]}	RG:Z:p
      """

    Examples:
      | bound              | constraint          | limit |
      | minFindsPerGroup=3 | min_finds_per_group | 3     |
      | maxFindsPerGroup=1 | max_finds_per_group | 1     |
      | minFindsPerTag=3   | min_finds_per_tag   | 3     |
      | maxFindsPerTag=1   | max_finds_per_tag   | 1     |

  Scenario: a remove rule whose group fails its own find bound does not remove the read
    Given a file "own.fq" containing:
      """
      @own
      GATCGATCGATCGATCAAAAAAAA
      +
      IIIIIIIIIIIIIIIIIIIIIIII
      """
    # Tag p matches [0,8), [4,12) and [8,16), three finds where phix allows one.
    When I run `unmux own.fq --pool p --group i7={s1=AAAAAAAA} --group i7::loc=0:16:24,minFindsPerGroup=1 --group phix={p=GATCGATC} --group phix::loc=0:0:16,maxFindsPerGroup=1 --sample S1=i7::s1 --out all.sam --remove phix=all.sam --unassigned all.sam --qc-tag`
    Then the exit code is 0
    And the records of "all.sam" are:
      """
      own	4	*	0	0	*	*	0	0	GATCGATCGATCGATCAAAAAAAA	IIIIIIIIIIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"phix","constraint":"max_finds_per_group","found":3,"limit":1,"groups":[{"g":"i7","tag":"AAAAAAAA","loc":"0:16:24","obs":"AAAAAAAA","sub":0,"ind":0}]}	RG:Z:p
      """

  Scenario: a remove rule whose group follows a group that fails its find bound does not remove the read
    Given a file "after.fq" containing:
      """
      @after
      GATCGATCTTTTTTTT
      +
      IIIIIIIIIIIIIIII
      """
    # phix is matched only after i7 matches, and i7 finds no tag in [8,16).
    When I run `unmux after.fq --pool p --group i7={s1=AAAAAAAA} --group i7::loc=0:8:16,minFindsPerGroup=1 --group phix={p=GATCGATC} --group phix::loc=0:0:8,prev=i7 --sample S1=i7::s1 --out all.sam --remove phix=all.sam --unassigned all.sam --qc-tag`
    Then the exit code is 0
    And the records of "all.sam" are:
      """
      after	4	*	0	0	*	*	0	0	GATCGATCTTTTTTTT	IIIIIIIIIIIIIIII	ZS:Z:{"v":1,"outcome":"unassigned","reason":"find_constraint","group":"i7","constraint":"min_finds_per_group","found":0,"limit":1,"groups":[]}	RG:Z:p
      """
