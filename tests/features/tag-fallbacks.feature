Feature: Tag stream fallbacks
  A stream of a `--tag` may list other streams after `|` to use, in order, when it is absent, so
  `CB=nuc|ba+bb` writes the sheet column when its cell is filled and the barcode `ba` otherwise. The
  `|` binds tighter than `+`. A tag is still left out when every choice for one of its streams is
  absent. A record's quality tag is written only when every chosen stream has base qualities, so a
  record whose tag takes column text has none.

  Background:
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAGGGGACGTACGT
      +
      ABCDEFGHIJKLMNOP
      @r2
      AAAATTTTACGTACGT
      +
      ABCDEFGHIJKLMNOP
      @r3
      CCCCGGGGACGTACGT
      +
      ABCDEFGHIJKLMNOP
      @r4
      CCCCTTTTACGTACGT
      +
      ABCDEFGHIJKLMNOP
      """
    And a file "sheet.tsv" containing:
      """
      sample	a	b	nucleus
      s1	a1	b1	N1
      s1	a1	b2	N1
      s2	a2	b1	
      s2	a2	b2	
      """

  Scenario Outline: a blank column falls back to the barcode
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract ba=@a --extract bb=@b --extract r=0:8:end --tag CB=nuc|ba+bb --tag BC=ba --template r --sample-sheet sheet.tsv --threads <threads> --out out.<ext>`
    Then the exit code is 0
    And the records of "out.<ext>" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:N1-GGGG	BC:Z:AAAA	QT:Z:ABCD	RG:Z:s1
      r2	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:N1-TTTT	BC:Z:AAAA	QT:Z:ABCD	RG:Z:s1
      r3	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:CCCC-GGGG	CY:Z:ABCD EFGH	BC:Z:CCCC	QT:Z:ABCD	RG:Z:s2
      r4	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:CCCC-TTTT	CY:Z:ABCD EFGH	BC:Z:CCCC	QT:Z:ABCD	RG:Z:s2
      """

    Examples:
      | ext | threads |
      | sam | 1       |
      | sam | 8       |
      | bam | 1       |
      | bam | 8       |

  Scenario: a stream with several fallbacks uses the first present one
    Given a file "sheet.tsv" containing:
      """
      sample	a	b	nucleus	other
      s1	a1	b1	N1	O1
      s2	a1	b2		O2
      s3	a2	b1		
      s4	a2	b2	N4	
      """
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract alt=%sample::other --extract ba=@a --extract bb=@b --extract r=0:8:end --tag CB=nuc|alt|ba+bb --template r --sample-sheet sheet.tsv --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:N1-GGGG	RG:Z:s1
      r2	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:O2-TTTT	RG:Z:s2
      r3	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:CCCC-GGGG	CY:Z:ABCD EFGH	RG:Z:s3
      r4	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:N4-TTTT	RG:Z:s4
      """

  Scenario: a tag is left out when every choice for one of its streams is absent
    Given a file "sheet.tsv" containing:
      """
      sample	a	b	nucleus	other
      s1	a1	b1	N1	
      s2	a2	b1		
      """
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract alt=%sample::other --extract bb=@b --extract r=0:8:end --tag xn=nuc|alt+bb --tag xb=bb --template r --sample-sheet sheet.tsv --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	xn:Z:N1-GGGG	xb:Z:GGGG	RG:Z:s1
      r3	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	xb:Z:GGGG	RG:Z:s2
      """

  Scenario Outline: a carried quality tag is dropped only for records whose tag takes column text
    Given a file "carried.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AAAAGGGGACGTACGT	ABCDEFGHIJKLMNOP	CB:Z:OLD	CY:Z:III
      r3	4	*	0	0	*	*	0	0	CCCCGGGGACGTACGT	ABCDEFGHIJKLMNOP	CB:Z:OLD	CY:Z:III
      """
    When I run `unmux carried.sam --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract ba=@a --extract bb=@b --tag CB=nuc|ba+bb <qual> --sample-sheet sheet.tsv --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	AAAAGGGGACGTACGT	ABCDEFGHIJKLMNOP	CB:Z:N1-GGGG	RG:Z:s1
      r3	4	*	0	255	*	*	0	0	CCCCGGGGACGTACGT	ABCDEFGHIJKLMNOP	CB:Z:CCCC-GGGG	CY:Z:ABCD EFGH	RG:Z:s2
      """

    Examples:
      | qual                |
      | --per-record        |
      | --tag CB::qual=CY   |

  Scenario: with qual=none a record whose tag takes bases keeps its carried quality tag
    Given a file "carried.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AAAAGGGGACGTACGT	ABCDEFGHIJKLMNOP	CB:Z:OLD	CY:Z:III
      r3	4	*	0	0	*	*	0	0	CCCCGGGGACGTACGT	ABCDEFGHIJKLMNOP	CB:Z:OLD	CY:Z:III
      """
    When I run `unmux carried.sam --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract ba=@a --extract bb=@b --tag CB=nuc|ba+bb --tag CB::qual=none --sample-sheet sheet.tsv --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	AAAAGGGGACGTACGT	ABCDEFGHIJKLMNOP	CB:Z:N1-GGGG	RG:Z:s1
      r3	4	*	0	255	*	*	0	0	CCCCGGGGACGTACGT	ABCDEFGHIJKLMNOP	CB:Z:CCCC-GGGG	CY:Z:III	RG:Z:s2
      """

  Scenario: a named quality tag is allowed when each part of the tag has a choice with bases
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract ba=@a --extract bb=@b --extract r=0:8:end --tag CB=nuc|ba+bb --tag CB::qual=xq --template r --sample-sheet sheet.tsv --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:N1-GGGG	RG:Z:s1
      r2	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:N1-TTTT	RG:Z:s1
      r3	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:CCCC-GGGG	xq:Z:ABCD EFGH	RG:Z:s2
      r4	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:CCCC-TTTT	xq:Z:ABCD EFGH	RG:Z:s2
      """
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract alt=%sample::a --extract ba=@a --extract bb=@b --extract r=0:8:end --tag CB=nuc|ba+alt --tag CB::qual=xq --template r --sample-sheet sheet.tsv --out out.sam`
    Then the exit code is 1
    And stderr contains "tag `CB` sets `qual=xq`, but its column stream `alt` has no base qualities"

  Scenario: the log says which tags fall back for a blank column
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract ba=@a --extract bb=@b --extract r=0:8:end --tag CB=nuc|ba+bb --tag xn=nuc --template r --sample-sheet sheet.tsv --out out.sam`
    Then the exit code is 0
    And stderr contains "column "nucleus" is blank in 2 of 4 sheet rows; records they route get no xn tag and fall back to the next stream for CB"
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract ba=@a --extract r=0:8:end --tag CB=nuc|ba+nuc --tag CB::qual=none --template r --sample-sheet sheet.tsv --out out.sam`
    Then the exit code is 0
    And stderr contains "column "nucleus" is blank in 2 of 4 sheet rows; records they route get no CB tag"

  Scenario: an empty choice around | is an error
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --extract ba=@a --extract bb=0:4:8 --extract r=0:8:end --template r --out out.sam --tag CB=ba|`
    Then the exit code is 1
    And stderr contains "tag `CB` has an empty choice around `|` in "ba|""
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --extract ba=@a --extract bb=0:4:8 --extract r=0:8:end --template r --out out.sam --tag CB=|ba`
    Then the exit code is 1
    And stderr contains "tag `CB` has an empty choice around `|` in "|ba""
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --extract ba=@a --extract bb=0:4:8 --extract r=0:8:end --template r --out out.sam --tag CB=ba||bb`
    Then the exit code is 1
    And stderr contains "tag `CB` has an empty choice around `|` in "ba||bb""
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --extract ba=@a --extract bb=0:4:8 --extract r=0:8:end --template r --out out.sam --tag CB=bb+ba|`
    Then the exit code is 1
    And stderr contains "tag `CB` has an empty choice around `|` in "ba|""

  Scenario: a fallback must be an extracted stream
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --extract ba=@a --extract r=0:8:end --tag CB=ba|missing --template r --out out.sam`
    Then the exit code is 1
    And stderr contains "tag `CB` references undefined extract stream `missing`"

  Scenario: a --template stream cannot fall back
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --extract ba=@a --extract r=0:8:end --out out.fq --template ba|r`
    Then the exit code is 1
    And stderr contains "a `--template` stream cannot fall back with `|`, which works in `--tag` only; a record whose template stream is absent is unassigned, and `--unassigned` writes it as its input reads. Got: "ba|r""
    And stderr does not contain "invalid characters"
    And a file "out.fq" does not exist
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --extract ba=@a --extract r=0:8:end --out out.fq --template r ~ba|r::raw=true`
    Then the exit code is 1
    And stderr contains "a `--template` stream cannot fall back with `|`, which works in `--tag` only; a record whose template stream is absent is unassigned, and `--unassigned` writes it as its input reads. Got: "~ba|r::raw=true""
