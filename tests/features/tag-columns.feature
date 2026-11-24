Feature: Tag-file column streams
  A tag file may carry columns beyond `id` and `seq`. The extract source `@grp::COLUMN` is the
  matched tag's text in that column, for use in `--tag` alone or joined with other streams. A blank
  value is the tag's `seq`, `@grp::id` is the tag ID, and `@grp::seq` is `@grp` itself. Column text
  has no base qualities, so it cannot go in a `--template`, a group `match=`, or under `~`, and a tag
  that uses it has no quality tag. Columns that no extract reads are ignored, whatever they hold.

  Background:
    Given a file "reads.fq" containing:
      """
      @r1
      AACCGGACGTAAAACCCC
      +
      ABCDEFGHIJKLMNOPQR
      @r2
      AACCGTACGTAAAACCCC
      +
      ABCDEFGHIJKLMNOPQR
      @r3
      TTGGCCTTTTAAAACCCC
      +
      ABCDEFGHIJKLMNOPQR
      @r4
      GAGAGATTTTAAAACCCC
      +
      ABCDEFGHIJKLMNOPQR
      """
    And a file "tags.tsv" containing:
      """
      id	seq	Cell	note
      w1	AACCGG	A+B	first well
      w2	TTGGCC		second well
      w3	CATCAT	C,D:1=x~y;@z	third well
      """

  Scenario Outline: a column alone, joined with streams, and beside the matched barcode
    When I run `unmux reads.fq --pool reads --group bc=tags.tsv --group bc::loc=0:0:6,dist=1 --extract well=@bc --extract cell=@bc::cell --extract bc2=0:6:10 --extract r=0:10:end --tag BC=well --tag CB=cell+bc2 --tag xc=cell+bc2+bc2 --tag xd=cell --template r --threads <threads> --out out.<ext>`
    Then the exit code is 0
    And the records of "out.<ext>" are:
      """
      r1	4	*	0	255	*	*	0	0	AAAACCCC	KLMNOPQR	BC:Z:AACCGG	QT:Z:ABCDEF	CB:Z:A+B-ACGT	xc:Z:A+B-ACGT-ACGT	xd:Z:A+B	RG:Z:reads
      r2	4	*	0	255	*	*	0	0	AAAACCCC	KLMNOPQR	BC:Z:AACCGG	QT:Z:ABCDEF	CB:Z:A+B-ACGT	xc:Z:A+B-ACGT-ACGT	xd:Z:A+B	RG:Z:reads
      r3	4	*	0	255	*	*	0	0	AAAACCCC	KLMNOPQR	BC:Z:TTGGCC	QT:Z:ABCDEF	CB:Z:TTGGCC-TTTT	xc:Z:TTGGCC-TTTT-TTTT	xd:Z:TTGGCC	RG:Z:reads
      r4	4	*	0	255	*	*	0	0	AAAACCCC	KLMNOPQR	CB:Z:TTTT	xc:Z:TTTT-TTTT	RG:Z:reads
      """

    Examples:
      | ext | threads |
      | sam | 1       |
      | sam | 2       |
      | sam | 8       |
      | bam | 1       |
      | bam | 2       |
      | bam | 8       |

  Scenario Outline: column text survives CRAM output
    When I run `unmux reads.fq --group bc=tags.tsv --group bc::loc=0:0:6,dist=1 --extract cell=@bc::cell --extract bc2=0:6:10 --extract r=0:10:end --tag xc=cell+bc2 --template r --threads <threads> --out out.cram`
    Then the exit code is 0
    When I run `unmux out.cram --out back.sam`
    Then the exit code is 0
    And the file "back.sam" contains "xc:Z:A+B-ACGT"
    And the file "back.sam" contains "xc:Z:TTGGCC-TTTT"

    Examples:
      | threads |
      | 1       |
      | 8       |

  Scenario Outline: column text goes in the FASTQ comment
    When I run `unmux reads.fq --group bc=tags.tsv --group bc::loc=0:0:6,dist=1 --extract cell=@bc::cell --extract bc2=0:6:10 --extract r=0:10:end --tag xc=cell+bc2 --template r --threads <threads> --out out.fq`
    Then the exit code is 0
    And the file "out.fq" is:
      """
      @r1 xc:Z:A+B-ACGT
      AAAACCCC
      +
      KLMNOPQR
      @r2 xc:Z:A+B-ACGT
      AAAACCCC
      +
      KLMNOPQR
      @r3 xc:Z:TTGGCC-TTTT
      AAAACCCC
      +
      KLMNOPQR
      @r4 xc:Z:TTTT
      AAAACCCC
      +
      KLMNOPQR
      """

    Examples:
      | threads |
      | 1       |
      | 2       |
      | 8       |

  Scenario: column text goes in the FASTA comment
    When I run `unmux reads.fq --group bc=tags.tsv --group bc::loc=0:0:6 --group bc::minFindsPerGroup=1 --extract cell=@bc::cell --extract r=0:10:end --tag xc=cell --template r --out out.fa`
    Then the exit code is 0
    And the file "out.fa" is:
      """
      >r1 xc:Z:A+B
      AAAACCCC
      >r3 xc:Z:TTGGCC
      AAAACCCC
      """

  Scenario: a value may hold any printable text
    Given a file "one.fq" containing:
      """
      @r1
      CATCATAA
      +
      IIIIIIII
      """
    When I run `unmux one.fq --group bc=tags.tsv --group bc::loc=0:0:6 --extract cell=@bc::cell --tag xc=cell`
    Then the exit code is 0
    And stdout contains "xc:Z:C,D:1=x~y;@z"

  Scenario: the tag ID of a file set, by any ID header
    Given a file "barcodes.tsv" containing:
      """
      name	barcode
      w1	AACCGG
      """
    When I run `unmux reads.fq --group bc=barcodes.tsv --group bc::loc=0:0:6 --extract a=@bc::id --extract b=@bc::name --extract c=@bc::barcode --tag xa=a --tag xb=b --tag BC=c`
    Then the exit code is 0
    And stdout contains "@r1 xa:Z:w1	xb:Z:w1	BC:Z:AACCGG	QT:Z:ABCDEF"

  Scenario: the tag ID of an inline set
    When I run `unmux reads.fq --group bc={w1=AACCGG,TTGGCC,=GAGAGA} --group bc::loc=0:0:6 --extract a=@bc::id --tag xa=a`
    Then the exit code is 0
    And stdout contains "@r1 xa:Z:w1"
    And stdout contains "@r3 xa:Z:TTGGCC"
    And stdout contains "@r4 xa:Z:GAGAGA"

  Scenario: an inline set has no other column
    When I run `unmux reads.fq --group bc={w1=AACCGG} --extract a=@bc::cell --tag xa=a --out out.sam`
    Then the exit code is 1
    And stderr contains "extract `a` reads column `cell` of group `bc`, but an inline set has only `id` and `seq`"
    And a file "out.sam" does not exist

  Scenario Outline: an inline tag ID must be printable ASCII to be read
    When I run `unmux reads.fq --group bc=<source> --extract a=@bc::id --tag xa=a`
    Then the exit code is <code>

    Examples:
      | source      | code |
      | {é=AACCGG}  | 1    |
      | {e=AACCGG}  | 0    |

  Scenario Outline: the seq column is the matched span with its qualities
    When I run `unmux reads.fq --group bc=tags.tsv --group bc::loc=0:0:6,dist=1 --extract s=<span> --tag BC=s --tag OX=s --tag OX::raw=true --tag xr=~s --template s --unassigned un.fq --out out.fq`
    Then the exit code is 0
    And the file "out.fq" is:
      """
      @r1 BC:Z:AACCGG	QT:Z:ABCDEF	OX:Z:AACCGG	BZ:Z:ABCDEF	xr:Z:CCGGTT
      AACCGG
      +
      ABCDEF
      @r2 BC:Z:AACCGG	QT:Z:ABCDEF	OX:Z:AACCGT	BZ:Z:ABCDEF	xr:Z:CCGGTT
      AACCGG
      +
      ABCDEF
      @r3 BC:Z:TTGGCC	QT:Z:ABCDEF	OX:Z:TTGGCC	BZ:Z:ABCDEF	xr:Z:GGCCAA
      TTGGCC
      +
      ABCDEF
      """

    Examples:
      | span         |
      | @bc          |
      | @bc::seq     |
      | @bc::Barcode |

  Scenario: a blank value is the corrected sequence of a reverse-strand match
    Given a file "one.fq" containing:
      """
      @r1
      GGCCAT
      +
      IIIIII
      """
    When I run `unmux one.fq --group bc=tags.tsv --group bc::bothStrands=true,dist=1 --extract b=@bc --extract cell=@bc::cell --tag xb=b --tag xc=cell --tag xr=cell --tag xr::raw=true`
    Then the exit code is 0
    And stdout contains "xb:Z:TTGGCC	xc:Z:TTGGCC	xr:Z:TTGGCC"

  Scenario: raw=true does not change column text
    When I run `unmux reads.fq --group bc=tags.tsv --group bc::loc=0:0:6,dist=1 --extract b=@bc --extract cell=@bc::cell --tag xc=cell+b --tag xc::raw=true`
    Then the exit code is 0
    And stdout contains "@r2 xc:Z:A+B-AACCGT"

  Scenario: a column tag has no default quality tag and accepts qual=none
    When I run `unmux reads.fq --group bc=tags.tsv --group bc::loc=0:0:6 --extract cell=@bc::cell --extract u=0:6:10 --tag CB=cell --tag RX=cell+u --tag BC=u+cell --tag BC::qual=none`
    Then the exit code is 0
    And stdout contains "@r1 CB:Z:A+B	RX:Z:A+B-ACGT	BC:Z:ACGT-A+B"
    And stdout contains "@r4 RX:Z:TTTT	BC:Z:TTTT"
    And stdout does not contain "CY:Z:"
    And stdout does not contain "QX:Z:"
    And stdout does not contain "QT:Z:"

  Scenario Outline: a column stream outside --tag, under ~, or with a quality tag is an error
    When I run `unmux reads.fq --group bc=tags.tsv --extract cell=@bc::cell --extract u=0:6:10 <use> --out out.sam`
    Then the exit code is 1
    And stderr contains "<message>"
    And a file "out.sam" does not exist

    Examples:
      | use                                                | message                                            |
      | --template cell                                    | template #1 uses column stream `cell`              |
      | --template u cell+u                                | template #2 uses column stream `cell`              |
      | --template cell::raw=true                          | template #1 uses column stream `cell`              |
      | --group g={ACGT} --group g::match=u+cell           | group `g` match= uses column stream `cell`         |
      | --tag xc=~cell                                     | tag `xc` reverse-complements column stream `cell`  |
      | --tag xc=u+~cell                                   | tag `xc` reverse-complements column stream `cell`  |
      | --tag xc=cell --tag xc::qual=xq                    | tag `xc` sets `qual=xq`, but its column stream     |
      | --tag CB=u+cell --tag CB::qual=CY                  | tag `CB` sets `qual=CY`, but its column stream     |

  Scenario: an unknown column is an error that lists the file's columns
    When I run `unmux reads.fq --group bc=tags.tsv --extract cell=@bc::nope --tag xc=cell --out out.sam --metrics-summary m.tsv`
    Then the exit code is 1
    And stderr contains "tags.tsv: no column `nope`; the columns are: id, seq, Cell, note"
    And a file "out.sam" does not exist
    And a file "m.tsv" does not exist

  Scenario: a value that is not printable ASCII is an error naming its line
    Given a file "spaced.tsv" containing:
      """
      # a comment
      id	seq	cell

      w1	AACCGG	ok
      w2	TTGGCC	é
      """
    When I run `unmux reads.fq --group bc=spaced.tsv --extract cell=@bc::cell --tag xc=cell --out out.sam`
    Then the exit code is 1
    And stderr contains "spaced.tsv: line 5: value "é" in column `cell`"
    And a file "out.sam" does not exist

  Scenario Outline: a malformed column span is an error
    When I run `unmux reads.fq --group bc=tags.tsv --extract x=<span> --tag xc=x`
    Then the exit code is 1
    And stderr contains "<message>"

    Examples:
      | span            | message                              |
      | @bc::           | a column span is `@group::COLUMN`    |
      | @bc::cell::x    | a column span is `@group::COLUMN`    |
      | @bc::cell+1:3   | a column span is `@group::COLUMN`    |
      | @bc::cell-1:3   | a column span is `@group::COLUMN`    |
      | @bc+1:3::cell   | a column span is `@group::COLUMN`    |
      | @bc::cell..@bc  | a column span is `@group::COLUMN`    |
      | @bc..@bc::cell  | a column span is `@group::COLUMN`    |
      | @nope::cell     | references undefined group `nope`    |

  Scenario: columns that no extract reads are ignored
    Given a file "extra.tsv" containing:
      """
      note	id	seq	note	
      two words	w1	AACCGG	x	y	z
      é	w2	TTGGCC
      """
    When I run `unmux reads.fq --group bc=extra.tsv --group bc::loc=0:0:6 --extract b=@bc --tag xb=b`
    Then the exit code is 0
    And stdout contains "@r1 xb:Z:AACCGG"
    And stdout contains "@r3 xb:Z:TTGGCC"

  Scenario: a column header that repeats is an error only when read
    Given a file "extra.tsv" containing:
      """
      id	seq	cell	CELL
      w1	AACCGG	x	y
      """
    When I run `unmux reads.fq --group bc=extra.tsv --extract cell=@bc::Cell --tag xc=cell`
    Then the exit code is 1
    And stderr contains "more than one column `Cell`"

  Scenario: headers that differ only by case are read by their exact spelling
    Given a file "extra.tsv" containing:
      """
      id	seq	cell	Cell
      w1	AACCGG	x	y
      """
    When I run `unmux reads.fq --group bc=extra.tsv --group bc::loc=0:0:6 --extract a=@bc::cell --extract b=@bc::Cell --tag xa=a --tag xb=b`
    Then the exit code is 0
    And stdout contains "@r1 xa:Z:x	xb:Z:y"

  Scenario: a header with spaces or punctuation is percent-escaped
    Given a file "extra.tsv" containing:
      """
      id	seq	cell type	a+b::c	µ
      w1	AACCGG	x	y	z
      """
    When I run `unmux reads.fq --group bc=extra.tsv --group bc::loc=0:0:6 --extract a=@bc::cell%20type --extract b=@bc::A%2Bb%3A%3Ac --extract c=@bc::%C2%B5 --tag xa=a+b+c`
    Then the exit code is 0
    And stdout contains "@r1 xa:Z:x-y-z"

  Scenario Outline: a value may hold interior spaces and surrounding spaces are trimmed
    Given a file "spaced.tsv" containing:
      """
      id	seq	cell	note
      w 1	AACCGG	  two  words  	x
      """
    When I run `unmux reads.fq --group bc=spaced.tsv --group bc::loc=0:0:6 --extract cell=@bc::cell --extract i=@bc::id --tag xc=cell --tag xi=i --out out.<ext>`
    Then the exit code is 0
    When I run `unmux out.<ext> --out back.sam`
    Then the exit code is 0
    And the file "back.sam" contains "xc:Z:two  words	xi:Z:w 1"

    Examples:
      | ext  |
      | sam  |
      | bam  |
      | cram |
      | fq   |
      | fa   |

  Scenario Outline: a match= group has columns
    Given a file "pairs.tsv" containing:
      """
      id	seq	cell
      p1	AACCACGT	A+B
      p2	TTGGTTTT	
      """
    When I run `unmux reads.fq --extract x=0:0:4 --extract y=0:6:10 --group pair=pairs.tsv --group pair::match=x+y --extract cell=@pair::cell --extract i=@pair::id --extract r=0:10:end --tag xc=cell+x --tag xi=i --template r --threads <threads> --out out.fa`
    Then the exit code is 0
    And the file "out.fa" is:
      """
      >r1 xc:Z:A+B-AACC	xi:Z:p1
      AAAACCCC
      >r2 xc:Z:A+B-AACC	xi:Z:p1
      AAAACCCC
      >r3 xc:Z:TTGGTTTT-TTGG	xi:Z:p2
      AAAACCCC
      >r4 xc:Z:GAGA
      AAAACCCC
      """

    Examples:
      | threads |
      | 1       |
      | 2       |
      | 8       |

  Scenario Outline: a carried default quality tag is dropped when a column tag replaces its tag
    Given a file "carried.sam" containing:
      """
      @HD	VN:1.6
      c1	4	*	0	0	*	*	0	0	AACCGGACGT	IIIIIIIIII	CB:Z:OLD	CY:Z:III	XX:Z:kept
      c2	4	*	0	0	*	*	0	0	AACCGGACGT	IIIIIIIIII	CR:Z:OLD	CB:Z:OLD	CY:Z:III
      c3	4	*	0	0	*	*	0	0	GAGAGAACGT	IIIIIIIIII	CB:Z:OLD	CY:Z:III
      """
    When I run `unmux carried.sam --pool p --group bc=tags.tsv --group bc::loc=0:0:6 --extract a=@bc::cell --extract u=0:6:10 --tag CB=a <qual> --threads <threads> --out out.<ext>`
    Then the exit code is 0
    And the records of "out.<ext>" are:
      """
      c1	4	*	0	255	*	*	0	0	AACCGGACGT	IIIIIIIIII	CB:Z:A+B	XX:Z:kept	RG:Z:p
      c2	4	*	0	255	*	*	0	0	AACCGGACGT	IIIIIIIIII	CR:Z:OLD	CB:Z:A+B	CY:Z:III	RG:Z:p
      c3	4	*	0	255	*	*	0	0	GAGAGAACGT	IIIIIIIIII	CB:Z:OLD	CY:Z:III	RG:Z:p
      """

    Examples:
      | qual                   | ext | threads |
      | --per-record           | sam | 1       |
      | --per-record           | bam | 8       |
      | --tag CB::qual=none    | sam | 8       |

  Scenario: a default quality tag left out for a column is warned about
    When I run `unmux reads.fq --group bc=tags.tsv --extract cell=@bc::cell --extract u=0:6:10 --tag CB=cell+u`
    Then the exit code is 0
    And stderr contains "tag `CB` has a column stream, which has no base qualities, so its default quality tag `CY` is left out"
    When I run `unmux reads.fq --group bc=tags.tsv --extract cell=@bc::cell --extract u=0:6:10 --tag CB=cell+u --tag CB::qual=none --tag xc=cell`
    Then the exit code is 0
    And stderr does not contain "is left out"

  Scenario: routing, bins, metrics and QC are the same with and without a column stream
    When I run `unmux reads.fq --pool p --group bc=tags.tsv --group bc::loc=0:0:6,dist=1 --extract b=@bc --extract r=0:10:end --tag BC=b --template r --sample-from-group bc --qc-tag --remove bc::w2=rm.fq --unassigned un.fq --out a.%sample.fq --metrics-summary a.summary.tsv --metrics-per-sample a.samples.tsv --metrics-unassigned a.unassigned.tsv`
    Then the exit code is 0
    When I run `unmux reads.fq --pool p --group bc=tags.tsv --group bc::loc=0:0:6,dist=1 --extract b=@bc --extract r=0:10:end --extract cell=@bc::cell --tag BC=b --tag xc=cell --template r --sample-from-group bc --qc-tag --remove bc::w2=rm2.fq --unassigned un2.fq --out b.%sample.fq --metrics-summary b.summary.tsv --metrics-per-sample b.samples.tsv --metrics-unassigned b.unassigned.tsv`
    Then the exit code is 0
    And the files "a.summary.tsv" and "b.summary.tsv" are identical
    And the files "a.samples.tsv" and "b.samples.tsv" are identical
    And the files "a.unassigned.tsv" and "b.unassigned.tsv" are identical
    And the files "rm.fq" and "rm2.fq" are identical
    And the files "un.fq" and "un2.fq" are identical
    And the file "b.w1.fq" contains "xc:Z:A+B"
    And the file "a.w1.fq" contains "ZS:Z:{"v":1,"outcome":"assigned","sample":"w1","sub_sample":null,"groups":[{"g":"bc","tag":"AACCGG","loc":"0:0:6","obs":"AACCGT","sub":1,"ind":0}]}"
    And the file "b.w1.fq" contains "ZS:Z:{"v":1,"outcome":"assigned","sample":"w1","sub_sample":null,"groups":[{"g":"bc","tag":"AACCGG","loc":"0:0:6","obs":"AACCGT","sub":1,"ind":0}]}"

  Scenario: the sub_sample column still routes and can also be read
    Given a file "subs.tsv" containing:
      """
      id	seq	sub_sample
      w1	AACCGG	lib1
      w2	TTGGCC
      """
    When I run `unmux reads.fq --group bc=subs.tsv --group bc::loc=0:0:6 --extract s=@bc::sub_sample --tag xs=s --sample-from-group bc --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "xs:Z:lib1	RG:Z:w1.lib1"
    And the file "out.sam" contains "xs:Z:TTGGCC	RG:Z:w2"

  Scenario: a read group has no tag columns
    Given the test file "aligned.sam"
    When I run `unmux aligned.sam --group rg=@RG --sample-from-group rg --extract a=@rg::id --tag xa=a --out %sample.sam`
    Then the exit code is 1
    And stderr contains "an `@RG` group has no tag columns"

  Scenario: an aligned read-group split refuses a column tag
    Given the test file "aligned.sam"
    When I run `unmux aligned.sam --group rg=@RG --sample-from-group rg --group bc={w1=AACCGG} --extract a=@bc::id --tag xa=a --out %sample.sam`
    Then the exit code is 1
    And stderr contains "group `bc` matches record sequences"
