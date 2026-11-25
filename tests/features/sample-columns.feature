Feature: Sample-sheet column streams
  The extract source `%sample::COLUMN` is the text in a `--sample-sheet` column of the row that
  routed the record, for use in `--tag` alone or joined with other streams. Each row has its own
  value, so two rows of one sample can give different tags. A blank cell leaves out every tag built
  from the stream, even one that joins it with barcodes, while the record's other tags are written.
  The reserved columns `sample` and `sub_sample` are the routed sample and sub-sample IDs with any
  sample flag; any other column needs `--sample-sheet`. Column text has no base qualities, so it
  cannot go in a `--template`, a group `match=`, or under `~`, and a tag that uses it has no quality
  tag. Unassigned and removed records get no `--tag`, as before.

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
      @r5
      GAGAGGGGACGTACGT
      +
      ABCDEFGHIJKLMNOP
      """
    And a file "sheet.tsv" containing:
      """
      sample	a	b	nucleus
      s1	a1	b1	N1
      s1	a1	b2	N2+N3
      s2	a2	b1	
      s2	a2	b2	N4
      """

  Scenario Outline: each sheet row gives its own value, and a blank cell leaves out the tags built from it
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract ba=@a --extract bb=@b --extract r=0:8:end --tag xn=nuc --tag nu=nuc+ba+bb --tag CB=ba+bb --template r --sample-sheet sheet.tsv --threads <threads> --out out.<ext>`
    Then the exit code is 0
    And the records of "out.<ext>" are:
      """
      r1	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	xn:Z:N1	nu:Z:N1-AAAA-GGGG	CB:Z:AAAA-GGGG	CY:Z:ABCD EFGH	RG:Z:s1
      r2	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	xn:Z:N2+N3	nu:Z:N2+N3-AAAA-TTTT	CB:Z:AAAA-TTTT	CY:Z:ABCD EFGH	RG:Z:s1
      r3	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	CB:Z:CCCC-GGGG	CY:Z:ABCD EFGH	RG:Z:s2
      r4	4	*	0	255	*	*	0	0	ACGTACGT	IJKLMNOP	xn:Z:N4	nu:Z:N4-CCCC-TTTT	CB:Z:CCCC-TTTT	CY:Z:ABCD EFGH	RG:Z:s2
      """
    And stderr contains "column "nucleus" is blank in 1 of 4 sheet rows; records they route get no xn or nu tag"

    Examples:
      | ext | threads |
      | sam | 1       |
      | sam | 2       |
      | sam | 8       |
      | bam | 1       |
      | bam | 2       |
      | bam | 8       |

  Scenario Outline: sample-sheet column text survives CRAM output
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract ba=@a --extract bb=@b --extract r=0:8:end --tag nu=nuc+ba+bb --template r --sample-sheet sheet.tsv --threads <threads> --out out.cram`
    Then the exit code is 0
    When I run `unmux out.cram --out back.sam`
    Then the exit code is 0
    And the file "back.sam" contains "nu:Z:N1-AAAA-GGGG"
    And the file "back.sam" contains "nu:Z:N2+N3-AAAA-TTTT"
    And the file "back.sam" contains "nu:Z:N4-CCCC-TTTT"
    And the file "back.sam" does not contain "nu:Z:CCCC-GGGG"

    Examples:
      | threads |
      | 1       |
      | 8       |

  Scenario Outline: sample-sheet column text goes in the FASTQ comment
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract bb=@b --extract r=0:8:end --tag nu=nuc+bb --tag xb=bb --template r --sample-sheet sheet.tsv --threads <threads> --out out.fq`
    Then the exit code is 0
    And the file "out.fq" is:
      """
      @r1 nu:Z:N1-GGGG	xb:Z:GGGG
      ACGTACGT
      +
      IJKLMNOP
      @r2 nu:Z:N2+N3-TTTT	xb:Z:TTTT
      ACGTACGT
      +
      IJKLMNOP
      @r3 xb:Z:GGGG
      ACGTACGT
      +
      IJKLMNOP
      @r4 nu:Z:N4-TTTT	xb:Z:TTTT
      ACGTACGT
      +
      IJKLMNOP
      """

    Examples:
      | threads |
      | 1       |
      | 2       |
      | 8       |

  Scenario: two extracts may read one column
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract n=%sample::nucleus --extract m=%sample::nucleus --extract s=%sample::sample --extract r=0:8:end --tag xn=n --tag xs=s --tag xm=m --template r --sample-sheet sheet.tsv --out out.fq`
    Then the exit code is 0
    And the file "out.fq" is:
      """
      @r1 xn:Z:N1	xs:Z:s1	xm:Z:N1
      ACGTACGT
      +
      IJKLMNOP
      @r2 xn:Z:N2+N3	xs:Z:s1	xm:Z:N2+N3
      ACGTACGT
      +
      IJKLMNOP
      @r3 xs:Z:s2
      ACGTACGT
      +
      IJKLMNOP
      @r4 xn:Z:N4	xs:Z:s2	xm:Z:N4
      ACGTACGT
      +
      IJKLMNOP
      """

  Scenario: unassigned and removed records get no tag
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract sm=%sample::sample --tag nu=nuc --tag xs=sm --sample-sheet sheet.tsv --remove a::a2=rm.fq --unassigned un.fq --out out.fq`
    Then the exit code is 0
    And the file "out.fq" is:
      """
      @r1 nu:Z:N1	xs:Z:s1
      AAAAGGGGACGTACGT
      +
      ABCDEFGHIJKLMNOP
      @r2 nu:Z:N2+N3	xs:Z:s1
      AAAATTTTACGTACGT
      +
      ABCDEFGHIJKLMNOP
      """
    And the file "rm.fq" is:
      """
      @r3
      CCCCGGGGACGTACGT
      +
      ABCDEFGHIJKLMNOP
      @r4
      CCCCTTTTACGTACGT
      +
      ABCDEFGHIJKLMNOP
      """
    And the file "un.fq" is:
      """
      @r5
      GAGAGGGGACGTACGT
      +
      ABCDEFGHIJKLMNOP
      """

  Scenario: the sample and sub_sample columns of a sample sheet
    Given a file "subs.tsv" containing:
      """
      sample	sub_sample	a
      s1	%pool	a1
      s2		a2
      """
    When I run `unmux reads.fq --pool p --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --extract sm=%sample::sample --extract ss=%sample::sub_sample --tag xs=sm --tag xb=ss --sample-sheet subs.tsv --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "xs:Z:s1	xb:Z:p	RG:Z:s1.p"
    And the file "out.sam" contains "xs:Z:s2	RG:Z:s2"
    And stderr contains "column "sub_sample" is blank in 1 of 2 sheet rows; records they route get no xb tag"

  Scenario: a sample sheet without a sub_sample column has a blank sub_sample
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract ss=%sample::sub_sample --extract sm=%sample::SAMPLE --tag xb=ss --tag xs=sm --sample-sheet sheet.tsv --sample-sheet-ignore nucleus --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "xs:Z:s1	RG:Z:s1"
    And the file "out.sam" does not contain "xb:Z:"

  Scenario: the sample and sub_sample columns of inline samples
    When I run `unmux reads.fq --pool p --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --extract sm=%sample::sample --extract ss=%sample::sub_sample --tag xs=sm --tag xb=ss --sample s1::%pool=a::a1 --sample s2=a::a2 --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "r1	4	*	0	255	*	*	0	0	AAAAGGGGACGTACGT	ABCDEFGHIJKLMNOP	xs:Z:s1	xb:Z:p	RG:Z:s1.p"
    And the file "out.sam" contains "r3	4	*	0	255	*	*	0	0	CCCCGGGGACGTACGT	ABCDEFGHIJKLMNOP	xs:Z:s2	RG:Z:s2"
    And stderr contains "column "sub_sample" is blank in 1 of 2 samples; records they route get no xb tag"

  Scenario: the sample and sub_sample columns with --sample-from-group
    Given a file "tags.tsv" containing:
      """
      id	seq	sub_sample
      w1	AAAA	lib1
      w2	CCCC
      """
    When I run `unmux reads.fq --group a=tags.tsv --group a::loc=0:0:4 --extract sm=%sample::sample --extract ss=%sample::sub_sample --tag xs=sm --tag xb=ss --sample-from-group a --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "xs:Z:w1	xb:Z:lib1	RG:Z:w1.lib1"
    And the file "out.sam" contains "xs:Z:w2	RG:Z:w2"
    And stderr contains "column "sub_sample" is blank in 1 of 2 tags of group `a`; records they route get no xb tag"

  Scenario Outline: a column other than sample or sub_sample needs a sample sheet
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --extract nuc=%sample::<column> --tag nu=nuc <samples> --out out.sam`
    Then the exit code is 1
    And stderr contains "<message>"
    And a file "out.sam" does not exist

    Examples:
      | column  | samples                | message                                                                                                                |
      | nucleus | --sample s1=a::a1      | extract `nuc` reads sample-sheet column "nucleus", but no `--sample-sheet` is given (only `sample` and `sub_sample` work with `--sample`) |
      | nucleus | --sample-from-group a  | extract `nuc` reads sample-sheet column "nucleus"; with `--sample-from-group a` the row is a tag, so use `@a::nucleus` |
      | a%20b   | --sample-from-group a  | so use `@a::a%20b`                                                                                                     |
      | nucleus |                        | `--extract` uses `%sample::nucleus` but no samples are defined; add `--sample-sheet` (or remove the extract)          |
      | sample  |                        | `--extract` uses `%sample::sample` but no samples are defined; add `--sample`, `--sample-sheet`, or `--sample-from-group` |

  Scenario Outline: a malformed sample-sheet column span is an error
    When I run `unmux reads.fq --group a={a1=AAAA} --sample-sheet sheet.tsv --extract x=<span> --tag xc=x`
    Then the exit code is 1
    And stderr contains "<message>"

    Examples:
      | span                  | message                                                                                    |
      | %sample               | `%sample` is not an `--extract` source; use `%sample::sample` for the sample ID            |
      | %sub_sample           | `%sub_sample` is not an `--extract` source; use `%sample::sub_sample`                      |
      | %pool                 | `%pool` is not an `--extract` source; use `%sample::COLUMN` (or `--rg-tag` for read-group fields) |
      | %ordinal              | `%ordinal` is not an `--extract` source; use `%sample::COLUMN`                             |
      | %source               | `%source` is not an `--extract` source; use `%sample::COLUMN`                              |
      | %sample::             | a sample-sheet column span is `%sample::COLUMN` with no offset                             |
      | %sample::nucleus+1:4  | a sample-sheet column span is `%sample::COLUMN` with no offset                             |
      | %sample::a..b         | a sample-sheet column span is `%sample::COLUMN` with no offset                             |
      | %sample::nucleus::x   | a sample-sheet column span is `%sample::COLUMN` with no offset                             |

  Scenario Outline: a sample-sheet column stream outside --tag, under ~, or with a quality tag is an error
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group b={b1=GGGG,b2=TTTT} --sample-sheet sheet.tsv --extract nuc=%sample::nucleus --extract u=0:8:12 <use> --out out.sam`
    Then the exit code is 1
    And stderr contains "<message>"
    And a file "out.sam" does not exist

    Examples:
      | use                                      | message                                           |
      | --template nuc                           | template #1 uses column stream `nuc`              |
      | --group g={ACGT} --group g::match=u+nuc  | group `g` match= uses column stream `nuc`         |
      | --tag xn=~nuc                            | tag `xn` reverse-complements column stream `nuc`  |
      | --tag CB=u+nuc --tag CB::qual=CY         | tag `CB` sets `qual=CY`, but its column stream    |

  Scenario: a tag built from a sample-sheet column has no default quality tag
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --extract u=0:8:12 --tag CB=nuc+u --sample-sheet sheet.tsv --out out.sam`
    Then the exit code is 0
    And stderr contains "tag `CB` has a column stream, which has no base qualities, so its default quality tag `CY` is left out"
    And the file "out.sam" contains "CB:Z:N1-ACGT	RG:Z:s1"
    And the file "out.sam" does not contain "CY:Z:"

  Scenario Outline: two rows that can route one record to one sample must agree on a read column
    Given a file "overlap.tsv" containing:
      """
      sample	a	b	nucleus
      <first>
      <second>
      <third>
      """
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --tag nu=nuc --sample-sheet overlap.tsv --out out.sam`
    Then the exit code is 1
    And stderr contains "<message>"
    And a file "out.sam" does not exist

    Examples:
      | first             | second            | third             | message                                                                                                    |
      | s1	a1	b1	N1  | s1	a1	b1	N2  |                   | sample sheet rows 1 and 2 can both route a record to `s1` but disagree on "nucleus" ("N1" vs "N2")          |
      | s2	a2	b2	N4  | s1	a1	b1	N1  | s1	a1		N2   | sample sheet rows 2 and 3 can both route a record to `s1` but disagree on "nucleus" ("N1" vs "N2")          |
      | s1	a1		N1   | s1	AAAA	b2    |                   | sample sheet rows 1 and 2 can both route a record to `s1` but disagree on "nucleus" ("N1" vs a blank cell)  |
      | s1		b1	N1   | s1	a1		N2   |                   | sample sheet rows 1 and 2 can both route a record to `s1` but disagree on "nucleus" ("N1" vs "N2")          |

  Scenario: rows that agree, route to different samples, or select different tags are not in conflict
    Given a file "agree.tsv" containing:
      """
      sample	a	b	nucleus
      s1	a1	b1	N1
      s1	a1	b1	N1
      s1	a1		N1
      s1	a2	b2	N9
      s2	a2	b1	N1
      """
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --tag nu=nuc --sample-sheet agree.tsv --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "nu:Z:N1	RG:Z:s1"
    And the file "out.sam" contains "nu:Z:N1	RG:Z:s2"
    And the file "out.sam" contains "nu:Z:N9	RG:Z:s1"

  Scenario Outline: a read column is found by its header
    Given a file "headers.tsv" containing:
      """
      sample	a	<header>
      s1	a1	x
      """
    When I run `unmux reads.fq --group a={a1=AAAA} --group a::loc=0:0:4 --extract v=%sample::<column> --tag xv=v --sample-sheet headers.tsv`
    Then the exit code is 0
    And stdout contains "@r1 xv:Z:x"

    Examples:
      | header     | column      |
      | nucleus    | nucleus     |
      | nucleus    | NUCLEUS     |
      | same cell  | same%20cell |
      | a+b        | a%2Bb       |
      | µ          | %C2%B5      |

  Scenario: headers that differ only by case are read by their exact spelling
    Given a file "headers.tsv" containing:
      """
      sample	a	Nucleus	nucleus
      s1	a1	x	y
      """
    When I run `unmux reads.fq --group a={a1=AAAA} --group a::loc=0:0:4 --extract u=%sample::Nucleus --extract v=%sample::nucleus --tag xu=u --tag xv=v --sample-sheet headers.tsv`
    Then the exit code is 0
    And stdout contains "@r1 xu:Z:x	xv:Z:y"

  Scenario Outline: a read column that is missing, ambiguous, or a group is an error
    Given a file "headers.tsv" containing:
      """
      sample	a	Nucleus	NUCLEUS
      s1	a1	x	y
      """
    When I run `unmux reads.fq --group a={a1=AAAA} --extract v=%sample::<column> --tag xv=v --sample-sheet headers.tsv --out out.sam`
    Then the exit code is 1
    And stderr contains "<message>"
    And a file "out.sam" does not exist

    Examples:
      | column   | message                                                                                                    |
      | nucleu   | headers.tsv: extract `v` reads sample-sheet column "nucleu"; the columns are: sample, a, Nucleus, NUCLEUS   |
      | nucleus  | extract `v` reads sample-sheet column "nucleus", but more than one column matches it regardless of case     |
      | a        | extract `v` reads sample-sheet column "a", a group column; use `@a::id` for the matched tag                |

  Scenario: a read column need not be a group, but other columns still must
    Given a file "extra.tsv" containing:
      """
      sample	a	nucleus	note
      s1	a1	x	y
      """
    When I run `unmux reads.fq --group a={a1=AAAA} --extract v=%sample::nucleus --tag xv=v --sample-sheet extra.tsv`
    Then the exit code is 1
    And stderr contains "sample sheet column `note` is not a declared --group"
    When I run `unmux reads.fq --group a={a1=AAAA} --group a::loc=0:0:4 --extract v=%sample::nucleus --tag xv=v --sample-sheet extra.tsv --sample-sheet-ignore note`
    Then the exit code is 0
    And stdout contains "@r1 xv:Z:x"

  Scenario: a read column that is also ignored is read, with a warning
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract nuc=%sample::nucleus --tag nu=nuc --sample-sheet sheet.tsv --sample-sheet-ignore nucleus`
    Then the exit code is 0
    And stderr contains "--sample-sheet-ignore skips "nucleus", which `%sample::nucleus` reads; the column is read"
    And stdout contains "@r1 nu:Z:N1"

  Scenario: values keep their case, leading zeros and punctuation, and are trimmed of spaces
    Given a file "values.tsv" containing:
      """
      sample	a	b	v
      s1	a1	b1	  007 Ab  
      s1	a1	b2	x,y;z+w:1
      s2	a2	b1	   
      s2	a2	b2	~!@#$%^&*()
      """
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --group a::loc=0:0:4 --group b={b1=GGGG,b2=TTTT} --group b::loc=0:4:8 --extract v=%sample::v --tag xv=v --sample-sheet values.tsv --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "xv:Z:007 Ab	RG:Z:s1"
    And the file "out.sam" contains "xv:Z:x,y;z+w:1	RG:Z:s1"
    And the file "out.sam" contains "xv:Z:~!@#$%^&*()	RG:Z:s2"
    And the file "out.sam" contains "r3	4	*	0	255	*	*	0	0	CCCCGGGGACGTACGT	ABCDEFGHIJKLMNOP	RG:Z:s2"

  Scenario: a value that is not printable ASCII is an error that gives its row
    Given a file "values.tsv" containing:
      """
      sample	a	v
      # a comment
      s1	a1	ok

      s2	a2	é
      """
    When I run `unmux reads.fq --group a={a1=AAAA,a2=CCCC} --extract v=%sample::v --tag xv=v --sample-sheet values.tsv --out out.sam`
    Then the exit code is 1
    And stderr contains "values.tsv: sample sheet row 2: value "é" in column "v" may hold only printable ASCII and spaces"
    And a file "out.sam" does not exist

  Scenario: a sample ID that is not printable ASCII cannot be read
    Given a file "tags.tsv" containing:
      """
      id	seq
      é	AAAA
      """
    When I run `unmux reads.fq --group a=tags.tsv --extract sm=%sample::sample --tag xs=sm --sample-from-group a --out out.sam`
    Then the exit code is 1
    And stderr contains "extract `sm` reads "sample", but value "é" may hold only printable ASCII and spaces"
