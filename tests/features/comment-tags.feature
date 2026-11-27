Feature: Lifting SAM tags from FASTX read-name comments
  A FASTX read-name comment may carry SAM auxiliary tags (the samtools `fastq -T`
  convention that unmux itself writes on FASTX output). On input those tags are
  lifted back into real SAM tags, the inverse of that output, so a tag-bearing
  FASTQ converts to a uBAM/SAM with real tags and no extraction configured.
  Alignment tags such as NM and MD are dropped, as they are from SAM/BAM/CRAM
  input, while Dorado's `pa:B` poly(A) array is kept on both.
  Every type unmux writes reads back: `A`, `f`, `Z`, `H`, `B` arrays of every
  subtype, and `i` from -2147483648 to 4294967295 in the narrowest BAM type that
  holds it. Lifting is per-field: a field that is not shaped like an
  `XX:T:VALUE` tag (a bare UMI, an Illumina CASAVA string, free text, or a value
  other than a `Z` string that holds a space) is skipped while the tags in the
  same comment still lift. A tag whose value is outside the SAM grammar of its
  type or a BAM record cannot hold, or a tag given twice in one comment, stops
  the run with an error that identifies the record and the tag.

  Scenario: a single SAM tag in a FASTQ comment becomes a real SAM tag
    Given a file "in.fq" containing:
      """
      @read1 RX:Z:ACGT
      AAAACCCC
      +
      IIIIIIII
      """
    When I run `unmux in.fq --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "RX:Z:ACGT"

  Scenario: a non-tag field is skipped while the valid tags in the same comment lift
    Given a file "in.fq" containing:
      """
      @read1	notatag	RX:Z:ACGT	BC:Z:GGGG
      AAAACCCC
      +
      IIIIIIII
      """
    When I run `unmux in.fq --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "RX:Z:ACGT"
    And the file "out.sam" contains "BC:Z:GGGG"
    And the file "out.sam" does not contain "notatag"

  Scenario: a bare UMI in the comment is left alone, not lifted into a tag
    Given a file "in.fq" containing:
      """
      @read1 GATTACAUMI
      AAAACCCC
      +
      IIIIIIII
      """
    When I run `unmux in.fq --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "AAAACCCC"
    And the file "out.sam" does not contain "GATTACAUMI"

  Scenario: a numeric or array value followed by a space is not a tag, so its field is skipped
    Given a file "in.fq" containing:
      """
      @r1 ns:i:5 BC:Z:ACGT
      AAAACCCC
      +
      IIIIIIII
      @r2 ch:i:12 extra
      AAAACCCC
      +
      IIIIIIII
      @r3 ML:B:C,1,2 MM:Z:C+m,0;
      AAAACCCC
      +
      IIIIIIII
      @r4	ns:i:5 	BC:Z:ACGT
      AAAACCCC
      +
      IIIIIIII
      """
    When I run `unmux in.fq --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      r1	4	*	0	0	*	*	0	0	AAAACCCC	IIIIIIII	RG:Z:in
      r2	4	*	0	0	*	*	0	0	AAAACCCC	IIIIIIII	RG:Z:in
      r3	4	*	0	0	*	*	0	0	AAAACCCC	IIIIIIII	RG:Z:in
      r4	4	*	0	0	*	*	0	0	AAAACCCC	IIIIIIII	BC:Z:ACGT	RG:Z:in
      """

  Scenario: a base-modification FASTQ keeps its MM, ML, and ns tags
    Given a file "mod.fq" containing:
      """
      @read1	MM:Z:C+m,0,1;	ML:B:C,200,100	ns:i:4000000000
      ACGTACGTCC
      +
      IIIIIIIIII
      """
    When I run `unmux mod.fq --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      read1	4	*	0	0	*	*	0	0	ACGTACGTCC	IIIIIIIIII	MM:Z:C+m,0,1;	ML:B:C,200,100	ns:i:4000000000	RG:Z:mod
      """

  Scenario Outline: the tags of <input> lose NM, MD and bwa's pa:f but keep Dorado's pa:B:i
    Given a file "in.fq" containing:
      """
      @r1	NM:i:10	MD:Z:16^T45	pa:f:0.95	RX:Z:ACGT
      ACGTACGT
      +
      IIIIIIII
      @r2	pa:B:i,1,2,3	AS:i:8	qs:i:10
      TTTTGGGG
      +
      JJJJJJJJ
      """
    And a file "in.fa" containing:
      """
      >r1	NM:i:10	MD:Z:16^T45	pa:f:0.95	RX:Z:ACGT
      ACGTACGT
      >r2	pa:B:i,1,2,3	AS:i:8	qs:i:10
      TTTTGGGG
      """
    And a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII	NM:i:10	MD:Z:16^T45	pa:f:0.95	RX:Z:ACGT
      r2	4	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ	pa:B:i,1,2,3	AS:i:8	qs:i:10
      """
    # bwa writes pa as a float, the score ratio to an ALT hit; Dorado writes its poly(A) positions as pa:B:i.
    When I run `unmux <input> --pool p --out whole.sam`
    Then the exit code is 0
    And the records of "whole.sam" are:
      """
      r1	4	*	0	0	*	*	0	0	ACGTACGT	<q1>	RX:Z:ACGT	RG:Z:p
      r2	4	*	0	0	*	*	0	0	TTTTGGGG	<q2>	pa:B:i,1,2,3	qs:i:10	RG:Z:p
      """
    When I run `unmux <input> --pool p --extract r=0:2:end --template r --out cut.sam`
    Then the exit code is 0
    And the records of "cut.sam" are:
      """
      r1	4	*	0	0	*	*	0	0	GTACGT	<c1>	RX:Z:ACGT	RG:Z:p
      r2	4	*	0	0	*	*	0	0	TTGGGG	<c2>	pa:B:i,1,2,3	qs:i:10	RG:Z:p
      """

    Examples:
      | input  | q1       | q2       | c1     | c2     |
      | in.fq  | IIIIIIII | JJJJJJJJ | IIIIII | JJJJJJ |
      | in.fa  | *        | *        | *      | *      |
      | in.sam | IIIIIIII | JJJJJJJJ | IIIIII | JJJJJJ |

  Scenario Outline: a FASTQ written from <input> input leaves alignment tags out of its comments
    Given a file "in.fq" containing:
      """
      @r1	NM:i:10	MD:Z:16^T45	pa:f:0.95	RX:Z:ACGT
      ACGTACGT
      +
      IIIIIIII
      @r2	pa:B:i,1,2,3	AS:i:8	qs:i:10
      TTTTGGGG
      +
      JJJJJJJJ
      """
    And a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	ACGTACGT	IIIIIIII	NM:i:10	MD:Z:16^T45	pa:f:0.95	RX:Z:ACGT
      r2	4	*	0	0	*	*	0	0	TTTTGGGG	JJJJJJJJ	pa:B:i,1,2,3	AS:i:8	qs:i:10
      """
    When I run `unmux <input> --out out.fq`
    Then the exit code is 0
    And the file "out.fq" is:
      """
      @r1 RX:Z:ACGT
      ACGTACGT
      +
      IIIIIIII
      @r2 pa:B:i,1,2,3	qs:i:10
      TTTTGGGG
      +
      JJJJJJJJ
      """

    Examples:
      | input  |
      | in.fq  |
      | in.sam |

  Scenario: every SAM tag type survives a SAM to FASTQ to SAM round trip
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      read1	4	*	0	0	*	*	0	0	ACGT	IIII	xa:A:Q	xn:i:-2147483648	xu:i:4294967295	xf:f:1.5	xz:Z:hello world	xh:H:CAFE	xc:B:c,-128,127	ML:B:C,200,100	xs:B:s,-32768,32767	xS:B:S,0,65535	xi:B:i,-2147483648,2147483647	xI:B:I,0,4294967295	xF:B:f,1.5,-2.25	xe:B:C
      """
    When I run `unmux in.sam --out mid.fq`
    Then the exit code is 0
    When I run `unmux mid.fq --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      read1	4	*	0	0	*	*	0	0	ACGT	IIII	xa:A:Q	xn:i:-2147483648	xu:i:4294967295	xf:f:1.5	xz:Z:hello world	xh:H:CAFE	xc:B:c,-128,127	ML:B:C,200,100	xs:B:s,-32768,32767	xS:B:S,0,65535	xi:B:i,-2147483648,2147483647	xI:B:I,0,4294967295	xF:B:f,1.5,-2.25	xe:B:C	RG:Z:mid
      """

  Scenario: a FASTA description lifts array and hex tags
    Given a file "in.fa" containing:
      """
      >read1 ML:B:C,200,100	xh:H:CAFE
      ACGTACGT
      """
    When I run `unmux in.fa --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      read1	4	*	0	0	*	*	0	0	ACGTACGT	*	ML:B:C,200,100	xh:H:CAFE	RG:Z:in
      """

  Scenario: a nameless record keeps its empty name and its tags through a FASTA round trip
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      *	4	*	0	0	*	*	0	0	ACGT	IIII	xa:A:Q	ML:B:C,1,2
      """
    When I run `unmux in.sam --out mid.fa`
    Then the exit code is 0
    When I run `unmux mid.fa --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      *	4	*	0	0	*	*	0	0	ACGT	*	xa:A:Q	ML:B:C,1,2	RG:Z:mid
      """

  Scenario: two nameless records read back from FASTA as single-end records, not a pair
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      *	4	*	0	0	*	*	0	0	ACGT	IIII	xa:A:Q
      *	4	*	0	0	*	*	0	0	CCGG	IIII	xa:A:R
      """
    When I run `unmux in.sam --out two.fa`
    Then the exit code is 0
    When I run `unmux two.fa --out back.sam`
    Then the exit code is 0
    And the records of "back.sam" are:
      """
      *	4	*	0	0	*	*	0	0	ACGT	*	xa:A:Q	RG:Z:two
      *	4	*	0	0	*	*	0	0	CCGG	*	xa:A:R	RG:Z:two
      """

  Scenario: three nameless records read back from FASTA as single-end records
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      *	4	*	0	0	*	*	0	0	ACGT	IIII	xa:A:Q
      *	4	*	0	0	*	*	0	0	CCGG	IIII	xa:A:R
      *	4	*	0	0	*	*	0	0	TTAA	IIII	xa:A:S
      """
    When I run `unmux in.sam --out three.fa`
    Then the exit code is 0
    When I run `unmux three.fa --out back.sam`
    Then the exit code is 0
    And the records of "back.sam" are:
      """
      *	4	*	0	0	*	*	0	0	ACGT	*	xa:A:Q	RG:Z:three
      *	4	*	0	0	*	*	0	0	CCGG	*	xa:A:R	RG:Z:three
      *	4	*	0	0	*	*	0	0	TTAA	*	xa:A:S	RG:Z:three
      """

  Scenario: a FASTA name after a space is kept unless it is shaped like a tag
    Given a file "in.fa" containing:
      """
      > read1 desc
      ACGT
      >	read2	BC:Z:GGGG
      ACGT
      > xa:A:Q	ML:B:C,1,2
      ACGT
      """
    When I run `unmux in.fa --out out.sam`
    Then the exit code is 0
    And the records of "out.sam" are:
      """
      read1	4	*	0	0	*	*	0	0	ACGT	*	RG:Z:in
      read2	4	*	0	0	*	*	0	0	ACGT	*	BC:Z:GGGG	RG:Z:in
      *	4	*	0	0	*	*	0	0	ACGT	*	xa:A:Q	ML:B:C,1,2	RG:Z:in
      """

  Scenario Outline: a comment tag that a BAM record cannot hold stops the run: <field>
    Given a file "in.fq" containing:
      """
      @read1	RX:Z:ACGT	<field>
      AAAACCCC
      +
      IIIIIIII
      """
    When I run `unmux in.fq --out out.sam`
    Then the exit code is 1
    And stderr contains "malformed FASTQ record `read1`: the comment tag `<tag>`"
    And the directory "." holds exactly "in.fq"

    Examples:
      | field            | tag  |
      | ns:i:4294967296  | ns:i |
      | ns:i:-2147483649 | ns:i |
      | ML:B:C,256       | ML:B |
      | ML:B:C,1,x       | ML:B |
      | ML:B:Q,1         | ML:B |
      | xh:H:CAF         | xh:H |
      | CO:Z:café        | CO:Z |
      | xf:f:nan         | xf:f |
      | xf:f:1e40        | xf:f |
      | RX:Z:GGGG        | RX   |

  Scenario: the help says that comment tags are lifted and a malformed one stops the run
    When I run `unmux --help`
    Then the exit code is 0
    And stdout contains "SAM tags in a FASTX record-name comment are lifted onto the record,"
    And stdout contains "and a malformed one stops the run."
