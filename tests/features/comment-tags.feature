Feature: Lifting SAM tags from FASTX read-name comments
  A FASTX read-name comment may carry SAM auxiliary tags (the samtools `fastq -T`
  convention that unmux itself writes on FASTX output). On input those tags are
  lifted back into real SAM tags, the inverse of that output, so a tag-bearing
  FASTQ converts to a uBAM/SAM with real tags and no extraction configured.
  Every type unmux writes reads back: `A`, `f`, `Z`, `H`, `B` arrays of every
  subtype, and `i` from -2147483648 to 4294967295 in the narrowest BAM type that
  holds it. Lifting is per-field: a field that is not shaped like an
  `XX:T:VALUE` tag (a bare UMI, an Illumina CASAVA string, free text) is skipped
  while the tags in the same comment still lift, but a tag whose value its type
  does not allow, or a BAM record cannot hold, stops the run with an error that
  identifies the record and the tag.

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
      r1	4	*	0	255	*	*	0	0	AAAACCCC	IIIIIIII	RG:Z:in
      r2	4	*	0	255	*	*	0	0	AAAACCCC	IIIIIIII	RG:Z:in
      r3	4	*	0	255	*	*	0	0	AAAACCCC	IIIIIIII	RG:Z:in
      r4	4	*	0	255	*	*	0	0	AAAACCCC	IIIIIIII	BC:Z:ACGT	RG:Z:in
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
      read1	4	*	0	255	*	*	0	0	ACGTACGTCC	IIIIIIIIII	MM:Z:C+m,0,1;	ML:B:C,200,100	ns:i:4000000000	RG:Z:mod
      """

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
      read1	4	*	0	255	*	*	0	0	ACGT	IIII	xa:A:Q	xn:i:-2147483648	xu:i:4294967295	xf:f:1.5	xz:Z:hello world	xh:H:CAFE	xc:B:c,-128,127	ML:B:C,200,100	xs:B:s,-32768,32767	xS:B:S,0,65535	xi:B:i,-2147483648,2147483647	xI:B:I,0,4294967295	xF:B:f,1.5,-2.25	xe:B:C	RG:Z:mid
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
      read1	4	*	0	255	*	*	0	0	ACGTACGT	*	ML:B:C,200,100	xh:H:CAFE	RG:Z:in
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
