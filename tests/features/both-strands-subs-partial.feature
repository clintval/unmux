Feature: Reverse-complement, substitution, and partial end matching
  Tag group attributes for matching on the reverse complement strand, error-correcting
  matches to canonical sequences, and truncated overhang matching. Each scenario runs the
  real binary on a small FASTQ fixture and asserts on the extracted span surfaced as a tag
  and template body in the FASTQ read-name comment.

  Scenario: Reverse-complement matching finds a tag on the antisense strand
    Given a file "read.fq" containing:
      """
      @r1
      GGGGGGGGGGCACGTTAAAAAA
      +
      IIIIIIIIIIIIIIIIIIIIII
      """
    # both_strands=true lets AACGTG match its rc (CACGTT) at [10,16). The corrected tag is the declared
    # canonical AACGTG (both_strands is matching-only); use `~` on the tag/template for the read strand.
    When I run `unmux read.fq --group grp={AACGTG} --group grp::loc=0:10:16,both_strands=true --extract bc=@grp --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:AACGTG"
    And stdout contains "AACGTG"

  Scenario: Without both_strands the antisense instance does not match
    Given a file "read.fq" containing:
      """
      @r1
      GGGGGGGGGGCACGTTAAAAAA
      +
      IIIIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux read.fq --group grp={AACGTG} --group grp::loc=0:10:16 --extract bc=@grp --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout does not contain "CACGTT"

  Scenario: A matched tag emits corrected bases to BC and raw bases to CR from one match
    Given a file "read.fq" containing:
      """
      @r1
      GAAGGGTTTTTT
      +
      IIIIIIIIIIII
      """
    # correct defaults to true, so BC carries the error-corrected canonical GAAGAG; CR::raw=true emits
    # the observed GAAGGG from the same match, and the templated read body also keeps the observed bases.
    When I run `unmux read.fq --group grp={GAAGAG,CCCCCC,TTACGT} --group grp::loc=0:0:6,dist=1 --extract bc=@grp --tag BC=bc --tag BC::qual=none --tag CR=bc --tag CR::raw=true --tag CR::qual=none --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:GAAGAG"
    And stdout contains "CR:Z:GAAGGG"
    And stdout contains "GAAGGG"

  Scenario: raw=true preserves the observed sequence
    Given a file "read.fq" containing:
      """
      @r1
      GGGTCCTTTTTT
      +
      IIIIIIIIIIII
      """
    When I run `unmux read.fq --group grp={GGATCC,AATTGG} --group grp::loc=0:0:6,dist=1 --extract bc=@grp --tag BC=bc --tag BC::raw=true --template bc::raw=true`
    Then the exit code is 0
    And stdout contains "BC:Z:GGGTCC"
    And stdout contains "GGGTCC"
    And stdout does not contain "GGATCC"

  Scenario: partial5 matches a tag truncated at the read 5-prime end
    Given a file "read.fq" containing:
      """
      @r1
      CGATCGTTTTTTTTTT
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux read.fq --group grp={ATCGATCG} --group grp::loc=0:0:8,partial5=4:0 --extract bc=@grp --tag BC=bc --tag BC::raw=true --template bc::raw=true`
    Then the exit code is 0
    And stdout contains "BC:Z:CGATCG"
    And stdout contains "CGATCG"

  Scenario: Without partial5 a truncated tag does not match
    Given a file "read.fq" containing:
      """
      @r1
      CGATCGTTTTTTTTTT
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux read.fq --group grp={ATCGATCG} --group grp::loc=0:0:8 --extract bc=@grp --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout does not contain "BC:Z:"

  Scenario: partial3 matches a tag truncated at the read 3-prime end
    Given a file "read.fq" containing:
      """
      @r1
      TTTTTTTTTTTACGTA
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux read.fq --group grp={ACGTACGT} --group grp::loc=0:11:end,partial3=4:0 --extract bc=@grp --tag BC=bc --tag BC::raw=true --template bc::raw=true`
    Then the exit code is 0
    And stdout contains "BC:Z:ACGTA"
    And stdout contains "ACGTA"

  Scenario: Reverse-complement error-correction emits the declared canonical, not its read-strand rc
    Given a file "read.fq" containing:
      """
      @r1
      GGGGGGGGGGCACGATAAAAAA
      +
      IIIIIIIIIIIIIIIIIIIIII
      """
    # The tag AACGTG matches antisense (its rc CACGTT) at [10:16]=CACGAT (1 mismatch). The corrected tag
    # is the declared canonical AACGTG (both_strands is matching-only), NOT its read-strand rc CACGTT; `~` on
    # the tag/template is how you would ask for the read strand.
    When I run `unmux read.fq --group grp={AACGTG,CCGGAA} --group grp::loc=0:10:16,both_strands=true,dist=1 --extract bc=@grp --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:AACGTG"
    And stdout contains "AACGTG"

  Scenario: A group loaded from a tag file matches on the antisense strand under both_strands=true
    Given a file "tags.tsv" containing:
      """
      id	seq
      t1	AACGTG
      """
    And a file "read.fq" containing:
      """
      @r1
      GGGGGGGGGGCACGTTAAAAAA
      +
      IIIIIIIIIIIIIIIIIIIIII
      """
    # The tag set comes from a headered TSV (id/seq columns) rather than an inline {...}; both_strands=true is
    # what lets the file-sourced AACGTG match its reverse complement CACGTT at [10,16). The corrected tag
    # is the declared canonical AACGTG. The companion "Without both_strands" scenario above (inline tag) drops
    # the same antisense instance, so both_strands=true remains load-bearing when the tags arrive from a file.
    When I run `unmux read.fq --group grp=tags.tsv --group grp::loc=0:10:16,both_strands=true --extract bc=@grp --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:AACGTG"

  Scenario: A multi-row tag file matches forward and antisense tags in one both_strands=true run
    Given a file "tags.tsv" containing:
      """
      id	seq
      t1	AACGTG
      t2	GGATCC
      """
    And a file "read.fq" containing:
      """
      @fwd
      AAAAAAAAAAGGATCCTTTTTT
      +
      IIIIIIIIIIIIIIIIIIIIII
      @rev
      GGGGGGGGGGCACGTTAAAAAA
      +
      IIIIIIIIIIIIIIIIIIIIII
      """
    # Two tags from the file: @fwd carries t2 (GGATCC) on the forward strand at [10,16), @rev carries the
    # reverse complement of t1 (AACGTG -> CACGTT) at [10,16). With both_strands=true both resolve in a single
    # run with no file-sourcing error, so the observed bases are GGATCC and CACGTT respectively.
    When I run `unmux read.fq --group grp=tags.tsv --group grp::loc=0:10:16,both_strands=true --extract bc=@grp --tag BC=bc --tag BC::raw=true --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:GGATCC"
    And stdout contains "BC:Z:CACGTT"

  Scenario: a ~ tag stream puts a reverse-complement-matched barcode on the read strand
    Given a file "read.fq" containing:
      """
      @r1
      GGGGGGGGGGCACGATAAAAAA
      +
      IIIIIIIIIIIIIIIIIIIIII
      """
    # grp's AACGTG matches antisense (its rc CACGTT) at [10:16]=CACGAT (1 mismatch). The corrected tag
    # is the declared canonical AACGTG by default; a `~` on the tag stream (BC=~bc) reverse-complements
    # it onto the read strand, so BC is CACGTT here, the same orientation as the read.
    When I run `unmux read.fq --group grp={AACGTG,CCGGAA} --group grp::loc=0:10:16,both_strands=true,dist=1 --extract bc=@grp --tag BC=~bc --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:CACGTT"

  Scenario: the same barcode forward and antisense yields one identical declared BC (identity-stable)
    Given a file "read.fq" containing:
      """
      @fwd
      AACGTGAAAA
      +
      IIIIIIIIII
      @rev
      CACGTTAAAA
      +
      IIIIIIIIII
      """
    # @fwd carries AACGTG forward; @rev carries its reverse complement CACGTT, found via both_strands=true.
    # Both are the SAME barcode, so the corrected BC is the declared canonical AACGTG for both reads,
    # not strand-dependent. The old behavior emitted CACGTT for the antisense read, splitting one
    # barcode into two values; this guards against that footgun.
    When I run `unmux read.fq --group grp={AACGTG} --group grp::loc=0:0:6,both_strands=true --extract bc=@grp --tag BC=bc --template bc`
    Then the exit code is 0
    And stdout contains "BC:Z:AACGTG"
    And stdout does not contain "BC:Z:CACGTT"
