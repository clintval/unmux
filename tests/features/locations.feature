Feature: Location windows bound barcode search and improve specificity for short tags
  Location windows restrict the search region for tag matching, preventing false positives in short
  tags (6-8 bp) that might match elsewhere. Supports file-scoped windows, start-anchored syntax, and
  negative coordinates from the read end. Each scenario is designed so the window is load-bearing:
  removing or changing it flips the outcome (assigned vs unassigned, or which tag wins), rather than
  merely re-finding a tag a whole-read search would find anyway. (The unassigned bin holds raw input
  segments and never carries SAM tags, so non-assignment is shown by the raw read landing in un.fq.)

  Scenario: A window that excludes the tag leaves the read unassigned
    Given a file "reads.fq" containing:
      """
      @r1
      GGGGGGGGGGATACTGCCCCCCCCCCCCCCATACTGTTTT
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group bc={s1=ATACTG} --group bc::loc=0:0:10 --extract rest=0:16:end --template rest --sample S1=bc::s1 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "un.fq" contains "GGGGGGGGGGATACTGCCCCCCCCCCCCCCATACTGTTTT"

  Scenario: A file-1 window ignores a tag that lives only in file 0
    Given a file "r1.fq" containing:
      """
      @r1
      GGAATTCCACGTACGT
      +
      IIIIIIIIIIIIIIII
      """
    And a file "i1.fq" containing:
      """
      @r1
      TTTTTTTT
      +
      IIIIIIII
      """
    # GGAATTCC is only in file 0; the loc=1 window searches file 1, so the read is unassigned and its
    # raw segment lands in un.fq. A whole-read (all-files) search would assign it to S1 instead.
    When I run `unmux r1.fq i1.fq --group bc={s1=GGAATTCC} --group bc::loc=1:0 --extract cb=@bc --extract body=0:0:end --template body --tag CB=cb --sample S1=bc::s1 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "un.fq" contains "GGAATTCC"

  Scenario: A tail window (negative coordinates) ignores a tag near the read start
    Given a file "reads.fq" containing:
      """
      @r1
      CCCGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group bc={s1=CCCGGG} --group bc::loc=0:-10:end --extract body=0:0:end --template body --sample S1=bc::s1 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "un.fq" contains "CCCGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

  Scenario: Moving the window to the second tag's region selects the other tag
    Given a file "reads.fq" containing:
      """
      @r1
      AAATTTNNNNNNNNNNNNNNNNNNNNNNNNGGGCCCNNNN
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # Without a window the leftmost tag (tagA at position 0) wins; windowing [30,36) makes tagB win.
    When I run `unmux reads.fq --group bc={tagA=AAATTT,tagB=GGGCCC} --group bc::loc=0:30:36 --extract cb=@bc --extract body=0:0:30 --template body --tag CB=cb --sample SA=bc::tagA --sample SB=bc::tagB --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And a file "SB.fq" exists
    And the file "SB.fq" contains "CB:Z:GGGCCC"
    And the file "SB.fq" does not contain "CB:Z:AAATTT"

  Scenario: A file-0 window misses a tag that lives only in another file
    Given a file "r1.fq" containing:
      """
      @r1
      ACGTACGTACGTACGT
      +
      IIIIIIIIIIIIIIII
      """
    And a file "i1.fq" containing:
      """
      @r1
      TTTTGGAATTCCTTTT
      +
      IIIIIIIIIIIIIIII
      """
    # GGAATTCC is only in file 1; a file-0 window cannot see it, so the read is unassigned (the
    # default whole-read search across all files would find it). File-scoping is load-bearing.
    When I run `unmux r1.fq i1.fq --group bc={s1=GGAATTCC} --group bc::loc=0:0:end --extract body=0:0:end --template body --sample S1=bc::s1 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "un.fq" contains "TTTTGGAATTCCTTTT"

  Scenario: A negative file index is rejected (there is no any-file wildcard; omit loc instead)
    Given a file "r1.fq" containing:
      """
      @r1
      GGAATTCCACGTACGT
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux r1.fq --group bc={s1=GGAATTCC} --group bc::loc=-1:0:end --extract body=0:0:end --template body`
    Then the exit code is 1
    And stderr contains "loc file index"
    And stderr contains "non-negative integer"
