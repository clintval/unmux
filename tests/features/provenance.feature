Feature: BAM/CRAM provenance headers and demux statistics
  Every output SAM/BAM/CRAM carries exactly one provenance program record with the full
  command line, and per-sample fan-out files carry their read group; a pass-through `--out` and the
  raw `--unassigned`/`--remove` bins carry a default `@RG` whose ID/SM/LB are the pool id, so every
  record belongs to a read group. A file that a bin shares with samples declares the read groups of
  both. Demux statistics are not
  embedded in output headers (that would force a full per-file rewrite once counts are tallied);
  they go to stderr and the --metrics-per-sample / --metrics-summary TSVs. Header internals are read
  back with `noodles` (the `the BAM header of "..."` steps decode the BGZF-framed header the same way
  the binary wrote it), so the @PG provenance record and per-sample @RG read groups are asserted
  directly rather than parked.

  Scenario: every output BAM carries a single @PG provenance record with unmux's id, version, and command line
    Given a file "reads.fq" containing:
      """
      @r1
      AAGGTTCCACGT
      +
      IIIIIIIIIIII
      @r2
      TTCCAAGGACGT
      +
      IIIIIIIIIIII
      """
    # Each fan-out BAM gets exactly one @PG program record stamped with ID/PN unmux, the running version,
    # and the verbatim command line. Asserting on the recorded flags (--pool lib01 and the exact --sample
    # token) proves the CL is the real invocation, not a placeholder.
    When I run `unmux reads.fq --pool lib01 --group cbt={c1=AAGGTTCC,c2=TTCCAAGG} --group cbt::loc=0:0:8,dist=0 --extract body=0:0:end --template body --sample dna01::%pool=cbt::c1 --sample dna02::%pool=cbt::c2 --out output/%sample.%sub_sample.bam`
    Then the exit code is 0
    And a file "output/dna01.lib01.bam" exists
    And the BAM header of "output/dna01.lib01.bam" has exactly 1 @PG line
    And the BAM header of "output/dna01.lib01.bam" contains "ID:unmux"
    And the BAM header of "output/dna01.lib01.bam" contains "PN:unmux"
    And the BAM header of "output/dna01.lib01.bam" records the running unmux version
    And the BAM header of "output/dna01.lib01.bam" contains "--pool lib01"
    And the BAM header of "output/dna01.lib01.bam" contains "--sample dna01::%pool=cbt::c1"

  Scenario: CRAM output is produced for each fan-out sample
    Given a file "reads.fq" containing:
      """
      @a
      AAAACCCCGGGGTTTT
      +
      IIIIIIIIIIIIIIII
      @b
      CCCCGGGGTTTTAAAA
      +
      IIIIIIIIIIIIIIII
      """
    # Files are eagerly created, so existence alone does not prove routing; the per-sample metric
    # (each sample = half the pool) proves the two reads actually fanned out into the two CRAMs.
    When I run `unmux reads.fq --group bc={s1=AAAA,s2=CCCC} --sample dna01=bc::s1 --sample dna02=bc::s2 --out output/%sample.cram --metrics-per-sample m.tsv`
    Then the exit code is 0
    And a file "output/dna01.cram" exists
    And a file "output/dna02.cram" exists
    And the file "m.tsv" contains "0.500000"

  Scenario: metrics TSVs are the primary, joinable output for demux stats
    Given a file "reads.fq" containing:
      """
      @a
      AAAACCCCGGGGTTTT
      +
      IIIIIIIIIIIIIIII
      @b
      CCCCGGGGTTTTAAAA
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --pool lib01 --group bc={s1=AAAA,s2=CCCC} --sample dna01=bc::s1 --sample dna02=bc::s2 --out output/%sample.bam --metrics-per-sample lib01.per_sample.tsv --metrics-summary lib01.summary.tsv`
    Then the exit code is 0
    And a file "lib01.per_sample.tsv" exists
    And a file "lib01.summary.tsv" exists
    And the file "lib01.per_sample.tsv" contains "pool"
    And the file "lib01.per_sample.tsv" contains "sample"
    And the file "lib01.per_sample.tsv" contains "sub_sample"
    And the file "lib01.per_sample.tsv" contains "frac_bases_unextracted"
    And the file "lib01.per_sample.tsv" contains "dna01"
    And the file "lib01.per_sample.tsv" contains "dna02"
    # Routing actually happened: each sample took half the pool and all reads were assigned (a
    # zero-read declared sample would show 0.000000 here and frac_assigned below would not be 1).
    And the file "lib01.per_sample.tsv" contains "0.500000"
    And the file "lib01.summary.tsv" contains "pool"
    And the file "lib01.summary.tsv" contains "total_reads"
    And the file "lib01.summary.tsv" contains "lib01"
    And the file "lib01.summary.tsv" contains "1.000000"

  Scenario: per-sample fan-out writes one BGZF BAM per sample
    Given a file "reads.fq" containing:
      """
      @a
      AAAACCCCGGGGTTTT
      +
      IIIIIIIIIIIIIIII
      @b
      CCCCGGGGTTTTAAAA
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --group bc={s1=AAAA,s2=CCCC} --sample dna01=bc::s1 --sample dna02=bc::s2 --out output/%sample.raw.unmapped.bam --metrics-per-sample m.tsv`
    Then the exit code is 0
    And a file "output/dna01.raw.unmapped.bam" exists
    And a file "output/dna02.raw.unmapped.bam" exists
    And the file "output/dna01.raw.unmapped.bam" is gzip-compressed
    And the file "output/dna02.raw.unmapped.bam" is gzip-compressed
    And the file "m.tsv" contains "0.500000"

  Scenario: a pass-through SAM carries a default pool read group on the header and every record
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGT
      +
      IIIIIIII
      """
    # No --sample, so the pool passes through undemultiplexed. The output still carries a read group
    # (ID/SM/LB = the pool id) so the uBAM is valid downstream, and the record references it via RG:Z.
    When I run `unmux reads.fq --pool lib01 --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "@RG"
    And the file "out.sam" contains "ID:lib01"
    And the file "out.sam" contains "SM:lib01"
    And the file "out.sam" contains "LB:lib01"
    And the file "out.sam" contains "RG:Z:lib01"

  Scenario: the --unassigned bin carries the default pool read group
    Given a file "reads.fq" containing:
      """
      @r1
      TTTTTTTT
      +
      IIIIIIII
      """
    # The read matches no sample, so it lands in --unassigned. That bin is not a per-sample target, so
    # it gets the default pool @RG (ID/SM/LB = pool) and the record references it via RG:Z.
    When I run `unmux reads.fq --pool lib01 --group bc={s1=AAAAAAAA} --group bc::loc=0:0:8,dist=0,minFindsPerGroup=1 --sample dna01=bc::s1 --unassigned un.sam --out assigned.sam`
    Then the exit code is 0
    And the file "un.sam" contains "@RG"
    And the file "un.sam" contains "ID:lib01"
    And the file "un.sam" contains "RG:Z:lib01"

  Scenario Outline: a <format> file that a bin shares with a sample declares the bin's read group
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      @r2
      TTTTTTTTGGCCATGC
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --pool p --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A --out out.<format> <bin>`
    Then the exit code is 0
    And the <what> "out.<format>" contains "@RG	ID:s1	SM:s1"
    And the <what> "out.<format>" contains "@RG	ID:p	SM:p	LB:p"
    And the records of "out.<format>" are:
      """
      r1	4	*	0	255	*	*	0	0	AAAAAAAACCCCGGGG	IIIIIIIIIIIIIIII	RG:Z:s1
      r2	4	*	0	255	*	*	0	0	TTTTTTTTGGCCATGC	IIIIIIIIIIIIIIII	RG:Z:p
      """

    Examples:
      | format | what          | bin                    |
      | sam    | file          | --unassigned out.sam   |
      | bam    | BAM header of | --remove bc::T=out.bam |

  Scenario: a CRAM file that a bin shares with a sample declares the bin's read group
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      @r2
      TTTTTTTTGGCCATGC
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --pool p --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A --out out.cram --unassigned out.cram`
    Then the exit code is 0
    # Splitting the CRAM by the read groups its header declares finds the bin's record under the pool's.
    When I run `unmux out.cram --group rg=@RG --sample-from-group rg --out split/%sample.sam`
    Then the exit code is 0
    And the records of "split/p.sam" are:
      """
      r2	4	*	0	255	*	*	0	0	TTTTTTTTGGCCATGC	IIIIIIIIIIIIIIII	RG:Z:p
      """

  Scenario Outline: SAM on stdout that a bin shares with a sample declares the bin's read group (<outputs>)
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AAAAAAAACCCCGGGG	IIIIIIIIIIIIIIII
      r2	4	*	0	0	*	*	0	0	TTTTTTTTGGCCATGC	IIIIIIIIIIIIIIII
      """
    When I run `unmux in.sam --pool p --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A <outputs>`
    Then the exit code is 0
    And stdout contains "@RG	ID:s1	SM:s1"
    And stdout contains "@RG	ID:p	SM:p	LB:p"
    And stdout contains "TTTTTTTTGGCCATGC	IIIIIIIIIIIIIIII	RG:Z:p"

    Examples:
      | outputs                          |
      | --out - --unassigned -           |
      | --out - --unassigned /dev/stdout |
      | --unassigned -                   |

  Scenario: SAM on the default stdout declares the samples' read groups
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AAAAAAAACCCCGGGG	IIIIIIIIIIIIIIII
      r2	4	*	0	0	*	*	0	0	TTTTTTTTGGCCATGC	IIIIIIIIIIIIIIII
      """
    When I run `unmux in.sam --pool p --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A`
    Then the exit code is 0
    And stdout contains "@RG	ID:s1	SM:s1"
    And stdout does not contain "@RG	ID:p"
    And stdout contains "AAAAAAAACCCCGGGG	IIIIIIIIIIIIIIII	RG:Z:s1"

  Scenario: a SAM file that a bin shares is refused when the pool read group is also a sample's
    Given a file "reads.fq" containing:
      """
      @r1
      AAAAAAAACCCCGGGG
      +
      IIIIIIIIIIIIIIII
      @r2
      TTTTTTTTGGCCATGC
      +
      IIIIIIIIIIIIIIII
      """
    When I run `unmux reads.fq --pool s1 --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A --out out.sam --unassigned out.sam`
    Then the exit code is 1
    And stderr contains "--unassigned records would share read group "s1" with a sample in "out.sam"; set --pool to another name"
    And a file "out.sam" does not exist
    # FASTQ has no read groups, so the same samples and bin may share a FASTQ file.
    When I run `unmux reads.fq --pool s1 --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A --out out.fq --unassigned out.fq`
    Then the exit code is 0
    And the file "out.fq" contains "@r2"

  Scenario: the pool read group is checked on stdout that a bin shares, but not on /dev/null
    Given a file "in.sam" containing:
      """
      @HD	VN:1.6
      r1	4	*	0	0	*	*	0	0	AAAAAAAACCCCGGGG	IIIIIIIIIIIIIIII
      r2	4	*	0	0	*	*	0	0	TTTTTTTTGGCCATGC	IIIIIIIIIIIIIIII
      """
    When I run `unmux in.sam --pool s1 --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A --out - --unassigned -`
    Then the exit code is 1
    And stderr contains "--unassigned records would share read group "s1" with a sample on stdout; set --pool to another name"
    # Without --out, records go to stdout as with --out -.
    When I run `unmux in.sam --pool s1 --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A --unassigned -`
    Then the exit code is 1
    And stderr contains "--unassigned records would share read group "s1" with a sample on stdout; set --pool to another name"
    # /dev/null keeps nothing, so its read groups are not checked.
    When I run `unmux in.sam --pool s1 --group bc={A=AAAAAAAA,T=TTTTTTTT} --sample s1=bc::A --out /dev/null --unassigned /dev/null`
    Then the exit code is 0

  Scenario: a pass-through FASTQ carries no read group (FASTX has no @RG concept)
    Given a file "reads.fq" containing:
      """
      @r1
      ACGTACGT
      +
      IIIIIIII
      """
    When I run `unmux reads.fq --pool lib01 --out out.fq`
    Then the exit code is 0
    And the file "out.fq" does not contain "RG:Z:"

  Scenario Outline: The default pool is the input file name without its format extension
    Given a file "<input>" containing:
      """
      @r
      ACGT
      +
      IIII
      """
    When I run `unmux <input> --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "@RG	ID:<pool>	SM:<pool>	LB:<pool>"
    And the file "out.sam" contains "RG:Z:<pool>"
    When I run `unmux <input> --out %pool.out.fq`
    Then the exit code is 0
    And a file "<pool>.out.fq" exists

    Examples:
      | input        | pool       |
      | lib1.fq      | lib1       |
      | lib1.fastq   | lib1       |
      | LIB1.FQ      | LIB1       |
      | a.b.fq       | a.b        |
      | reads.txt    | reads.txt  |
      | .fq          | .fq        |

  Scenario: Two inputs share a default pool without an extension and lone separators
    Given a file "lib1_R1.fq" containing:
      """
      @r
      ACGT
      +
      IIII
      """
    And a file "lib1_R2.fq" containing:
      """
      @r
      TTTT
      +
      JJJJ
      """
    When I run `unmux lib1_R1.fq lib1_R2.fq --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "@RG	ID:lib1_R	SM:lib1_R	LB:lib1_R"

  Scenario: Input names that differ inside one character do not stop the run
    Given a file "é.fq" containing:
      """
      @r
      ACGT
      +
      IIII
      """
    And a file "è.fq" containing:
      """
      @r
      TTTT
      +
      JJJJ
      """
    When I run `unmux é.fq è.fq --out out.sam`
    Then the exit code is 0
    And the file "out.sam" contains "@RG	ID:pool	SM:pool	LB:pool"

  Scenario: A name that a SAM header cannot hold is replaced in the read group only
    Given a file "ä.fq" containing:
      """
      @r
      ACGT
      +
      IIII
      """
    When I run `unmux ä.fq --out out.sam`
    Then the exit code is 0
    And stderr contains "--pool"
    And the file "out.sam" contains "@RG	ID:pool	SM:pool	LB:pool"
    And the file "out.sam" contains "RG:Z:pool"
    When I run `unmux ä.fq --out %pool.out.fq`
    Then the exit code is 0
    And a file "ä.out.fq" exists

  Scenario Outline: A --pool that a SAM header cannot hold is refused before any output for <format>
    Given a file "reads.fq" containing:
      """
      @r
      ACGT
      +
      IIII
      """
    When I run `unmux reads.fq --pool données --out out.<format>`
    Then the exit code is 1
    And stderr contains "--pool"
    And stderr contains "printable ASCII"
    And a file "out.<format>" does not exist
    When I run `unmux reads.fq --pool données --out %pool.fq`
    Then the exit code is 0
    And a file "données.fq" exists

    Examples:
      | format |
      | sam    |
      | bam    |
      | cram   |

  Scenario Outline: The command line is recorded in a header whatever its characters (<format>)
    Given a file "ä.fq" containing:
      """
      @r
      ACGT
      +
      IIII
      """
    When I run `unmux ä.fq --pool p --out out.<format>`
    Then the exit code is 0
    And the <what> "out.<format>" contains "\u{e4}.fq"
    And the <what> "out.<format>" contains "PN:unmux"

    Examples:
      | format | what          |
      | sam    | file          |
      | bam    | BAM header of |
