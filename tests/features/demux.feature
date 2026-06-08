Feature: Demux parsing and sample fan-out
  Sample demultiplexing routes each read to a fan-out target and writes one output per sample. These
  scenarios exercise the demux-specific machinery not covered by routing.feature (the --remove /
  --unassigned / AND-selector primitives) or end-to-end.feature / provenance.feature (per-sample
  FASTQ/BAM/CRAM fan-out and metrics): a 1:1 --sample-from-group with a concatenated
  match= target, and the %sub_sample placeholder that names the worked-example BAMs. Each
  scenario is built so the named feature decides the outcome; a decoy or a flipped parameter
  moves the read to a different file. (@PG / @RG header internals are asserted directly in
  provenance.feature / routing.feature via the `noodles`-backed header steps; these demux scenarios
  prove routing the simpler way - the file written in the right container, named by the right
  placeholders, with routing shown by per-sample FASTQ contents or the per-sample metrics.)

  Scenario: dual-index demux routes by the i7+i5 concatenation via sample-from-group
    Given a file "r1.fq" containing:
      """
      @q1
      TTAAGGCCTTAAGGCC
      +
      IIIIIIIIIIIIIIII
      """
    And a file "i1.fq" containing:
      """
      @q1
      AACCGGTT
      +
      IIIIIIII
      """
    And a file "i2.fq" containing:
      """
      @q1
      ACGTACGT
      +
      IIIIIIII
      """
    And a file "metadata.tsv" containing:
      """
      sample_id	barcode
      s1	AACCGGTTACGTACGT
      s2	ACGTACGTAACCGGTT
      """
    # The i7 read is AACCGGTT and the i5 read is ACGTACGT; match=i7+i5 joins them in that order to
    # AACCGGTTACGTACGT, which is s1's barcode. The decoy s2 carries the reversed concatenation
    # ACGTACGTAACCGGTT, so the join ORDER decides the sample: with match=i7+i5 the R1 body lands in
    # s1's file; flipping the group to match=i5+i7 would route the same body to s2 instead. The R1 body
    # TTAAGGCCTTAAGGCC is distinct from either index, so only the concatenated match can assign it.
    When I run `unmux --in 0=r1.fq --in 1=i1.fq --in 2=i2.fq --extract i7=1:0:8 --extract i5=2:0:8 --extract t.r1=0:0:end --group sample_bc=metadata.tsv --group sample_bc::match=i7+i5 --template t.r1 --sample-from-group sample_bc --out out/%sample.R%ordinal.fq`
    Then the exit code is 0
    And the file "out/s1.R1.fq" contains "TTAAGGCCTTAAGGCC"
    And the file "out/s2.R1.fq" does not contain "TTAAGGCCTTAAGGCC"

  Scenario: Per-index error budgets via two index groups and an AND selector (the BCL Convert model)
    Given a file "r1.fq" containing:
      """
      @r1
      GGGGGGGGGGGG
      +
      IIIIIIIIIIII
      @r2
      TTTTTTTTTTTT
      +
      IIIIIIIIIIII
      """
    And a file "i1.fq" containing:
      """
      @r1
      AAAAAAAT
      +
      IIIIIIII
      @r2
      AAAAAATT
      +
      IIIIIIII
      """
    And a file "i2.fq" containing:
      """
      @r1
      CCCCCCCG
      +
      IIIIIIII
      @r2
      CCCCCCCC
      +
      IIIIIIII
      """
    # Unlike match=i7+i5 (one mismatch budget over the 16 bp concatenation), each index is its own
    # group with its own dist=1. r1 has ONE mismatch in EACH index (i7 AAAAAAAT, i5 CCCCCCCG) and
    # assigns; r2 has TWO mismatches in i7 (AAAAAATT) and is rejected by g7's per-index cap, so it lands
    # unassigned. A single combined budget could not separate these (1+1 and 2+0 both total 2).
    When I run `unmux r1.fq i1.fq i2.fq --group g7={a=AAAAAAAA} --group g7::loc=1:0:8,dist=1 --group g5={b=CCCCCCCC} --group g5::loc=2:0:8,dist=1 --extract t=0:0:end --template t --sample S1=g7::a+g5::b --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "S1.fq" contains "GGGGGGGGGGGG"
    And the file "un.fq" contains "TTTTTTTTTTTT"
    And the file "S1.fq" does not contain "TTTTTTTTTTTT"

  Scenario: a pool fans out into per-sample BAMs named by %sample.%sub_sample
    Given a file "reads.fq" containing:
      """
      @r1
      AAGGTTCCACGTACGTACGTAAAA
      +
      IIIIIIIIIIIIIIIIIIIIIIII
      @r2
      TTCCAAGGTGCATGCATGCAGGGG
      +
      IIIIIIIIIIIIIIIIIIIIIIII
      """
    # The two reads carry distinct barcodes c1 (AAGGTTCC) and c2 (TTCCAAGG) routed to dna01 and dna02.
    # Each sample takes %pool (lib01) as its sub_sample, so %sample.%sub_sample names the files
    # dna01.lib01.raw.unmapped.bam and dna02.lib01.raw.unmapped.bam (the worked example).
    # Dropping ::%pool (a bare --sample dna01=...) expands %sub_sample to empty, naming the files
    # dna01..raw.unmapped.bam instead, so the .lib01. infix is load-bearing on the %sub_sample binding.
    # The .bam files are BGZF binary (gzip magic, and the packed bases never appear as text); routing is
    # proven by the per-sample metric showing each sample took half the pool.
    When I run `unmux reads.fq --pool lib01 --group cbt={c1=AAGGTTCC,c2=TTCCAAGG} --group cbt::loc=0:0:8,dist=0 --extract body=0:0:end --template body --sample dna01::%pool=cbt::c1 --sample dna02::%pool=cbt::c2 --out output/%sample.%sub_sample.raw.unmapped.bam --metrics-per-sample m.tsv`
    Then the exit code is 0
    And a file "output/dna01.lib01.raw.unmapped.bam" exists
    And a file "output/dna02.lib01.raw.unmapped.bam" exists
    And the file "output/dna01.lib01.raw.unmapped.bam" is gzip-compressed
    And the file "output/dna02.lib01.raw.unmapped.bam" is gzip-compressed
    And the file "output/dna01.lib01.raw.unmapped.bam" does not contain "ACGTACGTACGT"
    And the file "m.tsv" contains "0.500000"

  Scenario: a cbt matching a tag but no sample lands in the unassigned bin
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
    # The group knows both cbt tags c1 and c2, but the sample sheet covers only dna01 (=c1). @r1's cbt is
    # routed and lands in dna01.fq; @r2 matches the group's c2 tag yet no sample claims c2, so it is
    # unassigned and its raw segment lands in un.fq. Adding --sample dna02=cbt::c2 would pull @r2 out of
    # the unassigned bin, so its presence there is load-bearing on c2 being unclaimed.
    When I run `unmux reads.fq --group cbt={c1=AAGGTTCC,c2=TTCCAAGG} --group cbt::loc=0:0:8,dist=0 --extract body=0:0:end --template body --sample dna01=cbt::c1 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "dna01.fq" contains "AAGGTTCCACGT"
    And the file "un.fq" contains "TTCCAAGGACGT"
    And the file "dna01.fq" does not contain "TTCCAAGGACGT"

  Scenario: a variable-length cbt group resolves the anchored UMI offset before fan-out
    Given a file "reads.fq" containing:
      """
      @r7
      ATTTTTTTCACACACACACACACACACACGTACGTA
      +
      IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      """
    # The cbt group has a 6 nt and a 7 nt tag, 5'-anchored at loc 0:1:end. This read carries the 7 nt
    # TTTTTTT, so its match ends at position 8; @cbt+19:9 anchors the UMI 19 bp past that end at [27,36) =
    # ACGTACGTA, which the RX tag carries into the fanned-out S.fq. The +19:9 offset is resolved from the
    # matched tag's end, so the 7 nt (vs 6 nt) length does not shift the UMI window.
    When I run `unmux reads.fq --group cbt={c6=GGGGGG,c7=TTTTTTT} --group cbt::loc=0:1:end,dist=0 --extract umi=@cbt+19:9 --extract body=0:0:end --template body --tag RX=umi --tag RX::qual=none --sample S=cbt::c7 --out %sample.fq --unassigned un.fq`
    Then the exit code is 0
    And the file "S.fq" contains "RX:Z:ACGTACGTA"

  Scenario: every output-path placeholder expands across the fan-out
    Given a file "r1.fq" containing:
      """
      @frag1
      AAAAAAAACCCCCCCC
      +
      IIIIIIIIIIIIIIII
      @frag2
      TTTTTTTTGGGGGGGG
      +
      IIIIIIIIIIIIIIII
      """
    And a file "r2.fq" containing:
      """
      @frag1
      GGGGGGGG
      +
      IIIIIIII
      @frag2
      AAAAAAAA
      +
      IIIIIIII
      """
    # Exercises all five output-path placeholders in one run. frag1 matches sample dna01 (sub_sample
    # lib9) on r1's first 8 bp; its two templates (t1 from r1, t2 from r2) fan out to one FASTQ per
    # %ordinal, so --out resolves %pool/%sample/%sub_sample/%ordinal to POOL.dna01.lib9.R1.fq and
    # .R2.fq. frag2 matches no sample, so its raw mates land in --unassigned, where %source (the 0-based
    # input-file index) splits them into unassigned.POOL.0.fq (the r1 mate) and unassigned.POOL.1.fq
    # (the r2 mate). %source is valid only in the --unassigned/--remove bins; %ordinal/%sample/
    # %sub_sample only in --out.
    When I run `unmux r1.fq r2.fq --pool POOL --group bc={s1=AAAAAAAA} --group bc::loc=0:0:8 --extract t1=0:8:end --extract t2=1:0:end --template t1 --template t2 --sample dna01::lib9=bc::s1 --out %pool.%sample.%sub_sample.R%ordinal.fq --unassigned unassigned.%pool.%source.fq`
    Then the exit code is 0
    And a file "POOL.dna01.lib9.R1.fq" exists
    And a file "POOL.dna01.lib9.R2.fq" exists
    And a file "unassigned.POOL.0.fq" exists
    And a file "unassigned.POOL.1.fq" exists
    And the file "POOL.dna01.lib9.R1.fq" contains "CCCCCCCC"
    And the file "POOL.dna01.lib9.R2.fq" contains "GGGGGGGG"
    And the file "unassigned.POOL.0.fq" contains "TTTTTTTTGGGGGGGG"
    And the file "unassigned.POOL.1.fq" contains "AAAAAAAA"
