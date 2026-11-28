# unmux cookbook

Recipes for common unmux jobs. Each one says when to use it, the data it ran on, the exact command, and an excerpt of what unmux printed or wrote, so you can see the result before running it on your own data.

> [!NOTE]
> Every command on this page was run with unmux 0.11.0 and samtools 1.25, and every output is copied from that run. The script [`docs/recipes/generate.py`](recipes/generate.py) fetches the inputs, runs the recipes again, and rewrites this page.

Commands are bash and run in a directory that holds their inputs. Log lines appear without their timestamps, and an excerpt that skips lines marks the gap with …. Header excerpts come from `samtools view --no-PG`, which keeps samtools' own `@PG` line out.

## Contents

- [Concordance at a glance](#concordance-at-a-glance)
- [Conversions](#conversions)
  - [1. Paired FASTQ to an unmapped BAM](#1-paired-fastq-to-an-unmapped-bam)
  - [2. Unmapped BAM back to paired FASTQ](#2-unmapped-bam-back-to-paired-fastq)
  - [3. What plain SAM or BAM pass-through keeps and changes](#3-what-plain-sam-or-bam-pass-through-keeps-and-changes)
  - [4. Split an aligned BAM by read group](#4-split-an-aligned-bam-by-read-group)
  - [5. Join read pairs end to end for CHANGE-seq](#5-join-read-pairs-end-to-end-for-change-seq)
- [Illumina demultiplexing](#illumina-demultiplexing)
  - [6. Dual-index demultiplexing with a sample sheet](#6-dual-index-demultiplexing-with-a-sample-sheet)
  - [7. Undetermined reads and a PhiX bin](#7-undetermined-reads-and-a-phix-bin)
  - [8. Inline barcodes into unmapped BAMs](#8-inline-barcodes-into-unmapped-bams)
- [UMIs and duplex sequencing](#umis-and-duplex-sequencing)
  - [9. Random duplex UMIs into RX](#9-random-duplex-umis-into-rx)
  - [10. Fixed duplex UMIs corrected to a list](#10-fixed-duplex-umis-corrected-to-a-list)
- [Single cell](#single-cell)
  - [11. 10x Chromium cell barcodes and UMIs](#11-10x-chromium-cell-barcodes-and-umis)
  - [12. SPLiT-seq barcodes in three rounds](#12-split-seq-barcodes-in-three-rounds)
  - [13. Keep the reads of called cells with a joined match](#13-keep-the-reads-of-called-cells-with-a-joined-match)
  - [14. CITE-seq antibody barcodes](#14-cite-seq-antibody-barcodes)
- [Long reads](#long-reads)
  - [15. Nanopore native barcodes at both ends](#15-nanopore-native-barcodes-at-both-ends)
  - [16. PacBio barcode pairs, compared with lima](#16-pacbio-barcode-pairs-compared-with-lima)
  - [17. Direct RNA reads with U bases to BAM](#17-direct-rna-reads-with-u-bases-to-bam)
  - [18. Lift SAM tags from FASTQ comments back onto records](#18-lift-sam-tags-from-fastq-comments-back-onto-records)
- [Diagnostics](#diagnostics)
  - [19. Explain one read's assignment with --qc-tag](#19-explain-one-reads-assignment-with---qc-tag)
  - [20. Read the metrics files](#20-read-the-metrics-files)
  - [21. The tag-collision check](#21-the-tag-collision-check)
  - [22. Search windows wider than their tags](#22-search-windows-wider-than-their-tags)
- [Data sources](#data-sources)
- [Regenerating this page](#regenerating-this-page)

## Concordance at a glance

Each row counts, read by read, how unmux's calls compare with another tool's or with the truth, from the same results its recipe checks. A call is the sample, cell barcode or record each side gives a read, and a read pair counts once.

| Recipe | Compared with | Reads compared | Same call | Different call | unmux only | Other only | Neither | Agreement | Main reason for differences |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| [1](#1-paired-fastq-to-an-unmapped-bam) | samtools import | 2,000 | 2,000 | 0 | 0 | 0 | 0 | 100% | none; MAPQ (255 against 0) and the tags left out |
| [2](#2-unmapped-bam-back-to-paired-fastq) | samtools fastq -T RG | 2,000 | 2,000 | 0 | 0 | 0 | 0 | 100% | none; samtools' tab before the tags read as a space |
| [3](#3-what-plain-sam-or-bam-pass-through-keeps-and-changes) | samtools reset | 2 | 2 | 0 | 0 | 0 | 0 | 100% | none; MAPQ and the tags left out |
| [4](#4-split-an-aligned-bam-by-read-group) | samtools split | 17 | 17 | 0 | 0 | 0 | 0 | 100% | none; headers aside |
| [5](#5-join-read-pairs-end-to-end-for-change-seq) | mergeReads() of CHANGE-seq | 3,406 | 3,406 | 0 | 0 | 0 | 0 | 100% | none |
| [6](#6-dual-index-demultiplexing-with-a-sample-sheet) | UDiTaS's demultiplex() rules | 20,000 | 20,000 | 0 | 0 | 0 | 0 | 100% | none |
| [6](#6-dual-index-demultiplexing-with-a-sample-sheet) | fqtk 0.5.0 | 20,000 | 20,000 | 0 | 0 | 0 | 0 | 100% | none; with `--max-mismatches 2 --min-mismatch-delta 1` |
| [7](#7-undetermined-reads-and-a-phix-bin) | the run of each pair, and the UDiTaS rules within fig2b | 27,000 | 19,997 | 0 | 0 | 3 | 7,000 | 99.98% | fig2b pairs that match the PhiX k-mers (3) |
| [8](#8-inline-barcodes-into-unmapped-bams) | the run each read came from | 8,000 | 8,000 | 0 | 0 | 0 | 0 | 100% | none |
| [8](#8-inline-barcodes-into-unmapped-bams) | fqtk 0.5.0 | 8,000 | 8,000 | 0 | 0 | 0 | 0 | 100% | none; with `--max-mismatches 1 --min-mismatch-delta 1` |
| [11](#11-10x-chromium-cell-barcodes-and-umis) | the whole 10x keeplist | 100,000 | 97,291 | 0 | 0 | 0 | 2,709 | 100% | none |
| [12](#12-split-seq-barcodes-in-three-rounds) | splitcode 0.31.6, tutorial files | 10,000 | 3,473 | 6 | 3,637 | 211 | 2,673 | 47.40% | 43 of the tutorial's 192 round-1 barcodes are this library's |
| [12](#12-split-seq-barcodes-in-three-rounds) | splitcode 0.31.6, this library's lists | 10,000 | 7,116 | 0 | 0 | 0 | 2,884 | 100% | none |
| [15](#15-nanopore-native-barcodes-at-both-ends) | the run each read came from | 4,800 | 4,748 | 0 | 0 | 52 | 0 | 98.91% | reads left unclassified (52) |
| [16](#16-pacbio-barcode-pairs-compared-with-lima) | lima 26.2.1 | 1,922 | 1,253 | 0 | 11 | 8 | 650 | 98.50% | lima's score filter (11); one end beyond unmux's budget (8) |

Agreement is the share of same calls among the reads that at least one side calls; reads that neither calls are left out of it.

## Conversions

Format conversions with nothing to demultiplex: unmux as a converter between FASTQ, SAM and BAM, compared with samtools.

### 1. Paired FASTQ to an unmapped BAM

Use this to store raw read pairs as one unmapped BAM (uBAM), the input many alignment workflows such as GATK's expect, with a read group on every record.

**Data:** the first 2,000 read pairs of SRA run [SRR390728](https://www.ncbi.nlm.nih.gov/sra/SRR390728), fetched with fg-sra:

```bash
fg-sra fastq --min-spot-id 1 --max-spot-id 2000 --r1 SRR390728_1.fastq --r2 SRR390728_2.fastq SRR390728
```

```bash
unmux SRR390728_1.fastq.gz SRR390728_2.fastq.gz \
  --rg-tag PL=ILLUMINA \
  --out SRR390728.unmapped.bam
```

The log reports what unmux read and wrote:

```text
[INFO  unmux::demux] pool=None: 2 input(s), 1 output(s), 0 group(s), 0 extract(s), 0 template(s), 0 SAM tag(s)
[INFO  unmux::demux]   input  0: SRR390728_1.fastq.gz (FASTQ)
[INFO  unmux::demux]   input  1: SRR390728_2.fastq.gz (FASTQ)
[INFO  unmux::demux]   output 0: SRR390728.unmapped.bam (BAM)
[INFO  unmux::demux] threads=8: 8 worker(s) read, match, and compress; 0 decode BAM or CRAM input
[INFO  unmux::demux] Demultiplexed 2,000 records
[INFO  unmux::metrics] metrics pool=SRR390728: 2,000 reads (0 assigned, 2,000 pass-through, 0 unassigned, 0 removed) frac_assigned=0.0000 frac_bases_unextracted=0.0000 mean_read_len_by_input=36,36
```

```bash
samtools view --no-PG -h SRR390728.unmapped.bam | head -n 5
```

```text
@HD	VN:1.6
@RG	ID:SRR390728	SM:SRR390728	LB:SRR390728	PL:ILLUMINA
@PG	ID:unmux	PN:unmux	VN:0.11.0	CL:unmux SRR390728_1.fastq.gz SRR390728_2.fastq.gz --rg-tag PL=ILLUMINA --out SRR390728.unmapped.bam
SRR390728.1	77	*	0	255	*	*	0	0	CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGC	;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142	RG:Z:SRR390728
SRR390728.1	141	*	0	255	*	*	0	0	GATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(	RG:Z:SRR390728
```

Every record carries `RG:Z:SRR390728`. The pool, which becomes the read group's ID, SM and LB, defaults to the common stem of the input names; set it with `--pool`.

Each of the 4,000 records has the name, bases and qualities of its FASTQ read, flags 77 and 141 for the two reads of a pair, and no tag but RG.

Compared with `samtools import`, the records differ only in MAPQ and the read group:

```bash
samtools import -1 SRR390728_1.fastq.gz -2 SRR390728_2.fastq.gz -o samtools.unmapped.bam
diff <(samtools view samtools.unmapped.bam | cut -f 1-4,6-11) \
  <(samtools view SRR390728.unmapped.bam | cut -f 1-4,6-11) && echo identical
```

```text
identical
```

> [!NOTE]
> MAPQ is left out of the comparison: unmux writes 255 (unavailable) on unmapped records, where `samtools import` writes 0. All 4,000 records match otherwise.

### 2. Unmapped BAM back to paired FASTQ

Use this to hand reads stored in an unmapped BAM to a tool that wants FASTQ. The `%ordinal` placeholder writes one file per read of the pair.

**Data:** the unmapped BAM from [recipe 1](#1-paired-fastq-to-an-unmapped-bam).

```bash
unmux SRR390728.unmapped.bam --out "SRR390728.R%ordinal.fastq.gz"
```

```text
SRR390728.R1.fastq.gz
SRR390728.R2.fastq.gz
```

```bash
gzip -dc SRR390728.R1.fastq.gz | head -n 4
```

```text
@SRR390728.1 RG:Z:SRR390728
CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGC
+
;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142
```

The tags of each record, here the read group, go into the FASTQ comment. The files match `samtools fastq -T RG` once its tab before the tags is read as a space:

```bash
samtools fastq -T RG -1 samtools.R1.fastq.gz -2 samtools.R2.fastq.gz SRR390728.unmapped.bam
for r in 1 2; do
  diff <(gzip -dc samtools.R$r.fastq.gz | tr '\t' ' ') <(gzip -dc SRR390728.R$r.fastq.gz) && echo "R$r identical"
done
```

```text
R1 identical
R2 identical
```

Bases and qualities equal the original FASTQ for all 2,000 pairs.

### 3. What plain SAM or BAM pass-through keeps and changes

Use this to get an unmapped copy of aligned reads, to realign them or share them without the alignment, as `samtools reset` does.

**Data:** two aligned Nanopore reads with modified-base calls (MM/ML), one on each strand ([fwd_rev_modbase_records.bam](https://raw.githubusercontent.com/nanoporetech/modkit/5cecc3fb3a9336068d9e3c68d5c08d678153dd2c/tests/resources/fwd_rev_modbase_records.bam)).

```bash
unmux modbam.aligned.bam --out modbam.unmapped.bam
```

The input's header without its 34 `@SQ` lines, the first six columns of its records, and the tags of the first record:

```bash
samtools view --no-PG -H modbam.aligned.bam | grep -v '^@SQ' | cut -c 1-100
samtools view modbam.aligned.bam | cut -f 1-6
samtools view modbam.aligned.bam | head -n 1 | cut -f 12- | tr '\t' '\n' | cut -d : -f 1 | paste -s -d ' ' -
```

```text
@HD	VN:1.5	SO:coordinate	ob:0.0.2
@PG	ID:basecaller	PN:bonito	VN:0.6.1	CL:bonito basecaller dna_r10.4.1_e8.2_400bps_hac@v3.5.2 subsets
@PG	ID:aligner	PN:minimap2	VN:2.24	DS:mappy
@PG	ID:samtools	PN:samtools	PP:aligner	VN:1.16.1	CL:samtools sort bc_anchored_hac_synthetic_ground_t
@RG	ID:5598049b1b3264566b162bf035344e7ec610d608_dna_r10.4.1_e8.2_400bps_hac@v3.5.2	PL:ONT	DT:2022-08
068ce426-129e-4870-bd34-16cd78edaa43	0	oligo_1512_adapters	1	9	9S16M1D45M1D13M1D18M1I13M1D4M4S
e3a491b7-4e53-48f2-b75f-fa62a74ab729	16	oligo_1512_adapters	16	44	4S5M2D72M1D8M2D10M1D6M1I30M12S
NM MD RG qs ns ts mx ch st du rn f5 sm sd sv mv MM ML
```

The same for the output:

```bash
samtools view --no-PG -H modbam.unmapped.bam
samtools view modbam.unmapped.bam | cut -f 1-6
samtools view modbam.unmapped.bam | head -n 1 | cut -f 12- | tr '\t' '\n' | cut -d : -f 1 | paste -s -d ' ' -
```

```text
@HD	VN:1.6
@RG	ID:modbam.aligned	SM:modbam.aligned	LB:modbam.aligned
@PG	ID:unmux	PN:unmux	VN:0.11.0	CL:unmux modbam.aligned.bam --out modbam.unmapped.bam
068ce426-129e-4870-bd34-16cd78edaa43	4	*	0	255	*
e3a491b7-4e53-48f2-b75f-fa62a74ab729	4	*	0	255	*
RG qs ns ts mx ch st du rn f5 sm sd sv mv MM ML
```

What changes:

- Each record becomes unmapped: flag 4, no position or CIGAR, MAPQ 255. The read on the reverse strand (flag 16) gets back its sequenced orientation, bases and qualities.
- Tags that describe the alignment are dropped, here MD, NM; per-read tags such as MM, ML and the basecaller's tags are kept.
- The header keeps `@HD VN` and gets the pool's read group and an unmux `@PG` line. The input `@RG`, `@PG`, `@CO` and `@SQ` lines are not kept, and each record's RG:Z becomes the pool.

Apart from MAPQ and RG:Z, the records equal those of `samtools reset`:

```bash
samtools reset -o reset.bam modbam.aligned.bam
diff <(samtools view reset.bam | cut -f 1-4,6-11) <(samtools view modbam.unmapped.bam | cut -f 1-4,6-11) && echo identical
```

```text
identical
```

The generator also checks that MM and ML are unchanged on both records and that the reverse-strand read is the reverse complement of its aligned SEQ.

> [!TIP]
> To keep the alignments and split them into files instead, use a read-group split, as in [recipe 4](#4-split-an-aligned-bam-by-read-group).

### 4. Split an aligned BAM by read group

Use this to separate an aligned SAM or BAM that holds several read groups, such as lanes or samples, into one file per read group with the alignments kept whole.

**Data:** 34 aligned records in two read groups, some with no read group or one missing from the header ([split.sam](https://raw.githubusercontent.com/samtools/samtools/7162950d0988480faeddeb49f03520737e13a9c4/test/split/split.sam)).

```bash
unmux split.sam \
  --group "rg=@RG" \
  --sample-from-group rg \
  --out "by_rg/%sample.bam" \
  --unassigned by_rg/unknown.bam
```

```bash
for f in by_rg/*.bam; do echo "$f $(samtools view -c "$f")"; done
```

```text
by_rg/grp1.bam 14
by_rg/grp2.bam 14
by_rg/unknown.bam 6
```

```bash
samtools view --no-PG -H by_rg/grp1.bam | grep -v '^@CO'
```

```text
@HD	VN:1.4	SO:coordinate
@SQ	SN:ref1	LN:56	M5:08c04d512d4797d9ba2a156c1daba468
@SQ	SN:ref2	LN:60	M5:7c35feac7036c1cdef3bee0cc4b21437
@SQ	SN:ref3_unused	LN:70	M5:5fdd18c2c6ecac4838996d029bf395b5
@RG	ID:grp1	DS:Group 1	LB:Library 1	SM:Sample
@PG	ID:prog1	PN:emacs	CL:emacs	VN:23.1.1
@PG	ID:unmux	PN:unmux	VN:0.11.0	CL:unmux split.sam --group rg=@RG --sample-from-group rg --out by_rg/%sample.bam --unassigned by_rg/unknown.bam	PP:prog1
```

Each file has only its own `@RG` line, and the records are those `samtools split` writes:

```bash
samtools split --no-PG -u unknown.samtools.bam -f "%!.samtools.bam" split.sam
for rg in grp1 grp2 unknown; do
  diff <(samtools view by_rg/$rg.bam) <(samtools view $rg.samtools.bam) && echo "$rg identical"
done
```

```text
grp1 identical
grp2 identical
unknown identical
```

The headers differ only in where `@SQ` sits and the added unmux `@PG` line:

```bash
diff <(samtools view --no-PG -H by_rg/grp1.bam) <(samtools view --no-PG -H grp1.samtools.bam)
```

```diff
2,4d1
< @SQ	SN:ref1	LN:56	M5:08c04d512d4797d9ba2a156c1daba468
< @SQ	SN:ref2	LN:60	M5:7c35feac7036c1cdef3bee0cc4b21437
< @SQ	SN:ref3_unused	LN:70	M5:5fdd18c2c6ecac4838996d029bf395b5
7d3
< @PG	ID:unmux	PN:unmux	VN:0.11.0	CL:unmux split.sam --group rg=@RG --sample-from-group rg --out by_rg/%sample.bam --unassigned by_rg/unknown.bam	PP:prog1
28a25,27
> @SQ	SN:ref1	LN:56	M5:08c04d512d4797d9ba2a156c1daba468
> @SQ	SN:ref2	LN:60	M5:7c35feac7036c1cdef3bee0cc4b21437
> @SQ	SN:ref3_unused	LN:70	M5:5fdd18c2c6ecac4838996d029bf395b5
```

> [!TIP]
> With `--group "rg=@RG::SM"` the files are per sample, merging the read groups of one sample, and `@RG::SM::LB` adds a library level for `--out "%sample.%sub_sample.bam"`.

### 5. Join read pairs end to end for CHANGE-seq

Use this to turn read pairs into the single reads that CIRCLE-seq and CHANGE-seq align: the pipeline's `mergeReads()` writes the reverse complement of R1 followed by R2, with R1's qualities reversed, so both reads of a cleaved site line up on one strand.

**Data:** the treated (2,439 pairs) and control (967 pairs) test libraries of the CHANGE-seq pipeline ([TEST.r1.fastq.gz](https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf/test/data/input/TEST.r1.fastq.gz), [TEST.r2.fastq.gz](https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf/test/data/input/TEST.r2.fastq.gz), [TEST_control.r1.fastq.gz](https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf/test/data/input/TEST_control.r1.fastq.gz), [TEST_control.r2.fastq.gz](https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf/test/data/input/TEST_control.r2.fastq.gz)).

```bash
gzip -dc TEST.r1.fastq.gz | head -n 2 | cut -c 1-70
gzip -dc TEST.r2.fastq.gz | head -n 2 | cut -c 1-70
```

```text
@M01326:241:000000000-AHF83:1:1101:16662:2289 1:N:0:3
CCGTCTGGGACTCTGCATCTTTCAGACCCACCTGTCCTCTTGAACTTAGCCCTCCCTAACTTGGGAGCAA
@M01326:241:000000000-AHF83:1:1101:16662:2289 2:N:0:3
CTTGAAAAACTCTTTTTGTGGAATTTGCAAGTGGAGATTTCAGCCGCTTTGAGGTCCATGGTAGAAAAGG
```

```bash
unmux TEST.r1.fastq.gz TEST.r2.fastq.gz \
  --extract "r1=0:0:end" --extract "r2=1:0:end" \
  --template "~r1+r2" \
  --out TEST_merged.fastq.gz
unmux TEST_control.r1.fastq.gz TEST_control.r2.fastq.gz \
  --extract "r1=0:0:end" --extract "r2=1:0:end" \
  --template "~r1+r2" \
  --out control_TEST_merged.fastq.gz
```

```bash
gzip -dc TEST_merged.fastq.gz | head -n 2 | cut -c 1-70
```

```text
@M01326:241:000000000-AHF83:1:1101:16662:2289
GCACATTTAGAGTGTGCCTGTGTGTCTGTGTTGGAGAGGTCTGCTGCAGAGCCCAGGGCATCCCCCAACC
```

The `~` reverse-complements the R1 stream, qualities reversed with it, and `+` joins the two streams into one read. Every merged read of both libraries equals what [mergeReads.py](https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf/changeseq/mergeReads.py) writes for its pair, bases and qualities, and keeps the pair's name without its comment.

## Illumina demultiplexing

Short-read runs split by index reads or inline barcodes, with reads of known origin so each result can be scored.

### 6. Dual-index demultiplexing with a sample sheet

Use this to split a pooled Illumina run into per-sample unmapped BAMs when the index reads are in their own files, here with a UMI after the i5 index in I2, as in UDiTaS (Giannoukos et al. 2018), a method for measuring gene-editing outcomes.

**Data:** the first 20,000 reads of the UDiTaS fig2b MiSeq run, with its sample sheet ([Undetermined_S0_L001_R1_001.fastq](https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2b/Undetermined_S0_L001_R1_001.fastq.gz), [Undetermined_S0_L001_R2_001.fastq](https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2b/Undetermined_S0_L001_R2_001.fastq.gz), [Undetermined_S0_L001_I1_001.fastq](https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2b/Undetermined_S0_L001_I1_001.fastq.gz), [Undetermined_S0_L001_I2_001.fastq](https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2b/Undetermined_S0_L001_I2_001.fastq.gz), [sample_info.csv](https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2b/sample_info.csv)). Its 10 samples share the i5 index N504 and each has its own i7. I1 is the 8 bp i7 index; I2 is the 8 bp i5 index followed by a 10 bp UMI:

```bash
gzip -dc Undetermined_S0_L001_I1_001.fastq.gz | head -n 4
gzip -dc Undetermined_S0_L001_I2_001.fastq.gz | head -n 4
```

```text
@M03958:36:000000000-B85D5:1:1101:9410:1062 1:N:0:0
NATCTCCT
+
#8BCCGGG
@M03958:36:000000000-B85D5:1:1101:9410:1062 2:N:0:0
NGAGTAGAGGATTTTTGT
+
#8ACCGGGGGGGGGGGGG
```

Build a sample sheet and a tag file per index from the run's `sample_info.csv`:

```bash
awk -F , 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; print "sample\ti7\ti5"; next }
  { print $c["index_I1"] "_" $c["index_I2"] "\t" $c["index_I1"] "\t" $c["index_I2"] }' sample_info.csv > samples.tsv
awk -F , 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; print "id\tseq"; next }
  !seen[$c["index_I1"]]++ { print $c["index_I1"] "\t" $c["barcode_I1"] }' sample_info.csv > i7.tsv
awk -F , 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; print "id\tseq"; next }
  !seen[$c["index_I2"]]++ { print $c["index_I2"] "\t" $c["barcode_I2"] }' sample_info.csv > i5.tsv
```

```bash
head -n 4 samples.tsv
head -n 3 i7.tsv
cat i5.tsv
```

```tsv
sample	i7	i5
N701_N504	N701	N504
N702_N504	N702	N504
N703_N504	N703	N504
id	seq
N701	ATTCCGCT
N702	GCATGATC
id	seq
N504	AGAGTAGA
```

```bash
unmux Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz \
  Undetermined_S0_L001_I1_001.fastq.gz Undetermined_S0_L001_I2_001.fastq.gz \
  --pool fig2b \
  --group "i7=i7.tsv" --group "i7::loc=2:0:8,dist=1" \
  --group "i5=i5.tsv" --group "i5::loc=3:0:8,dist=1" \
  --extract "r1=0:0:end" --extract "r2=1:0:end" --extract "umi=3:8:end" \
  --template r1 r2 \
  --tag "RX=umi" \
  --sample-sheet samples.tsv \
  --rg-tag PL=ILLUMINA \
  --out "bam/%sample.unmapped.bam" \
  --metrics-per-sample metrics.per_sample.tsv
```

Each index group searches only its index read and corrects one substitution, and a read is assigned when its i7 and i5 match one row of the sheet. The UMI is I2 from base 8 on (counting from 0), written to RX with its qualities in QX.

```bash
for f in bam/*.bam; do echo "$f $(samtools view -c -f 64 "$f")"; done
```

```text
bam/N701_N504.unmapped.bam 1765
bam/N702_N504.unmapped.bam 3221
bam/N703_N504.unmapped.bam 2383
bam/N704_N504.unmapped.bam 2398
bam/N705_N504.unmapped.bam 1330
bam/N706_N504.unmapped.bam 1847
bam/N707_N504.unmapped.bam 1322
bam/N710_N504.unmapped.bam 1058
bam/N711_N504.unmapped.bam 2360
bam/N712_N504.unmapped.bam 2316
```

```bash
samtools view bam/N702_N504.unmapped.bam | head -n 2 | cut -f 1,2,12-
```

```text
M03958:36:000000000-B85D5:1:1101:9738:1064	77	RX:Z:GCTTCGGTTA	QX:Z:GGGGGGGGGG	RG:Z:N702_N504
M03958:36:000000000-B85D5:1:1101:9738:1064	141	RX:Z:GCTTCGGTTA	QX:Z:GGGGGGGGGG	RG:Z:N702_N504
```

20,000 of 20,000 reads were assigned. To check every call, the generator reimplements the rules of `demultiplex()` in UDiTaS's [uditas_helpers.py](https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/uditas/uditas_helpers.py): each index exact, or one mismatch from its closest barcode, and the pair listed in the sheet. unmux makes the same call on all 20,000 reads. UDiTaS also masks bases below quality 12 to N before matching, which unmux does not do; no read here has such a base that pushes an index past one mismatch, and on a run that has them, UDiTaS would leave those reads unmatched where unmux assigns them. RX and QX equal I2's bases and qualities from base 8 on, on every record.

> [!NOTE]
> UDiTaS writes the UMI with low-quality bases masked to N. unmux writes the bases as read and their qualities to QX, so a later step can weigh them.

**Compared with fqtk 0.5.0.** fqtk reads each read's barcodes from fixed positions given as read structures, and each sample's i7 and i5 joined in one `barcode` column:

```bash
awk 'BEGIN { print "sample_id\tbarcode" } FNR == 1 { next }
  FILENAME == "i7.tsv" { i7[$1] = $2; next } FILENAME == "i5.tsv" { i5[$1] = $2; next }
  { print $1 "\t" i7[$2] i5[$3] }' i7.tsv i5.tsv samples.tsv > fqtk.samples.tsv
head -n 3 fqtk.samples.tsv
```

```tsv
sample_id	barcode
N701_N504	ATTCCGCTAGAGTAGA
N702_N504	GCATGATCAGAGTAGA
```

```bash
fqtk demux \
  --inputs Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz \
    Undetermined_S0_L001_I1_001.fastq.gz Undetermined_S0_L001_I2_001.fastq.gz \
  --read-structures +T +T 8B 8B+M \
  --sample-metadata fqtk.samples.tsv \
  --max-mismatches 2 --min-mismatch-delta 1 \
  --umi-in-name false \
  --output fqtk
```

fqtk counts mismatches over the joined barcode, where unmux gives each index a budget of its own, so `--max-mismatches 2` admits one mismatch in each index, and `--min-mismatch-delta 1` takes the best sample when it beats the next by one; with `--umi-in-name false` the read names stay as they are, so each read can be matched with unmux's. Per read, against unmux:

| Reads | these settings | fqtk's defaults |
| --- | ---: | ---: |
| Same call from both | 20,000 | 19,828 |
| Different calls | 0 | 0 |
| A call from unmux only | 0 | 172 |
| A call from fqtk only | 0 | 0 |
| No call from either | 0 | 0 |

With these settings fqtk makes unmux's call on all 20,000 reads. They would also admit two mismatches in one index, which unmux does not; no read here has that. With fqtk's defaults, one mismatch and a lead of two, the 172 reads with a mismatch in each index are left unmatched.

How unmux expresses fqtk's settings:

| fqtk | unmux |
| --- | --- |
| `--read-structures`, with `B` for the barcode bases of each read | a group's `loc`, such as `i7::loc=2:0:8` |
| `--sample-metadata`, each sample's barcodes joined | `--sample-sheet`, a column per group, and a tag file per group |
| `--max-mismatches`, over all the barcode bases | `dist`, for each group |
| `--min-mismatch-delta` | `mode=nearest`, with `delta` |
| `M` segments, the UMI | `--extract` and `--tag "RX=umi"` |
| `T` segments, written to one FASTQ per sample and read | `--template` and `--out` |
| `--unmatched-prefix` | `--unassigned` |
| `demux-metrics.txt` | `--metrics-per-sample` and `--metrics-summary` |
| IUPAC codes in expected barcodes | IUPAC codes in tags |
| `--header-format illumina` and `--umi-in-name` | no equivalent |

What fqtk does that unmux does not: it can write FASTQ headers in the bcl2fastq style, with the sample barcode in the comment and the UMI in the read name. unmux writes tags such as BC and RX as SAM tags instead, in the record or in the FASTQ comment.

### 7. Undetermined reads and a PhiX bin

Use this to keep the reads that match no sample, the way `bcl-convert` writes Undetermined files, and to move spike-in reads such as PhiX into a bin of their own.

**Data:** a lane pooled from three real runs, so the origin of every read is known: the 20,000 fig2b reads of [recipe 6](#6-dual-index-demultiplexing-with-a-sample-sheet), the first 5,000 reads of UDiTaS's fig2c run, a library with an i5 that fig2b's sheet does not list ([Undetermined_S0_L001_R1_001.fastq](https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2c/Undetermined_S0_L001_R1_001.fastq.gz) and its other reads), and 2,000 read pairs of the PhiX control run [SRR27862880](https://www.ncbi.nlm.nih.gov/sra/SRR27862880), fetched with fg-sra:

```bash
fg-sra fastq --min-spot-id 1 --max-spot-id 2000 --r1 SRR27862880_1.fastq --r2 SRR27862880_2.fastq --defline '$sn $ri:N:0:$sg' SRR27862880
```

SRA kept the PhiX run's index only in its spot group, `GGGGGGGG+AGATCTCG` (the i7 reads as G with no index on a two-color instrument), so the pool's I1 and I2 for those pairs are built from it, with every quality F. That is the only constructed part.

The PhiX tags are 24-mers taken every 12 bases of the PhiX174 genome, [NC_001422.1](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text):

```bash
awk '!/^>/ { seq = seq $0 } END { print "id\tseq"
  for (i = 1; i + 23 <= length(seq); i += 12) print "phix_" i "\t" substr(seq, i, 24) }' NC_001422.1.fasta > phix.tsv
```

```bash
unmux R1.fastq.gz R2.fastq.gz I1.fastq.gz I2.fastq.gz \
  --pool lane1 \
  --group "i7=i7.tsv" --group "i7::loc=2:0:8,dist=1" \
  --group "i5=i5.tsv" --group "i5::loc=3:0:8,dist=1" \
  --group "phix=phix.tsv" --group "phix::loc=0:0:end,dist=1,bothStrands=true" \
  --extract "r1=0:0:end" --extract "r2=1:0:end" \
  --template r1 r2 \
  --sample-sheet samples.tsv \
  --out "fastq/%sample.R%ordinal.fastq.gz" \
  --remove "phix=removed/PhiX.%source.fastq.gz" \
  --unassigned "fastq/Undetermined.%source.fastq.gz" \
  --metrics-per-sample metrics.per_sample.tsv \
  --metrics-summary metrics.summary.tsv \
  --metrics-unassigned metrics.unassigned.tsv
```

The `phix` group looks for any of the 24-mers on either strand of R1. The `--remove` rule sends the reads it matches to their own files, and `--unassigned` keeps the rest that match no sample, one file per input since `%source` is the input's number from 0 (R1, R2, I1, I2).

```bash
for f in removed/*.fastq.gz fastq/Undetermined.*.fastq.gz; do
  echo "$f $(gzip -dc "$f" | wc -l | awk '{print $1 / 4}')"
done
```

```text
removed/PhiX.0.fastq.gz 2003
removed/PhiX.1.fastq.gz 2003
removed/PhiX.2.fastq.gz 2003
removed/PhiX.3.fastq.gz 2003
fastq/Undetermined.0.fastq.gz 5000
fastq/Undetermined.1.fastq.gz 5000
fastq/Undetermined.2.fastq.gz 5000
fastq/Undetermined.3.fastq.gz 5000
```

```bash
cat metrics.unassigned.tsv
```

```tsv
pool	reason	group	stream	reads	frac_of_pool
lane1	no_sample			5000	0.185185
```

Checked against the origin of every pair: the PhiX bin holds all 2,000 pairs of the PhiX run and 3 UDiTaS pairs, 1 whose R1 holds an exact PhiX 24-mer and 2 with one a substitution away from a tag. Undetermined holds the other 5,000 fig2c pairs and nothing else, and the samples got the other 19,997 fig2b pairs, each in the sample [recipe 6](#6-dual-index-demultiplexing-with-a-sample-sheet) gives it.

> [!TIP]
> With `dist=0` in the `phix` group, the bin takes 1,994 PhiX pairs and 0 others: an exact match spares those UDiTaS reads but misses PhiX reads in which every tag that fits has an error.

> [!WARNING]
> Do not add `minFindsPerGroup` to the index groups here. A read that fails a find constraint is unassigned before a `--remove` rule sees it, so with `minFindsPerGroup=1` on i7 and i5 all the PhiX run's pairs go to Undetermined and none to the PhiX bin.

> [!TIP]
> Pass the `--out` pattern itself to `--unassigned` to keep unassigned reads in the `--out` file instead of a file of their own.

### 8. Inline barcodes into unmapped BAMs

Use this when the sample barcode is the first bases of the read rather than an index read, as in RAD-seq: unmux checks the restriction site after it, demultiplexes on the barcode, trims it off, and keeps it in BC with its qualities in QT.

**Data:** 250 reads of each of the 32 RAD-seq runs of ground beetles in study [DRP001067](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB2148) (Takahashi et al. 2014), one run per sample, fetched with fg-sra and pooled one after another, so the sample of every read is known. Each 101 bp read starts with the sample's 5 bp barcode and the PstI site remnant TGCAG; the barcodes are the ones each run's BioSample record (SAMD00008878 to SAMD00008909) gives, at least 3 substitutions apart.

```bash
for run in DRR0219{59..90}; do
  fg-sra fastq --min-spot-id 1 --max-spot-id 250 --unpaired "$run.fastq" "$run"
done
```

```bash
head -n 4 barcodes.tsv
gzip -dc pool.fastq.gz | head -n 2
```

```text
id	seq
CR01	AAAAA
CR02	AACCC
CR03	AAGGG
@DRR021959.1
AAGAATGCAGGATATGGGATACAGACTCAGCCTTGTCGCAACCCGCACGGCACATCCGCTGCTCTGCTGGATTGCTGGGGAGTCCCCTGCATGGGAGGTTA
```

```bash
unmux pool.fastq.gz \
  --pool D0BTPACXX.8 \
  --group "bc=barcodes.tsv" --group "bc::loc=0:0:5,dist=1" \
  --group "site={PstI=TGCAG}" --group "site::loc=0:5:10,minFindsPerGroup=1" \
  --extract "barcode=0:0:5" --extract "insert=0:5:end" \
  --template insert \
  --tag "BC=barcode" \
  --sample-from-group bc \
  --rg-tag PL=ILLUMINA \
  --out "bam/%sample.unmapped.bam"
```

The `site` group requires the remnant right after the barcode, so a read without it would be unassigned; the insert keeps the remnant, as RAD pipelines expect.

A read whose barcode had one substitution:

```bash
samtools view bam/CR01.unmapped.bam | awk '$1 == "DRR021959.1"' | cut -f 1,2,12-
```

```text
DRR021959.1	4	BC:Z:AAGAA	QT:Z:CCCFD	RG:Z:CR01
```

All 8,000 reads went to their run's sample and none to a wrong one; 8,000 carry the remnant exactly. On every record BC and QT equal the read's first 5 bases and qualities, and the record holds the rest of the read.

> [!NOTE]
> The submitter deposited these runs already split by barcode, allowing one substitution, so no read here has more errors than `dist=1` corrects.

> [!TIP]
> With `dist=0`, 7,829 reads are assigned and the 171 with a substitution in the barcode are not.

**Compared with fqtk 0.5.0.** The same pool through fqtk, with one mismatch, as `dist=1` allows, and a lead of one over the next barcode; [recipe 6](#6-dual-index-demultiplexing-with-a-sample-sheet) gives how its settings map to unmux's:

```bash
awk 'BEGIN { print "sample_id\tbarcode" } NR > 1 { print $1 "\t" $2 }' barcodes.tsv > fqtk.samples.tsv
fqtk demux \
  --inputs pool.fastq.gz \
  --read-structures 5B+T \
  --sample-metadata fqtk.samples.tsv \
  --max-mismatches 1 --min-mismatch-delta 1 \
  --output fqtk
```

| Reads | fqtk |
| --- | ---: |
| Same call from both | 8,000 |
| Different calls | 0 |
| A call from unmux only | 0 |
| A call from fqtk only | 0 |
| No call from either | 0 |

fqtk makes unmux's call on all 8,000 reads. It reads the barcode from a fixed read structure and has no check like the `site` group's; every read here has the remnant, so that check changes no call.

## UMIs and duplex sequencing

Unique molecular identifiers (UMIs) into the RX tag for duplex consensus calling, on reads from the nf-core/fastquorum tests.

### 9. Random duplex UMIs into RX

Use this for duplex sequencing libraries with a random UMI at the start of each read: one RX tag gets both UMIs, R1's first, the UMIs and spacers come off the reads, and the unmapped BAM is ready for fgbio's duplex workflow.

**Data:** the 9,901 read pairs of SRR6109255 in the nf-core/fastquorum test set ([SRR6109255_1.fastq.gz](https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/fastqs/tiny/SRR6109255_1.fastq.gz), [SRR6109255_2.fastq.gz](https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/fastqs/tiny/SRR6109255_2.fastq.gz), [samplesheet.tiny.csv](https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/samplesheets/samplesheet.tiny.csv)). The sample sheet gives the read structure `10M1S+T 10M1S+T`: in each read a 10 bp UMI and a 1 bp spacer, then the template. The spacer is a T in 98.2% of reads.

```bash
gzip -dc SRR6109255_1.fastq.gz | head -n 4 | cut -c 1-60
gzip -dc SRR6109255_2.fastq.gz | head -n 4 | cut -c 1-60
```

```text
@SRR6109255.184 184/1
CGCCATACAATAATACTGAGGAAAAGGGAGCATTAAGATCGGAAGAGCACACGTCTGAAC
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
@SRR6109255.184 184/2
TAATGCTCCCTTTTCCTCAGTATTATTGTATGGCGAGATCGGAAGAGCGTCGTGTAGGGA
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGFGGGGG
```

```bash
unmux SRR6109255_1.fastq.gz SRR6109255_2.fastq.gz \
  --extract "umi1=0:0:10" --extract "umi2=1:0:10" \
  --extract "r1=0:11:end" --extract "r2=1:11:end" \
  --template r1 r2 \
  --tag "RX=umi1+umi2" \
  --out SRR6109255.unmapped.bam
```

```bash
samtools view SRR6109255.unmapped.bam | head -n 2 | cut -f 1,2,12-
samtools view SRR6109255.unmapped.bam | head -n 2 | cut -f 10 | cut -c 1-49
```

```text
SRR6109255.184	77	RX:Z:CGCCATACAA-TAATGCTCCC	QX:Z:CCCCCGGGGG CCCCCGGGGG	RG:Z:SRR6109255
SRR6109255.184	141	RX:Z:CGCCATACAA-TAATGCTCCC	QX:Z:CCCCCGGGGG CCCCCGGGGG	RG:Z:SRR6109255
AATACTGAGGAAAAGGGAGCATTAAGATCGGAAGAGCACACGTCTGAAC
TTTCCTCAGTATTATTGTATGGCGAGATCGGAAGAGCGTCGTGTAGGGA
```

On all 9,901 pairs, RX is R1's first 10 bases, a hyphen, then R2's, QX holds their qualities joined by a space, and both reads start after their spacer, at base 11 counting from 0. The pool, and so the read group, is the common stem of the input names, `SRR6109255`.

In this sample 870 of the 8,875 UMI pairs also appear swapped, as B-A for an A-B: 435 molecules read from both strands, in 1,042 read pairs. These are the pairs the paired strategy joins:

```bash
samtools view -f 64 SRR6109255.unmapped.bam | grep -o 'RX:Z:[^[:space:]]*' | cut -c 6- \
  | sort -u > umi_pairs.txt
awk 'BEGIN { FS = "-" } { print $2 "-" $1 }' umi_pairs.txt | sort | comm -12 - umi_pairs.txt \
  | awk 'END { print NR }'
```

```text
870
```

> [!IMPORTANT]
> Keep R1's UMI first. The read pairs from one strand of a molecule carry A-B and those from the other strand B-A, since their R1 starts from the other end. fgbio's paired strategy joins the two by that swap and tells the strands apart by it; sorting the two UMIs, or putting R2's first, would lose the strand.

> [!TIP]
> Next, align the reads, bring the tags back onto the alignments with fgbio ZipperBams, and run fgbio GroupReadsByUmi with `--strategy paired`, then CallDuplexConsensusReads.

### 10. Fixed duplex UMIs corrected to a list

Use this for duplex kits with a fixed set of UMIs, such as IDT xGen cfDNA: each read's UMI is corrected to the list, RX gets the corrected pair, and OX keeps the bases as read with their qualities in BZ.

**Data:** the nf-core/fastquorum test set for fixed UMIs and the UMI list ([synthetic_nonrandom_umi_1.fastq.gz](https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/fastqs/synthetic-nonrandom-umi/synthetic_nonrandom_umi_1.fastq.gz), [synthetic_nonrandom_umi_2.fastq.gz](https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/fastqs/synthetic-nonrandom-umi/synthetic_nonrandom_umi_2.fastq.gz), [idt-xgen-cfdna-umis.txt](https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/umis/idt-xgen-cfdna-umis.txt), [samplesheet.nonrandom.csv](https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/samplesheets/samplesheet.nonrandom.csv)). Its 9,901 read pairs are those of [recipe 9](#9-random-duplex-umis-into-rx) with the first 11 bases of each read, the random UMI and the base after it, replaced by 8 bases of quality I: one of the 32 IDT xGen cfDNA UMIs or, in 1,999 of the 19,802 reads, one substitution away from one. The rest of each read is unchanged. The read structure is `8M+T 8M+T`: an 8 bp UMI, then the template.

```bash
gzip -dc synthetic_nonrandom_umi_1.fastq.gz | head -n 4 | cut -c 1-60
```

```text
@SRR6109255.184 184/1
GATGTGTGAATACTGAGGAAAAGGGAGCATTAAGATCGGAAGAGCACACGTCTGAACTCC
+
IIIIIIIIGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
```

The list holds one UMI per line, so give it the `id` and `seq` columns unmux reads:

```bash
awk 'BEGIN { print "id\tseq" } { print $1 "\t" $1 }' idt-xgen-cfdna-umis.txt > umis.tsv
head -n 3 umis.tsv
```

```tsv
id	seq
GAGACGAT	GAGACGAT
TTCCAAGG	TTCCAAGG
```

```bash
unmux synthetic_nonrandom_umi_1.fastq.gz synthetic_nonrandom_umi_2.fastq.gz \
  --group "umi1=umis.tsv" --group "umi1::loc=0:0:8,dist=1" \
  --group "umi2=umis.tsv" --group "umi2::loc=1:0:8,dist=1" \
  --extract "fixed1=@umi1" --extract "fixed2=@umi2" \
  --extract "seen1=0:0:8" --extract "seen2=1:0:8" \
  --extract "r1=0:8:end" --extract "r2=1:8:end" \
  --template r1 r2 \
  --tag "RX=fixed1+fixed2" --tag "RX::qual=none" \
  --tag "OX=seen1+seen2" \
  --out umis.unmapped.bam
```

A pair whose R1 UMI had one substitution:

```bash
samtools view umis.unmapped.bam | awk '$1 == "SRR6109255.1031"' | cut -f 1,2,12-
```

```text
SRR6109255.1031	77	RX:Z:ACGACTTG-CGCATGAT	OX:Z:ACGACTAG-CGCCTGAT	BZ:Z:IIIIIIII IIIIIIII	RG:Z:synthetic_nonrandom_umi
SRR6109255.1031	141	RX:Z:ACGACTTG-CGCATGAT	OX:Z:ACGACTAG-CGCCTGAT	BZ:Z:IIIIIIII IIIIIIII	RG:Z:synthetic_nonrandom_umi
```

On all 9,901 pairs RX holds the two list UMIs and OX the bases as read, with BZ their qualities. 989 R1 UMIs and 1,010 R2 UMIs had one substitution and were corrected. Every UMI in this data is within one substitution of the list, so no pair is left uncorrected here.

To see what happens to a UMI that cannot be corrected, allow none (`dist=0`):

<details>
<summary>The command</summary>

```bash
unmux synthetic_nonrandom_umi_1.fastq.gz synthetic_nonrandom_umi_2.fastq.gz \
  --group "umi1=umis.tsv" --group "umi1::loc=0:0:8,dist=0" \
  --group "umi2=umis.tsv" --group "umi2::loc=1:0:8,dist=0" \
  --extract "fixed1=@umi1" --extract "fixed2=@umi2" \
  --extract "seen1=0:0:8" --extract "seen2=1:0:8" \
  --extract "r1=0:8:end" --extract "r2=1:8:end" \
  --template r1 r2 \
  --tag "RX=fixed1+fixed2" --tag "RX::qual=none" \
  --tag "OX=seen1+seen2" \
  --out umis.exact.bam
```

</details>

```bash
samtools view -c -f 64 umis.exact.bam
samtools view -c -f 64 -d RX umis.exact.bam
```

```text
9901
8000
```

The 1,901 pairs with a substitution in either UMI are still written, with OX and BZ but no RX: a tag is left out when any stream in it is missing.

> [!TIP]
> To set such pairs aside instead, add `minFindsPerGroup=1` to both groups and give `--unassigned`.

> [!NOTE]
> The 32 UMIs are at least 3 substitutions apart, so `dist=1` corrects without ties. With `dist=2` unmux refuses the set, since one read could then match two UMIs equally well.

## Single cell

Cell barcodes and UMIs into SAM tags, from 10x Chromium, split-pool and antibody-tag libraries.

### 11. 10x Chromium cell barcodes and UMIs

Use this to tag 10x Chromium reads with their cell barcode and UMI in an unmapped BAM, correcting each barcode to a list of known cells: CB holds the corrected barcode, CR and CY the barcode as read with its qualities, and UR, UY and UB the UMI.

**Data:** the 100,000 10x v2 read pairs of the nf-core/scrnaseq test set (R1: 16 bp barcode and 10 bp UMI; R2: cDNA) and the 10x v2 barcode keeplist ([S10_L001_R1_001.fastq.gz](https://raw.githubusercontent.com/nf-core/test-datasets/d934d6e8367fe2626184496b1889671cf2b02dab/testdata/S10_L001_R1_001.fastq.gz), [S10_L001_R2_001.fastq.gz](https://raw.githubusercontent.com/nf-core/test-datasets/d934d6e8367fe2626184496b1889671cf2b02dab/testdata/S10_L001_R2_001.fastq.gz), [10xv2_whitelist.txt](https://raw.githubusercontent.com/nf-core/test-datasets/d934d6e8367fe2626184496b1889671cf2b02dab/reference/kallisto/10xv2_whitelist.txt)).

The keeplist holds 737,280 barcodes, and this recipe uses a subset of it: the 12,518 within one substitution of one of the 14,897 barcodes seen in R1. A read can be corrected only to such a barcode, so the subset gives every read the call the whole list would. unmux reads tag files with `id` and `seq` columns, so the barcode is both:

```bash
gzip -dc S10_L001_R1_001.fastq.gz \
  | awk 'NR % 4 == 2 { print substr($0, 1, 16) }' | sort -u > seen.txt
awk 'NR == FNR { keep[$1]; next }
  { for (i = 1; i <= 16; i++) for (j = 1; j <= 4; j++) {
      s = substr($1, 1, i - 1) substr("ACGT", j, 1) substr($1, i + 1)
      if (s in keep) print s } }' 10xv2_whitelist.txt seen.txt \
  | sort -u | awk 'BEGIN { print "id\tseq" } { print $1 "\t" $1 }' > cells.tsv
```

```bash
unmux S10_L001_R1_001.fastq.gz S10_L001_R2_001.fastq.gz \
  --pool S10 \
  --group "cb=cells.tsv" \
  --group "cb::loc=0:0:16,dist=1,mode=nearest" \
  --extract "cell=@cb" \
  --extract "barcode=0:0:16" \
  --extract "umi=0:16:26" \
  --extract "cdna=1:0:end" \
  --template cdna \
  --tag "CB=cell" \
  --tag "CR=barcode" \
  --tag "UR=umi" --tag "UR::qual=UY" \
  --tag "UB=umi" \
  --out S10.unmapped.bam
```

A read whose barcode was corrected and one with no keeplist barcode within one substitution, in file order:

```bash
samtools view S10.unmapped.bam | grep -E '^(A00111:245:HGKNTDMXX:1:1101:13241:1000|A00111:245:HGKNTDMXX:1:1101:3893:1000)[[:space:]]' | cut -f 1,2,12-
```

```text
A00111:245:HGKNTDMXX:1:1101:3893:1000	4	CR:Z:TTGCCTCTGCATATTG	CY:Z:FFFFFFFFFFFFFFFF	UR:Z:TTCCTGATTG	UY:Z:FFFFFFFFFF	UB:Z:TTCCTGATTG	RG:Z:S10
A00111:245:HGKNTDMXX:1:1101:13241:1000	4	CB:Z:CGTCTACAGGCAATTA	CY:Z:#FFFFF:FFFFFFFFF	CR:Z:NGTCTACAGGCAATTA	UR:Z:TGAAACGCTA	UY:Z:FFFFFFFFFF	UB:Z:TGAAACGCTA	RG:Z:S10
```

Of 100,000 reads, 96,512 matched a keeplist barcode exactly, 779 were corrected, and 2,709 have no CB but keep CR, UR and UB: 2,162 with no keeplist barcode within one substitution and 547 with two or more. On every read, CB agrees with a lookup of the barcode and its one-substitution neighbors in the whole keeplist, CR, CY, UR, UB and UY equal the R1 bases and qualities they come from, and the record holds R2.

unmux warns about tags one read can match equally well, which `mode=nearest` leaves uncorrected, and that no reads are split by sample, which is expected here:

```text
[WARN  unmux::grammar] groups are defined but no `--sample`, `--sample-sheet`, or `--sample-from-group` was given; reads will not be demultiplexed (they pass through to `--out` undemultiplexed). Add a sample target to name a single sample or split by sample.
[WARN  unmux::demux] group `cb`: tags `AAACCTGGTAAATGTG` and `AACCATGGTAAATGTG` can both match one read within dist=1 (the bases `AACCCTGGTAAATGTG` match both); mode=nearest leaves ambiguous reads unassigned unless delta separates them
[WARN  unmux::demux] group `cb`: tags `AAACCTGTCAACGAAA` and `AATCCAGTCAACGAAA` can both match one read within dist=1 (the bases `AATCCTGTCAACGAAA` match both); mode=nearest leaves ambiguous reads unassigned unless delta separates them
…
[WARN  unmux::demux] group `cb`: at least 990 more pairs (counting stopped at 1000) of tags can both match one read within dist=1; mode=nearest leaves ambiguous reads unassigned unless delta separates them
```

> [!NOTE]
> unmux does not correct UMIs, so UB holds the same bases as UR; tools such as velocyto read UB.

> [!WARNING]
> Without `--pool S10` the pool, and with it the read group, would be `S10_L001_R`: the common stem of the two file names.

### 12. SPLiT-seq barcodes in three rounds

Use this for split-pool data such as SPLiT-seq, where the cell barcode is three 8 bp barcodes, one per round of barcoding, separated by fixed linkers in R2. This is the README's Quick Start command.

**Data:** 10,000 read pairs of SRA run [SRR13948565](https://www.ncbi.nlm.nih.gov/sra/SRR13948565), the C2C12_short_1k library of LR-Split-seq (Rebboah et al. 2021) made with the Parse Biosciences v1 kit, whose R2 has the chemistry v2 layout of the split-seq-pipeline: a 10 bp UMI, round 3, linker 1, round 2, linker 2 and round 1 (`NNNNNNNNNNIIIIIIIIGTGGCCGATGTTTCGCATCGGCGTACGACTIIIIIIIIATCCACGTGCTTGAGACTGTGGIIIIIIII`, where N is the UMI and I a barcode). The barcodes are the pipeline's published lists ([bc_8nt_v2.csv](https://raw.githubusercontent.com/yjzhang/split-seq-pipeline/c3923ea44a94a1fc0effe5ff8e2ca479c0bff934/split_seq/barcodes/bc_8nt_v2.csv), [bc_8nt_v1.csv](https://raw.githubusercontent.com/yjzhang/split-seq-pipeline/c3923ea44a94a1fc0effe5ff8e2ca479c0bff934/split_seq/barcodes/bc_8nt_v1.csv), [tools.py](https://raw.githubusercontent.com/yjzhang/split-seq-pipeline/c3923ea44a94a1fc0effe5ff8e2ca479c0bff934/split_seq/tools.py)).

```bash
fg-sra fastq --min-spot-id 1000001 --max-spot-id 1010000 --r1 SRR13948565_1.fastq --r2 SRR13948565_2.fastq SRR13948565
```

```text
@SRR13948565.1000001
TAATATCTTACTGGCATAGTGGCCGATGTTTCGCATCGGCGTACGACTATCATTCCATCCACGTGCTTGAGACTGTGGCGTGGTTG
+
66AAAEEEEEEEEEEEEEEAAE<AEAAAAEEE/EEEEEEEAEE6EEEEEEAEEEA<<EEEEAEE/AA/EEEEEEAEE/EEEEEEEA
```

Make the two tag files from the published lists:

```bash
awk -F , 'BEGIN { print "id\tseq" } { print "bc" $1 "\t" $2 }' bc_8nt_v2.csv > round1.tags.txt
awk -F , 'BEGIN { print "id\tseq" } { print "bc" $1 "\t" $2 }' bc_8nt_v1.csv > round2_3.tags.txt
head -n 3 round1.tags.txt
```

```tsv
id	seq
bc0	ACTCGTAA
bc1	AAACGATA
```

```bash
unmux R1.fastq.gz R2.fastq.gz \
  --group "round1=round1.tags.txt" \
  --group "round1::loc=1:78:86" \
  --group "round1::dist=1" \
  --group "round1::minFindsPerGroup=1" \
  --group "round2=round2_3.tags.txt" \
  --group "round2::loc=1:48:56" \
  --group "round2::dist=1" \
  --group "round2::minFindsPerGroup=1" \
  --group "round3=round2_3.tags.txt" \
  --group "round3::loc=1:10:18" \
  --group "round3::dist=1" \
  --group "round3::minFindsPerGroup=1" \
  --extract "cdna=0:0:end" \
  --extract "umi=1:0:10" \
  --extract "bc1=@round1" \
  --extract "bc2=@round2" \
  --extract "bc3=@round3" \
  --template "cdna" \
  --tag "CB=bc1+bc2+bc3" \
  --tag "CB::qual=CY" \
  --tag "CB::sep=-" \
  --tag "CB::qual-sep=%20" \
  --tag "RX=umi" \
  --tag "RX::qual=QX" \
  --out "splitseq.unmapped.bam"
```

```bash
samtools view splitseq.unmapped.bam | head -n 2 | cut -f 1,2,12-
```

```text
SRR13948565.1000001	4	CB:Z:CGTGGTTG-ATCATTCC-CTGGCATA	CY:Z:EEEEEEEA EEAEEEA< EEEEEEEE	RX:Z:TAATATCTTA	QX:Z:66AAAEEEEE	RG:Z:R
SRR13948565.1000002	4	CB:Z:ATCATGCA-GAGCTGAA-CGGATTGC	CY:Z:EEEEEEEE EEEEEEEE EEEEEEEE	RX:Z:GCTTGTACTG	QX:Z:AAAAAEEEEE	RG:Z:R
```

CB joins the three corrected barcodes with `-`, and CY their qualities with spaces; RX and QX hold the UMI. A read that misses any round fails `minFindsPerGroup=1` and is dropped.

7,116 of 10,000 reads were written. The generator looks up each read's three windows itself, a published barcode within one substitution in each, and every written read has the CB, CY, RX and QX that lookup gives, while every read it drops misses a round; each record holds its R1. 5,453 of the written reads belong to 546 of the 568 cells that split-pipe called in this library ([GSM5169184](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5169nnn/GSM5169184/suppl/GSM5169184_C2C12_short_1k_cell_metadata.csv.gz)).

```text
[WARN  unmux::grammar] groups are defined but no `--sample`, `--sample-sheet`, or `--sample-from-group` was given; reads will not be demultiplexed (they pass through to `--out` undemultiplexed). Add a sample target to name a single sample or split by sample.
[WARN  unmux::metrics] metrics pool=R: 28.8% of reads were unassigned (>= 20%)
```

The warning is expected: this command tags reads and does not split them by sample.

> [!NOTE]
> The Parse kit primes each round-1 well with an oligo-dT and a random-hexamer barcode, so one cell's reads carry either, and this command gives the cell two CBs; 404 called cells appear under both here. [recipe 13](#13-keep-the-reads-of-called-cells-with-a-joined-match) joins them through a cell list.

> [!NOTE]
> Without `--pool`, the pool is the common stem of the input names, here `R` from R1 and R2, so every record gets `RG:Z:R`; add `--pool` to set it.

**Compared with splitcode 0.31.6.** splitcode's SPLiT-seq tutorial ([tutorials_splitseq.rst](https://raw.githubusercontent.com/pachterlab/splitcode-tutorial/68387770995199a9ac4ccadc25ee2f8e96cd24be/docs/source/tutorials_splitseq.rst)) uses the same layout as this command: rounds 1, 2 and 3 at R2 bases 78 to 86, 48 to 56 and 10 to 18, one substitution, and the UMI in R2's first 10 bases. Its config.txt for barcode reformatting, with its barcode files ([r1_R.txt](https://raw.githubusercontent.com/pachterlab/splitcode-tutorial/68387770995199a9ac4ccadc25ee2f8e96cd24be/uploads/splitseq/r1_R.txt), [r1_T.txt](https://raw.githubusercontent.com/pachterlab/splitcode-tutorial/68387770995199a9ac4ccadc25ee2f8e96cd24be/uploads/splitseq/r1_T.txt), [r2_r3.txt](https://raw.githubusercontent.com/pachterlab/splitcode-tutorial/68387770995199a9ac4ccadc25ee2f8e96cd24be/uploads/splitseq/r2_r3.txt)):

```text
@extract <output_R2{*}>,1:0<umi[10]>
tags         distances  ids     groups    minFindsG	locations
r1_R.txt$    1          r1_R    round1    1          1,78,86
r1_T.txt$    1          r1_T    round1    1          1,78,86
r2_r3.txt$   1          r2      round2_3  2          1,48,56
r2_r3.txt$   1          r3      round2_3  2          1,10,18
```

```bash
splitcode -c config.txt --nFastqs=2 --select=0 \
--gzip -o output_R1.fastq.gz,unused_R2.fastq.gz R1.fastq.gz R2.fastq.gz
```

As the tutorial gives the command, with one `-o` file, splitcode 0.31.6 stops with `Error: Incorrect number of output files`, since it counts the `-o` files against `--nFastqs` even with `--select=0`. A second name, whose file stays empty, is the only change. The file output_R2.fastq.gz then holds each kept read's corrected rounds 3, 2 and 1.

The tutorial's round-1 files hold 192 barcodes, of which 43 are among the 96 of this library's list, while its round-2 and round-3 file is this library's list. A second run keeps the configuration and puts this library's oligo-dT and random-hexamer barcodes in r1_T.txt and r1_R.txt:

```bash
awk -F , '$1 < 48 { print $2 }' bc_8nt_v2.csv > r1_T.txt
awk -F , '$1 >= 48 { print $2 }' bc_8nt_v2.csv > r1_R.txt
```

Per read, against unmux:

| Reads | tutorial files | this library's lists |
| --- | ---: | ---: |
| Same call from both | 3,473 | 7,116 |
| Different calls | 6 | 0 |
| A call from unmux only | 3,637 | 0 |
| A call from splitcode only | 211 | 0 |
| No call from either | 2,673 | 2,884 |

With this library's lists, splitcode makes unmux's call on every read: 7,116 with a CB and 2,884 with none. With the tutorial's files, of the 3,637 reads only unmux calls, 3,628 have a round-1 barcode outside the tutorial's lists and 9 have round-1 bases within one substitution of two of its barcodes; the 211 only splitcode calls hold one of the tutorial's round-1 barcodes, 174 of them exactly, more than one substitution from any of this library's; and the 6 calls that differ agree on rounds 2 and 3, while round 1 is one substitution from this library's barcode and within one of a tutorial barcode, which splitcode takes.

How unmux expresses splitcode's settings:

| splitcode | unmux |
| --- | --- |
| tag files such as `r1_R.txt$`, one tag per line, in group `round1` | `--group "round1=round1.tags.txt"`, one tag file for both primers |
| `distances` 1 | `dist=1` |
| `locations` 1,78,86 | `loc=1:78:86` |
| `minFindsG` 1 for round 1, and 2 for rounds 2 and 3 together | `minFindsPerGroup=1` on each round |
| `@extract 1:0<umi[10]>` | `--extract "umi=1:0:10"` and `--tag "RX=umi"` |
| `@extract <output_R2{*}>`, the corrected tags stitched together | `--tag "CB=bc1+bc2+bc3"` |
| `--select=0`, writing R1 only | `--template "cdna"` |
| `--min-delta` | `mode=nearest`, with `delta` |
| `--keep-grp`, files of tag-group arrangements to keep | `--sample-sheet` |
| `--subs`, such as to turn R barcodes into T barcodes | a `match=` group with a column, as in [recipe 13](#13-keep-the-reads-of-called-cells-with-a-joined-match) |
| `--revcomp` | `bothStrands=true` |
| `next` and `previous`, with the bases between | `next=` |
| `--partial5` and `--partial3` | `partial5` and `partial3` |
| `--qtrim` and `--filter-len` | no equivalent |
| `--assign` and `--mapping`, an ID for each combination of tags | no equivalent |

What splitcode does that unmux does not: it trims bases by quality and filters reads by length, and it can give each combination of tags a short ID with a mapping file. It also nests configurations and remultiplexes reads.

### 13. Keep the reads of called cells with a joined match

Use this to keep only reads from cells you have called and write each cell's name: a `match=` group looks up the three corrected round barcodes, joined, in a list of cells, and a column of the list gives the name.

**Data:** the SPLiT-seq reads of [recipe 12](#12-split-seq-barcodes-in-three-rounds) and the 568 cells split-pipe called in that library ([GSM5169184](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5169nnn/GSM5169184/suppl/GSM5169184_C2C12_short_1k_cell_metadata.csv.gz)). Each cell's `cell_barcode` is its round-3 and round-2 barcodes, then its round-1 well, and the well's oligo-dT and random-hexamer barcodes are rows w and w + 48 of the round-1 list, so the list gets a row for each, with the same `cell`:

```bash
gzip -dc GSM5169184_C2C12_short_1k_cell_metadata.csv.gz | awk -F , 'NR == FNR { rt[$1] = $2; next }
  FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; print "id\tseq\tcell"; next }
  { cb = $c["cell_barcode"]; w = $c["rnd1_well"]; r3 = substr(cb, 1, 8); r2 = substr(cb, 9, 8)
    print cb "_dT\t" rt[w] r2 r3 "\t" cb; print cb "_hex\t" rt[w + 48] r2 r3 "\t" cb }' bc_8nt_v2.csv - > cells.tsv
head -n 3 cells.tsv
```

```tsv
id	seq	cell
CAAGGAGCATAGCGAC_22_dT	TACGGTTAATAGCGACCAAGGAGC	CAAGGAGCATAGCGAC_22
CAAGGAGCATAGCGAC_22_hex	CATAGCTAATAGCGACCAAGGAGC	CAAGGAGCATAGCGAC_22
```

```bash
unmux R1.fastq.gz R2.fastq.gz \
  --pool C2C12_short_1k \
  --group "round1=round1.tags.txt" --group "round1::loc=1:78:86,dist=1" \
  --group "round2=round2_3.tags.txt" --group "round2::loc=1:48:56,dist=1" \
  --group "round3=round2_3.tags.txt" --group "round3::loc=1:10:18,dist=1" \
  --group "cells=cells.tsv" \
  --group "cells::match=bc1+bc2+bc3,minFindsPerGroup=1" \
  --extract "bc1=@round1" --extract "bc2=@round2" --extract "bc3=@round3" \
  --extract "name=@cells::cell" \
  --extract "cdna=0:0:end" --extract "umi=1:0:10" \
  --template cdna \
  --tag "CB=name" \
  --tag "RX=umi" \
  --out cells.unmapped.bam \
  --unassigned other.unmapped.bam
```

```bash
samtools view cells.unmapped.bam | head -n 2 | cut -f 1,2,12-
```

```text
SRR13948565.1000001	4	CB:Z:CTGGCATAATCATTCC_42	RX:Z:TAATATCTTA	QX:Z:66AAAEEEEE	RG:Z:C2C12_short_1k
SRR13948565.1000007	4	CB:Z:CACTTCGAATCCTGTA_12	RX:Z:GTACGCTCTT	QX:Z:AAAAAEEEEE	RG:Z:C2C12_short_1k
```

5,453 reads were kept, of 546 cells, and 4,547 went to other.unmapped.bam, both records of each pair as they were read. Checked against the generator's own lookup of the three windows: the kept reads are exactly those whose joined barcode is a called cell's, each with that cell's name, both primers of a cell under one name. 381 of them needed a correction in some round and are kept because `match=` joins each round as corrected, the way `--tag` writes it.

unmux warns, as in [recipe 12](#12-split-seq-barcodes-in-three-rounds), that nothing is split by sample; that CB, a column of cells.tsv, has no base qualities to write to CY, which `CB::qual=none` silences; and that 45.5% of the reads were unassigned, the 4,547 in other.unmapped.bam:

```text
[WARN  unmux::grammar] groups are defined but no `--sample`, `--sample-sheet`, or `--sample-from-group` was given; reads will not be demultiplexed (they pass through to `--out` undemultiplexed). Add a sample target to name a single sample or split by sample.
[WARN  unmux::grammar] tag `CB` has a column stream, which has no base qualities, so its default quality tag `CY` is left out; set `CB::qual=none` to silence this.
[WARN  unmux::metrics] metrics pool=C2C12_short_1k: 45.5% of reads were unassigned (>= 20%)
```

> [!TIP]
> Without `minFindsPerGroup=1` on the `match=` group, reads of other barcodes are written too, without a CB.

### 14. CITE-seq antibody barcodes

Use this for antibody-derived tags (ADT) in CITE-seq or 10x feature barcoding: find the antibody barcode in R2, correct it, and write the antibody's name next to the cell barcode and UMI.

**Data:** the 10,000 antibody read pairs of the 10x 5k anti-CMV T cell set in nf-core/test-datasets, with their feature reference of 11 TotalSeq-C antibodies ([subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz](https://raw.githubusercontent.com/nf-core/test-datasets/e8d2f4f2af42cc757ed2448c2398e282e9990235/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/ab/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz), [subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz](https://raw.githubusercontent.com/nf-core/test-datasets/e8d2f4f2af42cc757ed2448c2398e282e9990235/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/ab/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz), [5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv](https://raw.githubusercontent.com/nf-core/test-datasets/e8d2f4f2af42cc757ed2448c2398e282e9990235/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv)). R1 is a 16 bp cell barcode and the UMI; in R2 the 15 bp antibody barcode follows 10 other bases, as the reference's pattern `^NNNNNNNNNN(BC)NNNNNNNNN` gives:

```bash
gzip -dc subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz | head -n 2
```

```text
@A00519:1254:HLVWHDSX3:4:1101:3296:1000 2:N:0:AAGGGTTTAC+NGCGTGAGTA
TTCATACTTGTCGCAACTTGATGATGGCAGAGTGAAGGTATAAGTTAGCATAGCGAACAAGCGCACCTTATTAAGATCTGAAGAGAGTCG
```

Give the feature reference the `id` and `seq` columns unmux reads:

```bash
awk -F , 'NR == 1 { print "id\tseq"; next } { print $1 "\t" $5 }' \
  5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv > antibodies.tsv
head -n 3 5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv
head -n 3 antibodies.tsv
```

```tsv
id,name,read,pattern,sequence,feature_type
CD3,CD3,R2,^NNNNNNNNNN(BC)NNNNNNNNN,CTCATTGTAACTCCT,Antibody Capture
CD4,CD4,R2,^NNNNNNNNNN(BC)NNNNNNNNN,TGTTCCCGCTCAACT,Antibody Capture
id	seq
CD3	CTCATTGTAACTCCT
CD4	TGTTCCCGCTCAACT
```

```bash
unmux subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz \
  subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz \
  --pool antiCMV_AB \
  --group "ab=antibodies.tsv" \
  --group "ab::loc=1:10:25,dist=1,minFindsPerGroup=1" \
  --extract "cell=0:0:16" \
  --extract "umi=0:16:end" \
  --extract "barcode=@ab" \
  --extract "antibody=@ab::id" \
  --extract "r2=1:0:end" \
  --template r2 \
  --tag "CR=cell" \
  --tag "UR=umi" --tag "UR::qual=UY" \
  --tag "fb=barcode" \
  --tag "fx=antibody" \
  --out adt.unmapped.bam
```

```bash
samtools view adt.unmapped.bam | head -n 2 | cut -f 1,2,12-
```

```text
A00519:1254:HLVWHDSX3:4:1101:3296:1000	4	CR:Z:NAATAAGGTGCGCTTG	CY:Z:#:FFFFFF:FFFFFFF	UR:Z:TTCGCTATGCTA	UY:Z:,FFFFFFFFF,,	fb:Z:GCGCAACTTGATGAT	fx:Z:CD8	RG:Z:antiCMV_AB
A00519:1254:HLVWHDSX3:4:1101:4182:1000	4	CR:Z:NGCGGTACATACTACG	CY:Z:#FF:FFFFFFFFFFF,	UR:Z:GATTATCAATTT	UY:Z:FFFFFFFFFFF,	fb:Z:TCTCAGACCTCCGTA	fx:Z:CD14	RG:Z:antiCMV_AB
```

The `@ab::id` stream is the matched antibody's `id` column, written as text to fx; fb holds the corrected barcode.

```bash
samtools view adt.unmapped.bam | grep -o 'fx:Z:[^[:space:]]*' \
  | sort | uniq -c | awk '{ print $1, $2 }' | sort -k1,1nr -k2,2 | head -n 5
```

```text
2591 fx:Z:CD3
1946 fx:Z:CD8
1524 fx:Z:CD14
646 fx:Z:CD19
620 fx:Z:CD4
```

9,460 of 10,000 reads were written: every read whose R2 bases 10 to 25 hold a panel barcode within one substitution (9,187 of them exactly), checked against the generator's own lookup. On every record fx and fb match that lookup, CR, CY, UR and UY equal the R1 bases and qualities, and the record holds R2. The panel's barcodes are at least 7 substitutions apart, so one substitution cannot reach a wrong antibody.

> [!NOTE]
> The tags fb and fx follow Cell Ranger's names for a feature barcode and its feature. Any two-letter tag works; lowercase ones are free for local use.

## Long reads

Nanopore and PacBio reads: barcodes at both ends, U bases, and tags carried in FASTQ comments.

### 15. Nanopore native barcodes at both ends

Use this to demultiplex Nanopore reads from a native barcoding kit, where each read starts with the barcode and ends with its reverse complement, trimming the barcodes off.

**Data:** 200 reads of each of 24 SRA runs from one MinION flow cell of SARS-CoV-2 ARTIC amplicons ([PRJNA730334](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA730334), batch VERD6, kits SQK-LSK109 with EXP-NBD104 and EXP-NBD114). Guppy split the reads into one run per native barcode, requiring the barcode at both ends and leaving it in the reads, so the barcode of every read is known. The 4,800 reads are pooled one run after another:

```bash
for run in SRR14560737 SRR14560738 SRR14560739 SRR14560652 SRR14560583 SRR14560584 \
  SRR14560585 SRR14560586 SRR14560587 SRR14560588 SRR14560763 SRR14560764 \
  SRR14560765 SRR14560523 SRR14560525 SRR14560631 SRR14560632 SRR14560633 \
  SRR14560634 SRR14560635 SRR14560769 SRR14560772 SRR14560526 SRR14560531; do
  fg-sra fastq --min-spot-id 1 --max-spot-id 200 --unpaired "$run.fastq" "$run"
done
```

Each tag is a native barcode between its flanks, from dorado ([barcode_kits.cpp](https://raw.githubusercontent.com/nanoporetech/dorado/651980f48c21eb93a5c46c65e1f004136f02c07a/dorado/utils/barcode_kits.cpp)), and the sheet asks for the same barcode at both ends:

```bash
head -n 3 native_barcodes.tsv
head -n 3 barcode_pairs.tsv
```

```tsv
id	seq
NB01	ATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT
NB02	ATTGCTAAGGTTAAACAGACGACTACAAACGGAATCGACAGCACCT
sample	front	rear
NB01	NB01	NB01
NB02	NB02	NB02
```

```bash
unmux reads.fastq.gz \
  --group "front=native_barcodes.tsv" \
  --group "front::loc=0:0:150,dist=12:12:12,mode=nearest" \
  --group "rear=native_barcodes.tsv" \
  --group "rear::loc=0:-150:end,dist=12:12:12,mode=nearest,bothStrands=true" \
  --extract "insert=@front..@rear" \
  --template insert \
  --sample-sheet barcode_pairs.tsv \
  --out "demux/%sample.fastq.gz" \
  --unassigned demux/unclassified.fastq.gz
```

The `front` group searches the first 150 bases and the `rear` group the last 150 on both strands, so it finds the reverse-complemented copy. A budget of 12 edits in a 46 bp tag suits the error rate of these R9.4.1 reads, and `mode=nearest` leaves a read that two barcodes fit equally well unclassified. The span `@front..@rear` keeps the bases between the two tags.

```bash
for f in demux/NB0[1-3].fastq.gz demux/unclassified.fastq.gz; do
  echo "$f $(gzip -dc "$f" | wc -l | awk '{print $1 / 4}')"
done
```

```text
demux/NB01.fastq.gz 199
demux/NB02.fastq.gz 199
demux/NB03.fastq.gz 200
demux/unclassified.fastq.gz 52
```

With a budget this large, unmux warns about pairs of tags that one read could match, listing 10 for each group and counting the rest:

```text
[WARN  unmux::demux] group `front`: tags `ATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT` and `ATTGCTAAGGTTAAACAGACGACTACAAACGGAATCGACAGCACCT` can both match one read within dist=12:12:12 (the bases `ATTGCTAAGGTTAAACAGACGACTACAACAACTTTCTTCAGCACCT` match both); mode=nearest leaves ambiguous reads unassigned unless delta separates them
…
[WARN  unmux::demux] group `front`: 266 more pairs of tags can both match one read within dist=12:12:12; mode=nearest leaves ambiguous reads unassigned unless delta separates them
[WARN  unmux::demux] group `rear`: tags `ATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT` and `ATTGCTAAGGTTAAACAGACGACTACAAACGGAATCGACAGCACCT` can both match one read within dist=12:12:12 (the bases `ATTGCTAAGGTTAAACAGACGACTACAACAACTTTCTTCAGCACCT` match both); mode=nearest leaves ambiguous reads unassigned unless delta separates them
…
[WARN  unmux::demux] group `rear`: 266 more pairs of tags can both match one read within dist=12:12:12; mode=nearest leaves ambiguous reads unassigned unless delta separates them
```

The windows are much wider than the tags, as long reads need, and unmux prints no window warning: it warns about one only where it estimates that over 0.3% of reads with an error could match another tag at a shifted position.

4,748 of 4,800 reads (98.9%) went to the barcode of their run and none to another; 52 are unclassified. Every trimmed read is part of its read, cut within the search windows: a median of 68 bases off the front and 54 off the rear.

> [!TIP]
> With `dist=8:8:8`, 4,351 reads go to the barcode of their run, none to another, and 449 are unclassified: a tighter budget costs reads here without saving any from a wrong barcode.

What a trim does to the basecaller's tags, on the ten modkit reads of [recipe 18](#18-lift-sam-tags-from-fastq-comments-back-onto-records), before and after a 20-base trim:

```bash
unmux modbam.unmapped.bam --out kept.bam
unmux modbam.unmapped.bam --extract "insert=0:20:end" --template insert --out trimmed.bam
samtools view kept.bam | head -n 1 | cut -f 12- | tr '\t' '\n' | cut -d : -f 1 | paste -s -d ' ' -
samtools view trimmed.bam | head -n 1 | cut -f 12- | tr '\t' '\n' | cut -d : -f 1 | paste -s -d ' ' -
```

```text
RG qs ns ts mx ch st du rn f5 sm sd sv mv MM ML
RG qs ns ts mx ch st du rn f5 sm sd sv
```

```text
[WARN  unmux::demux] --template left out per-base tags from 10 records that did not copy a read exactly: MM, ML, mv
```

> [!WARNING]
> A `--template` record keeps the tags of its read, but per-base tags such as MM, ML and the move table mv describe the read's bases one by one, so unmux leaves them out of a record that does not copy its read exactly, with one warning for the run. Demultiplex without a template to keep the modification calls; the barcodes then stay in the reads.

### 16. PacBio barcode pairs, compared with lima

Use this to demultiplex PacBio HiFi reads that still carry their barcodes, one at each end, as lima does, here with a different barcode at each end of a sample's amplicons.

**Data:** the first 1,922 HiFi reads of PacBio's CYP2D6 amplicon set, 22 samples each with a pair of barcodes, before lima ran on them, cut from the first 10 MB of the BAM by an HTTP range request ([prefix.bam](https://downloads.pacbcloud.com/public/dataset/pbAmpliconAnalysis_CYP2D6/workflow_inputs/m64013e_210223_231140.hifi_reads.bam), [biosamples.csv](https://downloads.pacbcloud.com/public/dataset/pbAmpliconAnalysis_CYP2D6/workflow_inputs/biosamples.csv), [Sequel_96_barcodes_v1.fasta](https://raw.githubusercontent.com/PacificBiosciences/hifi-amplicon-workflow/2964cd1fbddae04243d56914935881cca01f6d13/data/Sequel_96_barcodes_v1.fasta)). Each read starts with a pad of a few bases, a barcode and the M13 forward or reverse primer, and ends with the same for the other barcode, reverse-complemented.

```bash
curl -sS -A unmux-recipes/1.0 -r 0-9999999 -o prefix.bam https://downloads.pacbcloud.com/public/dataset/pbAmpliconAnalysis_CYP2D6/workflow_inputs/m64013e_210223_231140.hifi_reads.bam
```

Keep the complete records of the cut, and build a tag of each barcode with its M13 primer and a sheet with each pair in both orders:

```bash
samtools view -h prefix.bam 2> /dev/null | samtools view -b -o hifi_reads.bam -
awk -F , 'NR == FNR { if (/^>/) id = substr($0, 2); else seq[id] = $0; next }
  FNR == 1 { print "id\tseq" > "barcodes.tsv"; print "sample\tstart\tend" > "barcode_pairs.tsv"; next }
  { split($1, bc, "--")
    if (!(bc[1] in seen)) { seen[bc[1]]; print bc[1] "\t" seq[bc[1]] "GTAAAACGACGGCCAGT" > "barcodes.tsv" }
    if (!(bc[2] in seen)) { seen[bc[2]]; print bc[2] "\t" seq[bc[2]] "CAGGAAACAGCTATGAC" > "barcodes.tsv" }
    print $2 "\t" bc[1] "\t" bc[2] > "barcode_pairs.tsv"; print $2 "\t" bc[2] "\t" bc[1] > "barcode_pairs.tsv" }' \
  Sequel_96_barcodes_v1.fasta biosamples.csv
```

```bash
head -n 3 barcodes.tsv
head -n 3 barcode_pairs.tsv
```

```tsv
id	seq
bc1016	CATAGAGAGATAGTATGTAAAACGACGGCCAGT
bc1056	AGACACACACGCACATCAGGAAACAGCTATGAC
sample	start	end
NA02016	bc1016	bc1056
NA02016	bc1056	bc1016
```

lima 26.2.1, with its HiFi preset for a different barcode at each end:

```bash
lima --hifi-preset ASYMMETRIC --split-named -j 2 --biosample-csv biosamples.csv \
  hifi_reads.bam Sequel_96_barcodes_v1.fasta lima.bam
```

and unmux:

```bash
unmux hifi_reads.bam \
  --group "start=barcodes.tsv" \
  --group "start::loc=0:0:60,dist=2:2:3,mode=nearest" \
  --group "end=barcodes.tsv" \
  --group "end::loc=0:-60:end,dist=2:2:3,mode=nearest,bothStrands=true" \
  --extract "insert=@start..@end" \
  --template insert \
  --sample-sheet barcode_pairs.tsv \
  --out "unmux/%sample.bam" \
  --unassigned unmux/unassigned.bam
```

unmux warns about pairs of tags that one read could match within this budget, listing 10 for each group, and about the share of reads it leaves unassigned:

```text
[WARN  unmux::demux] group `start`: tags `CATAGAGAGATAGTATGTAAAACGACGGCCAGT` and `CGCATGACACGTGTGTGTAAAACGACGGCCAGT` can both match one read within dist=2:2:3 (the bases `CGCATAGAGAGTGTGTGTAAAACGACGGCCAGT` match both); mode=nearest leaves ambiguous reads unassigned unless delta separates them
…
[WARN  unmux::demux] group `start`: 51 more pairs of tags can both match one read within dist=2:2:3; mode=nearest leaves ambiguous reads unassigned unless delta separates them
[WARN  unmux::demux] group `end`: tags `CATAGAGAGATAGTATGTAAAACGACGGCCAGT` and `CGCATGACACGTGTGTGTAAAACGACGGCCAGT` can both match one read within dist=2:2:3 (the bases `CGCATAGAGAGTGTGTGTAAAACGACGGCCAGT` match both); mode=nearest leaves ambiguous reads unassigned unless delta separates them
…
[WARN  unmux::demux] group `end`: 51 more pairs of tags can both match one read within dist=2:2:3; mode=nearest leaves ambiguous reads unassigned unless delta separates them
[WARN  unmux::metrics] metrics pool=hifi_reads: 34.2% of reads were unassigned (>= 20%)
```

Per read, against lima's report (1,261 reads with a sample):

| Call | Reads |
| --- | ---: |
| Same sample from lima and unmux | 1,253 |
| A sample from lima only | 8 |
| A sample from unmux only | 11 |
| No sample from either | 650 |

The two never call different samples for one read. lima scores each barcode by alignment, and its ASYMMETRIC preset removes a read whose combined score is 80 or less: the 11 reads only unmux assigns score from 53 to 78, though both of their ends fit one sheet row within unmux's budget. unmux needs a tag within its budget at both ends, and in the 8 reads only lima assigns, its `--qc-tag` found a tag at one end only in all 8.

Of the 650 reads neither assigns, unmux found a tag at one end only in 504 and at neither end in 134, and lima passed fewer than two barcode regions in 623. The other 12 carry tags at both ends but not a sample's pair: the same barcode twice in 10 and a pair no sample has in 2, and lima reports the same barcodes. One read, of 47 bases, is missing from lima's report.

Every trimmed read is part of its read, cut within the search windows: a median of 38 bases off the start and 38 off the end, each a 33-base tag and its pad.

How unmux expresses lima's settings:

| Setting | lima | unmux |
| --- | --- | --- |
| The same barcode at both ends | `--hifi-preset SYMMETRIC`, which sets `--same` | a sheet row with that barcode in both groups, as in [recipe 15](#15-nanopore-native-barcodes-at-both-ends) |
| A different barcode at each end | `--hifi-preset ASYMMETRIC`, which sets `--different` | a sheet row for the pair in each order, since a read can come from either strand |
| The reverse-complemented barcode at the far end | searched by default | `bothStrands=true` on the end group |
| Where to search | `-W` in bases, or `-w` as a multiple of the barcode length | `loc=0:0:N` for the start and `loc=0:-N:end` for the end |
| How well a barcode must fit | `--min-score` and `--min-end-score`, alignment scores out of 100 | `dist`, a budget of substitutions, indels and edits in all |
| How far the best must lead | `--min-score-lead` | `mode=nearest`, with `delta` for a larger lead |
| Sample names and files | `--biosample-csv`, and `--split-named` for a file per barcode pair | `--sample-sheet`, and `--out "%sample.bam"` for a file per sample |
| Barcode removal | clipped by default, kept with `--no-clip` | cut with `--template` and `@start..@end`, kept with no template |
| Read filters | `--min-ref-span`, `--min-scoring-regions` and `--min-length` | no equivalent |
| Finding which barcodes are in use | `--peek-guess` | no equivalent: the sheet lists them |
| Per-read reports | `bc` and `bq` tags, and the `.lima.report` file | `--qc-tag` JSON, and the metrics files |

> [!NOTE]
> Both tools keep the tags PacBio wrote (ec, np, rq, sn, we, ws, zm) on the reads they clip, and lima adds its own. unmux sets RG to the sample's read group.

### 17. Direct RNA reads with U bases to BAM

Use this to store direct RNA reads, which basecallers can write with U for uracil, in an unmapped BAM.

**Data:** a dorado test read written with U, with SAM tags in its FASTQ comment ([fastq_with_us_and_tags.fq](https://raw.githubusercontent.com/nanoporetech/dorado/651980f48c21eb93a5c46c65e1f004136f02c07a/tests/data/bam_reader/fastq_with_us_and_tags.fq)).

```bash
head -n 2 direct_rna.fastq | cut -c 1-100
```

```text
@8623ac42-0956-4692-9a93-bdd99bf1a94e	st:Z:2023-06-22T07:17:48.308+00:00	RG:Z:6a94c5e38fbe36232d63fd
GUUAUGUUAGUCUACCUUGGUUGCAAGUUGUGGUGCUUGCUGGUUUUUGGCUUGAAGUCUGUGUUUGAAAGCCCUGGCCUUCACAACCCAAACUGAAACU
```

```bash
unmux direct_rna.fastq --out direct_rna.bam
```

```bash
samtools view direct_rna.bam | cut -f 10 | cut -c 1-60
samtools view direct_rna.bam | cut -f 1,2,12-
```

```text
GTTATGTTAGTCTACCTTGGTTGCAAGTTGTGGTGCTTGCTGGTTTTTGGCTTGAAGTCT
8623ac42-0956-4692-9a93-bdd99bf1a94e	4	st:Z:2023-06-22T07:17:48.308+00:00	RG:Z:direct_rna
```

BAM has no code for U, so all 1,473 U bases are stored as T, as `samtools import` stores them. The tags in the comment become SAM tags (more in [recipe 18](#18-lift-sam-tags-from-fastq-comments-back-onto-records)), and RG:Z is replaced by the pool, as for any input.

SAM and FASTQ output keep the U as read:

```bash
unmux direct_rna.fastq --out direct_rna.sam
grep -v '^@' direct_rna.sam | cut -f 10 | cut -c 1-40
unmux direct_rna.fastq --out copy.fastq
head -n 2 copy.fastq | tail -n 1 | cut -c 1-40
```

```text
GUUAUGUUAGUCUACCUUGGUUGCAAGUUGUGGUGCUUGC
GUUAUGUUAGUCUACCUUGGUUGCAAGUUGUGGUGCUUGC
```

> [!TIP]
> When matching, U in a read or a tag counts as T, so DNA barcodes find RNA reads written with U, on either strand. The generator checks this with a 24-mer of this read's reverse complement written with T.

### 18. Lift SAM tags from FASTQ comments back onto records

Use this when reads went through FASTQ with their SAM tags in the name comment, as `samtools fastq -T` and basecallers write them, and the tags, MM and ML among them, should be SAM tags again.

**Data:** ten unmapped Nanopore reads with modified-base calls ([bc_anchored_10_reads.unmapped.bam](https://raw.githubusercontent.com/nanoporetech/modkit/5cecc3fb3a9336068d9e3c68d5c08d678153dd2c/tests/resources/bc_anchored_10_reads.unmapped.bam)), written to FASTQ with every tag in the comment:

```bash
samtools fastq -T '*' modbam.unmapped.bam > modbam.fastq
head -n 1 modbam.fastq | tr '\t' '\n' | cut -c 1-60
```

```text
@068ce426-129e-4870-bd34-16cd78edaa43
NM:i:10
MD:Z:16^T45^C13^A0C9G1C2G2G12^G4
RG:Z:5598049b1b3264566b162bf035344e7ec610d608_dna_r10.4.1_e8
qs:i:10
ns:i:1542
ts:i:0
mx:i:4
ch:i:442
st:Z:2022-08-18T14:16:10.445+00:00
du:f:0.3855
rn:i:6335
f5:Z:merged-merged.pod5
sm:f:94.241
sd:f:26.357
sv:Z:quantile
mv:B:c,5,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0
MM:Z:C+h?,4,1,11,0,0,2,0,2;C+m?,4,1,11,0,0,2,0,2;
ML:B:C,78,34,67,52,38,53,81,15,27,17,188,203,217,112,163,239
```

```bash
unmux modbam.fastq --out modbam.roundtrip.bam
```

```bash
samtools view modbam.roundtrip.bam | head -n 1 | cut -f 1,2,12- | tr '\t' '\n' | cut -c 1-60
```

```text
068ce426-129e-4870-bd34-16cd78edaa43
4
NM:i:10
MD:Z:16^T45^C13^A0C9G1C2G2G12^G4
RG:Z:modbam
qs:i:10
ns:i:1542
ts:i:0
mx:i:4
ch:i:442
st:Z:2022-08-18T14:16:10.445+00:00
du:f:0.3855
rn:i:6335
f5:Z:merged-merged.pod5
sm:f:94.241
sd:f:26.357
sv:Z:quantile
mv:B:c,5,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0
MM:Z:C+h?,4,1,11,0,0,2,0,2;C+m?,4,1,11,0,0,2,0,2;
ML:B:C,78,34,67,52,38,53,81,15,27,17,188,203,217,112,163,239
```

On all 10 reads the bases, qualities and every tag, MM and ML included, equal the original BAM's, except RG:Z, which becomes the pool as for any input. Even NM and MD come back: tags lifted from a FASTQ comment are kept as they are, while the same tags on a SAM or BAM input record are dropped as alignment tags.

> [!NOTE]
> A comment that is not a SAM tag, such as Illumina's `1:N:0:ACGT`, is not kept, and a malformed tag stops the run.

## Diagnostics

Ways to see what a run did and why.

### 19. Explain one read's assignment with --qc-tag

Use this to see why a read went where it did. The `--qc-tag` option writes a JSON summary of the match into a tag, ZS unless another is given, in the FASTQ comment or on the SAM record.

**Data:** the pooled lane of [recipe 7](#7-undetermined-reads-and-a-phix-bin).

```bash
unmux R1.fastq.gz R2.fastq.gz I1.fastq.gz I2.fastq.gz \
  --pool lane1 \
  --group "i7=i7.tsv" --group "i7::loc=2:0:8,dist=1" \
  --group "i5=i5.tsv" --group "i5::loc=3:0:8,dist=1" \
  --group "phix=phix.tsv" --group "phix::loc=0:0:end,dist=1,bothStrands=true" \
  --extract "r1=0:0:end" --extract "r2=1:0:end" \
  --template r1 r2 \
  --sample-sheet samples.tsv \
  --out "fastq/%sample.R%ordinal.fastq.gz" \
  --remove "phix=removed/PhiX.%source.fastq.gz" \
  --unassigned "fastq/Undetermined.%source.fastq.gz" \
  --qc-tag
```

A fig2b read whose i7 had one mismatch, a fig2c read, a PhiX read, and the fig2b read the PhiX rule took with one substitution:

```bash
gzip -dc fastq/N701_N504.R1.fastq.gz | grep '^@M03958:36:000000000-B85D5:1:1101:9622:1073 '
gzip -dc fastq/Undetermined.0.fastq.gz | grep '^@M04856:31:000000000-B34WN:1:1101:17469:1000 '
gzip -dc removed/PhiX.0.fastq.gz | grep '^@A00488:97:HTCHGDSXY:3:1101:1145:1000 '
gzip -dc removed/PhiX.0.fastq.gz | grep '^@M03958:36:000000000-B85D5:1:1101:19685:2435 '
```

```text
@M03958:36:000000000-B85D5:1:1101:9622:1073 ZS:Z:{"v":1,"outcome":"assigned","sample":"N701_N504","sub_sample":null,"groups":[{"g":"i7","tag":"ATTCCGCT","loc":"2:0:8","obs":"NTTCCGCT","sub":1,"ind":0},{"g":"i5","tag":"AGAGTAGA","loc":"3:0:8","obs":"NGAGTAGA","sub":1,"ind":0}]}
@M04856:31:000000000-B34WN:1:1101:17469:1000 ZS:Z:{"v":1,"outcome":"unassigned","reason":"no_sample","groups":[]}
@A00488:97:HTCHGDSXY:3:1101:1145:1000 ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"phix","tag":"GCTGTTCAACCACTAATAGGTAAG","loc":"0:6:30","obs":"CTTACCTATTAGTGGTTGAACAGC","sub":0,"ind":0,"str":"-"}]}
@M03958:36:000000000-B85D5:1:1101:19685:2435 ZS:Z:{"v":1,"outcome":"removed","rule":"phix","group":"phix","groups":[{"g":"i7","tag":"TTCTCCGT","loc":"2:0:8","obs":"TTCTCCGT","sub":0,"ind":0},{"g":"i5","tag":"AGAGTAGA","loc":"3:0:8","obs":"AGAGTAGA","sub":0,"ind":0},{"g":"phix","tag":"TGCCGCGGATTGGTTTCGCTGAAT","loc":"0:199:223","obs":"ATTCAGCGAAAACAATCCGCGGCA","sub":1,"ind":0,"str":"-"}]}
```

Each matched group lists the tag it matched, where (`loc` is file, start and end), the observed bases, and the substitutions and indels it took. The fig2c read matched neither index, so no sample claims it and its reason is `no_sample`; a removed read gives the `--remove` rule that took it, and `str` marks a match on the reverse strand. The last record is how to find out why a read was removed by mistake. The generator checks each record's outcome against the file it is in.

### 20. Read the metrics files

Use the three metrics files to check a run: how many reads each sample got, how many matched no sample and why, and the totals.

**Data:** the pooled lane of [recipe 7](#7-undetermined-reads-and-a-phix-bin), whose command wrote all three.

```bash
cat metrics.summary.tsv
```

```tsv
pool	total_reads	reads_assigned	reads_pass_through	reads_unassigned	reads_removed	frac_assigned	frac_bases_unextracted	mean_read_len_by_input
lane1	27000	19997	0	5000	2003	0.740630	0.000000	289,289,8,17
```

A read is a template: one pair here. Every read is assigned, passed through (when no sample is declared), unassigned or removed, and `frac_assigned` is the assigned share. The `mean_read_len_by_input` column gives the mean length of each input in order, which catches an index read shorter than its tags.

```bash
head -n 4 metrics.per_sample.tsv
```

```tsv
pool	sample	sub_sample	reads	frac_of_pool	frac_bases_unextracted
lane1	N701_N504		1763	0.065296	0.000000
lane1	N702_N504		3221	0.119296	0.000000
lane1	N703_N504		2383	0.088259	0.000000
```

The `frac_bases_unextracted` column is the share of a sample's input bases left out of its output, here the index reads.

```bash
cat metrics.unassigned.tsv
```

```tsv
pool	reason	group	stream	reads	frac_of_pool
lane1	no_sample			5000	0.185185
```

The per-sample reads add up to `reads_assigned` (19,997), the reasons to `reads_unassigned` (5,000), and with the 2,003 removed they make up all 27,000 pairs of the pool; each sample's count, the removed count and the unassigned count also equal the pairs in its files. A reason of `no_sample` means no sample claimed what the groups matched; a group with `minFindsPerGroup` adds rows that give the group and the bound it missed.

### 21. The tag-collision check

Before it reads a record, unmux checks the tags of each group against each other. A set in which one read could match two tags equally well is refused, with the ways out.

**Data:** the 10x and SPLiT-seq inputs from [recipe 11](#11-10x-chromium-cell-barcodes-and-umis) and [recipe 12](#12-split-seq-barcodes-in-three-rounds).

The 10x command with `dist=1` and no mode. Cell barcodes are not designed to be far apart, so some differ in only two places, and a read with one error between them matches both:

```bash
unmux S10_L001_R1_001.fastq.gz S10_L001_R2_001.fastq.gz \
  --pool S10 \
  --group "cb=cells.tsv" \
  --group "cb::loc=0:0:16,dist=1" \
  --extract "cell=@cb" \
  --extract "barcode=0:0:16" \
  --extract "umi=0:16:26" \
  --extract "cdna=1:0:end" \
  --template cdna \
  --tag "CB=cell" \
  --tag "CR=barcode" \
  --tag "UR=umi" --tag "UR::qual=UY" \
  --tag "UB=umi" \
  --out S10.unmapped.bam
```

```text
[ERROR unmux] group `cb`: tags `AAACCTGGTAAATGTG` and `AACCATGGTAAATGTG` can both match one read within dist=1 (the bases `AACCCTGGTAAATGTG` match both); with mode=all this silently mis-assigns reads: separate the tags, lower dist, or use mode=nearest, which only warns and lists the colliding pairs
```

Under the default `mode=all` the first match would silently win. With `mode=nearest`, as in the 10x recipe, such reads are left uncorrected and the check only warns, listing up to ten pairs.

The SPLiT-seq command with one indel allowed (`dist=1:1:2`). Its barcodes are 4 substitutions apart, but an insertion or deletion can bring a read within reach of two:

<details>
<summary>The command</summary>

```bash
unmux R1.fastq.gz R2.fastq.gz \
  --group "round1=round1.tags.txt" \
  --group "round1::loc=1:78:86" \
  --group "round1::dist=1:1:2" \
  --group "round1::minFindsPerGroup=1" \
  --group "round2=round2_3.tags.txt" \
  --group "round2::loc=1:48:56" \
  --group "round2::dist=1:1:2" \
  --group "round2::minFindsPerGroup=1" \
  --group "round3=round2_3.tags.txt" \
  --group "round3::loc=1:10:18" \
  --group "round3::dist=1:1:2" \
  --group "round3::minFindsPerGroup=1" \
  --extract "cdna=0:0:end" \
  --extract "umi=1:0:10" \
  --extract "bc1=@round1" \
  --extract "bc2=@round2" \
  --extract "bc3=@round3" \
  --template "cdna" \
  --tag "CB=bc1+bc2+bc3" \
  --tag "CB::qual=CY" \
  --tag "CB::sep=-" \
  --tag "CB::qual-sep=%20" \
  --tag "RX=umi" \
  --tag "RX::qual=QX" \
  --out "splitseq.unmapped.bam"
```

</details>

```text
[ERROR unmux] group `round1`: 823 pairs of tags can both match one read within dist=1:1 only through an indel, such as `ACTCGTAA` and `AAACGATA` (the bases `AATCGTA` match both); set mode=nearest to leave such reads unassigned (recommended), or write mode=all to keep assigning each to one of the tags by tie-break
```

Here unmux asks for a mode. With `mode=nearest` such reads are left unassigned, and with `mode=all` each is still assigned by a tie-break, with warnings that list the pairs; [recipe 22](#22-search-windows-wider-than-their-tags) compares the two. For both commands the generator checks the exit status, the message, and that no file was written.

> [!TIP]
> For a set designed to be far apart, such as Illumina indexes, a refusal usually means `dist` is too large for the set.

### 22. Search windows wider than their tags

Widening a window lets a barcode that an insertion or deletion shifted still be found. But in a window wider than its tags, a read with an error in its barcode can also match another tag, shifted, as closely as its own. unmux estimates how often that happens and warns.

**Data:** the SPLiT-seq reads of [recipe 12](#12-split-seq-barcodes-in-three-rounds), of which 2,884 miss a round in the README windows.

The README command with each window 1 bp wider on both sides, one indel allowed, and `mode=all` as the collision check of [recipe 21](#21-the-tag-collision-check) asks:

```bash
unmux R1.fastq.gz R2.fastq.gz \
  --group "round1=round1.tags.txt" \
  --group "round1::loc=1:77:87" \
  --group "round1::dist=1:1:2,mode=all" \
  --group "round1::minFindsPerGroup=1" \
  --group "round2=round2_3.tags.txt" \
  --group "round2::loc=1:47:57" \
  --group "round2::dist=1:1:2,mode=all" \
  --group "round2::minFindsPerGroup=1" \
  --group "round3=round2_3.tags.txt" \
  --group "round3::loc=1:9:19" \
  --group "round3::dist=1:1:2,mode=all" \
  --group "round3::minFindsPerGroup=1" \
  --extract "cdna=0:0:end" \
  --extract "umi=1:0:10" \
  --extract "bc1=@round1" \
  --extract "bc2=@round2" \
  --extract "bc3=@round3" \
  --template "cdna" \
  --tag "CB=bc1+bc2+bc3" \
  --tag "CB::qual=CY" \
  --tag "CB::sep=-" \
  --tag "CB::qual-sep=%20" \
  --tag "RX=umi" \
  --tag "RX::qual=QX" \
  --out "splitseq.unmapped.bam"
```

```text
[WARN  unmux::demux] group `round1`: the loc window is 2 bp wider than the tags, so for an estimated 15.7% of reads with an error in the barcode another tag matches at a shifted position as closely as the true tag, and the leftmost match wins; use mode=nearest
[WARN  unmux::demux] group `round2`: the loc window is 2 bp wider than the tags, so for an estimated 15.7% of reads with an error in the barcode another tag matches at a shifted position as closely as the true tag, and the leftmost match wins; use mode=nearest
[WARN  unmux::demux] group `round3`: the loc window is 2 bp wider than the tags, so for an estimated 15.7% of reads with an error in the barcode another tag matches at a shifted position as closely as the true tag, and the leftmost match wins; use mode=nearest
```

The run also prints 33 warnings that list tag pairs one indel apart, as `mode=all` does. Following the advice, the same command with `mode=nearest` in each round:

<details>
<summary>The command</summary>

```bash
unmux R1.fastq.gz R2.fastq.gz \
  --group "round1=round1.tags.txt" \
  --group "round1::loc=1:77:87" \
  --group "round1::dist=1:1:2,mode=nearest" \
  --group "round1::minFindsPerGroup=1" \
  --group "round2=round2_3.tags.txt" \
  --group "round2::loc=1:47:57" \
  --group "round2::dist=1:1:2,mode=nearest" \
  --group "round2::minFindsPerGroup=1" \
  --group "round3=round2_3.tags.txt" \
  --group "round3::loc=1:9:19" \
  --group "round3::dist=1:1:2,mode=nearest" \
  --group "round3::minFindsPerGroup=1" \
  --extract "cdna=0:0:end" \
  --extract "umi=1:0:10" \
  --extract "bc1=@round1" \
  --extract "bc2=@round2" \
  --extract "bc3=@round3" \
  --template "cdna" \
  --tag "CB=bc1+bc2+bc3" \
  --tag "CB::qual=CY" \
  --tag "CB::sep=-" \
  --tag "CB::qual-sep=%20" \
  --tag "RX=umi" \
  --tag "RX::qual=QX" \
  --out "splitseq.unmapped.bam"
```

</details>

It prints no window warning: under `mode=nearest`, a window wider than the tags by twice the indel budget is as wide as unmux advises.

Real reads have no known truth, so the table scores each run two ways: against the 568 cells split-pipe called, whose barcodes a correct CB most likely belongs to, and against the README windows, which take no shifted match:

| Command | Written | With a CB | CB of a called cell | Other CB | No CB | CB changed from the README windows |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| README windows, `dist=1` | 7,116 | 7,116 | 5,453 | 1,663 | 0 |  |
| 1 bp wider, `dist=1:1:2,mode=all` | 9,238 | 9,238 | 6,266 | 2,972 | 0 | 106 |
| 1 bp wider, `dist=1:1:2,mode=nearest` | 9,238 | 8,168 | 6,159 | 2,009 | 1,070 | 1 |

The wider windows find 813 more reads of called cells, but under `mode=all` they also change the CB of 106 reads that the README windows assign, and add 1,309 reads whose CB is no called cell's. Under `mode=nearest` fewer CBs change (1), but 1,070 reads that tie in some round are written with no CB at all:

```bash
samtools view -c splitseq.unmapped.bam
samtools view -c -d CB splitseq.unmapped.bam
```

```text
9238
8168
```

> [!WARNING]
> Under `mode=nearest`, a read that ties between two tags of a round is still written, without the CB built from that round, though every round sets `minFindsPerGroup=1`; its `--qc-tag` lists only the rounds that matched. Keep only reads with a CB, for example with `samtools view -d CB`.

## Data sources

Every download is pinned: GitHub files by commit, SRA runs by accession and spot range, NCBI sequences by accession and version.

| Data | Pinned URL |
| --- | --- |
| SRA run SRR390728, spots 1 to 2,000 | <https://www.ncbi.nlm.nih.gov/sra/SRR390728> |
| UDiTaS fig2b run, R1 (first 20,000 records) | <https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2b/Undetermined_S0_L001_R1_001.fastq.gz> |
| UDiTaS fig2b run, R2 (first 20,000 records) | <https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2b/Undetermined_S0_L001_R2_001.fastq.gz> |
| UDiTaS fig2b run, I1 (first 20,000 records) | <https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2b/Undetermined_S0_L001_I1_001.fastq.gz> |
| UDiTaS fig2b run, I2 (first 20,000 records) | <https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2b/Undetermined_S0_L001_I2_001.fastq.gz> |
| UDiTaS fig2b sample sheet | <https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2b/sample_info.csv> |
| UDiTaS fig2c run, R1 (first 5,000 records) | <https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2c/Undetermined_S0_L001_R1_001.fastq.gz> |
| UDiTaS fig2c run, R2 (first 5,000 records) | <https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2c/Undetermined_S0_L001_R2_001.fastq.gz> |
| UDiTaS fig2c run, I1 (first 5,000 records) | <https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2c/Undetermined_S0_L001_I1_001.fastq.gz> |
| UDiTaS fig2c run, I2 (first 5,000 records) | <https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/data/fig2c/Undetermined_S0_L001_I2_001.fastq.gz> |
| CHANGE-seq test reads TEST, r1 | <https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf/test/data/input/TEST.r1.fastq.gz> |
| CHANGE-seq test reads TEST, r2 | <https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf/test/data/input/TEST.r2.fastq.gz> |
| CHANGE-seq test reads TEST_control, r1 | <https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf/test/data/input/TEST_control.r1.fastq.gz> |
| CHANGE-seq test reads TEST_control, r2 | <https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf/test/data/input/TEST_control.r2.fastq.gz> |
| CHANGE-seq read merging, mergeReads() | <https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf/changeseq/mergeReads.py> |
| SRA runs of BioProject PRJNA730334 batch VERD6, one per native barcode, spots 1 to 200 of each | <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA730334> |
| PacBio pbAmpliconAnalysis_CYP2D6 HiFi reads before lima, first 10 MB | <https://downloads.pacbcloud.com/public/dataset/pbAmpliconAnalysis_CYP2D6/workflow_inputs/m64013e_210223_231140.hifi_reads.bam> |
| PacBio CYP2D6 barcode pairs and samples | <https://downloads.pacbcloud.com/public/dataset/pbAmpliconAnalysis_CYP2D6/workflow_inputs/biosamples.csv> |
| PacBio Sequel 96 barcodes v1 | <https://raw.githubusercontent.com/PacificBiosciences/hifi-amplicon-workflow/2964cd1fbddae04243d56914935881cca01f6d13/data/Sequel_96_barcodes_v1.fasta> |
| splitcode tutorial: SPLiT-seq processing | <https://raw.githubusercontent.com/pachterlab/splitcode-tutorial/68387770995199a9ac4ccadc25ee2f8e96cd24be/docs/source/tutorials_splitseq.rst> |
| splitcode tutorial round-1 random-hexamer barcodes | <https://raw.githubusercontent.com/pachterlab/splitcode-tutorial/68387770995199a9ac4ccadc25ee2f8e96cd24be/uploads/splitseq/r1_R.txt> |
| splitcode tutorial round-1 oligo-dT barcodes | <https://raw.githubusercontent.com/pachterlab/splitcode-tutorial/68387770995199a9ac4ccadc25ee2f8e96cd24be/uploads/splitseq/r1_T.txt> |
| splitcode tutorial round-2 and round-3 barcodes | <https://raw.githubusercontent.com/pachterlab/splitcode-tutorial/68387770995199a9ac4ccadc25ee2f8e96cd24be/uploads/splitseq/r2_r3.txt> |
| SRA runs DRR021959 to DRR021990 of study DRP001067, spots 1 to 250 of each | <https://www.ncbi.nlm.nih.gov/bioproject/PRJDB2148> |
| SRA run SRR13948565, spots 1,000,001 to 1,010,000 | <https://www.ncbi.nlm.nih.gov/sra/SRR13948565> |
| GEO GSM5169184, the cells split-pipe called in SRR13948565 | <https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5169nnn/GSM5169184/suppl/GSM5169184_C2C12_short_1k_cell_metadata.csv.gz> |
| nf-core/test-datasets (modules) 10x antibody reads, read 1 | <https://raw.githubusercontent.com/nf-core/test-datasets/e8d2f4f2af42cc757ed2448c2398e282e9990235/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/ab/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz> |
| nf-core/test-datasets (modules) 10x antibody reads, read 2 | <https://raw.githubusercontent.com/nf-core/test-datasets/e8d2f4f2af42cc757ed2448c2398e282e9990235/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/fastqs/ab/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz> |
| nf-core/test-datasets (modules) antibody feature reference | <https://raw.githubusercontent.com/nf-core/test-datasets/e8d2f4f2af42cc757ed2448c2398e282e9990235/data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells/5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv> |
| SRA run SRR27862880, spots 1 to 2,000 | <https://www.ncbi.nlm.nih.gov/sra/SRR27862880> |
| NCBI RefSeq NC_001422.1, the PhiX174 genome | <https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text> |
| UDiTaS demultiplexing rules, demultiplex() | <https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a/uditas/uditas_helpers.py> |
| samtools 1.25 test file test/split/split.sam | <https://raw.githubusercontent.com/samtools/samtools/7162950d0988480faeddeb49f03520737e13a9c4/test/split/split.sam> |
| modkit test file tests/resources/fwd_rev_modbase_records.bam | <https://raw.githubusercontent.com/nanoporetech/modkit/5cecc3fb3a9336068d9e3c68d5c08d678153dd2c/tests/resources/fwd_rev_modbase_records.bam> |
| modkit test file tests/resources/bc_anchored_10_reads.unmapped.bam | <https://raw.githubusercontent.com/nanoporetech/modkit/5cecc3fb3a9336068d9e3c68d5c08d678153dd2c/tests/resources/bc_anchored_10_reads.unmapped.bam> |
| nf-core/test-datasets (scrnaseq) 10x v2 read 1 | <https://raw.githubusercontent.com/nf-core/test-datasets/d934d6e8367fe2626184496b1889671cf2b02dab/testdata/S10_L001_R1_001.fastq.gz> |
| nf-core/test-datasets (scrnaseq) 10x v2 read 2 | <https://raw.githubusercontent.com/nf-core/test-datasets/d934d6e8367fe2626184496b1889671cf2b02dab/testdata/S10_L001_R2_001.fastq.gz> |
| nf-core/test-datasets (scrnaseq) 10x v2 barcode keeplist | <https://raw.githubusercontent.com/nf-core/test-datasets/d934d6e8367fe2626184496b1889671cf2b02dab/reference/kallisto/10xv2_whitelist.txt> |
| nf-core/test-datasets (fastquorum) sample sheet for SRR6109255 | <https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/samplesheets/samplesheet.tiny.csv> |
| nf-core/test-datasets (fastquorum) duplex reads SRR6109255, read 1 | <https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/fastqs/tiny/SRR6109255_1.fastq.gz> |
| nf-core/test-datasets (fastquorum) duplex reads SRR6109255, read 2 | <https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/fastqs/tiny/SRR6109255_2.fastq.gz> |
| nf-core/test-datasets (fastquorum) sample sheet for the fixed UMIs | <https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/samplesheets/samplesheet.nonrandom.csv> |
| nf-core/test-datasets (fastquorum) reads with fixed UMIs, read 1 | <https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/fastqs/synthetic-nonrandom-umi/synthetic_nonrandom_umi_1.fastq.gz> |
| nf-core/test-datasets (fastquorum) reads with fixed UMIs, read 2 | <https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/fastqs/synthetic-nonrandom-umi/synthetic_nonrandom_umi_2.fastq.gz> |
| nf-core/test-datasets (fastquorum) IDT xGen cfDNA UMI list | <https://raw.githubusercontent.com/nf-core/test-datasets/2ca9f8176f1f40ec1124c5182c58524766c95411/testdata/umis/idt-xgen-cfdna-umis.txt> |
| split-seq-pipeline ligation barcodes (rounds 2 and 3) | <https://raw.githubusercontent.com/yjzhang/split-seq-pipeline/c3923ea44a94a1fc0effe5ff8e2ca479c0bff934/split_seq/barcodes/bc_8nt_v1.csv> |
| split-seq-pipeline reverse-transcription barcodes (round 1) | <https://raw.githubusercontent.com/yjzhang/split-seq-pipeline/c3923ea44a94a1fc0effe5ff8e2ca479c0bff934/split_seq/barcodes/bc_8nt_v2.csv> |
| split-seq-pipeline read layout (amp_seq, chemistry v2) | <https://raw.githubusercontent.com/yjzhang/split-seq-pipeline/c3923ea44a94a1fc0effe5ff8e2ca479c0bff934/split_seq/tools.py> |
| dorado native barcodes NB01 to NB24 and their flanks | <https://raw.githubusercontent.com/nanoporetech/dorado/651980f48c21eb93a5c46c65e1f004136f02c07a/dorado/utils/barcode_kits.cpp> |
| dorado test read written with U bases | <https://raw.githubusercontent.com/nanoporetech/dorado/651980f48c21eb93a5c46c65e1f004136f02c07a/tests/data/bam_reader/fastq_with_us_and_tags.fq> |

## Regenerating this page

One command fetches the inputs into `docs/recipes/cache/`, builds unmux, runs every recipe with its checks, and rewrites this page. It needs Python 3, samtools 1.16 or newer, fg-sra, fqtk, splitcode and, for the lima comparison, Docker.

```bash
make -C docs/recipes            # fetch, build, run, and write docs/recipes.md
make -C docs/recipes check      # fail if the page would change
make -C docs/recipes run UNMUX=/path/to/unmux
```

The targets are `fetch`, `run`, `page`, `check` and `clean`, and the variables `SAMTOOLS`, `FG_SRA`, `FQTK`, `SPLITCODE` and `LIMA` give other ways to run those tools, for example `SAMTOOLS="pixi exec -c conda-forge -c bioconda -s samtools samtools"`.
