# sasscode

[![Install with bioconda](https://img.shields.io/badge/Install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/sasscode/README.html)
<!--TODO: [![Anaconda Version](https://anaconda.org/bioconda/sasscode/badges/version.svg)](http://bioconda.github.io/recipes/sasscode/README.html)-->
[![Build Status](https://github.com/clintval/sasscode/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/clintval/sasscode/actions/workflows/ci.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/clintval/sasscode/badge.svg?branch=main)](https://coveralls.io/github/clintval/sasscode?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Flexible read parsing and demultiplexing to FASTX/SAM/BAM/CRAM, splitcode-style.

![sasscode](.github/img/cover.jpg)

Install with mamba, conda, or run directly with pixi:

```bash
pixi exec \
    -c conda-forge -c bioconda \
    sasscode --help
```

## Introduction

The tool `sasscode` reads multiple FASTX/SAM/BAM/CRAM inputs, identifies and extracts technical sequences (barcodes, UMIs, adapters) with error tolerance using the [`sassy`](https://github.com/RagnarGrootKoerkamp/sassy) approximate matcher, and writes FASTX/SAM/BAM/CRAM files with preserved per-read segment qualities, fanning a read pool out into per-sample (and optionally per-sample and per-library) files in a single pass.

It aims to do in one shot what otherwise would take a combination of separate tools, including:

- [`fgbio FastqToBam`](https://github.com/fulcrumgenomics/fgbio): stitch multiple FASTQs into one unmapped BAM via read structures and SAM tags
- [`fqtk`](https://github.com/fulcrumgenomics/fqtk): fast sample demultiplexing driven by per-read barcode structures
- [`Picard FastqToSam`](https://github.com/broadinstitute/picard): convert raw FASTQs into an unmapped BAM
- [`qualrepair`](https://bioconda.github.io/recipes/qualrepair/README.html): repair the base qualities that splitcode mangles during extraction
- [`samtools split`](https://github.com/samtools/samtools): split one BAM into per-read-group (per-sample) BAMs
- [`splitcode`](https://github.com/pachterlab/splitcode): identify, extract, and edit technical sequences from a declarative config
- [`UMI-tools extract`](https://github.com/CGATOxford/UMI-tools): pull UMIs out of reads and onto the read names

## Quick Start

Process a SPLiT-seq run (splitcode's [SPLiT-seq example](https://splitcode.readthedocs.io/en/latest/tutorials_splitseq.html)) in one call: three rounds of 8 bp cell barcodes plus a 10 bp UMI on R2, with the cDNA on R1, emitted as an unmapped BAM with `CB` and `RX` tags and their associated quality scores in tags `CY` and `QX`:

```bash
sasscode R1.fastq.gz R2.fastq.gz \
  --group round1=round1.tags.txt \
  --group round1::loc=1:78:86 \
  --group round1::dist=1 \
  --group round1::minFindsPerGroup=1 \
  --group round2=round2_3.tags.txt \
  --group round2::loc=1:48:56 \
  --group round2::dist=1 \
  --group round2::minFindsPerGroup=1 \
  --group round3=round2_3.tags.txt \
  --group round3::loc=1:10:18 \
  --group round3::dist=1 \
  --group round3::minFindsPerGroup=1 \
  --extract umi=1:0:10 \
  --extract bc1=@round1 \
  --extract bc2=@round2 \
  --extract bc3=@round3 \
  --extract cdna=0:0:end \
  --template cdna \
  --tag CB=bc1,bc2,bc3 \
  --tag CB::qual=CY \
  --tag CB::sep='-' \
  --tag CB::qual-sep=' ' \
  --tag RX=umi \
  --tag RX::qual=QX \
  --out splitseq.unmapped.bam
```

## Features

- **FASTX/SAM/BAM/CRAM in and out** (SAM/BAM/CRAM written unmapped), with per-segment qualities carried through extraction.
- **splitcode-style matching**: tag groups, variable-length tags, location windows, mismatch/indel tolerance, sequential `next`/`previous` anchoring, `@extract` spans, and error-correction to canonical barcodes.
- **Single-pass demultiplexing, optionally nested**: split a pool into per-sample outputs and, if you like, each sample into sub-samples (a library of samples, sub-libraries of lysates, and so on), each written as its own FASTX/SAM/BAM/CRAM with `SM`/`LB` read group identifiers.
- **Configured on the CLI or in sheets** with the same grammar: flags repeat and accumulate, so a spec can be built up piece by piece (`--group g1::loc=... --group g1::dist=...`), or written compactly with comma lists (`--group g1::loc=...,dist=...`).

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
