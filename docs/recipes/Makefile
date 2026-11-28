# Regenerate the unmux cookbook, docs/recipes.md, from pinned public data.
#
#   make -C docs/recipes            the same as `page`
#   make -C docs/recipes fetch      fetch the pinned inputs into cache/
#   make -C docs/recipes run        run every recipe and its checks
#   make -C docs/recipes page       run them and write docs/recipes.md
#   make -C docs/recipes check      fail if the regenerated page differs from docs/recipes.md
#   make -C docs/recipes clean      remove cache/
#
# Set UNMUX to test a given binary; otherwise `cargo build --release` builds one. The tool
# commands below can be overridden, for example
# SAMTOOLS="pixi exec -c conda-forge -c bioconda -s samtools samtools".

REPO := $(abspath $(CURDIR)/../..)
PYTHON ?= python3
UNMUX ?=
SAMTOOLS ?= samtools
FG_SRA ?= fg-sra
FQTK ?= fqtk
SPLITCODE ?= splitcode
LIMA ?= docker run --rm --platform linux/amd64 -v "$$PWD":/w -w /w quay.io/biocontainers/lima:26.2.1--h9ee0642_0 lima
CARGO_BUILD ?= cargo build --release --manifest-path $(REPO)/Cargo.toml

ifeq ($(strip $(UNMUX)),)
BINARY := $(REPO)/target/release/unmux
BUILD := $(CARGO_BUILD)
else
BINARY := $(UNMUX)
BUILD := true
endif

GENERATE = $(PYTHON) -I generate.py --unmux "$(BINARY)" --samtools '$(SAMTOOLS)' --fg-sra '$(FG_SRA)' --fqtk '$(FQTK)' --splitcode '$(SPLITCODE)' --lima '$(LIMA)'

.PHONY: page fetch run check clean build

page: build
	$(GENERATE) page

fetch:
	$(GENERATE) fetch

run: build
	$(GENERATE) run

check: build
	$(GENERATE) check

build:
	$(BUILD)

clean:
	rm -rf cache
