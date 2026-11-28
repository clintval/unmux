"""Shared helpers of the UMI recipes: the read structure an nf-core/fastquorum sample sheet gives a FASTQ."""

import csv


def read_structure(sheet, fastq_name):
    with open(sheet, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["fastq_1"].endswith("/" + fastq_name):
                return row["read_structure"]
    raise SystemExit(f"{sheet} has no row for {fastq_name}")
