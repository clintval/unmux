import illumina

UDITAS_INPUTS = [
    *[("uditas", f"Undetermined_S0_L001_{r}_001.fastq", f"Undetermined_S0_L001_{r}_001.fastq.gz") for r in ("R1", "R2", "I1", "I2")],
    ("uditas", "sample_info.csv", "sample_info.csv"),
    ("uditas", "uditas_helpers.py", None),
]
UDITAS_SHEET = """awk -F , 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; print "sample\\ti7\\ti5"; next }
  { print $c["index_I1"] "_" $c["index_I2"] "\\t" $c["index_I1"] "\\t" $c["index_I2"] }' sample_info.csv > samples.tsv
awk -F , 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; print "id\\tseq"; next }
  !seen[$c["index_I1"]]++ { print $c["index_I1"] "\\t" $c["barcode_I1"] }' sample_info.csv > i7.tsv
awk -F , 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; print "id\\tseq"; next }
  !seen[$c["index_I2"]]++ { print $c["index_I2"] "\\t" $c["barcode_I2"] }' sample_info.csv > i5.tsv"""
UDITAS_COMMAND = """unmux Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz \\
  Undetermined_S0_L001_I1_001.fastq.gz Undetermined_S0_L001_I2_001.fastq.gz \\
  --pool fig2b \\
  --group "i7=i7.tsv" --group "i7::loc=2:0:8,dist=1" \\
  --group "i5=i5.tsv" --group "i5::loc=3:0:8,dist=1" \\
  --extract "r1=0:0:end" --extract "r2=1:0:end" --extract "umi=3:8:end" \\
  --template r1 r2 \\
  --tag "RX=umi" \\
  --sample-sheet samples.tsv \\
  --rg-tag PL=ILLUMINA \\
  --out "bam/%sample.unmapped.bam" \\
  --metrics-per-sample metrics.per_sample.tsv"""
FQTK_UDITAS_SHEET = """awk 'BEGIN { print "sample_id\\tbarcode" } FNR == 1 { next }
  FILENAME == "i7.tsv" { i7[$1] = $2; next } FILENAME == "i5.tsv" { i5[$1] = $2; next }
  { print $1 "\\t" i7[$2] i5[$3] }' i7.tsv i5.tsv samples.tsv > fqtk.samples.tsv"""
FQTK_UDITAS = """fqtk demux \\
  --inputs Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz \\
    Undetermined_S0_L001_I1_001.fastq.gz Undetermined_S0_L001_I2_001.fastq.gz \\
  --read-structures +T +T 8B 8B+M \\
  --sample-metadata fqtk.samples.tsv \\
  --max-mismatches 2 --min-mismatch-delta 1 \\
  --umi-in-name false \\
  --output fqtk"""
UDITAS = {
    "index": "gzip -dc Undetermined_S0_L001_I1_001.fastq.gz | head -n 4\ngzip -dc Undetermined_S0_L001_I2_001.fastq.gz | head -n 4",
    "sheet": UDITAS_SHEET,
    "sheets": "head -n 4 samples.tsv\nhead -n 3 i7.tsv\ncat i5.tsv",
    "demux": UDITAS_COMMAND,
    "counts": 'for f in bam/*.bam; do echo "$f $(samtools view -c -f 64 "$f")"; done',
    "view": "samtools view bam/N702_N504.unmapped.bam | head -n 2 | cut -f 1,2,12-",
    "fqtk_version": "fqtk --version",
    "fqtk_sheet": FQTK_UDITAS_SHEET,
    "fqtk_head": "head -n 3 fqtk.samples.tsv",
    "fqtk": FQTK_UDITAS,
    "fqtk_defaults": FQTK_UDITAS.replace("  --max-mismatches 2 --min-mismatch-delta 1 \\\n", "").replace(
        "--output fqtk", "--output fqtk-defaults"
    ),
}


rule uditas_fetch:
    input:
        downloads(UDITAS_INPUTS),
    output:
        [work("uditas", path) for _, _, path in UDITAS_INPUTS if path],
    log:
        logfile("uditas.fetch"),
    params:
        inputs=UDITAS_INPUTS,
        recipe="uditas",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule uditas_run:
    input:
        rules.uditas_fetch.output,
        TOOLS,
    output:
        run=work("uditas", "run.json"),
        sheets=[work("uditas", name) for name in ("samples.tsv", "i7.tsv", "i5.tsv")],
    threads: 8
    log:
        logfile("uditas.run"),
    params:
        commands=UDITAS,
        used=UDITAS_INPUTS,
        recipe="uditas",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)


BINS_INPUTS = [
    *[(key, f"Undetermined_S0_L001_{r}_001.fastq", None) for key in ("uditas", "uditas-fig2c") for r in ("R1", "R2", "I1", "I2")],
    ("sra-SRR27862880", "SRR27862880_1.fastq", None),
    ("sra-SRR27862880", "SRR27862880_2.fastq", None),
    ("ncbi", "NC_001422.1.fasta", "NC_001422.1.fasta"),
]
BINS_FROM_UDITAS = ["samples.tsv", "i7.tsv", "i5.tsv", "sample_info.csv"]
BINS_COMMAND = """unmux R1.fastq.gz R2.fastq.gz I1.fastq.gz I2.fastq.gz \\
  --pool lane1 \\
  --group "i7=i7.tsv" --group "i7::loc=2:0:8,dist=1" \\
  --group "i5=i5.tsv" --group "i5::loc=3:0:8,dist=1" \\
  --group "phix=phix.tsv" --group "phix::loc=0:0:end,dist=1,bothStrands=true" \\
  --extract "r1=0:0:end" --extract "r2=1:0:end" \\
  --template r1 r2 \\
  --sample-sheet samples.tsv \\
  --out "fastq/%sample.R%ordinal.fastq.gz" \\
  --remove "phix=removed/PhiX.%source.fastq.gz" \\
  --unassigned "fastq/Undetermined.%source.fastq.gz" \\
  --metrics-per-sample metrics.per_sample.tsv \\
  --metrics-summary metrics.summary.tsv \\
  --metrics-unassigned metrics.unassigned.tsv"""
PHIX_KMERS = """awk '!/^>/ { seq = seq $0 } END { print "id\\tseq"
  for (i = 1; i + 23 <= length(seq); i += 12) print "phix_" i "\\t" substr(seq, i, 24) }' NC_001422.1.fasta > phix.tsv"""
BINS = {
    "kmers": PHIX_KMERS,
    "demux": BINS_COMMAND,
    "listing": (
        "for f in removed/*.fastq.gz fastq/Undetermined.*.fastq.gz; do\n"
        "  echo \"$f $(gzip -dc \"$f\" | wc -l | awk '{print $1 / 4}')\"\n"
        "done"
    ),
    "unassigned": "cat metrics.unassigned.tsv",
    "exact": (
        BINS_COMMAND.replace("phix::loc=0:0:end,dist=1", "phix::loc=0:0:end,dist=0")
        .replace("fastq/", "exact/")
        .replace("removed/", "exact/removed.")
        .replace("metrics.", "exact/metrics.")
    ),
    "strict": (
        BINS_COMMAND.replace('"i7::loc=2:0:8,dist=1"', '"i7::loc=2:0:8,dist=1,minFindsPerGroup=1"')
        .replace('"i5::loc=3:0:8,dist=1"', '"i5::loc=3:0:8,dist=1,minFindsPerGroup=1"')
        .replace("fastq/", "strict/")
        .replace("removed/", "strict/removed.")
        .replace("metrics.", "strict/metrics.")
    ),
}


rule bins_fetch:
    input:
        downloads(BINS_INPUTS),
        uditas=[work("uditas", name) for name in BINS_FROM_UDITAS],
    output:
        [work("bins", f"{r}.fastq.gz") for r in ("R1", "R2", "I1", "I2")],
        [work("bins", path) for _, _, path in BINS_INPUTS if path],
        [work("bins", name) for name in BINS_FROM_UDITAS],
    log:
        logfile("bins.fetch"),
    params:
        inputs=BINS_INPUTS,
        recipe="bins",
        copies=lambda wildcards, input: list(zip(input.uditas, BINS_FROM_UDITAS)),
        after=[("illumina.stage_bins", "")],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule bins_run:
    input:
        rules.bins_fetch.output,
        TOOLS,
    output:
        run=work("bins", "run.json"),
        phix=work("bins", "phix.tsv"),
        summary=work("bins", "metrics.summary.tsv"),
        per_sample=work("bins", "metrics.per_sample.tsv"),
        unassigned=work("bins", "metrics.unassigned.tsv"),
    threads: 8
    log:
        logfile("bins.run"),
    params:
        commands=BINS,
        used=BINS_INPUTS,
        recipe="bins",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)


INLINE_INPUTS = [("sra-carabus", f"{run}.fastq", None) for run, _, _ in illumina.CARABUS]
INLINE_COMMAND = """unmux pool.fastq.gz \\
  --pool D0BTPACXX.8 \\
  --group "bc=barcodes.tsv" --group "bc::loc=0:0:5,dist=1" \\
  --group "site={PstI=TGCAG}" --group "site::loc=0:5:10,minFindsPerGroup=1" \\
  --extract "barcode=0:0:5" --extract "insert=0:5:end" \\
  --template insert \\
  --tag "BC=barcode" \\
  --sample-from-group bc \\
  --rg-tag PL=ILLUMINA \\
  --out "bam/%sample.unmapped.bam\""""
FQTK_RAD_SHEET = """awk 'BEGIN { print "sample_id\\tbarcode" } NR > 1 { print $1 "\\t" $2 }' barcodes.tsv > fqtk.samples.tsv"""
FQTK_RAD = """fqtk demux \\
  --inputs pool.fastq.gz \\
  --read-structures 5B+T \\
  --sample-metadata fqtk.samples.tsv \\
  --max-mismatches 1 --min-mismatch-delta 1 \\
  --output fqtk"""
INLINE = {
    "head": "head -n 4 barcodes.tsv\ngzip -dc pool.fastq.gz | head -n 2",
    "demux": INLINE_COMMAND,
    "exact": INLINE_COMMAND.replace("bc::loc=0:0:5,dist=1", "bc::loc=0:0:5,dist=0").replace("bam/", "exact/"),
    "fqtk_version": "fqtk --version",
    "fqtk": f"{FQTK_RAD_SHEET}\n{FQTK_RAD}",
}
INLINE_TEMPLATES = {"view": "samtools view bam/{{sample}}.unmapped.bam | awk '$1 == \"{{read}}\"' | cut -f 1,2,12-"}


rule inline_fetch:
    input:
        downloads(INLINE_INPUTS),
    output:
        work("inline", "pool.fastq.gz"),
        work("inline", "barcodes.tsv"),
    log:
        logfile("inline.fetch"),
    params:
        inputs=INLINE_INPUTS,
        recipe="inline",
        copies=[],
        after=[("illumina.stage_inline", "")],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule inline_run:
    input:
        rules.inline_fetch.output,
        TOOLS,
    output:
        run=work("inline", "run.json"),
    threads: 8
    log:
        logfile("inline.run"),
    params:
        commands=INLINE,
        templates=INLINE_TEMPLATES,
        used=INLINE_INPUTS,
        recipe="inline",
    run:
        cookbook.record_rule(CACHE, output, params, log)
