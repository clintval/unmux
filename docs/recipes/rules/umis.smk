DUPLEX_INPUTS = [
    ("fastquorum", "SRR6109255_1.fastq.gz", "SRR6109255_1.fastq.gz"),
    ("fastquorum", "SRR6109255_2.fastq.gz", "SRR6109255_2.fastq.gz"),
    ("fastquorum", "samplesheet.tiny.csv", None),
]
DUPLEX_COMMAND = """unmux SRR6109255_1.fastq.gz SRR6109255_2.fastq.gz \\
  --extract "umi1=0:0:10" --extract "umi2=1:0:10" \\
  --extract "r1=0:11:end" --extract "r2=1:11:end" \\
  --template r1 r2 \\
  --tag "RX=umi1+umi2" \\
  --out SRR6109255.unmapped.bam"""
SWAP_COMMAND = """samtools view -f 64 SRR6109255.unmapped.bam | grep -o 'RX:Z:[^[:space:]]*' | cut -c 6- \\
  | sort -u > umi_pairs.txt
awk 'BEGIN { FS = "-" } { print $2 "-" $1 }' umi_pairs.txt | sort | comm -12 - umi_pairs.txt \\
  | awk 'END { print NR }'"""
DUPLEX = {
    "inputs": "gzip -dc SRR6109255_1.fastq.gz | head -n 4 | cut -c 1-60\ngzip -dc SRR6109255_2.fastq.gz | head -n 4 | cut -c 1-60",
    "duplex": DUPLEX_COMMAND,
    "view": "samtools view SRR6109255.unmapped.bam | head -n 2 | cut -f 1,2,12-\nsamtools view SRR6109255.unmapped.bam | head -n 2 | cut -f 10 | cut -c 1-49",
    "swap": SWAP_COMMAND,
}


rule duplex_fetch:
    input:
        downloads(DUPLEX_INPUTS),
    output:
        [work("duplex", path) for _, _, path in DUPLEX_INPUTS if path],
    log:
        logfile("duplex.fetch"),
    params:
        inputs=DUPLEX_INPUTS,
        recipe="duplex",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule duplex_run:
    input:
        rules.duplex_fetch.output,
        TOOLS,
    output:
        run=work("duplex", "run.json"),
    threads: 8
    log:
        logfile("duplex.run"),
    params:
        commands=DUPLEX,
        templates={},
        used=DUPLEX_INPUTS,
        recipe="duplex",
    run:
        cookbook.record_rule(CACHE, output, params, log)


DUPLEX_FIXED_INPUTS = [
    ("fastquorum", "synthetic_nonrandom_umi_1.fastq.gz", "synthetic_nonrandom_umi_1.fastq.gz"),
    ("fastquorum", "synthetic_nonrandom_umi_2.fastq.gz", "synthetic_nonrandom_umi_2.fastq.gz"),
    ("fastquorum", "idt-xgen-cfdna-umis.txt", "idt-xgen-cfdna-umis.txt"),
    ("fastquorum", "samplesheet.nonrandom.csv", None),
    ("fastquorum", "SRR6109255_1.fastq.gz", None),
    ("fastquorum", "SRR6109255_2.fastq.gz", None),
]
UMI_PREP = """awk 'BEGIN { print "id\\tseq" } { print $1 "\\t" $1 }' idt-xgen-cfdna-umis.txt > umis.tsv"""
FIXED_COMMAND = """unmux synthetic_nonrandom_umi_1.fastq.gz synthetic_nonrandom_umi_2.fastq.gz \\
  --group "umi1=umis.tsv" --group "umi1::loc=0:0:8,dist=1" \\
  --group "umi2=umis.tsv" --group "umi2::loc=1:0:8,dist=1" \\
  --extract "fixed1=@umi1" --extract "fixed2=@umi2" \\
  --extract "seen1=0:0:8" --extract "seen2=1:0:8" \\
  --extract "r1=0:8:end" --extract "r2=1:8:end" \\
  --template r1 r2 \\
  --tag "RX=fixed1+fixed2" --tag "RX::qual=none" \\
  --tag "OX=seen1+seen2" \\
  --out umis.unmapped.bam"""
DUPLEX_FIXED_EXACT = FIXED_COMMAND.replace("dist=1", "dist=0").replace("umis.unmapped.bam", "umis.exact.bam")
DUPLEX_FIXED = {
    "prep": f"{UMI_PREP}\nhead -n 3 umis.tsv",
    "inputs": "gzip -dc synthetic_nonrandom_umi_1.fastq.gz | head -n 4 | cut -c 1-60",
    "fixed": FIXED_COMMAND,
    "exact": DUPLEX_FIXED_EXACT,
    "counted": "samtools view -c -f 64 umis.exact.bam\nsamtools view -c -f 64 -d RX umis.exact.bam",
    "wide": {"command": FIXED_COMMAND.replace("dist=1", "dist=2"), "ok": [1]},
    "aside": DUPLEX_FIXED_EXACT.replace("dist=0", "dist=0,minFindsPerGroup=1").replace(
        "--out umis.exact.bam", "--out umis.kept.bam \\\n  --unassigned umis.aside.bam"
    ),
}
DUPLEX_FIXED_TEMPLATES = {"view": "samtools view umis.unmapped.bam | awk '$1 == \"{{pick}}\"' | cut -f 1,2,12-"}


rule duplex_fixed_fetch:
    input:
        downloads(DUPLEX_FIXED_INPUTS),
    output:
        [work("duplex_fixed", path) for _, _, path in DUPLEX_FIXED_INPUTS if path],
    log:
        logfile("duplex_fixed.fetch"),
    params:
        inputs=DUPLEX_FIXED_INPUTS,
        recipe="duplex_fixed",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule duplex_fixed_run:
    input:
        rules.duplex_fixed_fetch.output,
        TOOLS,
    output:
        run=work("duplex_fixed", "run.json"),
    threads: 8
    log:
        logfile("duplex_fixed.run"),
    params:
        commands=DUPLEX_FIXED,
        templates=DUPLEX_FIXED_TEMPLATES,
        used=DUPLEX_FIXED_INPUTS,
        recipe="duplex_fixed",
    run:
        cookbook.record_rule(CACHE, output, params, log)
