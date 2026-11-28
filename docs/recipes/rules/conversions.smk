FASTQ_TO_UBAM_INPUTS = [
    ("sra-SRR390728", "SRR390728_1.fastq", "SRR390728_1.fastq.gz"),
    ("sra-SRR390728", "SRR390728_2.fastq", "SRR390728_2.fastq.gz"),
]
FASTQ_TO_UBAM = {
    "convert": "unmux SRR390728_1.fastq.gz SRR390728_2.fastq.gz \\\n  --rg-tag PL=ILLUMINA \\\n  --out SRR390728.unmapped.bam",
    "view": "samtools view --no-PG -h SRR390728.unmapped.bam | head -n 5",
    "compare": (
        "samtools import -1 SRR390728_1.fastq.gz -2 SRR390728_2.fastq.gz -o samtools.unmapped.bam\n"
        "diff <(samtools view samtools.unmapped.bam | cut -f 1-4,6-11) \\\n"
        "  <(samtools view SRR390728.unmapped.bam | cut -f 1-4,6-11) && echo identical"
    ),
    "mapq": "samtools view samtools.unmapped.bam | cut -f 5 | sort -u",
}


rule fastq_to_ubam_fetch:
    input:
        downloads(FASTQ_TO_UBAM_INPUTS),
    output:
        [work("fastq_to_ubam", path) for _, _, path in FASTQ_TO_UBAM_INPUTS],
    log:
        logfile("fastq_to_ubam.fetch"),
    params:
        inputs=FASTQ_TO_UBAM_INPUTS,
        recipe="fastq_to_ubam",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule fastq_to_ubam_run:
    input:
        rules.fastq_to_ubam_fetch.output,
        TOOLS,
    output:
        run=work("fastq_to_ubam", "run.json"),
        bam=work("fastq_to_ubam", "SRR390728.unmapped.bam"),
    threads: 8
    log:
        logfile("fastq_to_ubam.run"),
    params:
        commands=FASTQ_TO_UBAM,
        used=FASTQ_TO_UBAM_INPUTS,
        recipe="fastq_to_ubam",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)


UBAM_TO_FASTQ = {
    "convert": 'unmux SRR390728.unmapped.bam --out "SRR390728.R%ordinal.fastq.gz"',
    "listing": "ls SRR390728.R*.fastq.gz",
    "view": "gzip -dc SRR390728.R1.fastq.gz | head -n 4",
    "compare": (
        "samtools fastq -T RG -1 samtools.R1.fastq.gz -2 samtools.R2.fastq.gz SRR390728.unmapped.bam\n"
        "for r in 1 2; do\n"
        "  diff <(gzip -dc samtools.R$r.fastq.gz | tr '\\t' ' ') <(gzip -dc SRR390728.R$r.fastq.gz) && echo \"R$r identical\"\n"
        "done"
    ),
}


rule ubam_to_fastq_fetch:
    input:
        bam=rules.fastq_to_ubam_run.output.bam,
        reads=rules.fastq_to_ubam_fetch.output,
    output:
        work("ubam_to_fastq", "SRR390728.unmapped.bam"),
    log:
        logfile("ubam_to_fastq.fetch"),
    params:
        recipe="ubam_to_fastq",
        inputs=[],
        copies=lambda wildcards, input: list([(input.bam, "SRR390728.unmapped.bam")]),
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule ubam_to_fastq_run:
    input:
        rules.ubam_to_fastq_fetch.output,
        TOOLS,
    output:
        run=work("ubam_to_fastq", "run.json"),
    threads: 8
    log:
        logfile("ubam_to_fastq.run"),
    params:
        commands=UBAM_TO_FASTQ,
        recipe="ubam_to_fastq",
        templates={},
        used=[],
    run:
        cookbook.record_rule(CACHE, output, params, log)


PASSTHROUGH_INPUTS = [("modkit", "fwd_rev_modbase_records.bam", "modbam.aligned.bam")]
PASSTHROUGH_BEFORE = (
    "samtools view --no-PG -H modbam.aligned.bam | grep -v '^@SQ' | cut -c 1-100\n"
    "samtools view modbam.aligned.bam | cut -f 1-6\n"
    "samtools view modbam.aligned.bam | head -n 1 | cut -f 12- | tr '\\t' '\\n' | cut -d : -f 1 | paste -s -d ' ' -"
)
PASSTHROUGH = {
    "convert": "unmux modbam.aligned.bam --out modbam.unmapped.bam",
    "sq": "samtools view --no-PG -H modbam.aligned.bam | grep -c '^@SQ'",
    "before": PASSTHROUGH_BEFORE,
    "after": PASSTHROUGH_BEFORE.replace("modbam.aligned.bam", "modbam.unmapped.bam").replace(
        " | grep -v '^@SQ' | cut -c 1-100", ""
    ),
    "compare": (
        "samtools reset -o reset.bam modbam.aligned.bam\n"
        "diff <(samtools view reset.bam | cut -f 1-4,6-11) <(samtools view modbam.unmapped.bam | cut -f 1-4,6-11) && echo identical"
    ),
}


rule passthrough_fetch:
    input:
        downloads(PASSTHROUGH_INPUTS),
    output:
        work("passthrough", "modbam.aligned.bam"),
    log:
        logfile("passthrough.fetch"),
    params:
        inputs=PASSTHROUGH_INPUTS,
        recipe="passthrough",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule passthrough_run:
    input:
        rules.passthrough_fetch.output,
        TOOLS,
    output:
        run=work("passthrough", "run.json"),
    threads: 8
    log:
        logfile("passthrough.run"),
    params:
        commands=PASSTHROUGH,
        used=PASSTHROUGH_INPUTS,
        recipe="passthrough",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)


SPLIT_RG_INPUTS = [("samtools", "split.sam", "split.sam")]
SPLIT_RG = {
    "split": 'unmux split.sam \\\n  --group "rg=@RG" \\\n  --sample-from-group rg \\\n  --out "by_rg/%sample.bam" \\\n  --unassigned by_rg/unknown.bam',
    "counts": 'for f in by_rg/*.bam; do echo "$f $(samtools view -c "$f")"; done',
    "header": "samtools view --no-PG -H by_rg/grp1.bam | grep -v '^@CO'",
    "compare": (
        'samtools split --no-PG -u unknown.samtools.bam -f "%!.samtools.bam" split.sam\n'
        "for rg in grp1 grp2 unknown; do\n"
        '  diff <(samtools view by_rg/$rg.bam) <(samtools view $rg.samtools.bam) && echo "$rg identical"\n'
        "done"
    ),
    "hdiff": {
        "command": "diff <(samtools view --no-PG -H by_rg/grp1.bam) <(samtools view --no-PG -H grp1.samtools.bam)",
        "ok": [0, 1],
    },
    "total": "grep -vc '^@' split.sam",
}


rule split_rg_fetch:
    input:
        downloads(SPLIT_RG_INPUTS),
    output:
        work("split_rg", "split.sam"),
    log:
        logfile("split_rg.fetch"),
    params:
        inputs=SPLIT_RG_INPUTS,
        recipe="split_rg",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule split_rg_run:
    input:
        rules.split_rg_fetch.output,
        TOOLS,
    output:
        run=work("split_rg", "run.json"),
    threads: 8
    log:
        logfile("split_rg.run"),
    params:
        commands=SPLIT_RG,
        used=SPLIT_RG_INPUTS,
        recipe="split_rg",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)


CHANGESEQ_INPUTS = [
    *[("changeseq", f"{s}.{r}.fastq.gz", f"{s}.{r}.fastq.gz") for s in ("TEST", "TEST_control") for r in ("r1", "r2")],
    ("changeseq", "mergeReads.py", None),
]
CHANGESEQ = {
    "inputs": "gzip -dc TEST.r1.fastq.gz | head -n 2 | cut -c 1-70\ngzip -dc TEST.r2.fastq.gz | head -n 2 | cut -c 1-70",
    "merge": """unmux TEST.r1.fastq.gz TEST.r2.fastq.gz \\
  --extract "r1=0:0:end" --extract "r2=1:0:end" \\
  --template "~r1+r2" \\
  --out TEST_merged.fastq.gz
unmux TEST_control.r1.fastq.gz TEST_control.r2.fastq.gz \\
  --extract "r1=0:0:end" --extract "r2=1:0:end" \\
  --template "~r1+r2" \\
  --out control_TEST_merged.fastq.gz""",
    "view": "gzip -dc TEST_merged.fastq.gz | head -n 2 | cut -c 1-70",
}


rule changeseq_fetch:
    input:
        downloads(CHANGESEQ_INPUTS),
    output:
        [work("changeseq", path) for _, _, path in CHANGESEQ_INPUTS if path],
    log:
        logfile("changeseq.fetch"),
    params:
        inputs=CHANGESEQ_INPUTS,
        recipe="changeseq",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule changeseq_run:
    input:
        rules.changeseq_fetch.output,
        TOOLS,
    output:
        run=work("changeseq", "run.json"),
    threads: 8
    log:
        logfile("changeseq.run"),
    params:
        commands=CHANGESEQ,
        used=CHANGESEQ_INPUTS,
        recipe="changeseq",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)
