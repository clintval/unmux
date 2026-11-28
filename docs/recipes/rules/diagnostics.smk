QC_TAG_FILES = ["R1.fastq.gz", "R2.fastq.gz", "I1.fastq.gz", "I2.fastq.gz", "samples.tsv", "i7.tsv", "i5.tsv", "phix.tsv"]
QC_TAG = {
    "qc": BINS_COMMAND.replace(
        "  --metrics-per-sample metrics.per_sample.tsv \\\n  --metrics-summary metrics.summary.tsv \\\n"
        "  --metrics-unassigned metrics.unassigned.tsv",
        "  --qc-tag",
    ),
}
QC_TAG_TEMPLATES = {
    "show": (
        "gzip -dc fastq/{{sample}}.R1.fastq.gz | grep '^@{{fixed}} '\n"
        "gzip -dc fastq/Undetermined.0.fastq.gz | grep '^@{{other}} '\n"
        "gzip -dc removed/PhiX.0.fastq.gz | grep '^@{{phix}} '\n"
        "gzip -dc removed/PhiX.0.fastq.gz | grep '^@{{stray}} '"
    ),
}


rule qc_tag_fetch:
    input:
        [work("bins", name) for name in QC_TAG_FILES],
    output:
        [work("qc_tag", name) for name in QC_TAG_FILES],
    log:
        logfile("qc_tag.fetch"),
    params:
        recipe="qc_tag",
        inputs=[],
        copies=lambda wildcards, input: list(zip(input, QC_TAG_FILES)),
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule qc_tag_run:
    input:
        rules.qc_tag_fetch.output,
        TOOLS,
    output:
        run=work("qc_tag", "run.json"),
    threads: 8
    log:
        logfile("qc_tag.run"),
    params:
        commands=QC_TAG,
        templates=QC_TAG_TEMPLATES,
        recipe="qc_tag",
        used=[],
    run:
        cookbook.record_rule(CACHE, output, params, log)


METRICS_FILES = ["metrics.summary.tsv", "metrics.per_sample.tsv", "metrics.unassigned.tsv"]
METRICS = {
    "summary": "cat metrics.summary.tsv",
    "per": "head -n 4 metrics.per_sample.tsv",
    "unassigned": "cat metrics.unassigned.tsv",
}


rule metrics_fetch:
    input:
        [work("bins", name) for name in METRICS_FILES],
    output:
        [work("metrics", name) for name in METRICS_FILES],
    log:
        logfile("metrics.fetch"),
    params:
        recipe="metrics",
        inputs=[],
        copies=lambda wildcards, input: list(zip(input, METRICS_FILES)),
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule metrics_run:
    input:
        rules.metrics_fetch.output,
        TOOLS,
    output:
        run=work("metrics", "run.json"),
    log:
        logfile("metrics.run"),
    params:
        commands=METRICS,
        recipe="metrics",
        templates={},
        used=[],
    run:
        cookbook.record_rule(CACHE, output, params, log)


INDEL_COMMAND = SPLITSEQ_COMMAND
WIDE_ALL_COMMAND = (
    SPLITSEQ_COMMAND.replace("loc=1:78:86", "loc=1:77:87").replace("loc=1:48:56", "loc=1:47:57").replace("loc=1:10:18", "loc=1:9:19")
)
WIDE_NEAREST_COMMAND = WIDE_ALL_COMMAND
for r in ("round1", "round2", "round3"):
    INDEL_COMMAND = INDEL_COMMAND.replace(f'"{r}::dist=1"', f'"{r}::dist=1:1:2"')
    WIDE_ALL_COMMAND = WIDE_ALL_COMMAND.replace(f'"{r}::dist=1"', f'"{r}::dist=1:1:2,mode=all"')
    WIDE_NEAREST_COMMAND = WIDE_NEAREST_COMMAND.replace(f'"{r}::dist=1"', f'"{r}::dist=1:1:2,mode=nearest"')

COLLISION_INPUTS = [
    *[(key, name, f"tenx/{path}") for key, name, path in TENX_INPUTS],
    *[(key, name, f"splitseq/{path}") for key, name, path in splitseq.INPUTS],
]
COLLISION_TEMPLATES = {"tenx": TENX_COMMAND.replace("dist=1,mode=nearest", "dist=1"), "indel": INDEL_COMMAND}


rule collision_fetch:
    input:
        downloads(COLLISION_INPUTS),
        cells=work("tenx", "cells.tsv"),
    output:
        [work("collision", path) for _, _, path in COLLISION_INPUTS],
        work("collision", "tenx/cells.tsv"),
        [work("collision", f"splitseq/{path}") for path in splitseq.STAGED[len(splitseq.INPUTS) :]],
    log:
        logfile("collision.fetch"),
    params:
        inputs=COLLISION_INPUTS,
        recipe="collision",
        copies=lambda wildcards, input: list([(input.cells, "tenx/cells.tsv")]),
        after=[("splitseq.splitseq_inputs", "splitseq")],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule collision_run:
    input:
        rules.collision_fetch.output,
        TOOLS,
    output:
        run=work("collision", "run.json"),
    log:
        logfile("collision.run"),
    params:
        templates=COLLISION_TEMPLATES,
        used=COLLISION_INPUTS,
        recipe="collision",
        commands={},
    run:
        cookbook.record_rule(CACHE, output, params, log)


WIDE_INPUTS = [*splitseq.INPUTS, *[(key, name, f"nearest/{path}") for key, name, path in splitseq.INPUTS]]
WIDE = {
    "wide_all": WIDE_ALL_COMMAND,
    "wide_nearest": {"command": WIDE_NEAREST_COMMAND, "cwd": "nearest"},
    "count": {
        "command": "samtools view -c splitseq.unmapped.bam\nsamtools view -c -d CB splitseq.unmapped.bam",
        "cwd": "nearest",
    },
    "qc": {
        "command": WIDE_NEAREST_COMMAND.replace('--out "splitseq.unmapped.bam"', '--qc-tag \\\n  --out "qc.unmapped.bam"'),
        "cwd": "nearest",
    },
}


rule wide_fetch:
    input:
        downloads(WIDE_INPUTS),
    output:
        [work("wide", path) for path in splitseq.STAGED],
        [work("wide", f"nearest/{path}") for path in splitseq.STAGED],
    log:
        logfile("wide.fetch"),
    params:
        inputs=WIDE_INPUTS,
        recipe="wide",
        copies=[],
        after=[("splitseq.splitseq_inputs", ""), ("splitseq.splitseq_inputs", "nearest")],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule wide_run:
    input:
        rules.wide_fetch.output,
        work("splitseq", "splitseq.unmapped.bam"),
        TOOLS,
    output:
        run=work("wide", "run.json"),
    threads: 8
    log:
        logfile("wide.run"),
    params:
        commands=WIDE,
        used=splitseq.INPUTS,
        recipe="wide",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)
