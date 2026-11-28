import nanopore

NANOPORE_INPUTS = (
    [("dorado", "barcode_kits.cpp", None)]
    + [("sra-verd6", f"{run}.fastq", None) for _, run in nanopore.VERD6]
    + [("modkit", "bc_anchored_10_reads.unmapped.bam", "mods/modbam.unmapped.bam")]
)
NANOPORE_COMMAND = """unmux reads.fastq.gz \\
  --group "front=native_barcodes.tsv" \\
  --group "front::loc=0:0:150,dist=12:12:12,mode=nearest" \\
  --group "rear=native_barcodes.tsv" \\
  --group "rear::loc=0:-150:end,dist=12:12:12,mode=nearest,bothStrands=true" \\
  --extract "insert=@front..@rear" \\
  --template insert \\
  --sample-sheet barcode_pairs.tsv \\
  --out "demux/%sample.fastq.gz" \\
  --unassigned demux/unclassified.fastq.gz"""
MODS_COMMAND = """unmux modbam.unmapped.bam --out kept.bam
unmux modbam.unmapped.bam --extract "insert=0:20:end" --template insert --out trimmed.bam
samtools view kept.bam | head -n 1 | cut -f 12- | tr '\\t' '\\n' | cut -d : -f 1 | paste -s -d ' ' -
samtools view trimmed.bam | head -n 1 | cut -f 12- | tr '\\t' '\\n' | cut -d : -f 1 | paste -s -d ' ' -"""
NANOPORE = {
    "tags": "head -n 3 native_barcodes.tsv\nhead -n 3 barcode_pairs.tsv",
    "demux": NANOPORE_COMMAND,
    "counts": 'for f in demux/NB0[1-3].fastq.gz demux/unclassified.fastq.gz; do\n  echo "$f $(gzip -dc "$f" | wc -l | awk \'{print $1 / 4}\')"\ndone',
    "tight": NANOPORE_COMMAND.replace("dist=12:12:12", "dist=8:8:8").replace("demux/", "demux8/"),
    "mods": {"command": MODS_COMMAND, "cwd": "mods"},
}


rule nanopore_fetch:
    input:
        downloads(NANOPORE_INPUTS),
    output:
        work("nanopore", "mods/modbam.unmapped.bam"),
        work("nanopore", "native_barcodes.tsv"),
        work("nanopore", "barcode_pairs.tsv"),
        work("nanopore", "reads.fastq.gz"),
    log:
        logfile("nanopore.fetch"),
    params:
        inputs=NANOPORE_INPUTS,
        recipe="nanopore",
        copies=[],
        after=[("nanopore.stage", "")],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule nanopore_run:
    input:
        rules.nanopore_fetch.output,
        TOOLS,
    output:
        run=work("nanopore", "run.json"),
    threads: 8
    log:
        logfile("nanopore.run"),
    params:
        commands=NANOPORE,
        used=NANOPORE_INPUTS,
        recipe="nanopore",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)


PACBIO_INPUTS = [
    ("pacbio-cyp2d6", name, name) for name in ("prefix.bam", "biosamples.csv", "Sequel_96_barcodes_v1.fasta")
]
PACBIO_PREP = """samtools view -h prefix.bam 2> /dev/null | samtools view -b -o hifi_reads.bam -
awk -F , 'NR == FNR { if (/^>/) id = substr($0, 2); else seq[id] = $0; next }
  FNR == 1 { print "id\\tseq" > "barcodes.tsv"; print "sample\\tstart\\tend" > "barcode_pairs.tsv"; next }
  { split($1, bc, "--")
    if (!(bc[1] in seen)) { seen[bc[1]]; print bc[1] "\\t" seq[bc[1]] "GTAAAACGACGGCCAGT" > "barcodes.tsv" }
    if (!(bc[2] in seen)) { seen[bc[2]]; print bc[2] "\\t" seq[bc[2]] "CAGGAAACAGCTATGAC" > "barcodes.tsv" }
    print $2 "\\t" bc[1] "\\t" bc[2] > "barcode_pairs.tsv"; print $2 "\\t" bc[2] "\\t" bc[1] > "barcode_pairs.tsv" }' \\
  Sequel_96_barcodes_v1.fasta biosamples.csv"""
LIMA_COMMAND = """lima --hifi-preset ASYMMETRIC --split-named -j 2 --biosample-csv biosamples.csv \\
  hifi_reads.bam Sequel_96_barcodes_v1.fasta lima.bam"""
PACBIO_COMMAND = """unmux hifi_reads.bam \\
  --group "start=barcodes.tsv" \\
  --group "start::loc=0:0:60,dist=2:2:3,mode=nearest" \\
  --group "end=barcodes.tsv" \\
  --group "end::loc=0:-60:end,dist=2:2:3,mode=nearest,bothStrands=true" \\
  --extract "insert=@start..@end" \\
  --template insert \\
  --sample-sheet barcode_pairs.tsv \\
  --out "unmux/%sample.bam" \\
  --unassigned unmux/unassigned.bam"""
PACBIO = {
    "prep": PACBIO_PREP,
    "heads": "head -n 3 barcodes.tsv\nhead -n 3 barcode_pairs.tsv",
    "version": "lima --version",
    "lima": LIMA_COMMAND,
    "demux": PACBIO_COMMAND,
    "qc": PACBIO_COMMAND.replace('--out "unmux/%sample.bam"', '--qc-tag \\\n  --out "qc/%sample.bam"').replace(
        "--unassigned unmux/unassigned.bam", "--unassigned qc/unassigned.bam"
    ),
}


rule pacbio_fetch:
    input:
        downloads(PACBIO_INPUTS),
    output:
        [work("pacbio", path) for _, _, path in PACBIO_INPUTS],
    log:
        logfile("pacbio.fetch"),
    params:
        inputs=PACBIO_INPUTS,
        recipe="pacbio",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule pacbio_run:
    input:
        rules.pacbio_fetch.output,
        TOOLS,
    output:
        run=work("pacbio", "run.json"),
    threads: 8
    log:
        logfile("pacbio.run"),
    params:
        commands=PACBIO,
        used=PACBIO_INPUTS,
        recipe="pacbio",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)


DIRECT_RNA_INPUTS = [("dorado", "fastq_with_us_and_tags.fq", "direct_rna.fastq")]
DIRECT_RNA = {
    "head": "head -n 2 direct_rna.fastq | cut -c 1-100",
    "convert": "unmux direct_rna.fastq --out direct_rna.bam",
    "view": "samtools view direct_rna.bam | cut -f 10 | cut -c 1-60\nsamtools view direct_rna.bam | cut -f 1,2,12-",
    "text": (
        "unmux direct_rna.fastq --out direct_rna.sam\n"
        "grep -v '^@' direct_rna.sam | cut -f 10 | cut -c 1-40\n"
        "unmux direct_rna.fastq --out copy.fastq\n"
        "head -n 2 copy.fastq | tail -n 1 | cut -c 1-40"
    ),
    "import": "samtools import -0 direct_rna.fastq -o samtools.bam",
}
DIRECT_RNA_TEMPLATES = {
    "probe": (
        'unmux direct_rna.fastq --group "x={dna={{probe}}}" --group "x::bothStrands=true" '
        '--sample-from-group x --out "probe.%sample.fastq"'
    ),
}


rule direct_rna_fetch:
    input:
        downloads(DIRECT_RNA_INPUTS),
    output:
        work("direct_rna", "direct_rna.fastq"),
    log:
        logfile("direct_rna.fetch"),
    params:
        inputs=DIRECT_RNA_INPUTS,
        recipe="direct_rna",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule direct_rna_run:
    input:
        rules.direct_rna_fetch.output,
        TOOLS,
    output:
        run=work("direct_rna", "run.json"),
    threads: 8
    log:
        logfile("direct_rna.run"),
    params:
        commands=DIRECT_RNA,
        templates=DIRECT_RNA_TEMPLATES,
        used=DIRECT_RNA_INPUTS,
        recipe="direct_rna",
    run:
        cookbook.record_rule(CACHE, output, params, log)


COMMENT_TAGS_INPUTS = [("modkit", "bc_anchored_10_reads.unmapped.bam", "modbam.unmapped.bam")]
COMMENT_PREP = "samtools fastq -T '*' modbam.unmapped.bam > modbam.fastq"
COMMENT_TAGS = {
    "prep": f"{COMMENT_PREP}\nhead -n 1 modbam.fastq | tr '\\t' '\\n' | cut -c 1-60",
    "convert": "unmux modbam.fastq --out modbam.roundtrip.bam",
    "view": "samtools view modbam.roundtrip.bam | head -n 1 | cut -f 1,2,12- | tr '\\t' '\\n' | cut -c 1-60",
    "direct": "unmux modbam.unmapped.bam --out direct.bam",
}


rule comment_tags_fetch:
    input:
        downloads(COMMENT_TAGS_INPUTS),
    output:
        work("comment_tags", "modbam.unmapped.bam"),
    log:
        logfile("comment_tags.fetch"),
    params:
        inputs=COMMENT_TAGS_INPUTS,
        recipe="comment_tags",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule comment_tags_run:
    input:
        rules.comment_tags_fetch.output,
        TOOLS,
    output:
        run=work("comment_tags", "run.json"),
    threads: 8
    log:
        logfile("comment_tags.run"),
    params:
        commands=COMMENT_TAGS,
        used=COMMENT_TAGS_INPUTS,
        recipe="comment_tags",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)
