import splitseq

TENX_INPUTS = [
    ("nfcore", name, name) for name in ("S10_L001_R1_001.fastq.gz", "S10_L001_R2_001.fastq.gz", "10xv2_whitelist.txt")
]
CELLS_COMMAND = """gzip -dc S10_L001_R1_001.fastq.gz \\
  | awk 'NR % 4 == 2 { print substr($0, 1, 16) }' | sort -u > seen.txt
awk 'NR == FNR { keep[$1]; next }
  { for (i = 1; i <= 16; i++) for (j = 1; j <= 4; j++) {
      s = substr($1, 1, i - 1) substr("ACGT", j, 1) substr($1, i + 1)
      if (s in keep) print s } }' 10xv2_whitelist.txt seen.txt \\
  | sort -u | awk 'BEGIN { print "id\\tseq" } { print $1 "\\t" $1 }' > cells.tsv"""
TENX_COMMAND = """unmux S10_L001_R1_001.fastq.gz S10_L001_R2_001.fastq.gz \\
  --pool S10 \\
  --group "cb=cells.tsv" \\
  --group "cb::loc=0:0:16,dist=1,mode=nearest" \\
  --extract "cell=@cb" \\
  --extract "barcode=0:0:16" \\
  --extract "umi=0:16:26" \\
  --extract "cdna=1:0:end" \\
  --template cdna \\
  --tag "CB=cell" \\
  --tag "CR=barcode" \\
  --tag "UR=umi" --tag "UR::qual=UY" \\
  --tag "UB=umi" \\
  --out S10.unmapped.bam"""
TENX = {"cells": CELLS_COMMAND, "tag": TENX_COMMAND}
TENX_TEMPLATES = {"view": "samtools view S10.unmapped.bam | grep -E '^({{fixed}}|{{none}})[[:space:]]' | cut -f 1,2,12-"}


rule tenx_fetch:
    input:
        downloads(TENX_INPUTS),
    output:
        [work("tenx", path) for _, _, path in TENX_INPUTS],
    log:
        logfile("tenx.fetch"),
    params:
        inputs=TENX_INPUTS,
        recipe="tenx",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule tenx_run:
    input:
        rules.tenx_fetch.output,
        TOOLS,
    output:
        run=work("tenx", "run.json"),
        cells=work("tenx", "cells.tsv"),
    threads: 8
    log:
        logfile("tenx.run"),
    params:
        commands=TENX,
        templates=TENX_TEMPLATES,
        used=TENX_INPUTS,
        recipe="tenx",
    run:
        cookbook.record_rule(CACHE, output, params, log)


SPLITSEQ_INPUTS = splitseq.INPUTS + [
    ("splitseq", "tools.py", None),
    ("splitcode", "tutorials_splitseq.rst", None),
    ("splitcode", "r1_R.txt", None),
    ("splitcode", "r1_T.txt", None),
    ("splitcode", "r2_r3.txt", None),
]
SPLITSEQ_COMMAND = """unmux R1.fastq.gz R2.fastq.gz \\
  --group "round1=round1.tags.txt" \\
  --group "round1::loc=1:78:86" \\
  --group "round1::dist=1" \\
  --group "round1::minFindsPerGroup=1" \\
  --group "round2=round2_3.tags.txt" \\
  --group "round2::loc=1:48:56" \\
  --group "round2::dist=1" \\
  --group "round2::minFindsPerGroup=1" \\
  --group "round3=round2_3.tags.txt" \\
  --group "round3::loc=1:10:18" \\
  --group "round3::dist=1" \\
  --group "round3::minFindsPerGroup=1" \\
  --extract "cdna=0:0:end" \\
  --extract "umi=1:0:10" \\
  --extract "bc1=@round1" \\
  --extract "bc2=@round2" \\
  --extract "bc3=@round3" \\
  --template "cdna" \\
  --tag "CB=bc1+bc2+bc3" \\
  --tag "CB::qual=CY" \\
  --tag "CB::sep=-" \\
  --tag "CB::qual-sep=%20" \\
  --tag "RX=umi" \\
  --tag "RX::qual=QX" \\
  --out "splitseq.unmapped.bam\""""
SPLITCODE_COMMAND = """splitcode -c config.txt --nFastqs=2 --select=0 \\
--gzip -o output_R1.fastq.gz,unused_R2.fastq.gz R1.fastq.gz R2.fastq.gz"""
SPLITCODE_LISTS = """awk -F , '$1 < 48 { print $2 }' bc_8nt_v2.csv > r1_T.txt
awk -F , '$1 >= 48 { print $2 }' bc_8nt_v2.csv > r1_R.txt"""
SPLITSEQ = {
    "tags": f"{splitseq.SPLITSEQ_TAGS}\nhead -n 3 round1.tags.txt",
    "demux": SPLITSEQ_COMMAND,
    "view": "samtools view splitseq.unmapped.bam | head -n 2 | cut -f 1,2,12-",
    "splitcode_version": "splitcode --version",
    "as_written": {
        "command": SPLITCODE_COMMAND.replace(",unused_R2.fastq.gz", ""),
        "ok": [1],
        "cwd": "splitcode-tutorial",
    },
    "tutorial": {"command": SPLITCODE_COMMAND, "cwd": "splitcode-tutorial"},
    "lists": {"command": SPLITCODE_LISTS, "cwd": "splitcode-library"},
    "library": {"command": SPLITCODE_COMMAND, "cwd": "splitcode-library"},
}


rule splitseq_fetch:
    input:
        downloads(SPLITSEQ_INPUTS),
    output:
        [work("splitseq", path) for path in [*splitseq.STAGED, *splitseq.SPLITCODE_STAGED]],
    log:
        logfile("splitseq.fetch"),
    params:
        inputs=SPLITSEQ_INPUTS,
        recipe="splitseq",
        copies=[],
        after=[("splitseq.splitseq_inputs", ""), ("splitseq.splitcode_inputs", "")],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule splitseq_run:
    input:
        rules.splitseq_fetch.output,
        TOOLS,
    output:
        run=work("splitseq", "run.json"),
        bam=work("splitseq", "splitseq.unmapped.bam"),
    threads: 8
    log:
        logfile("splitseq.run"),
    params:
        commands=SPLITSEQ,
        used=SPLITSEQ_INPUTS,
        recipe="splitseq",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)


MATCH_INPUTS = splitseq.INPUTS
MATCH_COMMAND = """unmux R1.fastq.gz R2.fastq.gz \\
  --pool C2C12_short_1k \\
  --group "round1=round1.tags.txt" --group "round1::loc=1:78:86,dist=1" \\
  --group "round2=round2_3.tags.txt" --group "round2::loc=1:48:56,dist=1" \\
  --group "round3=round2_3.tags.txt" --group "round3::loc=1:10:18,dist=1" \\
  --group "cells=cells.tsv" \\
  --group "cells::match=bc1+bc2+bc3,minFindsPerGroup=1" \\
  --extract "bc1=@round1" --extract "bc2=@round2" --extract "bc3=@round3" \\
  --extract "name=@cells::cell" \\
  --extract "cdna=0:0:end" --extract "umi=1:0:10" \\
  --template cdna \\
  --tag "CB=name" \\
  --tag "RX=umi" \\
  --out cells.unmapped.bam \\
  --unassigned other.unmapped.bam"""
MATCH = {
    "cells": f"{splitseq.SPLITSEQ_CELLS}\nhead -n 3 cells.tsv",
    "match": MATCH_COMMAND,
    "view": "samtools view cells.unmapped.bam | head -n 2 | cut -f 1,2,12-",
}


rule match_fetch:
    input:
        downloads(MATCH_INPUTS),
    output:
        [work("match", path) for path in splitseq.STAGED],
    log:
        logfile("match.fetch"),
    params:
        inputs=MATCH_INPUTS,
        recipe="match",
        copies=[],
        after=[("splitseq.splitseq_inputs", "")],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule match_run:
    input:
        rules.match_fetch.output,
        TOOLS,
    output:
        run=work("match", "run.json"),
    threads: 8
    log:
        logfile("match.run"),
    params:
        commands=MATCH,
        used=MATCH_INPUTS,
        recipe="match",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)


CITE_NAMES = (
    "subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz",
    "subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz",
    "5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv",
)
CITE_INPUTS = [("nfcore-modules", name, name) for name in CITE_NAMES]
CITE_PREP = """awk -F , 'NR == 1 { print "id\\tseq"; next } { print $1 "\\t" $5 }' \\
  5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv > antibodies.tsv"""
CITE_COMMAND = """unmux subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz \\
  subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz \\
  --pool antiCMV_AB \\
  --group "ab=antibodies.tsv" \\
  --group "ab::loc=1:10:25,dist=1,minFindsPerGroup=1" \\
  --extract "cell=0:0:16" \\
  --extract "umi=0:16:end" \\
  --extract "barcode=@ab" \\
  --extract "antibody=@ab::id" \\
  --extract "r2=1:0:end" \\
  --template r2 \\
  --tag "CR=cell" \\
  --tag "UR=umi" --tag "UR::qual=UY" \\
  --tag "fb=barcode" \\
  --tag "fx=antibody" \\
  --out adt.unmapped.bam"""
CITE = {
    "prep": f"{CITE_PREP}\nhead -n 3 5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv\nhead -n 3 antibodies.tsv",
    "input": f"gzip -dc {CITE_NAMES[1]} | head -n 2",
    "tag": CITE_COMMAND,
    "view": "samtools view adt.unmapped.bam | head -n 2 | cut -f 1,2,12-",
    "count": (
        "samtools view adt.unmapped.bam | grep -o 'fx:Z:[^[:space:]]*' \\\n"
        "  | sort | uniq -c | awk '{ print $1, $2 }' | sort -k1,1nr -k2,2 | head -n 5"
    ),
}


rule cite_fetch:
    input:
        downloads(CITE_INPUTS),
    output:
        [work("cite", path) for _, _, path in CITE_INPUTS],
    log:
        logfile("cite.fetch"),
    params:
        inputs=CITE_INPUTS,
        recipe="cite",
        copies=[],
        after=[],
    run:
        cookbook.stage_rule(CACHE, params, log)


rule cite_run:
    input:
        rules.cite_fetch.output,
        TOOLS,
    output:
        run=work("cite", "run.json"),
    threads: 8
    log:
        logfile("cite.run"),
    params:
        commands=CITE,
        used=CITE_INPUTS,
        recipe="cite",
        templates={},
    run:
        cookbook.record_rule(CACHE, output, params, log)
