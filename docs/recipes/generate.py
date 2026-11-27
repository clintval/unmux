#!/usr/bin/env python3
"""Regenerate docs/recipes.md, the unmux cookbook.

Fetch the pinned public inputs into a cache directory, run every recipe and its checks with
an unmux build, and write the page from what the runs print and write. The Makefile beside
this script runs it; by hand:

    cargo build --release
    python3 -I docs/recipes/generate.py fetch
    python3 -I docs/recipes/generate.py run
    python3 -I docs/recipes/generate.py page
    python3 -I docs/recipes/generate.py check

The options --samtools, --fg-sra and --lima give the commands that run those tools.
"""

import argparse
import collections
import csv
import difflib
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAGE = ROOT / "docs" / "recipes.md"
USER_AGENT = "unmux-recipes/1.0"

GITHUB = "https://raw.githubusercontent.com"
SAMTOOLS_SHA = "7162950d0988480faeddeb49f03520737e13a9c4"
NFCORE_SHA = "d934d6e8367fe2626184496b1889671cf2b02dab"
FASTQUORUM_SHA = "2ca9f8176f1f40ec1124c5182c58524766c95411"
DORADO_SHA = "651980f48c21eb93a5c46c65e1f004136f02c07a"
MODKIT_SHA = "5cecc3fb3a9336068d9e3c68d5c08d678153dd2c"
SPLITSEQ_SHA = "c3923ea44a94a1fc0effe5ff8e2ca479c0bff934"

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class Source:
    """A pinned public input: a whole file, the first records of a gzipped FASTQ, or an SRA spot range."""

    def __init__(self, key, name, url, sha256, label, records=None, sra=None):
        self.key = key
        self.name = name
        self.url = url
        self.sha256 = sha256
        self.label = label
        self.records = records
        self.sra = sra


def sra_slice(key, accession, first, last, outputs, label, extra=""):
    """Sources for the files one `fg-sra fastq` call writes for a spot range of a run."""
    flags = " ".join(f"--{flag} {name}" for flag, name in outputs)
    command = f"fg-sra fastq --min-spot-id {first} --max-spot-id {last} {flags} {extra}".strip() + f" {accession}"
    url = f"https://www.ncbi.nlm.nih.gov/sra/{accession}"
    return [
        Source(key, name, url, "", f"SRA run {accession}, spots {first:,} to {last:,}", sra=command)
        for _, name in outputs
    ]


UNPINNED_OK = False

PINS = {
    "sra-SRR390728/SRR390728_1.fastq": "93f6ea4eae9361fa6109c62caed201f8d462d11a2af3fdb2307eb3a597e1aabd",
    "sra-SRR390728/SRR390728_2.fastq": "479caf671e0766784f78ff7aedc7578ca6a1cee66f983a7c6fdcc1d8df99e134",
    "uditas/Undetermined_S0_L001_R1_001.fastq": "cc4c86e51b466a7c8f6885f3a0c6b7c832f9a5d80cdf5840d6ebe689f3abbfa1",
    "uditas/Undetermined_S0_L001_R2_001.fastq": "91bfe77bcf39b4e248a6e67814d01e3c9eb52af4e89bfe8fc08694b2d1b9d283",
    "uditas/Undetermined_S0_L001_I1_001.fastq": "d45f9944f59c394b5804d0706d17d3d84ee04b5ff58c60c1795a223d03d53317",
    "uditas/Undetermined_S0_L001_I2_001.fastq": "cba62a532e8a8eb3cbfa7d6faeaf032ab4edcc627443fe5eeda7143f541f1e0e",
    "uditas/sample_info.csv": "11807fb316b2e86a2b59bff67b765eaf4778ce67c38501f910856b25a9c892cc",
    "uditas-fig2c/Undetermined_S0_L001_R1_001.fastq": "8183bf5ab0a273466b4c17a700a34f48ba50418d7911c93abce6c10cbf392019",
    "uditas-fig2c/Undetermined_S0_L001_R2_001.fastq": "99ab74d3e3afc5c8d7851747981fdf52b3c600eced7fce8cc2c7cce5ca609d57",
    "uditas-fig2c/Undetermined_S0_L001_I1_001.fastq": "9b10600a26d8dd53abb77265b1e310e0a3987a0ff500051d8d496a106e18f20c",
    "uditas-fig2c/Undetermined_S0_L001_I2_001.fastq": "403d5b9c34908ca2544c46cf9ef9a307ff36f5555f3cf10e16a3b40209cb0c15",
    "changeseq/TEST.r1.fastq.gz": "a4843b9ff978958de6863718193c8380218e551717e659f3856b891539465e72",
    "changeseq/TEST.r2.fastq.gz": "ee0117a9f829f86aa942efc521a2e97d507fc84e51e6b95ffa483ccafa55ca1a",
    "changeseq/TEST_control.r1.fastq.gz": "c66c3a4f6791da9cadc6fe81e36337db5b6a8f9bbbfcfdb6a49d4b48f8eaccef",
    "changeseq/TEST_control.r2.fastq.gz": "d9b90ec361f8a31b34cfa4b9fbc8e40e139e4784ac57fb9369522731d6f70a10",
    "changeseq/mergeReads.py": "61854937da7e2a0be34f323298802a37ca99b0f82de4b0d409e2d5244a6c9467",
    "sra-verd6/SRR14560737.fastq": "3a692351b5d50e7387d468df2869eb19a701906c1013a94179ec64a80cacf6f1",
    "sra-verd6/SRR14560738.fastq": "54729bb12872bef5f8f2e5a96669e45752fd7d71a5e8586b640b54936a11d870",
    "sra-verd6/SRR14560739.fastq": "54cf8134c59224cc98746faa52febf2ef9e045a51faf7eb653f23dd4ebae7ae8",
    "sra-verd6/SRR14560652.fastq": "dcc6dae37460689c7f051705bbc747ccd7e8e1cca92cc4a3ac3527ff9cac8856",
    "sra-verd6/SRR14560583.fastq": "bfa90841686f0a6945f7552e7db7458cbcffb7d5a75ad41ef269a198ea57c105",
    "sra-verd6/SRR14560584.fastq": "ecbdb88f15476b226059aea8d8ee81a0b3ba94afe9e610a2e1883e01c41d1994",
    "sra-verd6/SRR14560585.fastq": "4fa1317adf546e55c818bb11a71af9e324f48e7128f88d2a1244fd3f1d36fc80",
    "sra-verd6/SRR14560586.fastq": "cc42fa60cda9831ef4c3c1b4e60b53a1348a500e6919340e64c8f2282198e268",
    "sra-verd6/SRR14560587.fastq": "d2ce4814a9cdeacbb2828962e4f57060a3cd6f7ecb70d5c409b99346502b18bd",
    "sra-verd6/SRR14560588.fastq": "01fc16f31974a478f0fa75b9002ab4f8d6edd95492c6d6ed8271bef58f05bda1",
    "sra-verd6/SRR14560763.fastq": "3150aa81144c401c18b850808f0d3ab43a1648a18c6e6257b7679cef315f03d6",
    "sra-verd6/SRR14560764.fastq": "f03ecee43a0ae83583034af37abf189df639f30682b99808a2f7540d32192826",
    "sra-verd6/SRR14560765.fastq": "503f7cca77c1934756fddb88e9df66ed60621b5c2da59af59c4adf7e0df24fa0",
    "sra-verd6/SRR14560523.fastq": "d2d13bfa942b7833add772dfe81719bd3ad90ba32e2af9b30a3c24a25dbeb542",
    "sra-verd6/SRR14560525.fastq": "66a9d6ff51b98a9143a2d7948e086cead3f4d24861c4b553a00e1573b0fe50c9",
    "sra-verd6/SRR14560631.fastq": "1463b94bc77bedb1c8e6c14c7f699abefef865f1a0ade5dcb10811001832ff07",
    "sra-verd6/SRR14560632.fastq": "87802c2392aceea734eebc1a7d589aa927508843f98336c436b1c178f6aae6d7",
    "sra-verd6/SRR14560633.fastq": "237105db8f009af9264abb482e2a50ea716881ffcd159315666bca6df094b49b",
    "sra-verd6/SRR14560634.fastq": "6c421ff0111e2dbc51ee08f175c8bf04499db5d4a879af9e47bc53e48cabd82c",
    "sra-verd6/SRR14560635.fastq": "0f9752c8f7574cc0145f9ef82d5cbfcfa4b2977ff3d3b00e9f83a0ec147c4798",
    "sra-verd6/SRR14560769.fastq": "f5afceed5aacc5a7330bc205486985910cfc721268dbd58641c9843ec98a2015",
    "sra-verd6/SRR14560772.fastq": "ab17c2abe95ca9935858003f337fdae2e03a1528f35fc5cd786a810919b43a0b",
    "sra-verd6/SRR14560526.fastq": "3b9459789f9a6b3c380ff7bc31bdc257b03296d125798b2cb4bb818089fc3010",
    "sra-verd6/SRR14560531.fastq": "5a79a546ca2b9297e49125437552b95ddc97276af2d3c36a7179f9ff1d9efa36",
    "pacbio-cyp2d6/prefix.bam": "055d6bba9b3296e09fb731203e5a5c66b5c888b618076eab4fa9fee027efd187",
    "pacbio-cyp2d6/biosamples.csv": "62a94113d99bd7394297b15587e320cfbd9913003cb2302121ead540673ed01e",
    "pacbio-cyp2d6/Sequel_96_barcodes_v1.fasta": "98d81e93dff2c7e08ef26833ecbb066fb0e4bc3de3a60f48a7ca3705d029d57d",
    "sra-carabus/DRR021959.fastq": "0c238893ff23c25c2931e3ae0d5fa81de27d85b4090a7a6b40b69f12cbaa7e58",
    "sra-carabus/DRR021960.fastq": "5b27a3556070beab59e9ef456c8f0903d912c97d22911ca764fcf2e8222494f0",
    "sra-carabus/DRR021961.fastq": "d18849fb726000845295874fc3dfc2542387999a9e5dbf1b810a2703df288ebd",
    "sra-carabus/DRR021962.fastq": "d046c33aee0788148486d4bdf4a15e3d74d5d0742d45587efa750893988df82c",
    "sra-carabus/DRR021963.fastq": "0e6b8e96145f52c6851a411ffadba91aa85beda083e845a4c4ebef293fb35a67",
    "sra-carabus/DRR021964.fastq": "71ae671058495a2a697a117d3fa91a320eae68d6c800d0f215a5bf7efc1631ea",
    "sra-carabus/DRR021965.fastq": "4382ab3114205f3fcb262952a5e80274b8b4a2228dfecde5ab476941136e5970",
    "sra-carabus/DRR021966.fastq": "b0272673caecfe1d62db048ad3a869d155890f7f5a661b889571bc676e934811",
    "sra-carabus/DRR021967.fastq": "ddd4c041514d5f6d6f822707703fa1035a57b0177f81c69df0d1a43feeefd51f",
    "sra-carabus/DRR021968.fastq": "8cd1479c4679da842eea250a77b482c94e8028a3a336563d9e953cdd41f2e3ef",
    "sra-carabus/DRR021969.fastq": "2693e8bd774c8d5b145e37260088656f5b8ec7ca3c27b7ae9d8a46ea8fdc1c8f",
    "sra-carabus/DRR021970.fastq": "49238791a24bbe435b3b172ab1957bf3ec2ae24e9cd5ec7269bae925a30f18aa",
    "sra-carabus/DRR021971.fastq": "9dc6623d30215a82e4b2d9fe87c40f3cc5bec010bdbda4182dd014457b615b1e",
    "sra-carabus/DRR021972.fastq": "1e82b6637e399e4a558b29e1c92d17583604606bddfede5f7e2dda1c8c985cf0",
    "sra-carabus/DRR021973.fastq": "b4f0583690c8f3a0c6075aa323673f2766bb139271fd2cff2ae634579d1dd33c",
    "sra-carabus/DRR021974.fastq": "34d12541515b0a02a210677b9a86e97d58cd33326167951724b0a88b83896909",
    "sra-carabus/DRR021975.fastq": "9bc1833659db71bd8b4f1e9fe1604ea3131fe8570bb7f6aecbfa1c39db34b24e",
    "sra-carabus/DRR021976.fastq": "804616e2328f35e7929251e9c710551a95d2862326c5308cfce3c69f7bd8e9e9",
    "sra-carabus/DRR021977.fastq": "c9f6e640faed7c00dfa70291069819742bf3cc902dee2e14fded35c6abe16aaa",
    "sra-carabus/DRR021978.fastq": "c7fa6248da59579b711c6bb06a9ae711015f08f0c68fb67c9811936bd66117d1",
    "sra-carabus/DRR021979.fastq": "f8332f39eea56c76bcc5c5ab077956b537d0094d1023bdfb0e4fc129a1817f52",
    "sra-carabus/DRR021980.fastq": "63b6f04302bf01addb15430dcaf5601292d93b4e2891d672c361b274d2229061",
    "sra-carabus/DRR021981.fastq": "fec38c24f73f8cf0dad59d72a96f7788d5da86294247c6b5fbeddf15821234cf",
    "sra-carabus/DRR021982.fastq": "79ddea5945e7cf09e233fe1a5c7f69bd49156414bcec8ebee86ae6326180ae0b",
    "sra-carabus/DRR021983.fastq": "98945a94011a38079280c2b45e9e57a33dccf37dd06985e33c01837283314cf2",
    "sra-carabus/DRR021984.fastq": "83823d52123cf0c957a14f93fc9bca8140315c1e8605cf565dc7f4d172550945",
    "sra-carabus/DRR021985.fastq": "0077438bad31dca17989ce979c60b11dda4c9200ce00746ec0f4f3203bf5dcc0",
    "sra-carabus/DRR021986.fastq": "f5b542fdf1e0b80649322f38d57b217f43b707bda988cf0397a590d63120761a",
    "sra-carabus/DRR021987.fastq": "ea8aea3cd4683cceb6b5bad273862e4d9e20b6e85580399ad3780a32627de0d2",
    "sra-carabus/DRR021988.fastq": "e0985fef474c6061a3dca3d691300847e5c20ddb60213366b7ef66aec0b31a3f",
    "sra-carabus/DRR021989.fastq": "fbde07c93e7fe23bd24efedaa93dac1226637a464e5d8158b9cad6294256996f",
    "sra-carabus/DRR021990.fastq": "27b6d3d91794876f75fb7eeada811ea9d970c63eb2df94229a7ed60c5699e910",
    "sra-SRR13948565/SRR13948565_1.fastq": "fcad4d5e5173a21e9db9accdecee6454cf123abe8128a08019a2719a9479c6d5",
    "sra-SRR13948565/SRR13948565_2.fastq": "5abd93388975a6fe2ef9b18785485935c215ec35f0a53ddb0f56ce7af7529238",
    "geo/GSM5169184_C2C12_short_1k_cell_metadata.csv.gz": "5c50c7cc1fe85f51fe9f4724ceb211d8e0091c9cfd2acec4b22c9ab957d511cc",
    "nfcore-modules/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz": "12de3133f682ea61d5e9158021d59c20d17372bdbe97ccb40a23ca6320a377b4",
    "nfcore-modules/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz": "e531f406e8f8273d949ff62015d73edbb90982e1064326f6208e46d1a8431daa",
    "nfcore-modules/5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv": "167c46e917b3879b7b5416eba04a33388d624a1dce481e211a8bdfff0022aaf7",
    "sra-SRR27862880/SRR27862880_1.fastq": "9bcb902a6ba0d5fc3030658dbffef3573129c51e4c4b3d361d399b040c55222d",
    "sra-SRR27862880/SRR27862880_2.fastq": "344783d6ff278b44559abb61d1487adacb426fde41a5e2592f4019a73a267d35",
    "ncbi/NC_001422.1.fasta": "caa155c83967cc080f4f188f14f1f28c18d77448e87c3d00d3579ed348f4de2a",
    "uditas/uditas_helpers.py": "9aafa10ce759e93ed59a15686e8039b2d9ef8af95ec0e39b102315dbee0a2ec8",
    "samtools/split.sam": "f937e3eebd2056046ece910f934725c17946fb3be16c2339ab03c8f5ddcf6794",
    "modkit/fwd_rev_modbase_records.bam": "246f175ecd169f4a39b30a5130a7fe64234490c920f5110fb33a49d7849f7d9b",
    "modkit/bc_anchored_10_reads.unmapped.bam": "e1a802203410fccc9f6b4258cef3d9ca2cf93bb3dfc487f63e2acfa57678e144",
    "nfcore/S10_L001_R1_001.fastq.gz": "a318e85709d690b25c763fb3f156fc98d0ce46360115f3de3d9b1b2ecc55a1f7",
    "nfcore/S10_L001_R2_001.fastq.gz": "1514b26ce8e972ac7ba1860d591fa72e8c1fb6f716825651c2eb24b59c1a38d2",
    "nfcore/10xv2_whitelist.txt": "b0dda4b114d8fc7ea8def94fb472e6205c1831a72e3ac245ed13cc332d9f20d1",
    "fastquorum/samplesheet.tiny.csv": "c6c8603f79cb9ee317700d0f0b055a48b0652879552dfa6dded1c31233a03a26",
    "fastquorum/SRR6109255_1.fastq.gz": "1506c177f3fb7d58c193604de4d9c12cf365f9fb9345e26ec9c6fd7fb9dd407f",
    "fastquorum/SRR6109255_2.fastq.gz": "feded00ee29cb946b9a387392d8f0c2a5b63dca2d57cfc3b13365313d7d33340",
    "fastquorum/samplesheet.nonrandom.csv": "b3245e6abb9df4721a2fa21de08db92e5a2605c8c0a63a375efbbc5e83a9c87d",
    "fastquorum/synthetic_nonrandom_umi_1.fastq.gz": "6da282703b2eb7671800f2a34621ff2a8a0ce05b48deff6da09755466c54d855",
    "fastquorum/synthetic_nonrandom_umi_2.fastq.gz": "11d50af9ee3bf94882c3c9001057fad3777840f0097cbbc442f639db2d196cd3",
    "fastquorum/idt-xgen-cfdna-umis.txt": "b1ebe64e3adf7107154b0f08143750d9d3ca03b118091ccb1dc1124ed09553e4",
    "splitseq/bc_8nt_v1.csv": "76dd0395f1f450176d668619db4d7fd1a1cd502b8870c7d15411f7c53d29cff3",
    "splitseq/bc_8nt_v2.csv": "6180689a6bf2b26d5d6578039d6f7150ee3e8165ee3b0f6ee7d2045884fda5e2",
    "splitseq/tools.py": "87223418472ffc3e214b7b902c6c7751462a69c62e65bb4daef75cc093525ba4",
    "dorado/barcode_kits.cpp": "6dd726dd7a3abe8eabbdfa39095269fe57393b82b411570b65cedc6d2b78fc44",
    "dorado/fastq_with_us_and_tags.fq": "932a2a32f845b02db19e1688bbdc129611d984998267ca8ac28b8dd4b7b68877",
}

CARABUS_BARCODES = (
    "AAAAA AACCC AAGGG AATTT ACACG ACCAT ACGTA ACTGC CAACT CACAG CAGTC CATGA CCAAC CCCCA CCGGT CCTTG "
    "GAAGC GACTA GAGAT GATCG GCATT GCCGG GCGCC GCTAA TAATG TACGT TAGCA TATAC TCAGA TCCTC TCGAG TCTCT"
).split()
CARABUS = [(f"DRR0219{59 + k}", f"CR{k + 1:02d}", barcode) for k, barcode in enumerate(CARABUS_BARCODES)]
CARABUS_FETCH = (
    "for run in DRR0219{59..90}; do\n"
    "  fg-sra fastq --min-spot-id 1 --max-spot-id 250 --unpaired \"$run.fastq\" \"$run\"\n"
    "done"
)

VERD6 = [
    ("NB01", "SRR14560737"), ("NB02", "SRR14560738"), ("NB03", "SRR14560739"), ("NB04", "SRR14560652"),
    ("NB05", "SRR14560583"), ("NB06", "SRR14560584"), ("NB07", "SRR14560585"), ("NB08", "SRR14560586"),
    ("NB09", "SRR14560587"), ("NB10", "SRR14560588"), ("NB11", "SRR14560763"), ("NB12", "SRR14560764"),
    ("NB13", "SRR14560765"), ("NB14", "SRR14560523"), ("NB15", "SRR14560525"), ("NB16", "SRR14560631"),
    ("NB17", "SRR14560632"), ("NB18", "SRR14560633"), ("NB19", "SRR14560634"), ("NB20", "SRR14560635"),
    ("NB21", "SRR14560769"), ("NB22", "SRR14560772"), ("NB23", "SRR14560526"), ("NB24", "SRR14560531"),
]
VERD6_FETCH = (
    "for run in " + " ".join(run for _, run in VERD6[:6]) + " \\\n"
    + "".join("  " + " ".join(run for _, run in VERD6[k:k + 6]) + (" \\\n" if k < 18 else "; do\n") for k in (6, 12, 18))
    + "  fg-sra fastq --min-spot-id 1 --max-spot-id 200 --unpaired \"$run.fastq\" \"$run\"\ndone"
)
CYP2D6 = "https://downloads.pacbcloud.com/public/dataset/pbAmpliconAnalysis_CYP2D6/workflow_inputs"
CYP2D6_FETCH = f"curl -sS -A unmux-recipes/1.0 -r 0-9999999 -o prefix.bam {CYP2D6}/m64013e_210223_231140.hifi_reads.bam"

UDITAS = "https://raw.githubusercontent.com/editasmedicine/uditas/658ad062acf2b073dd42523704fae9c71e14bd8a"
CHANGESEQ = "https://raw.githubusercontent.com/tsailabSJ/changeseq/742427fb1ee9e016d247fa263039f1cd349a42cf"
NFCORE_MODULES = ("https://raw.githubusercontent.com/nf-core/test-datasets/e8d2f4f2af42cc757ed2448c2398e282e9990235/"
                  "data/genomics/homo_sapiens/10xgenomics/cellranger/5k_cmvpos_tcells")

SOURCES = [
    *sra_slice("sra-SRR390728", "SRR390728", 1, 2000, [("r1", "SRR390728_1.fastq"), ("r2", "SRR390728_2.fastq")], ""),
    *[
        Source("uditas", f"Undetermined_S0_L001_{r}_001.fastq", f"{UDITAS}/data/fig2b/Undetermined_S0_L001_{r}_001.fastq.gz",
               "", f"UDiTaS fig2b run, {r} (first 20,000 records)", records=20000)
        for r in ("R1", "R2", "I1", "I2")
    ],
    Source("uditas", "sample_info.csv", f"{UDITAS}/data/fig2b/sample_info.csv", "", "UDiTaS fig2b sample sheet"),
    *[
        Source("uditas-fig2c", f"Undetermined_S0_L001_{r}_001.fastq", f"{UDITAS}/data/fig2c/Undetermined_S0_L001_{r}_001.fastq.gz",
               "", f"UDiTaS fig2c run, {r} (first 5,000 records)", records=5000)
        for r in ("R1", "R2", "I1", "I2")
    ],
    *[
        Source("changeseq", f"{s}.{r}.fastq.gz", f"{CHANGESEQ}/test/data/input/{s}.{r}.fastq.gz", "",
               f"CHANGE-seq test reads {s}, {r}")
        for s in ("TEST", "TEST_control") for r in ("r1", "r2")
    ],
    Source("changeseq", "mergeReads.py", f"{CHANGESEQ}/changeseq/mergeReads.py", "", "CHANGE-seq read merging, mergeReads()"),
    *[
        Source("sra-verd6", f"{run}.fastq", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA730334", "",
               "SRA runs of BioProject PRJNA730334 batch VERD6, one per native barcode, spots 1 to 200 of each",
               sra=VERD6_FETCH)
        for _, run in VERD6
    ],
    Source("pacbio-cyp2d6", "prefix.bam", f"{CYP2D6}/m64013e_210223_231140.hifi_reads.bam", "",
           "PacBio pbAmpliconAnalysis_CYP2D6 HiFi reads before lima, first 10 MB", sra=CYP2D6_FETCH),
    Source("pacbio-cyp2d6", "biosamples.csv", f"{CYP2D6}/biosamples.csv", "", "PacBio CYP2D6 barcode pairs and samples"),
    Source("pacbio-cyp2d6", "Sequel_96_barcodes_v1.fasta",
           f"{GITHUB}/PacificBiosciences/hifi-amplicon-workflow/2964cd1fbddae04243d56914935881cca01f6d13/data/Sequel_96_barcodes_v1.fasta",
           "", "PacBio Sequel 96 barcodes v1"),
    *[
        Source("sra-carabus", f"{run}.fastq", "https://www.ncbi.nlm.nih.gov/bioproject/PRJDB2148", "",
               "SRA runs DRR021959 to DRR021990 of study DRP001067, spots 1 to 250 of each", sra=CARABUS_FETCH)
        for run, _, _ in CARABUS
    ],
    *sra_slice("sra-SRR13948565", "SRR13948565", 1000001, 1010000,
               [("r1", "SRR13948565_1.fastq"), ("r2", "SRR13948565_2.fastq")], ""),
    Source("geo", "GSM5169184_C2C12_short_1k_cell_metadata.csv.gz",
           "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5169nnn/GSM5169184/suppl/GSM5169184_C2C12_short_1k_cell_metadata.csv.gz",
           "", "GEO GSM5169184, the cells split-pipe called in SRR13948565"),
    *[
        Source("nfcore-modules", path.split("/")[-1], f"{NFCORE_MODULES}/{path}", "", label)
        for path, label in (
            ("fastqs/ab/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz", "nf-core/test-datasets (modules) 10x antibody reads, read 1"),
            ("fastqs/ab/subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz", "nf-core/test-datasets (modules) 10x antibody reads, read 2"),
            ("5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv", "nf-core/test-datasets (modules) antibody feature reference"),
        )
    ],
    *sra_slice("sra-SRR27862880", "SRR27862880", 1, 2000, [("r1", "SRR27862880_1.fastq"), ("r2", "SRR27862880_2.fastq")],
               "", extra="--defline '$sn $ri:N:0:$sg'"),
    Source("ncbi", "NC_001422.1.fasta",
           "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text",
           "", "NCBI RefSeq NC_001422.1, the PhiX174 genome"),
    Source("uditas", "uditas_helpers.py", f"{UDITAS}/uditas/uditas_helpers.py", "", "UDiTaS demultiplexing rules, demultiplex()"),
    Source(
        "samtools", "split.sam",
        f"{GITHUB}/samtools/samtools/{SAMTOOLS_SHA}/test/split/split.sam",
        "", "samtools 1.25 test file test/split/split.sam",
    ),
    Source(
        "modkit", "fwd_rev_modbase_records.bam",
        f"{GITHUB}/nanoporetech/modkit/{MODKIT_SHA}/tests/resources/fwd_rev_modbase_records.bam",
        "", "modkit test file tests/resources/fwd_rev_modbase_records.bam",
    ),
    Source(
        "modkit", "bc_anchored_10_reads.unmapped.bam",
        f"{GITHUB}/nanoporetech/modkit/{MODKIT_SHA}/tests/resources/bc_anchored_10_reads.unmapped.bam",
        "", "modkit test file tests/resources/bc_anchored_10_reads.unmapped.bam",
    ),
    Source(
        "nfcore", "S10_L001_R1_001.fastq.gz",
        f"{GITHUB}/nf-core/test-datasets/{NFCORE_SHA}/testdata/S10_L001_R1_001.fastq.gz",
        "", "nf-core/test-datasets (scrnaseq) 10x v2 read 1",
    ),
    Source(
        "nfcore", "S10_L001_R2_001.fastq.gz",
        f"{GITHUB}/nf-core/test-datasets/{NFCORE_SHA}/testdata/S10_L001_R2_001.fastq.gz",
        "", "nf-core/test-datasets (scrnaseq) 10x v2 read 2",
    ),
    Source(
        "nfcore", "10xv2_whitelist.txt",
        f"{GITHUB}/nf-core/test-datasets/{NFCORE_SHA}/reference/kallisto/10xv2_whitelist.txt",
        "", "nf-core/test-datasets (scrnaseq) 10x v2 barcode keeplist",
    ),
    Source(
        "fastquorum", "samplesheet.tiny.csv",
        f"{GITHUB}/nf-core/test-datasets/{FASTQUORUM_SHA}/testdata/samplesheets/samplesheet.tiny.csv",
        "", "nf-core/test-datasets (fastquorum) sample sheet for SRR6109255",
    ),
    Source(
        "fastquorum", "SRR6109255_1.fastq.gz",
        f"{GITHUB}/nf-core/test-datasets/{FASTQUORUM_SHA}/testdata/fastqs/tiny/SRR6109255_1.fastq.gz",
        "", "nf-core/test-datasets (fastquorum) duplex reads SRR6109255, read 1",
    ),
    Source(
        "fastquorum", "SRR6109255_2.fastq.gz",
        f"{GITHUB}/nf-core/test-datasets/{FASTQUORUM_SHA}/testdata/fastqs/tiny/SRR6109255_2.fastq.gz",
        "", "nf-core/test-datasets (fastquorum) duplex reads SRR6109255, read 2",
    ),
    Source(
        "fastquorum", "samplesheet.nonrandom.csv",
        f"{GITHUB}/nf-core/test-datasets/{FASTQUORUM_SHA}/testdata/samplesheets/samplesheet.nonrandom.csv",
        "", "nf-core/test-datasets (fastquorum) sample sheet for the fixed UMIs",
    ),
    Source(
        "fastquorum", "synthetic_nonrandom_umi_1.fastq.gz",
        f"{GITHUB}/nf-core/test-datasets/{FASTQUORUM_SHA}/testdata/fastqs/synthetic-nonrandom-umi/synthetic_nonrandom_umi_1.fastq.gz",
        "", "nf-core/test-datasets (fastquorum) reads with fixed UMIs, read 1",
    ),
    Source(
        "fastquorum", "synthetic_nonrandom_umi_2.fastq.gz",
        f"{GITHUB}/nf-core/test-datasets/{FASTQUORUM_SHA}/testdata/fastqs/synthetic-nonrandom-umi/synthetic_nonrandom_umi_2.fastq.gz",
        "", "nf-core/test-datasets (fastquorum) reads with fixed UMIs, read 2",
    ),
    Source(
        "fastquorum", "idt-xgen-cfdna-umis.txt",
        f"{GITHUB}/nf-core/test-datasets/{FASTQUORUM_SHA}/testdata/umis/idt-xgen-cfdna-umis.txt",
        "", "nf-core/test-datasets (fastquorum) IDT xGen cfDNA UMI list",
    ),
    Source(
        "splitseq", "bc_8nt_v1.csv",
        f"{GITHUB}/yjzhang/split-seq-pipeline/{SPLITSEQ_SHA}/split_seq/barcodes/bc_8nt_v1.csv",
        "", "split-seq-pipeline ligation barcodes (rounds 2 and 3)",
    ),
    Source(
        "splitseq", "bc_8nt_v2.csv",
        f"{GITHUB}/yjzhang/split-seq-pipeline/{SPLITSEQ_SHA}/split_seq/barcodes/bc_8nt_v2.csv",
        "", "split-seq-pipeline reverse-transcription barcodes (round 1)",
    ),
    Source(
        "splitseq", "tools.py",
        f"{GITHUB}/yjzhang/split-seq-pipeline/{SPLITSEQ_SHA}/split_seq/tools.py",
        "", "split-seq-pipeline read layout (amp_seq, chemistry v2)",
    ),
    Source(
        "dorado", "barcode_kits.cpp",
        f"{GITHUB}/nanoporetech/dorado/{DORADO_SHA}/dorado/utils/barcode_kits.cpp",
        "", "dorado native barcodes NB01 to NB24 and their flanks",
    ),
    Source(
        "dorado", "fastq_with_us_and_tags.fq",
        f"{GITHUB}/nanoporetech/dorado/{DORADO_SHA}/tests/data/bam_reader/fastq_with_us_and_tags.fq",
        "", "dorado test read written with U bases",
    ),
]


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def write_text(path, text):
    data = text.encode()
    if path.suffix == ".gz":
        data = gzip.compress(data, mtime=0)
    path.write_bytes(data)


def write_fastq(path, records):
    write_text(path, "".join(f"@{n}\n{s}\n+\n{q}\n" for n, s, q in records))


def write_tsv(path, header, rows):
    write_text(path, "".join("\t".join(map(str, r)) + "\n" for r in [header] + list(rows)))


def read_fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            yield header[1:].rstrip("\n"), seq, qual


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def read_fasta(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">")).upper()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fmt(n):
    return f"{n:,}"


def pct(a, b):
    return f"{100 * a / b:.1f}%"


def slug(title):
    return re.sub(r"[^\w\- ]", "", title.lower()).replace(" ", "-")


def fence(text, lang):
    return f"```{lang}\n{text.rstrip(chr(10))}\n```"


def callout(kind, text):
    lines = [f"[!{kind}]"] + text.strip().split("\n")
    return "\n".join(f"> {line}".rstrip() for line in lines)


def link(text, url):
    return f"[{text}]({url})"


def details(summary, body):
    return f"<details>\n<summary>{summary}</summary>\n\n{body}\n\n</details>"


def ref(key):
    return f"\x00{key}\x00"


class Run:
    def __init__(self, command, code, out, err):
        self.command = command
        self.code = code
        self.out = out
        self.err = err


class Context:
    def __init__(self, cache, unmux, tools):
        self.cache = cache
        self.bin = cache / "bin"
        self.bin.mkdir(parents=True, exist_ok=True)
        shim = self.bin / "unmux"
        if unmux is not None:
            if shim.is_symlink() or shim.exists():
                shim.unlink()
            shim.symlink_to(unmux)
        for name, command in tools.items():
            wrapper = self.bin / name
            wrapper.write_text(f"#!/bin/sh\nexec {command} \"$@\"\n")
            wrapper.chmod(0o755)
        self.env = {
            "PATH": f"{self.bin}{os.pathsep}{os.environ.get('PATH', '')}",
            "HOME": os.environ.get("HOME", "/tmp"),
            "NO_COLOR": "1",
            "LC_ALL": "C",
            "TZ": "UTC",
        }
        self.sources = {(s.key, s.name): s for s in SOURCES}
        for source in SOURCES:
            source.sha256 = PINS.get(f"{source.key}/{source.name}", "")
        self.used = []
        self.unmux_version = self.samtools_version = None

    def versions(self):
        self.unmux_version = self.sh(self.cache, "unmux --version").out.strip()
        self.samtools_version = self.sh(self.cache, "samtools --version | head -n 1").out.strip()

    def sh(self, cwd, command, ok=(0,)):
        proc = subprocess.run(["bash", "-c", command], cwd=cwd, env=self.env, capture_output=True)
        out = proc.stdout.decode()
        err = re.sub(r"^\[\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ ", "[", proc.stderr.decode(), flags=re.M)
        if proc.returncode not in ok:
            raise SystemExit(f"command exited {proc.returncode} in {cwd}:\n{command}\n{err}{out[-2000:]}")
        return Run(command, proc.returncode, out, err)

    def fetch(self, key, name):
        source = self.sources[(key, name)]
        if source not in self.used:
            self.used.append(source)
        directory = self.cache / "downloads" / key
        path = directory / name
        if path.exists() and (not source.sha256 or sha256(path) == source.sha256):
            verify(source, path)
            return path
        directory.mkdir(parents=True, exist_ok=True)
        if source.sra:
            print(f"fetching {source.label}", file=sys.stderr)
            staging = directory / ".staging"
            if staging.exists():
                shutil.rmtree(staging)
            staging.mkdir()
            self.sh(staging, source.sra)
            for other in (s for s in SOURCES if s.sra == source.sra):
                verify(other, staging / other.name)
                (staging / other.name).replace(directory / other.name)
            shutil.rmtree(staging)
            return path
        partial = path.with_name(name + ".part")
        for attempt in range(1, 4):
            try:
                download(source, partial)
                break
            except OSError as error:
                if attempt == 3:
                    raise SystemExit(f"{source.url}: {error}")
                print(f"retrying {source.url} after: {error}", file=sys.stderr)
                time.sleep(5 * attempt)
        verify(source, partial)
        partial.replace(path)
        return path

    def workdir(self, name):
        path = self.cache / "work" / name
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True)
        return path


class CheckFailed(Exception):
    pass


def check(condition, message):
    if not condition:
        raise CheckFailed(message)


def fetch_command(ctx, key, name):
    return ctx.sources[(key, name)].sra


def verify(source, path):
    digest = sha256(path)
    if not source.sha256 and not UNPINNED_OK:
        raise SystemExit(f"{source.key}/{source.name} has no sha256 in PINS; add it (see --print-pins)")
    if source.sha256 and digest != source.sha256:
        raise SystemExit(f"{source.label} ({source.name}): sha256 {digest} does not match the pinned {source.sha256}")


def download(source, partial):
    request = urllib.request.Request(source.url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=300) as response, open(partial, "wb") as out:
        if source.records:
            with gzip.GzipFile(fileobj=response) as stream:
                for i, line in enumerate(stream):
                    if i == 4 * source.records:
                        break
                    out.write(line)
        else:
            shutil.copyfileobj(response, out)


def place(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if str(dst).endswith(".gz") and not str(src).endswith(".gz"):
        write_text(dst, Path(src).read_text())
    else:
        shutil.copyfile(src, dst)


def sam_records(ctx, cwd, path, options=""):
    out = ctx.sh(cwd, f"samtools view {options} {path}").out
    records = []
    for line in out.splitlines():
        fields = line.split("\t")
        tags = {t[:2]: t[5:] for t in fields[11:]}
        records.append((fields, tags))
    return records


def log_lines(text, keep):
    lines = [line for line in text.rstrip("\n").split("\n") if any(k in line for k in keep)]
    return "\n".join(lines)


def pair_warnings(text):
    """Keep the first tag-pair warning of each group, its count of the rest, and every other warning."""
    pair = re.compile(r"group `([^`]+)`: tags ")
    lines = [line for line in text.split("\n") if "WARN" in line]
    listed = collections.Counter(m.group(1) for m in map(pair.search, lines) if m)
    out, seen = [], set()
    for line in lines:
        m = pair.search(line)
        if not m:
            out.append(line)
        elif m.group(1) not in seen:
            seen.add(m.group(1))
            out += [line, "…"] if listed[m.group(1)] > 1 else [line]
    return "\n".join(out), max(listed.values(), default=0)


def tally(parts, total=None):
    parts = [(n, what) for n, what in parts if n]
    if len(parts) == 1 and parts[0][0] == total:
        return f"{parts[0][1]} in all {fmt(total)}"
    words = [f"{what} in {fmt(n)}" for n, what in parts]
    return words[0] if len(words) == 1 else ", ".join(words[:-1]) + " and " + words[-1]


def native_barcodes(ctx):
    text = ctx.fetch("dorado", "barcode_kits.cpp").read_text()
    front = re.search(r'NB_1st_FRONT = "([ACGT]+)"', text).group(1)
    rear = re.search(r'NB_1st_REAR = "([ACGT]+)"', text).group(1)
    codes = dict(re.findall(r'\{"(NB\d\d)", "([ACGT]+)"\}', text))
    return front, rear, [(f"NB{k:02d}", codes[f"NB{k:02d}"]) for k in range(1, 25)]


def splitseq_layout(ctx):
    tools = ctx.fetch("splitseq", "tools.py").read_text()
    block = tools.split("elif chemistry=='v2':", 1)[1]
    amp = re.search(r"amp_seq = '([NIACGT]+)'", block).group(1)
    parts = [p for p in re.split(r"(N+|I+)", amp) if p]
    linkers = [p for p in parts if p[0] not in "NI"]
    return linkers[0], linkers[1], amp


class Recipe:
    def __init__(self, section, title, blocks):
        self.section = section
        self.title = title
        self.blocks = blocks


def source_links(ctx, *keys):
    links = []
    for key, name in keys:
        source = ctx.sources[(key, name)]
        links.append(link(name, source.url))
    return ", ".join(links)


def recipe_fastq_to_ubam(ctx):
    w = ctx.workdir("01-fastq-to-ubam")
    for r in (1, 2):
        place(ctx.fetch("sra-SRR390728", f"SRR390728_{r}.fastq"), w / f"SRR390728_{r}.fastq.gz")
    command = "unmux SRR390728_1.fastq.gz SRR390728_2.fastq.gz \\\n  --rg-tag PL=ILLUMINA \\\n  --out SRR390728.unmapped.bam"
    run = ctx.sh(w, command)
    view_cmd = "samtools view --no-PG -h SRR390728.unmapped.bam | head -n 5"
    view = ctx.sh(w, view_cmd)
    compare = (
        "samtools import -1 SRR390728_1.fastq.gz -2 SRR390728_2.fastq.gz -o samtools.unmapped.bam\n"
        "diff <(samtools view samtools.unmapped.bam | cut -f 1-4,6-11) \\\n"
        "  <(samtools view SRR390728.unmapped.bam | cut -f 1-4,6-11) && echo identical"
    )
    same = ctx.sh(w, compare)
    mapq = ctx.sh(w, "samtools view samtools.unmapped.bam | cut -f 5 | sort -u").out.split()
    records = sam_records(ctx, w, "SRR390728.unmapped.bam")
    mates = [list(read_fastq(w / f"SRR390728_{r}.fastq.gz")) for r in (1, 2)]
    check(len(records) == 2 * len(mates[0]), "one record per read")
    for k, (fields, tags) in enumerate(records):
        name, seq, qual = mates[k % 2][k // 2]
        check((fields[0], fields[9], fields[10]) == (name, seq, qual), f"record {k} equals its FASTQ read")
        check(int(fields[1]) == (77 if k % 2 == 0 else 141) and tags == {"RG": "SRR390728"}, f"record {k} flags and tags")
    check(same.out.strip() == "identical" and mapq == ["0"], f"samtools import comparison: {same.out!r} {mapq}")
    return Recipe("Conversions", "Paired FASTQ to an unmapped BAM", [
        "Use this to store raw read pairs as one unmapped BAM (uBAM), the input many alignment "
        "workflows such as GATK's expect, with a read group on every record.",
        f"**Data:** the first {fmt(len(mates[0]))} read pairs of SRA run "
        f"{link('SRR390728', ctx.sources[('sra-SRR390728', 'SRR390728_1.fastq')].url)}, fetched with fg-sra:",
        fence(fetch_command(ctx, "sra-SRR390728", "SRR390728_1.fastq"), "bash"),
        fence(command, "bash"),
        "The log reports what unmux read and wrote:",
        fence(log_lines(run.err, ["input", "output", "Demultiplexed", "metrics"]), "text"),
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        "Every record carries `RG:Z:SRR390728`. The pool, which becomes the read group's ID, SM "
        "and LB, defaults to the common stem of the input names; set it with `--pool`.",
        f"Each of the {fmt(len(records))} records has the name, bases and qualities of its FASTQ read, "
        "flags 77 and 141 for the two reads of a pair, and no tag but RG.",
        "Compared with `samtools import`, the records differ only in MAPQ and the read group:",
        fence(compare, "bash"),
        fence(same.out, "text"),
        callout("NOTE", f"MAPQ is left out of the comparison: unmux writes 255 (unavailable) on unmapped "
                        f"records, where `samtools import` writes 0. All {fmt(len(records))} records match "
                        "otherwise."),
    ])


def recipe_ubam_to_fastq(ctx):
    w = ctx.workdir("02-ubam-to-fastq")
    up = ctx.cache / "work" / "01-fastq-to-ubam"
    shutil.copyfile(up / "SRR390728.unmapped.bam", w / "SRR390728.unmapped.bam")
    command = 'unmux SRR390728.unmapped.bam --out "SRR390728.R%ordinal.fastq.gz"'
    ctx.sh(w, command)
    listing = ctx.sh(w, "ls SRR390728.R*.fastq.gz")
    view = ctx.sh(w, "gzip -dc SRR390728.R1.fastq.gz | head -n 4")
    compare = (
        "samtools fastq -T RG -1 samtools.R1.fastq.gz -2 samtools.R2.fastq.gz SRR390728.unmapped.bam\n"
        "for r in 1 2; do\n"
        "  diff <(gzip -dc samtools.R$r.fastq.gz | tr '\\t' ' ') <(gzip -dc SRR390728.R$r.fastq.gz) && echo \"R$r identical\"\n"
        "done"
    )
    same = ctx.sh(w, compare)
    assert same.out.split("\n")[:2] == ["R1 identical", "R2 identical"], same.out
    original = ctx.cache / "work" / "01-fastq-to-ubam"
    back = [(s, q) for _, s, q in read_fastq(w / "SRR390728.R1.fastq.gz")]
    orig = [(s, q) for _, s, q in read_fastq(original / "SRR390728_1.fastq.gz")]
    assert back == orig
    return Recipe("Conversions", "Unmapped BAM back to paired FASTQ", [
        "Use this to hand reads stored in an unmapped BAM to a tool that wants FASTQ. The "
        "`%ordinal` placeholder writes one file per read of the pair.",
        f"**Data:** the unmapped BAM from {ref('fastq_to_ubam')}.",
        fence(command, "bash"),
        fence(listing.out, "text"),
        fence("gzip -dc SRR390728.R1.fastq.gz | head -n 4", "bash"),
        fence(view.out, "text"),
        "The tags of each record, here the read group, go into the FASTQ comment. The files match "
        "`samtools fastq -T RG` once its tab before the tags is read as a space:",
        fence(compare, "bash"),
        fence(same.out, "text"),
        f"Bases and qualities equal the original FASTQ for all {fmt(len(back))} pairs.",
    ])


def recipe_passthrough(ctx):
    w = ctx.workdir("03-passthrough")
    shutil.copyfile(ctx.fetch("modkit", "fwd_rev_modbase_records.bam"), w / "modbam.aligned.bam")
    command = "unmux modbam.aligned.bam --out modbam.unmapped.bam"
    ctx.sh(w, command)
    before_cmd = (
        "samtools view --no-PG -H modbam.aligned.bam | grep -v '^@SQ' | cut -c 1-100\n"
        "samtools view modbam.aligned.bam | cut -f 1-6\n"
        "samtools view modbam.aligned.bam | head -n 1 | cut -f 12- | tr '\\t' '\\n' | cut -d : -f 1 | paste -s -d ' ' -"
    )
    after_cmd = before_cmd.replace("modbam.aligned.bam", "modbam.unmapped.bam").replace(" | grep -v '^@SQ' | cut -c 1-100", "")
    sq = int(ctx.sh(w, "samtools view --no-PG -H modbam.aligned.bam | grep -c '^@SQ'").out)
    before = ctx.sh(w, before_cmd)
    after = ctx.sh(w, after_cmd)
    compare = (
        "samtools reset -o reset.bam modbam.aligned.bam\n"
        "diff <(samtools view reset.bam | cut -f 1-4,6-11) <(samtools view modbam.unmapped.bam | cut -f 1-4,6-11) && echo identical"
    )
    same = ctx.sh(w, compare)
    assert same.out.strip() == "identical", same.out
    aligned = sam_records(ctx, w, "modbam.aligned.bam")
    unmapped = sam_records(ctx, w, "modbam.unmapped.bam")
    assert all(a[1]["MM"] == u[1]["MM"] and a[1]["ML"] == u[1]["ML"] for a, u in zip(aligned, unmapped))
    assert all(u[0][9] == (revcomp(a[0][9]) if int(a[0][1]) & 16 else a[0][9]) for a, u in zip(aligned, unmapped))
    dropped = sorted(set(aligned[0][1]) - set(unmapped[0][1]))
    return Recipe("Conversions", "What plain SAM or BAM pass-through keeps and changes", [
        "Use this to get an unmapped copy of aligned reads, to realign them or share them without "
        "the alignment, as `samtools reset` does.",
        f"**Data:** two aligned Nanopore reads with modified-base calls (MM/ML), one on each strand "
        f"({source_links(ctx, ('modkit', 'fwd_rev_modbase_records.bam'))}).",
        fence(command, "bash"),
        f"The input's header without its {sq} `@SQ` lines, the first six columns of its records, "
        "and the tags of the first record:",
        fence(before_cmd, "bash"),
        fence(before.out, "text"),
        "The same for the output:",
        fence(after_cmd, "bash"),
        fence(after.out, "text"),
        "What changes:",
        "- Each record becomes unmapped: flag 4, no position or CIGAR, MAPQ 255. The read on the "
        "reverse strand (flag 16) gets back its sequenced orientation, bases and qualities.\n"
        f"- Tags that describe the alignment are dropped, here {', '.join(dropped)}; per-read tags "
        "such as MM, ML and the basecaller's tags are kept.\n"
        "- The header keeps `@HD VN` and gets the pool's read group and an unmux `@PG` line. The "
        "input `@RG`, `@PG`, `@CO` and `@SQ` lines are not kept, and each record's RG:Z becomes the pool.",
        "Apart from MAPQ and RG:Z, the records equal those of `samtools reset`:",
        fence(compare, "bash"),
        fence(same.out, "text"),
        "The generator also checks that MM and ML are unchanged on both records and that the "
        "reverse-strand read is the reverse complement of its aligned SEQ.",
        callout("TIP", "To keep the alignments and split them into files instead, use a read-group "
                       f"split, as in {ref('split_rg')}."),
    ])


CHANGESEQ_COMMAND = """unmux TEST.r1.fastq.gz TEST.r2.fastq.gz \\
  --extract "r1=0:0:end" --extract "r2=1:0:end" \\
  --template "~r1+r2" \\
  --out TEST_merged.fastq.gz
unmux TEST_control.r1.fastq.gz TEST_control.r2.fastq.gz \\
  --extract "r1=0:0:end" --extract "r2=1:0:end" \\
  --template "~r1+r2" \\
  --out control_TEST_merged.fastq.gz"""


def recipe_changeseq(ctx):
    w = ctx.workdir("05-changeseq")
    for s in ("TEST", "TEST_control"):
        for r in ("r1", "r2"):
            shutil.copyfile(ctx.fetch("changeseq", f"{s}.{r}.fastq.gz"), w / f"{s}.{r}.fastq.gz")
    ctx.fetch("changeseq", "mergeReads.py")
    input_cmd = "gzip -dc TEST.r1.fastq.gz | head -n 2 | cut -c 1-70\ngzip -dc TEST.r2.fastq.gz | head -n 2 | cut -c 1-70"
    inputs = ctx.sh(w, input_cmd)
    ctx.sh(w, CHANGESEQ_COMMAND)
    view_cmd = "gzip -dc TEST_merged.fastq.gz | head -n 2 | cut -c 1-70"
    view = ctx.sh(w, view_cmd)
    counts = {}
    for sample, merged in (("TEST", "TEST_merged"), ("TEST_control", "control_TEST_merged")):
        r1 = list(read_fastq(w / f"{sample}.r1.fastq.gz"))
        r2 = list(read_fastq(w / f"{sample}.r2.fastq.gz"))
        out = list(read_fastq(w / f"{merged}.fastq.gz"))
        check(len(out) == len(r1) == len(r2), f"{merged}: one read per pair")
        for (n1, s1, q1), (_, s2, q2), (name, seq, qual) in zip(r1, r2, out):
            check(name == n1.split()[0], f"{merged}: the name of {n1}")
            check(seq == revcomp(s1) + s2 and qual == q1[::-1] + q2, f"{merged}: {name} equals mergeReads()")
        counts[sample] = len(out)
    return Recipe("Conversions", "Join read pairs end to end for CHANGE-seq", [
        "Use this to turn read pairs into the single reads that CIRCLE-seq and CHANGE-seq align: the "
        "pipeline's `mergeReads()` writes the reverse complement of R1 followed by R2, with R1's "
        "qualities reversed, so both reads of a cleaved site line up on one strand.",
        f"**Data:** the treated ({fmt(counts['TEST'])} pairs) and control ({fmt(counts['TEST_control'])} "
        f"pairs) test libraries of the CHANGE-seq pipeline "
        f"({source_links(ctx, ('changeseq', 'TEST.r1.fastq.gz'), ('changeseq', 'TEST.r2.fastq.gz'), ('changeseq', 'TEST_control.r1.fastq.gz'), ('changeseq', 'TEST_control.r2.fastq.gz'))}).",
        fence(input_cmd, "bash"),
        fence(inputs.out, "text"),
        fence(CHANGESEQ_COMMAND, "bash"),
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        "The `~` reverse-complements the R1 stream, qualities reversed with it, and `+` joins the two "
        "streams into one read. Every merged read of both libraries equals what "
        f"{link('mergeReads.py', ctx.sources[('changeseq', 'mergeReads.py')].url)} writes for its "
        "pair, bases and qualities, and keeps the pair's name without its comment.",
    ])


def recipe_split_rg(ctx):
    w = ctx.workdir("04-split-by-read-group")
    shutil.copyfile(ctx.fetch("samtools", "split.sam"), w / "split.sam")
    command = 'unmux split.sam \\\n  --group "rg=@RG" \\\n  --sample-from-group rg \\\n  --out "by_rg/%sample.bam" \\\n  --unassigned by_rg/unknown.bam'
    ctx.sh(w, command)
    counts_cmd = 'for f in by_rg/*.bam; do echo "$f $(samtools view -c "$f")"; done'
    counts = ctx.sh(w, counts_cmd)
    header_cmd = "samtools view --no-PG -H by_rg/grp1.bam | grep -v '^@CO'"
    header = ctx.sh(w, header_cmd)
    compare = (
        'samtools split --no-PG -u unknown.samtools.bam -f "%!.samtools.bam" split.sam\n'
        "for rg in grp1 grp2 unknown; do\n"
        '  diff <(samtools view by_rg/$rg.bam) <(samtools view $rg.samtools.bam) && echo "$rg identical"\n'
        "done"
    )
    same = ctx.sh(w, compare)
    assert same.out.split() == ["grp1", "identical", "grp2", "identical", "unknown", "identical"], same.out
    hdiff_cmd = "diff <(samtools view --no-PG -H by_rg/grp1.bam) <(samtools view --no-PG -H grp1.samtools.bam)"
    hdiff = ctx.sh(w, hdiff_cmd, ok=(0, 1))
    total = int(ctx.sh(w, "grep -vc '^@' split.sam").out)
    return Recipe("Conversions", "Split an aligned BAM by read group", [
        "Use this to separate an aligned SAM or BAM that holds several read groups, such as lanes "
        "or samples, into one file per read group with the alignments kept whole.",
        f"**Data:** {fmt(total)} aligned records in two read groups, some with no read group or one "
        f"missing from the header ({source_links(ctx, ('samtools', 'split.sam'))}).",
        fence(command, "bash"),
        fence(counts_cmd, "bash"),
        fence(counts.out, "text"),
        fence(header_cmd, "bash"),
        fence(header.out, "text"),
        "Each file has only its own `@RG` line, and the records are those `samtools split` writes:",
        fence(compare, "bash"),
        fence(same.out, "text"),
        "The headers differ only in where `@SQ` sits and the added unmux `@PG` line:",
        fence(hdiff_cmd, "bash"),
        fence(hdiff.out, "diff"),
        callout("TIP", "With `--group \"rg=@RG::SM\"` the files are per sample, merging the read "
                       "groups of one sample, and `@RG::SM::LB` adds a library level for "
                       "`--out \"%sample.%sub_sample.bam\"`."),
    ])


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


def uditas_calls(i1, i2, rows, masked=True):
    """Each read's call under the rules of UDiTaS's demultiplex(): sample, or the bin it goes to."""
    bc1 = [r["barcode_I1"] for r in rows]
    bc2 = [r["barcode_I2"] for r in rows]
    name1 = dict(zip(bc1, (r["index_I1"] for r in rows)))
    name2 = dict(zip(bc2, (r["index_I2"] for r in rows)))
    pairs = {(r["index_I1"], r["index_I2"]) for r in rows}

    def mask(seq, qual):
        return "".join(b if ord(q) - 33 >= 12 else "N" for b, q in zip(seq, qual)) if masked else seq

    def closest(seq, barcodes):
        near = [(hamming(seq, b), b) for b in barcodes if hamming(seq, b) <= 1]
        best = min((d for d, _ in near), default=None)
        return [b for d, b in near if d == best]

    calls = {}
    for (name, s1, q1), (_, s2, q2) in zip(i1, i2):
        s1, s2 = mask(s1, q1), mask(s2[:8], q2[:8])
        if not (s1 in bc1 and s2 in bc2):
            m1, m2 = closest(s1, bc1), closest(s2, bc2)
            if not (m1 and m2):
                calls[name] = "mismatched_adapters"
                continue
            s1, s2 = m1[0], m2[0]
        pair = (name1[s1], name2[s2])
        calls[name] = f"{pair[0]}_{pair[1]}" if pair in pairs else "not_in_exp_list"
    return calls


def recipe_uditas(ctx):
    w = ctx.workdir("05-uditas")
    for r in ("R1", "R2", "I1", "I2"):
        place(ctx.fetch("uditas", f"Undetermined_S0_L001_{r}_001.fastq"), w / f"Undetermined_S0_L001_{r}_001.fastq.gz")
    shutil.copyfile(ctx.fetch("uditas", "sample_info.csv"), w / "sample_info.csv")
    ctx.fetch("uditas", "uditas_helpers.py")
    with open(w / "sample_info.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    index_cmd = "gzip -dc Undetermined_S0_L001_I1_001.fastq.gz | head -n 4\ngzip -dc Undetermined_S0_L001_I2_001.fastq.gz | head -n 4"
    index_reads = ctx.sh(w, index_cmd)
    ctx.sh(w, UDITAS_SHEET)
    sheets_cmd = "head -n 4 samples.tsv\nhead -n 3 i7.tsv\ncat i5.tsv"
    sheets = ctx.sh(w, sheets_cmd)
    ctx.sh(w, UDITAS_COMMAND)
    counts_cmd = 'for f in bam/*.bam; do echo "$f $(samtools view -c -f 64 "$f")"; done'
    counts = ctx.sh(w, counts_cmd)
    view_cmd = "samtools view bam/N702_N504.unmapped.bam | head -n 2 | cut -f 1,2,12-"
    view = ctx.sh(w, view_cmd)
    i1 = list(read_fastq(w / "Undetermined_S0_L001_I1_001.fastq.gz"))
    i2 = list(read_fastq(w / "Undetermined_S0_L001_I2_001.fastq.gz"))
    by_name = {name.split()[0]: (s, q) for name, s, q in i2}
    i1 = [(n.split()[0], s, q) for n, s, q in i1]
    i2 = [(n.split()[0], s, q) for n, s, q in i2]
    calls = {}
    for path in sorted((w / "bam").glob("*.bam")):
        sample = path.name.split(".")[0]
        for fields, tags in sam_records(ctx, w, f"bam/{path.name}"):
            seq, qual = by_name[fields[0]]
            check(tags["RX"] == seq[8:] and tags["QX"] == qual[8:], f"{fields[0]}: RX and QX equal I2 bases 8 to 18")
            if int(fields[1]) & 64:
                check(fields[0] not in calls, f"{fields[0]} is in one sample only")
                calls[fields[0]] = sample
    for name, _, _ in i1:
        calls.setdefault(name, "unassigned")
    plain = uditas_calls(i1, i2, rows, masked=False)
    masked = uditas_calls(i1, i2, rows)

    def same(a, b):
        return a == b or (a == "unassigned" and b in ("mismatched_adapters", "not_in_exp_list"))

    check(all(same(calls[n], plain[n]) for n in calls), "unmux equals the UDiTaS rules without masking")
    agree = sum(1 for n in calls if same(calls[n], masked[n]))
    differ = [n for n in calls if not same(calls[n], masked[n])]
    check(all(masked[n] == "mismatched_adapters" and calls[n] != "unassigned" for n in differ),
          "every difference is a read UDiTaS leaves unmatched and unmux assigns")
    quals = {n: (q1, q2[:8]) for (n, _, q1), (_, _, q2) in zip(i1, i2)}
    low = [sum(1 for n in differ if any(ord(c) - 33 < 12 for c in quals[n][k])) for k in (0, 1)]
    assigned = sum(1 for c in calls.values() if c != "unassigned")
    return Recipe("Illumina demultiplexing", "Dual-index demultiplexing with a sample sheet", [
        "Use this to split a pooled Illumina run into per-sample unmapped BAMs when the index reads "
        "are in their own files, here with a UMI after the i5 index in I2, as in UDiTaS "
        "(Giannoukos et al. 2018), a method for measuring gene-editing outcomes.",
        f"**Data:** the first {fmt(len(i1))} reads of the UDiTaS fig2b MiSeq run, with its sample "
        f"sheet ({source_links(ctx, ('uditas', 'Undetermined_S0_L001_R1_001.fastq'), ('uditas', 'Undetermined_S0_L001_R2_001.fastq'), ('uditas', 'Undetermined_S0_L001_I1_001.fastq'), ('uditas', 'Undetermined_S0_L001_I2_001.fastq'), ('uditas', 'sample_info.csv'))}). "
        f"Its {len(rows)} samples share the i5 index N504 and each has its own i7. I1 is the 8 bp i7 "
        "index; I2 is the 8 bp i5 index followed by a 10 bp UMI:",
        fence(index_cmd, "bash"),
        fence(index_reads.out, "text"),
        "Build a sample sheet and a tag file per index from the run's `sample_info.csv`:",
        fence(UDITAS_SHEET, "bash"),
        fence(sheets_cmd, "bash"),
        fence(sheets.out, "tsv"),
        fence(UDITAS_COMMAND, "bash"),
        "Each index group searches only its index read and corrects one substitution, and a read "
        "is assigned when its i7 and i5 match one row of the sheet. The UMI is I2 from base 8 on "
        "(counting from 0), written to RX with its qualities in QX.",
        fence(counts_cmd, "bash"),
        fence(counts.out, "text"),
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        f"{fmt(assigned)} of {fmt(len(i1))} reads were assigned. To check every call, the generator "
        "reimplements the rules of `demultiplex()` in UDiTaS's "
        f"{link('uditas_helpers.py', ctx.sources[('uditas', 'uditas_helpers.py')].url)}: each index "
        "exact, or one mismatch from its closest barcode, and the pair listed in the sheet. "
        + (f"unmux makes the same call on all {fmt(agree)} reads. UDiTaS also masks bases below "
           "quality 12 to N before matching, which unmux does not do; no read here has such a base "
           "that pushes an index past one mismatch, and on a run that has them, UDiTaS would leave "
           "those reads unmatched where unmux assigns them."
           if not differ else
           f"unmux makes the same call on {fmt(agree)} reads. UDiTaS also masks bases below quality "
           f"12 to N before matching, which unmux does not do, and the other {fmt(len(differ))} reads "
           "are all ones that masking pushes over one mismatch, so UDiTaS leaves them unmatched while "
           f"unmux assigns them ({fmt(low[0])} have a masked base in I1, {fmt(low[1])} in the i5 part "
           "of I2). With the masking left out of the reimplementation, every read agrees.")
        + " RX and QX equal I2's bases and qualities from base 8 on, on every record.",
        callout("NOTE", "UDiTaS writes the UMI with low-quality bases masked to N. unmux writes the bases as "
                        "read and their qualities to QX, so a later step can weigh them."),
    ])


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


def phix_reads(ctx, names_and_reads):
    """Names of read pairs with an exact PhiX 24-mer, either strand, in R1."""
    phix = read_fasta(ctx.fetch("ncbi", "NC_001422.1.fasta"))
    kmers = {s[i:i + 24] for s in (phix, revcomp(phix)) for i in range(len(s) - 23)}
    return {n for n, seq in names_and_reads if any(seq[i:i + 24] in kmers for i in range(len(seq) - 23))}


def bins_inputs(ctx, w):
    """The pool: UDiTaS fig2b and fig2c reads, then PhiX pairs with index reads built from their spot group."""
    origin = {}
    files = {r: [] for r in ("R1", "R2", "I1", "I2")}
    for key in ("uditas", "uditas-fig2c"):
        for r in files:
            for name, seq, qual in read_fastq(ctx.fetch(key, f"Undetermined_S0_L001_{r}_001.fastq")):
                files[r].append((name, seq, qual))
                origin[name.split()[0]] = key
    mates = [list(read_fastq(ctx.fetch("sra-SRR27862880", f"SRR27862880_{m}.fastq"))) for m in (1, 2)]
    for (n1, s1, q1), (n2, s2, q2) in zip(*mates):
        name, comment = n1.split(" ", 1)
        i7, i5 = comment.split(":")[3].split("+")
        files["R1"].append((n1, s1, q1))
        files["R2"].append((n2, s2, q2))
        files["I1"].append((name, i7, "F" * len(i7)))
        files["I2"].append((name, i5, "F" * len(i5)))
        origin[name] = "phix-run"
    check(len(origin) == len(files["R1"]), "read names are unique across the pooled runs")
    for r, records in files.items():
        write_fastq(w / f"{r}.fastq.gz", records)
    for name in ("samples.tsv", "i7.tsv", "i5.tsv"):
        shutil.copyfile(ctx.cache / "work" / "05-uditas" / name, w / name)
    shutil.copyfile(ctx.fetch("ncbi", "NC_001422.1.fasta"), w / "NC_001422.1.fasta")
    return origin, files


def recipe_bins(ctx):
    w = ctx.workdir("06-bins")
    origin, files = bins_inputs(ctx, w)
    ctx.sh(w, PHIX_KMERS)
    ctx.sh(w, BINS_COMMAND)
    listing_cmd = (
        "for f in removed/*.fastq.gz fastq/Undetermined.*.fastq.gz; do\n"
        "  echo \"$f $(gzip -dc \"$f\" | wc -l | awk '{print $1 / 4}')\"\n"
        "done"
    )
    listing = ctx.sh(w, listing_cmd)
    unassigned_cmd = "cat metrics.unassigned.tsv"
    unassigned = ctx.sh(w, unassigned_cmd)
    removed = {name.split()[0] for name, _, _ in read_fastq(w / "removed" / "PhiX.0.fastq.gz")}
    undetermined = {name.split()[0] for name, _, _ in read_fastq(w / "fastq" / "Undetermined.0.fastq.gz")}
    assigned = {}
    for path in sorted((w / "fastq").glob("N*.R1.fastq.gz")):
        for name, _, _ in read_fastq(path):
            assigned[name.split()[0]] = path.name.split(".")[0]
    phix = phix_reads(ctx, [(n.split()[0], s) for n, s, _ in files["R1"]])
    with open(ctx.cache / "work" / "05-uditas" / "sample_info.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    i1 = [(n.split()[0], s, q) for n, s, q in read_fastq(ctx.fetch("uditas", "Undetermined_S0_L001_I1_001.fastq"))]
    i2 = [(n.split()[0], s, q) for n, s, q in read_fastq(ctx.fetch("uditas", "Undetermined_S0_L001_I2_001.fastq"))]
    expected = uditas_calls(i1, i2, rows, masked=False)
    by_origin = collections.Counter(origin.values())
    removed_by = collections.Counter(origin[n] for n in removed)
    check(len(assigned) + len(removed) + len(undetermined) == len(origin), "each pair is in exactly one output")
    check(removed_by["phix-run"] == by_origin["phix-run"], "every pair of the PhiX run is in the PhiX bin")
    check(all(origin[n] == "uditas" and expected[n] == s for n, s in assigned.items()), "fig2b reads go to their sheet sample")
    check(all(origin[n] == "uditas-fig2c" for n in undetermined), "Undetermined holds only fig2c pairs")
    check(len(undetermined) == by_origin["uditas-fig2c"] - removed_by["uditas-fig2c"], "every fig2c pair not removed is undetermined")
    check(len(assigned) == by_origin["uditas"] - removed_by["uditas"], "every fig2b pair not removed is assigned")
    kept_out = sorted(n for n in removed if origin[n] != "phix-run")
    exact = sum(1 for n in kept_out if n in phix)
    loose = BINS_COMMAND.replace("phix::loc=0:0:end,dist=1", "phix::loc=0:0:end,dist=0").replace(
        "fastq/", "exact/").replace("removed/", "exact/removed.").replace("metrics.", "exact/metrics.")
    ctx.sh(w, loose)
    exact_removed = collections.Counter(origin[n.split()[0]] for n, _, _ in read_fastq(w / "exact" / "removed.PhiX.0.fastq.gz"))
    strict = (BINS_COMMAND
              .replace('"i7::loc=2:0:8,dist=1"', '"i7::loc=2:0:8,dist=1,minFindsPerGroup=1"')
              .replace('"i5::loc=3:0:8,dist=1"', '"i5::loc=3:0:8,dist=1,minFindsPerGroup=1"')
              .replace("fastq/", "strict/").replace("removed/", "strict/removed.").replace("metrics.", "strict/metrics."))
    ctx.sh(w, strict)
    strict_removed = w / "strict" / "removed.PhiX.0.fastq.gz"
    removed_strict = {n.split()[0] for n, _, _ in read_fastq(strict_removed)} if strict_removed.exists() else set()
    undetermined_strict = {n.split()[0] for n, _, _ in read_fastq(w / "strict" / "Undetermined.0.fastq.gz")}
    run_pairs = {n for n, o in origin.items() if o == "phix-run"}
    check(run_pairs <= undetermined_strict and not run_pairs & removed_strict,
          "with minFindsPerGroup on the index groups, the PhiX run's pairs go to Undetermined, not the PhiX bin")
    return Recipe("Illumina demultiplexing", "Undetermined reads and a PhiX bin", [
        "Use this to keep the reads that match no sample, the way `bcl-convert` writes "
        "Undetermined files, and to move spike-in reads such as PhiX into a bin of their own.",
        f"**Data:** a lane pooled from three real runs, so the origin of every read is known: the "
        f"{fmt(by_origin['uditas'])} fig2b reads of {ref('uditas')}, the first {fmt(by_origin['uditas-fig2c'])} "
        f"reads of UDiTaS's fig2c run, a library with an i5 that fig2b's sheet does not list "
        f"({source_links(ctx, ('uditas-fig2c', 'Undetermined_S0_L001_R1_001.fastq'))} and its other reads), "
        f"and {fmt(by_origin['phix-run'])} read pairs of the PhiX control run "
        f"{link('SRR27862880', ctx.sources[('sra-SRR27862880', 'SRR27862880_1.fastq')].url)}, fetched with fg-sra:",
        fence(fetch_command(ctx, "sra-SRR27862880", "SRR27862880_1.fastq"), "bash"),
        "SRA kept the PhiX run's index only in its spot group, `GGGGGGGG+AGATCTCG` (the i7 reads as "
        "G with no index on a two-color instrument), so the pool's I1 and I2 for those pairs are "
        "built from it, with every quality F. That is the only constructed part.",
        f"The PhiX tags are 24-mers taken every 12 bases of the PhiX174 genome, "
        f"{link('NC_001422.1', ctx.sources[('ncbi', 'NC_001422.1.fasta')].url)}:",
        fence(PHIX_KMERS, "bash"),
        fence(BINS_COMMAND, "bash"),
        "The `phix` group looks for any of the 24-mers on either strand of R1. The `--remove` rule "
        "sends the reads it matches to their own files, and `--unassigned` keeps the rest that "
        "match no sample, one file per input since `%source` is the input's number from 0 (R1, R2, "
        "I1, I2).",
        fence(listing_cmd, "bash"),
        fence(listing.out, "text"),
        fence(unassigned_cmd, "bash"),
        fence(unassigned.out, "tsv"),
        f"Checked against the origin of every pair: the PhiX bin holds all {fmt(by_origin['phix-run'])} "
        f"pairs of the PhiX run and {fmt(len(kept_out))} UDiTaS pairs, {fmt(exact)} whose R1 holds an "
        f"exact PhiX 24-mer and {fmt(len(kept_out) - exact)} with one a substitution away from a tag. "
        f"Undetermined holds the other {fmt(len(undetermined))} fig2c "
        f"pairs and nothing else, and the samples got the other {fmt(len(assigned))} fig2b pairs, each "
        f"in the sample {ref('uditas')} gives it.",
        callout("TIP", f"With `dist=0` in the `phix` group, the bin takes {fmt(exact_removed['phix-run'])} "
                       f"PhiX pairs and {fmt(sum(exact_removed.values()) - exact_removed['phix-run'])} others: "
                       "an exact match spares those UDiTaS reads but misses PhiX reads in which every "
                       "tag that fits has an error."),
        callout("WARNING", "Do not add `minFindsPerGroup` to the index groups here. A read that fails a find "
                           "constraint is unassigned before a `--remove` rule sees it, so with "
                           "`minFindsPerGroup=1` on i7 and i5 all the PhiX run's pairs go to "
                           "Undetermined and none to the PhiX bin."),
        callout("TIP", "Pass the `--out` pattern itself to `--unassigned` to keep unassigned reads in "
                       "the `--out` file instead of a file of their own."),
    ])


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


def recipe_inline(ctx):
    w = ctx.workdir("07-inline")
    origin, records = {}, []
    for run, sample, _ in CARABUS:
        for name, seq, qual in read_fastq(ctx.fetch("sra-carabus", f"{run}.fastq")):
            records.append((name, seq, qual))
            origin[name] = sample
    check(len(origin) == len(records), "read names are unique across the 32 runs")
    write_fastq(w / "pool.fastq.gz", records)
    write_tsv(w / "barcodes.tsv", ["id", "seq"], [(sample, barcode) for _, sample, barcode in CARABUS])
    head_cmd = "head -n 4 barcodes.tsv\ngzip -dc pool.fastq.gz | head -n 2"
    head = ctx.sh(w, head_cmd)
    ctx.sh(w, INLINE_COMMAND)
    reads = {name: (seq, qual) for name, seq, qual in records}
    correct = wrong = 0
    mismatched = None
    for path in sorted((w / "bam").glob("*.unmapped.bam")):
        sample = path.name.split(".")[0]
        for fields, tags in sam_records(ctx, w, f"bam/{path.name}"):
            seq, qual = reads[fields[0]]
            correct += origin[fields[0]] == sample
            wrong += origin[fields[0]] != sample
            check((tags["BC"], tags["QT"]) == (seq[:5], qual[:5]), f"{fields[0]}: BC and QT are the read's first 5 bases")
            check((fields[9], fields[10]) == (seq[5:], qual[5:]), f"{fields[0]}: the record is the read from base 5 on")
            if mismatched is None and tags["BC"] != dict((s, b) for _, s, b in CARABUS)[sample]:
                mismatched = (fields[0], sample)
    check(wrong == 0 and correct == len(records), f"every read goes to its run's sample ({correct}, {wrong})")
    view_cmd = f"samtools view bam/{mismatched[1]}.unmapped.bam | awk '$1 == \"{mismatched[0]}\"' | cut -f 1,2,12-"
    view = ctx.sh(w, view_cmd)
    exact = INLINE_COMMAND.replace("bc::loc=0:0:5,dist=1", "bc::loc=0:0:5,dist=0").replace("bam/", "exact/")
    ctx.sh(w, exact)
    exact_n = sum(len(sam_records(ctx, w, f"exact/{p.name}")) for p in (w / "exact").glob("*.bam"))
    site = sum(1 for s, _ in reads.values() if s[5:10] == "TGCAG")
    return Recipe("Illumina demultiplexing", "Inline barcodes into unmapped BAMs", [
        "Use this when the sample barcode is the first bases of the read rather than an index read, "
        "as in RAD-seq: unmux checks the restriction site after it, demultiplexes on the barcode, "
        "trims it off, and keeps it in BC with its qualities in QT.",
        f"**Data:** 250 reads of each of the 32 RAD-seq runs of ground beetles in study "
        f"{link('DRP001067', ctx.sources[('sra-carabus', 'DRR021959.fastq')].url)} (Takahashi et al. "
        "2014), one run per sample, fetched with fg-sra and pooled one after another, so the "
        "sample of every read is known. Each 101 bp read starts with the sample's 5 bp barcode and "
        "the PstI site remnant TGCAG; the barcodes are the ones each run's BioSample record "
        "(SAMD00008878 to SAMD00008909) gives, at least 3 substitutions apart.",
        fence(CARABUS_FETCH, "bash"),
        fence(head_cmd, "bash"),
        fence(head.out, "text"),
        fence(INLINE_COMMAND, "bash"),
        "The `site` group requires the remnant right after the barcode, so a read without it would "
        "be unassigned; the insert keeps the remnant, as RAD pipelines expect.",
        "A read whose barcode had one substitution:",
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        f"All {fmt(len(records))} reads went to their run's sample and none to a wrong one; "
        f"{fmt(site)} carry the remnant exactly. On every record BC and QT equal the read's first "
        "5 bases and qualities, and the record holds the rest of the read.",
        callout("NOTE", "The submitter deposited these runs already split by barcode, allowing one "
                        "substitution, so no read here has more errors than `dist=1` corrects."),
        callout("TIP", f"With `dist=0`, {fmt(exact_n)} reads are assigned and the "
                       f"{fmt(len(records) - exact_n)} with a substitution in the barcode are not."),
    ])


def read_structure(sheet, fastq_name):
    with open(sheet, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["fastq_1"].endswith("/" + fastq_name):
                return row["read_structure"]
    raise SystemExit(f"{sheet} has no row for {fastq_name}")


def swapped(umi_pair):
    a, b = umi_pair.split("-")
    return f"{b}-{a}"


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


def recipe_duplex(ctx):
    w = ctx.workdir("08-duplex-random")
    for name in ("SRR6109255_1.fastq.gz", "SRR6109255_2.fastq.gz"):
        shutil.copyfile(ctx.fetch("fastquorum", name), w / name)
    structure = read_structure(ctx.fetch("fastquorum", "samplesheet.tiny.csv"), "SRR6109255_1.fastq.gz")
    assert structure == "10M1S+T 10M1S+T", structure
    input_cmd = "gzip -dc SRR6109255_1.fastq.gz | head -n 4 | cut -c 1-60\ngzip -dc SRR6109255_2.fastq.gz | head -n 4 | cut -c 1-60"
    inputs = ctx.sh(w, input_cmd)
    ctx.sh(w, DUPLEX_COMMAND)
    view_cmd = "samtools view SRR6109255.unmapped.bam | head -n 2 | cut -f 1,2,12-\nsamtools view SRR6109255.unmapped.bam | head -n 2 | cut -f 10 | cut -c 1-49"
    view = ctx.sh(w, view_cmd)
    r1 = list(read_fastq(w / "SRR6109255_1.fastq.gz"))
    r2 = list(read_fastq(w / "SRR6109255_2.fastq.gz"))
    records = sam_records(ctx, w, "SRR6109255.unmapped.bam")
    assert len(records) == 2 * len(r1) == 2 * len(r2)
    pairs = collections.Counter()
    for k, (a, b) in enumerate(zip(r1, r2)):
        (f1, t1), (f2, t2) = records[2 * k], records[2 * k + 1]
        rx, qx = f"{a[1][:10]}-{b[1][:10]}", f"{a[2][:10]} {b[2][:10]}"
        assert f1[0] == f2[0] == a[0].split()[0] and int(f1[1]) & 64 and int(f2[1]) & 128
        assert t1["RX"] == t2["RX"] == rx and t1["QX"] == t2["QX"] == qx, (rx, qx, t1, t2)
        assert (f1[9], f1[10], f2[9], f2[10]) == (a[1][11:], a[2][11:], b[1][11:], b[2][11:])
        pairs[rx] += 1
    both = {p for p in pairs if swapped(p) in pairs and swapped(p) != p}
    count = ctx.sh(w, SWAP_COMMAND)
    assert int(count.out) == len(both), (count.out, len(both))
    spacers = collections.Counter(s[10] for _, s, _ in r1 + r2)
    return Recipe("UMIs and duplex sequencing", "Random duplex UMIs into RX", [
        "Use this for duplex sequencing libraries with a random UMI at the start of each read: one "
        "RX tag gets both UMIs, R1's first, the UMIs and spacers come off the reads, and the "
        "unmapped BAM is ready for fgbio's duplex workflow.",
        f"**Data:** the {fmt(len(r1))} read pairs of SRR6109255 in the nf-core/fastquorum test set "
        f"({source_links(ctx, ('fastquorum', 'SRR6109255_1.fastq.gz'), ('fastquorum', 'SRR6109255_2.fastq.gz'), ('fastquorum', 'samplesheet.tiny.csv'))}). "
        f"The sample sheet gives the read structure `{structure}`: in each read a 10 bp UMI and a "
        f"1 bp spacer, then the template. The spacer is a T in {pct(spacers['T'], len(r1) + len(r2))} "
        "of reads.",
        fence(input_cmd, "bash"),
        fence(inputs.out, "text"),
        fence(DUPLEX_COMMAND, "bash"),
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        f"On all {fmt(len(r1))} pairs, RX is R1's first 10 bases, a hyphen, then R2's, QX holds their "
        "qualities joined by a space, and both reads start after their spacer, at base 11 counting "
        "from 0. The pool, and so the read group, is the common stem of the input names, "
        "`SRR6109255`.",
        f"In this sample {fmt(len(both))} of the {fmt(len(pairs))} UMI pairs also appear swapped, "
        f"as B-A for an A-B: {fmt(len(both) // 2)} molecules read from both strands, in "
        f"{fmt(sum(pairs[p] for p in both))} read pairs. These are the pairs the paired strategy joins:",
        fence(SWAP_COMMAND, "bash"),
        fence(count.out, "text"),
        callout("IMPORTANT", "Keep R1's UMI first. The read pairs from one strand of a molecule carry A-B "
                             "and those from the other strand B-A, since their R1 starts from the "
                             "other end. fgbio's paired strategy joins the two by that swap and tells "
                             "the strands apart by it; sorting the two UMIs, or putting R2's first, "
                             "would lose the strand."),
        callout("TIP", "Next, align the reads, bring the tags back onto the alignments with fgbio "
                       "ZipperBams, and run fgbio GroupReadsByUmi with `--strategy paired`, then "
                       "CallDuplexConsensusReads."),
    ])


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


def recipe_duplex_fixed(ctx):
    w = ctx.workdir("09-duplex-fixed")
    for name in ("synthetic_nonrandom_umi_1.fastq.gz", "synthetic_nonrandom_umi_2.fastq.gz", "idt-xgen-cfdna-umis.txt"):
        shutil.copyfile(ctx.fetch("fastquorum", name), w / name)
    structure = read_structure(ctx.fetch("fastquorum", "samplesheet.nonrandom.csv"), "synthetic_nonrandom_umi_1.fastq.gz")
    assert structure == "8M+T 8M+T", structure
    prep_cmd = UMI_PREP + "\nhead -n 3 umis.tsv"
    prep = ctx.sh(w, prep_cmd)
    umis = [row["seq"] for row in read_tsv(w / "umis.tsv")]
    input_cmd = "gzip -dc synthetic_nonrandom_umi_1.fastq.gz | head -n 4 | cut -c 1-60"
    inputs = ctx.sh(w, input_cmd)
    ctx.sh(w, FIXED_COMMAND)
    r1 = list(read_fastq(w / "synthetic_nonrandom_umi_1.fastq.gz"))
    r2 = list(read_fastq(w / "synthetic_nonrandom_umi_2.fastq.gz"))
    for made, end in ((r1, 1), (r2, 2)):
        real = list(read_fastq(ctx.fetch("fastquorum", f"SRR6109255_{end}.fastq.gz")))
        check([m[0] for m in made] == [r[0] for r in real], f"the R{end} reads are those of SRR6109255, in order")
        check(all(m[1][8:] == r[1][11:] and m[2][8:] == r[2][11:] and m[2][:8] == "I" * 8 for m, r in zip(made, real)),
              f"each R{end} read is 8 bases of quality I and the SRR6109255 read after its first 11 bases")
    records = sam_records(ctx, w, "umis.unmapped.bam")
    assert len(records) == 2 * len(r1)

    def nearest(seen):
        best = sorted((hamming(seen, u), u) for u in umis)
        assert best[0][0] <= 1 and best[1][0] > best[0][0], (seen, best[:2])
        return best[0][1]

    fixed = [0, 0]
    pick = None
    for k, (a, b) in enumerate(zip(r1, r2)):
        (f1, t1), (f2, t2) = records[2 * k], records[2 * k + 1]
        n1, n2 = nearest(a[1][:8]), nearest(b[1][:8])
        fixed[0] += n1 != a[1][:8]
        fixed[1] += n2 != b[1][:8]
        for t in (t1, t2):
            assert t["RX"] == f"{n1}-{n2}" and t["OX"] == f"{a[1][:8]}-{b[1][:8]}" and "QX" not in t
            assert t["BZ"] == f"{a[2][:8]} {b[2][:8]}"
        assert (f1[9], f2[9]) == (a[1][8:], b[1][8:])
        if pick is None and n1 != a[1][:8]:
            pick = f1[0]
    view_cmd = f"samtools view umis.unmapped.bam | awk '$1 == \"{pick}\"' | cut -f 1,2,12-"
    view = ctx.sh(w, view_cmd)
    exact = FIXED_COMMAND.replace("dist=1", "dist=0").replace("umis.unmapped.bam", "umis.exact.bam")
    ctx.sh(w, exact)
    exact_cmd = "samtools view -c -f 64 umis.exact.bam\nsamtools view -c -f 64 -d RX umis.exact.bam"
    counted = ctx.sh(w, exact_cmd)
    either = sum(1 for a, b in zip(r1, r2) if a[1][:8] not in umis or b[1][:8] not in umis)
    exact_records = sam_records(ctx, w, "umis.exact.bam")
    assert all("OX" in t and "BZ" in t for _, t in exact_records)
    assert counted.out.split() == [str(len(r1)), str(len(r1) - either)], (counted.out, either)
    wide = ctx.sh(w, FIXED_COMMAND.replace("dist=1", "dist=2"), ok=(1,))
    check("can both match" in wide.err, "dist=2 is refused as a collision")
    aside = exact.replace("dist=0", "dist=0,minFindsPerGroup=1").replace(
        "--out umis.exact.bam", "--out umis.kept.bam \\\n  --unassigned umis.aside.bam")
    ctx.sh(w, aside)
    kept = sam_records(ctx, w, "umis.kept.bam", "-f 64")
    set_aside = sam_records(ctx, w, "umis.aside.bam", "-f 64")
    check(len(set_aside) == either and len(kept) == len(r1) - either and all("RX" in t for _, t in kept),
          "with minFindsPerGroup=1 the pairs without a corrected UMI go to --unassigned")
    return Recipe("UMIs and duplex sequencing", "Fixed duplex UMIs corrected to a list", [
        "Use this for duplex kits with a fixed set of UMIs, such as IDT xGen cfDNA: each read's UMI "
        "is corrected to the list, RX gets the corrected pair, and OX keeps the bases as read with "
        "their qualities in BZ.",
        "**Data:** the nf-core/fastquorum test set for fixed UMIs and the UMI list "
        f"({source_links(ctx, ('fastquorum', 'synthetic_nonrandom_umi_1.fastq.gz'), ('fastquorum', 'synthetic_nonrandom_umi_2.fastq.gz'), ('fastquorum', 'idt-xgen-cfdna-umis.txt'), ('fastquorum', 'samplesheet.nonrandom.csv'))}). "
        f"Its {fmt(len(r1))} read pairs are those of {ref('duplex')} with the first 11 bases of each "
        "read, the random UMI and the base after it, replaced by 8 bases of quality I: one of the "
        f"{len(umis)} IDT xGen cfDNA UMIs or, in {fmt(fixed[0] + fixed[1])} of the "
        f"{fmt(2 * len(r1))} reads, one substitution away from one. The rest of each read is "
        f"unchanged. The read structure is `{structure}`: an 8 bp UMI, then the template.",
        fence(input_cmd, "bash"),
        fence(inputs.out, "text"),
        "The list holds one UMI per line, so give it the `id` and `seq` columns unmux reads:",
        fence(prep_cmd, "bash"),
        fence(prep.out, "tsv"),
        fence(FIXED_COMMAND, "bash"),
        "A pair whose R1 UMI had one substitution:",
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        f"On all {fmt(len(r1))} pairs RX holds the two list UMIs and OX the bases as read, with BZ "
        f"their qualities. {fmt(fixed[0])} R1 UMIs and {fmt(fixed[1])} R2 UMIs had one substitution "
        "and were corrected. Every UMI in this data is within one substitution of the list, so no "
        "pair is left uncorrected here.",
        "To see what happens to a UMI that cannot be corrected, allow none (`dist=0`):",
        details("The command", fence(exact, "bash")),
        fence(exact_cmd, "bash"),
        fence(counted.out, "text"),
        f"The {fmt(either)} pairs with a substitution in either UMI are still written, with OX and BZ "
        "but no RX: a tag is left out when any stream in it is missing.",
        callout("TIP", "To set such pairs aside instead, add `minFindsPerGroup=1` to both groups and "
                       "give `--unassigned`."),
        callout("NOTE", "The 32 UMIs are at least 3 substitutions apart, so `dist=1` corrects without "
                        "ties. With `dist=2` unmux refuses the set, since one read could then match "
                        "two UMIs equally well."),
    ])


def tenx_inputs(ctx, w):
    for name in ("S10_L001_R1_001.fastq.gz", "S10_L001_R2_001.fastq.gz", "10xv2_whitelist.txt"):
        dst = w / name
        if not dst.exists():
            shutil.copyfile(ctx.fetch("nfcore", name), dst)


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


def recipe_tenx(ctx):
    w = ctx.workdir("08-10x")
    tenx_inputs(ctx, w)
    ctx.sh(w, CELLS_COMMAND)
    cells = [row["seq"] for row in read_tsv(w / "cells.tsv")]
    keep = {line.strip() for line in open(w / "10xv2_whitelist.txt") if line.strip()}
    run = ctx.sh(w, TENX_COMMAND)
    records = sam_records(ctx, w, "S10.unmapped.bam")
    r1 = {name.split()[0]: (seq, qual) for name, seq, qual in read_fastq(w / "S10_L001_R1_001.fastq.gz")}
    r2 = {name.split()[0]: (seq, qual) for name, seq, qual in read_fastq(w / "S10_L001_R2_001.fastq.gz")}
    check(len(records) == len(r1), "one record per read pair")

    def neighbors(observed):
        return {observed[:p] + b + observed[p + 1:] for p in range(16) for b in BASES if b != observed[p]}

    seen = {seq[:16] for seq, _ in r1.values()}
    reach = {s for o in seen for s in neighbors(o) | {o} if s in keep}
    check(len(cells) == len(set(cells)) and set(cells) == reach,
          "cells.tsv holds each keeplist barcode within one substitution of a barcode in R1, once")

    def oracle(observed):
        if observed in keep:
            return observed, "exact"
        hits = neighbors(observed) & keep
        return (hits.pop(), "corrected") if len(hits) == 1 else (None, "tied" if hits else "far")

    agree = exact = corrected = 0
    kinds = collections.Counter()
    pick_fixed = pick_none = None
    for fields, tags in records:
        seq, qual = r1[fields[0]]
        want, kind = oracle(seq[:16])
        kinds[kind] += 1
        got = tags.get("CB")
        agree += got == want
        check(tags["CR"] == seq[:16] and tags["CY"] == qual[:16], f"{fields[0]}: CR and CY are R1 bases 0 to 16")
        check(tags["UR"] == tags["UB"] == seq[16:26] and tags["UY"] == qual[16:26], f"{fields[0]}: UR, UB and UY are R1 bases 16 to 26")
        check((fields[9], fields[10]) == r2[fields[0]], f"{fields[0]}: the record is R2")
        if got is None:
            pick_none = pick_none or (fields[0] if kind == "far" else None)
        elif got == seq[:16]:
            exact += 1
        else:
            corrected += 1
            pick_fixed = pick_fixed or fields[0]
    check(agree == len(records), f"CB agrees with the whole keeplist on every read ({agree} of {len(records)})")
    none = kinds["far"] + kinds["tied"]
    view_cmd = f"samtools view S10.unmapped.bam | grep -E '^({pick_fixed}|{pick_none})[[:space:]]' | cut -f 1,2,12-"
    view = ctx.sh(w, view_cmd)
    shown = [line for line in run.err.split("\n") if "WARN" in line]
    excerpt = "\n".join(shown[:3] + (["…"] + shown[-1:] if len(shown) > 4 else shown[3:]))
    return Recipe("Single cell", "10x Chromium cell barcodes and UMIs", [
        "Use this to tag 10x Chromium reads with their cell barcode and UMI in an unmapped BAM, "
        "correcting each barcode to a list of known cells: CB holds the corrected barcode, CR and "
        "CY the barcode as read with its qualities, and UR, UY and UB the UMI.",
        f"**Data:** the {fmt(len(records))} 10x v2 read pairs of the nf-core/scrnaseq test set "
        f"(R1: 16 bp barcode and 10 bp UMI; R2: cDNA) and the 10x v2 barcode keeplist "
        f"({source_links(ctx, ('nfcore', 'S10_L001_R1_001.fastq.gz'), ('nfcore', 'S10_L001_R2_001.fastq.gz'), ('nfcore', '10xv2_whitelist.txt'))}).",
        f"The keeplist holds {fmt(len(keep))} barcodes, and this recipe uses a subset of it: the "
        f"{fmt(len(cells))} within one substitution of one of the {fmt(len(seen))} barcodes seen in R1. "
        "A read can be corrected only to such a barcode, so the subset gives every read the call the "
        "whole list would. unmux reads tag files with `id` and `seq` columns, so the barcode is both:",
        fence(CELLS_COMMAND, "bash"),
        fence(TENX_COMMAND, "bash"),
        "A read whose barcode was corrected and one with no keeplist barcode within one substitution, "
        "in file order:",
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        f"Of {fmt(len(records))} reads, {fmt(exact)} matched a keeplist barcode exactly, "
        f"{fmt(corrected)} were corrected, and {fmt(none)} have no CB but keep CR, UR and UB: "
        f"{fmt(kinds['far'])} with no keeplist barcode within one substitution"
        + (f" and {fmt(kinds['tied'])} with two or more." if kinds["tied"] else ".")
        + " On every read, CB agrees with a lookup of the barcode and its one-substitution neighbors "
        "in the whole keeplist, CR, CY, UR, UB and UY equal the R1 bases and qualities they come "
        "from, and the record holds R2.",
        "unmux warns about tags one read can match equally well, which `mode=nearest` leaves "
        "uncorrected, and that no reads are split by sample, which is expected here:",
        fence(excerpt, "text"),
        callout("NOTE", "unmux does not correct UMIs, so UB holds the same bases as UR; tools such as "
                        "velocyto read UB."),
        callout("WARNING", "Without `--pool S10` the pool, and with it the read group, would be "
                           "`S10_L001_R`: the common stem of the two file names."),
    ])


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


GEO_CELLS = "GSM5169184_C2C12_short_1k_cell_metadata.csv.gz"

SPLITSEQ_TAGS = """awk -F , 'BEGIN { print "id\\tseq" } { print "bc" $1 "\\t" $2 }' bc_8nt_v2.csv > round1.tags.txt
awk -F , 'BEGIN { print "id\\tseq" } { print "bc" $1 "\\t" $2 }' bc_8nt_v1.csv > round2_3.tags.txt"""

SPLITSEQ_CELLS = """gzip -dc GSM5169184_C2C12_short_1k_cell_metadata.csv.gz | awk -F , 'NR == FNR { rt[$1] = $2; next }
  FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; print "id\\tseq\\tcell"; next }
  { cb = $c["cell_barcode"]; w = $c["rnd1_well"]; r3 = substr(cb, 1, 8); r2 = substr(cb, 9, 8)
    print cb "_dT\\t" rt[w] r2 r3 "\\t" cb; print cb "_hex\\t" rt[w + 48] r2 r3 "\\t" cb }' bc_8nt_v2.csv - > cells.tsv"""


def splitseq_inputs(ctx, w):
    for m in (1, 2):
        place(ctx.fetch("sra-SRR13948565", f"SRR13948565_{m}.fastq"), w / f"R{m}.fastq.gz")
    for name in ("bc_8nt_v1.csv", "bc_8nt_v2.csv"):
        shutil.copyfile(ctx.fetch("splitseq", name), w / name)
    shutil.copyfile(ctx.fetch("geo", GEO_CELLS), w / GEO_CELLS)
    ctx.sh(w, SPLITSEQ_TAGS)
    ctx.sh(w, SPLITSEQ_CELLS)


def splitseq_oracle(w, windows=((78, 86), (48, 56), (10, 18))):
    """Each read's CB, CY, RX and QX under the README windows: a unique tag within one substitution per round."""
    round1 = [r["seq"] for r in read_tsv(w / "round1.tags.txt")]
    round23 = [r["seq"] for r in read_tsv(w / "round2_3.tags.txt")]

    def fix(observed, tags):
        hits = [t for t in tags if len(observed) == len(t) and hamming(observed, t) <= 1]
        return hits[0] if len(hits) == 1 else None

    calls = {}
    for name, seq, qual in read_fastq(w / "R2.fastq.gz"):
        found = [fix(seq[a:b], tags) for (a, b), tags in zip(windows, (round1, round23, round23))]
        calls[name] = None if None in found else (
            "-".join(found), " ".join(qual[a:b] for a, b in windows), seq[:10], qual[:10])
    return calls


def called_cells(ctx, w):
    """The cells split-pipe called, by joined round 1, 2 and 3 barcode, checked against the published wells."""
    rows = read_tsv(w / "cells.tsv")
    v1 = {int(line.split(",")[0]): line.strip().split(",")[1] for line in open(w / "bc_8nt_v1.csv") if line.strip()}
    with gzip.open(w / GEO_CELLS, "rt") as fh:
        for row in csv.DictReader(fh):
            cb = row["cell_barcode"]
            check(cb[:8] == v1[int(row["rnd3_well"]) + 1] and cb[8:16] == v1[int(row["rnd2_well"]) + 1],
                  f"{cb}: rounds 3 and 2 are the barcodes of its wells")
    return {r["seq"]: r["cell"] for r in rows}


def splitseq_scores(ctx, w, path, cells, reference=None):
    """Written reads, with a CB, with a called cell's CB, and with a CB other than the reference run's."""
    written = with_cb = in_cells = changed = 0
    for fields, tags in sam_records(ctx, w, path):
        written += 1
        cb = tags.get("CB")
        if cb is None:
            continue
        with_cb += 1
        in_cells += cb.replace("-", "") in cells
        if reference is not None and reference.get(fields[0]) not in (None, cb):
            changed += 1
    return written, with_cb, in_cells, changed


def recipe_splitseq(ctx):
    w = ctx.workdir("09-splitseq")
    splitseq_inputs(ctx, w)
    _, _, amp = splitseq_layout(ctx)
    tags_cmd = SPLITSEQ_TAGS + "\nhead -n 3 round1.tags.txt"
    tags_out = ctx.sh(w, tags_cmd)
    run = ctx.sh(w, SPLITSEQ_COMMAND)
    view_cmd = "samtools view splitseq.unmapped.bam | head -n 2 | cut -f 1,2,12-"
    view = ctx.sh(w, view_cmd)
    oracle = splitseq_oracle(w)
    r1 = {n: (s, q) for n, s, q in read_fastq(w / "R1.fastq.gz")}
    records = sam_records(ctx, w, "splitseq.unmapped.bam")
    for fields, tags in records:
        want = oracle[fields[0]]
        check(want is not None, f"{fields[0]}: written only when every round matches")
        check((tags["CB"], tags["CY"], tags["RX"], tags["QX"]) == want, f"{fields[0]}: CB, CY, RX and QX match the oracle")
        check((fields[9], fields[10]) == r1[fields[0]], f"{fields[0]}: the record is R1")
    check(len(records) == sum(1 for v in oracle.values() if v), "every read the oracle calls is written")
    cells = called_cells(ctx, w)
    kept_cells = collections.Counter(cells[t["CB"].replace("-", "")] for _, t in records if t["CB"].replace("-", "") in cells)
    by_primer = collections.defaultdict(set)
    for _, t in records:
        joined = t["CB"].replace("-", "")
        if joined in cells:
            by_primer[cells[joined]].add(joined)
    both = sum(1 for seqs in by_primer.values() if len(seqs) == 2)
    r2 = next(read_fastq(w / "R2.fastq.gz"))
    reads = len(oracle)
    return Recipe("Single cell", "SPLiT-seq barcodes in three rounds", [
        "Use this for split-pool data such as SPLiT-seq, where the cell barcode is three 8 bp "
        "barcodes, one per round of barcoding, separated by fixed linkers in R2. This is the "
        "README's Quick Start command.",
        f"**Data:** {fmt(reads)} read pairs of SRA run "
        f"{link('SRR13948565', ctx.sources[('sra-SRR13948565', 'SRR13948565_1.fastq')].url)}, the "
        "C2C12_short_1k library of LR-Split-seq (Rebboah et al. 2021) made with the Parse Biosciences "
        "v1 kit, whose R2 has the chemistry v2 layout of the split-seq-pipeline: a 10 bp UMI, round "
        f"3, linker 1, round 2, linker 2 and round 1 (`{amp}`, where N is the UMI and I a barcode). "
        f"The barcodes are the pipeline's published lists "
        f"({source_links(ctx, ('splitseq', 'bc_8nt_v2.csv'), ('splitseq', 'bc_8nt_v1.csv'), ('splitseq', 'tools.py'))}).",
        fence(fetch_command(ctx, "sra-SRR13948565", "SRR13948565_1.fastq"), "bash"),
        fence(f"@{r2[0]}\n{r2[1]}\n+\n{r2[2]}", "text"),
        "Make the two tag files from the published lists:",
        fence(tags_cmd, "bash"),
        fence(tags_out.out, "tsv"),
        fence(SPLITSEQ_COMMAND, "bash"),
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        "CB joins the three corrected barcodes with `-`, and CY their qualities with spaces; RX and "
        "QX hold the UMI. A read that misses any round fails `minFindsPerGroup=1` and is dropped.",
        f"{fmt(len(records))} of {fmt(reads)} reads were written. The generator looks up each read's "
        "three windows itself, a published barcode within one substitution in each, and every written "
        "read has the CB, CY, RX and QX that lookup gives, while every read it drops misses a round; "
        f"each record holds its R1. {fmt(sum(kept_cells.values()))} of the written reads belong to "
        f"{fmt(len(kept_cells))} of the {fmt(len(set(cells.values())))} cells that split-pipe called "
        f"in this library ({link('GSM5169184', ctx.sources[('geo', GEO_CELLS)].url)}).",
        fence(log_lines(run.err, ["WARN"]), "text"),
        "The warning is expected: this command tags reads and does not split them by sample.",
        callout("NOTE", "The Parse kit primes each round-1 well with an oligo-dT and a random-hexamer "
                        "barcode, so one cell's reads carry either, and this command gives the cell two "
                        f"CBs; {fmt(both)} called cells appear under both here. {ref('match')} joins "
                        "them through a cell list."),
        callout("NOTE", "Without `--pool`, the pool is the common stem of the input names, here `R` "
                        "from R1 and R2, so every record gets `RG:Z:R`; add `--pool` to set it."),
    ])


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


def recipe_match(ctx):
    w = ctx.workdir("10-match-cells")
    splitseq_inputs(ctx, w)
    cells_cmd = SPLITSEQ_CELLS + "\nhead -n 3 cells.tsv"
    cells_out = ctx.sh(w, cells_cmd)
    run = ctx.sh(w, MATCH_COMMAND)
    view_cmd = "samtools view cells.unmapped.bam | head -n 2 | cut -f 1,2,12-"
    view = ctx.sh(w, view_cmd)
    cells = called_cells(ctx, w)
    oracle = splitseq_oracle(w)
    kept = sam_records(ctx, w, "cells.unmapped.bam")
    other = {f[0] for f, _ in sam_records(ctx, w, "other.unmapped.bam")}
    want = {n: cells[c[0].replace("-", "")] for n, c in oracle.items() if c and c[0].replace("-", "") in cells}
    check({f[0]: t.get("CB") for f, t in kept} == want, "kept reads are the oracle's called-cell reads, each with its cell")
    check(len(kept) + len(other) == len(oracle) and not other & set(want), "every other read goes to other.unmapped.bam")
    r2 = {n: s for n, s, _ in read_fastq(w / "R2.fastq.gz")}
    observed = {n: "".join(s[a:b] for a, b in ((78, 86), (48, 56), (10, 18))) for n, s in r2.items()}
    corrected = sum(1 for f, _ in kept if oracle[f[0]][0].replace("-", "") != observed[f[0]])
    names = len(set(want.values()))
    warn = log_lines(run.err, ["WARN"])
    share = pct(len(other), len(oracle))
    check("no `--sample`" in warn and "`CY` is left out" in warn and f"{share} of reads were unassigned" in warn,
          "unmux warns about the missing sample target, the qualities of CB and the unassigned share")
    return Recipe("Single cell", "Keep the reads of called cells with a joined match", [
        "Use this to keep only reads from cells you have called and write each cell's name: a "
        "`match=` group looks up the three corrected round barcodes, joined, in a list of cells, "
        "and a column of the list gives the name.",
        f"**Data:** the SPLiT-seq reads of {ref('splitseq')} and the {fmt(len(set(cells.values())))} "
        f"cells split-pipe called in that library "
        f"({link('GSM5169184', ctx.sources[('geo', GEO_CELLS)].url)}). Each cell's `cell_barcode` is "
        "its round-3 and round-2 barcodes, then its round-1 well, and the well's oligo-dT and "
        "random-hexamer barcodes are rows w and w + 48 of the round-1 list, so the list gets a row "
        "for each, with the same `cell`:",
        fence(cells_cmd, "bash"),
        fence(cells_out.out, "tsv"),
        fence(MATCH_COMMAND, "bash"),
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        f"{fmt(len(kept))} reads were kept, of {fmt(names)} cells, and {fmt(len(other))} went to "
        "other.unmapped.bam, both records of each pair as they were read. Checked against the "
        "generator's own lookup of the three windows: the kept reads are exactly those whose joined "
        "barcode is a called cell's, each with that cell's name, both primers of a cell under one "
        f"name. {fmt(corrected)} of them needed a correction in some round and are kept because "
        "`match=` joins each round as corrected, the way `--tag` writes it.",
        f"unmux warns, as in {ref('splitseq')}, that nothing is split by sample; that CB, a column of "
        "cells.tsv, has no base qualities to write to CY, which `CB::qual=none` silences; and that "
        f"{share} of the reads were unassigned, the {fmt(len(other))} in other.unmapped.bam:",
        fence(warn, "text"),
        callout("TIP", "Without `minFindsPerGroup=1` on the `match=` group, reads of other barcodes "
                       "are written too, without a CB."),
    ])


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


def recipe_cite(ctx):
    w = ctx.workdir("11-cite-seq")
    names = ("subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz",
             "subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz",
             "5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv")
    for name in names:
        shutil.copyfile(ctx.fetch("nfcore-modules", name), w / name)
    prep_cmd = CITE_PREP + "\nhead -n 3 5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv\nhead -n 3 antibodies.tsv"
    features = ctx.sh(w, prep_cmd)
    input_cmd = f"gzip -dc {names[1]} | head -n 2"
    inputs = ctx.sh(w, input_cmd)
    ctx.sh(w, CITE_COMMAND)
    view_cmd = "samtools view adt.unmapped.bam | head -n 2 | cut -f 1,2,12-"
    view = ctx.sh(w, view_cmd)
    count_cmd = (
        "samtools view adt.unmapped.bam | grep -o 'fx:Z:[^[:space:]]*' \\\n"
        "  | sort | uniq -c | awk '{ print $1, $2 }' | sort -k1,1nr -k2,2 | head -n 5"
    )
    counts = ctx.sh(w, count_cmd)
    panel = {row["id"]: row["seq"] for row in read_tsv(w / "antibodies.tsv")}
    gap = min(hamming(a, b) for a in panel.values() for b in panel.values() if a != b)
    r1 = {n.split()[0]: (s, q) for n, s, q in read_fastq(w / names[0])}
    r2 = {n.split()[0]: (s, q) for n, s, q in read_fastq(w / names[1])}

    def lookup(window):
        hits = [i for i, s in panel.items() if hamming(window, s) <= 1 and len(window) == len(s)]
        return hits[0] if len(hits) == 1 else None

    want = {n: lookup(s[10:25]) for n, (s, _) in r2.items()}
    records = sam_records(ctx, w, "adt.unmapped.bam")
    for fields, tags in records:
        (s1, q1), (s2, q2) = r1[fields[0]], r2[fields[0]]
        check(tags["fx"] == want[fields[0]] and tags["fb"] == panel[tags["fx"]], f"{fields[0]}: fx and fb match the panel lookup")
        check((tags["CR"], tags["CY"], tags["UR"], tags["UY"]) == (s1[:16], q1[:16], s1[16:], q1[16:]), f"{fields[0]}: CR, CY, UR, UY are R1")
        check((fields[9], fields[10]) == (s2, q2), f"{fields[0]}: the record is R2")
    check(len(records) == sum(1 for v in want.values() if v), "every read whose window holds a panel barcode is written")
    exact = sum(1 for n, (s, _) in r2.items() if s[10:25] in panel.values())
    return Recipe("Single cell", "CITE-seq antibody barcodes", [
        "Use this for antibody-derived tags (ADT) in CITE-seq or 10x feature barcoding: find the "
        "antibody barcode in R2, correct it, and write the antibody's name next to the cell "
        "barcode and UMI.",
        f"**Data:** the {fmt(len(r2))} antibody read pairs of the 10x 5k anti-CMV T cell set in "
        f"nf-core/test-datasets, with their feature reference of {len(panel)} TotalSeq-C antibodies "
        f"({source_links(ctx, *(('nfcore-modules', n) for n in names))}). R1 is a 16 bp cell barcode "
        "and the UMI; in R2 the 15 bp antibody barcode follows 10 other bases, as the reference's "
        "pattern `^NNNNNNNNNN(BC)NNNNNNNNN` gives:",
        fence(input_cmd, "bash"),
        fence(inputs.out, "text"),
        "Give the feature reference the `id` and `seq` columns unmux reads:",
        fence(prep_cmd, "bash"),
        fence(features.out, "tsv"),
        fence(CITE_COMMAND, "bash"),
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        "The `@ab::id` stream is the matched antibody's `id` column, written as text to fx; fb "
        "holds the corrected barcode.",
        fence(count_cmd, "bash"),
        fence(counts.out, "text"),
        f"{fmt(len(records))} of {fmt(len(r2))} reads were written: every read whose R2 bases 10 to 25 "
        f"hold a panel barcode within one substitution ({fmt(exact)} of them exactly), checked against "
        "the generator's own lookup. On every record fx and fb match that lookup, CR, CY, UR and UY "
        f"equal the R1 bases and qualities, and the record holds R2. The panel's barcodes are at least "
        f"{gap} substitutions apart, so one substitution cannot reach a wrong antibody.",
        callout("NOTE", "The tags fb and fx follow Cell Ranger's names for a feature barcode and its "
                        "feature. Any two-letter tag works; lowercase ones are free for local use."),
    ])


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


def recipe_nanopore(ctx):
    w = ctx.workdir("15-nanopore")
    front, rear, barcodes = native_barcodes(ctx)
    write_tsv(w / "native_barcodes.tsv", ["id", "seq"], [(i, front + s + rear) for i, s in barcodes])
    write_tsv(w / "barcode_pairs.tsv", ["sample", "front", "rear"], [(i, i, i) for i, _ in barcodes])
    origin, records = {}, []
    for barcode, run in VERD6:
        for name, seq, qual in read_fastq(ctx.fetch("sra-verd6", f"{run}.fastq")):
            records.append((name, seq, qual))
            origin[name] = barcode
    check(len(origin) == len(records), "read names are unique across the 24 runs")
    write_fastq(w / "reads.fastq.gz", records)
    tags_cmd = "head -n 3 native_barcodes.tsv\nhead -n 3 barcode_pairs.tsv"
    tags_head = ctx.sh(w, tags_cmd)
    run = ctx.sh(w, NANOPORE_COMMAND)
    counts_cmd = 'for f in demux/NB0[1-3].fastq.gz demux/unclassified.fastq.gz; do\n  echo "$f $(gzip -dc "$f" | wc -l | awk \'{print $1 / 4}\')"\ndone'
    counts = ctx.sh(w, counts_cmd)
    reads = {name: seq for name, seq, _ in records}
    correct = wrong = 0
    cut = []
    for path in sorted((w / "demux").glob("NB*.fastq.gz")):
        sample = path.name.split(".")[0]
        for name, seq, _ in read_fastq(path):
            correct += origin[name] == sample
            wrong += origin[name] != sample
            start = reads[name].find(seq)
            check(start >= 0, f"{name}: the trimmed read is part of the read")
            cut.append((start, len(reads[name]) - start - len(seq)))
    unclassified = sum(1 for _ in read_fastq(w / "demux" / "unclassified.fastq.gz"))
    check(wrong == 0, f"no read goes to another barcode ({wrong})")
    check(correct + unclassified == len(records), "every read is assigned or unclassified")
    check(all(a <= 150 and b <= 150 for a, b in cut), "the cuts fall inside the search windows")
    front_cut = sorted(a for a, _ in cut)
    rear_cut = sorted(b for _, b in cut)
    ctx.sh(w, NANOPORE_COMMAND.replace("dist=12:12:12", "dist=8:8:8").replace("demux/", "demux8/"))
    tight = collections.Counter()
    for path in sorted((w / "demux8").glob("NB*.fastq.gz")):
        for name, _, _ in read_fastq(path):
            tight["right" if origin[name] == path.name.split(".")[0] else "wrong"] += 1
    tight["unclassified"] = sum(1 for _ in read_fastq(w / "demux8" / "unclassified.fastq.gz"))
    check(sum(tight.values()) == len(records), "with 8 edits every read is assigned or unclassified")
    check(tight["wrong"] == 0 and tight["right"] < correct, "8 edits assign fewer reads and none to a wrong barcode")
    w2 = ctx.workdir("15-nanopore-mods")
    shutil.copyfile(ctx.fetch("modkit", "bc_anchored_10_reads.unmapped.bam"), w2 / "modbam.unmapped.bam")
    mods = ctx.sh(w2, MODS_COMMAND)
    kept_tags, trimmed_tags = (line.split() for line in mods.out.strip().split("\n"))
    per_base = [tag for tag in kept_tags if tag not in trimmed_tags]
    check(set(per_base) == {"MM", "ML", "mv"} and trimmed_tags == [tag for tag in kept_tags if tag not in per_base],
          "a trimming template keeps the read's tags except MM, ML and mv")
    left_out = log_lines(mods.err, ["left out per-base tags"])
    check(bool(left_out) and all(tag in left_out for tag in per_base), "unmux warns which per-base tags it left out")
    warn, listed = pair_warnings(run.err)
    check(listed > 0 and "window" not in warn, "unmux warns about tag pairs and not about the windows")
    return Recipe("Long reads", "Nanopore native barcodes at both ends", [
        "Use this to demultiplex Nanopore reads from a native barcoding kit, where each read starts "
        "with the barcode and ends with its reverse complement, trimming the barcodes off.",
        f"**Data:** 200 reads of each of 24 SRA runs from one MinION flow cell of SARS-CoV-2 ARTIC "
        f"amplicons ({link('PRJNA730334', ctx.sources[('sra-verd6', 'SRR14560737.fastq')].url)}, batch "
        "VERD6, kits SQK-LSK109 with EXP-NBD104 and EXP-NBD114). Guppy split the reads into one run "
        "per native barcode, requiring the barcode at both ends and leaving it in the reads, so the "
        f"barcode of every read is known. The {fmt(len(records))} reads are pooled one run after another:",
        fence(VERD6_FETCH, "bash"),
        f"Each tag is a native barcode between its flanks, from dorado "
        f"({source_links(ctx, ('dorado', 'barcode_kits.cpp'))}), and the sheet asks for the same "
        "barcode at both ends:",
        fence(tags_cmd, "bash"),
        fence(tags_head.out, "tsv"),
        fence(NANOPORE_COMMAND, "bash"),
        "The `front` group searches the first 150 bases and the `rear` group the last 150 on both "
        "strands, so it finds the reverse-complemented copy. A budget of 12 edits in a 46 bp tag "
        "suits the error rate of these R9.4.1 reads, and `mode=nearest` leaves a read that two "
        "barcodes fit equally well unclassified. The span `@front..@rear` keeps the bases between "
        "the two tags.",
        fence(counts_cmd, "bash"),
        fence(counts.out, "text"),
        f"With a budget this large, unmux warns about pairs of tags that one read could match, "
        f"listing {listed} for each group and counting the rest:",
        fence(warn, "text"),
        "The windows are much wider than the tags, as long reads need, and unmux prints no window "
        "warning: it warns about one only where it estimates that over 0.3% of reads with an error "
        "could match another tag at a shifted position.",
        f"{fmt(correct)} of {fmt(len(records))} reads ({pct(correct, len(records))}) went to the barcode "
        f"of their run and {fmt(wrong) if wrong else 'none'} to another; {fmt(unclassified)} are "
        "unclassified. Every trimmed read is part of its read, cut within the search windows: a "
        f"median of {front_cut[len(front_cut) // 2]} bases off the front and "
        f"{rear_cut[len(rear_cut) // 2]} off the rear.",
        callout("TIP", f"With `dist=8:8:8`, {fmt(tight['right'])} reads go to the barcode of their run, "
                       f"{fmt(tight['wrong']) if tight['wrong'] else 'none'} to another, and "
                       f"{fmt(tight['unclassified'])} are unclassified: a tighter budget costs reads "
                       "here without saving any from a wrong barcode."),
        "What a trim does to the basecaller's tags, on the ten modkit reads of "
        f"{ref('comment_tags')}, before and after a 20-base trim:",
        fence(MODS_COMMAND, "bash"),
        fence(mods.out, "text"),
        fence(left_out, "text"),
        callout("WARNING", "A `--template` record keeps the tags of its read, but per-base tags such as "
                           "MM, ML and the move table mv describe the read's bases one by one, so unmux "
                           "leaves them out of a record that does not copy its read exactly, with one "
                           "warning for the run. Demultiplex without a template to keep the "
                           "modification calls; the barcodes then stay in the reads."),
    ])


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


def recipe_pacbio(ctx):
    w = ctx.workdir("16-pacbio")
    for name in ("prefix.bam", "biosamples.csv", "Sequel_96_barcodes_v1.fasta"):
        shutil.copyfile(ctx.fetch("pacbio-cyp2d6", name), w / name)
    ctx.sh(w, PACBIO_PREP)
    heads_cmd = "head -n 3 barcodes.tsv\nhead -n 3 barcode_pairs.tsv"
    heads = ctx.sh(w, heads_cmd)
    tag_ids = {row["seq"]: row["id"] for row in read_tsv(w / "barcodes.tsv")}
    tag_len = {len(seq) for seq in tag_ids}
    check(len(tag_len) == 1, "every tag has one length")
    tag_len = tag_len.pop()
    lima_version = ctx.sh(w, "lima --version").out.split("\n")[0].strip()
    check(re.fullmatch(r"lima \d+\.\d+\.\d+", lima_version) is not None, f"lima reports its version ({lima_version!r})")
    ctx.sh(w, LIMA_COMMAND)
    run = ctx.sh(w, PACBIO_COMMAND)
    with open(w / "biosamples.csv", newline="") as fh:
        pairs = {}
        for row in csv.DictReader(fh):
            a, b = row["Barcode"].split("--")
            pairs[(a, b)] = pairs[(b, a)] = row["Bio Sample Name"]
    samples = len(set(pairs.values()))
    report = {row["ZMW"].split("/")[1]: row for row in read_tsv(w / "lima.lima.report")}
    lima = {z: pairs.get((row["IdxFirstNamed"], row["IdxCombinedNamed"])) if row["PassedFilters"] == "1" else None
            for z, row in report.items()}
    inputs = sam_records(ctx, w, "hifi_reads.bam", "--no-PG")
    reads = {f[0].split("/")[1]: f[9] for f, _ in inputs}
    carried = {f[0].split("/")[1]: {k: v for k, v in tags.items() if k != "RG"} for f, tags in inputs}
    check(len(reads) == len(inputs), "one record per ZMW")
    unmux_calls, cuts = {}, []
    for path in sorted(p for p in (w / "unmux").glob("*.bam") if p.name != "unassigned.bam"):
        for fields, tags in sam_records(ctx, w, f"unmux/{path.name}"):
            z = fields[0].split("/")[1]
            unmux_calls[z] = path.name.split(".")[0]
            start = reads[z].find(fields[9])
            check(start >= 0, f"{fields[0]}: the trimmed read is part of the read")
            cuts.append((start, len(reads[z]) - start - len(fields[9])))
            check({k: v for k, v in tags.items() if k != "RG"} == carried[z] and tags["RG"] == unmux_calls[z],
                  f"{fields[0]}: the record keeps its read's tags, with the sample's read group")
    check(all(a <= 60 and b <= 60 for a, b in cuts), "the cuts fall inside the search windows")
    table = collections.Counter()
    for z in reads:
        a, b = lima.get(z), unmux_calls.get(z)
        table[("same" if a == b else "different") if a and b else "lima only" if a else "unmux only" if b else "neither"] += 1
    check(table["different"] == 0, f"lima and unmux never call different samples ({table['different']})")
    lima_n = sum(1 for z in reads if lima.get(z))
    unmux_only = [z for z in reads if unmux_calls.get(z) and not lima.get(z)]
    lima_only = [z for z in reads if lima.get(z) and not unmux_calls.get(z)]
    neither = [z for z in reads if not lima.get(z) and not unmux_calls.get(z)]
    check(all(z in report and report[z]["PassedFilters"] != "1" for z in unmux_only), "lima removes every read only unmux assigns")
    scores = sorted(int(report[z]["ScoreCombined"]) for z in unmux_only)
    check(scores[-1] <= 80, "every read only unmux assigns has a lima combined score of 80 or less")
    ctx.sh(w, PACBIO_COMMAND.replace('--out "unmux/%sample.bam"', '--qc-tag \\\n  --out "qc/%sample.bam"').replace(
        "--unassigned unmux/unassigned.bam", "--unassigned qc/unassigned.bam"))
    qc = {f[0].split("/")[1]: json.loads(t["ZS"]) for f, t in sam_records(ctx, w, "qc/unassigned.bam")}
    check(set(qc) == set(lima_only) | set(neither), "unmux leaves out exactly the reads it does not assign")

    def ends(z):
        found = {g["g"] for g in qc[z]["groups"]}
        return "both" if {"start", "end"} <= found else "one" if found else "none"

    sides = collections.Counter(ends(z) for z in lima_only)
    left = collections.Counter(ends(z) for z in neither)
    lima_short = sum(1 for z in neither if z in report and int(report[z]["NumBarcodedRegionsPassed"]) < 2)
    twice = unlisted = 0
    for z in (z for z in neither if ends(z) == "both"):
        found = {g["g"]: tag_ids[g["tag"]] for g in qc[z]["groups"]}
        row = report[z]
        check({found["start"], found["end"]} == {row["IdxFirstNamed"], row["IdxCombinedNamed"]},
              f"{z}: lima reports the barcodes unmux found")
        check((found["start"], found["end"]) not in pairs, f"{z}: the barcodes are not a sample's pair")
        twice += found["start"] == found["end"]
        unlisted += found["start"] != found["end"]
    missing = [z for z in neither if z not in report]
    check(all(z in report for z in reads if z not in missing), "lima's report has every other read")
    warn, listed = pair_warnings(run.err)
    check(listed > 0 and "unassigned" in warn and "per-base" not in run.err,
          "unmux warns about tag pairs and the unassigned share, and leaves out no per-base tag")
    lima_out = sorted(w.glob("lima.bc*.bam"))
    lima_tags = sam_records(ctx, w, lima_out[0].name)[0][1]
    pacbio_tags = sorted(t for t in inputs[0][1] if t != "RG")
    check(all(t in lima_tags for t in pacbio_tags), "lima keeps the input tags")
    front = sorted(a for a, _ in cuts)[len(cuts) // 2]
    rear = sorted(b for _, b in cuts)[len(cuts) // 2]
    rows = "\n".join([
        "| Call | Reads |",
        "| --- | ---: |",
        f"| Same sample from lima and unmux | {fmt(table['same'])} |",
        f"| A sample from lima only | {fmt(table['lima only'])} |",
        f"| A sample from unmux only | {fmt(table['unmux only'])} |",
        f"| No sample from either | {fmt(table['neither'])} |",
    ])
    modes = "\n".join([
        "| Setting | lima | unmux |",
        "| --- | --- | --- |",
        f"| The same barcode at both ends | `--hifi-preset SYMMETRIC`, which sets `--same` | a sheet row with that barcode in both groups, as in {ref('nanopore')} |",
        "| A different barcode at each end | `--hifi-preset ASYMMETRIC`, which sets `--different` | a sheet row for the pair in each order, since a read can come from either strand |",
        "| The reverse-complemented barcode at the far end | searched by default | `bothStrands=true` on the end group |",
        "| Where to search | `-W` in bases, or `-w` as a multiple of the barcode length | `loc=0:0:N` for the start and `loc=0:-N:end` for the end |",
        "| How well a barcode must fit | `--min-score` and `--min-end-score`, alignment scores out of 100 | `dist`, a budget of substitutions, indels and edits in all |",
        "| How far the best must lead | `--min-score-lead` | `mode=nearest`, with `delta` for a larger lead |",
        "| Sample names and files | `--biosample-csv`, and `--split-named` for a file per barcode pair | `--sample-sheet`, and `--out \"%sample.bam\"` for a file per sample |",
        "| Barcode removal | clipped by default, kept with `--no-clip` | cut with `--template` and `@start..@end`, kept with no template |",
        "| Read filters | `--min-ref-span`, `--min-scoring-regions` and `--min-length` | no equivalent |",
        "| Finding which barcodes are in use | `--peek-guess` | no equivalent: the sheet lists them |",
        "| Per-read reports | `bc` and `bq` tags, and the `.lima.report` file | `--qc-tag` JSON, and the metrics files |",
    ])
    lima_sources = source_links(ctx, ("pacbio-cyp2d6", "prefix.bam"), ("pacbio-cyp2d6", "biosamples.csv"),
                                ("pacbio-cyp2d6", "Sequel_96_barcodes_v1.fasta"))
    missing_text = ("" if not missing else
                    f"One read, of {len(reads[missing[0]])} bases, is missing from lima's report." if len(missing) == 1 else
                    f"{fmt(len(missing))} reads are missing from lima's report.")
    return Recipe("Long reads", "PacBio barcode pairs, compared with lima", [
        "Use this to demultiplex PacBio HiFi reads that still carry their barcodes, one at each end, "
        "as lima does, here with a different barcode at each end of a sample's amplicons.",
        f"**Data:** the first {fmt(len(reads))} HiFi reads of PacBio's CYP2D6 amplicon set, {samples} "
        "samples each with a pair of barcodes, before lima ran on them, cut from the first 10 MB of "
        f"the BAM by an HTTP range request ({lima_sources}). Each read starts with a pad of a few "
        "bases, a barcode and the M13 forward or reverse primer, and ends with the same for the "
        "other barcode, reverse-complemented.",
        fence(CYP2D6_FETCH, "bash"),
        "Keep the complete records of the cut, and build a tag of each barcode with its M13 primer "
        "and a sheet with each pair in both orders:",
        fence(PACBIO_PREP, "bash"),
        fence(heads_cmd, "bash"),
        fence(heads.out, "tsv"),
        f"{lima_version}, with its HiFi preset for a different barcode at each end:",
        fence(LIMA_COMMAND, "bash"),
        "and unmux:",
        fence(PACBIO_COMMAND, "bash"),
        f"unmux warns about pairs of tags that one read could match within this budget, listing "
        f"{listed} for each group, and about the share of reads it leaves unassigned:",
        fence(warn, "text"),
        f"Per read, against lima's report ({fmt(lima_n)} reads with a sample):",
        rows,
        "The two never call different samples for one read. lima scores each barcode by alignment, "
        "and its ASYMMETRIC preset removes a read whose combined score is 80 or less: the "
        f"{fmt(len(unmux_only))} reads only unmux assigns score from {scores[0]} to {scores[-1]}, though "
        "both of their ends fit one sheet row within unmux's budget. unmux needs a tag within its "
        f"budget at both ends, and in the {fmt(len(lima_only))} reads only lima assigns, its "
        f"`--qc-tag` found a tag {tally([(sides['one'], 'at one end only'), (sides['both'], 'at both ends'), (sides['none'], 'at neither end')], len(lima_only))}.",
        f"Of the {fmt(len(neither))} reads neither assigns, unmux found a tag "
        f"{tally([(left['one'], 'at one end only'), (left['none'], 'at neither end')])}, and lima "
        f"passed fewer than two barcode regions in {fmt(lima_short)}. The other {fmt(left['both'])} "
        "carry tags at both ends but not a sample's pair: "
        f"{tally([(twice, 'the same barcode twice'), (unlisted, 'a pair no sample has')])}, and lima "
        f"reports the same barcodes. {missing_text}",
        f"Every trimmed read is part of its read, cut within the search windows: a median of {front} "
        f"bases off the start and {rear} off the end, each a {tag_len}-base tag and its pad.",
        "How unmux expresses lima's settings:",
        modes,
        callout("NOTE", f"Both tools keep the tags PacBio wrote ({', '.join(pacbio_tags)}) on the reads "
                        "they clip, and lima adds its own. unmux sets RG to the sample's read group."),
    ])


def recipe_direct_rna(ctx):
    w = ctx.workdir("14-direct-rna")
    shutil.copyfile(ctx.fetch("dorado", "fastq_with_us_and_tags.fq"), w / "direct_rna.fastq")
    head_cmd = "head -n 2 direct_rna.fastq | cut -c 1-100"
    head = ctx.sh(w, head_cmd)
    command = "unmux direct_rna.fastq --out direct_rna.bam"
    ctx.sh(w, command)
    view_cmd = "samtools view direct_rna.bam | cut -f 10 | cut -c 1-60\nsamtools view direct_rna.bam | cut -f 1,2,12-"
    view = ctx.sh(w, view_cmd)
    text_cmd = (
        "unmux direct_rna.fastq --out direct_rna.sam\n"
        "grep -v '^@' direct_rna.sam | cut -f 10 | cut -c 1-40\n"
        "unmux direct_rna.fastq --out copy.fastq\n"
        "head -n 2 copy.fastq | tail -n 1 | cut -c 1-40"
    )
    text = ctx.sh(w, text_cmd)
    seq_in = next(read_fastq(w / "direct_rna.fastq"))[1]
    rec = sam_records(ctx, w, "direct_rna.bam")[0]
    assert rec[0][9] == seq_in.replace("U", "T")
    ctx.sh(w, "samtools import -0 direct_rna.fastq -o samtools.bam")
    check(sam_records(ctx, w, "samtools.bam")[0][0][9] == rec[0][9], "BAM bases equal samtools import's")
    probe = revcomp(seq_in[100:124].replace("U", "T"))
    ctx.sh(w, f'unmux direct_rna.fastq --group "x={{dna={probe}}}" --group "x::bothStrands=true" '
              '--sample-from-group x --out "probe.%sample.fastq"')
    check(len(list(read_fastq(w / "probe.dna.fastq"))) == 1, "a T-only tag from the other strand finds the U read")
    u = seq_in.count("U")
    return Recipe("Long reads", "Direct RNA reads with U bases to BAM", [
        "Use this to store direct RNA reads, which basecallers can write with U for uracil, in an "
        "unmapped BAM.",
        f"**Data:** a dorado test read written with U, with SAM tags in its FASTQ comment "
        f"({source_links(ctx, ('dorado', 'fastq_with_us_and_tags.fq'))}).",
        fence(head_cmd, "bash"),
        fence(head.out, "text"),
        fence(command, "bash"),
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        f"BAM has no code for U, so all {fmt(u)} U bases are stored as T, as `samtools import` stores "
        f"them. The tags in the comment become SAM tags (more in {ref('comment_tags')}), and RG:Z is "
        "replaced by the pool, as for any input.",
        "SAM and FASTQ output keep the U as read:",
        fence(text_cmd, "bash"),
        fence(text.out, "text"),
        callout("TIP", "When matching, U in a read or a tag counts as T, so DNA barcodes find RNA "
                       "reads written with U, on either strand. The generator checks this with a "
                       "24-mer of this read's reverse complement written with T."),
    ])


COMMENT_PREP = "samtools fastq -T '*' modbam.unmapped.bam > modbam.fastq"


def recipe_comment_tags(ctx):
    w = ctx.workdir("15-comment-tags")
    shutil.copyfile(ctx.fetch("modkit", "bc_anchored_10_reads.unmapped.bam"), w / "modbam.unmapped.bam")
    prep_cmd = COMMENT_PREP + "\nhead -n 1 modbam.fastq | tr '\\t' '\\n' | cut -c 1-60"
    header = ctx.sh(w, prep_cmd)
    command = "unmux modbam.fastq --out modbam.roundtrip.bam"
    ctx.sh(w, command)
    view_cmd = "samtools view modbam.roundtrip.bam | head -n 1 | cut -f 1,2,12- | tr '\\t' '\\n' | cut -c 1-60"
    view = ctx.sh(w, view_cmd)
    original = sam_records(ctx, w, "modbam.unmapped.bam")
    back = sam_records(ctx, w, "modbam.roundtrip.bam")
    assert len(original) == len(back)
    for (fa, ta), (fb, tb) in zip(original, back):
        assert fa[0] == fb[0] and fa[9] == fb[9] and fa[10] == fb[10]
        assert {k: v for k, v in ta.items() if k != "RG"} == {k: v for k, v in tb.items() if k != "RG"}
        assert tb["RG"] == "modbam" != ta["RG"]
    ctx.sh(w, "unmux modbam.unmapped.bam --out direct.bam")
    assert all("NM" not in t and "MD" not in t for _, t in sam_records(ctx, w, "direct.bam"))
    return Recipe("Long reads", "Lift SAM tags from FASTQ comments back onto records", [
        "Use this when reads went through FASTQ with their SAM tags in the name comment, as "
        "`samtools fastq -T` and basecallers write them, and the tags, MM and ML among them, should "
        "be SAM tags again.",
        f"**Data:** ten unmapped Nanopore reads with modified-base calls "
        f"({source_links(ctx, ('modkit', 'bc_anchored_10_reads.unmapped.bam'))}), written to FASTQ "
        "with every tag in the comment:",
        fence(prep_cmd, "bash"),
        fence(header.out, "text"),
        fence(command, "bash"),
        fence(view_cmd, "bash"),
        fence(view.out, "text"),
        f"On all {fmt(len(back))} reads the bases, qualities and every tag, MM and ML included, equal "
        "the original BAM's, except RG:Z, which becomes the pool as for any input. Even NM and MD "
        "come back: tags lifted from a FASTQ comment are kept as they are, while the same tags on a "
        "SAM or BAM input record are dropped as alignment tags.",
        callout("NOTE", "A comment that is not a SAM tag, such as Illumina's `1:N:0:ACGT`, is not "
                        "kept, and a malformed tag stops the run."),
    ])


def recipe_qc_tag(ctx):
    w = ctx.workdir("16-qc-tag")
    source = ctx.cache / "work" / "06-bins"
    for name in ("R1.fastq.gz", "R2.fastq.gz", "I1.fastq.gz", "I2.fastq.gz", "samples.tsv", "i7.tsv", "i5.tsv", "phix.tsv"):
        shutil.copyfile(source / name, w / name)
    command = BINS_COMMAND.replace("  --metrics-per-sample metrics.per_sample.tsv \\\n  --metrics-summary metrics.summary.tsv \\\n  --metrics-unassigned metrics.unassigned.tsv", "  --qc-tag")
    ctx.sh(w, command)
    barcodes = {row["id"]: row["seq"] for row in read_tsv(w / "i7.tsv")}
    fixed = sample = None
    i1 = {n.split()[0]: s for n, s, _ in read_fastq(w / "I1.fastq.gz")}
    for path in sorted((w / "fastq").glob("N*.R1.fastq.gz")):
        i7 = barcodes[path.name.split("_")[0]]
        hit = next((n.split()[0] for n, _, _ in read_fastq(path) if hamming(i1[n.split()[0]], i7) == 1), None)
        if hit:
            fixed, sample = hit, path.name.split(".")[0]
            break
    other = next(n.split()[0] for n, _, _ in read_fastq(w / "fastq" / "Undetermined.0.fastq.gz"))
    removed = [n.split()[0] for n, _, _ in read_fastq(w / "removed" / "PhiX.0.fastq.gz")]
    phix = next(n for n in removed if n.startswith("A00488"))
    qc = {h.split()[0]: json.loads(h.split("ZS:Z:", 1)[1]) for h, _, _ in read_fastq(w / "removed" / "PhiX.0.fastq.gz")}
    stray = next(n for n in removed if n.startswith("M03958") and all(g["sub"] > 0 for g in qc[n]["groups"] if g["g"] == "phix"))
    show = (
        f"gzip -dc fastq/{sample}.R1.fastq.gz | grep '^@{fixed} '\n"
        f"gzip -dc fastq/Undetermined.0.fastq.gz | grep '^@{other} '\n"
        f"gzip -dc removed/PhiX.0.fastq.gz | grep '^@{phix} '\n"
        f"gzip -dc removed/PhiX.0.fastq.gz | grep '^@{stray} '"
    )
    out = ctx.sh(w, show)
    lines = out.out.splitlines()
    check(len(lines) == 4, "one header line per read shown")
    outcomes = [json.loads(line.split("ZS:Z:", 1)[1])["outcome"] for line in lines]
    check(outcomes == ["assigned", "unassigned", "removed", "removed"], f"QC outcomes {outcomes}")
    check(json.loads(lines[0].split("ZS:Z:", 1)[1])["sample"] == sample, "the assigned read's QC sample")
    matched = [g["g"] for g in json.loads(lines[1].split("ZS:Z:", 1)[1])["groups"]]
    why = {(): "matched neither index", ("i7",): "matched an i7 of the sheet but not its i5"}[tuple(matched)]
    return Recipe("Diagnostics", "Explain one read's assignment with --qc-tag", [
        "Use this to see why a read went where it did. The `--qc-tag` option writes a JSON summary "
        "of the match into a tag, ZS unless another is given, in the FASTQ comment or on the SAM "
        "record.",
        f"**Data:** the pooled lane of {ref('bins')}.",
        fence(command, "bash"),
        "A fig2b read whose i7 had one mismatch, a fig2c read, a PhiX read, and the fig2b read the "
        "PhiX rule took with one substitution:",
        fence(show, "bash"),
        fence(out.out, "text"),
        "Each matched group lists the tag it matched, where (`loc` is file, start and end), the "
        f"observed bases, and the substitutions and indels it took. The fig2c read {why}, so no "
        "sample claims it and its reason is `no_sample`; a removed read gives the `--remove` rule "
        "that took it, and `str` marks a match on the reverse strand. The last record is how to "
        "find out why a read was removed by mistake. The generator checks each record's outcome "
        "against the file it is in.",
    ])


def recipe_metrics(ctx):
    w = ctx.cache / "work" / "06-bins"
    summary = ctx.sh(w, "cat metrics.summary.tsv")
    per = ctx.sh(w, "head -n 4 metrics.per_sample.tsv")
    unassigned = ctx.sh(w, "cat metrics.unassigned.tsv")
    s = read_tsv(w / "metrics.summary.tsv")[0]
    rows = read_tsv(w / "metrics.per_sample.tsv")
    reasons = read_tsv(w / "metrics.unassigned.tsv")
    assigned = sum(int(r["reads"]) for r in rows)
    un = sum(int(r["reads"]) for r in reasons)
    total = int(s["total_reads"])
    check(assigned == int(s["reads_assigned"]) and un == int(s["reads_unassigned"]), "per-sample and reason rows add up")
    check(assigned + un + int(s["reads_removed"]) + int(s["reads_pass_through"]) == total, "the four outcomes add up")
    pool = sum(1 for _ in read_fastq(w / "R1.fastq.gz"))
    check(total == pool, "total_reads is the number of pairs in the pool")
    for row in rows:
        written = sum(1 for _ in read_fastq(w / "fastq" / f"{row['sample']}.R1.fastq.gz"))
        check(int(row["reads"]) == written, f"{row['sample']}: per-sample reads equal the pairs in its files")
    removed = sum(1 for _ in read_fastq(w / "removed" / "PhiX.0.fastq.gz"))
    undetermined = sum(1 for _ in read_fastq(w / "fastq" / "Undetermined.0.fastq.gz"))
    check(int(s["reads_removed"]) == removed and int(s["reads_unassigned"]) == undetermined, "bins match the summary")
    return Recipe("Diagnostics", "Read the metrics files", [
        "Use the three metrics files to check a run: how many reads each sample got, how many "
        "matched no sample and why, and the totals.",
        f"**Data:** the pooled lane of {ref('bins')}, whose command wrote all three.",
        fence("cat metrics.summary.tsv", "bash"),
        fence(summary.out, "tsv"),
        "A read is a template: one pair here. Every read is assigned, passed through (when no "
        "sample is declared), unassigned or removed, and `frac_assigned` is the assigned share. "
        "The `mean_read_len_by_input` column gives the mean length of each input in order, which "
        "catches an index read shorter than its tags.",
        fence("head -n 4 metrics.per_sample.tsv", "bash"),
        fence(per.out, "tsv"),
        "The `frac_bases_unextracted` column is the share of a sample's input bases left out of its "
        "output, here the index reads.",
        fence("cat metrics.unassigned.tsv", "bash"),
        fence(unassigned.out, "tsv"),
        f"The per-sample reads add up to `reads_assigned` ({fmt(assigned)}), the reasons to "
        f"`reads_unassigned` ({fmt(un)}), and with the {fmt(int(s['reads_removed']))} removed they "
        f"make up all {fmt(total)} pairs of the pool; each sample's count, the removed count and the "
        "unassigned count also equal the pairs in its files. A reason of `no_sample` means no sample "
        "claimed what the groups matched; a group with `minFindsPerGroup` adds rows that give the "
        "group and the bound it missed.",
    ])


def rounds(command, old, new):
    for r in ("round1", "round2", "round3"):
        command = command.replace(f'"{r}::{old}', f'"{r}::{new}')
    return command


INDEL_COMMAND = rounds(SPLITSEQ_COMMAND, 'dist=1"', 'dist=1:1:2"')
WIDE_ALL_COMMAND = rounds(
    INDEL_COMMAND.replace("loc=1:78:86", "loc=1:77:87").replace("loc=1:48:56", "loc=1:47:57").replace("loc=1:10:18", "loc=1:9:19"),
    'dist=1:1:2"', 'dist=1:1:2,mode=all"',
)
WIDE_NEAREST_COMMAND = rounds(WIDE_ALL_COMMAND, "dist=1:1:2,mode=all", "dist=1:1:2,mode=nearest")


def recipe_collision(ctx):
    w = ctx.workdir("17-collision-10x")
    tenx_inputs(ctx, w)
    shutil.copyfile(ctx.cache / "work" / "08-10x" / "cells.tsv", w / "cells.tsv")
    command = TENX_COMMAND.replace("dist=1,mode=nearest", "dist=1")
    error = refused(ctx, w, command, "with mode=all this silently mis-assigns reads", snapshot(w))
    w2 = ctx.workdir("17-collision-splitseq")
    splitseq_inputs(ctx, w2)
    error2 = refused(ctx, w2, INDEL_COMMAND, "only through an indel", snapshot(w2))
    return Recipe("Diagnostics", "The tag-collision check", [
        "Before it reads a record, unmux checks the tags of each group against each other. A set in "
        "which one read could match two tags equally well is refused, with the ways out.",
        f"**Data:** the 10x and SPLiT-seq inputs from {ref('tenx')} and {ref('splitseq')}.",
        "The 10x command with `dist=1` and no mode. Cell barcodes are not designed to be far apart, "
        "so some differ in only two places, and a read with one error between them matches both:",
        fence(command, "bash"),
        fence(error, "text"),
        "Under the default `mode=all` the first match would silently win. With `mode=nearest`, as in "
        "the 10x recipe, such reads are left uncorrected and the check only warns, listing up to ten "
        "pairs.",
        "The SPLiT-seq command with one indel allowed (`dist=1:1:2`). Its barcodes are 4 "
        "substitutions apart, but an insertion or deletion can bring a read within reach of two:",
        details("The command", fence(INDEL_COMMAND, "bash")),
        fence(error2, "text"),
        "Here unmux asks for a mode. With `mode=nearest` such reads are left unassigned, and with "
        "`mode=all` each is still assigned by a tie-break, with warnings that list the pairs; "
        f"{ref('wide')} compares the two. For both commands the generator checks the exit status, "
        "the message, and that no file was written.",
        callout("TIP", "For a set designed to be far apart, such as Illumina indexes, a refusal "
                       "usually means `dist` is too large for the set."),
    ])


def recipe_wide(ctx):
    w = ctx.workdir("18-wide-windows")
    splitseq_inputs(ctx, w)
    cells = called_cells(ctx, w)
    base = ctx.cache / "work" / "09-splitseq"
    reference = {f[0]: t.get("CB") for f, t in sam_records(ctx, base, "splitseq.unmapped.bam")}
    rows = [("README windows, `dist=1`", splitseq_scores(ctx, base, "splitseq.unmapped.bam", cells))]
    run = ctx.sh(w, WIDE_ALL_COMMAND)
    rows.append(("1 bp wider, `dist=1:1:2,mode=all`", splitseq_scores(ctx, w, "splitseq.unmapped.bam", cells, reference)))
    run2 = ctx.sh(w, WIDE_NEAREST_COMMAND)
    rows.append(("1 bp wider, `dist=1:1:2,mode=nearest`", splitseq_scores(ctx, w, "splitseq.unmapped.bam", cells, reference)))
    table = "\n".join([
        "| Command | Written | With a CB | CB of a called cell | Other CB | No CB | CB changed from the README windows |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ] + [f"| {label} | {fmt(n)} | {fmt(cb)} | {fmt(inc)} | {fmt(cb - inc)} | {fmt(n - cb)} | "
         + ("" if label.startswith("README") else fmt(chg)) + " |" for label, (n, cb, inc, chg) in rows])
    (_, readme), (_, wide_all), (_, wide_nearest) = rows
    check(readme[0] == readme[1] and wide_all[0] == wide_all[1], "every written read has a CB except under mode=nearest")
    missing = wide_nearest[0] - wide_nearest[1]
    check(missing > 0, "mode=nearest writes reads without a CB")
    count_cmd = "samtools view -c splitseq.unmapped.bam\nsamtools view -c -d CB splitseq.unmapped.bam"
    counted = ctx.sh(w, count_cmd)
    check(counted.out.split() == [str(wide_nearest[0]), str(wide_nearest[1])], "samtools counts match")
    qc = WIDE_NEAREST_COMMAND.replace('--out "splitseq.unmapped.bam"', '--qc-tag \\\n  --out "qc.unmapped.bam"')
    ctx.sh(w, qc)
    tied = [json.loads(t["ZS"]) for _, t in sam_records(ctx, w, "qc.unmapped.bam") if "CB" not in t]
    check(len(tied) == missing and all(len(z["groups"]) < 3 for z in tied), "QC of a read without CB lists fewer than three rounds")
    warned = [line for line in run.err.split("\n") if "wider than" in line]
    pairs = len([line for line in run.err.split("\n") if "can both match" in line])
    quiet = "wider than" not in run2.err
    oracle = splitseq_oracle(w)
    return Recipe("Diagnostics", "Search windows wider than their tags", [
        "Widening a window lets a barcode that an insertion or deletion shifted still be found. But "
        "in a window wider than its tags, a read with an error in its barcode can also match another "
        "tag, shifted, as closely as its own. unmux estimates how often that happens and warns.",
        f"**Data:** the SPLiT-seq reads of {ref('splitseq')}, of which "
        f"{fmt(sum(1 for v in oracle.values() if not v))} miss a round in the README windows.",
        "The README command with each window 1 bp wider on both sides, one indel allowed, and "
        f"`mode=all` as the collision check of {ref('collision')} asks:",
        fence(WIDE_ALL_COMMAND, "bash"),
        fence("\n".join(warned), "text"),
        f"The run also prints {fmt(pairs)} warnings that list tag pairs one indel apart, as `mode=all` "
        "does. Following the advice, the same command with `mode=nearest` in each round:",
        details("The command", fence(WIDE_NEAREST_COMMAND, "bash")),
        ("It prints no window warning: under `mode=nearest`, a window wider than the tags by twice "
         "the indel budget is as wide as unmux advises."
         if quiet else fence(log_lines(run2.err, ["wider than"]), "text")),
        "Real reads have no known truth, so the table scores each run two ways: against the 568 "
        "cells split-pipe called, whose barcodes a correct CB most likely belongs to, and against "
        "the README windows, which take no shifted match:",
        table,
        f"The wider windows find {fmt(wide_all[2] - readme[2])} more reads of called cells, but under "
        f"`mode=all` they also change the CB of {fmt(wide_all[3])} reads that the README windows "
        f"assign, and add {fmt((wide_all[1] - wide_all[2]) - (readme[1] - readme[2]))} reads whose CB "
        f"is no called cell's. Under `mode=nearest` fewer CBs change ({fmt(wide_nearest[3])}), but "
        f"{fmt(missing)} reads that tie in some round are written with no CB at all:",
        fence(count_cmd, "bash"),
        fence(counted.out, "text"),
        callout("WARNING", "Under `mode=nearest`, a read that ties between two tags of a round is still "
                           "written, without the CB built from that round, though every round sets "
                           "`minFindsPerGroup=1`; its `--qc-tag` lists only the rounds that matched. "
                           "Keep only reads with a CB, for example with `samtools view -d CB`."),
    ])


def recipe_refuse_input(ctx):
    w = ctx.workdir("19-refuse-input")
    place(ctx.fetch("sra-SRR390728", "SRR390728_1.fastq"), w / "SRR390728_1.fastq.gz")
    command = 'unmux SRR390728_1.fastq.gz \\\n  --extract "trimmed=0:10:end" --template trimmed \\\n  --out "%pool.fastq.gz"'
    error = refused(ctx, w, command, 'output "SRR390728_1.fastq.gz" (--out) is input 0', snapshot(w))
    return Recipe("Refusals that protect data", "An output over an input", [
        f"The pool placeholder expands to the input's stem, so this command on the reads of "
        f"{ref('fastq_to_ubam')}, meant to write a trimmed copy, would write over its own input:",
        fence(command, "bash"),
        fence(error, "text"),
        "unmux compares every output, bins and metrics files included, with every file it reads "
        "under any name (inputs, tag files, the sample sheet, the reference), before it writes "
        "anything. The generator checks the exit status, the message, and that the input is "
        "unchanged.",
    ])


def refusal_inputs(ctx, w):
    source = ctx.cache / "work" / "05-uditas"
    for name in ("Undetermined_S0_L001_R1_001.fastq.gz", "Undetermined_S0_L001_R2_001.fastq.gz",
                 "Undetermined_S0_L001_I1_001.fastq.gz", "Undetermined_S0_L001_I2_001.fastq.gz",
                 "samples.tsv", "i7.tsv", "i5.tsv"):
        shutil.copyfile(source / name, w / name)
    return snapshot(w)


def snapshot(w):
    return {p.relative_to(w).as_posix(): sha256(p) if p.is_file() else "directory" for p in sorted(w.rglob("*"))}


def refused(ctx, w, command, expected, before):
    run = ctx.sh(w, command, ok=(1,))
    error = log_lines(run.err, ["ERROR"])
    check(expected in error, f"the refusal says: {expected}")
    check(snapshot(w) == before, "the refused run wrote and changed no file")
    return error


def recipe_refuse_tag(ctx):
    w = ctx.workdir("20-refuse-tag")
    before = refusal_inputs(ctx, w)
    command = UDITAS_COMMAND.replace(
        '--extract "umi=3:8:end" \\\n', '--extract "umi=3:8:end" \\\n  --extract "b7=2:0:8" --extract "b5=3:0:8" \\\n').replace(
        '  --tag "RX=umi" \\\n', '  --tag "BC=b7+b5" \\\n  --tag "RX=umi" --tag "RX::qual=QT" \\\n')
    error = refused(ctx, w, command, "tag `QT` would be written twice", before)
    return Recipe("Refusals that protect data", "Two writers for one SAM tag", [
        f"The UDiTaS command of {ref('uditas')}, with the index bases added as BC and the UMI "
        "qualities sent to QT by mistake:",
        fence(command, "bash"),
        fence(error, "text"),
        "QT is already the default quality tag of BC, so one would overwrite the other. Each "
        "binding, its quality tag and `--qc-tag` must write a different tag, and RG is reserved "
        "for the read group. The one exception is CB and CR, which share their default CY as the "
        "SAM spec pairs them. The generator checks the exit status, the message, and that no file "
        "changed.",
    ])


REFUSE_BIN_COMMAND = """unmux Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz \\
  Undetermined_S0_L001_I1_001.fastq.gz Undetermined_S0_L001_I2_001.fastq.gz \\
  --group "i7=i7.tsv" --group "i7::loc=2:0:8,dist=1" \\
  --group "i5=i5.tsv" --group "i5::loc=3:0:8,dist=1" \\
  --extract "r1=0:0:end" --extract "r2=1:0:end" \\
  --template r1 r2 \\
  --sample "N701_N504=i7::N701+i5::N504" \\
  --out "N701_N504.R%ordinal.fastq.gz" \\
  --unassigned "N701_N504.R%source.fastq.gz\""""


def recipe_refuse_bin(ctx):
    w = ctx.workdir("21-refuse-bin")
    before = refusal_inputs(ctx, w)
    error = refused(ctx, w, REFUSE_BIN_COMMAND, 'outputs "N701_N504.R1.fastq.gz" (sample N701_N504) and '
                                                '"N701_N504.R1.fastq.gz" (--unassigned) are the same file', before)
    return Recipe("Refusals that protect data", "A bin that shares a sample's file", [
        f"Pulling one sample out of the UDiTaS run of {ref('uditas')} and setting the rest beside it, "
        "with names that look alike:",
        fence(REFUSE_BIN_COMMAND, "bash"),
        fence(error, "text"),
        "The `%ordinal` placeholder counts reads of a template from 1, while `%source` counts input "
        "files from 0, so the unassigned reads of input 1 (R2) would be written into the sample's R1 "
        "file. unmux expands every pattern for every sample and refuses the run before writing; the "
        "generator checks the exit status, the message, and that no file changed.",
        callout("TIP", "To keep unassigned reads in the sample's own files on purpose, pass the "
                       "`--out` pattern itself: `--unassigned \"N701_N504.R%ordinal.fastq.gz\"`."),
    ])


SECTIONS = {
    "Conversions": "Format conversions with nothing to demultiplex: unmux as a converter between FASTQ, SAM and BAM, compared with samtools.",
    "Illumina demultiplexing": "Short-read runs split by index reads or inline barcodes, with reads of known origin so each result can be scored.",
    "UMIs and duplex sequencing": "Unique molecular identifiers (UMIs) into the RX tag for duplex consensus calling, on reads from the nf-core/fastquorum tests.",
    "Single cell": "Cell barcodes and UMIs into SAM tags, from 10x Chromium, split-pool and antibody-tag libraries.",
    "Long reads": "Nanopore and PacBio reads: barcodes at both ends, U bases, and tags carried in FASTQ comments.",
    "Diagnostics": "Ways to see what a run did and why.",
    "Refusals that protect data": "",
}


def render(ctx, recipes):
    lines = [
        "# unmux cookbook",
        "",
        "Recipes for common unmux jobs. Each one says when to use it, the data it ran on, the exact "
        "command, and an excerpt of what unmux printed or wrote, so you can see the result before "
        "running it on your own data.",
        "",
        callout("NOTE", f"Every command on this page was run with {ctx.unmux_version} and "
                        f"{ctx.samtools_version}, and every output is copied from that run. The script "
                        "[`docs/recipes/generate.py`](recipes/generate.py) fetches the inputs, runs the "
                        "recipes again, and rewrites this page."),
        "",
        "Commands are bash and run in a directory that holds their inputs. Log lines appear "
        "without their timestamps, and an excerpt that skips lines marks the gap with …. Header "
        "excerpts come from `samtools view --no-PG`, which keeps samtools' own `@PG` line out.",
        "",
        "## Contents",
        "",
    ]
    number = 0
    sections = []
    anchors = {}
    for recipe in recipes:
        if not sections or sections[-1][0] != recipe.section:
            sections.append((recipe.section, []))
        number += 1
        sections[-1][1].append((number, recipe))
        anchors[recipe.key] = f"[recipe {number}](#{slug(f'{number}. {recipe.title}')})"
    for section, items in sections:
        lines.append(f"- [{section}](#{slug(section)})")
        for n, recipe in items:
            title = f"{n}. {recipe.title}"
            lines.append(f"  - [{title}](#{slug(title)})")
    lines.append("- [Data sources](#data-sources)")
    lines.append("- [Regenerating this page](#regenerating-this-page)")
    for section, items in sections:
        lines += ["", f"## {section}", ""]
        if SECTIONS.get(section):
            lines += [SECTIONS[section], ""]
        if section == "Refusals that protect data":
            lines += [callout("CAUTION", "Each command in this section is a mistake. unmux refuses it before "
                                         "reading a record or writing a file, exits with status 1, and "
                                         "leaves every file as it was."), ""]
        for n, recipe in items:
            lines += [f"### {n}. {recipe.title}", ""]
            for block in recipe.blocks:
                if block:
                    lines += [block, ""]
    lines += ["## Data sources", "",
              "Every download is pinned: GitHub files by commit, SRA runs by accession and spot range, "
              "NCBI sequences by accession and version.", "",
              "| Data | Pinned URL |", "| --- | --- |"]
    rows = []
    for source in SOURCES:
        row = f"| {source.label} | <{source.url}> |"
        if source in ctx.used and row not in rows:
            rows.append(row)
    lines += rows
    lines += ["", "## Regenerating this page", "",
              "One command fetches the inputs into `docs/recipes/cache/`, builds unmux, runs every "
              "recipe with its checks, and rewrites this page. It needs Python 3, samtools 1.16 or "
              "newer, fg-sra and, for the lima comparison, Docker.", "",
              fence("make -C docs/recipes            # fetch, build, run, and write docs/recipes.md\n"
                    "make -C docs/recipes check      # fail if the page would change\n"
                    "make -C docs/recipes run UNMUX=/path/to/unmux", "bash"),
              "", "The targets are `fetch`, `run`, `page`, `check` and `clean`, and the variables "
              "`SAMTOOLS`, `FG_SRA` and `LIMA` give other ways to run those tools, for example "
              "`SAMTOOLS=\"pixi exec -c conda-forge -c bioconda -s samtools samtools\"`.", ""]
    page = re.sub(r"\n{3,}", "\n\n", "\n".join(lines))
    return re.sub("\x00([a-z_]+)\x00", lambda m: anchors.get(m.group(1), f"recipe {m.group(1)}"), page)


RECIPES = [
    recipe_fastq_to_ubam,
    recipe_ubam_to_fastq,
    recipe_passthrough,
    recipe_split_rg,
    recipe_changeseq,
    recipe_uditas,
    recipe_bins,
    recipe_inline,
    recipe_duplex,
    recipe_duplex_fixed,
    recipe_tenx,
    recipe_splitseq,
    recipe_match,
    recipe_cite,
    recipe_nanopore,
    recipe_pacbio,
    recipe_direct_rna,
    recipe_comment_tags,
    recipe_qc_tag,
    recipe_metrics,
    recipe_collision,
    recipe_wide,
    recipe_refuse_input,
    recipe_refuse_tag,
    recipe_refuse_bin,
]



def run_recipes(ctx, only=None):
    recipes = []
    for recipe in RECIPES:
        if only and recipe.__name__[len("recipe_"):] not in only:
            continue
        key = recipe.__name__[len("recipe_"):]
        print(f"running {key}", file=sys.stderr)
        try:
            result = recipe(ctx)
        except (CheckFailed, AssertionError) as error:
            frame = traceback.extract_tb(error.__traceback__)[-1]
            detail = str(error) or frame.line
            raise SystemExit(f"recipe {key} failed a check at generate.py:{frame.lineno}: {detail}")
        result.key = key
        recipes.append(result)
    return render(ctx, recipes)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("action", choices=["fetch", "run", "page", "check"], nargs="?", default="page",
                        help="fetch the inputs; run every recipe and its checks; also write the page; "
                             "or compare the page with docs/recipes.md")
    parser.add_argument("--unmux", type=Path, default=ROOT / "target" / "release" / "unmux")
    parser.add_argument("--samtools", default=os.environ.get("SAMTOOLS", "samtools"))
    parser.add_argument("--fg-sra", default=os.environ.get("FG_SRA", "fg-sra"))
    parser.add_argument("--lima", default=os.environ.get("LIMA", "lima"))
    parser.add_argument("--cache", type=Path, default=ROOT / "docs" / "recipes" / "cache")
    parser.add_argument("--only", default="", help="with run: a comma-separated list of recipes to run")
    parser.add_argument("--print-pins", action="store_true",
                        help="with fetch: accept unpinned inputs and print a PINS table of what was fetched")
    args = parser.parse_args()
    if args.print_pins:
        global UNPINNED_OK
        UNPINNED_OK = True
    only = [k for k in args.only.split(",") if k]
    if only and args.action != "run":
        raise SystemExit("--only works with run, since a page needs every recipe")
    if not __debug__:
        raise SystemExit("run without -O: the recipe checks are assert statements")
    args.cache.mkdir(parents=True, exist_ok=True)
    unmux = None
    if args.action != "fetch":
        unmux = args.unmux.resolve()
        if not unmux.exists():
            raise SystemExit(f"{unmux} does not exist; build it with `cargo build --release`")
    tools = {"samtools": args.samtools, "fg-sra": args.fg_sra, "lima": args.lima}
    ctx = Context(args.cache.resolve(), unmux, tools)
    for source in SOURCES:
        ctx.fetch(source.key, source.name)
    ctx.used = []
    if args.print_pins:
        print("PINS = {")
        for source in SOURCES:
            print(f'    "{source.key}/{source.name}": "{sha256(args.cache / "downloads" / source.key / source.name)}",')
        print("}")
    if args.action == "fetch":
        print(f"fetched {len(SOURCES)} inputs into {args.cache / 'downloads'}", file=sys.stderr)
        return
    ctx.versions()
    page = run_recipes(ctx, only)
    rendered = ctx.cache / "recipes.md"
    rendered.write_text(page)
    if args.action == "run":
        print(f"every recipe passed its checks; the page is at {rendered}", file=sys.stderr)
    elif args.action == "check":
        current = PAGE.read_text() if PAGE.exists() else ""
        if current != page:
            diff = difflib.unified_diff(current.splitlines(), page.splitlines(), "docs/recipes.md", "rendered", lineterm="")
            print("\n".join(list(diff)[:60]), file=sys.stderr)
            raise SystemExit(f"{PAGE} differs from the rendered page; run `make -C docs/recipes page` to update it")
        print(f"{PAGE} is up to date", file=sys.stderr)
    else:
        PAGE.write_text(page)
        print(f"wrote {PAGE}", file=sys.stderr)


if __name__ == "__main__":
    main()
