//! The decoded outputs, metrics, and errors of representative runs, pinned at every thread count.

use std::fs;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::OnceLock;

use assert_cmd::Command;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use md5::{Digest, Md5};
use noodles::sam::alignment::io::Write as _;
use noodles::{bam, bgzf, cram, sam};
use rstest::rstest;
use tempfile::TempDir;

const THREADS: [&str; 6] = ["1", "2", "3", "4", "8", "16"];

const BARCODES: [&str; 12] = [
    "AACCGGTT", "ACGTACGT", "AGTCAGTC", "ATATCGCG", "CAGTTGAC", "CCAATTGG", "CGATATCG", "CTGACTGA",
    "GACTGTCA", "GGTTAACC", "TACGCATG", "TTGGCCAA",
];

const LINKER: &str = "GTGGCCGATGTTTCGC";

/// A SplitMix64 generator.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn one_in(&mut self, n: usize) -> bool {
        self.below(n) == 0
    }

    fn bases(&mut self, n: usize) -> String {
        (0..n).map(|_| char::from(b"ACGT"[self.below(4)])).collect()
    }

    fn quals(&mut self, n: usize) -> String {
        (0..n)
            .map(|_| char::from(b'#' + self.below(40) as u8))
            .collect()
    }

    /// A barcode, or rarely random bases, with rare substitutions, insertions, and deletions.
    fn barcode(&mut self) -> String {
        let mut bases = match self.one_in(20) {
            true => self.bases(8),
            false => BARCODES[self.below(12)].to_string(),
        };
        if self.one_in(6) {
            let at = self.below(8);
            let base = ["A", "C", "G", "T", "N", "a", "R"][self.below(7)];
            bases.replace_range(at..at + 1, base);
        }
        if self.one_in(15) {
            bases.remove(self.below(8));
        }
        if self.one_in(15) {
            let at = self.below(7);
            bases.insert_str(at, &self.bases(1));
        }
        bases
    }
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::new(6));
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

fn bgzip(bytes: &[u8], block: usize) -> Vec<u8> {
    let mut writer = bgzf::io::Writer::new(Vec::new());
    for chunk in bytes.chunks(block) {
        writer.write_all(chunk).unwrap();
        writer.flush().unwrap();
    }
    writer.finish().unwrap()
}

/// One FASTQ record, with `comment` holding its leading separator.
fn fastq(name: &str, comment: &str, bases: &str, quals: &str) -> String {
    format!("@{name}{comment}\n{bases}\n+\n{quals}\n")
}

/// Encode SAM text as the format that `name` ends with.
fn write_alignments(dir: &Path, name: &str, text: &str) {
    let mut reader = sam::io::Reader::new(text.as_bytes());
    let header = reader.read_header().unwrap();
    let records: Vec<_> = reader
        .record_bufs(&header)
        .collect::<Result<_, _>>()
        .unwrap();
    let file = fs::File::create(dir.join(name)).unwrap();
    if name.ends_with(".bam") {
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        for record in &records {
            writer.write_alignment_record(&header, record).unwrap();
        }
        writer.try_finish().unwrap();
    } else {
        let mut writer = cram::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        for record in &records {
            writer.write_alignment_record(&header, record).unwrap();
        }
        writer.try_finish(&header).unwrap();
    }
}

/// Every input file, written once per test process.
fn inputs() -> &'static Path {
    static DIR: OnceLock<TempDir> = OnceLock::new();
    DIR.get_or_init(|| {
        let dir = TempDir::new().unwrap();
        write_inputs(dir.path());
        dir
    })
    .path()
}

fn write_inputs(dir: &Path) {
    let put = |name: &str, bytes: &[u8]| fs::write(dir.join(name), bytes).unwrap();
    let mut rng = Rng(7);

    let rows = |columns: &str, cell: &dyn Fn(usize) -> String| -> String {
        (0..12).fold(format!("{columns}\n"), |text, i| {
            format!("{text}{}\n", cell(i))
        })
    };
    let tags = |i: usize| format!("t{i:02}\t{}\tcell {}", BARCODES[i], i % 5);
    put("bc.tsv", rows("id\tseq\tcell", &tags).as_bytes());
    let sheet = |i: usize| format!("s{}\tlib{}\tt{i:02}\tn{i}", i % 7, i % 7 % 3);
    put(
        "sheet.tsv",
        rows("sample\tsub_sample\tbc\tnote", &sheet).as_bytes(),
    );

    let [mut r1, mut r2, mut i1, mut i2, mut third] = [const { String::new() }; 5];
    let [mut inter, mut any, mut fasta, mut sc1, mut sc2] = [const { String::new() }; 5];
    let mut flat = String::from("@HD\tVN:1.6\n@RG\tID:a\tSM:one\tLB:x\n@RG\tID:b\tSM:two\tLB:y\n");
    let mut pairs = flat.clone();
    let mut mates = String::from("@HD\tVN:1.6\n");
    let mut aligned = String::from(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n\
         @RG\tID:a\tSM:one\tLB:x\tPL:ILLUMINA\n@RG\tID:b\tSM:two\tLB:y\n",
    );
    for i in 0..3000 {
        let name = format!("read:{i}");
        let comment = match i % 6 {
            0 => " 1:N:0:ACGT",
            1 => "\tXX:Z:two words\tXI:i:-5",
            2 => " free text\tNM:i:3\tXF:f:1.5\tXA:A:q",
            _ => "",
        };
        let long = rng.below(60);
        let umi = rng.bases(6);
        let one = format!("{}{umi}{}", rng.barcode(), rng.bases(long));
        let quals = match i % 500 {
            3 => " ~!".repeat(one.len())[..one.len()].to_string(),
            _ => rng.quals(one.len()),
        };
        r1 += &fastq(&name, comment, &one, &quals);
        let mut two = if i % 251 == 9 {
            String::new()
        } else {
            let long = 30 + rng.below(40);
            rng.bases(long)
        };
        if rng.one_in(12) && two.len() > 20 {
            two.replace_range(4..16, "GATCGGAAGAGC");
        }
        let two_quals = rng.quals(two.len());
        r2 += &fastq(&name, "", &two, &two_quals);
        inter += &fastq(&format!("{name}/1"), comment, &one, &quals);
        inter += &fastq(&format!("{name}/2"), "", &two, &two_quals);
        let (a, b) = (rng.barcode(), rng.barcode());
        i1 += &fastq(&name, "", &a, &rng.quals(a.len()));
        i2 += &fastq(&name, "", &b, &rng.quals(b.len()));

        let stagger = rng.below(4);
        let cell = format!(
            "{}{}{LINKER}{}{}{}",
            rng.bases(stagger),
            rng.barcode(),
            rng.barcode(),
            rng.bases(6),
            rng.bases(20)
        );
        let cell = match rng.one_in(30) {
            true => cell.replace("CCGATG", "CCATG"),
            false => cell,
        };
        sc1 += &fastq(&name, "", &cell, &rng.quals(cell.len()));
        sc2 += &fastq(&name, "", &two, &two_quals);
        third += &fastq(&name, "", "ACGTACGTAC", "IIIIIIIIII");

        let (left, right) = (rng.below(40), rng.below(40));
        let mut anywhere = format!("{}{}{}", rng.bases(left), rng.barcode(), rng.bases(right));
        if rng.one_in(8) {
            anywhere += &rng.barcode();
        }
        if rng.one_in(9) {
            anywhere = anywhere[4.min(anywhere.len())..].to_string();
        }
        if rng.one_in(3) {
            let complement = |base: char| match base {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                other => other,
            };
            anywhere = anywhere.chars().rev().map(complement).collect();
        }
        any += &fastq(&name, "", &anywhere, &rng.quals(anywhere.len()));
        fasta += &format!(">{name}{comment}\n{}\n{}\n", &one[..10], &one[10..]);

        let letters: String = one.to_uppercase().replace('R', "N");
        let letter_quals = quals.replace(' ', "!");
        let flag = if i % 5 == 1 { 16 } else { 4 };
        let group = ["\tRG:Z:a", "\tRG:Z:b", ""][i % 3];
        flat += &format!(
            "{name}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{letters}\t{letter_quals}\t\
             CB:Z:old\tCY:Z:III\tXB:B:s,1,-2,{}\tAS:i:{}{group}\n",
            i % 300,
            i % 90
        );
        mates += &format!("{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{letters}\t{letter_quals}\n");
        let (second, second_quals) = match two.is_empty() {
            true => ("ACGTN", "IIII#"),
            false => (two.as_str(), two_quals.as_str()),
        };
        let mate =
            format!("{name}\t141\t*\t0\t0\t*\t*\t0\t0\t{second}\t{second_quals}\tXM:i:2{group}\n");
        let first =
            format!("{name}\t77\t*\t0\t0\t*\t*\t0\t0\t{letters}\t{letter_quals}\tXM:i:1{group}\n");
        pairs += &if i % 4 == 0 {
            mate + &first
        } else {
            first + &mate
        };
        let flag = [0, 16, 256, 2048, 4, 1024][i % 6];
        aligned += &format!(
            "{name}\t{flag}\tchr1\t{}\t{}\t{}M\t*\t0\t0\t{letters}\t{letter_quals}\t\
             NM:i:{}\tXQ:Z:q{group}\n",
            1 + i * 7,
            i % 61,
            letters.len(),
            i % 4
        );
    }
    put("r1.fq", r1.as_bytes());
    put("r2.fq", r2.as_bytes());
    put("i1.fq", i1.as_bytes());
    put("i2.fq", i2.as_bytes());
    put("third.fq", third.as_bytes());
    put("sc1.fq", sc1.as_bytes());
    put("sc2.fq", sc2.as_bytes());
    put("any.fq", any.as_bytes());
    put("inter.fq", inter.as_bytes());
    put(
        "inter.odd.fq",
        (inter.clone() + &fastq("z/1", "", "ACGT", "IIII")).as_bytes(),
    );
    put("reads.fa", fasta.as_bytes());
    put("flat.sam", flat.as_bytes());
    put("crlf.r1.fq", r1.replace('\n', "\r\n").as_bytes());
    put("noeol.r2.fq", r2.trim_end_matches('\n').as_bytes());
    put("r1.fq.gz", &gzip(r1.as_bytes()));
    put("r2.bgzf.fq.gz", &bgzip(r2.as_bytes(), 4099));
    let (head, tail) = r1.split_at(r1.len() / 3);
    put(
        "r1.members.fq.gz",
        &[gzip(head.as_bytes()), gzip(b""), gzip(tail.as_bytes())].concat(),
    );
    put(
        "r1.mixed.fq.gz",
        &[bgzip(head.as_bytes(), 60_000), gzip(tail.as_bytes())].concat(),
    );
    write_alignments(dir, "flat.bam", &flat);
    write_alignments(dir, "flat.cram", &flat.replace("\t16\t*\t", "\t4\t*\t"));
    write_alignments(dir, "pairs.bam", &pairs);
    write_alignments(dir, "mates.bam", &mates);
    write_alignments(dir, "aligned.bam", &aligned);

    let mut big = String::new();
    for i in 0..70_000 {
        let long = 4 + rng.below(12);
        let bases = format!("{}{}", rng.barcode(), rng.bases(long));
        big += &fastq(&format!("b{i}"), "", &bases, &rng.quals(bases.len()));
    }
    put("big.fq", big.as_bytes());
    put("big.fq.gz", &gzip(big.as_bytes()));
    put("big.bgzf.fq.gz", &bgzip(big.as_bytes(), 65_000));

    let record_at = |n: usize| big.match_indices("\n@b").nth(n - 1).unwrap().0 + 1;
    let (middle, next) = (record_at(40_000), record_at(40_001));
    let record = &big[middle..next];
    let broken = |name: &str, text: &str| {
        put(
            name,
            [&big[..middle], text, &big[next..]].concat().as_bytes(),
        );
    };
    broken("bad-prefix.fq", &record[1..]);
    broken("bad-plus.fq", &record.replace("\n+\n", "\nX\n"));
    broken("bad-length.fq", &record.replace("\n+\n", "A\n+\n"));
    broken("bad-tag.fq", &record.replacen('\n', "\tML:B:C,256\n", 1));
    broken("blank-line.fq", &format!("\n{record}"));
    broken("renamed.fq", &record.replace("@b", "@c"));
    put("no-quals.fq", &big.as_bytes()[..middle - 30]);
    put("half.fq", &big.as_bytes()[..middle]);
    let whole = gzip(big.as_bytes());
    put("cut.fq.gz", &whole[..whole.len() / 2]);
    let mut flipped = whole;
    let checksum = flipped.len() - 8;
    flipped[checksum] ^= 0xff;
    put("bad-checksum.fq.gz", &flipped);
    let blocks = bgzip(big.as_bytes(), 65_000);
    put("cut.bgzf.fq.gz", &blocks[..blocks.len() / 2]);
    let mut flipped = blocks;
    let mut end = 0;
    while end < flipped.len() / 2 {
        end += usize::from(u16::from_le_bytes([flipped[end + 16], flipped[end + 17]])) + 1;
    }
    flipped[end - 8] ^= 0xff;
    put("bad-checksum.bgzf.fq.gz", &flipped);
}

/// Text without its `@PG` lines, which hold the command line.
fn without_programs(text: &[u8]) -> Vec<u8> {
    text.split_inclusive(|&b| b == b'\n')
        .filter(|line| !line.starts_with(b"@PG"))
        .flatten()
        .copied()
        .collect()
}

/// The decoded content of an output, by its name or else its first bytes.
fn decoded(name: &str, bytes: Vec<u8>) -> Vec<u8> {
    if name.ends_with(".cram") {
        let mut reader = cram::io::Reader::new(&bytes[..]);
        let header = reader.read_header().unwrap();
        let mut writer = sam::io::Writer::new(Vec::new());
        writer.write_header(&header).unwrap();
        for record in reader.records(&header) {
            writer
                .write_alignment_record(&header, &record.unwrap())
                .unwrap();
        }
        return without_programs(writer.get_ref());
    }
    let mut content = bytes;
    if content.starts_with(&[0x1f, 0x8b]) {
        let mut inflated = Vec::new();
        MultiGzDecoder::new(&content[..])
            .read_to_end(&mut inflated)
            .unwrap();
        content = inflated;
    }
    if !content.starts_with(b"BAM\x01") {
        return without_programs(&content);
    }
    let text = u32::from_le_bytes(content[4..8].try_into().unwrap()) as usize;
    let (header, records) = content[8..].split_at(text);
    [without_programs(header), records.to_vec()].concat()
}

/// Every file under `dir`, as sorted relative paths.
fn tree(dir: &Path) -> Vec<String> {
    let mut paths = Vec::new();
    let mut pending = vec![dir.to_path_buf()];
    while let Some(next) = pending.pop() {
        for entry in fs::read_dir(next).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                pending.push(path);
            } else {
                let relative = path.strip_prefix(dir).unwrap();
                paths.push(relative.to_string_lossy().into_owned());
            }
        }
    }
    paths.sort();
    paths
}

/// Run `unmux` in a fresh directory, with `{in}` standing for the directory of the inputs.
fn unmux(args: &str, threads: &str) -> (TempDir, std::process::Output) {
    let dir = TempDir::new().unwrap();
    let args = args.replace("{in}", &inputs().to_string_lossy());
    let output = Command::new(env!("CARGO_BIN_EXE_unmux"))
        .current_dir(dir.path())
        .env("RUST_LOG", "error")
        .env("NO_COLOR", "1")
        .args(args.split_whitespace())
        .args(["--threads", threads])
        .output()
        .unwrap();
    (dir, output)
}

/// A digest of every output of a run, stdout included.
fn digest(args: &str, threads: &str) -> String {
    let (dir, output) = unmux(args, threads);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "--threads {threads}: {stderr}");
    let mut hasher = Md5::new();
    let mut add = |name: &str, bytes: Vec<u8>| {
        let content = decoded(name, bytes);
        hasher.update(format!("{name} {}\n", content.len()));
        hasher.update(&content);
    };
    for name in tree(dir.path()) {
        add(&name, fs::read(dir.path().join(&name)).unwrap());
    }
    add("stdout", output.stdout);
    format!("{:x}", hasher.finalize())
}

const BC: &str = "--group bc={in}/bc.tsv --group bc::loc=0:0:8 --group bc::dist=1";
const METRICS: &str =
    "--metrics-per-sample ps.tsv --metrics-summary sum.tsv --metrics-unassigned un.tsv";
const CUT: &str = "--extract a=0:14:end --extract b=1:0:end --extract umi=0:8:14 --extract s=@bc";
const BINS: &str = "--unassigned un.R%source.fq.gz";
const CELLS: &str = "--group b1={in}/bc.tsv --group b1::loc=0:0:13 --group b1::dist=1:1:1 \
     --group b1::mode=all \
     --group l={GTGGCCGATGTTTCGC} --group l::loc=0:6:32 --group l::dist=2:1:3 \
     --group b2={in}/bc.tsv --group b2::dist=1:1:1 --group b2::mode=all \
     --group l::next=b2:0-9 \
     --extract c1=@b1 --extract c2=@b2 --extract umi=@b2+0:6 --extract cdna=1:0:end \
     --extract gap=@b1..@l --extract name=@b1::id";

#[rstest]
#[case::paired_fastq(
    "6e2d7a8606af3af998aea22a8be04d97",
    "{in}/r1.fq {in}/r2.fq {BC} {CUT} \
    --template a b --tag RX=umi --tag BC=s --sample-from-group bc --out %sample.R%ordinal.fq.gz \
    {BINS} {METRICS}"
)]
#[case::gzip_to_fastq(
    "d78d615bf817217085361d8ebcf42a67",
    "{in}/r1.fq.gz {in}/r2.bgzf.fq.gz {BC} \
    {CUT} --template a b --tag RX=umi --sample-from-group bc --out %sample.R%ordinal.fq \
    --unassigned un.fq {METRICS}"
)]
#[case::gzip_members_to_bam(
    "88a7fe6ab0d6bdd0e319e5c3cf5258e4",
    "{in}/r1.members.fq.gz {in}/r2.fq \
    {BC} {CUT} --template a b --tag BC=s --tag RX=umi --sample-from-group bc --out %sample.bam \
    {BINS} {METRICS}"
)]
#[case::mixed_gzip_to_one_bam(
    "a4a8c7d879c73548ae18236d18c1e746",
    "{in}/r1.mixed.fq.gz \
    {in}/r2.bgzf.fq.gz {BC} {CUT} --template ~b a+s --tag BC=~s --tag OX=umi+s \
    --sample-from-group bc --out all.bam --unassigned un.bam"
)]
#[case::crlf_and_no_final_newline(
    "46153248a8ad077d766c6ed6f9a8666f",
    "{in}/crlf.r1.fq \
    {in}/noeol.r2.fq {BC} {CUT} --template a b --tag BC=s --tag BC::raw=true \
    --sample-from-group bc --out %sample.R%ordinal.fa.gz --unassigned un.R%source.fa {METRICS}"
)]
#[case::interleaved_to_sam(
    "8e4faf8141addeae099c1d7ecd5d9275",
    "{in}/inter.fq {BC} \
    --extract a=0:8:end --extract b=1:0:end --template a b --sample-from-group bc \
    --out %sample.sam --unassigned un.sam --qc-tag"
)]
#[case::interleaved_per_record(
    "d2e828b787627914ea23bf6d6146f461",
    "{in}/inter.fq --per-record \
    {BC} --extract a=0:8:end --template a --sample-from-group bc --out %sample.fq.gz {BINS} \
    {METRICS}"
)]
#[case::index_reads_matched_joined(
    "adb5e05cf4cbe48e7cfabde19ad6ee39",
    "{in}/r1.fq {in}/r2.fq \
    {in}/i1.fq {in}/i2.fq --group bc={in}/bc.tsv --group bc::match=x+y --group bc::dist=1 \
    --group bc::mode=nearest --extract x=2:0:4 --extract y=2:4:8 --extract a=0:0:end \
    --extract b=1:0:end --template a b --tag BC=x+y --sample-sheet {in}/sheet.tsv \
    --sample-sheet-ignore note --out %sample.%sub_sample.R%ordinal.fq.gz {BINS} {METRICS}"
)]
#[case::dual_index_sheet_to_cram(
    "418bc5c79d4aa89f2bf07347bada68a8",
    "{in}/r1.fq {in}/r2.fq \
    {in}/i1.fq {in}/i2.fq --group bc={in}/bc.tsv --group bc::loc=2:0:8 --group bc::dist=1 \
    --group i5={in}/bc.tsv --group i5::loc=3:0:9 --group i5::dist=1:1:1 \
    --group i5::mode=all \
    --extract a=0:0:end --extract b=1:0:end --extract x=@bc --extract y=@i5 \
    --extract n=%sample::note --extract c=@i5::cell --template a b --tag BC=x+y --tag xn=n \
    --tag xc=c+x --sample-sheet {in}/sheet.tsv --out %sample.cram --rg-tag PL=ILLUMINA {METRICS}"
)]
#[case::cells_to_bam(
    "d1c03d2ec28be5a03f34e1f0371bfdcd",
    "{in}/sc1.fq {in}/sc2.fq {CELLS} \
    --group b1::minFindsPerGroup=1 --template cdna --tag CB=c1+c2 --tag CR=c1+c2 \
    --tag CR::raw=true --tag RX=umi --tag xg=gap --tag xn=name --out out.bam --unassigned un.bam \
    --qc-tag=ZQ {METRICS}"
)]
#[case::cells_split_with_removal(
    "198650b175068301447cd77d9feea281",
    "{in}/sc1.fq {in}/sc2.fq \
    {CELLS} --group ad={GATCGGAAGAGC} --group ad::loc=1:0:end --group ad::dist=1 \
    --remove ad=ad.%source.fq.gz --remove b1::t03 --template cdna+umi --tag CB=c1+c2 --tag RX=umi \
    --sample-from-group b1 --out %sample.fq.gz --unassigned un.%source.fq --qc-tag {METRICS}"
)]
#[case::anywhere_on_both_strands(
    "38e250fe91c23bacfa95b118126a3be7",
    "{in}/any.fq \
    --group bc={in}/bc.tsv --group bc::dist=1 --group bc::bothStrands=true \
    --group bc::mode=nearest --extract s=@bc --extract r=0:0:end --extract l=@bc-0:9 --template r \
    --tag CB=s --tag xl=l --sample-from-group bc --out %sample.fq.gz --unassigned un.fq.gz \
    {METRICS}"
)]
#[case::anywhere_with_indels(
    "b935414b7cc2b0a4b1a01471b39e2623",
    "{in}/any.fq \
    --group bc={in}/bc.tsv --group bc::dist=1:1:2 --group bc::mode=nearest --group bc::delta=1 \
    --extract s=@bc --extract r=@bc+0:end --template r --tag CB=s --tag CR=s --tag CR::raw=true \
    --sample-from-group bc --out %sample.bam --unassigned un.bam --qc-tag {METRICS}"
)]
#[case::anywhere_on_both_strands_with_indels(
    "dde8d56f319c14c10d25dc1b3cdbdbbd",
    "{in}/any.fq \
    --group bc={in}/bc.tsv --group bc::dist=1:1:1 --group bc::bothStrands=true \
    --group bc::mode=nearest \
    --group bc::maxFindsPerTag=1 --extract s=@bc --extract r=0:0:end --template r --tag CB=s \
    --sample-from-group bc --out %sample.bam --unassigned un.bam --qc-tag {METRICS}"
)]
#[case::anywhere_on_both_strands_with_indels_apart(
    "b33e26e12be76b5f7a3fe6420ea9d3cf",
    "{in}/any.fq \
    --group bc={a=AACCGGTT,b=ACGTACGT,c=AGTCAGTC,d=ATATCGCG,e=CAGTTGAC,f=CCAATTGG,g=GGTTAACC,h=TACGCATG,i=TTGGCCAA} \
    --group bc::dist=0:1 --group bc::bothStrands=true --group bc::maxFindsPerTag=1 \
    --extract s=@bc --extract r=0:0:end --template r --tag CB=s \
    --sample-from-group bc --out %sample.bam --unassigned un.bam --qc-tag {METRICS}"
)]
#[case::anywhere_exact_and_unique(
    "14e0f61dc6b46800628b3fe995adbda7",
    "{in}/any.fq \
    --group bc={in}/bc.tsv --group bc::findOne --extract r=0:0:end --template r \
    --sample-from-group bc --out %sample.fq --unassigned un.fq {METRICS}"
)]
#[case::few_tags_anywhere(
    "37a715883cb0433ff7b76f2fb0edfeae",
    "{in}/any.fq \
    --group x={a=AACCGGTT,b=ACGTNCGT,c=TTGGCCAAG} --group x::dist=2 --group x::bothStrands=true \
    --group x::mode=nearest --extract r=0:0:end --extract s=@x --template r --tag CB=s \
    --sample-from-group x --out %sample.sam --unassigned un.sam --qc-tag {METRICS}"
)]
#[case::partial_ends(
    "1b1297829ced425558ce06aaf44164ea",
    "{in}/any.fq --group bc={in}/bc.tsv \
    --group bc::dist=1 --group bc::partial5=5 --group bc::partial3=5:0.2 \
    --group bc::bothStrands=true --group bc::mode=nearest --extract s=@bc --extract r=0:0:end \
    --template r --tag CB=s --sample-from-group bc --out %sample.bam --unassigned un.bam --qc-tag \
    {METRICS}"
)]
#[case::anchored(
    "19ed39d8603b82ed91b916abdcce7cfb",
    "{in}/r1.fq --group bc={in}/bc.tsv \
    --group bc::loc=0:0:9 --group bc::dist=1:1:1 --group bc::anchor=5p --group bc::mode=nearest \
    --group e={x=A,c=C,g=G} --group e::loc=0:5:-1 --group e::anchor=3p --extract s=@bc \
    --extract r=@bc+0:end --extract t=@e --template r --tag CB=s --tag xt=t --sample-from-group bc \
    --out %sample.fq.gz --unassigned un.fq.gz --qc-tag {METRICS}"
)]
#[case::missing_template_stream(
    "b750bd750071f95f9f73217706e24663",
    "{in}/r1.fq {in}/r2.fq \
    --group bc={in}/bc.tsv --group bc::loc=0:0:8 --extract s=@bc --extract b=1:0:end \
    --template s b --tag RX=s --out out.R%ordinal.fq.gz --unassigned un.%source.fq.gz --qc-tag \
    {METRICS}"
)]
#[case::fasta(
    "9874ee22a96a0a291a7c739cb18bd242",
    "{in}/reads.fa {BC} --extract a=0:8:end \
    --extract s=@bc --template a --tag BC=s --sample-from-group bc --out %sample.fa.gz \
    --unassigned un.fa {METRICS}"
)]
#[case::sam_to_bam(
    "34e48fefa48da8dedce930f5dc273d52",
    "{in}/flat.sam {BC} --extract a=0:8:end \
    --extract s=@bc --template a --tag CB=s --sample-from-group bc --out %sample.bam \
    --unassigned un.bam {METRICS}"
)]
#[case::bam_tagged_whole(
    "001fba68e1eea1954fc2cf7815dde6fc",
    "{in}/flat.bam {BC} --extract s=@bc \
    --tag CB=s --tag CB::qual=none --tag xc=c --extract c=@bc::cell --out out.bam --qc-tag \
    {METRICS}"
)]
#[case::paired_bam_to_bam(
    "1cf50f172fcf437f49f43edb1f82b405",
    "{in}/pairs.bam {BC} \
    --extract a=0:8:end --extract b=1:0:end --extract s=@bc --template a b --tag BC=s \
    --sample-from-group bc --out %sample.bam --unassigned un.bam {METRICS}"
)]
#[case::paired_bam_to_fastq(
    "e39593f936c9b2f2b02c8cb788116683",
    "{in}/pairs.bam \
    --out out.R%ordinal.fq.gz"
)]
#[case::cram_to_bam(
    "9cacebfa1c105e8cb63cd034ca9ca3cc",
    "{in}/flat.cram {BC} --extract a=0:8:end \
    --extract s=@bc --template a --tag CB=s --sample-from-group bc --out %sample.bam \
    --unassigned un.fq.gz {METRICS}"
)]
#[case::bam_split_by_read_group(
    "c58f0daafb88466f600065616b3e377f",
    "{in}/flat.bam \
    --group rg=@RG::SM::LB --sample-from-group rg --out %sample.%sub_sample.bam \
    --unassigned un.sam {METRICS}"
)]
#[case::aligned_split(
    "d5296696cef46049cfc76fab626dc3a7",
    "{in}/aligned.bam --group rg=@RG \
    --sample-from-group rg --out %sample.bam --unassigned un.bam {METRICS}"
)]
#[case::aligned_split_to_sam(
    "97d887fe7aa4a1774c7c1f221868c268",
    "{in}/aligned.bam \
    --group rg=@RG --sample-from-group rg --out %sample.sam --unassigned un.bam --qc-tag \
    {METRICS}"
)]
#[case::passed_through_to_bam(
    "22ff70356631c3da24c66e929150de1b",
    "{in}/r1.fq {in}/r2.fq \
    --out out.bam --metrics-summary sum.tsv"
)]
#[case::passed_through_to_cram(
    "08116100f0ee24d34c26c00ec7d8dc4e",
    "{in}/r1.fq {in}/r2.fq \
    --out out.cram"
)]
#[case::passed_through_to_stdout("c3bce5e725ff04c1b33d6aa1571e6115", "{in}/r1.fq.gz")]
#[case::blend_to_stdout("8ddb9e1c9c7773fc0b40693e6faf4503", "{in}/r1.fq {in}/mates.bam")]
#[case::stored_blocks(
    "62786fa6478896ae3f54170584685ef1",
    "{in}/r1.fq {BC} \
    --sample-from-group bc --out %sample.fq.gz --unassigned un.bam --compression 0"
)]
#[case::best_compression(
    "3bfe5a72f868c41fd11665555984acfc",
    "{in}/r1.fq {BC} \
    --sample-from-group bc --out %sample.bam --unassigned un.fq.gz --compression 9"
)]
#[case::many_batches_to_fastq(
    "0371f4d9647519802ab0155117348cd2",
    "{in}/big.fq.gz {BC} \
    --extract a=0:8:end --extract s=0:0:8 --template a --tag BC=s --sample-from-group bc \
    --out %sample.fq.gz {BINS} {METRICS}"
)]
#[case::many_batches_to_bam(
    "dd6b21d311d07f4bd028344fb98b8f76",
    "{in}/big.bgzf.fq.gz {BC} \
    --extract s=@bc --tag BC=s --sample-from-group bc --out all.bam --unassigned un.fq.gz \
    {METRICS}"
)]
fn outputs_are_pinned_at_every_thread_count(#[case] expected: &str, #[case] args: &str) {
    let args = args
        .replace("{BC}", BC)
        .replace("{METRICS}", METRICS)
        .replace("{CUT}", CUT)
        .replace("{BINS}", BINS)
        .replace("{CELLS}", CELLS);
    for threads in THREADS {
        assert_eq!(digest(&args, threads), expected, "--threads {threads}");
    }
}

#[rstest]
#[case::name_prefix(
    "failed to read FASTQ record: invalid name prefix",
    "{in}/bad-prefix.fq"
)]
#[case::plus_line(
    "failed to read FASTQ record: invalid description prefix",
    "{in}/bad-plus.fq"
)]
#[case::lengths("malformed FASTQ record `b40000`", "{in}/bad-length.fq")]
#[case::comment_tag(
    "malformed FASTQ record `b40000`: the comment tag `ML:B`",
    "{in}/bad-tag.fq"
)]
#[case::blank_line(
    "failed to read FASTQ record: invalid name prefix",
    "{in}/blank-line.fq"
)]
#[case::truncated_text("FASTQ record", "{in}/no-quals.fq")]
#[case::record_counts("unequal record counts", "{in}/big.fq {in}/half.fq")]
#[case::record_counts_first_shorter("unequal record counts", "{in}/half.fq {in}/big.fq")]
#[case::truncated_second_input("FASTQ record", "{in}/big.fq {in}/no-quals.fq")]
#[case::interleaved_odd("an unpaired mate (odd number of records)", "{in}/inter.odd.fq")]
#[case::record_names(
    "mismatched read names (`b40000` vs `c40000`)",
    "{in}/big.fq {in}/renamed.fq"
)]
#[case::truncated_gzip("failed to read FASTQ record", "{in}/cut.fq.gz")]
#[case::gzip_checksum("failed to read FASTQ record", "{in}/bad-checksum.fq.gz")]
#[case::truncated_bgzf("failed to read FASTQ record", "{in}/cut.bgzf.fq.gz")]
#[case::bgzf_checksum("failed to read FASTQ record", "{in}/bad-checksum.bgzf.fq.gz")]
#[case::no_qualities("no qualities", "{in}/reads.fa --template a --extract a=0:0:4")]
#[case::too_many_mates(
    "at most two reads per record (a pair); got 3",
    "{in}/r1.fq {in}/r2.fq {in}/third.fq --unassigned un.bam"
)]
fn errors_leave_no_output_at_every_thread_count(#[case] expected: &str, #[case] args: &str) {
    for threads in ["1", "2", "8"] {
        let (dir, output) = unmux(
            &format!(
                "{args} --group x={{GGGGGGGG}} --group x::loc=0:0:8 --group x::dist=8 \
                 --group x::minFindsPerGroup={finds} --out out.R%ordinal.fq.gz",
                finds = if args.contains("un.bam") { 2 } else { 1 }
            ),
            threads,
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        let code = output.status.code();
        assert_eq!(code, Some(1), "--threads {threads}: {stderr}");
        assert!(stderr.contains(expected), "--threads {threads}: {stderr}");
        assert_eq!(tree(dir.path()), Vec::<String>::new());
    }
}
