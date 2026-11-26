//! A gzip input is sniffed by its first decoded bytes however late in the stream they start, and
//! gzip data that fails to decode fails the run, from a file of any name or from stdin.

use std::io::Write;
use std::path::Path;

use assert_cmd::assert::Assert;
use assert_cmd::Command;
use flate2::{Compression, GzBuilder};
use predicates::str::contains;

/// A FASTQ of `records` reads.
fn fastq(records: usize) -> String {
    (0..records)
        .map(|i| format!("@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n"))
        .collect()
}

/// One gzip member of `content`, its header set by `builder`.
fn gzip_with(builder: GzBuilder, content: &str) -> Vec<u8> {
    let mut encoder = builder.write(Vec::new(), Compression::default());
    encoder.write_all(content.as_bytes()).unwrap();
    encoder.finish().unwrap()
}

/// One gzip member of `content`.
fn gzip(content: &str) -> Vec<u8> {
    gzip_with(GzBuilder::new(), content)
}

/// Run unmux in `dir` on `bytes`, written to the file `input` or piped to stdin when it is `-`.
fn unmux(dir: &Path, input: &str, bytes: &[u8], out: &str) -> Assert {
    let mut command = Command::cargo_bin("unmux").unwrap();
    command.current_dir(dir).args([input, "--out", out]);
    if input == "-" {
        command.write_stdin(bytes);
    } else {
        std::fs::write(dir.join(input), bytes).unwrap();
    }
    command.assert()
}

#[test]
fn a_corrupt_gzip_fails_the_run_whatever_its_name_or_source() {
    let mut corrupt = gzip(&fastq(20_000));
    corrupt[10..14].fill(0xff);
    for input in ["corrupt.dat", "-", "corrupt.fq.gz"] {
        let dir = tempfile::tempdir().unwrap();
        unmux(dir.path(), input, &corrupt, "out.bam")
            .code(1)
            .stderr(contains(format!("failed to open input file 0 ({input})")))
            .stderr(contains("corrupt deflate stream"));
        assert!(!dir.path().join("out.bam").exists(), "{input}");
    }
}

#[test]
fn a_gzip_whose_content_starts_past_the_sniff_window_is_read() {
    let reads = fastq(3);
    let long_header = gzip_with(GzBuilder::new().filename(vec![b'a'; 2_000]), &reads);
    let mut late_member = gzip("").repeat(100);
    late_member.extend(gzip(&reads));
    for (fixture, bytes) in [("long header", long_header), ("late member", late_member)] {
        for input in ["reads.dat", "-", "reads.fq.gz"] {
            let dir = tempfile::tempdir().unwrap();
            unmux(dir.path(), input, &bytes, "out.fq").success();
            let out = std::fs::read_to_string(dir.path().join("out.fq")).unwrap();
            assert_eq!(out, reads, "{fixture} from {input}");
        }
    }
}

#[test]
fn a_gzip_of_empty_content_is_read_as_zero_records() {
    for bytes in [gzip(""), gzip("").repeat(100)] {
        for input in ["empty.dat", "-", "empty.fq.gz"] {
            let dir = tempfile::tempdir().unwrap();
            unmux(dir.path(), input, &bytes, "out.fq").success();
            let out = std::fs::read_to_string(dir.path().join("out.fq")).unwrap();
            assert_eq!(out, "", "{} bytes from {input}", bytes.len());
        }
    }
}
