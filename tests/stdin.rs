//! Stdin input (`--in 0=-`): the only path that needs a real subprocess with a
//! piped stdin, so it lives here rather than in the in-process `run_demux` unit
//! tests. Covers a plain and a gzipped FASTQ streamed through stdin, and that
//! format is still auto-detected over the non-seekable pipe.

use std::io::Read;

use assert_cmd::Command;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;

/// Read a FASTQ output file back as a single string.
fn read_to_string(path: &std::path::Path) -> String {
    let mut s = String::new();
    std::fs::File::open(path)
        .unwrap()
        .read_to_string(&mut s)
        .unwrap();
    s
}

#[test]
fn stdin_fastq_passes_through_to_out() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.fq");
    Command::cargo_bin("unmux")
        .unwrap()
        .args(["--in", "0=-", "--out", out.to_str().unwrap()])
        .write_stdin("@a\nAAAACCCC\n+\nIIIIIIII\n@b\nTTTTGGGG\n+\nIIIIIIII\n")
        .assert()
        .success();
    let content = read_to_string(&out);
    assert!(
        content.contains("AAAACCCC"),
        "read a passed through: {content}"
    );
    assert!(
        content.contains("TTTTGGGG"),
        "read b passed through: {content}"
    );
}

#[test]
fn bare_invocation_defaults_to_stdin_and_stdout() {
    // No `--in` and no `--out`: file 0 defaults to stdin and output goes to
    // stdout, so a bare `unmux` is a stdin->stdout filter. Format (FASTQ) is
    // sniffed off the pipe and mirrored out.
    let assert = Command::cargo_bin("unmux")
        .unwrap()
        .write_stdin("@a\nAAAACCCC\n+\nIIIIIIII\n@b\nTTTTGGGG\n+\nIIIIIIII\n")
        .assert()
        .success();
    let stdout = String::from_utf8_lossy(&assert.get_output().stdout);
    assert!(
        stdout.contains("AAAACCCC") && stdout.contains("TTTTGGGG"),
        "both reads stream stdin->stdout: {stdout}"
    );
}

#[test]
fn stdin_gzipped_fastq_is_auto_detected() {
    // Format + gzip are sniffed by peeking the pipe (no seek). A gzipped FASTQ
    // on stdin decodes.
    let mut encoder = GzEncoder::new(Vec::new(), Compression::new(6));
    encoder.write_all(b"@a\nAAAACCCC\n+\nIIIIIIII\n").unwrap();
    let gzipped = encoder.finish().unwrap();

    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.fq");
    Command::cargo_bin("unmux")
        .unwrap()
        .args(["--in", "0=-", "--out", out.to_str().unwrap()])
        .write_stdin(gzipped)
        .assert()
        .success();
    assert!(read_to_string(&out).contains("AAAACCCC"));
}

#[test]
fn two_stdin_inputs_are_rejected() {
    Command::cargo_bin("unmux")
        .unwrap()
        .args(["--in", "0=-", "--in", "1=-", "--out", "/dev/null"])
        .write_stdin("@a\nACGT\n+\nIIII\n")
        .assert()
        .failure()
        .stderr(predicates::str::contains(
            "stdin (`-`) may be used for at most one input",
        ));
}

#[test]
fn stdin_short_first_write_is_still_sniffed_as_sam() {
    // Over a pipe a single fill_buf can return only the first (tiny) write, so
    // a producer that emits `@` then (after a flush) the rest of a SAM header
    // must still be sniffed as SAM, not misread as FASTQ. assert_cmd's
    // write_stdin delivers atomically and cannot reproduce this, so drive a raw
    // subprocess with two writes separated by a delay.
    use std::process::{Command as StdCommand, Stdio};
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.sam");
    let mut child = StdCommand::new(env!("CARGO_BIN_EXE_unmux"))
        .args(["--in", "0=-", "--out", out.to_str().unwrap()])
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    {
        let mut stdin = child.stdin.take().unwrap();
        stdin.write_all(b"@").unwrap(); // a lone `@` under-fills a naive single-read sniff
        stdin.flush().unwrap();
        std::thread::sleep(std::time::Duration::from_millis(200));
        stdin
            .write_all(b"HD\tVN:1.6\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n")
            .unwrap();
    } // stdin dropped -> EOF
    let status = child.wait().unwrap();
    assert!(
        status.success(),
        "a short first pipe write must still sniff as SAM"
    );
    assert!(
        read_to_string(&out).contains("r1"),
        "the SAM record round-trips"
    );
}

#[test]
fn stdin_aligned_bam_splits_by_read_group_keeping_alignments() {
    use noodles::sam::alignment::io::Write as _;
    use noodles::{bam, sam};

    let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
                 @RG\tID:rg1\tSM:sA\n@RG\tID:rg2\tSM:sB\n\
                 q1\t99\tchr1\t100\t60\t8M\t=\t150\t58\tACGTACGT\tIIIIIIII\tMC:Z:8M\tRG:Z:rg1\n\
                 q2\t0\tchr1\t120\t30\t8M\t*\t0\t0\tCCCCAAAA\tIIIIIIII\tRG:Z:rg2\n\
                 q1\t147\tchr1\t150\t60\t8M\t=\t100\t-58\tGGGGTTTT\tIIIIIIII\tMC:Z:8M\tRG:Z:rg1\n";
    let mut reader = sam::io::Reader::new(&text[..]);
    let header = reader.read_header().unwrap();
    let records: Vec<_> = reader
        .record_bufs(&header)
        .collect::<Result<_, _>>()
        .unwrap();
    let mut bam_bytes = Vec::new();
    {
        let mut writer = bam::io::Writer::new(&mut bam_bytes);
        writer.write_header(&header).unwrap();
        for record in &records {
            writer.write_alignment_record(&header, record).unwrap();
        }
        writer.try_finish().unwrap();
    }

    let dir = tempfile::tempdir().unwrap();
    let pattern = dir.path().join("%sample.bam");
    Command::cargo_bin("unmux")
        .unwrap()
        .args([
            "-",
            "--group",
            "rg=@RG",
            "--sample-from-group",
            "rg",
            "--out",
            pattern.to_str().unwrap(),
        ])
        .write_stdin(bam_bytes)
        .assert()
        .success();

    let mut reader = bam::io::Reader::new(std::fs::File::open(dir.path().join("rg1.bam")).unwrap());
    let out_header = reader.read_header().unwrap();
    assert_eq!(out_header.reference_sequences().len(), 1);
    assert_eq!(out_header.read_groups().len(), 1);
    let kept: Vec<_> = reader
        .record_bufs(&out_header)
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(kept, vec![records[0].clone(), records[2].clone()]);
}
