//! A BAM cut inside a BGZF block or its header fails at every thread count, from a file or from
//! stdin, and a BAM cut at a block boundary or missing its EOF marker block warns.

use std::io::Write;

use assert_cmd::Command;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, bgzf, sam};

/// A BAM whose BGZF blocks each end on a record boundary, cut `into` bytes past the start of a
/// middle block.
fn truncated_bam(into: usize) -> Vec<u8> {
    let (mut bam, block_starts) = bam_with_block_starts();
    bam.truncate(block_starts[block_starts.len() / 2] + into);
    bam
}

/// A BAM of 20,000 records whose BGZF blocks each end on a record boundary, with the offset each
/// block before its EOF marker block starts at.
fn bam_with_block_starts() -> (Vec<u8>, Vec<usize>) {
    let header = sam::Header::default();
    let mut writer = bam::io::Writer::from(Vec::new());
    writer.write_header(&header).unwrap();
    let mut boundaries = vec![writer.get_ref().len()];
    for i in 0..20_000 {
        let record = RecordBuf::builder()
            .set_name(format!("r{i}").into_bytes())
            .set_flags(Flags::UNMAPPED)
            .set_sequence(Sequence::from(b"ACGTACGTAC".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30; 10]))
            .build();
        writer.write_alignment_record(&header, &record).unwrap();
        if i % 1_000 == 999 {
            boundaries.push(writer.get_ref().len());
        }
    }
    let raw = writer.into_inner();

    let mut bgzf = bgzf::io::Writer::new(Vec::new());
    let mut block_starts = Vec::new();
    let mut start = 0;
    for end in boundaries {
        block_starts.push(bgzf.get_ref().len());
        bgzf.write_all(&raw[start..end]).unwrap();
        bgzf.flush().unwrap();
        start = end;
    }
    (bgzf.finish().unwrap(), block_starts)
}

#[test]
fn truncated_bam_fails_at_every_thread_count() {
    let bam = truncated_bam(30);
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.bam");
    std::fs::write(&input, &bam).unwrap();
    let out = dir.path().join("out.bam");
    for threads in ["1", "4", "8", "16"] {
        for stdin in [false, true] {
            let mut command = Command::cargo_bin("unmux").unwrap();
            let source = if stdin { "-" } else { input.to_str().unwrap() };
            command.args([source, "--out", out.to_str().unwrap(), "--threads", threads]);
            if stdin {
                command.write_stdin(bam.clone());
            }
            command
                .assert()
                .failure()
                .stderr(predicates::str::contains("failed to read BAM"));
        }
    }
}

#[test]
fn bam_cut_inside_a_block_header_fails_at_every_thread_count() {
    let bam = truncated_bam(10);
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.bam");
    std::fs::write(&input, &bam).unwrap();
    let out = dir.path().join("out.bam");
    for threads in ["1", "8", "16"] {
        for stdin in [false, true] {
            let mut command = Command::cargo_bin("unmux").unwrap();
            let source = if stdin { "-" } else { input.to_str().unwrap() };
            command.args([source, "--out", out.to_str().unwrap(), "--threads", threads]);
            if stdin {
                command.write_stdin(bam.clone());
            }
            command
                .assert()
                .failure()
                .stderr(predicates::str::contains("is truncated"));
            assert!(!out.exists(), "a failed run commits no output");
        }
    }
}

#[test]
fn bam_cut_at_a_block_boundary_warns_at_every_thread_count() {
    let (intact, block_starts) = bam_with_block_starts();
    let boundary = block_starts[block_starts.len() / 2];
    let eof_block = intact.len() - 28;
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.sam");
    for (name, bam, warns) in [
        ("intact.bam", &intact[..], false),
        ("boundary.bam", &intact[..boundary], true),
        ("no-eof-block.bam", &intact[..eof_block], true),
    ] {
        let input = dir.path().join(name);
        std::fs::write(&input, bam).unwrap();
        for threads in ["1", "8", "16"] {
            for stdin in [false, true] {
                let mut command = Command::cargo_bin("unmux").unwrap();
                let source = if stdin { "-" } else { input.to_str().unwrap() };
                command.args([source, "--out", out.to_str().unwrap(), "--threads", threads]);
                if stdin {
                    command.write_stdin(bam.to_vec());
                }
                let output = command.assert().success().get_output().clone();
                let stderr = String::from_utf8_lossy(&output.stderr);
                let context = format!("{name} at {threads} threads, stdin {stdin}: {stderr}");
                let warning = stderr
                    .lines()
                    .find(|line| line.contains("BGZF EOF marker is absent"));
                assert_eq!(warning.is_some(), warns, "{context}");
                assert!(
                    warning.is_none_or(|line| line.contains(source)),
                    "{context}"
                );
            }
        }
    }
}
