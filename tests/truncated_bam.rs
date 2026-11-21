//! A BAM cut inside a BGZF block fails at every thread count, from a file or from stdin.

use std::io::Write;

use assert_cmd::Command;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;
use noodles::{bam, bgzf, sam};

/// A BAM whose BGZF blocks each end on a record boundary, cut partway into a middle block.
fn truncated_bam() -> Vec<u8> {
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
    let mut compressed = bgzf.finish().unwrap();
    compressed.truncate(block_starts[block_starts.len() / 2] + 30);
    compressed
}

#[test]
fn truncated_bam_fails_at_every_thread_count() {
    let bam = truncated_bam();
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
