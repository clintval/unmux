//! A gzipped FASTQ, FASTA, or SAM cut short fails at every thread count, from a file or from
//! stdin, even when the cut leaves only whole records. Intact and concatenated gzip members read
//! in full, and a bgzipped input cut at a block boundary or missing its EOF marker block warns.

use std::io::Write;
use std::ops::Range;
use std::path::Path;

use assert_cmd::assert::Assert;
use assert_cmd::Command;
use flate2::write::GzEncoder;
use flate2::Compression;
use noodles::bgzf;

/// The text formats, by file extension.
const FORMATS: [&str; 3] = ["fq", "fa", "sam"];

/// The records `range` of `format`, as text.
fn records(format: &str, range: Range<usize>) -> Vec<u8> {
    range
        .map(|i| match format {
            "fq" => format!("@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n"),
            "fa" => format!(">r{i}\nACGTACGTAC\n"),
            _ => format!("r{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"),
        })
        .collect::<String>()
        .into_bytes()
}

/// The records `range` of `format`, as a whole file: a SAM file starts with its header.
fn file(format: &str, range: Range<usize>) -> Vec<u8> {
    let header: &[u8] = if format == "sam" {
        b"@HD\tVN:1.6\n"
    } else {
        b""
    };
    [header, &records(format, range)].concat()
}

/// One gzip member that holds `bytes`.
fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::new(6));
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

/// The number of records in the `format` output at `path`.
fn count(format: &str, path: &Path) -> usize {
    let text = std::fs::read_to_string(path).unwrap();
    match format {
        "fq" => text.lines().count() / 4,
        "fa" => text.lines().filter(|line| line.starts_with('>')).count(),
        _ => text.lines().filter(|line| !line.starts_with('@')).count(),
    }
}

/// Run unmux over `bytes` in `format`, as a file and as stdin at 1 and 8 threads, checking each
/// run with `check` and the path of its output.
fn run_each(format: &str, bytes: &[u8], check: impl Fn(Assert, &Path)) {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join(format!("in.{format}.gz"));
    std::fs::write(&input, bytes).unwrap();
    for threads in ["1", "8"] {
        for stdin in [false, true] {
            let out = dir.path().join(format!("out.{threads}.{stdin}.{format}"));
            let mut command = Command::cargo_bin("unmux").unwrap();
            let source = if stdin { "-" } else { input.to_str().unwrap() };
            command.args([source, "--out", out.to_str().unwrap(), "--threads", threads]);
            if stdin {
                command.write_stdin(bytes.to_vec());
            }
            check(command.assert(), &out);
        }
    }
}

/// Check that a run failed on its truncated input without committing its output.
fn fails_as_truncated(assert: Assert, out: &Path) {
    assert
        .failure()
        .stderr(predicates::str::contains("is truncated"));
    assert!(!out.exists(), "a failed run commits no output");
}

/// Check that a run read all `records` of `format`, warning of an absent EOF marker block only
/// when it `warns`.
fn reads_in_full(format: &str, records: usize, warns: bool) -> impl Fn(Assert, &Path) + '_ {
    move |assert, out| {
        let output = assert.success().get_output().clone();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(
            stderr.contains("BGZF EOF marker is absent"),
            warns,
            "{stderr}"
        );
        assert_eq!(count(format, out), records, "{stderr}");
    }
}

#[test]
fn gzip_cut_at_a_sync_flush_fails_at_every_thread_count() {
    for format in FORMATS {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::new(6));
        encoder.write_all(&file(format, 0..10_000)).unwrap();
        encoder.flush().unwrap();
        run_each(format, encoder.get_ref(), fails_as_truncated);
    }
}

#[test]
fn gzip_cut_inside_its_trailer_fails_at_every_thread_count() {
    for format in FORMATS {
        let gzip = gzip(&file(format, 0..20_000));
        for missing in [1, 4, 8] {
            run_each(format, &gzip[..gzip.len() - missing], fails_as_truncated);
        }
    }
}

#[test]
fn one_record_gzip_fastq_cut_anywhere_fails() {
    let gzip = gzip(&file("fq", 0..1));
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fq.gz");
    let out = dir.path().join("out.fq");
    for len in 1..gzip.len() {
        std::fs::write(&input, &gzip[..len]).unwrap();
        Command::cargo_bin("unmux")
            .unwrap()
            .args([input.to_str().unwrap(), "--out", out.to_str().unwrap()])
            .assert()
            .failure();
        assert!(!out.exists(), "a failed run commits no output");
    }
}

#[test]
fn intact_and_concatenated_gzip_read_in_full_at_every_thread_count() {
    for format in FORMATS {
        let intact = gzip(&file(format, 0..20_000));
        let concatenated = [
            gzip(&file(format, 0..10_000)),
            gzip(&records(format, 10_000..20_000)),
        ]
        .concat();
        for bytes in [intact, concatenated] {
            run_each(format, &bytes, reads_in_full(format, 20_000, false));
        }
    }
}

/// A bgzipped FASTQ of the records `range` whose BGZF blocks each hold 1,000 whole records, with
/// the offset each block before its EOF marker block starts at.
fn bgzf_fastq(range: Range<usize>) -> (Vec<u8>, Vec<usize>) {
    let mut writer = bgzf::io::Writer::new(Vec::new());
    let mut block_starts = Vec::new();
    for start in range.clone().step_by(1_000) {
        block_starts.push(writer.get_ref().len());
        let end = (start + 1_000).min(range.end);
        writer.write_all(&records("fq", start..end)).unwrap();
        writer.flush().unwrap();
    }
    (writer.finish().unwrap(), block_starts)
}

#[test]
fn bgzf_fastq_cut_inside_a_block_fails_at_every_thread_count() {
    let (intact, block_starts) = bgzf_fastq(0..20_000);
    let middle = block_starts[block_starts.len() / 2];
    for into in [10, 30] {
        run_each("fq", &intact[..middle + into], fails_as_truncated);
    }
}

#[test]
fn bgzf_fastq_cut_at_a_block_boundary_warns_at_every_thread_count() {
    let (intact, block_starts) = bgzf_fastq(0..20_000);
    let boundary = block_starts[block_starts.len() / 2];
    let eof_block = intact.len() - 28;
    run_each("fq", &intact, reads_in_full("fq", 20_000, false));
    run_each("fq", &intact[..boundary], reads_in_full("fq", 10_000, true));
    run_each(
        "fq",
        &intact[..eof_block],
        reads_in_full("fq", 20_000, true),
    );
}

#[test]
fn concatenated_bgzf_and_gzip_fastq_read_in_full_at_every_thread_count() {
    let (first, _) = bgzf_fastq(0..10_000);
    let (second, _) = bgzf_fastq(10_000..20_000);
    let plain = gzip(&records("fq", 10_000..20_000));
    for bytes in [
        [&first[..], &second[..]].concat(),
        [&first[..], &plain[..]].concat(),
    ] {
        run_each("fq", &bytes, reads_in_full("fq", 20_000, false));
    }
}
