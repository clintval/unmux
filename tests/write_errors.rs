//! A write error on an output file fails the run and leaves no file behind, at every thread count.

#![cfg(unix)]

use std::path::Path;
use std::time::Duration;

use assert_cmd::Command;

/// Every path under `dir`, relative to it, sorted.
fn tree(dir: &Path) -> Vec<String> {
    let mut paths = Vec::new();
    let mut pending = vec![dir.to_path_buf()];
    while let Some(next) = pending.pop() {
        for entry in std::fs::read_dir(next).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                pending.push(path.clone());
            }
            let relative = path.strip_prefix(dir).unwrap();
            paths.push(relative.to_string_lossy().into_owned());
        }
    }
    paths.sort();
    paths
}

/// A poorly compressible FASTQ of `records` reads, half per barcode.
fn reads(records: usize) -> String {
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let mut random = |below: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % below) as u8
    };
    let mut fastq = String::new();
    for record in 0..records {
        let barcode = if record % 2 == 0 {
            "AAAAAAAA"
        } else {
            "TTTTTTTT"
        };
        let bases: String = (0..142)
            .map(|_| char::from(b"ACGT"[usize::from(random(4))]))
            .collect();
        let quals: String = (0..150).map(|_| char::from(b'!' + random(41))).collect();
        fastq.push_str(&format!("@r{record}\n{barcode}{bases}\n+\n{quals}\n"));
    }
    fastq
}

/// The exit code and stderr of a run whose files cannot grow past 64 of the shell's blocks.
fn unmux_under_a_file_size_limit(dir: &Path, args: &[&str]) -> (Option<i32>, String) {
    let output = Command::new("sh")
        .current_dir(dir)
        .args(["-c", "trap '' XFSZ; ulimit -f 64; exec \"$@\"", "sh"])
        .arg(env!("CARGO_BIN_EXE_unmux"))
        .args(["reads.fq", "--group", "bc={a=AAAAAAAA,b=TTTTTTTT}"])
        .args(["--group", "bc::loc=0:0:8", "--sample-from-group", "bc"])
        .args(args)
        .timeout(Duration::from_secs(60))
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    (output.status.code(), stderr)
}

/// A level 0 write tends to fail while records are still arriving, and a level 5 write after.
#[test]
fn a_write_error_fails_the_run_and_leaves_no_output_at_every_thread_count() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("reads.fq"), reads(5_000)).unwrap();
    for out in ["out/%sample.fq.gz", "out/%sample.bam"] {
        for level in ["5", "0"] {
            for threads in ["1", "4", "8"] {
                let args = ["--out", out, "--compression", level, "--threads", threads];
                let (code, stderr) = unmux_under_a_file_size_limit(dir.path(), &args);
                let case = args.join(" ");
                assert_eq!(code, Some(1), "{case}: {stderr}");
                assert!(stderr.contains("File too large"), "{case}: {stderr}");
                assert_eq!(tree(dir.path()), ["reads.fq"], "{case}");
            }
        }
    }
}
