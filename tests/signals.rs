//! A termination signal mid-run removes the staged outputs and the directories made for them.

use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

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

#[test]
fn a_termination_signal_removes_staged_outputs_and_exits_128_plus_the_signal() {
    for (signal, status) in [("TERM", 143), ("INT", 130), ("HUP", 129)] {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("keep.fq"), "kept\n").unwrap();
        let mut child = Command::new(env!("CARGO_BIN_EXE_unmux"))
            .current_dir(dir.path())
            .args(["--out", "new/deep/out.fq.gz"])
            .stdin(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        let mut stdin = child.stdin.take().unwrap();
        for i in 0..100 {
            writeln!(
                stdin,
                "@r{i}\n{}\n+\n{}",
                "ACGT".repeat(25),
                "I".repeat(100)
            )
            .unwrap();
        }
        stdin.flush().unwrap();
        let started = Instant::now();
        while !tree(dir.path()).iter().any(|path| path.ends_with(".tmp")) {
            assert!(
                started.elapsed() < Duration::from_secs(20),
                "no staged output"
            );
            std::thread::sleep(Duration::from_millis(20));
        }
        let killed = Command::new("kill")
            .args([format!("-{signal}"), child.id().to_string()])
            .status()
            .unwrap();
        assert!(killed.success());
        assert_eq!(child.wait().unwrap().code(), Some(status), "SIG{signal}");
        drop(stdin);
        assert_eq!(tree(dir.path()), ["keep.fq"], "SIG{signal}");
    }
}
