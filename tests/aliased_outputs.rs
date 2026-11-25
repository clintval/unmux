//! Outputs that are one file under two names are refused before any record is written.

#![cfg(unix)]

use std::path::Path;

use assert_cmd::Command;

const READS: &str =
    "@r1\nAAAAAAAACCCCGGGG\n+\nIIIIIIIIIIIIIIII\n@r2\nTTTTTTTTGGCCATGC\n+\nIIIIIIIIIIIIIIII\n";

fn unmux(dir: &Path, args: &[&str]) -> (i32, String) {
    std::fs::write(dir.join("reads.fq"), READS).unwrap();
    let output = Command::cargo_bin("unmux")
        .unwrap()
        .current_dir(dir)
        .arg("reads.fq")
        .args(["--group", "bc={a=AAAAAAAA,b=TTTTTTTT}"])
        .args(args)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    (output.status.code().unwrap(), stderr)
}

fn entries(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

#[test]
fn hard_linked_outputs_are_refused_and_left_unchanged() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("o");
    std::fs::create_dir(&out).unwrap();
    std::fs::write(out.join("a.fq"), "keep").unwrap();
    std::fs::hard_link(out.join("a.fq"), out.join("b.fq")).unwrap();
    let (code, stderr) = unmux(
        dir.path(),
        &["--sample-from-group", "bc", "--out", "o/%sample.fq"],
    );
    assert_eq!(code, 1, "{stderr}");
    assert!(
        stderr
            .contains("outputs \"o/a.fq\" (sample a) and \"o/b.fq\" (sample b) are the same file"),
        "{stderr}"
    );
    assert_eq!(entries(&out), ["a.fq", "b.fq"]);
    assert_eq!(std::fs::read_to_string(out.join("b.fq")).unwrap(), "keep");
}

#[test]
fn an_output_through_a_dangling_directory_link_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("o");
    std::fs::create_dir(&out).unwrap();
    std::os::unix::fs::symlink("new", out.join("L")).unwrap();
    let args = [
        "--sample",
        "a=bc::a",
        "--out",
        "o/new/%sample.fq",
        "--unassigned",
        "o/L/a.fq",
    ];
    let (code, stderr) = unmux(dir.path(), &args);
    assert_eq!(code, 1, "{stderr}");
    assert!(
        stderr.contains("(--unassigned) are the same file"),
        "{stderr}"
    );
    assert_eq!(entries(&out), ["L"]);
}

/// On a file system that ignores Unicode normalization, NFC and NFD spellings are one file.
#[cfg(target_os = "macos")]
#[test]
fn outputs_one_file_only_on_this_file_system_are_refused_and_cleaned_up() {
    let (nfc, nfd) = ("o/\u{e9}.fq", "o/e\u{301}.fq");
    let remove = format!("bc::b={nfd}");
    for args in [
        ["--unassigned", nfc, "--remove", remove.as_str()],
        ["--metrics-summary", nfc, "--metrics-per-sample", nfd],
    ] {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("o");
        std::fs::create_dir(&out).unwrap();
        let mut all = vec!["--sample", "a=bc::a", "--out", "o/%sample.fq"];
        all.extend(args);
        let (code, stderr) = unmux(dir.path(), &all);
        assert_eq!(code, 1, "{stderr}");
        assert!(
            stderr.contains("are one file on this file system"),
            "{stderr}"
        );
        assert!(entries(&out).is_empty(), "{:?}", entries(&out));
    }
}
