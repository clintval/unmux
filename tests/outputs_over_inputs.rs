//! An output that is one of the input files under another name is refused before any output exists.

#![cfg(unix)]

use std::path::Path;

use assert_cmd::Command;

const READS: &str =
    "@r1\nAAAAAAAACCCCGGGG\n+\nIIIIIIIIIIIIIIII\n@r2\nTTTTTTTTGGCCATGC\n+\nIIIIIIIIIIIIIIII\n";

const TRIMMED: &str = "@r1\nAAAA\n+\nIIII\n@r2\nTTTT\n+\nIIII\n";

/// Trim the reads of `in.fq` in `dir` to four bases into `out`, returning the exit code and stderr.
fn trim(dir: &Path, out: &str) -> (i32, String) {
    let output = Command::cargo_bin("unmux")
        .unwrap()
        .current_dir(dir)
        .args(["in.fq", "--extract", "a=0:0:4", "--template", "a"])
        .args(["--out", out])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    (output.status.code().unwrap(), stderr)
}

fn read(path: impl AsRef<Path>) -> String {
    std::fs::read_to_string(path).unwrap()
}

fn refusal(out: &str) -> String {
    format!("output \"{out}\" (--out) is input 0 (in.fq); give it a distinct path")
}

#[test]
fn an_absolute_path_to_the_input_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("in.fq"), READS).unwrap();
    let out = dir.path().join("in.fq").to_string_lossy().into_owned();
    let (code, stderr) = trim(dir.path(), &out);
    assert_eq!(code, 1, "{stderr}");
    assert!(stderr.contains(&refusal(&out)), "{stderr}");
    assert_eq!(read(dir.path().join("in.fq")), READS);
}

#[test]
fn a_symbolic_link_to_the_input_or_its_directory_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("in.fq"), READS).unwrap();
    std::os::unix::fs::symlink("in.fq", dir.path().join("link.fq")).unwrap();
    std::os::unix::fs::symlink(".", dir.path().join("here")).unwrap();
    for out in ["link.fq", "here/in.fq"] {
        let (code, stderr) = trim(dir.path(), out);
        assert_eq!(code, 1, "{out}: {stderr}");
        assert!(stderr.contains(&refusal(out)), "{stderr}");
        assert_eq!(read(dir.path().join("in.fq")), READS, "{out}");
    }
    let link = std::fs::symlink_metadata(dir.path().join("link.fq")).unwrap();
    assert!(link.file_type().is_symlink());
}

#[test]
fn a_hard_link_to_the_input_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("in.fq"), READS).unwrap();
    std::fs::hard_link(dir.path().join("in.fq"), dir.path().join("hard.fq")).unwrap();
    let (code, stderr) = trim(dir.path(), "hard.fq");
    assert_eq!(code, 1, "{stderr}");
    assert!(stderr.contains(&refusal("hard.fq")), "{stderr}");
    assert_eq!(read(dir.path().join("in.fq")), READS);
    assert_eq!(read(dir.path().join("hard.fq")), READS);
}

#[test]
fn a_case_variant_of_the_input_is_refused_only_where_it_is_the_same_file() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("in.fq"), READS).unwrap();
    let folds_case = dir.path().join("IN.FQ").exists();
    let (code, stderr) = trim(dir.path(), "IN.FQ");
    assert_eq!(read(dir.path().join("in.fq")), READS, "{stderr}");
    if folds_case {
        assert_eq!(code, 1, "{stderr}");
        assert!(stderr.contains(&refusal("IN.FQ")), "{stderr}");
    } else {
        assert_eq!(code, 0, "{stderr}");
        assert_eq!(read(dir.path().join("IN.FQ")), TRIMMED);
    }
}

#[test]
fn a_process_substitution_input_is_not_refused() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("in.fq"), READS).unwrap();
    let output = std::process::Command::new("bash")
        .current_dir(dir.path())
        .args([
            "-c",
            r#""$0" <(cat in.fq) --out out.fq --metrics-summary m.tsv"#,
        ])
        .arg(assert_cmd::cargo::cargo_bin("unmux"))
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert_eq!(read(dir.path().join("out.fq")), READS);
    assert_eq!(read(dir.path().join("in.fq")), READS);
}
