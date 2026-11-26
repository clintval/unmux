//! An output that is one of the input files under another name is refused before any output exists.

#![cfg(unix)]

use std::path::Path;

use assert_cmd::Command;

const READS: &str =
    "@r1\nAAAAAAAACCCCGGGG\n+\nIIIIIIIIIIIIIIII\n@r2\nTTTTTTTTGGCCATGC\n+\nIIIIIIIIIIIIIIII\n";

const TRIMMED: &str = "@r1\nAAAA\n+\nIIII\n@r2\nTTTT\n+\nIIII\n";

/// Run unmux with `args` in `dir`, returning the exit code and stderr.
fn unmux(dir: &Path, args: &[&str]) -> (i32, String) {
    let output = Command::cargo_bin("unmux")
        .unwrap()
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    (output.status.code().unwrap(), stderr)
}

/// Trim the reads of `in.fq` in `dir` to four bases into `out`, returning the exit code and stderr.
fn trim(dir: &Path, out: &str) -> (i32, String) {
    unmux(
        dir,
        &[
            "in.fq",
            "--extract",
            "a=0:0:4",
            "--template",
            "a",
            "--out",
            out,
        ],
    )
}

/// Run a bash `script` in `dir` with the unmux binary as `$0`, returning the exit code and stderr.
fn bash(dir: &Path, script: &str) -> (i32, String) {
    let output = std::process::Command::new("bash")
        .current_dir(dir)
        .args(["-c", script])
        .arg(assert_cmd::cargo::cargo_bin("unmux"))
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    (output.status.code().unwrap(), stderr)
}

fn read(path: impl AsRef<Path>) -> String {
    std::fs::read_to_string(path).unwrap()
}

fn refused(out: &str, option: &str, input: &str) -> String {
    format!("output \"{out}\" ({option}) is {input}; give it a distinct path")
}

fn refusal(out: &str) -> String {
    refused(out, "--out", "input 0 (in.fq)")
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
    for out in ["hard.fq", "missing/../hard.fq"] {
        let (code, stderr) = trim(dir.path(), out);
        assert_eq!(code, 1, "{out}: {stderr}");
        assert!(stderr.contains(&refusal(out)), "{stderr}");
        assert_eq!(read(dir.path().join("in.fq")), READS, "{out}");
        assert_eq!(read(dir.path().join("hard.fq")), READS, "{out}");
        assert!(!dir.path().join("missing").exists(), "{out}");
    }
}

#[test]
fn a_case_variant_of_the_input_is_refused_only_where_it_is_the_same_file() {
    for (option, out) in [
        ("--out", "IN.FQ"),
        ("--out", "missing/../IN.FQ"),
        ("--metrics-summary", "missing/../IN.FQ"),
    ] {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("in.fq"), READS).unwrap();
        let folds_case = dir.path().join("IN.FQ").exists();
        let (code, stderr) = match option {
            "--out" => trim(dir.path(), out),
            _ => unmux(dir.path(), &["in.fq", "--out", "x.fq", option, out]),
        };
        assert_eq!(
            read(dir.path().join("in.fq")),
            READS,
            "{option} {out}: {stderr}"
        );
        if folds_case {
            assert_eq!(code, 1, "{option} {out}: {stderr}");
            let refusal = refused(out, option, "input 0 (in.fq)");
            assert!(stderr.contains(&refusal), "{stderr}");
            assert!(!dir.path().join("missing").exists(), "{option} {out}");
        } else {
            assert_eq!(code, 0, "{option} {out}: {stderr}");
            let written = read(dir.path().join("IN.FQ"));
            match option {
                "--out" => assert_eq!(written, TRIMMED),
                _ => assert!(written.contains("reads_assigned"), "{written}"),
            }
        }
    }
}

#[test]
fn a_regular_file_read_as_standard_input_is_refused() {
    for (input, shown) in [
        (None, "<stdin>"),
        (Some("-"), "<stdin>"),
        (Some("/dev/stdin"), "/dev/stdin"),
        (Some("/dev/fd/0"), "/dev/fd/0"),
    ] {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("in.fq"), READS).unwrap();
        let output = std::process::Command::new(assert_cmd::cargo::cargo_bin("unmux"))
            .current_dir(dir.path())
            .args(input)
            .args(["--extract", "a=0:0:4", "--template", "a", "--out", "in.fq"])
            .stdin(std::fs::File::open(dir.path().join("in.fq")).unwrap())
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{input:?}: {stderr}");
        let refusal = refused("in.fq", "--out", &format!("input 0 ({shown})"));
        assert!(stderr.contains(&refusal), "{stderr}");
        assert_eq!(read(dir.path().join("in.fq")), READS, "{input:?}");
    }
}

#[test]
fn an_input_read_through_another_descriptor_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("in.fq"), READS).unwrap();
    let script = r#""$0" /dev/fd/3 --extract a=0:0:4 --template a --out in.fq 3<in.fq"#;
    let (code, stderr) = bash(dir.path(), script);
    assert_eq!(code, 1, "{stderr}");
    let refusal = refused("in.fq", "--out", "input 0 (/dev/fd/3)");
    assert!(stderr.contains(&refusal), "{stderr}");
    assert_eq!(read(dir.path().join("in.fq")), READS);
}

#[test]
fn a_process_substitution_input_is_not_refused() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("in.fq"), READS).unwrap();
    let script = r#""$0" <(cat in.fq) --out out.fq --metrics-summary m.tsv"#;
    let (code, stderr) = bash(dir.path(), script);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(read(dir.path().join("out.fq")), READS);
    assert_eq!(read(dir.path().join("in.fq")), READS);
}
