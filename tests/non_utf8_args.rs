//! An argument that is not UTF-8 is an ordinary error, not a panic.

#![cfg(unix)]

use std::ffi::OsStr;
use std::os::unix::ffi::OsStrExt;

use assert_cmd::Command;

#[test]
fn a_non_utf8_input_path_is_a_clean_error() {
    let output = Command::new(env!("CARGO_BIN_EXE_unmux"))
        .arg(OsStr::from_bytes(b"missing\xff.fq"))
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(!stderr.contains("panicked"), "{stderr}");
    assert!(stderr.contains("missing"), "{stderr}");
}
