//! Cucumber behavior specs for unmux.
//!
//! Every feature file under `tests/features/` is an executable acceptance spec:
//! it drives the real `unmux` binary on small fixtures and asserts on its exit
//! code, stdout/stderr, output files, and (via the `noodles`-backed `the BAM
//! header of "..."` steps) BAM/CRAM header internals such as the `@PG`
//! provenance record and per-sample `@RG` read groups. The engine itself is
//! also validated by the inline unit tests in `src/lib`. The `@narrative`
//! filter in `main` is retained so any future prose-only spec can be parked as
//! documentation, but at present nothing is tagged.
//!
//! Each scenario runs in its own temporary working directory: `Given a file
//! "..." containing:` writes a fixture there, the `I run` step executes the
//! binary with that directory as the cwd, and the file-oriented `Then` steps
//! inspect what it wrote.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Output;

use assert_cmd::Command;
use cucumber::gherkin::Step;
use cucumber::{given, then, when, World};
use noodles::{bam, sam};
use tempfile::TempDir;

/// Decode a BAM file's header and re-render it as SAM header text
/// (`@HD`/`@RG`/`@PG` lines). This is what the header-assertion steps inspect;
/// it is the one thing the plain-text harness cannot do (BAM headers are
/// BGZF-framed binary), so `noodles` reads them the same way the binary wrote
/// them.
fn bam_header_text(path: &Path) -> String {
    let file = std::fs::File::open(path).unwrap_or_else(|_| panic!("open BAM file {path:?}"));
    let mut reader = bam::io::Reader::new(file);
    let header = reader
        .read_header()
        .unwrap_or_else(|e| panic!("read BAM header {path:?}: {e}"));
    let mut buf = Vec::new();
    let mut writer = sam::io::Writer::new(&mut buf);
    writer
        .write_header(&header)
        .unwrap_or_else(|e| panic!("serialize BAM header {path:?}: {e}"));
    String::from_utf8_lossy(&buf).into_owned()
}

/// Shared state across steps: the scenario's working directory and the most
/// recent invocation.
#[derive(Debug, Default, World)]
struct SassWorld {
    /// Per-scenario working directory, created lazily on first use and removed
    /// at scenario end.
    dir: Option<TempDir>,
    /// Arguments of the most recent invocation.
    args: Vec<String>,
    /// Exit code of the most recent invocation.
    code: Option<i32>,
    /// Captured stdout (lossy UTF-8).
    stdout: String,
    /// Captured stderr (lossy UTF-8).
    stderr: String,
}

impl SassWorld {
    /// The scenario's working directory, created on first use.
    fn work_dir(&mut self) -> PathBuf {
        self.dir
            .get_or_insert_with(|| TempDir::new().expect("create a temp working directory"))
            .path()
            .to_path_buf()
    }
}

fn run_unmux(dir: &Path, args: &[String]) -> Output {
    Command::cargo_bin("unmux")
        .expect("binary `unmux` builds")
        .current_dir(dir)
        .args(args)
        .output()
        .expect("unmux runs")
}

#[given(regex = r#"^a file "([^"]+)" containing:$"#)]
async fn given_file(world: &mut SassWorld, name: String, step: &Step) {
    let dir = world.work_dir();
    // Gherkin dedents the docstring but keeps the newline after the opening
    // `"""` (a leading blank line) and the indentation before the closing one;
    // trim that framing and end with one newline so the fixture is a clean
    // FASTX/SAM file.
    let content = format!("{}\n", step.docstring.clone().unwrap_or_default().trim());
    let path = dir.join(&name);
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).expect("create fixture parent directory");
    }
    fs::write(&path, content).expect("write fixture file");
}

#[given(regex = r#"^an empty file "([^"]+)"$"#)]
async fn given_empty_file(world: &mut SassWorld, name: String) {
    let dir = world.work_dir();
    let path = dir.join(&name);
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).expect("create fixture parent directory");
    }
    fs::write(&path, b"").expect("write empty fixture file");
}

#[when(regex = r"^I run `unmux ?(.*)`$")]
async fn i_run(world: &mut SassWorld, arg_line: String) {
    let dir = world.work_dir();
    let args: Vec<String> = arg_line.split_whitespace().map(str::to_owned).collect();
    let out = run_unmux(&dir, &args);
    world.code = out.status.code();
    world.stdout = String::from_utf8_lossy(&out.stdout).into_owned();
    world.stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    world.args = args;
}

#[then(regex = r"^the exit code is (\d+)$")]
async fn exit_code_is(world: &mut SassWorld, expected: i32) {
    assert_eq!(
        world.code,
        Some(expected),
        "args={:?}\nstderr:\n{}",
        world.args,
        world.stderr
    );
}

#[then(regex = r#"^stdout contains "(.*)"$"#)]
async fn stdout_contains(world: &mut SassWorld, needle: String) {
    assert!(
        world.stdout.contains(&needle),
        "stdout did not contain {needle:?}\nstdout:\n{}",
        world.stdout
    );
}

#[then(regex = r#"^stderr contains "(.*)"$"#)]
async fn stderr_contains(world: &mut SassWorld, needle: String) {
    assert!(
        world.stderr.contains(&needle),
        "stderr did not contain {needle:?}\nstderr:\n{}",
        world.stderr
    );
}

#[then(regex = r#"^a file "([^"]+)" exists$"#)]
async fn file_exists(world: &mut SassWorld, name: String) {
    let path = world.work_dir().join(&name);
    assert!(path.exists(), "expected output file {name:?} to exist");
}

#[then(regex = r#"^the file "([^"]+)" is gzip-compressed$"#)]
async fn file_is_gzip(world: &mut SassWorld, name: String) {
    let path = world.work_dir().join(&name);
    let bytes = fs::read(&path).unwrap_or_else(|_| panic!("read output file {name:?}"));
    // gzip and BGZF (the BAM/`.gz` container) both begin with the gzip magic 1f
    // 8b.
    assert!(
        bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b,
        "{name:?} is not gzip/BGZF (first bytes: {:02x?})",
        bytes.get(..2)
    );
}

#[then(regex = r#"^the file "([^"]+)" contains "(.*)"$"#)]
async fn file_contains(world: &mut SassWorld, name: String, needle: String) {
    let path = world.work_dir().join(&name);
    let bytes = fs::read(&path).unwrap_or_else(|_| panic!("read output file {name:?}"));
    let content = String::from_utf8_lossy(&bytes);
    assert!(
        content.contains(&needle),
        "{name:?} did not contain {needle:?}\ncontent:\n{content}"
    );
}

#[then(regex = r#"^stdout does not contain "(.*)"$"#)]
async fn stdout_lacks(world: &mut SassWorld, needle: String) {
    assert!(
        !world.stdout.contains(&needle),
        "stdout unexpectedly contained {needle:?}\nstdout:\n{}",
        world.stdout
    );
}

#[then(regex = r#"^the file "([^"]+)" does not contain "(.*)"$"#)]
async fn file_lacks(world: &mut SassWorld, name: String, needle: String) {
    let path = world.work_dir().join(&name);
    let bytes = fs::read(&path).unwrap_or_else(|_| panic!("read output file {name:?}"));
    let content = String::from_utf8_lossy(&bytes);
    assert!(
        !content.contains(&needle),
        "{name:?} unexpectedly contained {needle:?}\ncontent:\n{content}"
    );
}

#[then(regex = r#"^a file "([^"]+)" does not exist$"#)]
async fn file_absent(world: &mut SassWorld, name: String) {
    let path = world.work_dir().join(&name);
    assert!(!path.exists(), "file {name:?} unexpectedly exists");
}

#[then(regex = r#"^the file "([^"]+)" is larger than the file "([^"]+)"$"#)]
async fn file_larger_than(world: &mut SassWorld, bigger: String, smaller: String) {
    let dir = world.work_dir();
    let big = fs::metadata(dir.join(&bigger))
        .unwrap_or_else(|_| panic!("read output file {bigger:?}"))
        .len();
    let small = fs::metadata(dir.join(&smaller))
        .unwrap_or_else(|_| panic!("read output file {smaller:?}"))
        .len();
    assert!(
        big > small,
        "expected {bigger:?} ({big} bytes) to be larger than {smaller:?} ({small} bytes)"
    );
}

#[then(regex = r#"^the BAM header of "([^"]+)" contains "(.*)"$"#)]
async fn bam_header_contains(world: &mut SassWorld, name: String, needle: String) {
    let header = bam_header_text(&world.work_dir().join(&name));
    assert!(
        header.contains(&needle),
        "BAM header of {name:?} did not contain {needle:?}\nheader:\n{header}"
    );
}

#[then(regex = r#"^the BAM header of "([^"]+)" does not contain "(.*)"$"#)]
async fn bam_header_lacks(world: &mut SassWorld, name: String, needle: String) {
    let header = bam_header_text(&world.work_dir().join(&name));
    assert!(
        !header.contains(&needle),
        "BAM header of {name:?} unexpectedly contained {needle:?}\nheader:\n{header}"
    );
}

#[then(regex = r#"^the BAM header of "([^"]+)" has exactly (\d+) @PG lines?$"#)]
async fn bam_header_pg_count(world: &mut SassWorld, name: String, expected: usize) {
    let header = bam_header_text(&world.work_dir().join(&name));
    let count = header.lines().filter(|l| l.starts_with("@PG")).count();
    assert_eq!(
        count, expected,
        "expected {expected} @PG line(s) in {name:?}, found {count}\nheader:\n{header}"
    );
}

#[then(regex = r#"^the BAM header of "([^"]+)" records the running unmux version$"#)]
async fn bam_header_version(world: &mut SassWorld, name: String) {
    let header = bam_header_text(&world.work_dir().join(&name));
    // The @PG VN field must be the binary's own version. Bin and test share one
    // Cargo package, so the test's CARGO_PKG_VERSION is the running unmux
    // version; this keeps the assertion correct across bumps.
    let want = format!("VN:{}", env!("CARGO_PKG_VERSION"));
    let pg = header
        .lines()
        .find(|l| l.starts_with("@PG"))
        .unwrap_or_else(|| panic!("{name:?} has no @PG line\nheader:\n{header}"));
    assert!(
        pg.split('\t').any(|f| f == want),
        "expected @PG to carry {want:?} in {name:?}\n@PG: {pg}"
    );
}

fn main() {
    futures::executor::block_on(SassWorld::cucumber().filter_run_and_exit(
        "tests/features",
        |feature, _rule, scenario| {
            let parked = |tags: &[String]| tags.iter().any(|t| t == "narrative");
            !(parked(&feature.tags) || parked(&scenario.tags))
        },
    ));
}
