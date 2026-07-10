//! The colored `--help`/`-h` header (tool name + version) and the tertiary
//! color applied to indented code examples. These are styled on the
//! `clap::Command` at startup (`src/main.rs`), so the test drives the real
//! binary and inspects the raw ANSI in its output. Escape sequences are
//! reconstructed from clap's re-exported `anstyle` (same version as the
//! binary) rather than hard-coded, so the assertions do not depend on the
//! exact SGR byte layout.

use assert_cmd::Command;
use clap::builder::styling::{AnsiColor, Effects, Style};

/// Primary color: the tool name (green, bold). Mirrors `NAME` in `src/main.rs`.
const NAME: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
/// Secondary color: the version (cyan). Mirrors `VERSION` in `src/main.rs`.
const VERSION: Style = AnsiColor::Cyan.on_default();
/// Tertiary color: indented code examples (yellow). Mirrors `CODE` in `src/main.rs`.
const CODE: Style = AnsiColor::Yellow.on_default();

/// Run the binary with `flag` and return its captured stdout (with ANSI).
fn help_stdout(flag: &str) -> String {
    let output = Command::cargo_bin("unmux")
        .unwrap()
        .arg(flag)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    String::from_utf8(output).unwrap()
}

#[test]
fn long_help_opens_with_colored_name_and_version_header() {
    let stdout = help_stdout("--help");
    let header = format!(
        "{}unmux{} {}{}{}\n\n",
        NAME.render(),
        NAME.render_reset(),
        VERSION.render(),
        env!("CARGO_PKG_VERSION"),
        VERSION.render_reset(),
    );
    assert!(
        stdout.starts_with(&header),
        "--help should open with the green name + cyan version header; got:\n{stdout}"
    );
    // The pre-existing about text follows the blank line.
    assert!(stdout.contains("Flexible record parsing and demultiplexing"));
}

#[test]
fn short_help_also_shows_the_colored_header() {
    let stdout = help_stdout("-h");
    let header = format!(
        "{}unmux{} {}{}{}\n\n",
        NAME.render(),
        NAME.render_reset(),
        VERSION.render(),
        env!("CARGO_PKG_VERSION"),
        VERSION.render_reset(),
    );
    assert!(
        stdout.starts_with(&header),
        "-h should also open with the colored header; got:\n{stdout}"
    );
}

#[test]
fn indented_code_examples_are_painted_in_the_tertiary_color() {
    let stdout = help_stdout("--help");
    let code = CODE.render().to_string();
    // A top-level quick-start example line.
    assert!(
        stdout.contains(&format!("{code}  unmux in.fq --out out.bam")),
        "top-level code examples should be painted in the tertiary color"
    );
    // A per-option example line (from --extract).
    assert!(
        stdout.contains(&format!("{code}  r=0:9:end")),
        "per-option code examples should be painted in the tertiary color"
    );
}

#[test]
fn prose_lines_are_not_painted_in_the_tertiary_color() {
    let stdout = help_stdout("--help");
    let code = CODE.render().to_string();
    // Section labels and prose start at the left margin and stay default-colored.
    assert!(
        !stdout.contains(&format!("{code}Quick start:")),
        "section labels should not be painted in the tertiary color"
    );
    assert!(
        !stdout.contains(&format!("{code}Mental model:")),
        "prose should not be painted in the tertiary color"
    );
}

#[test]
fn backtick_terms_are_stripped_of_their_backticks() {
    // Every backtick in the help text is markup around a code term; none should
    // survive into the rendered output, in either the long or short help.
    assert!(
        !help_stdout("--help").contains('`'),
        "no backticks in --help"
    );
    assert!(!help_stdout("-h").contains('`'), "no backticks in -h");
}

#[test]
fn inline_backtick_terms_are_painted_in_the_tertiary_color() {
    let code = CODE.render().to_string();
    let reset = CODE.render_reset().to_string();
    // A backtick term inside default-colored prose is painted, then the color is
    // reset so the surrounding prose stays default. `NAME=SOURCE` opens the
    // --group long help; `N=PATH` opens the --in short help.
    assert!(
        help_stdout("--help").contains(&format!("{code}NAME=SOURCE{reset}")),
        "prose backtick terms should be painted in --help"
    );
    assert!(
        help_stdout("-h").contains(&format!("{code}N=PATH{reset}")),
        "prose backtick terms should be painted in -h"
    );
}
