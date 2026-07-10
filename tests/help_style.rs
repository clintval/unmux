//! The green one-line description at the top of `--help`/`-h` and the tertiary
//! color applied to code examples and backtick-wrapped terms. These are styled
//! on the `clap::Command` at startup (`src/main.rs`), so the test drives the
//! real binary and inspects the raw ANSI in its output. Escape sequences are
//! reconstructed from clap's re-exported `anstyle` (same version as the
//! binary) rather than hard-coded, so the assertions do not depend on the
//! exact SGR byte layout.

use assert_cmd::Command;
use clap::builder::styling::{AnsiColor, Effects, Style};

/// Primary color: the one-line description (green, bold). Mirrors `TITLE` in `src/main.rs`.
const TITLE: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
/// Tertiary color: code examples and backtick terms (yellow). Mirrors `CODE` in `src/main.rs`.
const CODE: Style = AnsiColor::Yellow.on_default();
/// Footer color: the license/attribution line (magenta). Mirrors `FOOTER` in `src/main.rs`.
const FOOTER: Style = AnsiColor::Magenta.on_default();

/// The tool's one-line description; the first line of the help text.
const DESCRIPTION: &str = "Flexible record parsing and demultiplexing to FASTX/SAM/BAM/CRAM.";

/// Run the binary with `flag` and return its captured stdout (with ANSI).
/// `NO_COLOR` is cleared so these color assertions are deterministic even when
/// the test host sets it.
fn help_stdout(flag: &str) -> String {
    let output = Command::cargo_bin("unmux")
        .unwrap()
        .env_remove("NO_COLOR")
        .arg(flag)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    String::from_utf8(output).unwrap()
}

#[test]
fn long_help_opens_with_the_green_description_line() {
    let stdout = help_stdout("--help");
    let title = format!("{}{DESCRIPTION}{}\n", TITLE.render(), TITLE.render_reset());
    assert!(
        stdout.starts_with(&title),
        "--help should open with the green one-line description; got:\n{stdout}"
    );
}

#[test]
fn short_help_opens_with_the_green_description_line() {
    let stdout = help_stdout("-h");
    let title = format!("{}{DESCRIPTION}{}\n", TITLE.render(), TITLE.render_reset());
    assert!(
        stdout.starts_with(&title),
        "-h should open with the green one-line description; got:\n{stdout}"
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

#[test]
fn both_help_flags_end_with_the_colored_license_footer() {
    let footer = format!(
        "{}MIT License 2026 · Clint Valentine{}",
        FOOTER.render(),
        FOOTER.render_reset(),
    );
    for flag in ["--help", "-h"] {
        assert!(
            help_stdout(flag).contains(&footer),
            "the colored license footer should appear in {flag}"
        );
    }
}

/// Visible width of a rendered line: characters minus the ANSI SGR escapes.
fn visible_width(line: &str) -> usize {
    let mut n = 0;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\u{1b}' {
            while let Some(&nc) = chars.peek() {
                chars.next();
                if nc.is_ascii_alphabetic() {
                    break;
                }
            }
        } else {
            n += 1;
        }
    }
    n
}

#[test]
fn both_help_flags_wrap_within_eighty_columns() {
    // `--help` is hand-wrapped; `-h` uses clap's next-line layout. Neither the
    // colored nor the NO_COLOR rendering may exceed 80 rendered columns.
    for flag in ["--help", "-h"] {
        let widest = help_stdout(flag)
            .lines()
            .map(visible_width)
            .max()
            .unwrap_or(0);
        assert!(
            widest <= 80,
            "{flag} has a line of {widest} columns; must be <= 80"
        );
    }
}

/// Run the binary with `NO_COLOR` set to `value` and return its stdout.
fn no_color_stdout(flag: &str, value: &str) -> String {
    let output = Command::cargo_bin("unmux")
        .unwrap()
        .env("NO_COLOR", value)
        .arg(flag)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    String::from_utf8(output).unwrap()
}

#[test]
fn no_color_env_suppresses_all_ansi() {
    // Per https://no-color.org, any value (including empty) disables color.
    for value in ["1", "", "anything"] {
        for flag in ["--help", "-h"] {
            let out = no_color_stdout(flag, value);
            assert!(
                !out.contains('\u{1b}'),
                "NO_COLOR={value:?} {flag} must emit no ANSI escapes"
            );
            // Backticks are still stripped and the content is intact.
            assert!(
                !out.contains('`'),
                "NO_COLOR={value:?} {flag} strips backticks"
            );
            assert!(
                out.contains("Flexible record parsing and demultiplexing"),
                "NO_COLOR={value:?} {flag} keeps the description"
            );
            assert!(
                out.contains("MIT License 2026 · Clint Valentine"),
                "NO_COLOR={value:?} {flag} keeps the footer"
            );
        }
    }
}
