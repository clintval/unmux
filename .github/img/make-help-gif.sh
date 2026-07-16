#!/usr/bin/env bash
# Render .github/img/help.gif: run `unmux --help` and gently scroll the whole thing.
#
# The GIF is recorded with VHS (https://github.com/charmbracelet/vhs). The help
# text is far taller than the terminal, so help.tape defines a small shell wrapper
# that pipes the real help through a printer that advances a couple of lines at a
# time; the terminal then eases downward on its own instead of dumping every line
# at once.
#
# Requires: vhs, gifsicle, cargo. Override the binary under test with UNMUX_BIN=/path/to/unmux.

set -euo pipefail

cd "$(dirname "$0")/../.."

for tool in vhs gifsicle; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "$tool not found; install it (e.g. brew install $tool)" >&2
    exit 1
  fi
done

if [[ -z "${UNMUX_BIN:-}" ]]; then
  cargo build --release
  UNMUX_BIN="$PWD/target/release/unmux"
fi
export UNMUX_BIN

gif=".github/img/help.gif"

echo "Recording $gif with UNMUX_BIN=$UNMUX_BIN"
vhs .github/img/help.tape

# VHS emits a huge, unoptimized GIF. Shrink it without visible loss (no lossy
# pass; it speckles the low-contrast gray comments). The median-cut color method
# matters: it weights the palette by pixel frequency, so the abundant dim gray
# annotation text keeps a palette entry. gifsicle's default (diversity) chases
# distinct hues and discards that desaturated near-background gray, which erases
# the inline annotations.
echo "Optimizing $gif with gifsicle"
gifsicle -O3 --color-method median-cut --colors 12 "$gif" -o "$gif.opt"
mv "$gif.opt" "$gif"
echo "Wrote $gif ($(du -h "$gif" | cut -f1))"
