#!/usr/bin/env bash
# Regenerate the animated README demo of `unmux --help`.
#
# Pipeline:
#   1. VHS renders the ENTIRE help into one very tall terminal image (nothing
#      scrolls off during capture).
#   2. ffmpeg pans a viewport down that image with motion blur for a smooth scroll.
#   3. libsvtav1 encodes the pan as an animated AVIF.
#
# AVIF uses the AV1 video codec, so the smooth pan stays tiny (~3MB) where a GIF or
# WebP would be tens to hundreds of MB (those formats have no motion compensation).
# Unlike an MP4, an AVIF renders and autoplays inline in a GitHub README straight
# from a committed file.
#
# Also writes help.txt: a plain-text snapshot of the help that CI diffs to decide
# whether --help changed and the demo needs re-rendering.
#
# Requires: cargo, vhs, ttyd, ffmpeg (with libsvtav1), ImageMagick.
# Override the binary under test with UNMUX_BIN=/path/to/unmux (skips cargo build).

set -euo pipefail
cd "$(dirname "$0")/../.."

command -v vhs >/dev/null || { echo "error: 'vhs' not found on PATH" >&2; exit 1; }
command -v ffmpeg >/dev/null || { echo "error: 'ffmpeg' not found on PATH" >&2; exit 1; }
magick=$(command -v magick || command -v convert || true)
[[ -n "$magick" ]] || { echo "error: ImageMagick (magick or convert) not found on PATH" >&2; exit 1; }
ffmpeg -hide_banner -encoders 2>/dev/null | grep -q libsvtav1 \
  || { echo "error: this ffmpeg has no libsvtav1 (AV1) encoder" >&2; exit 1; }

if [[ -z "${UNMUX_BIN:-}" ]]; then
  cargo build --release
  UNMUX_BIN="$PWD/target/release/unmux"
fi
export UNMUX_BIN

img=".github/img"
mkdir -p "$img"

# Plain-text help snapshot (NO_COLOR keeps it stable); CI diffs this to detect changes.
NO_COLOR=1 "$UNMUX_BIN" --help > "$img/help.txt"

# Geometry and motion knobs.
width=1200; font_size=21; padding=20; viewport_h=860
speed=180; fps=30; hold_top=1.5; hold_bottom=2.0

# Size the capture terminal to the help so the whole thing fits without scrolling.
help_lines=$(wc -l < "$img/help.txt")
tall_h=$(( (help_lines + 4) * 26 + 2 * padding ))

tape=$(mktemp -t help-tape.XXXXXX)
# No `Set FontFamily`: VHS ships its own JetBrains Mono and uses it by default, which
# renders identically everywhere. Naming a font explicitly makes VHS look for a system
# install and fall back with broken cell metrics (wide gaps, wrapping) when absent.
# The wrapper just dumps the real help into the tall terminal; $UNMUX_BIN and $@ are
# escaped so the shell inside VHS expands them at run time.
cat > "$tape" <<TAPE
Output ${img}/.tall.gif
Set Shell "zsh"
Set Theme "Dracula"
Set FontSize ${font_size}
Set Width ${width}
Set Height ${tall_h}
Set Padding ${padding}
Hide
Type \`unmux() { "\$UNMUX_BIN" "\$@"; }\`
Enter
Type "clear"
Enter
Show
Type "unmux --help"
Sleep 500ms
Enter
Sleep 3s
TAPE

echo "Rendering ${help_lines}-line help in a ${width}x${tall_h} terminal..."
vhs "$tape"
rm -f "$tape"

# The final frame holds the whole help; pull it out as a still.
ffmpeg -y -i "${img}/.tall.gif" -vf reverse -frames:v 1 "${img}/.tall.png" 2>/dev/null

# Find the last row containing text so the pan ends exactly at the footer, whatever
# the help length. Collapse each row to one pixel of average brightness after a
# threshold, then take the last row that still has ink.
bottom=$("$magick" "${img}/.tall.png" -colorspace Gray -threshold 18% -resize "1x${tall_h}!" -depth 8 txt:- \
  | sed -nE 's/^0,([0-9]+): \(([0-9]+).*/\1 \2/p' \
  | awk '$2 > 1 { last = $1 } END { print last }')
content_bottom=$(( bottom + 30 ))
maxy=$(( content_bottom - viewport_h ))
(( maxy < 0 )) && maxy=0

pan_dur=$(awk "BEGIN { print ${maxy} / ${speed} }")
total=$(awk "BEGIN { print ${hold_top} + ${pan_dur} + ${hold_bottom} }")

echo "Panning ${width}x${viewport_h} down ${maxy}px over ${pan_dur}s; encoding AVIF..."
# 240fps oversample -> tmix averages 3 sub-frames into motion blur -> 30fps output.
ffmpeg -y -framerate 240 -loop 1 -i "${img}/.tall.png" -t "$total" \
  -vf "crop=${width}:${viewport_h}:0:max(0\,min(${maxy}\,(t-${hold_top})*${speed})),tmix=frames=3:weights='1 1 1',fps=${fps},format=yuv420p" \
  -c:v libsvtav1 -crf 32 -preset 6 -g 300 -loop 0 "${img}/help.avif"

rm -f "${img}/.tall.gif" "${img}/.tall.png"
echo "Wrote ${img}/help.avif ($(du -h "${img}/help.avif" | cut -f1)) and ${img}/help.txt"
