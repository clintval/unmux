Output .github/img/help.gif

Set Shell "zsh"
Set Theme "Dracula"
Set FontSize 21
Set Width 1200
Set Height 1074
Set Padding 20
Set Framerate 24

Hide
Type `unmux() { "$UNMUX_BIN" "$@" | { rows=$(stty size </dev/tty 2>/dev/null | awk '{print $1}'); [ -z "$rows" ] && rows=41; fill=$((rows - 2)); i=0; while IFS= read -r line || [ -n "$line" ]; do printf '%s\n' "$line"; i=$((i + 1)); if [ "$i" -eq "$fill" ]; then sleep 1.2; elif [ "$i" -gt "$fill" ]; then sleep 0.09; fi; done; }; }`
Enter
Type `clear`
Enter
Show

Sleep 800ms
Type `unmux --help`
Sleep 900ms
Enter

Sleep 33s
Sleep 1500ms
