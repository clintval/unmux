"""The tag-collision check: the 10x and SPLiT-seq commands refused before they read a record or write a file."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cookbook as cb


def snapshot(w):
    return {p.relative_to(w).as_posix(): cb.sha256(p) if p.is_file() else "directory" for p in sorted(w.rglob("*"))}


def refused(ctx, w, command, expected):
    before = snapshot(w)
    run = ctx.sh(w, command, ok=(1,))
    error = cb.log_lines(run.err, ["ERROR"])
    cb.check(expected in error, f"the refusal says: {expected}")
    cb.check(snapshot(w) == before, "the refused run wrote and changed no file")
    return error


def body(ctx, w, runs, templates):
    error = refused(ctx, w / "tenx", templates["tenx"], "with mode=all this silently mis-assigns reads")
    error2 = refused(ctx, w / "splitseq", templates["indel"], "only through an indel")
    values = {
        "command": cb.fence(templates["tenx"], "bash"),
        "error": cb.fence(error, "text"),
        "indel": cb.details("The command", cb.fence(templates["indel"], "bash")),
        "error2": cb.fence(error2, "text"),
    }
    return values, []


cb.main("collision", body)
