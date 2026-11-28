"""Shared helpers of the cookbook workflow: sources, running commands, checks, and Markdown pieces."""

import collections
import difflib
import gzip
import hashlib
import importlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
import urllib.request
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
USER_AGENT = "unmux-recipes/1.0"
BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CONCORDANCE = ("same", "different", "unmux", "other", "neither")
NO_MAPQ = (0, 1, 2, 3, 5, 6, 7, 8, 9, 10)


class Source:
    """A pinned public input: a whole file, the first records of a gzipped FASTQ, or what a fetch command writes."""

    def __init__(self, entry):
        self.key = entry["key"]
        self.name = entry["name"]
        self.url = entry["url"]
        self.label = entry["label"]
        self.sha256 = entry.get("sha256", "")
        self.records = entry.get("records")
        self.sra = entry.get("fetch")
        self.entry = entry


def load_sources(path=ROOT / "config" / "sources.yaml"):
    return [Source(entry) for entry in yaml.safe_load(Path(path).read_text())]


SOURCES = load_sources()


class CheckFailed(Exception):
    pass


class RuleFailed(Exception):
    pass


def check(condition, message):
    if not condition:
        raise CheckFailed(message)


class Run:
    def __init__(self, command, code, out, err):
        self.command = command
        self.code = code
        self.out = out
        self.err = err


class Context:
    """Runs commands the way every recipe does, and resolves the inputs a recipe uses."""

    def __init__(self, cache):
        self.cache = Path(cache).resolve()
        self.bin = self.cache / "bin"
        self.env = {
            "PATH": f"{self.bin}{os.pathsep}{os.environ.get('PATH', '')}",
            "HOME": os.environ.get("HOME", "/tmp"),
            "NO_COLOR": "1",
            "LC_ALL": "C",
            "TZ": "UTC",
        }
        self.sources = {(s.key, s.name): s for s in SOURCES}
        self.used = []

    def sh(self, cwd, command, ok=(0,)):
        proc = subprocess.run(["bash", "-c", command], cwd=cwd, env=self.env, capture_output=True, check=False)
        out = proc.stdout.decode()
        err = re.sub(r"^\[\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ ", "[", proc.stderr.decode(), flags=re.MULTILINE)
        if proc.returncode not in ok:
            raise SystemExit(f"command exited {proc.returncode} in {cwd}:\n{command}\n{err}{out[-2000:]}")
        return Run(command, proc.returncode, out, err)

    def fetch(self, key, name):
        check((key, name) in self.sources, f"{key}/{name} is not in config/sources.yaml")
        if (key, name) not in self.used:
            self.used.append((key, name))
        path = self.cache / "downloads" / key / name
        if not path.exists():
            raise SystemExit(f"{path} is missing; its fetch rule should have run first")
        return path

    def work(self, name):
        return self.cache / "work" / name

    def workdir(self, name):
        path = self.work(name)
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True)
        return path


def stage(cache, recipe, files=(), copies=()):
    """Make a recipe's work directory afresh with its inputs: (key, name, path) downloads and (source, path) copies."""
    ctx = Context(cache)
    w = ctx.workdir(recipe)
    for key, name, path in files:
        source = ctx.fetch(key, name)
        if path:
            place(source, w / path)
    for source, path in copies:
        place(source, w / path)
    return ctx, w


def record(cache, recipe, commands, templates, used, out):
    """Run a recipe's commands in order in its work directory, and write what ran and what it printed."""
    ctx = Context(cache)
    runs = {}
    for name, spec in commands.items():
        spec = {"command": spec} if isinstance(spec, str) else spec
        cwd = ctx.work(recipe) / spec.get("cwd", "")
        run = ctx.sh(cwd, spec["command"], ok=tuple(spec.get("ok", (0,))))
        runs[name] = {"command": run.command, "code": run.code, "out": run.out, "err": run.err}
    used = [list(u[:2]) for u in used]
    Path(out).write_text(json.dumps({"runs": runs, "templates": templates, "used": used}, indent=1) + "\n")


def load_runs(path):
    data = json.loads(Path(path).read_text())
    runs = {name: Run(r["command"], r["code"], r["out"], r["err"]) for name, r in data["runs"].items()}
    return runs, data["templates"], [tuple(u) for u in data["used"]]


def main(recipe, body):
    """Run a recipe's checks from the command line: run.json in, results JSON out; a failed check stops with its place."""
    run_json, out = sys.argv[1], sys.argv[2]
    cache = Path(run_json).resolve().parents[2]
    ctx = Context(cache)
    runs, templates, used = load_runs(run_json)
    ctx.used = list(used)
    try:
        values, concordance = body(ctx, ctx.work(recipe), runs, templates)
    except CheckFailed as error:
        frames = traceback.extract_tb(error.__traceback__)
        frame = next((f for f in reversed(frames) if Path(f.filename).name != "cookbook.py"), frames[-1])
        raise SystemExit(f"recipe {recipe} failed a check at {Path(frame.filename).name}:{frame.lineno}: {error}")
    Path(out).write_text(json.dumps({"values": values, "concordance": concordance, "used": ctx.used}, indent=1) + "\n")


def fetch_groups():
    """Sources fetched together: each fetch command with all the files it writes, and each plain download alone."""
    groups = {}
    for source in SOURCES:
        groups.setdefault(source.sra or (source.key, source.name), []).append(source)
    return list(groups.values())


def rule_name(group):
    return "fetch_" + re.sub(r"\W+", "_", f"{group[0].key}_{group[0].name}").strip("_").lower()


def verify(source, path):
    digest = sha256(path)
    if not source.sha256:
        raise SystemExit(f"{source.key}/{source.name} has no sha256 in config/sources.yaml; its sha256 is {digest}")
    if digest != source.sha256:
        raise SystemExit(f"{source.label} ({source.name}): sha256 {digest} does not match the pinned {source.sha256}")


def manifest_agrees(source):
    """The accession, spot range and byte range given for a source must be the ones its fetch command uses."""
    entry = source.entry
    if "accession" in entry and entry["accession"] not in source.sra and "{" not in source.sra:
        raise SystemExit(f"{source.key}/{source.name}: the fetch command does not name {entry['accession']}")
    if "spots" in entry and "--min-spot-id {} --max-spot-id {}".format(*entry["spots"]) not in source.sra:
        raise SystemExit(f"{source.key}/{source.name}: the fetch command does not take spots {entry['spots']}")
    if "bytes" in entry and " -r {}-{} ".format(*entry["bytes"]) not in source.sra:
        raise SystemExit(f"{source.key}/{source.name}: the fetch command does not take bytes {entry['bytes']}")


def download(source, partial):
    request = urllib.request.Request(source.url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=300) as response, open(partial, "wb") as out:
        if source.records:
            with gzip.GzipFile(fileobj=response) as stream:
                for i, line in enumerate(stream):
                    if i == 4 * source.records:
                        break
                    out.write(line)
        else:
            shutil.copyfileobj(response, out)


def fetch_group(cache, group):
    """Fetch the files of one group into the cache, each checked against its pinned sha256."""
    ctx = Context(cache)
    first = group[0]
    directory = ctx.cache / "downloads" / first.key
    directory.mkdir(parents=True, exist_ok=True)
    if first.sra:
        for source in group:
            manifest_agrees(source)
        staging = directory / f".staging-{rule_name(group)}"
        if staging.exists():
            shutil.rmtree(staging)
        staging.mkdir()
        ctx.sh(staging, first.sra)
        for source in group:
            verify(source, staging / source.name)
            (staging / source.name).replace(directory / source.name)
        shutil.rmtree(staging)
        return
    partial = directory / f"{first.name}.part"
    for attempt in range(1, 4):
        try:
            download(first, partial)
            break
        except OSError as error:
            if attempt == 3:
                raise SystemExit(f"{first.url}: {error}")
            print(f"retrying {first.url} after: {error}", file=sys.stderr)
            time.sleep(5 * attempt)
    verify(first, partial)
    partial.replace(directory / first.name)


def write_tools(cache, unmux, lima):
    """Put the unmux build and the lima command on the PATH that every recipe command runs with."""
    directory = Path(cache) / "bin"
    directory.mkdir(parents=True, exist_ok=True)
    shim = directory / "unmux"
    if shim.is_symlink() or shim.exists():
        shim.unlink()
    shim.symlink_to(Path(unmux).resolve())
    wrapper = directory / "lima"
    wrapper.write_text(f"#!/bin/sh\nexec {lima} \"$@\"\n")
    wrapper.chmod(0o755)


def write_versions(cache, out):
    ctx = Context(cache)
    versions = {
        "unmux": ctx.sh(ctx.cache, "unmux --version").out.strip(),
        "samtools": ctx.sh(ctx.cache, "samtools --version | head -n 1").out.strip(),
    }
    Path(out).write_text(json.dumps(versions, indent=1) + "\n")


def logged(log, step):
    """Run one step of a rule with its messages in the rule's log, and print them as well if it fails."""
    log = Path(log if isinstance(log, (str, Path)) else log[0])
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as handle:
        try:
            step(handle)
        except SystemExit as error:
            handle.write(f"{error}\n")
            handle.flush()
            print(log.read_text(), file=sys.stderr, end="")
            raise RuleFailed(str(error)) from None


def fetch_rule(cache, params, log):
    sources = {(s.key, s.name): s for s in SOURCES}
    logged(log, lambda _: fetch_group(cache, [sources[tuple(name)] for name in params.names]))


def tools_rule(cache, input, params, log):
    logged(log, lambda _: write_tools(cache, input[0], params.lima))


def versions_rule(cache, output, log):
    logged(log, lambda _: write_versions(cache, output[0]))


def stage_rule(cache, params, log):
    """Stage a recipe's inputs, then run its staging functions, given as ("module.function", subdirectory)."""

    def step(_):
        ctx, w = stage(cache, params.recipe, params.inputs, params.copies)
        for function, subdirectory in params.after:
            module, name = function.split(".")
            getattr(importlib.import_module(module), name)(ctx, w / subdirectory)

    logged(log, step)


def record_rule(cache, output, params, log):
    logged(log, lambda _: record(cache, params.recipe, params.commands, params.templates, params.used, output.run))


def python_rule(arguments, log):
    """Run a script of the workflow in an isolated Python, with its messages in the log."""

    def step(handle):
        proc = subprocess.run([sys.executable, "-I", *map(str, arguments)], stderr=subprocess.PIPE, text=True, check=False)
        handle.write(proc.stderr)
        if proc.returncode:
            raise SystemExit(proc.stderr.strip().split("\n")[-1] if proc.stderr.strip() else f"exit {proc.returncode}")

    logged(log, step)


def check_rule(input, output, log):
    python_rule([input.script, input.run, output[0]], log)


def render_rule(input, output, log):
    python_rule([ROOT / "scripts" / "render.py", "fragment", input.template, input.results, output[0]], log)


def assemble_rule(input, output, log):
    python_rule([ROOT / "scripts" / "render.py", "page", input.order, Path(output[0]).parent, input.versions, output[0]], log)


def page_rule(input, output, log):
    logged(log, lambda _: shutil.copyfile(input[0], output[0]))


def check_page_rule(input, log):
    logged(log, lambda _: compare_page(input[0], ROOT.parent / "recipes.md"))


def compare_page(rendered, committed):
    new = Path(rendered).read_text()
    old = Path(committed).read_text() if Path(committed).exists() else ""
    if new != old:
        diff = difflib.unified_diff(old.splitlines(), new.splitlines(), "docs/recipes.md", "rendered", lineterm="")
        print("\n".join(list(diff)[:60]), file=sys.stderr)
        raise SystemExit(f"{committed} differs from the rendered page; run `make -C docs/recipes` to update it")
    print(f"{committed} is up to date", file=sys.stderr)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def write_text(path, text):
    data = text.encode()
    if Path(path).suffix == ".gz":
        data = gzip.compress(data, mtime=0)
    Path(path).write_bytes(data)


def write_fastq(path, records):
    write_text(path, "".join(f"@{n}\n{s}\n+\n{q}\n" for n, s, q in records))


def write_tsv(path, header, rows):
    write_text(path, "".join("\t".join(map(str, r)) + "\n" for r in [header] + list(rows)))


def read_fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            yield header[1:].rstrip("\n"), seq, qual


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def read_fasta(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">")).upper()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def place(src, dst):
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if str(dst).endswith(".gz") and not str(src).endswith(".gz"):
        write_text(dst, Path(src).read_text())
    else:
        shutil.copyfile(src, dst)


def fmt(n):
    return f"{n:,}"


def pct(a, b):
    return f"{100 * a / b:.1f}%"


def slug(title):
    return re.sub(r"[^\w\- ]", "", title.lower()).replace(" ", "-")


def fence(text, lang):
    return f"```{lang}\n{text.rstrip(chr(10))}\n```"


def callout(kind, text):
    lines = [f"[!{kind}]"] + text.strip().split("\n")
    return "\n".join(f"> {line}".rstrip() for line in lines)


def link(text, url):
    return f"[{text}]({url})"


def details(summary, body):
    return f"<details>\n<summary>{summary}</summary>\n\n{body}\n\n</details>"


def ref(key):
    return "{{ref:" + key + "}}"


def source_links(ctx, *keys):
    return ", ".join(link(name, ctx.sources[(key, name)].url) for key, name in keys)


def fetch_command(ctx, key, name):
    return ctx.sources[(key, name)].sra


def sam_records(ctx, cwd, path, options=""):
    out = ctx.sh(cwd, f"samtools view {options} {path}").out
    records = []
    for line in out.splitlines():
        fields = line.split("\t")
        tags = {t[:2]: t[5:] for t in fields[11:]}
        records.append((fields, tags))
    return records


def by_name(ctx, cwd, path, columns=None):
    """The records of a SAM or BAM file by read name, each as the tuple of the given columns (all by default)."""
    grouped = collections.defaultdict(list)
    for fields, _ in sam_records(ctx, cwd, path):
        grouped[fields[0]].append(tuple(fields[i] for i in columns) if columns else tuple(fields))
    return {name: tuple(rows) for name, rows in grouped.items()}


def log_lines(text, keep):
    lines = [line for line in text.rstrip("\n").split("\n") if any(k in line for k in keep)]
    return "\n".join(lines)


def pair_warnings(text):
    """Keep the first tag-pair warning of each group, its count of the rest, and every other warning."""
    pair = re.compile(r"group `([^`]+)`: tags ")
    lines = [line for line in text.split("\n") if "WARN" in line]
    listed = collections.Counter(m.group(1) for m in map(pair.search, lines) if m)
    out, seen = [], set()
    for line in lines:
        m = pair.search(line)
        if not m:
            out.append(line)
        elif m.group(1) not in seen:
            seen.add(m.group(1))
            out += [line, "…"] if listed[m.group(1)] > 1 else [line]
    return "\n".join(out), max(listed.values(), default=0)


def tally(parts, total=None):
    parts = [(n, what) for n, what in parts if n]
    if len(parts) == 1 and parts[0][0] == total:
        return f"{parts[0][1]} in all {fmt(total)}"
    words = [f"{what} in {fmt(n)}" for n, what in parts]
    return words[0] if len(words) == 1 else ", ".join(words[:-1]) + " and " + words[-1]


def concord(compared_with, mine, theirs, note):
    """Count, over the reads of both sides, where unmux's call and another's agree; a call of None is no call."""
    reads = set(mine) | set(theirs)
    counts = collections.Counter()
    for read in reads:
        a, b = mine.get(read), theirs.get(read)
        counts["same" if a is not None and a == b else "different" if a is not None and b is not None
               else "unmux" if a is not None else "other" if b is not None else "neither"] += 1
    row = {"with": compared_with, "reads": len(reads), "note": note, **{kind: counts[kind] for kind in CONCORDANCE}}
    check(sum(row[kind] for kind in CONCORDANCE) == row["reads"], f"{compared_with}: the counts add up to the reads compared")
    return row


def agreement(row):
    called = row["reads"] - row["neither"]
    if not called:
        return "n/a"
    return "100%" if row["same"] == called else f"{10000 * row['same'] // called / 100:.2f}%"


def agreement_table(other, columns):
    """Per-read calls of unmux beside another tool's, one column per run of that tool."""
    labels = [("same", "Same call from both"), ("different", "Different calls"), ("unmux", "A call from unmux only"),
              ("other", f"A call from {other} only"), ("neither", "No call from either")]
    lines = ["| Reads | " + " | ".join(name for name, _ in columns) + " |",
             "| --- | " + " | ".join("---:" for _ in columns) + " |"]
    lines += [f"| {label} | " + " | ".join(fmt(row[kind]) for _, row in columns) + " |" for kind, label in labels]
    return "\n".join(lines)


def tool_version(run, tool):
    found = re.search(r"\d+\.\d+\.\d+", run.out)
    check(found is not None, f"{tool} reports its version ({run.out.strip()!r})")
    return f"{tool} {found.group(0)}"
