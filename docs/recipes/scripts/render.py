"""Fill a recipe's template from its results, or assemble docs/recipes.md from every recipe."""

import json
import re
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cookbook as cb

PLACEHOLDER = re.compile(r"\{\{([a-z0-9_]+)\}\}")


def fragment(template, results, out):
    text = Path(template).read_text()
    values = json.loads(Path(results).read_text())["values"]
    wanted = set(PLACEHOLDER.findall(text))
    if wanted - set(values):
        raise SystemExit(f"{template}: no value for {sorted(wanted - set(values))}")
    if set(values) - wanted:
        raise SystemExit(f"{template}: no place for {sorted(set(values) - wanted)}")
    Path(out).write_text(PLACEHOLDER.sub(lambda m: values[m.group(1)], text).strip("\n") + "\n")


def page(order, results, versions, out):
    sections = yaml.safe_load(Path(order).read_text())
    versions = json.loads(Path(versions).read_text())
    results = Path(results)
    lines = [
        "# unmux cookbook",
        "",
        ("Recipes for common unmux jobs. Each one says when to use it, the data it ran on, the exact "
        "command, and an excerpt of what unmux printed or wrote, so you can see the result before "
        "running it on your own data."),
        "",
        cb.callout("NOTE", f"Every command on this page was run with {versions['unmux']} and "
                           f"{versions['samtools']}, and every output is copied from that run. The workflow in "
                           "[`docs/recipes/`](recipes/) fetches the inputs, runs the recipes again, and "
                           "rewrites this page."),
        "",
        ("Commands are bash and run in a directory that holds their inputs. Log lines appear "
        "without their timestamps, and an excerpt that skips lines marks the gap with …. Header "
        "excerpts come from `samtools view --no-PG`, which keeps samtools' own `@PG` line out."),
        "",
        "## Contents",
        "",
    ]
    numbered, anchors, n = [], {}, 0
    for section in sections:
        items = []
        for recipe in section["recipes"]:
            n += 1
            items.append((n, recipe, json.loads((results / f"{recipe['key']}.json").read_text())))
            anchors[recipe["key"]] = f"[recipe {n}](#{cb.slug(f'{n}. {recipe['title']}')})"
        numbered.append((section, items))
    lines.append("- [Concordance at a glance](#concordance-at-a-glance)")
    for section, items in numbered:
        lines.append(f"- [{section['section']}](#{cb.slug(section['section'])})")
        for n, recipe, _ in items:
            title = f"{n}. {recipe['title']}"
            lines.append(f"  - [{title}](#{cb.slug(title)})")
    lines.append("- [Data sources](#data-sources)")
    lines.append("- [Regenerating this page](#regenerating-this-page)")
    lines += ["", "## Concordance at a glance", "",
              ("Each row counts, read by read, how unmux's calls compare with another tool's or with the "
              "truth, from the same results its recipe checks. A call is the sample, cell barcode or record "
              "each side gives a read, and a read pair counts once."), "",
              ("| Recipe | Compared with | Reads compared | Same call | Different call | unmux only | Other only "
              "| Neither | Agreement | Main reason for differences |"),
              "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
    for _, items in numbered:
        for n, recipe, result in items:
            for row in result["concordance"]:
                clean = row["different"] == row["unmux"] == row["other"] == 0
                note = ("none" + (f"; {row['note']}" if row["note"] else "")) if clean else row["note"]
                cells = [f"[{n}](#{cb.slug(f'{n}. {recipe['title']}')})", row["with"], cb.fmt(row["reads"]),
                         *(cb.fmt(row[kind]) for kind in cb.CONCORDANCE), cb.agreement(row), note]
                lines.append("| " + " | ".join(cells) + " |")
    lines += ["", ("Agreement is the share of same calls among the reads that at least one side calls; reads "
              "that neither calls are left out of it."), ""]
    used = set()
    for section, items in numbered:
        lines += ["", f"## {section['section']}", ""]
        if section.get("intro"):
            lines += [section["intro"], ""]
        for n, recipe, result in items:
            used |= {tuple(u) for u in result["used"]}
            lines += [f"### {n}. {recipe['title']}", ""]
            lines += [(results / f"{recipe['key']}.md").read_text().rstrip("\n"), ""]
    lines += ["## Data sources", "",
              ("Every download is pinned: GitHub files by commit, SRA runs by accession and spot range, "
              "NCBI sequences by accession and version."), "",
              "| Data | Pinned URL |", "| --- | --- |"]
    rows = []
    for source in cb.SOURCES:
        row = f"| {source.label} | <{source.url}> |"
        if (source.key, source.name) in used and row not in rows:
            rows.append(row)
    lines += rows
    lines += ["", "## Regenerating this page", "",
              ("A Snakemake workflow in `docs/recipes/` fetches the inputs into `docs/recipes/cache/`, runs every "
              "recipe with its checks, and rewrites this page. pixi installs the tools it runs, pinned in "
              "`pixi.toml` and `pixi.lock`; the lima comparison also needs Docker."), "",
              cb.fence("make -C docs/recipes            # build unmux, run every recipe, and write docs/recipes.md\n"
                       "make -C docs/recipes check      # fail if the page would change\n"
                       "make -C docs/recipes check CACHE=/tmp/cookbook   # the same, from an empty cache", "bash"),
              "", ("Each recipe has its own target, such as `results/uditas.json` for its checks, so "
              "`pixi run snakemake --cores 8 results/uditas.json` runs one recipe and what it depends on. "
              "Its commands are in `rules/`, its checks in `scripts/checks/`, and its prose in "
              "`templates/`; a failed check stops the run and gives the recipe."), ""]
    text = re.sub(r"\n{3,}", "\n\n", "\n".join(lines))
    text = re.sub(r"\{\{ref:([a-z_]+)\}\}", lambda m: anchors.get(m.group(1), f"recipe {m.group(1)}"), text)
    Path(out).write_text(text)


if __name__ == "__main__":
    {"fragment": fragment, "page": page}[sys.argv[1]](*sys.argv[2:])
