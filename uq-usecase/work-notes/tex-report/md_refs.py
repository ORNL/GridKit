#!/usr/bin/env python3
"""Citation bookkeeping for datacenter_pilot_report.md.

The markdown's "## References" list is the single source of truth. Each entry
carries its BibTeX key in a hidden tag at the end of the line:

    23. M. B. Giles, "Multilevel Monte Carlo Methods," ... <!-- bib:Giles2015 -->

convert_pilot_report.py reads the list on every run, so [N] in the text always
maps to the key of entry N. The full BibTeX entry lives in uq-usecase_refs.bib.

Usage:
    python md_refs.py check                      # validate only
    python md_refs.py renumber --dry-run         # preview renumbering
    python md_refs.py renumber                   # renumber by first appearance
    python md_refs.py renumber --prune           # ...and drop uncited entries
"""

import argparse
import re
import shutil
import sys
from pathlib import Path

from merge_bib import parse_entries

HERE = Path(__file__).resolve().parent
DEFAULT_BIB = HERE / "uq-usecase_refs.bib"

REF_HEADING_RE = re.compile(r"^## References[ \t]*$", re.MULTILINE)
ENTRY_RE = re.compile(r"^(\d+)\.\s+(.*)$")
TAG_RE = re.compile(r"<!--\s*bib:\s*([^\s>]+)\s*-->")
# Code and math spans are matched first and left alone, so "[0]" inside them is
# never treated as a citation. Group 1 holds "N" or "N, M, ..." for a citation.
CITE_RE = re.compile(
    r"```.*?```|`[^`\n]*`|\$\$.*?\$\$|\$[^$\n]*\$|\[(\d+(?:\s*,\s*\d+)*)\](?!\()",
    re.DOTALL,
)


def split_references(md):
    """Return (body, ref_section, tail): text before/inside/after '## References'."""
    m = REF_HEADING_RE.search(md)
    if not m:
        return md, "", ""
    nxt = re.search(r"^## ", md[m.end():], re.MULTILINE)
    end = m.end() + nxt.start() if nxt else len(md)
    return md[: m.start()], md[m.start(): end], md[end:]


def parse_reference_entries(ref_section):
    """Return (preamble_lines, entries); entry = dict(num, key, lines)."""
    preamble, entries = [], []
    for line in ref_section.splitlines():
        m = ENTRY_RE.match(line)
        if m:
            tag = TAG_RE.search(line)
            entries.append({"num": m.group(1), "key": tag.group(1) if tag else None,
                            "lines": [m.group(2)]})
        elif entries and line.strip():
            entries[-1]["lines"].append(line)
            tag = TAG_RE.search(line)
            if tag and not entries[-1]["key"]:
                entries[-1]["key"] = tag.group(1)
        elif not entries:
            preamble.append(line)
    return preamble, entries


def cited_numbers(body):
    """Citation numbers in order of first appearance (code/math spans ignored)."""
    seen = []
    for m in CITE_RE.finditer(body):
        if m.group(1):
            for n in re.split(r"\s*,\s*", m.group(1)):
                if n not in seen:
                    seen.append(n)
    return seen


def bib_keys(path):
    return {k for k, _ in parse_entries(Path(path).read_text(encoding="utf-8")) if k}


def check(md, bib_path=DEFAULT_BIB, main_bib_path=None):
    """Validate citations. Returns (errors, warnings, number->key map)."""
    errors, warnings = [], []
    body, ref_section, _ = split_references(md)
    if not ref_section:
        return ["no '## References' section in the markdown"], [], {}

    _, entries = parse_reference_entries(ref_section)
    num_to_key, key_to_num = {}, {}
    for e in entries:
        n, k = e["num"], e["key"]
        if n in num_to_key:
            errors.append(f"References: number {n} used twice")
        if not k:
            errors.append(f"References: entry {n} has no <!-- bib:Key --> tag")
            continue
        if k in key_to_num:
            errors.append(f"References: key {k} used by both {key_to_num[k]} and {n}")
        num_to_key[n], key_to_num[k] = k, n

    cited = cited_numbers(body)
    for n in cited:
        if n not in num_to_key and not any(e["num"] == n for e in entries):
            errors.append(f"text cites [{n}] but References has no entry {n}")

    ours = bib_keys(bib_path)
    for n, k in num_to_key.items():
        if k not in ours:
            errors.append(f"entry {n} key '{k}' is not in {Path(bib_path).name}")

    for e in entries:
        if e["num"] not in cited:
            warnings.append(f"entry {e['num']} ({e['key']}) is never cited "
                            "(not in the PDF; remove it or run 'md_refs.py renumber --prune')")

    if main_bib_path and Path(main_bib_path).exists():
        main_keys = {k.lower() for k in bib_keys(main_bib_path)}
        missing = sorted({num_to_key[n] for n in cited if n in num_to_key
                          and num_to_key[n].lower() not in main_keys})
        if missing:
            warnings.append(f"cited keys missing from {main_bib_path} (run merge_bib.py, "
                            f"then upload references.bib to Overleaf): {', '.join(missing)}")
    return errors, warnings, num_to_key


def renumber(md, prune=False):
    """Renumber references by first citation in the text; returns (new_md, old->new)."""
    body, ref_section, tail = split_references(md)
    preamble, entries = parse_reference_entries(ref_section)
    by_num = {e["num"]: e for e in entries}
    cited = [n for n in cited_numbers(body) if n in by_num]
    order = cited + ([] if prune else [e["num"] for e in entries if e["num"] not in cited])
    mapping = {old: str(i) for i, old in enumerate(order, 1)}

    def _sub(m):
        if not m.group(1):
            return m.group(0)
        nums = re.split(r"\s*,\s*", m.group(1))
        return "[" + ", ".join(mapping.get(n, n) for n in nums) + "]"

    new_body = CITE_RE.sub(_sub, body)
    lines = list(preamble)
    while lines and not lines[-1].strip():
        lines.pop()
    lines.append("")
    for old in order:
        e = by_num[old]
        lines.append(f"{mapping[old]}. {e['lines'][0]}")
        lines.extend(e["lines"][1:])
    new_ref = "\n".join(lines) + "\n"
    if tail:
        new_ref += "\n"
    return new_body + new_ref + tail, mapping


def main():
    import convert_pilot_report as cpr

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["check", "renumber"])
    ap.add_argument("--md", default=cpr.SOURCE_MD)
    ap.add_argument("--bib", default=str(DEFAULT_BIB))
    ap.add_argument("--prune", action="store_true", help="drop uncited entries (renumber)")
    ap.add_argument("--dry-run", action="store_true", help="show changes; do not write")
    args = ap.parse_args()

    md_path = Path(args.md)
    md = md_path.read_text(encoding="utf-8")
    errors, warnings, _ = check(md, args.bib, cpr.REPORT_MAIN_BIB)
    for w in warnings:
        print(f"WARNING: {w}")
    for e in errors:
        print(f"ERROR: {e}")
    if errors:
        sys.exit(1)
    if args.command == "check":
        print("citations OK")
        return

    new_md, mapping = renumber(md, prune=args.prune)
    changed = {o: n for o, n in mapping.items() if o != n}
    _, old_entries = parse_reference_entries(split_references(md)[1])
    dropped = [e["num"] for e in old_entries if e["num"] not in mapping]
    print(f"renumbered: {changed or 'nothing'}")
    print(f"dropped:    {dropped or 'nothing'}")
    if args.dry_run or new_md == md:
        print("no file written" + (" (dry run)" if args.dry_run else ""))
        return
    shutil.copy(md_path, str(md_path) + ".bak")
    md_path.write_text(new_md, encoding="utf-8")
    print(f"wrote {md_path} (backup: {md_path}.bak)")


if __name__ == "__main__":
    main()
