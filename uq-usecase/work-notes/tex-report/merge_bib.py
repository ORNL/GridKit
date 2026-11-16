#!/usr/bin/env python3
"""Safely merge the pilot's references into the main SciDAC references.bib.

Appends only the entries whose key is NOT already present (case-insensitive),
so nothing existing in the main bibliography is touched or duplicated. Writes a
timestamped-free backup (<references>.bak) first, then verifies there are no
duplicate keys afterwards.

Usage:
    python merge_bib.py /path/to/downloaded/references.bib
    python merge_bib.py /path/to/references.bib --ours uq-usecase_refs.bib

Then upload the updated references.bib back to Overleaf.
"""

import argparse
import re
import shutil
import sys
from pathlib import Path

DEFAULT_OURS = Path(__file__).resolve().parent / "uq-usecase_refs.bib"


def parse_entries(text):
    """Return [(key, full_entry_text)] using brace matching (handles nested braces)."""
    entries, pos = [], 0
    for m in re.finditer(r"@\w+\s*\{", text):
        if m.start() < pos:
            continue  # inside an already-parsed entry
        km = re.match(r"@\w+\s*\{\s*([^,\s}]+)", text[m.start():])
        key = km.group(1) if km else None
        depth, j = 0, m.end() - 1
        while j < len(text):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        entries.append((key, text[m.start():j + 1]))
        pos = j + 1
    return entries


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("references", help="path to the main references.bib to augment")
    ap.add_argument("--ours", default=str(DEFAULT_OURS), help="pilot bib (default: uq-usecase_refs.bib)")
    ap.add_argument("--dry-run", action="store_true", help="report only; do not write")
    args = ap.parse_args()

    main_path = Path(args.references)
    if not main_path.exists():
        sys.exit(f"ERROR: not found: {main_path}")

    main_text = main_path.read_text(encoding="utf-8")
    ours = parse_entries(Path(args.ours).read_text(encoding="utf-8"))
    existing = {k.lower() for k, _ in parse_entries(main_text)}

    to_add = [(k, e) for k, e in ours if k.lower() not in existing]
    skipped = [k for k, _ in ours if k.lower() in existing]

    print(f"main bib entries: {len(existing)}")
    print(f"pilot entries:    {len(ours)}")
    print(f"to append:        {len(to_add)}")
    print(f"skipped (present):{skipped}")

    if args.dry_run:
        print("dry run; nothing written.")
        return

    shutil.copy(str(main_path), str(main_path) + ".bak")
    with main_path.open("a", encoding="utf-8") as f:
        f.write("\n\n% ===== Pilot UQ use-case references (datacenter_pilot_report) =====\n\n")
        for _, e in to_add:
            f.write(e.strip() + "\n\n")

    after = parse_entries(main_path.read_text(encoding="utf-8"))
    all_keys = [k.lower() for k, _ in after]
    dups = sorted({k for k in all_keys if all_keys.count(k) > 1})
    print(f"backup written:   {main_path}.bak")
    print(f"total entries now:{len(after)}")
    print(f"duplicate keys (MUST be empty): {dups}")
    if dups:
        sys.exit("ERROR: duplicate keys introduced; restore from .bak")


if __name__ == "__main__":
    main()
