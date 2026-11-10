import re
import sys
from pathlib import Path

D = Path("/home/isatkaus/gridkit/uq-usecase/work-notes")


def slug(h):
    s = h.lower()
    s = re.sub(r"[^\w\s-]", "", s)
    return re.sub(r"\s+", "-", s.strip())


def anchors(text):
    heads = [slug(m.group(2)) for m in re.finditer(r"^(#{2,6})\s+(.*)$", text, re.M)]
    seen, out = {}, set()
    for h in heads:
        n = seen.get(h, 0)
        out.add(h if n == 0 else f"{h}-{n}")
        seen[h] = n + 1
    return out


def tables_ok(text):
    # \| is a valid escaped pipe inside a cell and must not be counted
    L = [l.replace(r"\|", "") for l in text.splitlines()]
    bad, i = [], 0
    while i < len(L):
        if (
            L[i].strip().startswith("|")
            and i + 1 < len(L)
            and re.fullmatch(r"\|[\s:|-]+\|", L[i + 1].strip())
        ):
            n, j = L[i].count("|"), i
            while j < len(L) and L[j].strip().startswith("|"):
                if L[j].count("|") != n:
                    bad.append(j + 1)
                j += 1
            i = j
        else:
            i += 1
    return bad


files = {p.name: p.read_text() for p in D.glob("*.md")}
anc = {n: anchors(t) for n, t in files.items()}
fail = 0
for name, text in files.items():
    bad = tables_ok(text)
    fences = sum(1 for l in text.splitlines() if l.strip().startswith("```"))
    status = "ok" if not bad and fences % 2 == 0 else "PROBLEM"
    if status != "ok":
        fail += 1
    print(
        f"{name:34s} {len(text.splitlines()):5d} lines  tables {bad or 'ok'}  "
        f"fences {'ok' if fences%2==0 else 'UNBALANCED'}"
    )
    for tgt, frag in re.findall(r"\]\((\S*?\.md)#([^\s)]+)\)", text):
        if tgt in anc and frag not in anc[tgt]:
            print(f"    BROKEN -> {tgt}#{frag}")
            fail += 1
    for frag in re.findall(r"\]\(#([^\s)]+)\)", text):
        if frag not in anc[name]:
            print(f"    BROKEN internal -> #{frag}")
            fail += 1
    for tgt in re.findall(r"\]\((\S*?\.md)\)", text):
        if not (D / tgt).exists():
            print(f"    MISSING FILE -> {tgt}")
            fail += 1
    for img in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text):
        if not (D / img).exists():
            print(f"    MISSING IMAGE -> {img}")
            fail += 1
print("\nALL OK" if not fail else f"\n{fail} problem(s)")
sys.exit(1 if fail else 0)
