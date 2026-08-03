#!/usr/bin/env python3
"""Validate the generated line schema and every example document.

The schema is authored as code in docs/_ext/gridkit/line_schema.py and
generated at documentation build time; this script imports the same
builder, so the tests and the published schema cannot drift.

Checks, in order:

1. The built schema is valid JSON Schema 2020-12.
2. Every examples/EMT/Lines/*.catalog.yaml validates against the catalog
   document definition.
3. Every examples/EMT/Lines/*.line.json validates against the line schema
   and resolves: includes exist and are valid catalogs, names of each kind
   are unique across a document and its includes, the tower type and every
   conductor type and attachment reference resolve, no attachment point is
   used twice, and supplied tensions leave a positive minimum height.
4. Resolution reproduces the expected flat per-conductor coordinates
   (golden arrays, in document order).
5. Known-invalid documents are rejected (negative cases).

Exits nonzero on any failure. Requires jsonschema >= 4.18 and PyYAML.
"""

import importlib.util
import json
import math
import sys
from pathlib import Path

import yaml
from jsonschema import Draft202012Validator

REPO = Path(__file__).resolve().parents[2]
EXAMPLES = REPO / "examples" / "EMT" / "Lines"

# Load the builder module directly so the tests do not import the gridkit
# Sphinx extension package, which needs a documentation environment.
_spec = importlib.util.spec_from_file_location(
    "line_schema", REPO / "docs" / "_ext" / "gridkit" / "line_schema.py"
)
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
build_schema = _module.build_schema

failures = []


def check(label, ok, detail=""):
    print(f"{'ok  ' if ok else 'FAIL'} {label}" + (f": {detail}" if detail else ""))
    if not ok:
        failures.append(label)


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


schema = build_schema()
try:
    Draft202012Validator.check_schema(schema)
    check("generated line schema is valid 2020-12", True)
except Exception as e:
    check("generated line schema is valid 2020-12", False, str(e))

line_validator = Draft202012Validator(schema)
# Catalog files validate against the catalog document definition inside the
# same schema.
catalog_validator = Draft202012Validator(
    {
        "$schema": schema["$schema"],
        "$ref": "#/$defs/catalog-document",
        "$defs": schema["$defs"],
    }
)


def schema_errors(validator, doc):
    return [e.json_path + ": " + e.message for e in validator.iter_errors(doc)]


# ---- reference resolver -------------------------------------------------


MU0 = 4e-7 * math.pi


def resolve(doc, doc_path):
    """Resolve a line document against its catalogs.

    Returns (rows, errors). Rows are flat per-conductor records in document
    order: x, h, phase, circuit, tension, and the conductor type data with
    relative material properties scaled to absolute SI. These are the
    quantities the parameter models consume; type names do not survive
    resolution.
    """
    errs = []
    tables = {"conductors": {}, "towers": {}}

    def merge(sections, origin):
        for kind, table in tables.items():
            for name, data in (sections.get(kind) or {}).items():
                if name in table:
                    errs.append(f"{kind} name collision: {name} ({origin})")
                else:
                    table[name] = data

    merge(doc.get("catalog", {}), "document catalog")
    for inc in doc.get("include", []):
        inc_path = doc_path.parent / inc
        if not inc_path.is_file():
            errs.append(f"include not found: {inc}")
            continue
        catalog = yaml.safe_load(inc_path.read_text(encoding="utf-8"))
        inc_errs = schema_errors(catalog_validator, catalog)
        if inc_errs:
            errs.append(f"include {inc} is not a valid catalog: "
                        + "; ".join(inc_errs))
            continue
        merge(catalog, inc)

    for name, ctype in tables["conductors"].items():
        radius = ctype["radius"]
        if radius.get("inner", 0.0) >= radius["outer"]:
            errs.append(f"conductor type {name}: inner radius must be "
                        f"below the outer radius")

    tower = tables["towers"].get(doc["tower"])
    if tower is None:
        errs.append(f"unknown tower type: {doc['tower']}")

    rows = []
    used = set()
    for i, entry in enumerate(doc["conductors"]):
        ctype = tables["conductors"].get(entry["type"])
        if ctype is None:
            errs.append(f"conductors[{i}] references unknown conductor "
                        f"type: {entry['type']}")
        if tower is not None:
            point = tower["attachments"].get(entry["at"])
            if point is None:
                errs.append(f"conductors[{i}] references unknown "
                            f"attachment: {entry['at']}")
            elif entry["at"] in used:
                errs.append(f"attachment used more than once: {entry['at']}")
            else:
                used.add(entry["at"])
        else:
            point = None
        if ctype is None or point is None:
            continue
        rows.append({
            "x": point["x"],
            "h": point["h"],
            "phase": entry["phase"],
            "circuit": entry.get("circuit", 1),
            "tension": entry.get("tension"),
            "radius": {
                "outer": ctype["radius"]["outer"],
                "inner": ctype["radius"].get("inner", 0.0),
            },
            "conductivity": ctype["conductivity"],
            "permeability": ctype.get("permeability", 1.0) * MU0,
            "weight": ctype["weight"],
        })

    span = doc["path"]["span"]
    for row in rows:
        if row["tension"] is None:
            continue
        a = row["tension"] / row["weight"]
        try:
            sag = a * (math.cosh(span / (2.0 * a)) - 1.0)
        except OverflowError:
            sag = math.inf
        if row["h"] - sag <= 0:
            errs.append(f"tension leaves nonpositive minimum height at "
                        f"x={row['x']}")
    return rows, errs


# ---- catalogs and examples ----------------------------------------------

catalogs = sorted(EXAMPLES.glob("*.catalog.yaml"))
check("found a catalog file", bool(catalogs))
for f in catalogs:
    doc = yaml.safe_load(f.read_text(encoding="utf-8"))
    check(f.name, not schema_errors(catalog_validator, doc))

# Expected (x, h) per conductor, in document order.
GOLDEN = {
    "69kv-wood-pole.line.json": (
        [-1.0668, 1.0668, -1.0668, 0.0],
        [14.9352, 13.7160, 12.4968, 10.3632],
    ),
    "138kv-delta.line.json": (
        [-3.0, 3.0, 0.0, 0.0],
        [17.0, 17.0, 21.5, 25.0],
    ),
    "345kv-horizontal.line.json": (
        [-7.4285, -6.9715, -0.2285, 0.2285, 6.9715, 7.4285, -4.3, 4.3],
        [24.0, 24.0, 24.0, 24.0, 24.0, 24.0, 29.5, 29.5],
    ),
    "500kv-double-circuit.line.json": (
        [-8.7285, -8.2715, -10.0285, -9.5715, -8.7285, -8.2715,
         8.2715, 8.7285, 9.5715, 10.0285, 8.2715, 8.7285, -5.5, 5.5],
        [30.0, 30.0, 38.0, 38.0, 46.0, 46.0,
         30.0, 30.0, 38.0, 38.0, 46.0, 46.0, 53.0, 53.0],
    ),
    "765kv-horizontal.line.json": (
        [-13.9285, -13.4715, -13.9285, -13.4715, -0.2285, 0.2285,
         -0.2285, 0.2285, 13.4715, 13.9285, 13.4715, 13.9285, -9.8, 9.8],
        [31.7715, 31.7715, 32.2285, 32.2285, 31.7715, 31.7715,
         32.2285, 32.2285, 31.7715, 31.7715, 32.2285, 32.2285, 40.0, 40.0],
    ),
}

examples = sorted(EXAMPLES.glob("*.line.json"))
check("found the example documents", len(examples) == len(GOLDEN))
for f in examples:
    doc = load(f)
    errs = schema_errors(line_validator, doc)
    rows, res_errs = resolve(doc, f)
    errs += res_errs
    check(f.name, not errs, "; ".join(errs))
    if errs or f.name not in GOLDEN:
        continue
    x, h = GOLDEN[f.name]
    ok = ([row["x"] for row in rows] == x
          and [row["h"] for row in rows] == h)
    check(f"{f.name} resolves to the expected coordinates", ok)

# ---- negative cases -----------------------------------------------------


def mutated(base, fn):
    doc = json.loads(json.dumps(base))
    fn(doc)
    return doc


local_doc = load(EXAMPLES / "69kv-wood-pole.line.json")
include_doc = load(EXAMPLES / "345kv-horizontal.line.json")

NEGATIVE = [
    ("path with both length and points", local_doc, lambda d: d["path"].update(
        {"points": [{"latitude": 0, "longitude": 0},
                    {"latitude": 1, "longitude": 1}]})),
    ("path with neither length nor points", local_doc,
     lambda d: d["path"].pop("length")),
    ("path without span", local_doc, lambda d: d["path"].pop("span")),
    ("legacy phase label s", local_doc,
     lambda d: d["conductors"][0].update({"phase": "s"})),
    ("inline conductor data on an entry", include_doc,
     lambda d: d["conductors"][0].update({"radius": {"outer": 0.01}})),
    ("entry without an attachment point", local_doc,
     lambda d: d["conductors"][0].pop("at")),
    ("entry without a type", include_doc,
     lambda d: d["conductors"][0].pop("type")),
    ("tower as inline geometry object", local_doc, lambda d: d.update(
        {"tower": {"x": [0.0], "height": [10.0]}})),
    ("zero tension", include_doc,
     lambda d: d["conductors"][0].update({"tension": 0.0})),
    ("uppercase type name in the catalog", local_doc, lambda d: d["catalog"][
        "conductors"].update({"Drake": d["catalog"]["conductors"]["phase-acsr"]})),
    ("empty catalog section", local_doc,
     lambda d: d["catalog"].update({"towers": {}})),
    ("radius without outer", local_doc, lambda d: d["catalog"]["conductors"][
        "phase-acsr"].update({"radius": {"inner": 0.005}})),
    ("absolute conductor permeability", local_doc,
     lambda d: d["catalog"]["conductors"]["phase-acsr"].update(
         {"permeability": 1.2566370614e-6})),
    ("absolute earth permittivity", local_doc,
     lambda d: d["earth"].update({"permittivity": 8.8541878128e-12})),
]

for label, base, fn in NEGATIVE:
    check(f"rejects {label}", not line_validator.is_valid(mutated(base, fn)))

BAD_CATALOGS = [
    ("catalog file without a header", {"conductors": {
        "a1": {"radius": {"outer": 0.01}, "conductivity": 3.5e7,
               "weight": 16.0}}}),
    ("catalog file without sections", {"catalog": "1.0", "name": "empty"}),
    ("catalog file with a stray section", {"catalog": "1.0", "earths": {}}),
]

for label, doc in BAD_CATALOGS:
    check(f"rejects {label}", not catalog_validator.is_valid(doc))

# Loader-level negatives, resolved against the real example paths.
SEMANTIC = [
    ("unknown conductor type", include_doc,
     lambda d: d["conductors"][0].update({"type": "nope"}), "unknown conductor"),
    ("unknown tower type", include_doc,
     lambda d: d.update({"tower": "nope"}), "unknown tower"),
    ("unknown attachment point", include_doc,
     lambda d: d["conductors"][0].update({"at": "nowhere"}),
     "unknown attachment"),
    ("attachment used twice", include_doc,
     lambda d: d["conductors"][1].update({"at": "left-1"}), "more than once"),
    ("local/include name collision", include_doc, lambda d: d.update(
        {"catalog": {"conductors": {"drake-acsr": {
            "radius": {"outer": 0.01}, "conductivity": 3.5e7,
            "weight": 16.0}}}}), "collision"),
    ("missing include", include_doc,
     lambda d: d.update({"include": ["missing.catalog.yaml"]}), "not found"),
    ("tension exceeding the attachment height", include_doc,
     lambda d: d["conductors"][0].update({"tension": 2000.0}),
     "minimum height"),
    ("inner radius above the outer radius", local_doc,
     lambda d: d["catalog"]["conductors"]["phase-acsr"].update(
         {"radius": {"outer": 0.00915, "inner": 0.01}}), "inner radius"),
]

for label, base, fn, needle in SEMANTIC:
    _, errs = resolve(mutated(base, fn), EXAMPLES / "345kv-horizontal.line.json")
    check(f"rejects {label}", any(needle in e for e in errs))

print()
if failures:
    print(f"{len(failures)} failure(s)")
    sys.exit(1)
print("all checks passed")
