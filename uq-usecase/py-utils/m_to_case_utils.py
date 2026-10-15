"""Build GridKit PhasorDynamics `case.json` files from TAMU source data.

Provenance (verified 2026-09-29, see work-notes/m_to_case_helper.md):

    ACTIVSg200.aux           PowerWorld power-flow case  ->  bus init, gen
                             dispatch, loads, shunts
    ACTIVSg200_dynamics.dyr  PSS/E dynamics              ->  Genrou, SexsPti,
                                                             Tgov1 parameters
    case_ACTIVSg200.m        MATPOWER export of a DIFFERENT operating point;
                             NOT the source of illinois.case.json

`Vnom` is deliberately never written: it is not a `LoadZIP` JSON parameter.
GridKit derives it from the bus initial voltage and logs an error for any
unrecognised key.
"""

from __future__ import annotations

import copy
import json
import math
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------
# case configuration
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class CaseSpec:
    """Everything that differs between GridKit cases."""

    name: str
    case_dir: Path
    case_json: str
    solver_json: str
    aux: Path | None = None
    dyr: Path | None = None
    template_m: Path | None = None
    load_class: str = "LoadZIP"
    exciter_class: str = "SexsPti"
    governor_class: str = "Tgov1"
    companion_suffixes: dict[str, str] = field(
        default_factory=lambda: {"Tgov1": "_tgov1", "SexsPti": "_sexs_pti"}
    )

    @property
    def case_json_path(self) -> Path:
        return self.case_dir / self.case_json

    @property
    def solver_json_path(self) -> Path:
        return self.case_dir / self.solver_json


_GRIDKIT = Path.home() / "gridkit"
_TAMU = Path("/kfs2/projects/scidac/scidac-data")

# Source tree, not build/: the build copy is a CMake artifact of this file.
ILLINOIS = CaseSpec(
    name="illinois",
    case_dir=_GRIDKIT / "examples/PhasorDynamics/Large/Illinois",
    case_json="illinois.case.json",
    solver_json="illinois.solver.json",
    aux=_TAMU / "ACTIVSg200/raw-tamu-data/ACTIVSg200.aux",
    dyr=_TAMU / "ACTIVSg200/raw-tamu-data/ACTIVSg200_dynamics.dyr",
    template_m=_TAMU / "ACTIVSg200/raw-tamu-data/case_ACTIVSg200.m",
    load_class="LoadZIP",
    exciter_class="SexsPti",
)

HAWAII40 = CaseSpec(
    name="hawaii40",
    case_dir=_GRIDKIT / "examples/PhasorDynamics/Medium/Hawaii",
    case_json="hawaii.case.json",
    solver_json="hawaii.solver.json",
    template_m=_TAMU / "Hawaii40/raw-tamu-data/Hawaii40_20231026.m",
    load_class="LoadZ",
    exciter_class="Ieeet1",
    companion_suffixes={"Tgov1": "_tgov1", "Ieeet1": "_ieeet1"},
)

CASES = {c.name: c for c in (ILLINOIS, HAWAII40)}


# --------------------------------------------------------------------------
# PowerWorld .aux parser
# --------------------------------------------------------------------------

_AUX_SECTION = r"DATA \(\s*{name}\s*,\s*\[(?P<hdr>.*?)\]\s*\)\s*\{{(?P<body>.*?)\n\}}"


def parse_aux_section(raw: str, name: str) -> pd.DataFrame:
    """Parse every `DATA (<name>, [...]) { ... }` block into one DataFrame.

    Some object types (e.g. Branch) appear in more than one block; all are
    concatenated. Rows whose field count does not match the header are
    dropped and reported via the `.attrs["dropped"]` counter.
    """
    pat = re.compile(_AUX_SECTION.format(name=re.escape(name)), re.S)
    frames, dropped = [], 0
    for m in pat.finditer(raw):
        header = [
            c.strip() for c in m.group("hdr").replace("\n", "").split(",") if c.strip()
        ]
        rows = []
        for line in m.group("body").splitlines():
            line = line.strip()
            if not line or line.startswith("//"):
                continue
            fields = re.findall(r'"[^"]*"|\S+', line)
            if len(fields) != len(header):
                dropped += 1
                continue
            rows.append([f[1:-1].strip() if f.startswith('"') else f for f in fields])
        if rows:
            frames.append(pd.DataFrame(rows, columns=header))
    if not frames:
        raise ValueError(f"no DATA ({name}, ...) block found in .aux")
    out = pd.concat(frames, ignore_index=True)
    out.attrs["dropped"] = dropped
    return out


def parse_aux(
    path: str | Path, sections=("Bus", "Gen", "Load", "Shunt", "Branch")
) -> dict[str, pd.DataFrame]:
    """Parse the requested object tables out of a PowerWorld `.aux` file."""
    raw = Path(path).read_text(errors="replace")
    return {s: parse_aux_section(raw, s) for s in sections}


def _num(df: pd.DataFrame, col: str) -> pd.Series:
    return pd.to_numeric(df[col], errors="coerce")


# --------------------------------------------------------------------------
# PSS/E .dyr parser
# --------------------------------------------------------------------------

# PSS/E GENROU record order. GENROU carries a single subtransient reactance,
# so Xqpp is set equal to Xdpp when building the JSON device.
GENROU_ORDER = [
    "Tdop",
    "Tdopp",
    "Tqop",
    "Tqopp",
    "H",
    "D",
    "Xd",
    "Xq",
    "Xdp",
    "Xqp",
    "Xdpp",
    "Xl",
    "S10",
    "S12",
]
_DYR_RECORD = re.compile(r"\s*(\d+)\s+'([A-Z0-9 ]+)'\s+(\S+)\s*(.*)", re.S)


def parse_dyr(
    path: str | Path, allow_multi_unit: bool = False
) -> dict[tuple[int, str], list[float]]:
    """Return `{(bus, model): [values...]}` for every record in a `.dyr`.

    KNOWN LIMITATION: the key ignores the PSS/E unit id, so a bus carrying two
    machines of the same model would silently keep only the last one. ACTIVSg200
    is one unit per bus, but the 2000- and 10k-bus cases are not. Rather than
    fail quietly, this raises; pass `allow_multi_unit=True` to override once
    proper `(bus, model, unit)` keying is implemented.
    """
    out: dict[tuple[int, str], list[float]] = {}
    seen: dict[tuple[int, str], list[str]] = {}
    for chunk in Path(path).read_text(errors="replace").split("/"):
        chunk = chunk.strip()
        if not chunk:
            continue
        m = _DYR_RECORD.match(chunk)
        if not m:
            continue
        bus, model = int(m.group(1)), m.group(2).strip()
        unit = m.group(3).strip().strip("'")
        seen.setdefault((bus, model), []).append(unit)
        out[(bus, model)] = [float(x) for x in m.group(4).split()]

    collisions = {k: v for k, v in seen.items() if len(v) > 1}
    if collisions and not allow_multi_unit:
        head = list(collisions.items())[:5]
        raise ValueError(
            f"{len(collisions)} (bus, model) collision(s) in {Path(path).name}: "
            f"{head}. Multiple units per bus are not yet supported -- the dyr key "
            f"drops the PSS/E unit id. Implement (bus, model, unit) keying, or "
            f"pass allow_multi_unit=True to accept last-record-wins."
        )
    return out


def genrou_params(dyr_values: list[float], mva: float, p0: float, q0: float) -> dict:
    """PSS/E GENROU record -> GridKit `Genrou.params`."""
    v = dict(zip(GENROU_ORDER, dyr_values))
    return {
        "p0": p0,
        "q0": q0,
        "H": v["H"],
        "D": v["D"],
        "Ra": 0.0,
        "Tdop": v["Tdop"],
        "Tdopp": v["Tdopp"],
        "Tqop": v["Tqop"],
        "Tqopp": v["Tqopp"],
        "Xd": v["Xd"],
        "Xdp": v["Xdp"],
        "Xdpp": v["Xdpp"],
        "Xq": v["Xq"],
        "Xqp": v["Xqp"],
        "Xqpp": v["Xdpp"],
        "Xl": v["Xl"],
        "S10": v["S10"],
        "S12": v["S12"],
        "mva": mva,
    }


def sexs_params(dyr_values: list[float]) -> dict:
    """PSS/E SEXS `[TA/TB, TB, K, TE, Emin, Emax]` -> GridKit `SexsPti.params`."""
    ta_over_tb, tb, k, te, emin, emax = dyr_values[:6]
    return {
        "Ta": ta_over_tb * tb,
        "Tb": tb,
        "Te": te,
        "K": k,
        "Efdmax": emax,
        "Efdmin": emin,
    }


def tgov1_params(dyr_values: list[float], mva: float) -> dict:
    """PSS/E TGOV1 `[R, T1, Vmax, Vmin, T2, T3, Dt]` -> GridKit `Tgov1.params`."""
    r, t1, vmax, vmin, t2, t3, dt = dyr_values[:7]
    return {
        "Trate": mva,
        "R": r,
        "T1": t1,
        "T2": t2,
        "T3": t3,
        "Pvmin": vmin,
        "Pvmax": vmax,
        "Dt": dt,
    }


# --------------------------------------------------------------------------
# case.json inspection
# --------------------------------------------------------------------------


def load_case(path: str | Path) -> dict:
    with open(path) as fh:
        return json.load(fh)


def devices(case: dict, cls: str | None = None) -> list[dict]:
    devs = case.get("devices", [])
    return [d for d in devs if cls is None or d.get("class") == cls]


def is_shunt(dev: dict) -> bool:
    """Shunt/capacitor rather than real load."""
    p = dev.get("params", {})
    return str(dev.get("id", "")).startswith("shunt_") or (
        float(p.get("Pnom", 0.0)) == 0.0 and float(p.get("Qnom", 0.0)) < 0.0
    )


def audit_case_json(case: dict) -> pd.DataFrame:
    """Counts that the docs assert, computed rather than hardcoded."""
    from collections import Counter

    cls = Counter(d.get("class") for d in case.get("devices", []))
    lz = devices(case, "LoadZIP")
    real = [d for d in lz if not is_shunt(d)]
    per_bus = Counter(d["ports"]["bus"] for d in real)
    multiplicity = Counter(per_bus.values())

    rows = [("buses", len(case.get("buses", [])))]
    rows += [(f"class:{c}", n) for c, n in sorted(cls.items())]
    rows += [
        ("LoadZIP total", len(lz)),
        ("LoadZIP real", len(real)),
        ("LoadZIP shunt", len(lz) - len(real)),
        ("load buses", len(per_bus)),
        ("multi-ZIP buses", sum(n for k, n in multiplicity.items() if k > 1)),
        ('"Vnom" keys (must be 0)', sum("Vnom" in d.get("params", {}) for d in lz)),
    ]
    rows += [(f"buses with {k} ZIP", n) for k, n in sorted(multiplicity.items())]
    return pd.DataFrame(rows, columns=["quantity", "value"])


def companion_report(case: dict, spec: CaseSpec) -> pd.DataFrame:
    """One row per Genrou: which companion devices exist for it."""
    by_cls = {c: {d["id"] for d in devices(case, c)} for c in spec.companion_suffixes}
    rows = []
    for g in devices(case, "Genrou"):
        gid = g["id"]
        row = {"genrou": gid}
        for cls, suffix in spec.companion_suffixes.items():
            row[cls] = (gid + suffix) in by_cls[cls] or gid in by_cls[cls]
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# structural / numerical diff
# --------------------------------------------------------------------------


def diff_case_json(ref: dict, new: dict, tol: float = 1e-7) -> pd.DataFrame:
    """Compare two cases device-by-device and bus-by-bus.

    Returns one row per differing quantity. An empty frame means the two
    cases agree to `tol`.
    """
    rows = []

    ref_bus = {b["number"]: b for b in ref.get("buses", [])}
    new_bus = {b["number"]: b for b in new.get("buses", [])}
    for num in sorted(set(ref_bus) | set(new_bus)):
        if num not in ref_bus or num not in new_bus:
            rows.append(
                {
                    "kind": "bus",
                    "id": num,
                    "param": "<presence>",
                    "ref": num in ref_bus,
                    "new": num in new_bus,
                    "delta": None,
                }
            )
            continue
        for k in ("Vr", "Vi"):
            a = ref_bus[num]["init"][k]
            b = new_bus[num]["init"][k]
            if abs(a - b) > tol:
                rows.append(
                    {
                        "kind": "bus",
                        "id": num,
                        "param": k,
                        "ref": a,
                        "new": b,
                        "delta": b - a,
                    }
                )

    ref_dev = {(d["class"], d["id"]): d for d in ref.get("devices", [])}
    new_dev = {(d["class"], d["id"]): d for d in new.get("devices", [])}
    for key in sorted(set(ref_dev) | set(new_dev)):
        cls, did = key
        if key not in ref_dev or key not in new_dev:
            rows.append(
                {
                    "kind": cls,
                    "id": did,
                    "param": "<presence>",
                    "ref": key in ref_dev,
                    "new": key in new_dev,
                    "delta": None,
                }
            )
            continue
        pa = ref_dev[key].get("params", {})
        pb = new_dev[key].get("params", {})
        for p in sorted(set(pa) | set(pb)):
            a, b = pa.get(p), pb.get(p)
            if a is None or b is None:
                rows.append(
                    {
                        "kind": cls,
                        "id": did,
                        "param": p,
                        "ref": a,
                        "new": b,
                        "delta": None,
                    }
                )
            elif isinstance(a, bool) or isinstance(b, bool):
                if a != b:
                    rows.append(
                        {
                            "kind": cls,
                            "id": did,
                            "param": p,
                            "ref": a,
                            "new": b,
                            "delta": None,
                        }
                    )
            elif abs(float(a) - float(b)) > tol:
                rows.append(
                    {
                        "kind": cls,
                        "id": did,
                        "param": p,
                        "ref": a,
                        "new": b,
                        "delta": float(b) - float(a),
                    }
                )

    return pd.DataFrame(rows, columns=["kind", "id", "param", "ref", "new", "delta"])


def make_creator_overlay(shipped: dict, rebuilt: dict, tol: float = 1e-7) -> dict:
    """Deviations of the shipped case from a pure source-data rebuild.

    Applying this overlay to `rebuilt` reproduces `shipped`. It is generated,
    never hand-written, so it stays in sync if GridKit ships a new case.
    """
    d = diff_case_json(rebuilt, shipped, tol=tol)
    d = d[d["param"] != "<presence>"]
    overlay: dict[str, dict] = {}
    for _, r in d.iterrows():
        overlay.setdefault(r["id"], {})[r["param"]] = r["new"]
    return overlay


def apply_overlay(case: dict, overlay: dict) -> dict:
    """Apply a creator overlay in place and return the case."""
    by_id = {d["id"]: d for d in case.get("devices", [])}
    for did, params in overlay.items():
        if did in by_id:
            by_id[did].setdefault("params", {}).update(params)
    return case


# --------------------------------------------------------------------------
# rebuild case.json from .aux + .dyr
# --------------------------------------------------------------------------


def build_case_from_aux_dyr(
    spec: CaseSpec,
    base_case: dict,
    aux: dict[str, pd.DataFrame],
    dyr: dict[tuple[int, str], list[float]],
    creator_overlay: dict | None = None,
    baseMVA: float = 100.0,
) -> tuple[dict, dict]:
    """Rebuild a case from source data, keeping the base file's topology.

    Bus init, generator dispatch and dynamics, loads and shunts are all taken
    from `aux`/`dyr`. Branch, BusFault, signals and monitors are inherited
    from `base_case` unchanged.

    With `creator_overlay=None` the result is pure source data. Pass an
    overlay from `make_creator_overlay` to reproduce a shipped case exactly.
    """
    case = copy.deepcopy(base_case)
    counts = {
        "bus": 0,
        "genrou": 0,
        "sexs": 0,
        "tgov1": 0,
        "load": 0,
        "shunt": 0,
        "missing_dyr": [],
    }

    bus = aux["Bus"]
    vmag = dict(zip(_num(bus, "BusNum").astype(int), _num(bus, "BusPUVolt")))
    vang = dict(zip(_num(bus, "BusNum").astype(int), _num(bus, "BusAngle")))
    for b in case.get("buses", []):
        n = b["number"]
        if n in vmag:
            rad = math.radians(vang[n])
            b["init"]["Vr"] = vmag[n] * math.cos(rad)
            b["init"]["Vi"] = vmag[n] * math.sin(rad)
            counts["bus"] += 1

    gen = aux["Gen"].copy()
    gen["_bus"] = _num(gen, "BusNum").astype(int)
    gen["_rank"] = gen.groupby("_bus").cumcount() + 1
    gen["_gid"] = gen["_bus"].astype(str) + "_" + gen["_rank"].astype(str)
    gen["_online"] = gen["GenStatus"].str.strip().str.lower().eq("closed")
    gen_by_id = {r["_gid"]: r for _, r in gen.iterrows()}

    exc = spec.exciter_class
    suffixes = spec.companion_suffixes
    kept: list[dict] = []
    for d in case["devices"]:
        cls, did = d.get("class"), d.get("id", "")
        gid = did
        for s in suffixes.values():
            if did.endswith(s):
                gid = did[: -len(s)]
        if cls not in ("Genrou", exc, spec.governor_class):
            kept.append(d)
            continue
        row = gen_by_id.get(gid)
        if row is None or not row["_online"]:
            continue  # offline or absent upstream -> drop the whole triplet
        b, mva = int(row["_bus"]), float(row["GenMVABase"])
        if cls == "Genrou":
            rec = dyr.get((b, "GENROU"))
            if rec is None:
                counts["missing_dyr"].append(gid)
            else:
                d["params"] = genrou_params(
                    rec,
                    mva,
                    float(row["GenMWSetPoint"]) / baseMVA,
                    float(row["GenMvrSetPoint"]) / baseMVA,
                )
                counts["genrou"] += 1
        elif cls == exc and (rec := dyr.get((b, "SEXS"))) is not None:
            d["params"] = sexs_params(rec)
            counts["sexs"] += 1
        elif cls == spec.governor_class and (rec := dyr.get((b, "TGOV1"))) is not None:
            d["params"] = tgov1_params(rec, mva)
            counts["tgov1"] += 1
        kept.append(d)
    case["devices"] = kept

    load = aux["Load"].copy()
    load["_bus"] = _num(load, "BusNum").astype(int)
    load["_rank"] = load.groupby("_bus").cumcount() + 1
    load_by_key = {(int(r["_bus"]), int(r["_rank"])): r for _, r in load.iterrows()}

    shunt = aux["Shunt"].copy()
    shunt["_bus"] = _num(shunt, "BusNum").astype(int)
    shunt_mvr = dict(zip(shunt["_bus"], _num(shunt, "SSNMVR")))

    for d in devices(case, "LoadZIP"):
        b = d["ports"]["bus"]
        p = d["params"]
        if is_shunt(d):
            if b in shunt_mvr and b in vmag:
                p["Qnom"] = -(shunt_mvr[b] / baseMVA) * vmag[b] ** 2
                p["Pnom"] = 0.0
                counts["shunt"] += 1
        else:
            rank = int(str(d["id"]).rsplit("_", 1)[1])
            row = load_by_key.get((b, rank))
            if row is not None:
                p["Pnom"] = float(row["LoadSMW"]) / baseMVA
                p["Qnom"] = float(row["LoadSMVR"]) / baseMVA
                counts["load"] += 1
        p.pop("Vnom", None)  # never a GridKit parameter

    if creator_overlay:
        apply_overlay(case, creator_overlay)
    return case, counts


# --------------------------------------------------------------------------
# build case.json from a (solved) MATPOWER .m  -- the aleatoric-track workhorse
# --------------------------------------------------------------------------


def read_m(path: str | Path):
    """Parse a MATPOWER `.m`. This package exports `CaseFrames`, not `Case`."""
    from matpowercaseframes import CaseFrames

    return CaseFrames(str(path))


def gen_ids(gen_df: pd.DataFrame, warn_multi_unit: bool = True) -> pd.Series:
    """`"{bus}_{rank}"` with rank counted over all rows at a bus, online or not.

    Rank is positional in `mpc.gen`. That is stable for a single machine per
    bus, but if an upstream tool (Sienna, a PCM export) reorders the rows of a
    multi-unit bus, ranks shift and dynamics would be attached to the wrong
    machine. Warn so it cannot happen unnoticed.
    """
    bus = gen_df["GEN_BUS"].astype(int)
    if warn_multi_unit:
        multi = bus.value_counts()
        multi = multi[multi > 1]
        if len(multi):
            warnings.warn(
                f"{len(multi)} bus(es) carry multiple generators "
                f"{multi.head().to_dict()}; ids are positional so the .m row order "
                f"must match the base case. Verify before trusting the mapping.",
                stacklevel=2,
            )
    return bus.astype(str) + "_" + (bus.groupby(bus).cumcount() + 1).astype(str)


def build_case_from_m(
    spec: CaseSpec,
    base_case: dict,
    m_path: str | Path,
    dyr: dict | None = None,
    creator_overlay: dict | None = None,
    patch_bus: bool = True,
    patch_gen: bool = True,
    patch_load: bool = True,
    remove_offline: bool = True,
    strict: bool = True,
) -> tuple[dict, dict]:
    """Build a `case.json` from a solved MATPOWER `.m`.

    This is the scenario path: perturb a `.m`, solve its power flow, then map
    the solution onto the base case's device list.

    Loads at multi-device buses are split in proportion to the base case's
    existing shares, so a bus total is distributed the same way the original
    case distributed it. Shunts are rebuilt from `mpc.bus` `BS`, never `QD`.

    Generators offline in the `.m` have their whole device triplet removed.
    Generators online in the `.m` but absent from the base case need a `dyr`
    to be added; without one this raises when `strict`.
    """
    case = copy.deepcopy(base_case)
    mp = read_m(m_path)
    baseMVA = float(mp.baseMVA)
    bus, gen = mp.bus.copy(), mp.gen.copy()
    counts = {
        "bus": 0,
        "gen": 0,
        "gen_removed": 0,
        "gen_added": 0,
        "load": 0,
        "shunt": 0,
        "gen_unmatched": [],
    }

    bus_i = bus["BUS_I"].astype(int)
    vm = dict(zip(bus_i, bus["VM"].astype(float)))
    va = dict(zip(bus_i, bus["VA"].astype(float)))
    pd_ = dict(zip(bus_i, bus["PD"].astype(float)))
    qd_ = dict(zip(bus_i, bus["QD"].astype(float)))
    bs_ = dict(zip(bus_i, bus["BS"].astype(float)))

    if patch_bus:
        for b in case.get("buses", []):
            n = b["number"]
            if n in vm:
                rad = math.radians(va[n])
                b["init"]["Vr"] = vm[n] * math.cos(rad)
                b["init"]["Vi"] = vm[n] * math.sin(rad)
                counts["bus"] += 1

    gen["_gid"] = gen_ids(gen)
    gen["_online"] = gen["GEN_STATUS"].astype(float) != 0
    grow = {r["_gid"]: r for _, r in gen.iterrows()}
    online_ids = {g for g, r in grow.items() if r["_online"]}

    json_genrou = {d["id"] for d in devices(case, "Genrou")}
    missing = sorted(online_ids - json_genrou)
    if missing and dyr is None:
        msg = (
            f"{len(missing)} generator(s) online in the .m but absent from the "
            f"base case and no .dyr supplied: {missing[:5]}"
        )
        if strict:
            raise ValueError(msg)
        counts["gen_unmatched"] = missing

    suffixes = spec.companion_suffixes
    if remove_offline:
        drop = json_genrou - online_ids
        if drop:
            keep = []
            for d in case["devices"]:
                did, cls = d.get("id", ""), d.get("class")
                gid = did
                for s in suffixes.values():
                    if did.endswith(s):
                        gid = did[: -len(s)]
                if (
                    cls in ("Genrou", spec.exciter_class, spec.governor_class)
                    and gid in drop
                ):
                    continue
                keep.append(d)
            counts["gen_removed"] = len(case["devices"]) - len(keep)
            case["devices"] = keep

    if patch_gen:
        for d in devices(case, "Genrou"):
            r = grow.get(d["id"])
            if r is None or not r["_online"]:
                continue
            d["params"]["p0"] = float(r["PG"]) / baseMVA
            d["params"]["q0"] = float(r["QG"]) / baseMVA
            d["params"]["mva"] = float(r["MBASE"])
            counts["gen"] += 1

    if dyr is not None and missing:
        counts["gen_added"] = _add_gens_from_dyr(
            case, spec, missing, grow, dyr, baseMVA
        )

    if patch_load and spec.load_class == "LoadZIP":
        lz = devices(case, "LoadZIP")
        by_bus: dict[int, list[dict]] = {}
        for d in lz:
            by_bus.setdefault(d["ports"]["bus"], []).append(d)
        for b, devs in by_bus.items():
            real = [d for d in devs if not is_shunt(d)]
            for d in devs:
                if is_shunt(d) and b in bs_ and b in vm:
                    d["params"]["Pnom"] = 0.0
                    d["params"]["Qnom"] = -(bs_[b] / baseMVA) * vm[b] ** 2
                    counts["shunt"] += 1
            if not real or b not in pd_:
                continue
            p_t, q_t = pd_[b] / baseMVA, qd_[b] / baseMVA
            sp = sum(float(d["params"].get("Pnom", 0.0)) for d in real)
            sq = sum(float(d["params"].get("Qnom", 0.0)) for d in real)
            for d in real:
                p = d["params"]
                p["Pnom"] = (
                    (float(p.get("Pnom", 0.0)) * p_t / sp) if sp else p_t / len(real)
                )
                p["Qnom"] = (
                    (float(p.get("Qnom", 0.0)) * q_t / sq) if sq else q_t / len(real)
                )
                p.pop("Vnom", None)  # not a GridKit parameter
                counts["load"] += 1

    if creator_overlay:
        apply_overlay(case, creator_overlay)
    return case, counts


def _add_gens_from_dyr(case, spec, missing, grow, dyr, baseMVA) -> int:
    """Append Genrou + exciter + governor triplets for newly online units."""
    sig = {s["name"]: s["signal_id"] for s in case.get("signals", [])}
    next_sig = (max(sig.values()) + 1) if sig else 0
    added = 0
    for gid in missing:
        r = grow[gid]
        b, mva = int(r["GEN_BUS"]), float(r["MBASE"])
        rec = dyr.get((b, "GENROU"))
        if rec is None:
            continue
        ports = {}
        for nm in ("speed", "pmech", "efd"):
            case.setdefault("signals", []).append(
                {"signal_id": next_sig, "name": f"{gid}_{nm}"}
            )
            ports[nm] = next_sig
            next_sig += 1
        case["devices"].append(
            {
                "class": "Genrou",
                "ports": {"bus": b, **ports},
                "id": gid,
                "params": genrou_params(
                    rec, mva, float(r["PG"]) / baseMVA, float(r["QG"]) / baseMVA
                ),
                "mon": ["delta", "omega"],
            }
        )
        if (g := dyr.get((b, "TGOV1"))) is not None:
            case["devices"].append(
                {
                    "class": spec.governor_class,
                    "ports": {"speed": ports["speed"], "pmech": ports["pmech"]},
                    "id": gid + spec.companion_suffixes[spec.governor_class],
                    "params": tgov1_params(g, mva),
                }
            )
        if (e := dyr.get((b, "SEXS"))) is not None:
            case["devices"].append(
                {
                    "class": spec.exciter_class,
                    "ports": {"bus": b, "efd": ports["efd"]},
                    "id": gid + spec.companion_suffixes[spec.exciter_class],
                    "params": sexs_params(e),
                }
            )
        added += 1
    return added


def write_case(case: dict, path: str | Path, indent: int = 4) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(case, indent=indent))
    return path
