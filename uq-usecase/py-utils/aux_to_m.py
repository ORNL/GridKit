"""Write a MATPOWER `.m` case directly from a PowerWorld `.aux` file.

No template `.m` is read. The MATPOWER table layout is encoded here as a
documented constant, so the only inputs are the TAMU source files.

Why not reuse `case_ACTIVSg200.m`: it is an export of a *different* operating
point, it is missing branch (64,82), and its values are rounded to ~7
significant figures where the `.aux` carries ~15. Verified over the 245 shared
branches: max |dR| = 5.0e-07, max |dX| = 5.0e-07.

Mappings verified against ACTIVSg200 (see work-notes/m_to_case_helper.md):

    bus      Bus + aggregated Load + aggregated Shunt
    gen      Gen setpoints / limits / MVA base / status
    branch   lines use LineR/X/C, transformers use LineR:1/X:1/C:1;
             RATE_A = LineAMVA, TAP = LineTap:1, SHIFT = LinePhase
    gencost  c0 = GenFixedCost, c1 = GenIOB * GenFuelCost,
             c2 = GenIOC * GenFuelCost, c3 = GenIOD * GenFuelCost
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from m_to_case_utils import _num

# MATPOWER table layouts (MATPOWER manual, case format v2)
BUS_COLS = [
    "BUS_I",
    "BUS_TYPE",
    "PD",
    "QD",
    "GS",
    "BS",
    "BUS_AREA",
    "VM",
    "VA",
    "BASE_KV",
    "ZONE",
    "VMAX",
    "VMIN",
]
GEN_COLS = [
    "GEN_BUS",
    "PG",
    "QG",
    "QMAX",
    "QMIN",
    "VG",
    "MBASE",
    "GEN_STATUS",
    "PMAX",
    "PMIN",
    "PC1",
    "PC2",
    "QC1MIN",
    "QC1MAX",
    "QC2MIN",
    "QC2MAX",
    "RAMP_AGC",
    "RAMP_10",
    "RAMP_30",
    "RAMP_Q",
    "APF",
]
BRANCH_COLS = [
    "F_BUS",
    "T_BUS",
    "BR_R",
    "BR_X",
    "BR_B",
    "RATE_A",
    "RATE_B",
    "RATE_C",
    "TAP",
    "SHIFT",
    "BR_STATUS",
    "ANGMIN",
    "ANGMAX",
]
GENCOST_COLS = ["MODEL", "STARTUP", "SHUTDOWN", "NCOST", "c2", "c1", "c0"]

PQ, PV, SLACK = 1, 2, 3

# Not represented in a PowerWorld .aux; MATPOWER defaults are used and reported.
AUX_ABSENT = {
    "gen": [
        "PC1",
        "PC2",
        "QC1MIN",
        "QC1MAX",
        "QC2MIN",
        "QC2MAX",
        "RAMP_AGC",
        "RAMP_10",
        "RAMP_30",
        "RAMP_Q",
        "APF",
    ],
    "branch": ["RATE_B", "RATE_C", "ANGMIN", "ANGMAX"],
    "gencost": ["STARTUP", "SHUTDOWN"],
}


def _agg(df: pd.DataFrame, cols: dict[str, str]) -> pd.DataFrame:
    out = df.copy()
    out["_bus"] = _num(out, "BusNum").astype(int)
    for c in cols:
        out[c] = _num(out, c)
    return out.groupby("_bus")[list(cols)].sum().rename(columns=cols)


def aux_to_frames(
    aux: dict[str, pd.DataFrame],
    vmax: float = 1.1,
    vmin: float = 0.9,
) -> dict[str, pd.DataFrame]:
    """Build MATPOWER `bus`, `gen`, `branch`, `gencost` frames from a `.aux`.

    `vmax`/`vmin` are defaults: PowerWorld's per-bus voltage-limit fields are
    empty in the ACTIVSg exports.
    """
    ax_bus, ax_gen, ax_br = aux["Bus"], aux["Gen"], aux["Branch"]

    # --- bus ---
    bus_num = _num(ax_bus, "BusNum").astype(int)
    load = _agg(aux["Load"], {"LoadSMW": "pd", "LoadSMVR": "qd"})
    shunt = _agg(aux["Shunt"], {"SSNMW": "gs", "SSNMVR": "bs"})
    gen_buses = set(_num(ax_gen, "BusNum").astype(int))
    is_slack = ax_bus["BusSlack"].astype(str).str.strip().str.upper().eq("YES")

    bus = pd.DataFrame(
        {
            "BUS_I": bus_num.values,
            "BUS_TYPE": [
                SLACK if s else (PV if n in gen_buses else PQ)
                for n, s in zip(bus_num, is_slack)
            ],
            "PD": bus_num.map(load["pd"]).fillna(0.0).values,
            "QD": bus_num.map(load["qd"]).fillna(0.0).values,
            "GS": bus_num.map(shunt["gs"]).fillna(0.0).values,
            "BS": bus_num.map(shunt["bs"]).fillna(0.0).values,
            "BUS_AREA": _num(ax_bus, "AreaNum").fillna(1).astype(int).values,
            "VM": _num(ax_bus, "BusPUVolt").values,
            "VA": _num(ax_bus, "BusAngle").values,
            "BASE_KV": _num(ax_bus, "BusNomVolt").values,
            "ZONE": _num(ax_bus, "ZoneNum").fillna(1).astype(int).values,
            "VMAX": vmax,
            "VMIN": vmin,
        }
    )[BUS_COLS]

    # --- gen ---
    online = ax_gen["GenStatus"].astype(str).str.strip().str.lower().eq("closed")
    gen = pd.DataFrame(np.zeros((len(ax_gen), len(GEN_COLS))), columns=GEN_COLS)
    gen["GEN_BUS"] = _num(ax_gen, "BusNum").astype(int).values
    gen["PG"] = np.where(online, _num(ax_gen, "GenMWSetPoint"), 0.0)
    gen["QG"] = np.where(online, _num(ax_gen, "GenMvrSetPoint"), 0.0)
    gen["QMAX"] = _num(ax_gen, "GenMVRMax").values
    gen["QMIN"] = _num(ax_gen, "GenMVRMin").values
    gen["VG"] = _num(ax_gen, "GenVoltSet").values
    gen["MBASE"] = _num(ax_gen, "GenMVABase").values
    gen["GEN_STATUS"] = online.astype(int).values
    gen["PMAX"] = _num(ax_gen, "GenMWMax").values
    gen["PMIN"] = _num(ax_gen, "GenMWMin").values

    # --- branch: lines and transformers carry impedance in different columns ---
    is_xf = ax_br["BranchDeviceType"].astype(str).eq("Transformer").values

    def pick(line_col, xf_col):
        return np.where(is_xf, _num(ax_br, xf_col).values, _num(ax_br, line_col).values)

    tap = _num(ax_br, "LineTap:1").fillna(0.0).values
    branch = pd.DataFrame(np.zeros((len(ax_br), len(BRANCH_COLS))), columns=BRANCH_COLS)
    branch["F_BUS"] = _num(ax_br, "BusNum").astype(int).values
    branch["T_BUS"] = _num(ax_br, "BusNum:1").astype(int).values
    branch["BR_R"] = pick("LineR", "LineR:1")
    branch["BR_X"] = pick("LineX", "LineX:1")
    branch["BR_B"] = pick("LineC", "LineC:1")
    branch["RATE_A"] = _num(ax_br, "LineAMVA").fillna(0.0).values
    branch["TAP"] = np.where(is_xf, tap, 0.0)  # MATPOWER: 0 means unity for lines
    branch["SHIFT"] = _num(ax_br, "LinePhase").fillna(0.0).values
    branch["BR_STATUS"] = (
        ax_br["LineStatus"]
        .astype(str)
        .str.strip()
        .str.lower()
        .eq("closed")
        .astype(int)
        .values
    )

    # --- gencost: PowerWorld cubic IO curve scaled by fuel price ---
    fuel = _num(ax_gen, "GenFuelCost").fillna(1.0)
    c3 = (_num(ax_gen, "GenIOD").fillna(0.0) * fuel).values
    gencost = pd.DataFrame(
        {
            "MODEL": 2,
            "STARTUP": 0.0,
            "SHUTDOWN": 0.0,
            "NCOST": 4 if np.any(c3 != 0) else 3,
            "c2": (_num(ax_gen, "GenIOC").fillna(0.0) * fuel).values,
            "c1": (_num(ax_gen, "GenIOB").fillna(0.0) * fuel).values,
            "c0": _num(ax_gen, "GenFixedCost").fillna(0.0).values,
        }
    )[GENCOST_COLS]
    if np.any(c3 != 0):
        gencost.insert(4, "c3", c3)

    return {"bus": bus, "gen": gen, "branch": branch, "gencost": gencost}


def _parse_cell_array(text: str, name: str) -> list[str] | None:
    """Read `mpc.<name> = { 'a'; 'b'; };`. matpowercaseframes skips these."""
    import re

    m = re.search(rf"mpc\.{name}\s*=\s*\{{(.*?)\}}\s*;", text, re.S)
    if not m:
        return None
    return re.findall(r"'([^']*)'", m.group(1))


def m_to_frames(m_path: str | Path) -> dict[str, pd.DataFrame]:
    """Read any MATPOWER `.m` into the same frame dict `aux_to_frames` returns.

    This is what makes the perturbation layer source-agnostic: frames can come
    from a PowerWorld `.aux` (bootstrap) or from a `.m` produced by Sienna, a
    PCM run, or an earlier stage of this pipeline.
    """
    from matpowercaseframes import CaseFrames

    m_path = Path(m_path)
    cf = CaseFrames(str(m_path))
    out = {"bus": cf.bus.copy(), "gen": cf.gen.copy(), "branch": cf.branch.copy()}
    if hasattr(cf, "gencost"):
        out["gencost"] = cf.gencost.copy()
    if hasattr(cf, "bus_name"):
        out["bus_name"] = list(cf.bus_name)
    text = m_path.read_text()
    for name in ("genfuel", "gentype"):
        vals = _parse_cell_array(text, name)
        if vals is not None:
            out[name] = vals
    out["baseMVA"] = float(cf.baseMVA)
    return out


def _block(name: str, df: pd.DataFrame, int_cols: dict[str, str]) -> str:
    # Column-wise, because iterrows() would upcast integer columns to float.
    cols = list(df.columns)
    data = {c: df[c].to_numpy() for c in cols}
    lines = [f"mpc.{name} = ["]
    for i in range(len(df)):
        cells = [
            f"{int(data[c][i]):d}" if c in int_cols else f"{data[c][i]:.10g}"
            for c in cols
        ]
        lines.append("\t" + "\t".join(cells) + ";")
    lines.append("];")
    return "\n".join(lines)


def _cell_block(name: str, values) -> str:
    rows = "\n".join(f"\t'{v}';" for v in values)
    return f"mpc.{name} = {{\n{rows}\n}};"


def write_matpower_m(
    source: dict,
    out_m: str | Path,
    case_name: str = "case_pw",
    baseMVA: float = 100.0,
    include_gencost: bool = True,
    **frame_kwargs,
) -> dict:
    """Write a MATPOWER case.

    `source` is either a `.aux` table dict (as returned by `parse_aux`, detected
    by the presence of a `Bus` key) or an already-built frame dict from
    `aux_to_frames` / `m_to_frames` / a perturbation step.
    """
    out_m = Path(out_m)
    if "Bus" in source:  # raw .aux tables
        aux = source
        f = aux_to_frames(aux, **frame_kwargs)
        names = {
            "bus_name": aux["Bus"].get("BusName"),
            "genfuel": (
                aux["Gen"]["GenFuelType"].str.lower()
                if "GenFuelType" in aux["Gen"]
                else None
            ),
        }
    else:
        f = source
        baseMVA = float(source.get("baseMVA", baseMVA))
        names = {"bus_name": source.get("bus_name"), "genfuel": source.get("genfuel")}
    ints = {
        "BUS_I": "d",
        "BUS_TYPE": "d",
        "BUS_AREA": "d",
        "ZONE": "d",
        "GEN_BUS": "d",
        "GEN_STATUS": "d",
        "F_BUS": "d",
        "T_BUS": "d",
        "BR_STATUS": "d",
        "MODEL": "d",
        "NCOST": "d",
    }
    tab = chr(9)

    parts = [
        f"function mpc = {case_name}",
        f"%{case_name.upper()}  Built from a PowerWorld .aux by write_matpower_m().",
        "%   No template .m was used; see py-utils/aux_to_m.py for field mappings.",
        "",
        "%% MATPOWER Case Format : Version 2",
        "mpc.version = '2';",
        "",
        "%%-----  Power Flow Data  -----%%",
        f"mpc.baseMVA = {baseMVA:g};",
        "",
        "%% bus data",
        f"%{tab}{tab.join(BUS_COLS)}",
        _block(
            "bus",
            f["bus"].astype(
                {"BUS_I": int, "BUS_TYPE": int, "BUS_AREA": int, "ZONE": int}
            ),
            ints,
        ),
        "",
        "%% generator data",
        f"%{tab}{tab.join(GEN_COLS)}",
        _block("gen", f["gen"].astype({"GEN_BUS": int, "GEN_STATUS": int}), ints),
        "",
        "%% branch data",
        f"%{tab}{tab.join(BRANCH_COLS)}",
        _block(
            "branch",
            f["branch"].astype({"F_BUS": int, "T_BUS": int, "BR_STATUS": int}),
            ints,
        ),
    ]
    if include_gencost and "gencost" in f:
        parts += [
            "",
            "%%-----  OPF Data  -----%%",
            "%% generator cost data",
            f"%{tab}{tab.join(f['gencost'].columns)}",
            _block("gencost", f["gencost"].astype({"MODEL": int, "NCOST": int}), ints),
        ]
    if names.get("bus_name") is not None:
        parts += ["", "%% bus names", _cell_block("bus_name", names["bus_name"])]
    if names.get("genfuel") is not None:
        parts += [
            "",
            "%% generator fuel type",
            _cell_block("genfuel", names["genfuel"]),
        ]

    out_m.parent.mkdir(parents=True, exist_ok=True)
    out_m.write_text("\n".join(parts) + "\n")

    on = f["gen"]["GEN_STATUS"] == 1
    return {
        "path": out_m,
        "buses": len(f["bus"]),
        "gens": len(f["gen"]),
        "gens_online": int(on.sum()),
        "branches": len(f["branch"]),
        "transformers": int((f["branch"]["TAP"] != 0).sum()),
        "slack_bus": f["bus"].loc[f["bus"]["BUS_TYPE"] == SLACK, "BUS_I"].tolist(),
        "total_PD_MW": float(f["bus"]["PD"].sum()),
        "total_PG_MW": float(f["gen"].loc[on, "PG"].sum()),
        "defaulted_fields": AUX_ABSENT,
    }
