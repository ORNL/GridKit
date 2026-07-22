"""
m_to_case.py
------------
Convert a MATPOWER .m power-flow solution file into a patched GridKit case.json.

Main entry point: patch_case_from_m()

What gets patched:
  - bus  init: Vr = Vm * cos(Va_rad), Vi = Vm * sin(Va_rad)
  - Genrou params: p0 = PG / baseMVA, q0 = QG / baseMVA  (online gens only)

What is NOT changed:
  - Machine dynamic parameters (H, D, Xd, ...)
  - Load parameters (ZIP load, constant power)
  - Solver configuration, monitor settings

Requires: matpowercaseframes>=2.1.0  (already in the gridkit conda env)
"""

import copy
import json
import math
import os

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _load_m(m_path):
    """Parse a MATPOWER .m file and return a MatpowerCaseData object."""
    try:
        from matpowercaseframes import Case
    except ImportError as e:
        raise ImportError(
            "matpowercaseframes is required: pip install matpowercaseframes"
        ) from e
    return Case(m_path)


def _bus_vm_va(case_data):
    """
    Return dict: bus_i (int) -> (Vm, Va_deg).

    Handles extended MATPOWER bus columns produced by PF solvers (columns
    beyond the standard 13 are ignored; only Vm and Va from cols 8-9 are used).
    """
    bus_df = case_data.bus
    vm_col = "VM" if "VM" in bus_df.columns else bus_df.columns[7]
    va_col = "VA" if "VA" in bus_df.columns else bus_df.columns[8]
    result = {}
    for row in bus_df.itertuples():
        bus_i = int(row.BUS_I)
        vm = float(getattr(row, vm_col))
        va = float(getattr(row, va_col))
        result[bus_i] = (vm, va)
    return result


def _gen_dispatch(case_data, baseMVA):
    """
    Return dict: json_gen_id (str, e.g. "2_1") -> (p0_pu, q0_pu).

    Only online generators (GEN_STATUS != 0) are included.
    The json_gen_id is derived by the same rank-within-bus rule used
    in gridkit_utils.attach_json_ids().

    p0, q0 are in per-unit on the system baseMVA base.
    MATPOWER PG, QG are in MW / MVAr; divide by baseMVA to get p.u.
    """
    gen_df = case_data.gen

    pg_col = "PG" if "PG" in gen_df.columns else gen_df.columns[1]
    qg_col = "QG" if "QG" in gen_df.columns else gen_df.columns[2]
    status_col = "GEN_STATUS" if "GEN_STATUS" in gen_df.columns else None

    rank = {}
    result = {}
    for row in gen_df.itertuples():
        bus = int(row.GEN_BUS)
        rank[bus] = rank.get(bus, 0) + 1
        status = int(getattr(row, status_col)) if status_col else 1
        if status == 0:
            continue
        gid = f"{bus}_{rank[bus]}"
        pg = float(getattr(row, pg_col))
        qg = float(getattr(row, qg_col))
        result[gid] = (pg / baseMVA, qg / baseMVA)
    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def patch_case_from_m(
    base_case_json_path,
    m_path,
    output_path=None,
    patch_bus_init=True,
    patch_gen_dispatch=True,
    verbose=True,
):
    """
    Produce a GridKit case.json with bus init conditions and generator dispatch
    taken from a MATPOWER PF solution .m file.

    Parameters
    ----------
    base_case_json_path : str
        Path to the base GridKit case JSON (e.g. hawaii.json).
    m_path : str
        Path to the MATPOWER .m PF solution file.
    output_path : str or None
        Where to write the patched JSON. If None, does not write.
    patch_bus_init : bool
        If True, update bus init Vr/Vi from Vm/Va in the .m solution.
    patch_gen_dispatch : bool
        If True, update Genrou p0/q0 from PG/QG in the .m solution.
    verbose : bool
        Print a summary of what was patched.

    Returns
    -------
    dict
        The patched case data (same structure as the GridKit case JSON).

    Notes
    -----
    MATPOWER Va is in degrees; GridKit expects Vr/Vi (rectangular) on a
    per-unit basis relative to v_base. Vm is already in per-unit.

    Generator p0/q0 in GridKit are in per-unit on the system MVA base
    (va_base field in the JSON header). MATPOWER PG/QG are in MW/MVAr;
    dividing by baseMVA converts to per-unit.

    Bus number correspondence: MATPOWER BUS_I == GridKit bus "number" field
    (confirmed by attach_json_ids() in gridkit_utils.py).

    Generator id correspondence: GridKit id "{bus}_{rank}" where rank is
    the 1-based position of the generator in the MATPOWER gen matrix among
    all generators at that bus (same convention as attach_json_ids()).
    """
    # --- load MATPOWER solution ---
    case_data = _load_m(m_path)
    baseMVA = float(case_data.baseMVA)

    bus_vm_va = _bus_vm_va(case_data) if patch_bus_init else {}
    gen_disp = _gen_dispatch(case_data, baseMVA) if patch_gen_dispatch else {}

    # --- load GridKit case ---
    with open(base_case_json_path) as fh:
        case = copy.deepcopy(json.load(fh))

    n_bus_patched = 0
    n_gen_patched = 0
    n_gen_missing = []

    # --- patch buses ---
    for bus in case.get("buses", []):
        bus_num = bus.get("number")
        if bus_num in bus_vm_va:
            vm, va_deg = bus_vm_va[bus_num]
            va_rad = math.radians(va_deg)
            bus["init"] = {
                "Vr": vm * math.cos(va_rad),
                "Vi": vm * math.sin(va_rad),
            }
            n_bus_patched += 1

    # --- patch Genrou dispatch ---
    for dev in case.get("devices", []):
        if dev.get("class") != "Genrou":
            continue
        gid = dev.get("id")
        if gid in gen_disp:
            p0, q0 = gen_disp[gid]
            dev.setdefault("params", {})
            dev["params"]["p0"] = float(p0)
            dev["params"]["q0"] = float(q0)
            n_gen_patched += 1
        else:
            n_gen_missing.append(gid)

    if verbose:
        print(f"patch_case_from_m: {os.path.basename(m_path)}")
        print(f"  baseMVA = {baseMVA}")
        if patch_bus_init:
            print(f"  bus init patched: {n_bus_patched}/{len(case.get('buses', []))}")
        if patch_gen_dispatch:
            print(f"  Genrou dispatch patched: {n_gen_patched}")
            if n_gen_missing:
                print(f"  Genrou ids not found in .m solution: {n_gen_missing}")

    if output_path is not None:
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        with open(output_path, "w") as fh:
            json.dump(case, fh, indent=2)
        if verbose:
            print(f"  written -> {output_path}")

    return case


def patch_cases_from_m_list(
    base_case_json_path,
    m_paths,
    output_dir,
    prefix="scenario",
    **kwargs,
):
    """
    Batch version of patch_case_from_m for a list of .m solution files.

    Each .m file produces one patched case JSON written to:
        output_dir/{prefix}_{i:04d}.json

    Parameters
    ----------
    base_case_json_path : str
    m_paths : list of str
        Ordered list of .m solution file paths (e.g. sorted hour_1 .. hour_8760).
    output_dir : str
    prefix : str
    **kwargs : passed to patch_case_from_m (patch_bus_init, patch_gen_dispatch, verbose)

    Returns
    -------
    list of str : paths to written scenario JSON files
    """
    os.makedirs(output_dir, exist_ok=True)
    written = []
    for i, m_path in enumerate(m_paths):
        out = os.path.join(output_dir, f"{prefix}_{i:04d}.json")
        patch_case_from_m(base_case_json_path, m_path, output_path=out, **kwargs)
        written.append(out)
    return written


def scenario_summary(scenario_json_paths, field="devices", device_class="Genrou"):
    """
    Load a list of scenario case JSON files and return a DataFrame summarising
    the dispatch (p0, q0) for each device_class across scenarios.

    Useful for a quick sanity check that the .m -> case.json mapping produced
    meaningful variation across scenarios.

    Returns
    -------
    pd.DataFrame with columns: scenario_idx, gen_id, p0, q0
    """
    rows = []
    for i, path in enumerate(scenario_json_paths):
        with open(path) as fh:
            case = json.load(fh)
        for dev in case.get("devices", []):
            if dev.get("class") != device_class:
                continue
            params = dev.get("params", {})
            rows.append(
                {
                    "scenario_idx": i,
                    "gen_id": dev.get("id"),
                    "p0": params.get("p0", float("nan")),
                    "q0": params.get("q0", float("nan")),
                }
            )
    return pd.DataFrame(rows)
