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


# ---------------------------------------------------------------------------
# build_case_from_solved_m: superset of patch_case_from_m
# ---------------------------------------------------------------------------
#
# Extends patch_case_from_m with:
#   - LoadZIP patching (Pnom, Qnom, Vnom) with multi-ZIP proportional scaling
#     and shunt/capacitor exclusion.
#   - Offline-generator device-triplet removal ({Genrou, Tgov1, SexsPti} or
#     {Genrou, Tgov1, Ieeet1}) when the .m marks a generator GEN_STATUS=0.
#
# NOT handled (documented limitation, deferred to production):
#   - Introducing a NEW Genrou (and companion devices) for a generator that
#     is online in the .m but absent from the base JSON. This requires a
#     full device-list rebuild, out of scope here. The function will raise
#     if this case is encountered so the caller notices immediately.
#
# See work-notes/m_to_case_helper.md for the full spec and worked examples.

# Genrou companion device classes that share the same {bus}_{rank} id prefix.
# Illinois uses SexsPti; Hawaii uses Ieeet1. Tgov1 is common.
_GENROU_COMPANION_CLASSES = ("Tgov1", "SexsPti", "Ieeet1", "Gast", "Hygov")


def _gen_status_by_id(case_data):
    """
    Return dict: json_gen_id (str, e.g. "2_1") -> GEN_STATUS (int, 0 or 1).

    Uses the same rank-within-bus convention as _gen_dispatch(): all rows in
    mpc.gen (online or offline) are ranked by their order at each bus.
    """
    gen_df = case_data.gen
    status_col = "GEN_STATUS" if "GEN_STATUS" in gen_df.columns else None
    rank = {}
    result = {}
    for row in gen_df.itertuples():
        bus = int(row.GEN_BUS)
        rank[bus] = rank.get(bus, 0) + 1
        status = int(getattr(row, status_col)) if status_col else 1
        gid = f"{bus}_{rank[bus]}"
        result[gid] = status
    return result


def _bus_pd_qd(case_data, baseMVA):
    """
    Return dict: bus_i (int) -> (Pnom_pu, Qnom_pu, Vnom_pu).

    Pnom_pu = PD / baseMVA, Qnom_pu = QD / baseMVA, Vnom_pu = VM.
    Vnom is included so callers can always set the ZIP reference voltage
    even when PD/QD are zero.
    """
    bus_df = case_data.bus
    pd_col = "PD" if "PD" in bus_df.columns else bus_df.columns[2]
    qd_col = "QD" if "QD" in bus_df.columns else bus_df.columns[3]
    vm_col = "VM" if "VM" in bus_df.columns else bus_df.columns[7]
    result = {}
    for row in bus_df.itertuples():
        bus_i = int(row.BUS_I)
        pd_mw = float(getattr(row, pd_col))
        qd_mvar = float(getattr(row, qd_col))
        vm = float(getattr(row, vm_col))
        result[bus_i] = (pd_mw / baseMVA, qd_mvar / baseMVA, vm)
    return result


def _is_shunt_loadzip(dev):
    """
    Recognise a reactive-compensation (capacitor / shunt) LoadZIP.

    Convention (illinois.md, uq_plan.md): a LoadZIP with Pnom == 0 AND
    Qnom < 0 represents reactive compensation and must NOT be updated with
    PD/QD from the .m file. The id may also start with "shunt_" in some
    cases; both patterns are accepted.
    """
    params = dev.get("params", {})
    pnom = float(params.get("Pnom", 0.0))
    qnom = float(params.get("Qnom", 0.0))
    dev_id = str(dev.get("id", ""))
    if dev_id.startswith("shunt_"):
        return True
    if pnom == 0.0 and qnom < 0.0:
        return True
    return False


def build_case_from_solved_m(
    base_case_json_path,
    m_path,
    output_path=None,
    patch_bus_init=True,
    patch_gen_dispatch=True,
    patch_load=True,
    remove_offline_gens=True,
    strict=True,
    verbose=True,
):
    """
    Build a GridKit case.json from a PM.jl-solved MATPOWER .m file, updating
    bus init conditions, generator dispatch, LoadZIP demand, and pruning
    offline-generator device triplets.

    This is a superset of patch_case_from_m(); use build_case_from_solved_m
    when the .m file has been re-solved after modifying load, wind, or
    generator status. Use patch_case_from_m for the simpler p0/q0 + Vr/Vi
    patch when only PG/QG/Vm/Va have changed.

    Parameters
    ----------
    base_case_json_path : str
        Path to the base GridKit case JSON (e.g. illinois.json).
    m_path : str
        Path to the MATPOWER .m solved PF file.
    output_path : str or None
        Where to write the resulting JSON. If None, does not write.
    patch_bus_init : bool
        Update bus init Vr/Vi from Vm/Va.
    patch_gen_dispatch : bool
        Update Genrou p0/q0 from PG/QG (online gens only).
    patch_load : bool
        Update LoadZIP Pnom/Qnom/Vnom from PD/QD/Vm.
    remove_offline_gens : bool
        If True, remove Genrou + companion (Tgov1, SexsPti/Ieeet1) devices
        for gens whose GEN_STATUS=0 in the .m.
    strict : bool
        If True, raise on structural mismatches: (a) an online gen in the .m
        that has no matching Genrou in the JSON (device rebuild required,
        deferred). If False, log a warning and continue.
    verbose : bool
        Print a per-category summary of what was modified.

    Returns
    -------
    (case_dict, counts_dict) : tuple
        case_dict: the resulting case data.
        counts_dict: {"bus_updated", "gen_updated", "gen_removed",
                     "load_updated", "shunt_skipped", "gen_unmatched"}

    Notes
    -----
    LoadZIP proportional scaling for multi-ZIP buses:
        For a bus with N LoadZIP devices (excluding shunts), let
        S_p = sum of existing Pnom across the N devices,
        S_q = sum of existing Qnom.
        For each device i, new Pnom_i = Pnom_i * (PD_pu / S_p),
                          new Qnom_i = Qnom_i * (QD_pu / S_q).
        If S_p == 0 and PD_pu > 0, distribute PD_pu equally across the N
        devices (documented degenerate case).
        Vnom is set to Vm at that bus for every device (all devices share
        the same bus voltage).
    """
    case_data = _load_m(m_path)
    baseMVA = float(case_data.baseMVA)

    bus_vm_va = _bus_vm_va(case_data) if patch_bus_init else {}
    gen_disp = _gen_dispatch(case_data, baseMVA) if patch_gen_dispatch else {}
    gen_status = _gen_status_by_id(case_data)
    bus_load = _bus_pd_qd(case_data, baseMVA) if patch_load else {}

    with open(base_case_json_path) as fh:
        case = copy.deepcopy(json.load(fh))

    counts = {
        "bus_updated": 0,
        "gen_updated": 0,
        "gen_removed": 0,
        "load_updated": 0,
        "shunt_skipped": 0,
        "gen_unmatched": 0,
    }
    unmatched_gen_ids = []
    warnings = []

    # ----- 1. bus init -----
    if patch_bus_init:
        for bus in case.get("buses", []):
            bus_num = bus.get("number")
            if bus_num in bus_vm_va:
                vm, va_deg = bus_vm_va[bus_num]
                va_rad = math.radians(va_deg)
                bus["init"] = {
                    "Vr": vm * math.cos(va_rad),
                    "Vi": vm * math.sin(va_rad),
                }
                counts["bus_updated"] += 1

    # ----- 2. Genrou dispatch update + offline flagging -----
    genrou_to_remove = set()  # set of Genrou ids to be removed at end
    for dev in case.get("devices", []):
        if dev.get("class") != "Genrou":
            continue
        gid = dev.get("id")
        status = gen_status.get(gid)
        if status is None:
            # Genrou in JSON but no corresponding row in .m: not solvable here.
            counts["gen_unmatched"] += 1
            unmatched_gen_ids.append(gid)
            continue
        if status == 0 and remove_offline_gens:
            genrou_to_remove.add(gid)
            continue
        # online: update dispatch
        if patch_gen_dispatch and gid in gen_disp:
            p0, q0 = gen_disp[gid]
            dev.setdefault("params", {})
            dev["params"]["p0"] = float(p0)
            dev["params"]["q0"] = float(q0)
            counts["gen_updated"] += 1

    # ----- 2b. detect online .m gens missing from JSON (rebuild required) -----
    json_genrou_ids = {
        dev.get("id") for dev in case.get("devices", []) if dev.get("class") == "Genrou"
    }
    for gid, status in gen_status.items():
        if status == 1 and gid not in json_genrou_ids:
            msg = (
                f"generator {gid} is online in .m but has no Genrou in base JSON; "
                "full device-list rebuild is required (out of scope for "
                "build_case_from_solved_m). Constrain gen-off perturbations to "
                "gens already in the base JSON, or extend this function."
            )
            if strict:
                raise ValueError(msg)
            warnings.append(msg)

    # ----- 3. remove offline-gen device triplets -----
    if remove_offline_gens and genrou_to_remove:
        new_devices = []
        for dev in case.get("devices", []):
            cls = dev.get("class")
            dev_id = str(dev.get("id", ""))
            # Genrou itself: drop if in the removal set.
            if cls == "Genrou" and dev_id in genrou_to_remove:
                counts["gen_removed"] += 1
                continue
            # Companion devices: drop if the id matches "<genrou_id>" or
            # "<genrou_id>_..." for a genrou in the removal set.
            if cls in _GENROU_COMPANION_CLASSES:
                match = False
                for gid in genrou_to_remove:
                    if dev_id == gid or dev_id.startswith(gid + "_"):
                        match = True
                        break
                if match:
                    counts["gen_removed"] += 1
                    continue
            new_devices.append(dev)
        case["devices"] = new_devices

    # ----- 4. LoadZIP patching -----
    if patch_load:
        # group LoadZIP devices by bus
        loadzip_by_bus = {}
        for dev in case.get("devices", []):
            if dev.get("class") != "LoadZIP":
                continue
            bus_num = dev.get("ports", {}).get("bus")
            loadzip_by_bus.setdefault(bus_num, []).append(dev)

        for bus_num, devs in loadzip_by_bus.items():
            if bus_num not in bus_load:
                continue
            pd_pu, qd_pu, vm = bus_load[bus_num]

            # split shunts from real loads
            real_loads = []
            for d in devs:
                if _is_shunt_loadzip(d):
                    counts["shunt_skipped"] += 1
                    # ensure Vnom is still refreshed for capacitor/shunt so its
                    # per-unit reference tracks the new bus voltage
                    d.setdefault("params", {})["Vnom"] = vm
                    continue
                real_loads.append(d)

            if not real_loads:
                continue

            # zero-load bus: set all real loads to Pnom=Qnom=0, Vnom=Vm
            if pd_pu == 0.0 and qd_pu == 0.0:
                for d in real_loads:
                    p = d.setdefault("params", {})
                    p["Pnom"] = 0.0
                    p["Qnom"] = 0.0
                    p["Vnom"] = vm
                    counts["load_updated"] += 1
                continue

            # sum existing Pnom, Qnom for proportional scaling
            s_p = sum(float(d.get("params", {}).get("Pnom", 0.0)) for d in real_loads)
            s_q = sum(float(d.get("params", {}).get("Qnom", 0.0)) for d in real_loads)

            # degenerate case: existing Pnom sum is zero but .m has PD > 0.
            # Distribute PD/QD equally across the N devices.
            if s_p == 0.0 and pd_pu != 0.0:
                warnings.append(
                    f"bus {bus_num}: base JSON sum(Pnom)==0 but .m has PD>0; "
                    "distributing PD equally across LoadZIP devices"
                )
                n = len(real_loads)
                for d in real_loads:
                    p = d.setdefault("params", {})
                    p["Pnom"] = pd_pu / n
                    p["Qnom"] = qd_pu / n if s_q == 0.0 else float(p.get("Qnom", 0.0))
                    p["Vnom"] = vm
                    counts["load_updated"] += 1
                continue

            # standard proportional scaling
            for d in real_loads:
                p = d.setdefault("params", {})
                pnom_old = float(p.get("Pnom", 0.0))
                qnom_old = float(p.get("Qnom", 0.0))
                p["Pnom"] = pnom_old * (pd_pu / s_p) if s_p != 0.0 else 0.0
                p["Qnom"] = qnom_old * (qd_pu / s_q) if s_q != 0.0 else 0.0
                p["Vnom"] = vm
                counts["load_updated"] += 1

    # ----- report -----
    if verbose:
        print(f"build_case_from_solved_m: {os.path.basename(m_path)}")
        print(f"  baseMVA = {baseMVA}")
        n_json_buses = len(case.get("buses", []))
        print(f"  buses updated: {counts['bus_updated']}/{n_json_buses}")
        print(f"  Genrous updated: {counts['gen_updated']}")
        if remove_offline_gens:
            print(f"  device triplets removed (offline gens): {counts['gen_removed']}")
        print(f"  LoadZIP updated: {counts['load_updated']}")
        print(f"  shunt/capacitor LoadZIP skipped: {counts['shunt_skipped']}")
        if counts["gen_unmatched"]:
            print(
                f"  Genrou in JSON not found in .m: "
                f"{counts['gen_unmatched']} ({unmatched_gen_ids})"
            )
        for w in warnings:
            print(f"  WARNING: {w}")

    if output_path is not None:
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        with open(output_path, "w") as fh:
            json.dump(case, fh, indent=2)
        if verbose:
            print(f"  written -> {output_path}")

    return case, counts


def build_cases_from_solved_m_list(
    base_case_json_path,
    m_paths,
    output_dir,
    prefix="scenario",
    **kwargs,
):
    """
    Batch version of build_case_from_solved_m for a list of solved .m files.

    Each .m file produces one case JSON written to:
        output_dir/{prefix}_{i:04d}.json

    Parameters
    ----------
    base_case_json_path : str
    m_paths : list of str
    output_dir : str
    prefix : str
    **kwargs : passed to build_case_from_solved_m.

    Returns
    -------
    (paths, counts_list) : tuple
        paths: list of output JSON paths.
        counts_list: list of counts dicts, one per input .m.
    """
    os.makedirs(output_dir, exist_ok=True)
    written = []
    counts_list = []
    for i, m_path in enumerate(m_paths):
        out = os.path.join(output_dir, f"{prefix}_{i:04d}.json")
        _, counts = build_case_from_solved_m(
            base_case_json_path, m_path, output_path=out, **kwargs
        )
        written.append(out)
        counts_list.append(counts)
    return written, counts_list
