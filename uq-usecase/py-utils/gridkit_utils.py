"""
GridKit helper utilities for use in gridkit_helper.ipynb and related notebooks.
"""

import copy
import json
import os
import shutil
import subprocess
from numbers import Number

import numpy as np
import pandas as pd
from scipy.stats import qmc, norm as sp_norm

# ---------------------------------------------------------------------------
# Monitorable variables by element class
# Source: GridKit docs INPUT_FORMAT.md
#   - Bus classes table   (init vars + other variables available to monitor)
#   - Device classes table (Variables available to monitor column)
# ---------------------------------------------------------------------------
MONITORABLE_VARS_BY_ELEMENT = {
    # --- Bus classes ---
    "Bus": ["Vr", "Vi", "Vm", "Va"],
    # --- Device classes ---
    "Branch": ["ir1", "ii1", "im1", "p1", "q1", "ir2", "ii2", "im2", "p2", "q2"],
    "Load": ["p", "q"],
    "Genrou": ["ir", "ii", "p", "q", "delta", "omega", "speed"],
    "Gensal": [
        "ir",
        "ii",
        "p",
        "q",
        "delta",
        "omega",
        "speed",
        "Eqp",
        "psidp",
        "psiqpp",
        "psidpp",
        "vd",
        "vq",
        "te",
        "id",
        "iq",
    ],
    "GenClassical": ["ir", "ii", "p", "q", "delta", "omega"],
    "Tgov1": [],
    "Ieeet1": ["efd", "ksat"],
    "SexsPti": ["efd"],
    "Ieeest": ["vss"],
    "BusFault": ["state", "ir", "ii"],
}


# ---------------------------------------------------------------------------
# UQ workflow helpers
# ---------------------------------------------------------------------------


def generate_samples(param_specs, N, seed=None, method="lhs"):
    """
    Sample parameters according to per-spec distributions.

    Parameters
    ----------
    param_specs : list of dict, each with keys:
        "class"  - device class string (e.g. "Genrou")
        "id"     - device id string    (e.g. "genrou_2_1")
        "param"  - parameter name      (e.g. "H")
        "dist"   - one of "uniform" (default) or "normal"

        For dist="uniform":
            Option A — percent:  "nominal" + "pct"
                lo = nominal * (1 - pct),  hi = nominal * (1 + pct)
            Option B — fixed:    "lo" + "hi"

        For dist="normal":
            "mean" and "std"

    N      : int, number of samples
    seed   : int or None, for reproducibility
    method : "lhs"    — Latin Hypercube (better space-filling, recommended)
             "random" — independent random draws per parameter

    Returns
    -------
    pd.DataFrame shape (N, len(param_specs))
        columns named "{class}_{id}_{param}", e.g. "Genrou_genrou_2_1_H"

    Notes
    -----
    LHS works by drawing N uniform samples stratified across [0,1] per
    dimension, then mapping each through the inverse CDF (ppf) of the
    requested distribution. This gives correct marginals for both uniform
    and normal distributions while maximising coverage of the input space.
    """
    cols = [f"{s['class']}_{s['id']}_{s['param']}" for s in param_specs]
    d = len(param_specs)

    if method == "lhs":
        sampler = qmc.LatinHypercube(d=d, seed=seed)
        unit = sampler.random(N)  # N × d uniform in (0, 1)
    else:
        rng = np.random.default_rng(seed)
        unit = rng.uniform(0, 1, size=(N, d))

    data = {}
    for j, spec in enumerate(param_specs):
        dist = spec.get("dist", "uniform")
        u = unit[:, j]

        if dist == "uniform":
            if "lo" in spec and "hi" in spec:
                lo, hi = spec["lo"], spec["hi"]
            else:
                lo = spec["nominal"] * (1 - spec["pct"])
                hi = spec["nominal"] * (1 + spec["pct"])
            # uniform ppf: lo + u*(hi-lo)
            data[cols[j]] = lo + u * (hi - lo)

        elif dist == "normal":
            mean, std = spec["mean"], spec["std"]
            # normal ppf via scipy
            data[cols[j]] = sp_norm.ppf(u, loc=mean, scale=std)

        else:
            raise ValueError(f"Unknown dist '{dist}' in param_spec for {cols[j]}")

    return pd.DataFrame(data)


def make_run_dir(
    base_case_dir,
    run_root,
    i,
    sample_row,
    param_specs,
    monitors_by_class,
    case_fn=None,
    solver_fn=None,
    solver_overrides=None,
):
    """
    Create run_root/run_{i:03d}/, copy case+solver JSON, patch params and monitors.

    Parameters
    ----------
    base_case_dir     : str, directory containing the base case files
    run_root          : str, parent directory for all runs
    i                 : int, run index
    sample_row        : pd.Series, one row from generate_samples() DataFrame
    param_specs       : list of dict (same structure as generate_samples)
    monitors_by_class : dict mapping device/bus class (lowercase ok) -> list of var names
                        e.g. {"bus": ["Vm","Va"], "infinite_bus": ["Vm","Va"],
                               "genrou": ["delta","omega"]}
                        Classes not listed get their mon arrays removed.
    case_fn           : str or None; if None, auto-detected (single .case.json in dir)
    solver_fn         : str or None; if None, auto-detected (single .solver.json in dir)
    solver_overrides  : dict or None; key-value pairs to patch into the solver JSON.
                        Supports top-level keys (e.g. "tmax") and "events" as a list.
                        For "events", each entry must have "index" (0-based position in
                        the events array) plus whichever fields to overwrite, e.g.:
                          {"tmax": 20.0,
                           "events": [
                               {"index": 0, "time": 2.0},
                               {"index": 1, "time": 2.1},
                           ]}

    Returns
    -------
    str : path to the new run directory
    """
    if case_fn is None:
        candidates = [f for f in os.listdir(base_case_dir) if f.endswith(".case.json")]
        if len(candidates) == 0:
            # fall back: any .json that is not a .solver.json (e.g. hawaii.json)
            candidates = [
                f
                for f in os.listdir(base_case_dir)
                if f.endswith(".json") and not f.endswith(".solver.json")
            ]
        if len(candidates) != 1:
            raise ValueError(
                f"Expected 1 case .json in {base_case_dir}, found: {candidates}"
            )
        case_fn = candidates[0]
    if solver_fn is None:
        candidates = [
            f for f in os.listdir(base_case_dir) if f.endswith(".solver.json")
        ]
        if len(candidates) != 1:
            raise ValueError(
                f"Expected 1 .solver.json in {base_case_dir}, found: {candidates}"
            )
        solver_fn = candidates[0]

    run_dir = os.path.join(run_root, f"run_{i:03d}")
    os.makedirs(run_dir, exist_ok=True)

    # load solver, strip testing keys; remove output_file so the simulator
    # uses its default (always writes mon.csv; output_file only controls a symlink alias)
    with open(os.path.join(base_case_dir, solver_fn)) as f:
        solver = json.load(f)
    solver.pop("reference_file", None)
    solver.pop("error_tolerance", None)
    solver.pop("output_file", None)

    # apply solver_overrides (tmax, events, dt, ...)
    if solver_overrides:
        for key, val in solver_overrides.items():
            if key == "events":
                # val is a list of {index, ...fields}; patch each event in-place
                for patch in val:
                    idx = patch["index"]
                    for k, v in patch.items():
                        if k != "index":
                            solver["events"][idx][k] = v
            else:
                solver[key] = val

    with open(os.path.join(run_dir, solver_fn), "w") as f:
        json.dump(solver, f, indent=2)

    # load and deep-copy case
    with open(os.path.join(base_case_dir, case_fn)) as f:
        case = copy.deepcopy(json.load(f))

    # --- patch device params ---
    spec_lookup = {(s["class"].lower(), s["id"]): s for s in param_specs}
    col_lookup = {
        (s["class"].lower(), s["id"]): f"{s['class']}_{s['id']}_{s['param']}"
        for s in param_specs
    }
    for dev in case.get("devices", []):
        key = (dev.get("class", "").lower(), dev.get("id", ""))
        if key in spec_lookup:
            spec = spec_lookup[key]
            dev.setdefault("params", {})[spec["param"]] = float(
                sample_row[col_lookup[key]]
            )

    # --- patch monitors ---
    mon_lower = {k.lower(): v for k, v in monitors_by_class.items()}
    for bus in case.get("buses", []):
        cls = bus.get("class", "").lower()
        mon_vars = mon_lower.get(cls, [])
        if mon_vars:
            bus["mon"] = mon_vars
        elif "mon" in bus:
            del bus["mon"]
    for dev in case.get("devices", []):
        cls = dev.get("class", "").lower()
        mon_vars = mon_lower.get(cls, [])
        if mon_vars:
            dev["mon"] = mon_vars
        elif "mon" in dev:
            del dev["mon"]

    with open(os.path.join(run_dir, case_fn), "w") as f:
        json.dump(case, f, indent=2)

    return run_dir


def run_sample(run_dir, runner, solver_fn=None, timeout=300):
    """
    Run DynamicSimulation in run_dir.

    Parameters
    ----------
    run_dir   : str
    runner    : str, full path to DynamicSimulation executable
    solver_fn : str or None; if None auto-detected
    timeout   : int, seconds

    Returns
    -------
    subprocess.CompletedProcess
    """
    if solver_fn is None:
        candidates = [f for f in os.listdir(run_dir) if f.endswith(".solver.json")]
        if len(candidates) != 1:
            raise ValueError(
                f"Expected 1 .solver.json in {run_dir}, found: {candidates}"
            )
        solver_fn = candidates[0]
    return subprocess.run(
        [runner, solver_fn],
        cwd=run_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def collect_and_save(run_root, samples_df, out_path, mon_fn="mon.csv", mode="stacked"):
    """
    Read all run_i/mon.csv and save as Parquet in one of two formats.

    Parameters
    ----------
    run_root   : str, parent dir containing run_000/, run_001/, ...
    samples_df : pd.DataFrame from generate_samples(); index = run index 0..N-1
    out_path   : str
        mode="stacked" → path to a single output .parquet file
        mode="per_run" → path to an output *directory*; writes run_NNN.parquet per run
    mon_fn     : str, monitor output filename (default "mon.csv")
    mode       : "stacked" (default) or "per_run"
        "stacked"  — single flat file with run_id + time + signals + param cols
        "per_run"  — one file per run, wide format (rows=timesteps, cols=signals only);
                     param values are NOT duplicated — they live in samples.csv

    Returns
    -------
    mode="stacked" → pd.DataFrame : run_id | time | <mon vars> | <param cols>
    mode="per_run" → list of str  : paths to written run_NNN.parquet files
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    if mode == "stacked":
        frames = []
        for i, row in samples_df.iterrows():
            mon_path = os.path.join(run_root, f"run_{i:03d}", mon_fn)
            if not os.path.exists(mon_path) or os.path.islink(mon_path):
                print(f"  WARNING: missing {mon_path}, skipping run {i}")
                continue
            df = pd.read_csv(mon_path)
            df = df.rename(columns={df.columns[0]: "time"})
            df.insert(0, "run_id", i)
            for col, val in row.items():
                df[col] = val
            frames.append(df)

        if not frames:
            raise RuntimeError("No mon.csv files found — did the runs complete?")

        combined = pd.concat(frames, ignore_index=True)
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        pq.write_table(pa.Table.from_pandas(combined, preserve_index=False), out_path)
        print(f"Saved {len(combined)} rows ({len(frames)} runs) -> {out_path}")
        return combined

    elif mode == "per_run":
        os.makedirs(out_path, exist_ok=True)
        written = []
        for i, _row in samples_df.iterrows():
            mon_path = os.path.join(run_root, f"run_{i:03d}", mon_fn)
            if not os.path.exists(mon_path) or os.path.islink(mon_path):
                print(f"  WARNING: missing {mon_path}, skipping run {i}")
                continue
            df = pd.read_csv(mon_path)
            df = df.rename(columns={df.columns[0]: "time"})
            out_file = os.path.join(out_path, f"run_{i:03d}.parquet")
            pq.write_table(pa.Table.from_pandas(df, preserve_index=False), out_file)
            written.append(out_file)

        if not written:
            raise RuntimeError("No mon.csv files found — did the runs complete?")

        print(f"Saved {len(written)} run files -> {out_path}/run_NNN.parquet")
        return written

    else:
        raise ValueError(f"Unknown mode '{mode}': must be 'stacked' or 'per_run'")


# ---------------------------------------------------------------------------
# Dispatch editing helpers (aleatoric UQ prep)
# ---------------------------------------------------------------------------


def read_genrou_dispatch(case_json_path):
    """
    Read all Genrou default dispatch values from a case JSON.

    Parameters
    ----------
    case_json_path : str
        Path to case JSON file (e.g. hawaii.json).

    Returns
    -------
    pd.DataFrame with columns:
        gen_id, bus, p0, q0
    """
    with open(case_json_path, "r") as f:
        case = json.load(f)

    rows = []
    for dev in case.get("devices", []):
        if dev.get("class") != "Genrou":
            continue
        params = dev.get("params", {})
        ports = dev.get("ports", {})
        rows.append(
            {
                "gen_id": dev.get("id"),
                "bus": ports.get("bus"),
                "p0": params.get("p0", np.nan),
                "q0": params.get("q0", np.nan),
            }
        )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["bus", "gen_id"]).reset_index(drop=True)
    return df


def plot_genrou_dispatch(dispatch_df, top_n_labels=None, case_name=None):
    """
    Plot Genrou dispatch defaults (p0 and q0).

    Parameters
    ----------
    dispatch_df : pd.DataFrame
        Output of read_genrou_dispatch().
    top_n_labels : int or None
        Deprecated compatibility argument. Labels are now shown by default.
    case_name : str or None
        Optional case label (e.g. "hawaii", "illinois") shown in plot titles.
    """
    import plotly.express as px

    if dispatch_df.empty:
        print("No Genrou devices found.")
        return

    case_prefix = f"[{case_name}] " if case_name else ""

    fig1 = px.bar(
        dispatch_df,
        x="gen_id",
        y="p0",
        title=f"{case_prefix}Genrou default dispatch p0 ({len(dispatch_df)} generators)",
        hover_data=["bus", "q0"],
    )
    _ = fig1.update_layout(xaxis_title="Genrou id", yaxis_title="p0")
    _ = fig1.update_xaxes(type="category", tickangle=70)
    _ = fig1.show()

    fig2 = px.bar(
        dispatch_df,
        x="gen_id",
        y="q0",
        title=f"{case_prefix}Genrou default dispatch q0 ({len(dispatch_df)} generators)",
        hover_data=["bus", "p0"],
    )
    _ = fig2.update_layout(xaxis_title="Genrou id", yaxis_title="q0")
    _ = fig2.update_xaxes(type="category", tickangle=70)
    _ = fig2.show()

    fig3 = px.scatter(
        dispatch_df,
        x="p0",
        y="q0",
        text="gen_id",
        title=f"{case_prefix}Genrou dispatch operating points (p0, q0)",
        hover_data=["bus"],
    )
    _ = fig3.update_traces(textposition="top center")
    _ = fig3.show()


def patch_genrou_dispatch(case_json_path, updates, output_case_path=None):
    """
    Patch selected Genrou ids with new p0/q0 and write updated case JSON.

    Parameters
    ----------
    case_json_path : str
        Input case file path.
    updates : dict or list[dict]
        Either mapping by id:
            {"2_1": {"p0": 0.03, "q0": 0.01}, ...}
        or list form:
            [{"id": "2_1", "p0": 0.03, "q0": 0.01}, ...]
    output_case_path : str or None
        Output path. If None, overwrite case_json_path.

    Returns
    -------
    pd.DataFrame
        before/after table with columns:
        gen_id, p0_before, p0_after, q0_before, q0_after
    """
    with open(case_json_path, "r") as f:
        case = copy.deepcopy(json.load(f))

    if isinstance(updates, dict):
        normalized = []
        for gid, vals in updates.items():
            row = {"id": gid}
            row.update(vals or {})
            normalized.append(row)
    else:
        normalized = list(updates)

    seen = set()
    duplicates = set()
    for row in normalized:
        gid = row.get("id")
        if gid in seen:
            duplicates.add(gid)
        seen.add(gid)
    if duplicates:
        raise ValueError(f"Duplicate ids in updates: {sorted(duplicates)}")

    by_id = {
        d.get("id"): d for d in case.get("devices", []) if d.get("class") == "Genrou"
    }

    change_rows = []
    missing_ids = []
    for row in normalized:
        gid = row.get("id")
        if gid not in by_id:
            missing_ids.append(gid)
            continue

        dev = by_id[gid]
        dev.setdefault("params", {})

        p0_new = row.get("p0", dev["params"].get("p0"))
        q0_new = row.get("q0", dev["params"].get("q0"))

        if p0_new is not None and not isinstance(p0_new, Number):
            raise TypeError(f"p0 for {gid} must be numeric, got {type(p0_new)}")
        if q0_new is not None and not isinstance(q0_new, Number):
            raise TypeError(f"q0 for {gid} must be numeric, got {type(q0_new)}")

        old_p0 = dev["params"].get("p0", np.nan)
        old_q0 = dev["params"].get("q0", np.nan)

        if "p0" in row:
            dev["params"]["p0"] = float(p0_new)
        if "q0" in row:
            dev["params"]["q0"] = float(q0_new)

        change_rows.append(
            {
                "gen_id": gid,
                "p0_before": old_p0,
                "p0_after": dev["params"].get("p0", np.nan),
                "q0_before": old_q0,
                "q0_after": dev["params"].get("q0", np.nan),
            }
        )

    out_path = output_case_path or case_json_path
    with open(out_path, "w") as f:
        json.dump(case, f, indent=2)

    changes_df = pd.DataFrame(change_rows)
    if not changes_df.empty:
        changes_df = changes_df.sort_values("gen_id").reset_index(drop=True)

    print(f"Patched {len(change_rows)} Genrou devices in {out_path}")
    if missing_ids:
        print(f"WARNING: ids not found: {sorted(missing_ids)}")

    return changes_df


def attach_json_ids(case_data, json_path) -> None:
    """Attach GridKit JSON identifiers to case_data dataframes in-place.

    Adds columns:
      - bus_df:    ``json_bus_number`` (int, same as BUS_I — confirms correspondence)
      - gen_df:    ``json_gen_id`` (str like "2_1"; None for offline/absent gens)
      - branch_df: ``json_branch_id`` (str like "BR_1_2_1")

    Parameters
    ----------
    case_data : MatpowerCaseData
        Object whose ``.bus``, ``.gen``, ``.branch`` DataFrames will be mutated.
    json_path : str or Path
        Path to a GridKit JSON case file (hawaii.json style).
    """
    with open(json_path) as fh:
        case = json.load(fh)

    # --- bus: json number == BUS_I (trivially identical, added for hover clarity) ---
    case_data.bus["json_bus_number"] = case_data.bus["BUS_I"].astype(int)

    # --- generators: derive json_gen_id from rank within GEN_BUS group ---
    # Build set of ids that actually exist in the JSON (for validation / offline marking)
    json_gen_ids = {d["id"] for d in case["devices"] if d.get("class") == "Genrou"}

    gen_df = case_data.gen
    has_status = "GEN_STATUS" in gen_df.columns

    json_gen_id_col = []
    # rank counter per bus
    rank: dict[int, int] = {}
    for row in gen_df.itertuples():
        bus = int(row.GEN_BUS)
        rank[bus] = rank.get(bus, 0) + 1
        candidate = f"{bus}_{rank[bus]}"
        # Mark as None if offline in .m (GEN_STATUS=0) or absent from JSON
        if has_status and int(getattr(row, "GEN_STATUS", 1)) == 0:
            json_gen_id_col.append(None)
        elif candidate not in json_gen_ids:
            # absent from JSON even though status=1 — mark None
            json_gen_id_col.append(None)
        else:
            json_gen_id_col.append(candidate)

    case_data.gen["json_gen_id"] = json_gen_id_col

    # --- branches: derive json_branch_id from (F_BUS, T_BUS) rank ---
    parallel_count: dict[tuple, int] = {}
    json_branch_id_col = []
    for row in case_data.branch.itertuples():
        key = (int(row.F_BUS), int(row.T_BUS))
        parallel_count[key] = parallel_count.get(key, 0) + 1
        json_branch_id_col.append(f"BR_{key[0]}_{key[1]}_{parallel_count[key]}")

    case_data.branch["json_branch_id"] = json_branch_id_col


def get_case_path_for_editing(test_run_dir, example_dir, case_name):
    """
    Resolve case path for edits.

    Prefer the copied case file in test_run_dir. If it does not exist,
    fall back to the source/example directory.
    """
    run_case = os.path.join(test_run_dir, case_name)
    if os.path.exists(run_case):
        return run_case
    return os.path.join(example_dir, case_name)
