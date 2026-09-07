"""Run-directory management and result collection for UQ sampling sweeps."""

from __future__ import annotations

import copy
import json
import os
import subprocess

import pandas as pd

# ---------------------------------------------------------------------------
# Monitorable variables by element class
# Source: GridKit docs INPUT_FORMAT.md
#   - Bus classes table   (init vars + other variables available to monitor)
#   - Device classes table (Variables available to monitor column)
#
# Reference for building the `monitors_by_class` argument to make_run_dir().
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
