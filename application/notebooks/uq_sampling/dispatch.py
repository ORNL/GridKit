"""Genrou dispatch inspection and patching (aleatoric UQ prep)."""

from __future__ import annotations

import copy
import json
from numbers import Number

import numpy as np
import pandas as pd


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
