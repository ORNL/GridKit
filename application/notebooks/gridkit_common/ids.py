"""Helpers linking MATPOWER case tables back to GridKit JSON case files."""

from __future__ import annotations

import json
import os


def attach_json_ids(case_data, json_path) -> None:
    """Attach GridKit JSON identifiers to case_data dataframes in-place.

    Adds columns:
      - bus_df:    ``json_bus_number`` (int, same as BUS_I — confirms correspondence)
      - gen_df:    ``json_gen_id`` (str like "2_1"; None for offline/absent gens)
      - branch_df: ``json_branch_id`` (str like "BR_1_2_1")

    Parameters
    ----------
    case_data : m_viz.MatpowerCaseData
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

    # --- gen: also store the JSON bus "number" for each gen (== GEN_BUS == BUS_I) ---
    # Only set when the gen is present in JSON; None for offline/absent gens.
    case_data.gen["json_gen_bus_number"] = [
        int(row.GEN_BUS) if gid is not None else None
        for row, gid in zip(gen_df.itertuples(), json_gen_id_col)
    ]

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
