"""Backfill a finished campaign's dispatch from scratch, then verify it against the records.

usage: python extract_dispatch.py <campaign>
"""

import argparse
import json
from pathlib import Path

import numpy as np

import campaign as cp


def backfill(scratch: Path, store: Path) -> None:
    tags = sorted(p.stem for p in (store / "records").glob("*.json"))
    for tag in tags:
        cp.save_dispatch(scratch, store, tag)
    n = len(list((store / "dispatch").glob("*.parquet")))
    print(f"records {len(tags)}; dispatch files {n}")


def verify(store: Path) -> None:
    import pandas as pd

    recs = [json.loads(p.read_text()) for p in sorted((store / "records").glob("*.json"))]
    by_op = {}
    for r in recs:
        if r["status"] == "ok":
            by_op.setdefault(r["tag"].split("__")[0], r)
    missing, unsolved, mismatch, resid = [], 0, 0, []
    for op, r in by_op.items():
        f = store / "dispatch" / f"{op}.parquet"
        if not f.exists():
            missing.append(op)
            continue
        df = pd.read_parquet(f)
        unsolved += int(not df["ac_solved"].all())
        mismatch += int((df["status_ac"] != df["status_dc"]).sum())
        on = df["status_ac"] == 1
        resid.append(float(df.loc[on, "pg_ac"].sum() - r["demand_mw"] - r["losses_mw"]))
    failed_sim = {r["tag"].split("__")[0] for r in recs if r["status"] == "sim_failed"}
    failed_sim_missing = [op for op in failed_sim if not (store / "dispatch" / f"{op}.parquet").exists()]
    print(f"ok operating points {len(by_op)}; missing dispatch {len(missing)}")
    print(f"sim_failed operating points {len(failed_sim)}; without dispatch {len(failed_sim_missing)}")
    print(f"ok points with AC unsettled {unsolved}; gen status dc vs ac mismatches {mismatch}")
    r = np.array(resid)
    print(f"PG_ac - demand - losses: max |resid| {np.abs(r).max():.3f} MW, mean {r.mean():.3f} MW")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("campaign")
    a = ap.parse_args()
    scratch = cp.SCRATCH_ROOT / a.campaign
    store = cp.STORE_ROOT / a.campaign
    cp.prepare_storage(scratch, store)
    backfill(scratch, store)
    verify(store)


if __name__ == "__main__":
    main()
