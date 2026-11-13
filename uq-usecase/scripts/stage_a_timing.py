"""Time stage A (batched DC OPF + AC PF) for a slice of unique dispatches.

Counts every Julia launch and its wall time, to test where the per-worker
fixed cost of sweep_worker.py goes. Stage A only: no dynamics, no CCT.

    python stage_a_timing.py --conditions C.csv --start 0 --n 4 --out run.json

--start and --n index the unique (hour, n_commit, dc_bus, dc_mw) dispatches
of --conditions in file order, wrapping modulo the number of dispatches.
"""

import argparse
import json
import socket
import sys
import time
from pathlib import Path

import pandas as pd

T_IMPORT0 = time.time()
sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402
import pm_utils  # noqa: E402

T_IMPORT = time.time() - T_IMPORT0

DATA = "/kfs2/projects/scidac/scidac-data/ACTIVSg200/grid-v1"
KEY = ["hour", "n_commit", "dc_bus", "dc_mw"]


def profiles(cfg):
    dem = pd.read_csv(f"{DATA}/demand.csv")
    wnd = pd.read_csv(f"{DATA}/wind.csv")
    demc = [c for c in dem.columns if c != "DateTime"]
    wndc = [c for c in wnd.columns if c != "DateTime"]
    wg = {c: cfg.gid[int(c.replace("gen", "")) - 1] for c in wndc}

    def at(h):
        return (
            {int(x.replace("bus", "")): float(dem.loc[h, x]) for x in demc},
            {wg[x]: float(wnd.loc[h, x]) for x in wndc},
        )

    return at


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--conditions", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--root", default="/scratch/isatkaus/gridkit-timing")
    args = p.parse_args()

    # every Julia launch goes through pm_utils.run_pm_solve_batch; time each one
    launches = []
    orig = pm_utils.run_pm_solve_batch

    def timed(*a, **k):
        t = time.time()
        out = orig(*a, **k)
        cases = a[3] if len(a) > 3 else k["cases"]
        launches.append(
            {"cases": len(cases), "pf_type": cases[0].get("pf_type"), "wall_s": round(time.time() - t, 2), "rc": out[2]}
        )
        return out

    pm_utils.run_pm_solve_batch = timed

    tag = Path(args.out).stem
    root = Path(args.root) / tag
    cfg = cp.CampaignConfig(
        name=tag,
        dc_alpha_p=0.8,
        enforce_q_lims=True,
        scratch=root / "scratch",
        store=root / "store",
    )

    C = pd.read_csv(args.conditions)
    keys = C[KEY].drop_duplicates().reset_index(drop=True)
    idx = [(args.start + i) % len(keys) for i in range(args.n)]
    chosen = keys.iloc[idx]

    at = profiles(cfg)
    specs = []
    for _, r in chosen.iterrows():
        demand, wind = at(int(r.hour))
        op = cp.OperatingPoint(
            hour=int(r.hour),
            n_commit=int(r.n_commit),
            dc_bus=None if pd.isna(r.dc_bus) else int(r.dc_bus),
            dc_mw=0.0 if pd.isna(r.dc_mw) else float(r.dc_mw),
        )
        specs.append(
            {
                "demand": demand,
                "wind": wind,
                "n_commit": int(r.n_commit),
                "workdir": cfg.scratch / op.tag,
                "dc_bus": op.dc_bus,
                "dc_mw": op.dc_mw,
            }
        )

    t0 = time.time()
    solved = cp.solve_points_batch(cfg, specs)
    t_stage_a = time.time() - t0
    n_ok = sum(1 for case, _ in solved if case is not None)

    result = {
        "host": socket.gethostname(),
        "conditions": args.conditions,
        "start": args.start,
        "n_dispatch": len(specs),
        "n_ok": n_ok,
        "t_import_s": round(T_IMPORT, 2),
        "t_stage_a_s": round(t_stage_a, 2),
        "n_launches": len(launches),
        "launches": launches,
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(result, indent=1))
    print(
        f"{tag}: {len(specs)} dispatches, {n_ok} ok, stage A {t_stage_a:.1f}s, "
        f"{len(launches)} launches, import {T_IMPORT:.1f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
