"""Pre-fault flatness check for Q-limit switched cases.

Runs the campaign's default disturbance on one sited operating point twice
(Q limits off, then on) and reports how far each monitored state drifts
during the pre-fault window [0, t_on). A switched case whose initial state
is consistent should be flat there, the same as the unlimited case.

    python qlim_dyn_check.py --dc-alpha-p 0.8
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
import campaign as cp  # noqa: E402
from qlim_pf_check import SCRATCH, pick_points  # noqa: E402
from sweep_worker import profiles  # noqa: E402


def prefault_drift(csv: Path, t_on: float):
    df = pd.read_csv(csv)
    tcol = next(c for c in df.columns if c.lower() in ("time", "t"))
    pre = df[df[tcol] < t_on]
    num = pre.drop(columns=[tcol]).select_dtypes(include=[np.number])
    drift = (num.iloc[-1] - num.iloc[0]).abs().sort_values(ascending=False)
    return drift, len(pre), num.iloc[0].abs().max()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--design", default="faultscreen1")
    p.add_argument("--dc-alpha-p", type=float, default=0.8)
    args = p.parse_args()

    cfg = cp.CampaignConfig(
        name="qlim_check", dc_alpha_p=args.dc_alpha_p, enforce_q_lims=False
    )
    at = profiles(cfg)
    r = pick_points(args.design, 1).iloc[0]
    demand, wind = at(int(r["hour"]))
    dc_bus, dc_mw = int(r["dc_bus"]), float(r["dc_mw"])

    for tag, flag in (("off", False), ("on", True)):
        cfg.enforce_q_lims = flag
        wd = SCRATCH / "dyn" / tag
        spec = {
            "demand": demand,
            "wind": wind,
            "n_commit": int(r["n_commit"]),
            "workdir": wd,
            "dc_bus": dc_bus,
            "dc_mw": dc_mw,
        }
        ((case, rec),) = cp.solve_points_batch(cfg, [spec])
        if case is None:
            print(tag, "solve failed:", rec.get("fail_reason"))
            continue
        out = cp.run_disturbance(cfg, case, cfg.disturbance, wd, with_cct=False)
        print(f"\n{tag}: status={out['status']} sim_rc={out['sim_rc']} "
              f"sim_errors={out['sim_errors']}")
        if out["status"] != "ok":
            continue
        drift, npre, scale = prefault_drift(out["_csv"], cfg.disturbance.t_on)
        print(f"  pre-fault samples={npre}  max|x| at t0={scale:.4g}")
        print("  largest pre-fault drift (first->last sample before t_on):")
        for col, val in drift.head(6).items():
            print(f"    {col:<28} {val:.3e}")


if __name__ == "__main__":
    main()
