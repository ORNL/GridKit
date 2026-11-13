"""Probe the dispatch-stage failures in the fault sweep.

    python dispatch_pf_probe.py

Solves four operating points from faultscreen1 that give solve_failed in the
campaign records, once with Q limits off and once with PV->PQ switching on.
Every Julia batch call is logged (stage, pass number, return code, first error
line), so a switched failure can be placed in a specific pass. Writes nothing
to the campaign store; cases go to the scratch `dispatch_probe/` tree.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
import campaign as cp  # noqa: E402
from sweep_worker import profiles  # noqa: E402

CONDS = "/kfs2/projects/scidac/scidac-data/gridkit-campaigns"
SCRATCH = Path("/scratch/isatkaus/gridkit-campaigns/dispatch_probe")
# (hour, n_commit, dc_bus, dc_mw): the solve_failed dispatches in
# faultscreen1_ap08_qlim, one per failure class.
POINTS = [
    (7536, 18, 148, 200.0),  # DC OPF, capacity short by 10.2 MW
    (6483, 21, 108, 200.0),  # AC PF, fails with and without switching
    (6483, 21, 148, 200.0),  # AC PF, fails only with switching
    (7536, 23, 108, 200.0),  # AC PF, fails only with switching
]

_orig = cp.pm_utils.run_pm_solve_batch
calls = []


def logged(jb, proj, script, cases, batch_dir):
    recs, err, rc = _orig(jb, proj, script, cases, batch_dir)
    pf = cases[0]["pf_type"] if cases else "?"
    first = next((ln for ln in (err or "").splitlines() if "ERROR" in ln), "")
    calls.append((pf, len(cases), rc, sum(bool(r.get("ok")) for r in recs), first[:140]))
    return recs, err, rc


cp.pm_utils.run_pm_solve_batch = logged


def main():
    cfg = cp.CampaignConfig(name="dispatch_probe", dc_alpha_p=0.8)
    at = profiles(cfg)
    for flag in (False, True):
        cfg.enforce_q_lims = flag
        label = "switched" if flag else "unlimited"
        specs = []
        for h, n, b, mw in POINTS:
            demand, wind = at(h)
            specs.append(
                {
                    "demand": demand,
                    "wind": wind,
                    "n_commit": n,
                    "workdir": SCRATCH / label / f"h{h}_n{n}_dc{b}",
                    "dc_bus": b,
                    "dc_mw": mw,
                }
            )
        calls.clear()
        out = cp.solve_points_batch(cfg, specs)
        print(f"\n== {label}")
        for (h, n, b, mw), (case, rec) in zip(POINTS, out):
            print(
                f"  h{h} n{n} dc{b}: {rec.get('status', 'ok')} "
                f"{rec.get('fail_stage', '')} {rec.get('fail_reason', '')} "
                f"net={rec.get('net_mw')} q_solves={rec.get('q_solves')}"
            )
        print("  julia batch calls (stage, cases, rc, ok, first error):")
        for c in calls:
            print("   ", c)


if __name__ == "__main__":
    main()
