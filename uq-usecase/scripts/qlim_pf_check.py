"""Check PV->PQ Q-limit switching in the campaign AC PF on real operating points.

Solves a few points of an existing campaign design twice (Q limits off, then
on) with `campaign.solve_points_batch`, and reports for each point:
  - units fixed at a limit, Julia solves used, losses, vmin/vmax
  - whether every in-service non-slack gen is inside [QMIN, QMAX] after
    switching (units fixed at a limit are checked against that limit)
  - the Q-limit violations the unlimited solution had

    python qlim_pf_check.py --design faultscreen1 --points 4 --dc-alpha-p 0.8

Writes per-point rows to stdout and the solved cases to the scratch
`qlim_check/` tree (not the campaign store).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
import aux_to_m  # noqa: E402
import campaign as cp  # noqa: E402
from sweep_worker import profiles  # noqa: E402

CONDS = "/kfs2/projects/scidac/scidac-data/gridkit-campaigns"
SCRATCH = Path("/scratch/isatkaus/gridkit-campaigns/qlim_check")


def pick_points(design: str, n: int) -> pd.DataFrame:
    """First `n` operating points that carry a datacenter, plus one without."""
    C = pd.read_csv(f"{CONDS}/{design}/conditions.csv")
    keys = ["hour", "n_commit", "dc_bus", "dc_mw"]
    pts = C.drop_duplicates(keys)[keys]
    sited = pts[pts["dc_bus"].notna()].head(n)
    base = pts[pts["dc_bus"].isna()].head(1)
    return pd.concat([sited, base], ignore_index=True)


def gen_violations(frames, switched_qg=None):
    """In-service non-slack gens whose QG is outside [QMIN, QMAX]."""
    bus, gen = frames["bus"], frames["gen"]
    slack = set(bus.loc[bus["BUS_TYPE"] == 3, "BUS_I"].astype(int))
    bad = []
    for i in range(len(gen)):
        g = gen.iloc[i]
        if g["GEN_STATUS"] != 1 or int(g["GEN_BUS"]) in slack:
            continue
        q = float(g["QG"])
        if q > g["QMAX"] + 0.1 or q < g["QMIN"] - 0.1:
            bad.append((i, int(g["GEN_BUS"]), round(q, 2), g["QMIN"], g["QMAX"]))
    return bad


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--design", default="faultscreen1")
    p.add_argument("--points", type=int, default=4)
    p.add_argument("--dc-alpha-p", type=float, default=0.8)
    args = p.parse_args()

    cfg = cp.CampaignConfig(
        name="qlim_check", dc_alpha_p=args.dc_alpha_p, enforce_q_lims=False
    )
    at = profiles(cfg)
    pts = pick_points(args.design, args.points)
    print(pts.to_string(), flush=True)

    def specs_for(tag):
        out = []
        for _, r in pts.iterrows():
            demand, wind = at(int(r["hour"]))
            dc_bus = None if pd.isna(r["dc_bus"]) else int(r["dc_bus"])
            dc_mw = 0.0 if pd.isna(r["dc_mw"]) else float(r["dc_mw"])
            name = f"h{int(r['hour'])}_n{int(r['n_commit'])}_dc{dc_bus}_{dc_mw:g}"
            out.append(
                {
                    "demand": demand,
                    "wind": wind,
                    "n_commit": int(r["n_commit"]),
                    "workdir": SCRATCH / tag / name,
                    "dc_bus": dc_bus,
                    "dc_mw": dc_mw,
                }
            )
        return out

    results = {}
    for tag, flag in (("off", False), ("on", True)):
        cfg.enforce_q_lims = flag
        results[tag] = cp.solve_points_batch(cfg, specs_for(tag))

    for j, (case_off, rec_off) in enumerate(results["off"]):
        case_on, rec_on = results["on"][j]
        wd_off = specs_for("off")[j]["workdir"]
        wd_on = specs_for("on")[j]["workdir"]
        print(f"\npoint {j}: dc_bus={rec_on['dc_bus']} dc_mw={rec_on['dc_mw']}")
        for lbl, rec in (("off", rec_off), ("on", rec_on)):
            print(
                f"  {lbl:>3}: status={rec.get('status', 'ok')} "
                f"q_lims={rec.get('q_lims')} q_fixed={rec.get('q_fixed')} "
                f"q_solves={rec.get('q_solves')} losses={rec.get('losses_mw')} "
                f"vmin={rec.get('vmin')} vmax={rec.get('vmax')}"
            )
        if rec_on.get("status") == "solve_failed":
            print("  switched solve FAILED:", rec_on.get("fail_reason"))
            continue
        off_frames = aux_to_m.m_to_frames(wd_off / "ac.m")
        on_frames = aux_to_m.m_to_frames(wd_on / "ac.m")
        print(f"  unlimited: Q-limit violations = {len(gen_violations(off_frames))}")
        print(f"  switched : Q-limit violations = {len(gen_violations(on_frames))}")
        v = gen_violations(on_frames)
        for row in v[:5]:
            print("    ", row)


if __name__ == "__main__":
    main()
