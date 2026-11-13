"""Pilot: PM.jl AC OPF against DC OPF + PV->PQ switching, on the qlim_pf_check points.

Two arms on the same operating points (same load, wind, commitment):
  sw   DC OPF dispatch, then AC PF with Q-limit switching (campaign default)
  opf  PM.jl AC OPF (`ac-opf`): re-dispatches P and V, Q limits as constraints

For each point and arm it reports convergence, losses, Q-limit violations,
vmax, slack headroom, dispatch change against the switched case, warm solve
time (the same case repeated in one Julia process), and optionally the
pre-fault drift of the dynamic initial state (`--dyn`).

    python acopf_pilot.py --design faultscreen1 --points 4 --dc-alpha-p 0.8 --dyn

Writes to the scratch `qlim_check/acopf_pilot/` tree, not the campaign store.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
import aux_to_m  # noqa: E402
import campaign as cp  # noqa: E402
import m_to_case_utils as mu  # noqa: E402
import pm_utils  # noqa: E402
from qlim_dyn_check import prefault_drift  # noqa: E402
from qlim_pf_check import SCRATCH, gen_violations, pick_points  # noqa: E402
from sweep_worker import profiles  # noqa: E402

ROOT = SCRATCH / "acopf_pilot"
WARM_COPIES = 5


def name_of(r):
    dc_bus = None if pd.isna(r["dc_bus"]) else int(r["dc_bus"])
    dc_mw = 0.0 if pd.isna(r["dc_mw"]) else float(r["dc_mw"])
    return f"h{int(r['hour'])}_n{int(r['n_commit'])}_dc{dc_bus}_{dc_mw:g}"


def arm_metrics(cfg, frames, wd_ac, dc_pg=None):
    """Common metrics for one solved .m (`frames` from m_to_frames)."""
    bus, gen = frames["bus"], frames["gen"]
    on = gen[gen["GEN_STATUS"] == 1]
    slack = bus.loc[bus["BUS_TYPE"] == 3, "BUS_I"].astype(int).tolist()
    s_gen = on[on["GEN_BUS"].isin(slack)]
    out = {
        "losses_mw": round(float(on["PG"].sum() - bus["PD"].sum()), 2),
        "q_viol": len(gen_violations(frames)),
        "vmax": round(float(bus["VM"].max()), 4),
        "vmin": round(float(bus["VM"].min()), 4),
        "n_vm_ge_1099": int((bus["VM"] >= 1.099).sum()),
        "slack_pg": round(float(s_gen["PG"].sum()), 2),
        "slack_headroom_mw": round(float((s_gen["PMAX"] - s_gen["PG"]).sum()), 2),
        "max_over_pmax": round(float((on["PG"] - on["PMAX"]).clip(lower=0).max()), 3),
    }
    if dc_pg is not None and len(dc_pg) == len(on):
        d = np.abs(on["PG"].values - dc_pg)
        out["pg_change_sum_mw"] = round(float(d.sum()), 2)
        out["pg_change_max_mw"] = round(float(d.max()), 2)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--design", default="faultscreen1")
    p.add_argument("--points", type=int, default=4)
    p.add_argument("--dc-alpha-p", type=float, default=0.8)
    p.add_argument("--dyn", action="store_true", help="run pre-fault drift per arm")
    args = p.parse_args()

    cfg = cp.CampaignConfig(
        name="qlim_check", dc_alpha_p=args.dc_alpha_p, enforce_q_lims=True
    )
    at = profiles(cfg)
    pts = pick_points(args.design, args.points)
    order = cp.merit_order(cfg)

    # Arm sw: campaign path, switching on.
    specs, prepped = [], []
    for _, r in pts.iterrows():
        demand, wind = at(int(r["hour"]))
        dc_bus = None if pd.isna(r["dc_bus"]) else int(r["dc_bus"])
        dc_mw = 0.0 if pd.isna(r["dc_mw"]) else float(r["dc_mw"])
        nm = name_of(r)
        frames = cp.apply_hour(cfg, demand, wind)
        if dc_bus is not None and dc_mw > 0:
            frames = cp.add_datacenter(cfg, frames, dc_bus, dc_mw)
        frames = cp.commit_n(cfg, frames, int(r["n_commit"]), order)
        prepped.append((nm, frames, dc_bus, dc_mw))
        specs.append({
            "demand": demand, "wind": wind, "n_commit": int(r["n_commit"]),
            "workdir": ROOT / "sw" / nm, "dc_bus": dc_bus, "dc_mw": dc_mw,
        })
    sw = cp.solve_points_batch(cfg, specs)

    # Arm opf: AC OPF on the same frames (no DC loss inflation: AC OPF has losses).
    cases, opf_dirs = [], []
    for nm, frames, _, _ in prepped:
        wd = ROOT / "opf" / nm
        wd.mkdir(parents=True, exist_ok=True)
        aux_to_m.write_matpower_m(frames, wd / "opfin.m", case_name="c")
        cases.append({"input": str(wd / "opfin.m"), "output": str(wd / "opf.m"),
                      "pf_type": "ac-opf"})
        opf_dirs.append(wd)
    # Same cases again: the repeats run warm (Julia already compiled).
    warm_dir = ROOT / "opf" / "_warm"
    warm_dir.mkdir(parents=True, exist_ok=True)
    for k, (nm, frames, _, _) in enumerate(prepped[:WARM_COPIES]):
        aux_to_m.write_matpower_m(frames, warm_dir / f"w{k}_in.m", case_name="c")
        cases.append({"input": str(warm_dir / f"w{k}_in.m"),
                      "output": str(warm_dir / f"w{k}_out.m"), "pf_type": "ac-opf"})
    jb = pm_utils.resolve_julia()
    proj = pm_utils.PF_BACKENDS["pm"]["project"]
    script = pm_utils.PF_BACKENDS["pm"]["script"]
    recs, err, rc = pm_utils.run_pm_solve_batch(jb, proj, script, cases, ROOT / "opf" / "_batch")
    if rc != 0 or len(recs) != len(cases):
        raise SystemExit(f"AC OPF batch failed rc={rc}\n{err[-2000:]}")
    cold, warm = recs[: len(prepped)], recs[len(prepped):]

    rows, dyn_rows = [], []
    for j, (nm, frames, dc_bus, dc_mw) in enumerate(prepped):
        sw_case, sw_rec = sw[j]
        sw_ac = aux_to_m.m_to_frames(ROOT / "sw" / nm / "ac.m")
        sw_dc_pg = aux_to_m.m_to_frames(ROOT / "sw" / nm / "dc.m")["gen"]
        sw_on = sw_ac["gen"][sw_ac["gen"]["GEN_STATUS"] == 1]
        base = {"point": nm, "dc_bus": dc_bus, "dc_mw": dc_mw}
        rows.append({**base, "arm": "sw", "ok": sw_case is not None,
                     "solve_time_s": None, "warm_time_s": None,
                     **arm_metrics(cfg, sw_ac, None)})
        rec = cold[j]
        if not rec.get("ok"):
            rows.append({**base, "arm": "opf", "ok": False,
                         "status": rec.get("status"), "solve_time_s": rec.get("solve_time"),
                         "warm_time_s": None})
            continue
        opf = aux_to_m.m_to_frames(opf_dirs[j] / "opf.m")
        m = arm_metrics(cfg, opf, None, dc_pg=sw_on["PG"].values)
        w = warm[j] if j < WARM_COPIES else {}
        rows.append({**base, "arm": "opf", "ok": True, "status": rec.get("status"),
                     "solve_time_s": round(float(rec["solve_time"]), 3),
                     "warm_time_s": round(float(w["solve_time"]), 3) if w.get("ok") else None,
                     **m})

        if args.dyn:
            for arm, fr_dir in (("sw", ROOT / "sw" / nm), ("opf", opf_dirs[j])):
                ac_m = fr_dir / ("ac.m" if arm == "sw" else "opf.m")
                solved = aux_to_m.m_to_frames(ac_m)
                case, _ = mu.build_case_from_m(cfg.spec, cfg.base_case, ac_m, dyr=cfg.dyr)
                case, _ = cp.wind_to_negative_load(cfg, case, solved)
                if dc_bus is not None and dc_mw > 0:
                    case = cp.split_out_datacenter(
                        cfg, case, dc_bus, dc_mw,
                        alpha_p=cfg.dc_alpha_p, alpha_i=1.0 - cfg.dc_alpha_p,
                    )
                dwd = ROOT / "dyn" / arm / nm
                out = cp.run_disturbance(cfg, case, cfg.disturbance, dwd, with_cct=False)
                if out["status"] != "ok":
                    dyn_rows.append({"point": nm, "arm": arm, "status": out["status"]})
                    continue
                drift, npre, scale = prefault_drift(out["_csv"], cfg.disturbance.t_on)
                dyn_rows.append({"point": nm, "arm": arm, "status": "ok",
                                 "npre": npre, "max_drift": float(drift.iloc[0]),
                                 "worst_col": str(drift.index[0])})

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    df = pd.DataFrame(rows)
    print(df.to_string(), flush=True)
    ROOT.mkdir(parents=True, exist_ok=True)
    df.to_csv(ROOT / f"summary_ap{args.dc_alpha_p:g}.csv", index=False)
    if dyn_rows:
        ddf = pd.DataFrame(dyn_rows)
        print("\npre-fault drift", flush=True)
        print(ddf.to_string(), flush=True)
        ddf.to_csv(ROOT / f"dyn_ap{args.dc_alpha_p:g}.csv", index=False)


if __name__ == "__main__":
    main()
