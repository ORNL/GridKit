"""Zoom-in on IDACalcIC failures for one failing dispatch.

Takes one scratch case (a dispatch that fails IDACalcIC in a campaign) and
re-runs it with only the fault description changed. Three fault modes:

  open    BusFault present, state0=false, no events. Fault never closes.
  closed  BusFault present, state0=true, no events. Fault closed from t=0.
  event   BusFault present, state0=false, fault_on at t=1.0 (campaign setup).
          tmax is 1.2 so the event actually fires.

`--stage bus` sweeps the BusFault over every bus at each alphaP (plus a
no-BusFault control per alphaP). `--stage solver` runs integrator options
(consistent IC type, tolerances, max order) at the case's own fault bus in
event mode.

Each run goes to /scratch/isatkaus/gridkit-campaigns/ida_zoom/<name>/<label>/.
results.csv in <name>/ has one row per run. status is ok / IDACalcIC /
IDASolve / other, from exit code, output csv and log.

    python ida_zoom.py --src <scratch case dir> --name ic181 --stage bus
"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")

import campaign as cp  # noqa: E402
import solver_utils as su  # noqa: E402

RUNNER = Path.home() / "gridkit/build/application/PhasorDynamics/DynamicSimulation"
BASE_SOLVER = Path(
    "/home/isatkaus/gridkit/examples/PhasorDynamics/Validation/ACTIVSg200/"
    "ACTIVSg200.solver.json"
)
OUT_ROOT = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom")
FAULT_ON_T = 1.0
MW_LEVELS = [25, 50, 100, 150, 200, 250, 300, 400, 600]

# mode -> (state0, fires_event, tmax)
MODES = {
    "open": (False, False, 0.2),
    "closed": (True, False, 0.2),
    "event": (False, True, 1.2),
}


def set_alpha(case, alpha):
    """Datacenter ZIP to alphaP = alpha, alphaI = 1 - alpha (no alphaZ, as campaign.py)."""
    for d in case["devices"]:
        if d["class"] == "LoadZIP" and d.get("id", "").startswith("datacenter_"):
            d["params"]["alphaP"] = float(alpha)
            d["params"]["alphaI"] = 1.0 - float(alpha)
            return
    raise RuntimeError("no datacenter LoadZIP in case")


def set_fault(case, bus, state0, x=None):
    """Point the single BusFault at `bus` (or remove it when bus is None). x overrides X."""
    idx = [i for i, d in enumerate(case["devices"]) if d["class"] == "BusFault"]
    if len(idx) != 1:
        raise RuntimeError(f"expected one BusFault, found {len(idx)}")
    if bus is None:
        del case["devices"][idx[0]]
        return
    dev = case["devices"][idx[0]]
    dev["ports"]["bus"] = int(bus)
    dev["params"]["state0"] = bool(state0)
    if x is not None:
        dev["params"]["X"] = float(x)


def run_one(job):
    """One dynamic run. job keys: label, src, outdir, fault_bus, alpha_p, mode, solver_mods."""
    label, src, outdir = job["label"], Path(job["src"]), Path(job["outdir"])
    wd = outdir / label
    wd.mkdir(parents=True, exist_ok=True)
    for f in wd.glob("*.csv"):
        f.unlink()

    case = json.loads((src / "case.json").read_text())
    set_alpha(case, job["alpha_p"])
    fb = job["fault_bus"]
    if fb == "none":
        set_fault(case, None, None)
        events, tmax = [], 0.2
    else:
        state0, fires, tmax = MODES[job["mode"]]
        set_fault(case, fb, state0, job.get("fault_x"))
        eid = su.fault_element_id(case, fb)
        events = (
            [{"time": FAULT_ON_T, "type": "fault_on", "element_id": eid}] if fires else []
        )
    if job.get("dc_mw") is not None:
        for d in case["devices"]:
            if d["class"] == "LoadZIP" and d.get("id", "").startswith("datacenter_"):
                f = (job["dc_mw"] / 100.0) / d["params"]["Pnom"]
                d["params"]["Pnom"] *= f
                d["params"]["Qnom"] *= f
    if job.get("dc_zero"):
        for d in case["devices"]:
            if d["class"] == "LoadZIP" and d.get("id", "").startswith("datacenter_"):
                d["params"]["Pnom"] = 0.0
                d["params"]["Qnom"] = 0.0
    (wd / "case.json").write_text(json.dumps(case))

    solver = su.apply_solver_overrides(
        json.loads(BASE_SOLVER.read_text()),
        {
            "system_model_file": "case.json",
            "tmax": tmax,
            "events": events,
            **job["solver_mods"],
        },
    )
    (wd / "solver.json").write_text(json.dumps(solver, indent=2))

    t0 = time.time()
    r = subprocess.run(
        [str(RUNNER), "solver.json"],
        cwd=wd,
        capture_output=True,
        text=True,
        env=cp.runner_env(),
    )
    log = r.stdout + r.stderr
    (wd / "sim.log").write_text(log)
    csv_ok = next(iter(wd.glob("*.csv")), None) is not None
    if r.returncode == 0 and csv_ok:
        status = "ok"
    elif "IDACalcIC" in log:
        status = "IDACalcIC"
    elif "IDASolve" in log:
        status = "IDASolve"
    else:
        status = "other"
    sm = job["solver_mods"]
    return {
        "label": label,
        "mode": job["mode"],
        "fault_bus": fb,
        "alpha_p": job["alpha_p"],
        "ic_type": sm.get("consistent_ic_type", "ya_ydp"),
        "rel_tol": sm.get("rel_tol", ""),
        "abs_tol": sm.get("abs_tol", ""),
        "max_order": sm.get("max_order", ""),
        "status": status,
        "rc": r.returncode,
        "seconds": round(time.time() - t0, 2),
    }


def build_jobs(stage, src, name, alphas, fault_bus0, buses=None):
    jobs = []

    def job(label, fb, alpha, mode, sm, **extra):
        return dict(
            label=label,
            src=str(src),
            outdir=str(OUT_ROOT / name),
            fault_bus=fb,
            alpha_p=alpha,
            mode=mode,
            solver_mods=sm,
            **extra,
        )

    if stage == "bus":
        c0 = json.loads((src / "case.json").read_text())
        buses = sorted(b["number"] for b in c0["buses"])
        for a in alphas:
            jobs.append(job(f"nofault_a{a:g}", "none", a, "open", {}))
            for b in buses:
                for mode in MODES:
                    jobs.append(job(f"{mode}_fb{b}_a{a:g}", b, a, mode, {}))
    elif stage == "mw":
        # datacenter size sweep at fixed fault; dc_mw in MW, P and Q scaled together
        for b in buses or [fault_bus0]:
            for a in alphas:
                for mw in MW_LEVELS:
                    jobs.append(job(f"mw{mw:g}_fb{b}_a{a:g}", b, a, "event", {}, dc_mw=mw))
    elif stage == "solver":
        # (label, solver overrides, extra job fields)
        variants = [
            ("default", {}, {}),
            ("ic_y", {"consistent_ic_type": "y"}, {}),
            ("rel1e-4", {"rel_tol": 1e-4}, {}),
            ("rel1e-8_abs1e-10", {"rel_tol": 1e-8, "abs_tol": 1e-10}, {}),
            ("max_order_1", {"max_order": 1}, {}),
            ("max_order_2", {"max_order": 2}, {}),
            ("ic_y_rel1e-4", {"consistent_ic_type": "y", "rel_tol": 1e-4}, {}),
            ("X0.05", {}, {"fault_x": 0.05}),
            ("X0.1", {}, {"fault_x": 0.1}),
            ("dc_zero", {}, {"dc_zero": True}),
        ]
        for b in buses or [fault_bus0]:
            for a in alphas:
                for label, sm, extra in variants:
                    jobs.append(job(f"{label}_fb{b}_a{a:g}", b, a, "event", sm, **extra))
                # alphaP = 0 is constant current (0/100/0), a load-type control.
                jobs.append(job(f"ci_fb{b}_a0", b, 0.0, "event", {}))
    else:
        raise ValueError(stage)
    return jobs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--src", required=True, help="scratch case dir of a failing dispatch")
    p.add_argument("--name", required=True)
    p.add_argument("--stage", choices=["bus", "solver", "mw"], required=True)
    p.add_argument("--alphas", default="0.2,1.0")
    p.add_argument("--fault-bus", type=int, default=None,
                   help="fault bus for --stage solver (default: the case's own BusFault)")
    p.add_argument("--buses", default=None,
                   help="comma list of fault buses for --stage solver (default: --fault-bus)")
    p.add_argument("--workers", type=int, default=16)
    args = p.parse_args()

    src = Path(args.src)
    alphas = [float(x) for x in args.alphas.split(",")]
    fb0 = args.fault_bus
    if fb0 is None:
        c = json.loads((src / "case.json").read_text())
        fb0 = next(d["ports"]["bus"] for d in c["devices"] if d["class"] == "BusFault")

    outdir = OUT_ROOT / args.name
    outdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(BASE_SOLVER, outdir / "ACTIVSg200.solver.json")

    buses = [int(x) for x in args.buses.split(",")] if args.buses else None
    jobs = build_jobs(args.stage, src, args.name, alphas, fb0, buses)
    print(f"{len(jobs)} runs", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(run_one, jobs))

    with open(outdir / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote", outdir / "results.csv", flush=True)


if __name__ == "__main__":
    main()
