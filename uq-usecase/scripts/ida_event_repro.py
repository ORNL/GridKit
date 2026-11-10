"""Event-time reproducer for the IDA IC failure (ida_ic_issue.md, test 4).

Writes self-contained GridKit case directories (case.json + solver.json) and runs
the DynamicSimulation binary on each. The case is the fault-181 dispatch with a
200 MW datacenter at bus 148, alphaP = 1.0. Only the BusFault bus and the event
time change between cases.

  ev<t>_fb<b>   BusFault at bus b, state0 false (open at t = 0), fault_on at t,
                tmax = t + 0.2. t in {0.5, 1.0, 2.0}, b in {1, 94}.
  noev_fb<b>    same open fault, no events, tmax 2.2. Control, b in {1, 94}.

If the failing solve is the fault-on re-initialisation at the event, bus 1 fails
at every t and bus 94 passes at every t. If the failure is tied to t = 0 or to
tmax, the pattern changes. Writes results to
  /scratch/isatkaus/gridkit-campaigns/ida_zoom/event_repro/results.csv
and leaves ida_zoom.py unchanged.

    python ida_event_repro.py
"""

import csv
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")

import ida_zoom as iz  # noqa: E402
import solver_utils as su  # noqa: E402

SRC = Path("/scratch/isatkaus/gridkit-campaigns/faultscreen1_ap1/h6483_n21_dc148-200/b181_d0.1")
OUT = iz.OUT_ROOT / "event_repro"
BUSES = [1, 94]
EVENT_TIMES = [0.5, 1.0, 2.0]


def make_job(label, bus, event_t):
    return {"label": label, "bus": bus, "event_t": event_t}


def build_jobs():
    jobs = []
    for b in BUSES:
        for t in EVENT_TIMES:
            jobs.append(make_job(f"ev{t:g}_fb{b}", b, t))
        jobs.append(make_job(f"noev_fb{b}", b, None))
    return jobs


def run_job(job):
    wd = OUT / job["label"]
    wd.mkdir(parents=True, exist_ok=True)
    case = json.loads((SRC / "case.json").read_text())
    iz.set_alpha(case, 1.0)
    iz.set_fault(case, job["bus"], state0=False)
    if job["event_t"] is None:
        events, tmax = [], 2.2
    else:
        eid = su.fault_element_id(case, job["bus"])
        events = [{"time": job["event_t"], "type": "fault_on", "element_id": eid}]
        tmax = job["event_t"] + 0.2
    (wd / "case.json").write_text(json.dumps(case))
    solver = su.apply_solver_overrides(
        json.loads(iz.BASE_SOLVER.read_text()),
        {"system_model_file": "case.json", "tmax": tmax, "events": events},
    )
    (wd / "solver.json").write_text(json.dumps(solver, indent=2))

    r = subprocess.run(
        [str(iz.RUNNER), "solver.json"], cwd=wd, capture_output=True, text=True,
        env=iz.cp.runner_env(),
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
    return {
        "label": job["label"],
        "fault_bus": job["bus"],
        "event_t": "" if job["event_t"] is None else job["event_t"],
        "tmax": tmax,
        "status": status,
        "rc": r.returncode,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = build_jobs()
    with ProcessPoolExecutor(max_workers=8) as ex:
        rows = list(ex.map(run_job, jobs))
    with open(OUT / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['label']:16s} t={r['event_t']!s:5s} tmax={r['tmax']:<4} {r['status']}", flush=True)


if __name__ == "__main__":
    main()
