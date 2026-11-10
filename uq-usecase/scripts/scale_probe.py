"""Scale probe: what does ACTIVSg2000 cost compared with ACTIVSg200?

    python scale_probe.py

Measures the three numbers a campaign projection needs -- one dynamic run,
monitor output size, and case load -- on both cases under identical settings,
so the ratios are meaningful rather than anecdotal. Read-only with respect to
the shipped cases; everything is written to a scratch directory.
"""

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402
import m_to_case_utils as mu  # noqa: E402
import solver_utils as su  # noqa: E402

CASES = {
    "ACTIVSg200": "/home/isatkaus/gridkit/cases/PhasorDynamics/ACTIVSg200/ACTIVSg200.case.json",
    "ACTIVSg2000": "/home/isatkaus/gridkit/cases/PhasorDynamics/ACTIVSg2000/ACTIVSg2000.case.json",
}
RUNNER = Path(
    "/home/isatkaus/gridkit/build/application/PhasorDynamics/DynamicSimulation"
)
SCRATCH = Path(f"/scratch/{os.environ['USER']}/scale-probe")


def run(case_path, name, monitors, tmax=10.0, duration=0.15):
    wd = SCRATCH / f"{name}_{'full' if monitors else 'omega'}"
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True)

    t0 = time.time()
    case = mu.load_case(case_path)
    t_load = time.time() - t0

    if monitors:
        case = su.enable_monitors(case, bus=("Vr", "Vi"), genrou=("delta", "omega"))
    eid = su.fault_element_id(case, None) if False else 0
    mu.write_case(case, wd / "case.json")

    solver = {
        "system_model_file": "case.json",
        "output_file": "out.csv",
        "dt_monitor": 1.0 / 240.0,
        "tmax": tmax,
        "events": [
            {"time": 1.0, "type": "fault_on", "element_id": eid},
            {"time": 1.0 + duration, "type": "fault_off", "element_id": eid},
        ],
    }
    (wd / "solver.json").write_text(json.dumps(solver, indent=2))

    t0 = time.time()
    r = subprocess.run(
        [str(RUNNER), "solver.json"],
        cwd=wd,
        capture_output=True,
        text=True,
        env=cp.runner_env(),
    )
    t_sim = time.time() - t0
    (wd / "sim.log").write_text(r.stdout + r.stderr)

    csvs = sorted(wd.glob("*.csv"))
    nbytes = sum(f.stat().st_size for f in csvs)
    ncols = 0
    if csvs:
        with open(csvs[0]) as fh:
            ncols = len(fh.readline().split(","))

    return {
        "case": name,
        "monitors": "full" if monitors else "default",
        "rc": r.returncode,
        "t_load_s": round(t_load, 3),
        "t_sim_s": round(t_sim, 2),
        "csv_MB": round(nbytes / 1e6, 2),
        "csv_cols": ncols,
        "case_json_MB": round(Path(case_path).stat().st_size / 1e6, 2),
        "n_bus": len(case["buses"]),
        "n_dev": len(case["devices"]),
    }


def main():
    SCRATCH.mkdir(parents=True, exist_ok=True)
    missing = [str(p) for p in [RUNNER] if not p.exists()]
    if missing:
        sys.exit(f"missing runner: {missing}")

    rows = []
    for name, path in CASES.items():
        for mon in (False, True):
            rows.append(run(path, name, mon))
            print(rows[-1], flush=True)

    print()
    hdr = [
        "case",
        "monitors",
        "rc",
        "n_bus",
        "n_dev",
        "case_json_MB",
        "t_load_s",
        "t_sim_s",
        "csv_cols",
        "csv_MB",
    ]
    print(" | ".join(h.rjust(12) for h in hdr))
    for r in rows:
        print(" | ".join(str(r[h]).rjust(12) for h in hdr))

    small = {r["monitors"]: r for r in rows if r["case"] == "ACTIVSg200"}
    big = {r["monitors"]: r for r in rows if r["case"] == "ACTIVSg2000"}
    print("\nratios (2000 / 200):")
    for k in ("full", "default"):
        if k in small and k in big and small[k]["t_sim_s"] > 0:
            print(
                f"  monitors={k:8s} sim {big[k]['t_sim_s']/small[k]['t_sim_s']:5.1f}x"
                f"   csv {big[k]['csv_MB']/max(small[k]['csv_MB'],1e-9):5.1f}x"
                f"   cols {big[k]['csv_cols']/max(small[k]['csv_cols'],1):5.1f}x"
            )


if __name__ == "__main__":
    main()
