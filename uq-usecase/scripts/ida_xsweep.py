"""Fault-strength sweep for the IDA IC failure (ida_ic_issue.md, test 3).

Closed BusFault (state0 true, no events) at buses 1, 148 and a passing control
(94). alphaP = 1.0, datacenter 200 MW at bus 148, dispatch of fault 181.
Only the fault reactance X changes, from 1.0 pu (near open) down to 0.01 pu.

Reads the case through ida_zoom (run_one, set_fault, set_alpha) and leaves
ida_zoom.py unchanged. Output:
  /scratch/isatkaus/gridkit-campaigns/ida_zoom/<name>/results.csv  (adds fault_x)

    python ida_xsweep.py --name xsweep181 --workers 16
"""

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")

import ida_zoom as iz  # noqa: E402

SRC = Path("/scratch/isatkaus/gridkit-campaigns/faultscreen1_ap1/h6483_n21_dc148-200/b181_d0.1")
BUSES = [1, 148, 94]
XS = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]


def build_jobs(outdir):
    jobs = []
    for b in BUSES:
        for x in XS:
            jobs.append(
                {
                    "label": f"closed_fb{b}_x{x:g}_a1",
                    "src": str(SRC),
                    "outdir": str(outdir),
                    "fault_bus": b,
                    "alpha_p": 1.0,
                    "mode": "closed",
                    "solver_mods": {},
                    "fault_x": x,
                    "dc_mw": None,
                    "dc_zero": False,
                }
            )
    return jobs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--name", required=True)
    p.add_argument("--workers", type=int, default=16)
    args = p.parse_args()

    outdir = iz.OUT_ROOT / args.name
    outdir.mkdir(parents=True, exist_ok=True)
    jobs = build_jobs(outdir)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(iz.run_one, jobs))
    for row, job in zip(rows, jobs):
        row["fault_x"] = job["fault_x"]

    with open(outdir / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['label']:24s} {r['status']}", flush=True)


if __name__ == "__main__":
    main()
