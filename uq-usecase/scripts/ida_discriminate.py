"""Discriminate what in an `open` vs `event` fault run changes the IC solve.

ida_zoom's `open` mode (state0 false, no events, tmax 0.2) passes at every bus,
while `event` mode (state0 false, fault_on at t = 1.0, tmax 1.2) fails at the
same buses as `closed`. The two differ in the event schedule and in tmax, so
they can not be separated from the existing modes. This script adds two modes,
in-process only, and leaves ida_zoom.py unchanged:

  open_ev_t02    state0 false, fault_on at 1.0, tmax 0.2  (event schedule, short run)
  open_noev_t12  state0 false, no events,        tmax 1.2  (long run, no event)

plus the existing `open`, `closed`, `event`, and a no-BusFault control (`none`).

usage (run as a SLURM job, 16 workers, about 1 h limit):
  python ida_discriminate.py --src <scratch case dir> --name <name> --buses 1,148,94,181 --alpha 1.0
"""

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")

import ida_zoom as iz  # noqa: E402

iz.MODES["open_ev_t02"] = (False, True, 0.2)
iz.MODES["open_noev_t12"] = (False, False, 1.2)
MODES_RUN = ["none", "open", "open_ev_t02", "open_noev_t12", "event", "closed"]


def build_jobs(src, outdir, buses, alpha):
    jobs = []
    for b in buses:
        for mode in MODES_RUN:
            fb = "none" if mode == "none" else b
            jobs.append(
                {
                    "label": f"{mode}_fb{b}_a{alpha:g}" if mode != "none" else f"none_a{alpha:g}",
                    "src": str(src),
                    "outdir": str(outdir),
                    "fault_bus": fb,
                    "alpha_p": alpha,
                    "mode": mode,
                    "solver_mods": {},
                    "fault_x": None,
                    "dc_mw": None,
                    "dc_zero": False,
                }
            )
    return jobs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--src", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--buses", required=True)
    p.add_argument("--alpha", type=float, default=1.0)
    p.add_argument("--workers", type=int, default=16)
    args = p.parse_args()

    src = Path(args.src)
    outdir = iz.OUT_ROOT / args.name
    outdir.mkdir(parents=True, exist_ok=True)
    buses = [int(x) for x in args.buses.split(",")]
    jobs = build_jobs(src, outdir, buses, args.alpha)

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(iz.run_one, jobs))

    with open(outdir / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['label']:28s} {r['status']}", flush=True)


if __name__ == "__main__":
    main()
