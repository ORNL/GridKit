"""IDA failure mechanism on the h1879 dispatch (ida_ic_issue.md §7, follow-up questions).

Three experiments on the `h1879_n26_dc148-200` dispatch at fault event mode
(fault closed at t = 1.0), reusing ida_zoom.run_one:

  nodc  : datacenter zeroed (dc_zero), all 200 fault buses. Tests whether the
          failures need the datacenter at all.
  mw    : datacenter MW sweep at alphaP 0.2, 0.8, 1.0 for the buses that fail at
          alphaP 1.0 (the 0.2 failures are a subset). Tests whether the 0.8 extra
          failures are a loadability threshold in MW.
  xsweep: fault reactance X sweep at alphaP 0.2 and 0.8 for the same buses. Tests
          whether alphaP moves the threshold in fault strength.

Output: /scratch/isatkaus/gridkit-campaigns/ida_zoom/mechanism/results.csv
    python ida_mechanism.py --workers 16
"""

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import ida_zoom as iz

SRC = Path("/scratch/isatkaus/gridkit-campaigns/faultscreen1_ap08/h1879_n26_dc148-200/b1_d0.1")
AP08 = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/ap08_bus_h1879/results.csv")
OUT = iz.OUT_ROOT / "mechanism"
MW_LEVELS = [25, 50, 100, 150, 200]
X_LEVELS = [0.5, 0.3, 0.2, 0.15, 0.1, 0.05, 0.02, 0.01]


def failing_buses_at(alpha):
    out = set()
    with AP08.open() as f:
        for row in csv.DictReader(f):
            if row["mode"] == "event" and row["alpha_p"] == str(alpha) and row["fault_bus"] != "none":
                if row["status"] != "ok":
                    out.add(int(row["fault_bus"]))
    return sorted(out)


def job(exp, label, fb, alpha, **extra):
    return dict(label=label, src=str(SRC), outdir=str(OUT), fault_bus=fb,
                alpha_p=alpha, mode="event", solver_mods={}, exp=exp, **extra)


def build_jobs():
    jobs = []
    for fb in range(1, 201):
        jobs.append(job("nodc", f"nodc_fb{fb}", fb, 1.0, dc_zero=True))
    buses = failing_buses_at(1.0)
    for a in (0.2, 0.8, 1.0):
        for mw in MW_LEVELS:
            for fb in buses:
                jobs.append(job("mw", f"mw{mw}_a{a}_fb{fb}", fb, a, dc_mw=mw))
    for a in (0.2, 0.8):
        for x in X_LEVELS:
            for fb in buses:
                jobs.append(job("xsweep", f"x{x}_a{a}_fb{fb}", fb, a, fault_x=x))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = build_jobs()
    print(f"{len(jobs)} runs; buses failing at alphaP 1.0: {failing_buses_at(1.0)}")
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        results = list(ex.map(iz.run_one, jobs))
    for j, r in zip(jobs, results):
        r["exp"] = j["exp"]
        r["dc_mw"] = j.get("dc_mw", "")
        r["dc_zero"] = bool(j.get("dc_zero", False))
        r["fault_x"] = j.get("fault_x", "")
    fields = ["exp", "label", "mode", "fault_bus", "alpha_p", "dc_mw", "dc_zero", "fault_x",
              "status", "rc", "seconds"]
    with (OUT / "results.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    print(f"wrote {OUT / 'results.csv'}")


if __name__ == "__main__":
    main()
