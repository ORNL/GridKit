"""Fold location vs alphaP (ida_ic_issue.md §7 item 8, hunting code).

Runs the fixed-state continuation (ida_fixed_state_test.continuation) for all 200 fault
buses at alphaP = 0.2, 0.8, 1.0 (dc200, dispatch h1879). For each bus records whether a
fault-on solution exists at X = 0.01 and, if not, X_crit = X / s_stop. Compares the
result with the event-mode IDA status at the same alphaP (mechanism/results.csv, 200 MW).

Run with SLURM (ida_fold_alpha.sh), 16 workers. Output:
/scratch/isatkaus/gridkit-campaigns/ida_zoom/fixed_state/alpha_sweep.csv
"""
import csv
import json
import os
from multiprocessing import Pool
from pathlib import Path

import ida_fixed_state_test as t

ALPHAS = (0.2, 0.8, 1.0)
OUT = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/fixed_state/alpha_sweep.csv")


def load_ida_by_alpha():
    """{(alpha, bus): IDA status} for event mode, 200 MW, default fault X."""
    out = {}
    with t.MECH.open() as f:
        for r in csv.DictReader(f):
            if r["mode"] != "event" or r["exp"] != "mw" or r["dc_mw"] != "200":
                continue
            if r["fault_x"]:
                continue
            out[(float(r["alpha_p"]), int(r["fault_bus"]))] = r["status"]
    return out


def job(args):
    alpha, bus = args
    t.DC_ALPHA = alpha
    case = json.loads(t.CASE.read_text())
    net = t.Net(case, dc_on=True)
    net.set_fault_bus(net.num.index(bus))
    s, _, _ = t.continuation(net)
    reached = s >= 1.0 - 1e-12
    return {
        "alpha_p": alpha,
        "fault_bus": bus,
        "reached": reached,
        "s_stop": round(s, 6),
        "x_crit": "" if reached else round(net.X / s, 6),
    }


def main():
    case = json.loads(t.CASE.read_text())
    buses = [b["number"] for b in case["buses"]]
    ida = load_ida_by_alpha()
    jobs = [(a, b) for a in ALPHAS for b in buses]
    fields = ["alpha_p", "fault_bus", "reached", "s_stop", "x_crit", "ida_status"]
    rows = []
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f, Pool(min(16, os.cpu_count() or 1)) as pool:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in pool.imap_unordered(job, jobs):
            r["ida_status"] = ida.get((r["alpha_p"], r["fault_bus"]), "")
            rows.append(r)
            w.writerow(r)
            f.flush()
    for a in ALPHAS:
        sub = [r for r in rows if r["alpha_p"] == a]
        n_fold = sum(not r["reached"] for r in sub)
        xs = sorted(r["x_crit"] for r in sub if not r["reached"])
        ida_n = sum(1 for r in sub if r["ida_status"])
        print(f"alphaP={a}: fold (no solution at X=0.01) {n_fold} of {len(sub)}; "
              f"IDA events with status {ida_n}")
        if xs:
            print(f"  X_crit range {xs[0]} to {xs[-1]} (median {xs[len(xs)//2]})")
        cross = {}
        for r in sub:
            if r["ida_status"]:
                k = (r["ida_status"], "fold" if not r["reached"] else "solution")
                cross[k] = cross.get(k, 0) + 1
        for k, v in sorted(cross.items()):
            print(f"  IDA={k[0]:<11} fixed-state={k[1]:<9} count={v}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
