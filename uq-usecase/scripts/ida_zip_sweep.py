"""Fold location vs datacenter constant-impedance share alphaZ (hunting code).

Dispatch h1879_n26_dc148-200 (200 MW at bus 148). The datacenter ZIP is set to
alphaP = 1 - alphaZ, alphaI = 0. The pre-fault injection is unchanged for any alphaZ
because |V0| = Vn at the load bus, so f = aZ + aI + aP = 1 at the pre-fault point.
For each fault bus, runs the fixed-state continuation (ida_fixed_state_test.continuation)
and records whether a fault-on solution exists at X = 0.01 and X_crit if not.

alphaZ = 0 must reproduce the 21 folds of alpha_sweep.csv at alphaP = 1.0 (check).

Run with SLURM (ida_zip_sweep.sh), 16 workers. Output:
/scratch/isatkaus/gridkit-campaigns/ida_zoom/fixed_state/zip_sweep.csv
"""
import csv
import numpy as np
import json
import os
from multiprocessing import Pool
from pathlib import Path

import ida_fixed_state_test as t

ALPHA_Z = (0.0, 0.25, 0.5, 0.75)
OUT = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/fixed_state/zip_sweep.csv")


def build(az, bus):
    case = json.loads(t.CASE.read_text())
    net = t.Net(case, dc_on=True)
    k = net.dc_idx
    net.aZ[k] = az
    net.aI[k] = 0.0
    net.aP[k] = 1.0 - az
    net.set_fault_bus(net.num.index(bus))
    return net


def job(args):
    az, bus = args
    net = build(az, bus)
    s, _, _ = t.continuation(net)
    reached = s >= 1.0 - 1e-12
    return {
        "alpha_z": az,
        "fault_bus": bus,
        "reached": reached,
        "s_stop": round(s, 6),
        "x_crit": "" if reached else round(net.X / s, 6),
    }


def main():
    case = json.loads(t.CASE.read_text())
    buses = [b["number"] for b in case["buses"]]
    for z in ALPHA_Z:
        net = build(z, buses[0])
        res = float(np.max(np.abs(net.F(net.V0, 0.0))))
        print(f"alphaZ={z}: pre-fault max |F| = {res:.2e}")
    jobs = [(z, b) for z in ALPHA_Z for b in buses]
    fields = ["alpha_z", "fault_bus", "reached", "s_stop", "x_crit"]
    rows = []
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f, Pool(min(16, os.cpu_count() or 1)) as pool:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in pool.imap_unordered(job, jobs):
            rows.append(r)
            w.writerow(r)
            f.flush()
    for z in ALPHA_Z:
        sub = [r for r in rows if r["alpha_z"] == z]
        n_fold = sum(not r["reached"] for r in sub)
        xs = sorted(r["x_crit"] for r in sub if not r["reached"])
        line = f"alphaZ={z}: folds {n_fold} of {len(sub)}"
        if xs:
            line += f"; X_crit range {xs[0]} to {xs[-1]} (median {xs[len(xs)//2]})"
        print(line)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
