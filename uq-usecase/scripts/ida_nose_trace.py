"""Fault-on nose trace at alphaP 0.8 (ida_ic_issue.md section 9; hunting code).

Same fixed-state algebra and Newton solver as ida_fixed_state_test.py, but records the
fault-bus voltage at every accepted continuation step, so the upper branch can be plotted.
Two buses: 148 (folds at s = 0.367, X_crit = 0.027) and 3 (no fold, reaches s = 1).
Dispatch h1879_n26_dc148-200, datacenter at bus 148, 200 MW. Upper branch only: the lower
branch is not computed.

Run with SLURM (ida_nose_trace.sh). Output:
/scratch/isatkaus/gridkit-campaigns/ida_zoom/nose_trace/trace.csv
"""
import csv
import json
from pathlib import Path

import numpy as np

import ida_fixed_state_test as t

OUT = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/nose_trace/trace.csv")
ALPHA = 0.8
BUSES = (148, 3)


def traced_continuation(net):
    """Same step control as ida_fixed_state_test.continuation, recording each accepted step."""
    s, V = 0.0, net.V0.copy()
    ds = t.DS0
    path = [(0.0, abs(net.V0[net.fault_bus]))]
    while s < 1.0 - 1e-12:
        s_try = min(1.0, s + ds)
        Vt, ok = t.newton(net, V, s_try)
        if ok and np.max(np.abs(Vt - V)) < t.MAX_DV:
            s, V = s_try, Vt
            ds = min(ds * 1.5, t.DS_MAX)
            path.append((s, abs(V[net.fault_bus])))
        else:
            ds /= 2.0
            if ds < t.DS_MIN:
                break
    return s, path


def main():
    t.DC_ALPHA = ALPHA
    case = json.loads(t.CASE.read_text())
    rows = []
    for bus in BUSES:
        net = t.Net(case, dc_on=True)
        net.set_fault_bus(net.num.index(bus))
        s_end, path = traced_continuation(net)
        reached = s_end >= 1.0 - 1e-12
        print(f"bus {bus}: s_stop={s_end:.6f} reached={reached} steps={len(path) - 1}")
        for s, v in path:
            rows.append({
                "bus": bus,
                "s": round(s, 6),
                "x_eff": "" if s == 0.0 else round(net.X / s, 6),
                "v_fault": round(float(v), 6),
                "s_stop": round(s_end, 6),
                "reached": reached,
            })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
