"""Fold check for the fixed-state test (ida_ic_issue.md §7 item 8, hunting code).

For selected fault buses, records the continuation path in s (fault shunt scale) and
the relative smallest singular value of the fault-on Jacobian at each accepted step.
A true fold shows the relative sigma_min falling toward zero as s approaches the stop.
A step-size stop shows sigma_min still finite, and Newton converging from the stop
point at a tiny s increment.

    python ida_fold_check.py
Output: /scratch/isatkaus/gridkit-campaigns/ida_zoom/fixed_state/fold_check.csv
"""
import csv
import json
from pathlib import Path

import numpy as np

import ida_fixed_state_test as t

BUSES = [146, 1, 64, 134, 2, 116, 131, 94]
TINY = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
OUT = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/fixed_state/fold_check.csv")


def jac(net, V, s):
    n = net.n
    g0 = t._split(net.F(V, s))
    x = t._split(V)
    J = np.zeros((2 * n, 2 * n))
    for j in range(2 * n):
        xp = x.copy()
        xp[j] += t.H_FD
        J[:, j] = (t._split(net.F(t._join(xp, n), s)) - g0) / t.H_FD
    return J


def rel_sigma_min(net, V, s):
    sv = np.linalg.svd(jac(net, V, s), compute_uv=False)
    return sv[-1] / sv[0]


def path(net):
    """Same stepping as t.continuation, but records rel sigma_min at each accepted step."""
    s, V = 0.0, net.V0.copy()
    ds = t.DS0
    rec = []
    while s < 1.0 - 1e-12:
        s_try = min(1.0, s + ds)
        Vt, ok = t.newton(net, V, s_try)
        if ok and np.max(np.abs(Vt - V)) < t.MAX_DV:
            s, V = s_try, Vt
            ds = min(ds * 1.5, t.DS_MAX)
            rec.append((s, rel_sigma_min(net, V, s)))
        else:
            ds /= 2.0
            if ds < t.DS_MIN:
                break
    return s, V, rec


def main():
    case = json.loads(t.CASE.read_text())
    rows = []
    for bus in BUSES:
        net = t.Net(case, dc_on=True)
        net.set_fault_bus(net.num.index(bus))
        s, V, rec = path(net)
        reached = s >= 1.0 - 1e-12
        print(f"\nbus {bus}: reached={reached} s_stop={s:.6f}")
        print("  last accepted steps (s, rel sigma_min):")
        for ss, sg in rec[-6:]:
            print(f"    {ss:.6f}  {sg:.3e}")
        tiny = {}
        if not reached:
            for d in TINY:
                _, ok = t.newton(net, V, s + d, maxit=50)
                tiny[d] = ok
            print("  Newton from stop point at s+d converges:",
                  {d: ok for d, ok in tiny.items()})
            print(f"  rel sigma_min at stop: {rel_sigma_min(net, V, s):.3e}")
        rows.append({
            "fault_bus": bus, "reached": reached, "s_stop": round(s, 8),
            "rel_sigma_min_stop": f"{rel_sigma_min(net, V, s):.3e}",
            "rel_sigma_min_prev": f"{rec[-2][1]:.3e}" if len(rec) > 1 else "",
            "newton_ok_at_s_plus": " ".join(f"{d:g}:{int(ok)}" for d, ok in tiny.items()),
        })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
