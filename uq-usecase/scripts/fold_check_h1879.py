"""Compare sim-failed rows at dispatch h1879 with the fixed-state fold sets.

    python fold_check_h1879.py

Reads the alphaP sweep of fixed-state folds (ida_zoom/fixed_state/alpha_sweep.csv,
rows with ida_status IDACalcIC, 200 MW) and the fault sweep records. For the
h1879 rows it prints each sim-failed run (n_commit, dc_bus, fault_bus, kind, and
whether it also fails at alphaP 0.2 switched) and checks two things:

  - every sim-failed run that is sim-failed only at alphaP 0.8 has its fault bus in
    the alphaP 0.8 fold set;
  - no ok datacenter run (dc_mw > 0) at alphaP 0.8 has its fault bus in that fold set
    (a falsifier of H1). Baseline runs (dc_mw = 0) are excluded, since H1 is about the
    datacenter.

The sweep CSV does not record its datacenter bus, so this assumes the sweep used the
same datacenter placement as the rows. That is not checked here.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import failure_audit as fa  # noqa: E402

SWEEP = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/fixed_state/alpha_sweep.csv")
AP08 = "faultscreen1_ap08_qlim"
AP02 = "faultscreen1_qlim"
HOUR = 1879


def fold_set(alpha):
    d = pd.read_csv(SWEEP)
    s = d[(d.alpha_p == alpha) & (d.ida_status == "IDACalcIC")]
    return sorted(int(b) for b in s.fault_bus)


def main():
    a = fa.load(AP08)
    b = fa.load(AP02)
    f08 = set(fold_set(0.8))
    print(f"h1879 fold sets (IDACalcIC, 200 MW): alphaP 0.2 {fold_set(0.2)}")
    print(f"                                      alphaP 0.8 {sorted(f08)}")

    only08_in_fold = 0
    only08_total = 0
    print("sim-failed runs at h1879:")
    sim = [t for t in a if a[t]["hour"] == HOUR and a[t]["status"] == "sim_failed"]
    for t in sorted(sim, key=lambda t: (a[t]["n_commit"], a[t]["dc_bus"], a[t]["fault_bus"])):
        r = a[t]
        kind, tt = fa.sim_kind(AP08, t)
        both = b[t]["status"] == "sim_failed"
        tag = "both" if both else "only0.8"
        extra = f" t={tt:.4f}" if tt else ""
        print(f"  n{r['n_commit']} dc{r['dc_bus']} fault{r['fault_bus']} {kind} {tag}{extra}")
        if not both:
            only08_total += 1
            only08_in_fold += r["fault_bus"] in f08

    ok08_in_fold = sorted({
        r["fault_bus"] for r in a.values()
        if r["hour"] == HOUR and r["dc_mw"] > 0 and r["status"] == "ok" and r["fault_bus"] in f08
    })
    print(f"alphaP 0.8 only-failing runs with fault in the 0.8 fold set: "
          f"{only08_in_fold} of {only08_total}")
    print(f"ok datacenter runs at alphaP 0.8 with fault in the 0.8 fold set (falsifiers): {ok08_in_fold}")


if __name__ == "__main__":
    main()
