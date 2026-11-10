"""Buses 116 and 131 in closed mode under changed IDA settings (ida_ic_issue.md §7 item 8).

Baseline is the default solver (dt_monitor 1/240 s, rel_tol default). Each variant changes
one setting, to tell a solver-setting problem from a real fast dynamic instability:
finer output, tighter and looser rel_tol, lower max_order, tighter abs_tol.

    python ida_closed_variants.py
Output: /scratch/isatkaus/gridkit-campaigns/ida_zoom/closed_variants/
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
import ida_zoom as iz  # noqa: E402

DISPATCH = Path("/scratch/isatkaus/gridkit-campaigns/faultscreen1_ap08/h1879_n26_dc148-200")
OUT = iz.OUT_ROOT / "closed_variants"
BUSES = [116, 131]
VARIANTS = {
    "base": {},
    "fine_out": {"dt_monitor": 1e-4},
    "rel_tight": {"rel_tol": 1e-8},
    "rel_loose": {"rel_tol": 1e-4},
    "order2": {"max_order": 2},
    "abs_tight": {"abs_tol": 1e-10},
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, mods in VARIANTS.items():
        for fb in BUSES:
            label = f"{name}_b{fb}"
            job = dict(label=label, src=str(DISPATCH / "b1_d0.1"), kind="closed",
                       fault_bus=fb, alpha_p=1.0, mode="closed", solver_mods=mods,
                       dc_zero=False, outdir=str(OUT))
            r = iz.run_one(job)
            log = (OUT / label / "sim.log").read_text()
            t_fail = re.findall(r"At t = ([0-9.eE+-]+) and h = ([0-9.eE+-]+)", log)
            t_str = f"t={float(t_fail[-1][0]):.5f} h={float(t_fail[-1][1]):.1e}" if t_fail else "no failure log"
            rows = sum(1 for _ in (OUT / label / "ACTIVSg200.omega.csv").open()) - 1 \
                if (OUT / label / "ACTIVSg200.omega.csv").exists() else 0
            print(f"{name:>10} bus {fb}: status={r['status']:<10} {t_str} output_rows={rows}")


if __name__ == "__main__":
    main()
