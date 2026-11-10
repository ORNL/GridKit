"""Closed-mode IDA rerun for the two IDASolve buses (ida_ic_issue.md §7 item 8, hunting code).

Buses 116 and 131 pass the fixed-state test (a fault-on solution exists at X = 0.01)
but fail IDA in event mode at t = 1.0. This runs them in closed mode (fault on at t = 0)
at the same dispatch, alphaP 1.0, 200 MW datacenter, fault X = 0.01, and reports the
IDA status and the failing time step from sim.log.

    python ida_closed_check.py
Output: /scratch/isatkaus/gridkit-campaigns/ida_zoom/closed_check/
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
import ida_zoom as iz  # noqa: E402

DISPATCH = Path("/scratch/isatkaus/gridkit-campaigns/faultscreen1_ap08/h1879_n26_dc148-200")
OUT = iz.OUT_ROOT / "closed_check"
BUSES = [116, 131]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for fb in BUSES:
        job = dict(label=f"closed_b{fb}", src=str(DISPATCH / "b1_d0.1"), kind="closed",
                   fault_bus=fb, alpha_p=1.0, mode="closed", solver_mods={},
                   dc_zero=False, outdir=str(OUT))
        r = iz.run_one(job)
        log = (OUT / job["label"] / "sim.log").read_text()
        t_fail = re.findall(r"At t = ([0-9.eE+-]+) and h = ([0-9.eE+-]+)", log)
        print(f"bus {fb}: status={r['status']} t/h at failure={t_fail[-1] if t_fail else 'n/a'}")


if __name__ == "__main__":
    main()
