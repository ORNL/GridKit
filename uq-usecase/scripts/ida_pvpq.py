"""AC power flow with PV->PQ Q-limit switching at the IDA test points (ida_ic_issue.md, test 2).

Uses the fault-181 dispatch with the datacenter at bus 148 (the same case as
pm_fault_check.py). For fault buses 1, 148, 94, 181 (the last two are controls)
and no fault, sweeps the datacenter MW on the pm_fault_check grid, from a warm
and a flat start. Cases are built with pm_fault_check.build_case. The solver is
ida_pvpq.jl (PowerModels.jl, Ipopt), not pm_solve.jl.

Compares convergence with the IDA IC status from pm_fault_check.ida_lookup
(ida_zoom/ic181_mw, alphaP 1.0). Output:
  /scratch/isatkaus/gridkit-campaigns/ida_zoom/pvpq/results.csv

    module load julia
    python ida_pvpq.py
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")

import pm_fault_check as pfc  # noqa: E402
import pm_utils  # noqa: E402

JL = Path("/home/isatkaus/gridkit/uq-usecase/scripts/ida_pvpq.jl")
OUT = pfc.OUT.parent / "pvpq"
BUSES = [1, 148, 94, 181, None]


def build_manifest(dc, fault):
    manifest = []
    for start, flat in (("warm", False), ("flat", True)):
        wd = OUT / start
        wd.mkdir(parents=True, exist_ok=True)
        for fb in BUSES:
            for mw in pfc.MW_GRID:
                name = "none" if fb is None else str(fb)
                tag = f"{start}_b{name}_mw{mw}"
                inp = wd / f"b{name}_mw{mw}.m"
                pfc.build_case(mw, fb, dc, fault, inp)
                manifest.append({"tag": tag, "input": str(inp), "flat": flat,
                                 "start": start, "fault_bus": fb, "dc_mw": mw})
    return manifest


def main():
    dc, fault = pfc.dc_and_fault_params()
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest(dc, fault)
    man_path = OUT / "manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2))

    jb = pm_utils.resolve_julia()
    if jb is None:
        raise RuntimeError("no Julia found; see pm_utils.JULIA_CANDIDATES")
    res_path = OUT / "results.json"
    proc = subprocess.run(
        [str(jb), f"--project={pfc.PROJ}", str(JL), str(man_path), str(res_path)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr[-2000:])
    recs = json.loads(res_path.read_text())
    meta = {r["tag"]: m for r, m in zip(recs, manifest)}

    ida = pfc.ida_lookup()
    rows = []
    for r in recs:
        m = meta[r["tag"]]
        fb = m["fault_bus"]
        ida_status = "" if fb is None else ida.get((fb, float(m["dc_mw"]), 1.0), "")
        rows.append({
            "start": m["start"], "fault_bus": "none" if fb is None else fb,
            "dc_mw": m["dc_mw"], "plain_converged": r["plain_converged"],
            "plain_qviol": r["plain_qviol"], "converged": r["converged"],
            "status": r["status"], "n_switched": r["n_switched"],
            "iterations": r["iterations"], "final_qviol": r["final_qviol"],
            "vmin": r["vmin"], "vmax": r["vmax"], "ida_ic_alpha1": ida_status,
        })

    with open(OUT / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    for start in ("warm", "flat"):
        print(f"\n== {start} start: P = plain PF ok, S = PV-PQ switched ok, X = fail; "
              f"IDA status at alphaP 1 in last column ==")
        print(f"{'fault':>6} " + " ".join(f"{mw:>5}" for mw in pfc.MW_GRID))
        for fb in BUSES:
            name = "none" if fb is None else str(fb)
            cells = []
            for mw in pfc.MW_GRID:
                row = next(x for x in rows if x["start"] == start and str(x["fault_bus"]) == name
                           and x["dc_mw"] == mw)
                cells.append("  " + ("P" if row["plain_converged"] else
                                     ("S" if row["converged"] else "X")))
            print(f"{name:>6} " + " ".join(f"{c:>5}" for c in cells))
    print("\nIDA IC (alphaP 1.0; ok / IC = IDACalcIC / IS = IDASolve; none = no fault):")
    print(f"{'fault':>6} " + " ".join(f"{mw:>5}" for mw in pfc.MW_GRID))
    for fb in BUSES[:-1]:
        cells = []
        for mw in pfc.MW_GRID:
            st = ida.get((fb, float(mw), 1.0), "-")
            cells.append({"ok": "  ok", "IDACalcIC": "  IC", "IDASolve": "  IS"}.get(st, "   -"))
        print(f"{fb:>6} " + " ".join(f"{c:>5}" for c in cells))


if __name__ == "__main__":
    main()
