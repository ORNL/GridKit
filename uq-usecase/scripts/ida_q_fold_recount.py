"""Fold recount on the Q-limited dispatch at 200 and 400 MW (hunting code).

Takes the converged outer-loop PV->PQ solution from q_limit_pf.py (last iteration per MW),
rebuilds a balanced GridKit case.json with ida_mw_balanced.build_case_json, checks the
pre-fault balance, and runs the same fixed-state fold scan (alphaP 0.8, X = 0.01) for every bus.
Compares the per-bus fold sets with the unlimited-Q sweep (fold_mw.csv).

usage: sbatch ida_q_fold_recount.sh
"""
import csv
import glob
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
import ida_mw_balanced as M  # noqa: E402

QLIM_DIR = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/q_limit")
UNLIM_CSV = M.OUT / "fold_mw.csv"
OUT = QLIM_DIR / "fold_recount"
MW_LIST = [200, 400]


def latest_solved(mw):
    """Last outer-loop iteration's solved .m for one MW (highest it<NN>.solved.m)."""
    files = sorted(glob.glob(str(QLIM_DIR / f"mw{mw}" / "it*.solved.m")))
    if not files:
        raise RuntimeError(f"no solved outer-loop output for {mw} MW")
    return Path(files[-1])


def main():
    buses = [b["number"] for b in json.loads(M.REF_CASE.read_text())["buses"]]
    cases = {}
    for mw in MW_LIST:
        src = latest_solved(mw)
        j = M.build_case_json(src, mw, OUT / "cases" / f"mw{mw}.json")
        f0 = M.prefault_balance(j)
        print(f"mw{mw}: from {src.name}, pre-fault |F| = {f0:.3e}")
        if f0 > M.BALANCE_TOL:
            raise RuntimeError(f"mw{mw} Q-limited case is not balanced (|F| = {f0:.3e})")
        cases[mw] = j

    jobs = [(mw, b, str(p)) for mw, p in cases.items() for b in buses]
    fields = ["mw", "fault_bus", "reached", "s_stop", "x_crit"]
    OUT.mkdir(parents=True, exist_ok=True)
    out_csv = OUT / "fold_qlim_mw.csv"
    with out_csv.open("w", newline="") as f, Pool(min(16, os.cpu_count() or 1)) as pool:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in pool.imap_unordered(M.job, jobs):
            w.writerow(r)
            f.flush()

    with out_csv.open() as f:
        qrows = list(csv.DictReader(f))
    with UNLIM_CSV.open() as f:
        urows = list(csv.DictReader(f))

    print("per-MW fold counts (unlimited-Q -> Q-limited):")
    for mw in MW_LIST:
        uf = {int(r["fault_bus"]) for r in urows if int(r["mw"]) == mw and r["reached"] == "False"}
        qf = {int(r["fault_bus"]) for r in qrows if int(r["mw"]) == mw and r["reached"] == "False"}
        print(f"  {mw} MW: {len(uf)} -> {len(qf)} folded buses")
        print(f"    added by Q limits: {sorted(qf - uf)}")
        print(f"    removed by Q limits: {sorted(uf - qf)}")
    print(f"wrote {out_csv}")


if __name__ == "__main__":
    main()
