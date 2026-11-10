"""Hunt matrix: fold sets vs datacenter MW, Q limits, load model, and dispatch (hunting code).

For one dispatch, builds balanced GridKit cases per MW in two Q modes and runs the fixed-state
fold scan (X = 0.01, all buses) for each load variant:
  unlim     PF as solved (no Q limits), as in ida_mw_balanced.py
  qlim      outer PV->PQ loop (q_limit_pf.run), last iteration's solution
  campaign  datacenter alphaP 0.8, alphaI 0.2 (campaign model)
  zconst    datacenter constant impedance (alphaZ 1)

usage: ida_hunt_matrix.py --dispatch h1879 --mw 100 150 --modes qlim --loads campaign --tag NAME
"""
import argparse
import csv
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
import ida_mw_balanced as M  # noqa: E402  (also sets up py-utils on sys.path)
import ida_fixed_state_test as t  # noqa: E402
import pm_utils  # noqa: E402
import q_limit_pf as Q  # noqa: E402

ROOT = Path("/scratch/isatkaus/gridkit-campaigns")
DISPATCHES = {
    "h1879": ROOT / "faultscreen1_ap08" / "h1879_n26_dc148-200",
    "h2151": ROOT / "faultscreen1_ap08" / "h2151_n27_dc148-200",
    "h4770": ROOT / "faultscreen1_ap08" / "h4770_n30_dc148-200",
}
OUT_ROOT = ROOT / "ida_zoom" / "hunt"
DC_BUS = 148
REF_MW = 200.0
ALPHA_P = 0.8
FIELDS = ["dispatch", "mw", "mode", "load", "fault_bus", "reached", "s_stop", "x_crit"]


def set_dispatch(name):
    """Point the balanced-sweep module at one dispatch and its own output directory."""
    d = DISPATCHES[name]
    M.DISPATCH = d
    M.BASE_M = d / "ac.m"
    M.REF_CASE = d / "b1_d0.1" / "case.json"
    M.OUT = OUT_ROOT / name
    return M.OUT


def ref_dc():
    return next(d for d in json.loads(M.REF_CASE.read_text())["devices"]
                if d["class"] == "LoadZIP" and d["id"].startswith("datacenter"))


def make_cases(jb, mws, modes):
    """Balanced case.json per (mode, MW). Returns {(mode, mw): path}."""
    cases = {}
    # Both modes start from the warm-start PF solve; the Q-limit loop reads its .m file.
    dc = {"bus": DC_BUS, "pnom": REF_MW, "qnom": ref_dc()["params"]["Qnom"] * 100.0}
    solved = M.solve_pf(jb, dc, mws)
    for mw in mws:
        if "unlim" in modes:
            src = solved[("warm", mw)] or solved[("flat", mw)]
            if src is None:
                raise RuntimeError(f"{mw} MW: no unlimited PF solution")
            j = M.build_case_json(src, mw, M.OUT / "cases" / f"unlim_mw{mw}.json")
            cases[("unlim", mw)] = j
            f0 = M.prefault_balance(j)
            print(f"unlim {mw} MW: pre-fault |F| = {f0:.3e}", flush=True)
            if f0 > M.BALANCE_TOL:
                raise RuntimeError(f"unlim {mw} MW not balanced")
        if "qlim" in modes:
            qdir = M.OUT / "q_limit" / f"mw{mw}"
            res = Q.run(M.OUT / "pf" / "warm" / f"mw{mw}.m", qdir, log=lambda *a: None)
            if not res["converged"]:
                raise RuntimeError(f"{mw} MW: Q-limit loop did not converge")
            last = sorted(qdir.glob("it*.solved.m"))[-1]
            j = M.build_case_json(last, mw, M.OUT / "cases" / f"qlim_mw{mw}.json")
            cases[("qlim", mw)] = j
            f0 = M.prefault_balance(j)
            print(f"qlim {mw} MW: {len(res['state'])} fixed units, pre-fault |F| = {f0:.3e}",
                  flush=True)
            if f0 > M.BALANCE_TOL:
                raise RuntimeError(f"qlim {mw} MW not balanced")
    return cases


def job(args):
    dispatch, mode, mw, load, bus, case_path = args
    t.DC_ALPHA = ALPHA_P
    case = json.loads(Path(case_path).read_text())
    net = t.Net(case, dc_on=True)
    if load == "zconst":
        k = net.dc_idx
        net.aZ[k], net.aI[k], net.aP[k] = 1.0, 0.0, 0.0
    net.set_fault_bus(net.num.index(bus))
    s, _, _ = t.continuation(net)
    reached = s >= 1.0 - 1e-12
    return {
        "dispatch": dispatch,
        "mw": mw,
        "mode": mode,
        "load": load,
        "fault_bus": bus,
        "reached": reached,
        "s_stop": round(s, 6),
        "x_crit": "" if reached else round(net.X / s, 6),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dispatch", choices=sorted(DISPATCHES), default="h1879")
    ap.add_argument("--mw", type=int, nargs="+", required=True)
    ap.add_argument("--modes", nargs="+", choices=["unlim", "qlim"], required=True)
    ap.add_argument("--loads", nargs="+", choices=["campaign", "zconst"], default=["campaign"])
    ap.add_argument("--buses", type=int, nargs="+", default=None,
                    help="restrict to these fault buses (smoke tests only)")
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()

    jb = pm_utils.resolve_julia()
    if jb is None:
        raise RuntimeError("no Julia found")
    out = set_dispatch(a.dispatch)
    out.mkdir(parents=True, exist_ok=True)
    cases = make_cases(jb, a.mw, a.modes)

    buses = a.buses or [b["number"] for b in json.loads(M.REF_CASE.read_text())["buses"]]
    jobs = [(a.dispatch, mode, mw, load, b, str(path))
            for (mode, mw), path in cases.items() for load in a.loads for b in buses]
    csv_path = out / f"{a.tag}.csv"
    with csv_path.open("w", newline="") as f, Pool(min(16, os.cpu_count() or 1)) as pool:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for r in pool.imap_unordered(job, jobs):
            w.writerow(r)
            f.flush()

    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    print(f"{a.dispatch} folds per (mode, load, MW), {len(buses)} buses:")
    for mode in a.modes:
        for load in a.loads:
            for mw in a.mw:
                sub = [r for r in rows if r["mode"] == mode and r["load"] == load
                       and int(r["mw"]) == mw]
                n = sum(r["reached"] == "False" for r in sub)
                print(f"  {mode:>5} {load:>8} {mw:>4} MW: {n} of {len(sub)}")
    print(f"wrote {csv_path}")


if __name__ == "__main__":
    main()
