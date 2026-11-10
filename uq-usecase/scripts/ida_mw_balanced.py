"""Balanced datacenter-MW sweep at bus 148 (ida_ic_issue.md §7 item 9, hunting code).

The earlier ic181_mw cases scaled only the datacenter Pnom; the generators were not re-solved,
so their pre-fault point is unbalanced (|F| 0.6 to 5 pu) except at 200 MW. This script:

  1. re-solves the no-fault AC power flow on the h1879 base (alphaP 0.8 campaign dispatch) at
     each MW, warm and flat start (pm_utils.run_pm_solve_batch, PowerModels.jl);
  2. writes a GridKit case.json per converged MW, from the reference case.json with bus
     Vr/Vi, Genrou p0/q0 and the datacenter LoadZIP Pnom/Qnom taken from the solved PF;
  3. checks the pre-fault balance |F| of each case, and checks that the 200 MW rebuild
     reproduces the reference case.json;
  4. runs the fixed-state fold scan (ida_fixed_state_test.continuation) for all buses at
     each converged MW, at alphaP = ALPHA_P.

Run with SLURM (ida_mw_balanced.sh). Output under
/scratch/isatkaus/gridkit-campaigns/ida_zoom/mw_balanced/:
  pf_status.csv, cases/mw{MW}.json, fold_mw.csv
"""
import csv
import json
import math
import os
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
import pm_utils  # noqa: E402
import pm_fault_check as pfc  # noqa: E402
import ida_fixed_state_test as t  # noqa: E402

DISPATCH = Path("/scratch/isatkaus/gridkit-campaigns/faultscreen1_ap08/h1879_n26_dc148-200")
BASE_M = DISPATCH / "ac.m"
REF_CASE = DISPATCH / "b1_d0.1" / "case.json"
OUT = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/mw_balanced")
DC_BUS = 148
REF_MW = 200.0
MW_GRID = [100, 150, 200, 250, 300, 350, 400]
ALPHA_P = 0.8
BALANCE_TOL = 1e-3


def read_rows(path, name):
    lines = Path(path).read_text().splitlines()
    s = next(i for i, l in enumerate(lines) if l.startswith(f"mpc.{name} = ["))
    rows = []
    for l in lines[s + 1:]:
        if l.strip().startswith("];"):
            break
        rows.append(l.rstrip().rstrip(";").split())
    return rows


def solve_pf(jb, dc, mw_list):
    """Warm and flat PF per MW. Returns {(start, mw): solved .m path or None}."""
    pfc.BASE_M = BASE_M
    fault = {"R": 0.0, "X": 1.0}
    solved = {}
    for start, flat in (("warm", False), ("flat", True)):
        wd = OUT / "pf" / start
        wd.mkdir(parents=True, exist_ok=True)
        cases, keys = [], []
        for mw in mw_list:
            inp = wd / f"mw{mw}.m"
            outp = wd / f"mw{mw}.solved.m"
            pfc.build_case(mw, None, dc, fault, inp)
            cases.append({"input": inp, "output": outp, "pf_type": "ac"})
            keys.append((start, mw, outp))
        recs, err, rc = pm_utils.run_pm_solve_batch(jb, pfc.PROJ, pfc.PROJ / "pm_solve.jl",
                                                    cases, wd, flat_start=flat)
        if rc != 0 or len(recs) != len(cases):
            raise RuntimeError(f"PF batch ({start}) failed rc={rc}\n{err[-2000:]}")
        for (st, mw, outp), r in zip(keys, recs):
            solved[(st, mw)] = outp if r.get("ok") else None
    return solved


def build_case_json(solved_m, mw, out_json):
    """Reference case.json with PF-solved bus voltages, gen P/Q and datacenter load."""
    case = json.loads(REF_CASE.read_text())
    buses = {int(r[0]): (float(r[7]), math.radians(float(r[8]))) for r in read_rows(solved_m, "bus")}
    # Fixed bus shunts (MATPOWER Bs, MVAr) are folded into the load Qnom as -Bs*V^2 in the
    # GridKit case (see m_to_case_utils: q_shunt), so Qnom must follow the solved V.
    bs = {int(r[0]): float(r[5]) for r in read_rows(solved_m, "bus")}
    v2_ref = {b["number"]: b["init"]["Vr"] ** 2 + b["init"]["Vi"] ** 2 for b in case["buses"]}
    gens = {}
    for r in read_rows(solved_m, "gen"):
        b = int(r[0])
        if b in gens:
            raise RuntimeError(f"two gen rows at bus {b}")
        gens[b] = (float(r[1]), float(r[2]))
    ref_dc = next(d for d in case["devices"]
                  if d["class"] == "LoadZIP" and d["id"].startswith("datacenter"))
    q_per_p = ref_dc["params"]["Qnom"] / (ref_dc["params"]["Pnom"])
    for b in case["buses"]:
        vm, va = buses[b["number"]]
        b["init"]["Vr"] = vm * math.cos(va)
        b["init"]["Vi"] = vm * math.sin(va)
    for d in case["devices"]:
        if d["class"] == "Genrou":
            pg, qg = gens[d["ports"]["bus"]]
            d["params"]["p0"] = pg / 100.0
            d["params"]["q0"] = qg / 100.0
        if d["class"] == "LoadZIP" and d["id"].startswith("load_"):
            b = d["ports"]["bus"]
            v2 = buses[b][0] ** 2
            d["params"]["Qnom"] -= bs[b] / 100.0 * (v2 - v2_ref[b])
        # IBR generators are negative constant-current LoadZIP devices (*_ibr_stand_in).
        if d["class"] == "LoadZIP" and d["id"].endswith("_ibr_stand_in"):
            pg, qg = gens[d["ports"]["bus"]]
            d["params"]["Pnom"] = -pg / 100.0
            d["params"]["Qnom"] = -qg / 100.0
        if d is ref_dc:
            d["params"]["Pnom"] = mw / 100.0
            d["params"]["Qnom"] = q_per_p * mw / 100.0
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(case, indent=4))
    return out_json


def check_reference(json_path):
    """Max differences between a rebuilt 200 MW case and the reference case.json."""
    ref = json.loads(REF_CASE.read_text())
    new = json.loads(Path(json_path).read_text())
    dv = max(max(abs(a["init"]["Vr"] - b["init"]["Vr"]), abs(a["init"]["Vi"] - b["init"]["Vi"]))
             for a, b in zip(ref["buses"], new["buses"]))
    dgen = 0.0
    for a, b in zip(ref["devices"], new["devices"]):
        if a["class"] == "Genrou":
            dgen = max(dgen, abs(a["params"]["p0"] - b["params"]["p0"]),
                       abs(a["params"]["q0"] - b["params"]["q0"]))
    return dv, dgen


def prefault_balance(json_path):
    t.DC_ALPHA = ALPHA_P
    t.CASE = Path(json_path)
    f0, _, _, _ = t.prefault_check("dc200")
    return float(f0)


def job(args):
    mw, bus, case_path = args
    t.DC_ALPHA = ALPHA_P
    case = json.loads(Path(case_path).read_text())
    net = t.Net(case, dc_on=True)
    net.set_fault_bus(net.num.index(bus))
    s, _, _ = t.continuation(net)
    reached = s >= 1.0 - 1e-12
    return {
        "mw": mw,
        "fault_bus": bus,
        "reached": reached,
        "s_stop": round(s, 6),
        "x_crit": "" if reached else round(net.X / s, 6),
    }


def main():
    jb = pm_utils.resolve_julia()
    if jb is None:
        raise RuntimeError("no Julia found; see pm_utils.JULIA_CANDIDATES")
    ref_dc = next(d for d in json.loads(REF_CASE.read_text())["devices"]
                  if d["class"] == "LoadZIP" and d["id"].startswith("datacenter"))
    dc = {"bus": DC_BUS, "pnom": REF_MW, "qnom": ref_dc["params"]["Qnom"] * 100.0}

    solved = solve_pf(jb, dc, MW_GRID)
    pf_rows, cases = [], {}
    for mw in MW_GRID:
        src = None
        for start in ("warm", "flat"):
            if solved[(start, mw)] is not None:
                src = (start, solved[(start, mw)])
                break
        pf_rows.append({"mw": mw,
                        "warm_ok": solved[("warm", mw)] is not None,
                        "flat_ok": solved[("flat", mw)] is not None,
                        "used": src[0] if src else ""})
        if src is None:
            print(f"mw{mw}: no PF solution (warm and flat); skipped")
            continue
        j = build_case_json(src[1], mw, OUT / "cases" / f"mw{mw}.json")
        f0 = prefault_balance(j)
        print(f"mw{mw}: PF from {src[0]} start, pre-fault |F| = {f0:.3e}")
        if f0 > BALANCE_TOL:
            raise RuntimeError(f"mw{mw} case is not balanced (|F| = {f0:.3e})")
        cases[mw] = j
        if mw == int(REF_MW):
            dv, dg = check_reference(j)
            print(f"  rebuilt 200 MW vs reference case.json: max |dVr/dVi| {dv:.3e}, "
                  f"max |dp0/dq0| {dg:.3e}")

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "pf_status.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mw", "warm_ok", "flat_ok", "used"])
        w.writeheader()
        w.writerows(pf_rows)

    buses = [b["number"] for b in json.loads(REF_CASE.read_text())["buses"]]
    jobs = [(mw, b, str(p)) for mw, p in cases.items() for b in buses]
    fields = ["mw", "fault_bus", "reached", "s_stop", "x_crit"]
    with (OUT / "fold_mw.csv").open("w", newline="") as f, \
            Pool(min(16, os.cpu_count() or 1)) as pool:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in pool.imap_unordered(job, jobs):
            w.writerow(r)
            f.flush()

    with (OUT / "fold_mw.csv").open() as f:
        rows = list(csv.DictReader(f))
    print(f"alphaP={ALPHA_P}: folds (no fault-on solution at X=0.01) per datacenter MW")
    for mw in MW_GRID:
        sub = [r for r in rows if int(r["mw"]) == mw]
        if sub:
            n_fold = sum(r["reached"] == "False" for r in sub)
            print(f"  {mw:>4} MW: {n_fold} of {len(sub)} buses folded")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
