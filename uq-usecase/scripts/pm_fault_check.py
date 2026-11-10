"""AC power-flow check of the faulted, datacenter-sized steady state.

    python pm_fault_check.py

For each fault bus (plus a no-fault control) and each datacenter size, writes a
MATPOWER .m from the solved base case (ac.m) with:
  - the datacenter P/Q at bus 148 set to base load + MW (Q scaled with P, as
    the GridKit LoadZIP does), and
  - the BusFault shunt at the fault bus (G, B from R, X as in the BusFault
    README, converted to MVAr at V = 1).

Solves all cases with pm_utils.run_pm_solve_batch (PowerModels.jl, Ipopt), once
from the warm start in the .m and once flat start. Counts generator Q-limit
violations in each converged solution. Compares convergence with the IDA
initial-condition results in ida_zoom/ic181_mw/results.csv.

PowerModels solves constant-power loads, i.e. the alphaP = 1 limit. A ZIP load
with alphaP = 0.2 is not represented here.
"""

import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import pm_utils  # noqa: E402

CASE = Path("/scratch/isatkaus/gridkit-campaigns/faultscreen1_ap1/h6483_n21_dc148-200")
BASE_M = CASE / "ac.m"
CASE_JSON = CASE / "b181_d0.1" / "case.json"
OUT = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/pf_mw")
IDA_CSV = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/ic181_mw/results.csv")
PROJ = Path("/home/isatkaus/gridkit/uq-usecase/pm-solver")

FAULT_BUSES = [1, 66, 94, 108, 148, 181]  # None = no-fault control
MW_GRID = [25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600]


def dc_and_fault_params():
    """Datacenter (bus, Pnom, Qnom) and BusFault (R, X) from the scratch case.json."""
    case = json.loads(CASE_JSON.read_text())
    dc = fault = None
    for d in case["devices"]:
        if d["class"] == "LoadZIP" and d.get("id", "").startswith("datacenter_"):
            dc = {"bus": int(d["ports"]["bus"]), "pnom": d["params"]["Pnom"] * 100.0,
                  "qnom": d["params"]["Qnom"] * 100.0}
        if d["class"] == "BusFault":
            fault = {"R": d["params"]["R"], "X": d["params"]["X"]}
    if dc is None or fault is None:
        raise RuntimeError("case.json needs one datacenter LoadZIP and one BusFault")
    return dc, fault


def read_block(lines, name):
    """Return (start, end) line indices of the rows of mpc.<name> = [ ... ];"""
    start = next(i for i, l in enumerate(lines) if l.startswith(f"mpc.{name} = ["))
    end = next(i for i in range(start + 1, len(lines)) if lines[i].strip().startswith("];"))
    return start + 1, end


def build_case(mw, fault_bus, dc, fault, out_path):
    lines = BASE_M.read_text().splitlines()
    s, e = read_block(lines, "bus")
    g_pu = fault["R"] / (fault["R"] ** 2 + fault["X"] ** 2)
    b_pu = -fault["X"] / (fault["R"] ** 2 + fault["X"] ** 2)
    for i in range(s, e):
        tok = lines[i].rstrip().rstrip(";").split()
        bus = int(tok[0])
        if bus == dc["bus"]:
            # Base load = solved PD/QD minus the 200 MW datacenter (Pnom, Qnom in MW).
            pd_base = float(tok[2]) - dc["pnom"]
            qd_base = float(tok[3]) - dc["qnom"]
            tok[2] = f"{pd_base + mw:.10g}"
            tok[3] = f"{qd_base + mw * dc['qnom'] / dc['pnom']:.10g}"
        if fault_bus is not None and bus == fault_bus:
            tok[4] = f"{float(tok[4]) + g_pu * 100.0:.10g}"
            tok[5] = f"{float(tok[5]) + b_pu * 100.0:.10g}"
        lines[i] = "\t" + "\t".join(tok) + ";"
    # The base ac.m has ';;' at the end of gen rows; PowerModels' MATPOWER
    # parser rejects it ("inconsistent number of items in each row").
    out_path.write_text("\n".join(lines).replace(";;", ";") + "\n")


def qlim_stats(m_path):
    """(number of generator Q-limit violations, min bus voltage) of a solved .m."""
    lines = Path(m_path).read_text().splitlines()
    s, e = read_block(lines, "gen")
    nviol = 0
    for i in range(s, e):
        tok = lines[i].rstrip().rstrip(";").split()
        qg, qmax, qmin, status = float(tok[2]), float(tok[3]), float(tok[4]), float(tok[7])
        if status == 1 and (qg > qmax + 1e-3 or qg < qmin - 1e-3):
            nviol += 1
    s, e = read_block(lines, "bus")
    vmin = min(float(l.rstrip().rstrip(";").split()[7]) for l in lines[s:e])
    return nviol, vmin


def run_all():
    dc, fault = dc_and_fault_params()
    OUT.mkdir(parents=True, exist_ok=True)
    jb = pm_utils.resolve_julia()
    if jb is None:
        raise RuntimeError("no Julia found; see pm_utils.JULIA_CANDIDATES")

    fb_list = FAULT_BUSES + [None]
    cases_all = {}
    for start_name, flat in (("warm", False), ("flat", True)):
        wd = OUT / start_name
        wd.mkdir(parents=True, exist_ok=True)
        cases, keys = [], []
        for fb in fb_list:
            for mw in MW_GRID:
                tag = f"b{fb if fb is not None else 'none'}_mw{mw}"
                inp = wd / f"{tag}.m"
                outp = wd / f"{tag}.solved.m"
                build_case(mw, fb, dc, fault, inp)
                cases.append({"input": inp, "output": outp, "pf_type": "ac"})
                keys.append((fb, mw, outp))
        recs, err, rc = pm_utils.run_pm_solve_batch(jb, PROJ, PROJ / "pm_solve.jl",
                                                    cases, wd, flat_start=flat)
        if rc != 0 or len(recs) != len(cases):
            raise RuntimeError(f"batch failed rc={rc}\n{err[-2000:]}")
        for (fb, mw, outp), r in zip(keys, recs):
            ok = bool(r.get("ok"))
            nviol, vmin = qlim_stats(outp) if ok else (None, None)
            cases_all[(start_name, fb, mw)] = {
                "ok": ok, "status": r.get("status"), "solve_time": r.get("solve_time"),
                "qviol": nviol, "vmin": vmin,
            }
    return cases_all


def ida_lookup():
    """(fault_bus, dc_mw, alpha_p) -> 'ok' | 'IDACalcIC' | 'IDASolve' | ... from the event-mode MW sweep."""
    res = {}
    if not IDA_CSV.exists():
        return res
    pat = re.compile(r"mw([\d.]+)_fb(\d+)_a([\d.]+)")
    with IDA_CSV.open() as f:
        for row in csv.DictReader(f):
            m = pat.fullmatch(row["label"])
            if not m or row["mode"] != "event":
                continue
            key = (int(m.group(2)), float(m.group(1)), float(m.group(3)))
            res[key] = row["status"]
    return res


def main():
    results = run_all()
    OUT.joinpath("pf_results.csv").write_text(
        "start,fault_bus,dc_mw,ok,status,solve_time,qviol,vmin\n"
        + "\n".join(
            f"{s},{fb if fb is not None else 'none'},{mw},{r['ok']},{r['status']},"
            f"{r['solve_time']},{r['qviol']},{r['vmin']}"
            for (s, fb, mw), r in sorted(results.items(), key=lambda x: (x[0][0], str(x[0][1]), x[0][2]))
        )
        + "\n"
    )

    ida = ida_lookup()
    for start in ("warm", "flat"):
        print(f"\n== PowerModels AC PF, {start} start (ok = converged) ==")
        print(f"{'fault':>6} " + " ".join(f"{mw:>5}" for mw in MW_GRID))
        for fb in FAULT_BUSES + [None]:
            row = []
            for mw in MW_GRID:
                r = results[(start, fb, mw)]
                row.append("  ok" if r["ok"] else "   X")
            name = "none" if fb is None else str(fb)
            print(f"{name:>6} " + " ".join(f"{c:>5}" for c in row))

    print("\n== IDA event-mode IC status, alphaP 1.0 (ok / IDACalcIC / IDASolve) ==")
    print(f"{'fault':>6} " + " ".join(f"{mw:>5}" for mw in MW_GRID))
    for fb in FAULT_BUSES:
        cells = []
        for mw in MW_GRID:
            st = ida.get((fb, float(mw), 1.0))
            cells.append({"ok": "  ok", "IDACalcIC": "  IC", "IDASolve": "  IS"}.get(st, "   -"))
        print(f"{fb:>6} " + " ".join(f"{c:>5}" for c in cells))


if __name__ == "__main__":
    main()
