"""IDA seed test (ida_ic_issue.md §7 item 7): does a PF-solved fault-on bus state
let IDACalcIC converge where the default guess fails?

The IDA guess at t = 0 is the bus `init` (Vr, Vi) in case.json (Bus::initialize).
For each fault bus we solve the fault-on AC PF with PowerModels (pm_solve.jl, via
pm_fault_check.build_case), write the solved bus voltages into a copy of case.json,
and run `closed` mode (fault on at t = 0, tmax 0.2) once with the original guess and
once with the seed. Datacenter 0 MW (dc_zero) and 200 MW are both tested.

Controls: buses 1 and 94 (ok at dc_zero in ida_mechanism.py).

    python ida_seed_test.py
Output: /scratch/isatkaus/gridkit-campaigns/ida_zoom/seed_test/results.csv
"""

import csv
import json
import math
import shutil
import sys
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import ida_zoom as iz  # noqa: E402
import pm_fault_check as pfc  # noqa: E402
import pm_utils  # noqa: E402

DISPATCH = Path("/scratch/isatkaus/gridkit-campaigns/faultscreen1_ap08/h1879_n26_dc148-200")
OUT = iz.OUT_ROOT / "seed_test"
FAULT_BUSES = [64, 146, 1, 94]
MW_CASES = {"dc0": 0, "dc200": 200}


def dc_and_fault():
    case = json.loads((DISPATCH / "b1_d0.1" / "case.json").read_text())
    dc = fault = None
    for d in case["devices"]:
        if d["class"] == "LoadZIP" and d.get("id", "").startswith("datacenter_"):
            dc = {"bus": int(d["ports"]["bus"]), "pnom": d["params"]["Pnom"] * 100.0,
                  "qnom": d["params"]["Qnom"] * 100.0}
        if d["class"] == "BusFault":
            fault = {"R": d["params"]["R"], "X": d["params"]["X"]}
    return dc, fault


def solved_voltages(m_path):
    """{bus: (Vr, Vi)} from a solved MATPOWER .m (Vm pu, Va deg)."""
    lines = Path(m_path).read_text().splitlines()
    s, e = pfc.read_block(lines, "bus")
    out = {}
    for l in lines[s:e]:
        tok = l.rstrip().rstrip(";").split()
        vm, va = float(tok[7]), math.radians(float(tok[8]))
        out[int(tok[0])] = (vm * math.cos(va), vm * math.sin(va))
    return out


def seeded_src(label, volts):
    """Copy the dispatch case dir with bus init replaced by `volts`. Returns the src dir."""
    src = OUT / "_src" / label
    src.mkdir(parents=True, exist_ok=True)
    case = json.loads((DISPATCH / "b1_d0.1" / "case.json").read_text())
    for b in case["buses"]:
        vr, vi = volts[b["number"]]
        b["init"] = {"Vr": vr, "Vi": vi}
    (src / "case.json").write_text(json.dumps(case))
    return src


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pfc.BASE_M = DISPATCH / "ac.m"
    dc, fault = dc_and_fault()
    jb = pm_utils.resolve_julia()
    if jb is None:
        raise RuntimeError("no Julia found")

    wd = OUT / "pf"
    wd.mkdir(parents=True, exist_ok=True)
    cases, keys = [], []
    for name, mw in MW_CASES.items():
        for fb in FAULT_BUSES:
            inp, outp = wd / f"{name}_b{fb}.m", wd / f"{name}_b{fb}.solved.m"
            pfc.build_case(mw, fb, dc, fault, inp)
            cases.append({"input": inp, "output": outp, "pf_type": "ac"})
            keys.append((name, fb))
    recs, err, rc = pm_utils.run_pm_solve_batch(jb, pfc.PROJ, pfc.PROJ / "pm_solve.jl",
                                                cases, wd, flat_start=False)
    if rc != 0 or len(recs) != len(cases):
        raise RuntimeError(f"PF batch failed rc={rc}\n{err[-2000:]}")
    pf_ok = {k: bool(r.get("ok")) for k, r in zip(keys, recs)}
    print("PF fault-on converged:", pf_ok)

    jobs = []
    for name, mw in MW_CASES.items():
        dc_zero = mw == 0
        for fb in FAULT_BUSES:
            base = dict(fault_bus=fb, alpha_p=1.0, mode="closed", solver_mods={},
                        dc_zero=dc_zero, outdir=str(OUT / "runs"))
            jobs.append(dict(base, label=f"orig_{name}_b{fb}", src=str(DISPATCH / "b1_d0.1"),
                             kind="orig", dc=name))
            if not pf_ok[(name, fb)]:
                continue  # no PF seed: the seeded run is skipped
            volts = solved_voltages(wd / f"{name}_b{fb}.solved.m")
            seed = seeded_src(f"{name}_b{fb}", volts)
            jobs.append(dict(base, label=f"seed_{name}_b{fb}", src=str(seed),
                             kind="seed", dc=name))

    results = [iz.run_one(j) for j in jobs]
    with (OUT / "results.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["kind", "dc", "fault_bus", "status"])
        for j, r in zip(jobs, results):
            w.writerow([j["kind"], j["dc"], j["fault_bus"], r["status"]])
            print(f"{j['kind']:>5} {j['dc']:>6} bus {j['fault_bus']:>4}: {r['status']}")
    shutil.rmtree(OUT / "_src", ignore_errors=True)
    print(f"wrote {OUT / 'results.csv'}")


if __name__ == "__main__":
    main()
