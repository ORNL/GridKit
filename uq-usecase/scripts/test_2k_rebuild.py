"""End-to-end check that the re-dispatch path works on ACTIVSg2000.

    python test_2k_rebuild.py

Rebuilds a case.json from the solved .m via the same code path a campaign
uses, converts the 103 units that carry no dynamic model into negative load,
and runs the dynamics. Proves the prerequisite for any 2000-bus campaign.
"""

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402
import m_to_case_utils as mu  # noqa: E402
import solver_utils as su  # noqa: E402

BASEMVA = 100.0
SPEC = mu.ACTIVSG2000
RUNNER = Path(
    "/home/isatkaus/gridkit/build/application/PhasorDynamics/DynamicSimulation"
)
W = Path(f"/scratch/{os.environ['USER']}/scale-probe/rebuild2k")


def static_to_negative_load(case, spec, gen_rows, static_ids):
    """Units online in the .m with no dynamic model become negative LoadZIP.

    Same stand-in used for wind on ACTIVSg200: alphaI=1, because a
    constant-power negative injection makes the algebraic system much harder.
    """
    by_bus = {}
    for d in case["devices"]:
        if d["class"] == "LoadZIP":
            by_bus.setdefault(d["ports"]["bus"], []).append(d)

    added = 0
    for gid in static_ids:
        r = gen_rows.get(gid)
        if r is None:
            continue
        bus = int(r["GEN_BUS"])
        case["devices"].append(
            {
                "class": "LoadZIP",
                "ports": {"bus": bus},
                "id": f"{gid}_ibr_stand_in",
                "params": {
                    "Pnom": -float(r["PG"]) / BASEMVA,
                    "Qnom": -float(r["QG"]) / BASEMVA,
                    "alphaI": 1.0,
                    "alphaP": 0.0,
                },
            }
        )
        added += 1
    return added


def main():
    if W.exists():
        shutil.rmtree(W)
    W.mkdir(parents=True)

    base = mu.load_case(SPEC.case_json_path)
    m = SPEC.template_m
    print(f"template .m : {m}")
    print(f"base case   : {len(base['buses'])} buses, {len(base['devices'])} devices")

    t0 = time.time()
    case, counts = mu.build_case_from_m(SPEC, base, m, strict=False)
    t_build = time.time() - t0
    print(f"\nbuild_case_from_m: {t_build:.2f} s")
    for k, v in counts.items():
        print(f"  {k:14s} {len(v) if isinstance(v, list) else v}")

    mp = mu.read_m(m)
    g = mp.gen.copy()
    g["_gid"] = mu.gen_ids(g)
    rows = {r["_gid"]: r for _, r in g.iterrows() if float(r["GEN_STATUS"]) != 0}
    n_ibr = static_to_negative_load(case, SPEC, rows, counts["gen_static"])
    print(f"\nstatic units -> negative load: {n_ibr}")

    import collections

    print(
        "device mix after rebuild:",
        dict(collections.Counter(d["class"] for d in case["devices"])),
    )

    # power balance: dynamic load must equal dispatched generation minus losses
    pl = (
        sum(d["params"]["Pnom"] for d in case["devices"] if d["class"] == "LoadZIP")
        * BASEMVA
    )
    pg = (
        sum(
            d["params"]["p0"]
            for d in case["devices"]
            if d["class"] in SPEC.machine_classes
        )
        * BASEMVA
    )
    print(
        f"net dynamic load {pl:9.1f} MW   machine p0 {pg:9.1f} MW   gap {pg-pl:7.1f} MW"
    )

    case = su.enable_monitors(case, bus=(), genrou=("delta", "omega"))
    mu.write_case(case, W / "case.json")
    solver = {
        "system_model_file": "case.json",
        "output_file": "out.csv",
        "dt_monitor": 1 / 240,
        "tmax": 5.0,
        "events": [
            {"time": 1.0, "type": "fault_on", "element_id": 0},
            {"time": 1.1, "type": "fault_off", "element_id": 0},
        ],
    }
    (W / "solver.json").write_text(json.dumps(solver, indent=2))

    t0 = time.time()
    r = subprocess.run(
        [str(RUNNER), "solver.json"],
        cwd=W,
        capture_output=True,
        text=True,
        env=cp.runner_env(),
    )
    t_sim = time.time() - t0
    (W / "sim.log").write_text(r.stdout + r.stderr)
    csv = next(iter(W.glob("*.csv")), None)
    print(
        f"\ndynamics rc={r.returncode}  {t_sim:.2f} s  csv={csv.name if csv else None}"
    )
    if r.returncode != 0:
        print("--- log tail ---")
        print("\n".join((r.stdout + r.stderr).splitlines()[-15:]))
        sys.exit(1)

    import pandas as pd

    df = pd.read_csv(csv)
    om = [c for c in df.columns if "omega" in c.lower()]
    print(f"monitored machines: {len(om)}   rows: {len(df)}")
    print(f"max |omega-1| = {(df[om] - 1.0).abs().to_numpy().max():.4e} pu")
    print("\nREBUILD PATH OK")


if __name__ == "__main__":
    main()
