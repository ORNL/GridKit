"""faultscreen-shaped micro-campaign on ACTIVSg2000.

    python test_2k_screen.py [n_faults]

One dispatch, many fault locations -- the stage-A/stage-B split that made the
200-bus campaigns cheap. Measures whether the lever still pays at 2000 buses
and what a contingency screen would cost there.

omega in the monitor output is a DEVIATION from synchronous speed, not absolute
pu, so a healthy machine sits near 0.0 and not near 1.0.
"""

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402
import m_to_case_utils as mu  # noqa: E402
import solver_utils as su  # noqa: E402

SPEC = mu.ACTIVSG2000
BASEMVA = 100.0
RUNNER = Path(
    "/home/isatkaus/gridkit/build/application/PhasorDynamics/DynamicSimulation"
)
W = Path(f"/scratch/{os.environ['USER']}/scale-probe/screen2k")


def rebuild():
    base = mu.load_case(SPEC.case_json_path)
    t0 = time.time()
    case, counts = mu.build_case_from_m(SPEC, base, SPEC.template_m, strict=False)
    t_build = time.time() - t0

    mp = mu.read_m(SPEC.template_m)
    g = mp.gen.copy()
    g["_gid"] = mu.gen_ids(g)
    rows = {r["_gid"]: r for _, r in g.iterrows() if float(r["GEN_STATUS"]) != 0}
    for gid in counts["gen_static"]:
        r = rows.get(gid)
        if r is None:
            continue
        case["devices"].append(
            {
                "class": "LoadZIP",
                "ports": {"bus": int(r["GEN_BUS"])},
                "id": f"{gid}_ibr_stand_in",
                "params": {
                    "Pnom": -float(r["PG"]) / BASEMVA,
                    "Qnom": -float(r["QG"]) / BASEMVA,
                    "alphaI": 1.0,
                    "alphaP": 0.0,
                },
            }
        )
    return case, counts, t_build


def run_fault(case, bus, tag, duration=0.1, tmax=5.0):
    wd = W / tag
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True)
    c = su.set_fault_bus(case, bus)
    c = su.enable_monitors(c, bus=(), genrou=("delta", "omega"))
    mu.write_case(c, wd / "case.json")
    (wd / "solver.json").write_text(
        json.dumps(
            {
                "system_model_file": "case.json",
                "output_file": "out.csv",
                "dt_monitor": 1 / 240,
                "tmax": tmax,
                "events": [
                    {"time": 1.0, "type": "fault_on", "element_id": 0},
                    {"time": 1.0 + duration, "type": "fault_off", "element_id": 0},
                ],
            }
        )
    )
    t0 = time.time()
    r = subprocess.run(
        [str(RUNNER), "solver.json"],
        cwd=wd,
        capture_output=True,
        text=True,
        env=cp.runner_env(),
    )
    dt = time.time() - t0
    (wd / "sim.log").write_text(r.stdout + r.stderr)
    csv = next(iter(wd.glob("*.csv")), None)

    rec = {
        "fault_bus": bus,
        "rc": r.returncode,
        "t_s": round(dt, 2),
        "MB": round(csv.stat().st_size / 1e6, 1) if csv else 0.0,
    }
    if r.returncode == 0 and csv is not None:
        d = pd.read_csv(csv)
        dl = d[[c for c in d.columns if c.endswith("_delta")]].to_numpy()
        om = d[[c for c in d.columns if c.endswith("_omega")]].to_numpy()
        spread = np.ptp(dl, axis=1).max() if dl.size else np.nan
        rec["max_delta_spread_deg"] = round(float(np.degrees(spread)), 1)
        rec["max_domega_pu"] = round(float(np.abs(om).max()), 5)
    return rec


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    W.mkdir(parents=True, exist_ok=True)

    case, counts, t_build = rebuild()
    print(
        f"rebuild: {t_build:.2f} s, {counts['gen']} machines patched, "
        f"{len(counts['gen_static'])} static units -> negative load"
    )

    # highest-degree buses: the ACTIVSg200 study found fault severity tracks
    # connectivity, so screen there first rather than uniformly
    mp = mu.read_m(SPEC.template_m)
    deg = {}
    for f, t in zip(mp.branch.F_BUS.astype(int), mp.branch.T_BUS.astype(int)):
        deg[f] = deg.get(f, 0) + 1
        deg[t] = deg.get(t, 0) + 1
    buses = [b for b, _ in sorted(deg.items(), key=lambda kv: -kv[1])[:n]]
    print(f"screening {len(buses)} highest-degree buses: {buses}\n")

    t0 = time.time()
    recs = [run_fault(case, b, f"f{b}") for b in buses]
    t_all = time.time() - t0

    df = pd.DataFrame(recs)
    print(df.to_string(index=False))
    ok = df[df.rc == 0]
    print(
        f"\n{len(ok)}/{len(df)} ran. fan-out {t_all:.1f} s total, "
        f"{t_all/len(df):.2f} s per fault"
    )
    print(
        f"rebuild {t_build:.2f} s is amortised over all faults: "
        f"stage-A share {100*t_build/(t_build+t_all):.1f}%"
    )
    print(f"monitor output {df.MB.sum():.0f} MB for {len(df)} faults")
    df.to_csv(W / "screen.csv", index=False)
    print(f"wrote {W/'screen.csv'}")


if __name__ == "__main__":
    main()
