"""Bisect which rebuild step breaks ACTIVSg2000 initialisation."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402
import m_to_case_utils as mu  # noqa: E402
import solver_utils as su  # noqa: E402

SPEC = mu.ACTIVSG2000
BASEMVA = 100.0
RUNNER = Path(
    "/home/isatkaus/gridkit/build/application/PhasorDynamics/DynamicSimulation"
)
W = Path(f"/scratch/{os.environ['USER']}/scale-probe/bisect2k")


def add_static(case, rows, ids):
    for gid in ids:
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


def run(case, tag):
    wd = W / tag
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True)
    case = su.enable_monitors(case, bus=(), genrou=("omega",))
    mu.write_case(case, wd / "case.json")
    (wd / "solver.json").write_text(
        json.dumps(
            {
                "system_model_file": "case.json",
                "output_file": "out.csv",
                "dt_monitor": 1 / 240,
                "tmax": 2.0,
                "events": [
                    {"time": 1.0, "type": "fault_on", "element_id": 0},
                    {"time": 1.1, "type": "fault_off", "element_id": 0},
                ],
            }
        )
    )
    r = subprocess.run(
        [str(RUNNER), "solver.json"],
        cwd=wd,
        capture_output=True,
        text=True,
        env=cp.runner_env(),
    )
    (wd / "sim.log").write_text(r.stdout + r.stderr)
    errs = [l for l in (r.stdout + r.stderr).splitlines() if "ERROR" in l]
    print(f"{tag:34s} rc={r.returncode}  {errs[0][:70] if errs else 'ok'}")
    return r.returncode


def main():
    W.mkdir(parents=True, exist_ok=True)
    base = mu.load_case(SPEC.case_json_path)
    mp = mu.read_m(SPEC.template_m)
    g = mp.gen.copy()
    g["_gid"] = mu.gen_ids(g)
    rows = {r["_gid"]: r for _, r in g.iterrows() if float(r["GEN_STATUS"]) != 0}

    run(mu.load_case(SPEC.case_json_path), "0_shipped")

    c = mu.load_case(SPEC.case_json_path)
    _, counts = mu.build_case_from_m(SPEC, base, SPEC.template_m, strict=False)
    static = counts["gen_static"]

    c = mu.load_case(SPEC.case_json_path)
    add_static(c, rows, static)
    run(c, "1_shipped+static_as_load")

    c, _ = mu.build_case_from_m(
        SPEC,
        base,
        SPEC.template_m,
        strict=False,
        remove_offline=False,
        patch_load=False,
        patch_gen=False,
    )
    run(c, "2_rebuild_busonly")

    c, _ = mu.build_case_from_m(
        SPEC,
        base,
        SPEC.template_m,
        strict=False,
        remove_offline=False,
        patch_load=False,
    )
    run(c, "3_rebuild_bus+gen")

    c, _ = mu.build_case_from_m(
        SPEC, base, SPEC.template_m, strict=False, remove_offline=False
    )
    run(c, "4_rebuild_bus+gen+load")

    c, _ = mu.build_case_from_m(SPEC, base, SPEC.template_m, strict=False)
    run(c, "5_rebuild_all(remove_offline)")

    c, _ = mu.build_case_from_m(SPEC, base, SPEC.template_m, strict=False)
    add_static(c, rows, static)
    run(c, "6_rebuild_all+static_as_load")


if __name__ == "__main__":
    main()
