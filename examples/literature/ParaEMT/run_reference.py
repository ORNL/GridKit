"""Run a pinned ParaEMT 9-bus disturbance case; this does not run GridKit.

The stepping sequence follows upstream/main_step1_simulation.py (BSD-3-Clause).
Model equations and imported upstream source files are unchanged.
"""

import argparse
from importlib.metadata import version
import json
from pathlib import Path
import platform
import shutil
import sys
import tempfile
import time
import warnings
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "upstream"))

# xlrd 1.2.0 supports XLSX but calls this removed Python ElementTree alias.
if not hasattr(ET.ElementTree, "getiterator"):
    ET.ElementTree.getiterator = ET.ElementTree.iter

from psutils import initialize_emt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dt-us", type=float, default=50.0)
    parser.add_argument("--duration", type=float, default=3.0)
    parser.add_argument("--event", choices=("trip", "governor-step"), default="trip")
    parser.add_argument("--save-states", action="store_true",
                        help="also save all upstream state vectors as states.npz")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ts = args.dt_us / 1e6
    # Keep a common 0.5 ms output grid when refining the integration time step.
    dsrate = round(0.0005 / ts)
    nsteps = round(args.duration / ts)
    if (ts <= 0 or dsrate < 1 or args.duration <= 1.0
            or not np.isclose(dsrate * ts, 0.0005, rtol=0, atol=1e-12)
            or not np.isclose(nsteps * ts, args.duration, rtol=0, atol=1e-12)
            or nsteps % dsrate):
        parser.error("dt must divide 500 us and duration must exceed 1 s and end on the output grid")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    # NumPy warns about one upstream scalar cast on every post-trip step.
    warnings.filterwarnings("ignore", category=DeprecationWarning, module="Lib_BW")
    start = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="paraemt-9bus-") as temporary:
        work = Path(temporary)
        (work / "cases").mkdir()
        (work / "models/9bus_psse").mkdir(parents=True)
        shutil.copyfile(ROOT / "pfd_9_1_1.json", work / "cases/pfd_9_1_1.json")
        shutil.copyfile(ROOT / "9bus.xlsx", work / "models/9bus_psse/9bus.xlsx")
        pfd, ini, dyd, emt = initialize_emt(
            str(work), 2, 1, 1, ts, args.duration, mode="lu", nparts=2)

    assert (dyd.gen_genrou_n, dyd.exc_sexs_n, dyd.gov_gast_n, dyd.pss_ieeest_n) == (3, 3, 3, 3)
    assert len(pfd.ibr_bus) == 0 and len(pfd.bus_num) == 9
    # Same events and load/PLL settings as the original upstream driver.
    emt.t_sc, emt.i_gen_sc, emt.flag_exc_gov = 100, 0, 1
    emt.dsp, emt.flag_sc = -0.02, 1
    emt.t_gentrip, emt.i_gentrip = 1, 0
    emt.flag_gentrip, emt.flag_reinit = 1, 1
    emt.t_release_f, emt.loadmodel_option = 0.0, 1
    if args.event == "governor-step":
        emt.t_sc = 1.0
        emt.t_gentrip = 1000.0
    init_end = time.perf_counter()
    event_applied_at = None

    # Preserve upstream prediction history, save cadence, and reinitialization.
    for tn in range(1, nsteps + 1):
        previous_step_flag = emt.flag_sc
        emt.StepChange(dyd, ini, tn)
        if previous_step_flag == 1 and emt.flag_sc == 0:
            event_applied_at = tn * ts
        previous_trip_flag = emt.flag_gentrip
        emt.GenTrip(pfd, dyd, ini, tn, "lu")
        if previous_trip_flag == 1 and emt.flag_gentrip == 0:
            event_applied_at = tn * ts
        emt.predictX(pfd, dyd, emt.ts)
        emt.Igs = emt.Igs * 0
        emt.updateIg(pfd, dyd, ini)
        emt.Igi = emt.Igi * 0
        emt.Iibr = emt.Iibr * 0
        emt.updateIibr(pfd, dyd, ini)
        emt.solveV(ini)
        emt.BusMea(pfd, dyd, tn)
        emt.updateX(pfd, dyd, ini, tn)
        emt.updateXibr(pfd, dyd, ini, ts)
        emt.x_pred = {0: emt.x_pred[1], 1: emt.x_pred[2], 2: emt.x_pv_1}
        if tn % dsrate == 0:
            k = tn // dsrate
            emt.t.append(tn * ts)
            emt.x[k] = emt.x_pv_1.copy()
            emt.x_bus[k] = emt.x_bus_pv_1.copy()
            emt.x_load[k] = emt.x_load_pv_1.copy()
            emt.v[k] = emt.Vsol.copy()
            if not all(np.isfinite(value).all() for value in
                       (emt.x[k], emt.x_bus[k], emt.v[k])):
                raise RuntimeError(f"Nonfinite output at t={tn * ts}")
        if emt.flag_gentrip == 0 and emt.flag_reinit == 1:
            emt.Re_Init(pfd, dyd, ini)
        else:
            emt.updateIhis(ini)
        if tn % round(0.5 / ts) == 0:
            print(f"t={tn * ts:.4f} s", flush=True)
    stop = time.perf_counter()
    if event_applied_at is None:
        raise RuntimeError("Disturbance was not applied")

    t = np.array(emt.t)
    x = np.stack(list(emt.x.values()))
    v = np.stack(list(emt.v.values()))
    xb = np.stack(list(emt.x_bus.values()))
    data = {"time_s": t}
    schema = {"time_s": "Simulation time, seconds"}
    for phase_idx, phase in enumerate("abc"):
        for bus_idx, bus in enumerate(pfd.bus_num):
            key = f"bus{bus}_v{phase}_pu"
            data[key] = v[:, phase_idx * len(pfd.bus_num) + bus_idx]
            schema[key] = "Instantaneous phase voltage / (sqrt(2/3) * nominal line-line RMS voltage)"
    for i, bus in enumerate(pfd.bus_num):
        key = f"bus{bus}_vm_pu"
        data[key] = xb[:, i * dyd.bus_odr + 3]
        schema[key] = "ParaEMT instantaneous three-phase magnitude sqrt(2/3*(va^2+vb^2+vc^2)); not windowed RMS"
    for i, bus in enumerate(pfd.gen_bus):
        columns = {
            "speed_pu": (dyd.gen_genrou_xi_st + i * dyd.gen_genrou_odr + 1, pfd.ws,
                         "Rotor electrical speed / synchronous electrical speed"),
            "efd_pu": (dyd.exc_sexs_xi_st + i * dyd.exc_sexs_odr + 1, 1,
                       "SEXS field voltage on ParaEMT exciter base"),
            "pm_pu": (dyd.gov_gast_xi_st + i * dyd.gov_gast_odr + 3, 1,
                      "GAST mechanical power on machine MVA base"),
            "pss_vs_pu": (dyd.pss_ieeest_xi_st + i * dyd.pss_ieeest_odr + 9, 1,
                          "IEEEST stabilizer output into SEXS"),
        }
        for name, (index, scale, description) in columns.items():
            key = f"gen{bus}_{name}"
            data[key] = x[:, index] / scale
            schema[key] = description
    pd.DataFrame(data).to_csv(output / "reference.csv.gz", index=False, float_format="%.12e",
                             compression={"method": "gzip", "mtime": 0})
    # All saved state vectors, including controller and bus-measurement states.
    if args.save_states:
        np.savez_compressed(output / "states.npz", time_s=t, x=x, v=v, x_bus=xb,
                            x_load=np.stack(list(emt.x_load.values())))
    (output / "columns.json").write_text(json.dumps(schema, indent=2) + "\n")
    machine_parameters = {name: value.tolist() for name, value in vars(dyd).items()
                          if name.startswith("ec_") and isinstance(value, np.ndarray)}
    (output / "machine_parameters.json").write_text(json.dumps(machine_parameters, indent=2) + "\n")
    metadata = {
        "simulator": "ParaEMT (GridKit not run)",
        "upstream_revision": "d79d735a4a587d56c5b88187d1a499195b6b2b84",
        "systemN": 2, "duration_s": args.duration, "dt_s": ts,
        "output_interval_s": dsrate * ts, "steps": nsteps, "samples": len(t),
        "network_solver": "serial sparse LU", "loadmodel_option": 1,
        "models": {"GENROU": 3, "SEXS": 3, "GAST": 3, "IEEEST": 3},
        "generator_bus": pfd.gen_bus.tolist(), "generator_MVA_base": pfd.gen_MVA_base.tolist(),
        "bus_number": pfd.bus_num.tolist(), "bus_base_kV_LL_RMS": pfd.bus_basekV.tolist(),
        "nominal_frequency_Hz": pfd.ws / (2 * np.pi),
        "event": {"type": args.event, "generator_bus": 1, "time_s": 1.0,
                  "governor_reference_increment_pu": -0.02 if args.event == "governor-step" else None,
                  "applied_time_s": event_applied_at, "flag_reinit_at_end": int(emt.flag_reinit)},
        "python": platform.python_version(),
        "packages": {p: version(p) for p in ("numpy", "scipy", "numba", "llvmlite", "pandas", "xlrd", "matplotlib")},
        "initialization_wall_s": init_end - start, "loop_wall_s": stop - init_end,
        "all_saved_values_finite": True,
        "state_layout_source": "upstream/Lib_BW.py Initialize.CombineX",
    }
    (output / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
