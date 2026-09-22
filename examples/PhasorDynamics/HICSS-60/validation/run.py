#!/usr/bin/env python3
"""Collect validation measurements for the table and Texas signals for Figure 7."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import plot as paper
from common import signals
from common.simulation import (ROOT, check_settings, monitor_signals, read_json, simulate,
                   simulation_result, write_json)
import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build/paper/application/PhasorDynamics/DynamicSimulation")
    parser.add_argument("--workdir", type=Path, default=ROOT / "build/HICSS-60/validation")
    parser.add_argument("--data", type=Path, default=HERE / "data/measurements.json")
    parser.add_argument("--cases", nargs="+", choices=tuple(paper.CASE_STYLES), default=tuple(paper.CASE_STYLES))
    parser.add_argument("--from-runs", action="store_true", help="Read existing signals.csv and validation.log without simulation")
    args = parser.parse_args()
    if args.data.exists() and not args.from_runs:
        raise FileExistsError(f"{args.data}: choose a new --data output")
    results = {}
    for name in args.cases:
        source = HERE / "solvers" / f"{name}.solver.json"
        config = read_json(source)
        reference = (source.parent / config["reference_file"]).resolve()
        work = (args.workdir / name).resolve()
        if args.from_runs:
            check_settings(read_json(work / "validation.solver.json"), config)
            counts, _ = simulation_result((work / "validation.log").read_text())
        else:
            case = read_json((source.parent / config["system_model_file"]).resolve())
            monitor_signals(case)
            write_json(work / "case.json", case)
            run = dict(config, system_model_file="case.json", output_file="signals.csv")
            run.pop("reference_file")
            counts, _ = simulate(args.binary.resolve(), work, run, "validation")
        columns, output = signals.load_csv(work / "signals.csv")
        errors = {}
        for signal, suffix in (("omega", "_omega"), ("vmag", "_Vm")):
            errors[signal], _, values = signals.compare(
                columns, output, reference.parent / f"{name}.{signal}.ref.csv", suffix)
            if errors[signal][1] >= config["error_tolerance"]:
                raise ValueError(f"{name} {signal}: relative max error exceeds validation tolerance")
            if name == "ACTIVSg2000" and signal == "omega":
                args.data.parent.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(args.data.parent / "ACTIVSg2000.omega.npz", time=output[:, 0], values=values)
        results[name] = {"counts": counts, "errors": errors, "status": "ok"}
    write_json(args.data, results)


if __name__ == "__main__":
    main()
