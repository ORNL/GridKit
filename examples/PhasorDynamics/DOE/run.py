#!/usr/bin/env python3
"""Collect monitor-free GridKit CPU timings for the DOE runtime table."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "HICSS-60"))
from common.simulation import ROOT, disable_monitors, pin_cpu, read_json, simulate, write_json

HERE = Path(__file__).resolve().parent
SOLVERS = HERE / "solvers"


def main():
    names = sorted(path.name.removesuffix(".solver.json") for path in SOLVERS.glob("*.solver.json"))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build/application/PhasorDynamics/DynamicSimulation")
    parser.add_argument("--workdir", type=Path, default=ROOT / "build/DOE")
    parser.add_argument("--data", type=Path, default=HERE / "data/measurements.json")
    parser.add_argument("--cases", nargs="+", choices=names, default=names)
    parser.add_argument("--trials", type=int, default=1)
    parser.add_argument("--cpu", type=int)
    args = parser.parse_args()
    print(f"CPU affinity: {pin_cpu(args.cpu)}", flush=True)
    results = {}
    if args.data.exists():
        results = read_json(args.data)
    for name in args.cases:
        source = SOLVERS / f"{name}.solver.json"
        config = read_json(source)
        case = read_json((source.parent / config["system_model_file"]).resolve())
        disable_monitors(case)
        work = (args.workdir / name).resolve()
        write_json(work / "case.json", case)
        config["system_model_file"] = "case.json"
        trials = [simulate(args.binary.resolve(), work, config, f"doe-{trial}") for trial in range(1, args.trials + 1)]
        if any(counts != trials[0][0] for counts, _ in trials):
            raise ValueError(f"{name}: variable counts changed between timing trials")
        results[name] = {"counts": trials[0][0], "cpu_seconds": [seconds for _, seconds in trials]}
        write_json(args.data, results)
        print(f"{name}: {results[name]['cpu_seconds']}", flush=True)


if __name__ == "__main__":
    main()
