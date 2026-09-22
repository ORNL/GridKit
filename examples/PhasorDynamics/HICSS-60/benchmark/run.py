#!/usr/bin/env python3
"""Collect monitor-free CPU timings, or reuse the matching runtime-sweep trials."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import plot as paper
from common.simulation import ROOT, disable_monitors, measure, pin_cpu, read_json, write_json

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build/paper/application/PhasorDynamics/DynamicSimulation")
    parser.add_argument("--workdir", type=Path, default=ROOT / "build/HICSS-60/benchmark")
    parser.add_argument("--data", type=Path, default=HERE / "data/measurements.json")
    parser.add_argument("--cases", nargs="+", choices=tuple(paper.CASE_STYLES), default=tuple(paper.CASE_STYLES))
    parser.add_argument("--cpu", type=int)
    reuse = parser.add_mutually_exclusive_group()
    reuse.add_argument("--from-runs", action="store_true", help="Read existing benchmark-N.log files without simulation")
    reuse.add_argument("--from-sweep", type=Path, help="Read nominal-mu trials from this runtime-sweep work directory")
    args = parser.parse_args()
    if args.data.exists() and not (args.from_runs or args.from_sweep):
        raise FileExistsError(f"{args.data}: choose a new --data output")
    if not (args.from_runs or args.from_sweep):
        print(f"CPU affinity: {pin_cpu(args.cpu)}", flush=True)
    results = {}
    for name in args.cases:
        source = HERE / "solvers" / f"{name}.solver.json"
        config = read_json(source)
        work = ((args.from_sweep or args.workdir) / name).resolve()
        case = read_json((source.parent / config["system_model_file"]).resolve())
        disable_monitors(case)
        if args.from_runs or args.from_sweep:
            if read_json(work / "case.json") != case:
                raise ValueError(f"{name}: recorded case differs from the benchmark case")
        else:
            write_json(work / "case.json", case)
        config["system_model_file"] = "case.json"
        label = f"mu-{config['mu']:g}" if args.from_sweep else "benchmark"
        results[name] = measure(args.binary.resolve(), work, config, label,
                                reuse=bool(args.from_runs or args.from_sweep))
    write_json(args.data, results)


if __name__ == "__main__":
    main()
