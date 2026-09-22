#!/usr/bin/env python3
"""Run the five paper cases and plot accepted IDA step sizes over 20 seconds."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import plot as paper
from common import trace as ida
from common.simulation import ROOT, read_json


HERE = Path(__file__).resolve().parent


def collect(binary, workdir):
    for name, (label, _) in paper.CASE_STYLES.items():
        elapsed = ida.collect(HERE / "solvers" / f"{name}.solver.json", binary,
                              workdir / name)
        print(f"{label}: {elapsed:.6f} s CPU", flush=True)


def plot():
    figure, axis = paper.subplots(legend=True)
    common = None
    for name, (label, color) in paper.CASE_STYLES.items():
        path = HERE / "solvers" / f"{name}.solver.json"
        study = read_json(path)
        case = read_json((path.parent / study["system_model_file"]).resolve())
        events = [float(event["time"]) for event in study["events"]]
        frequency = float(case.get("params", {}).get("freq_base", 60.0))
        settings = (study["tmax"], events, frequency)
        if common is not None and common != settings:
            raise ValueError("Cases must share duration, event times, and nominal frequency")
        common = settings
        ends = ida.boundaries(study)
        trace = ida.load(path.parent / study["solver_trace_file"], ends)
        for segment, (times, sizes) in enumerate(ida.step_segments(trace, ends)):
            axis.step(times, sizes, where="pre", color=color,
                      label=label if segment == 0 else None)

    paper.step_axis(figure, axis, *common)
    paper.legend(figure)
    paper.save_figure(figure, HERE / "figures/adaptive_step.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path,
                        default=ROOT / "build/paper/application/PhasorDynamics/DynamicSimulation")
    parser.add_argument("--workdir", type=Path, default=ROOT / "build/HICSS-60/adaptive-step")
    parser.add_argument("--run", action="store_true", help="Collect new data before plotting")
    args = parser.parse_args()
    if args.run:
        collect(args.binary.resolve(strict=True), args.workdir)
    plot()


if __name__ == "__main__":
    main()
