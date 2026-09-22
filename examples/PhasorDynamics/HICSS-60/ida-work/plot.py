#!/usr/bin/env python3
"""Average per-fault cumulative IDA fractions and generate Figure 5."""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import plot as paper
from common import trace as ida
from common.simulation import ROOT, check_settings, read_json
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, PercentFormatter
import numpy as np


HERE = Path(__file__).resolve().parent
CASES = tuple(paper.CASE_STYLES)
METRICS = (
    ("accepted_steps", "Accepted IDA steps"),
    ("residual_evals", "Residual evaluations"),
    ("jacobian_evals", "Jacobian evaluations"),
    ("error_test_failures", "Error-test failures"),
)


def average(rows, workdir, cases, allow_partial=False):
    rows = [row for row in rows if row["case"] in cases]
    failed = [row for row in rows if row["status"] != "ok"]
    if failed and not allow_partial:
        raise ValueError("Incomplete contingency sweep; inspect the statistics before averaging")
    if failed:
        print(f"Excluding {len(failed)} failed runs from the average", flush=True)
    keys = [(row["case"], int(row["bus"])) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate case/bus measurements")
    for name in cases:
        study = read_json(HERE.parent / "ctg-analysis" / "solvers" / f"{name}.solver.json")
        ends = ida.boundaries(study)
        grid = np.unique(np.r_[np.linspace(0, study["tmax"], 1001), ends])
        sums = np.zeros((len(grid), len(METRICS)))
        contributing = np.zeros(len(METRICS), dtype=int)
        faults = [row for row in rows if row["case"] == name and row["status"] == "ok"]
        if not faults:
            raise ValueError(f"No measurements for {name}")
        for row in faults:
            path = workdir / name / f"bus-{row['bus']}.trace.csv"
            check_settings(read_json(workdir / name / f"bus-{row['bus']}.solver.json"), study)
            trace = ida.load(path, ends)
            times, counts = ida.cumulative_counts(trace, ends)
            expected = [int(row[key]) for key, _ in METRICS]
            if not np.array_equal(counts[-1], expected):
                raise ValueError(f"{path}: trace totals disagree with contingency statistics")
            sampled = counts[np.searchsorted(times, grid, side="right") - 1]
            nonzero = counts[-1] > 0
            sums[:, nonzero] += sampled[:, nonzero] / counts[-1, nonzero]
            contributing += nonzero
        means = np.divide(sums, contributing, out=np.full_like(sums, np.nan),
                          where=contributing > 0)
        path = HERE / "data" / f"{name}.mean.csv"
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(["t", *ida.COUNTERS,
                             *[f"{key}_faults" for key in ida.COUNTERS]])
            for t, values in zip(grid, means):
                writer.writerow([f"{t:.9g}", *[f"{value:.12g}" for value in values], *contributing])
        print(f"{name}: averaged {len(faults)} faults; contributing counts {contributing.tolist()}", flush=True)


def plot():
    figure, axes = paper.subplots(2, 2, sharex=True, sharey=True, legend=True)
    common = None
    cases = CASES
    for name in cases:
        _, color = paper.CASE_STYLES[name]
        source = HERE.parent / "ctg-analysis" / "solvers" / f"{name}.solver.json"
        study = read_json(source)
        ends = ida.boundaries(study)
        if common is not None and ends != common:
            raise ValueError("Cases must share duration and event times")
        common = ends
        data = np.genfromtxt(HERE / "data" / f"{name}.mean.csv", delimiter=",", names=True)
        for axis, (key, _) in zip(axes.flat, METRICS):
            axis.step(data["t"], 100 * data[key], where="post", color=color)
    for index, (axis, (_, label)) in enumerate(zip(axes.flat, METRICS)):
        axis.set_xlim(0, common[-1])
        axis.set_ylim(0, 100)
        axis.xaxis.set_major_locator(MultipleLocator(2))
        axis.yaxis.set_major_locator(MultipleLocator(25))
        axis.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
        paper.event_markers(axis, common[1:-1])
        paper.panel_label(axis, index, label, bottom=True)
    for value, label in zip(axes[0, 0].get_yticks(), axes[0, 0].get_yticklabels()):
        if value == 0:
            label.set_visible(False)
    handles = [Line2D([], [], color=color, label=label)
               for label, color in (paper.CASE_STYLES[name] for name in cases)]
    paper.legend(figure, handles)
    paper.labels(figure, x=r"$t$ – Time [sec]", y="Mean cumulative fraction [%]")
    paper.save_figure(figure, HERE / "figures/ida_work.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, default=ROOT / "build/HICSS-60/contingency")
    parser.add_argument("--stats", type=Path, help="Average existing contingency results without simulation")
    parser.add_argument("--cases", nargs="+", choices=tuple(paper.CASE_STYLES))
    parser.add_argument("--allow-partial", action="store_true", help="Explicitly exclude failed faults from the average")
    args = parser.parse_args()
    if (args.cases or args.allow_partial) and not args.stats:
        parser.error("--cases and --allow-partial require --stats")
    if args.stats:
        with args.stats.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        average(rows, args.workdir, args.cases or list(dict.fromkeys(row["case"] for row in rows)),
                allow_partial=args.allow_partial)
    plot()


if __name__ == "__main__":
    main()
