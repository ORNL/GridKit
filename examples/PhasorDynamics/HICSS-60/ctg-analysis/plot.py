#!/usr/bin/env python3
"""Render Figure 4 from saved contingency statistics; intentionally four cases."""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import plot as paper
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).resolve().parent
EFFORT_CASES = ("Hawaii", "IEEE39", "WECC240", "ACTIVSg200")


def plot_results(rows, path):
    keys = [(row["case"], int(row["bus"])) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate case/bus measurements")
    figure, axes = paper.subplots(2, 2, sharey=True)
    ratios = []
    for panel, name in enumerate(EFFORT_CASES):
        axis = axes.flat[panel]
        good = [row for row in rows if row["case"] == name and row["status"] == "ok"]
        steps = [int(row["accepted_steps"]) for row in good]
        jacobians = [int(row["jacobian_evals"]) for row in good]
        if any(value <= 0 for value in steps) or any(value < 0 for value in jacobians):
            raise ValueError(f"Invalid successful measurements for {name}")
        label, color = paper.CASE_STYLES[name]
        effort = [jac / step for jac, step in zip(jacobians, steps)]
        ratios.extend(effort)
        axis.scatter(steps, effort, s=28, color=color, alpha=0.8,
                     edgecolors="white", linewidths=0.3)
        axis.tick_params(which="both", top=True, right=True,
                         labelleft=False, labelright=panel % 2 == 1)
        axis.xaxis.set_major_locator(MaxNLocator(3, integer=True))
        axis.margins(x=0.08)
        if not good:
            axis.set_xticks([])
        paper.panel_label(axis, panel, label)
    if ratios:
        lower, upper = min(ratios), max(ratios)
        pad = 0.20 * (upper - lower) or max(0.01, 0.20 * upper)
        axes[0, 0].set_ylim(max(0.0, lower - pad), upper + pad)
    else:
        axes[0, 0].set_ylim(0.0, 1.0)
    axes[0, 0].yaxis.set_major_locator(MaxNLocator(4))
    paper.labels(figure, x="Accepted IDA steps",
                 y="Jacobian evaluation per accepted step")
    paper.save_figure(figure, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=HERE / "data/contingency_stats.csv")
    parser.add_argument("--figure", type=Path, default=HERE / "figures/contingency_effort.png")
    args = parser.parse_args()
    with args.data.open(newline="") as stream:
        plot_results(list(csv.DictReader(stream)), args.figure)


if __name__ == "__main__":
    main()
