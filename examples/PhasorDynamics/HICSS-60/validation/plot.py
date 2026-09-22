#!/usr/bin/env python3
"""Two-panel GridKit vs PowerWorld signal figure in the publication style."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import plot as paper
from common import signals
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.ticker import AutoMinorLocator

HERE = Path(__file__).resolve().parent
SIGNAL_COLOR = "#08306b"


def add_traces(ax, time, values):
    segments = [np.column_stack([time, values[:, index]]) for index in range(values.shape[1])]
    ax.add_collection(
        LineCollection(
            segments,
            colors=SIGNAL_COLOR,
            linewidths=0.28,
            alpha=0.10,
            capstyle="round",
            rasterized=True,
        )
    )
    return float(np.nanmin(values)), float(np.nanmax(values))


def limits(lo, hi, floor, ceil):
    span = hi - lo
    pad = 0.04 * span if span > 0.0 else max(1.0e-3, 0.04 * abs(hi))
    lower, upper = lo - pad, hi + pad
    if floor is not None:
        lower = min(floor, lower)
    if ceil is not None:
        upper = max(ceil, upper)
    return lower, upper


def configure_axis(ax, index, name, boundaries, window, ylim):
    ax.set_xlim(0.0, window)
    ax.set_ylim(*ylim)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    paper.event_markers(ax, boundaries)
    paper.panel_label(ax, index, name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=HERE / "data/ACTIVSg2000.omega.npz")
    parser.add_argument("--figure", type=Path, default=HERE / "figures/ACTIVSg2000.omega.signals.png")
    args = parser.parse_args()
    source = HERE / "solvers/ACTIVSg2000.solver.json"
    study = json.loads(source.read_text())
    reference = (source.parent / study["reference_file"]).resolve()
    case = json.loads((source.parent / study["system_model_file"]).read_text())
    frequency = float(case["params"]["freq_base"])
    with np.load(args.data) as data:
        time_out, values_out = data["time"], data["values"] * frequency
    _, ref = signals.load_csv(reference)
    time_ref, values_ref = ref[:, 0], ref[:, 1:] * frequency
    signals.check_grid(time_out, values_out, time_ref, values_ref)
    window = study["tmax"]
    boundaries = [event["time"] for event in study["events"]]
    figure, (ax_top, ax_bot) = paper.subplots(
        2, 1, sharex=True, sharey=True, figsize=(3.5, 2.4))
    lo_out, hi_out = add_traces(ax_top, time_out, values_out)
    lo_ref, hi_ref = add_traces(ax_bot, time_ref, values_ref)
    ylim = (-0.15, 0.15)
    if min(lo_out, lo_ref) < ylim[0] or max(hi_out, hi_ref) > ylim[1]:
        ylim = limits(min(lo_out, lo_ref), max(hi_out, hi_ref), *ylim)
    configure_axis(ax_top, 0, "GridKit", boundaries, window, ylim)
    configure_axis(ax_bot, 1, "PowerWorld", boundaries, window, ylim)
    paper.labels(figure, x=r"$t$ – Time [sec]", y=r"$\omega$ – Speed deviation [Hz]")
    paper.save_figure(figure, args.figure, tight=False)


if __name__ == "__main__":
    main()
