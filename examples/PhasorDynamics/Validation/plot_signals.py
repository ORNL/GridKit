#!/usr/bin/env python3
"""Two-panel GridKit vs PowerWorld signal figure in the publication style."""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.ticker import AutoMinorLocator

SIGNAL_COLOR = "#08306b"
SPINE_COLOR = "#222222"
DPI = 600

RC_PARAMS = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.labelsize": 8.5,
    "axes.linewidth": 0.6,
    "axes.edgecolor": SPINE_COLOR,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.7,
    "ytick.minor.size": 1.7,
    "lines.linewidth": 1.05,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
}

FIELDS = {
    "omega": {"label": r"$\omega$ - Speed deviation [Hz]", "scale": "freq_base",
              "ylim": (-0.15, 0.15)},
    "vmag": {"label": r"$|V|$ - Voltage [p.u.]", "floor": 0.9},
}


def load(path, scale=1.0):
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data[:, 0], data[:, 1:] * scale


def variable(reference):
    parts = Path(reference).name.split(".")
    return parts[-3] if parts[-2:] == ["ref", "csv"] else Path(reference).stem.rsplit(".", 1)[-1]


def resolve(base, path):
    path = Path(path)
    return path if path.is_absolute() else base / path


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


def style_axis(ax, letter, name, boundaries, ylabel, window, ylim):
    ax.set_xlim(0.0, window)
    ax.set_ylim(*ylim)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_ylabel(ylabel)
    ax.grid(True, which="major", axis="y", color="#D4D8DD", linewidth=0.42, alpha=0.75)
    ax.grid(True, which="major", axis="x", color="#E6E8EB", linewidth=0.36, alpha=0.5)
    for boundary in boundaries:
        ax.axvline(boundary, color="#2F3437", linestyle=":", linewidth=0.55, alpha=0.30, zorder=0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color(SPINE_COLOR)
    ax.tick_params(which="both", top=True, right=True)
    ax.text(0.025, 0.95, f"({letter})", transform=ax.transAxes, fontsize=8.5,
            fontweight="bold", va="top", ha="left")
    ax.text(0.975, 0.93, name, transform=ax.transAxes, fontsize=8.0,
            fontweight="bold", va="top", ha="right")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("study", type=Path)
    parser.add_argument("--rundir", type=Path)
    parser.add_argument("--outdir", type=Path)
    args = parser.parse_args()
    study_path = args.study.resolve()
    study_dir = study_path.parent
    run_dir = args.rundir.resolve() if args.rundir else study_dir
    out_dir = args.outdir.resolve() if args.outdir else study_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    study = json.loads(study_path.read_text())
    output = resolve(run_dir, study["output_file"])
    reference = resolve(study_dir, study["reference_file"])
    field = variable(reference)
    spec = FIELDS.get(field, {})
    ylabel = spec.get("label", f"{field} [p.u.]")

    case_path = resolve(study_dir, study["system_model_file"])
    if not case_path.exists():
        case_path = resolve(run_dir, study["system_model_file"])
    params = json.loads(case_path.read_text()).get("params", {})
    scale = params[spec["scale"]] if "scale" in spec else 1.0

    plt.rcParams.update(RC_PARAMS)
    time_out, values_out = load(output, scale)
    time_ref, values_ref = load(reference, scale)
    window = min(time_out[-1], time_ref[-1])
    boundaries = [event["time"] for event in study["events"] if 0.0 < event["time"] < window]

    figure, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(3.5, 3.0), sharex=True, sharey=True,
        gridspec_kw={"hspace": 0.12}, constrained_layout=True)

    lo_out, hi_out = add_traces(ax_top, time_out, values_out)
    lo_ref, hi_ref = add_traces(ax_bot, time_ref, values_ref)
    ylim = spec.get("ylim") or limits(
        min(lo_out, lo_ref), max(hi_out, hi_ref), spec.get("floor"), spec.get("ceil"))

    style_axis(ax_top, "a", "GridKit", boundaries, ylabel, window, ylim)
    style_axis(ax_bot, "b", "PowerWorld", boundaries, ylabel, window, ylim)
    ax_bot.set_xlabel(r"$t$ - Time [sec]")

    png = out_dir / f"{output.stem}.signals.png"
    figure.savefig(png, dpi=DPI)
    plt.close(figure)
    print(f"Wrote {png}")


if __name__ == "__main__":
    main()
