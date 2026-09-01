#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e3e2df"
COLORS = ("#2a78d6", "#eb6834", "#1baf7a")
FIELDS = {
    "omega": ("Δω  (pu)", "generator speed deviation"),
    "speed": ("ω  (pu)", "generator speed"),
    "delta": ("δ  (rad)", "generator rotor angle"),
    "p": ("P  (pu)", "generator real power"),
    "q": ("Q  (pu)", "generator reactive power"),
    "vmag": ("|V|  (pu)", "bus voltage magnitude"),
    "Va": ("∠V  (rad)", "bus voltage angle"),
}


def sci(value):
    exponent = int(np.floor(np.log10(abs(value)))) if value else 0
    return f"{value / 10**exponent:.3f} \\times 10^{{{exponent}}}"


def column_key(column):
    parts = column.split("_")
    if len(parts) < 3:
        return column.strip()
    body = parts[1:-1]
    if body and body[-1].lower() == parts[0].lower():
        body.pop()
    return " ".join(body) or column.strip()


def load(path):
    with path.open() as stream:
        labels = [column_key(column) for column in stream.readline().strip().split(",")[1:]]
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return labels, data[:, 0], data[:, 1:]


def resolve(base, path):
    path = Path(path)
    return path if path.is_absolute() else base / path


def find_case(study, study_dir, run_dir):
    path = Path(study["system_model_file"])
    if path.is_absolute():
        return path
    source = study_dir / path
    return source if source.exists() else run_dir / path


def fault(study, case):
    on = next(event for event in study["events"] if event["type"] == "fault_on")
    off = next(event for event in study["events"] if event["type"] == "fault_off")
    devices = [device for device in case["devices"] if device["class"] == "BusFault"]
    device = devices[on["element_id"]]
    bus = next(bus for bus in case["buses"] if bus["number"] == device["ports"]["bus"])
    resistance = device["params"]["R"]
    reactance = device["params"]["X"]
    return on["time"], off["time"], (
        f"bus {bus['number']}  {bus['name'].strip()}  ({bus['params']['kv']:g} kV)",
        f"{on['time']:g}–{off['time']:g} s  ({(off['time'] - on['time']) * 1e3:g} ms)",
        f"Z = {resistance:g} + j{reactance:g} pu",
    )


def variable(reference):
    parts = Path(reference).name.split(".")
    return parts[-3] if parts[-2:] == ["ref", "csv"] else Path(reference).stem.rsplit(".", 1)[-1]


def plot(output, reference, png, title, ylabel, fault_data):
    labels, output_time, values = load(output)
    reference_labels, reference_time, reference_values = load(reference)
    if len(output_time) != len(reference_time):
        raise ValueError(f"Sample-count mismatch: {len(output_time)} != {len(reference_time)}")
    if not np.allclose(output_time, reference_time, rtol=0.0, atol=1e-6):
        raise ValueError("Output and reference time grids differ")
    if labels != reference_labels and sorted(labels) == sorted(reference_labels):
        values = values[:, [labels.index(label) for label in reference_labels]]

    samples = len(output_time)
    error = values - reference_values
    baseline = reference_values
    steps = np.diff(reference_time)
    dt = np.median(steps[steps > 0])
    duration = reference_time[-1] - reference_time[0]
    series = error.shape[1]
    rmse = np.sqrt(dt / (series * duration)) * np.linalg.norm(error)
    inf = np.abs(error).max()
    rel_rmse = np.linalg.norm(error) / np.linalg.norm(baseline)
    rel_inf = inf / np.abs(baseline).max()
    fault_on, fault_off, fault_lines = fault_data

    figure, axes = plt.subplots(3, 1, figsize=(11, 9.5), sharex=True, dpi=150)
    figure.patch.set_facecolor(SURFACE)
    panels = (
        (output_time, values, COLORS[0], "GridKit"),
        (reference_time, reference_values, COLORS[1], "PowerWorld"),
        (reference_time, error, COLORS[2], "Difference"),
    )
    for axis, (time, data, color, panel_title) in zip(axes, panels):
        axis.set_facecolor(SURFACE)
        axis.axvspan(fault_on, fault_off, color=INK, alpha=0.05, lw=0)
        axis.plot(time, data, color=color, lw=0.6, alpha=0.35, solid_joinstyle="round")
        axis.set_title(panel_title, loc="left", color=INK, fontsize=12, fontweight="semibold", pad=8)
        axis.grid(True, color=GRID, lw=0.6)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_color(GRID)
        axis.spines["bottom"].set_color(GRID)
        axis.tick_params(colors=MUTED, labelsize=9)
        axis.set_ylabel(ylabel, color=MUTED, fontsize=9)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
        axis.yaxis.get_offset_text().set(color=MUTED, fontsize=9)

    lower = min(values.min(), reference_values.min())
    upper = max(values.max(), reference_values.max())
    padding = 0.08 * (upper - lower)
    lower, upper = lower - padding, upper + padding
    axes[0].set_ylim(lower, upper)
    axes[1].set_ylim(lower, upper)
    axes[2].set_ylim((lower, upper) if lower < 0 < upper else (-0.5 * (upper - lower), 0.5 * (upper - lower)))
    axes[2].set_xlim(output_time[0], output_time[-1])
    axes[2].set_xlabel("time (s)", color=MUTED, fontsize=9)
    axes[2].add_patch(
        Rectangle(
            (0.558, 0.700),
            0.440,
            0.290,
            transform=axes[2].transAxes,
            facecolor=SURFACE,
            edgecolor=GRID,
            lw=0.8,
            zorder=3,
        )
    )

    cells = (
        (0.945, "", "$\\epsilon_{\\mathrm{RMSE}}$", "$\\epsilon_{\\infty}$", MUTED),
        (0.845, "Absolute", f"${sci(rmse)}$ pu", f"${sci(inf)}$ pu", INK),
        (0.755, "Relative", f"{rel_rmse * 100:.3f} %", f"{rel_inf * 100:.3f} %", INK),
    )
    for y, row_label, rmse_cell, inf_cell, color in cells:
        axes[2].text(0.578, y, row_label, transform=axes[2].transAxes, ha="left", va="center", color=MUTED, fontsize=10, zorder=4)
        axes[2].text(0.800, y, rmse_cell, transform=axes[2].transAxes, ha="right", va="center", color=color, fontsize=10, zorder=4)
        axes[2].text(0.985, y, inf_cell, transform=axes[2].transAxes, ha="right", va="center", color=color, fontsize=10, zorder=4)

    row, column = np.unravel_index(np.abs(error).argmax(), error.shape)
    x_min, x_max = axes[2].get_xlim()
    x_fraction = (reference_time[row] - x_min) / (x_max - x_min)
    axes[2].annotate(
        f"max at {reference_labels[column]},  t = {reference_time[row]:.3f} s",
        xy=(reference_time[row], error[row, column]),
        xytext=(0.32 if x_fraction < 0.25 else 0.03, 0.58),
        textcoords=axes[2].transAxes,
        color=INK,
        fontsize=9,
        ha="left",
        va="center",
        bbox={"facecolor": SURFACE, "edgecolor": "none", "alpha": 0.85, "pad": 2},
        arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8},
    )
    figure.text(0.055, 0.978, title, color=INK, fontsize=15, fontweight="semibold", va="top")
    figure.text(0.055, 0.947, f"{series} series · {samples} samples", color=MUTED, fontsize=9.5, va="top")
    figure.text(0.985, 0.978, "Fault", color=INK, fontsize=10, fontweight="semibold", ha="right", va="top")
    for index, line in enumerate(fault_lines):
        figure.text(0.985, 0.951 - index * 0.0225, line, color=MUTED, fontsize=9.5, ha="right", va="top")
    figure.tight_layout(rect=(0, 0, 1, 0.915))
    figure.savefig(png, facecolor=SURFACE)
    plt.close(figure)


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
    case = json.loads(find_case(study, study_dir, run_dir).read_text())
    case_name = case.get("header", {}).get("case_name", study_path.stem)
    fault_data = fault(study, case)
    output = resolve(run_dir, study["output_file"])
    reference = resolve(study_dir, study["reference_file"])
    field = variable(reference)
    ylabel, description = FIELDS.get(field, (f"{field}  (pu)", field))
    plot(output, reference, out_dir / f"{output.stem}.png", f"{case_name} · {description}", ylabel, fault_data)


if __name__ == "__main__":
    main()
