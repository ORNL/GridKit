#!/usr/bin/env python3
"""Plot every monitored field of the ACTIVSg10k fault study.

There is no reference solution for this case, so each figure shows the run on
its own terms: every monitored series as one trace, and beneath it the spread
those traces occupy, which is the only way thousands of series stay readable.
The series that swings furthest from its pre-fault value is drawn on top of
both panels and named, since a population plot hides individuals by design.

One figure per entry in the study's `monitors` list, written as
<output-file-stem>.png next to the CSV.

Usage:
    plot_wecc10k.py [study.json]
"""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#e3e2df"
POPULATION = "#2a78d6"
ACCENT = "#eb6834"

# Axis label and what the panel is showing, per monitored variable.
FIELDS = {
    "omega": ("$\\Delta\\omega$  (pu)", "generator speed deviation", "machine"),
    "speed": ("$\\omega$  (pu)", "generator speed", "machine"),
    "delta": ("$\\delta$  (rad)", "generator rotor angle", "machine"),
    "p": ("$P$  (pu)", "generator real power", "machine"),
    "q": ("$Q$  (pu)", "generator reactive power", "machine"),
    "vmag": ("$|V|$  (pu)", "bus voltage magnitude", "bus"),
    "vang": ("$\\angle V$  (rad)", "bus voltage angle", "bus"),
}


def column_key(column):
    """Reduce a GridKit column to the bare machine or bus it names.

    `Genrou_10737_1_genrou_omega` -> `10737 1`,  `Bus_30399_Vm` -> `30399`.
    """
    parts = column.split("_")
    if len(parts) < 3:
        return column.strip()
    device_class, body = parts[0], parts[1:-1]
    # Component ids repeat the class as a suffix; bus ids do not
    if body and body[-1].lower() == device_class.lower():
        body = body[:-1]
    return " ".join(body) if body else column.strip()


def load(path):
    with open(path) as f:
        labels = [column_key(c) for c in f.readline().strip().split(",")[1:]]
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return labels, data[:, 0], data[:, 1:]


def fault_details(study, study_dir):
    """Where the fault is, how long it lasts, and its impedance."""
    case = json.loads((study_dir / study["system_model_file"]).read_text())
    events = {e["type"]: e for e in study.get("events", [])}
    on = events.get("fault_on", {}).get("time", 0.0)
    off = events.get("fault_off", {}).get("time", 0.0)

    faults = [d for d in case["devices"] if d["class"] == "BusFault"]
    fault = faults[events.get("fault_on", {}).get("element_id", 0)]
    number = fault["ports"]["bus"]
    bus = next(b for b in case["buses"] if b["number"] == number)
    r = fault["params"].get("R", 0.0)
    x = fault["params"].get("X", 0.0)

    return on, off, [
        f"bus {number} {bus['name']} ({bus['params']['kv']:g} kV)",
        f"{on:g}-{off:g} s ({(off - on) * 1e3:.0f} ms)",
        f"Z = {r:g} + j{x:g} pu",
    ]


def plot_field(csv, png, field, fault):
    labels, t, y = load(csv)
    on, off, fault_lines = fault
    axis_label, description, unit_noun = FIELDS.get(
        field, (f"{field}  (pu)", field, "series"))

    # A population plot hides individuals, so lift out the one that moves most
    excursion = np.abs(y - y[0]).max(axis=0)
    worst = int(excursion.argmax())

    fig, axes = plt.subplots(2, 1, figsize=(11, 7.6), sharex=True, dpi=150)
    fig.patch.set_facecolor(SURFACE)

    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.axvspan(on, off, color=INK, alpha=0.05, lw=0)
        ax.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=INK_MUTED, labelsize=9)
        ax.set_ylabel(axis_label, color=INK_MUTED, fontsize=9)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0),
                            useMathText=True)
        ax.yaxis.get_offset_text().set(color=INK_MUTED, fontsize=9)

    # One LineCollection draws every series in a single artist; thousands of
    # Line2D objects would take minutes and render identically.
    segments = np.empty((y.shape[1], len(t), 2))
    segments[:, :, 0] = t
    segments[:, :, 1] = y.T
    # Enough traces overlap that per-trace alpha has to fall as the count rises
    alpha = float(np.clip(12.0 / y.shape[1], 0.03, 0.5))
    axes[0].add_collection(LineCollection(segments, colors=POPULATION, lw=0.4,
                                          alpha=alpha))
    axes[0].set_title(f"Every monitored {unit_noun}", loc="left", color=INK,
                      fontsize=12, fontweight="semibold", pad=8)

    lo, hi = y.min(), y.max()
    pad = 0.08 * (hi - lo) or 1e-12
    axes[0].set_ylim(lo - pad, hi + pad)

    # Where the series sit at wildly different levels -- machine ratings span
    # two orders of magnitude -- the spread of levels swamps the motion and
    # every trace reads flat. Re-reference to each series' own pre-fault value
    # so the dynamics share a scale. Where levels are already comparable
    # (speed deviation rests at zero) that plot would just repeat the first,
    # so show how the population spreads instead.
    levels = y[0]
    rebase = (levels.max() - levels.min()) > 3.0 * excursion.max()

    if rebase:
        d = y - levels
        segments[:, :, 1] = d.T
        axes[1].add_collection(LineCollection(segments, colors=POPULATION,
                                              lw=0.4, alpha=alpha))
        axes[1].plot(t, np.median(d, axis=1), color=POPULATION, lw=2.0,
                     solid_joinstyle="round", label="median")
        axes[1].set_title("Deviation from pre-fault value", loc="left",
                          color=INK, fontsize=12, fontweight="semibold", pad=8)
        axes[1].set_ylabel(axis_label.replace("$", "$\\Delta ", 1),
                           color=INK_MUTED, fontsize=9)
        span = np.abs(d).max()
        axes[1].set_ylim(-1.08 * span, 1.08 * span)
        highlight = d[:, worst]
    else:
        lo_b, hi_b = np.percentile(y, [5, 95], axis=1)
        axes[1].fill_between(t, y.min(axis=1), y.max(axis=1), color=POPULATION,
                             alpha=0.13, lw=0, label="min-max")
        axes[1].fill_between(t, lo_b, hi_b, color=POPULATION, alpha=0.30, lw=0,
                             label="5-95%")
        axes[1].plot(t, np.median(y, axis=1), color=POPULATION, lw=2.0,
                     solid_joinstyle="round", label="median")
        axes[1].set_title("Spread across the population", loc="left", color=INK,
                          fontsize=12, fontweight="semibold", pad=8)
        axes[1].set_ylim(lo - pad, hi + pad)
        highlight = y[:, worst]

    label = f"{unit_noun} {labels[worst]} (largest swing)"
    axes[0].plot(t, y[:, worst], color=ACCENT, lw=2.0, solid_joinstyle="round",
                 label=label)
    axes[1].plot(t, highlight, color=ACCENT, lw=2.0, solid_joinstyle="round",
                 label=label)

    row = int(np.abs(y[:, worst] - y[0, worst]).argmax())
    below = (y[row, worst] - axes[0].get_ylim()[0]) / (2 * pad + hi - lo) < 0.5
    axes[0].annotate(
        f"{unit_noun} {labels[worst]}",
        xy=(t[row], y[row, worst]),
        xytext=(30, 34 if below else -34), textcoords="offset points",
        color=INK, fontsize=9, va="bottom" if below else "top",
        bbox=dict(facecolor=SURFACE, edgecolor="none", alpha=0.85, pad=2),
        arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.8),
    )

    legend = axes[1].legend(loc="lower right", frameon=True, fontsize=9,
                            facecolor=SURFACE, edgecolor=GRID, ncol=2)
    for text in legend.get_texts():
        text.set_color(INK_MUTED)

    axes[1].set_xlim(t[0], t[-1])
    axes[1].set_xlabel("time (s)", color=INK_MUTED, fontsize=9)

    fig.text(0.055, 0.978, f"ACTIVSg10k  ·  {description}", color=INK,
             fontsize=15, fontweight="semibold", va="top")
    fig.text(0.055, 0.947, f"{y.shape[1]} {unit_noun}s · {len(t)} samples",
             color=INK_MUTED, fontsize=9.5, va="top")

    fig.text(0.985, 0.978, "Fault", color=INK, fontsize=10,
             fontweight="semibold", ha="right", va="top")
    for i, line in enumerate(fault_lines):
        fig.text(0.985, 0.951 - i * 0.0275, line, color=INK_MUTED, fontsize=9.5,
                 ha="right", va="top")

    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(png, facecolor=SURFACE)
    plt.close(fig)

    print(f"wrote {png}")
    print(f"  {y.shape[1]} {unit_noun}s · range [{y.min():.6g}, {y.max():.6g}]")
    print(f"  largest swing : {unit_noun} {labels[worst]} "
          f"({excursion[worst]:.6g} from rest, at t = {t[row]:.3f} s)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("study", nargs="?", type=Path,
                        default=Path(__file__).with_name("wecc10k.solver.json"))
    args = parser.parse_args()

    study_dir = args.study.parent
    study = json.loads(args.study.read_text())
    fault = fault_details(study, study_dir)

    for monitor in study.get("monitors", []):
        csv = study_dir / monitor["file_name"]
        if not csv.exists():
            print(f"skipping {monitor['file_name']}: not found")
            continue
        field = Path(monitor["file_name"]).stem.split(".")[-1]
        plot_field(csv, csv.with_suffix(".png"), field, fault)


if __name__ == "__main__":
    main()
