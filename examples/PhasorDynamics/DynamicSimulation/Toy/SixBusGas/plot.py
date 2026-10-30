"""Plot the SixBusGas response: gas pressures and flows, machine power and speed.

Usage: python3 plot.py SixBusGas.csv

The figure is written next to this script.
"""

import csv
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3de"
# Categorical slots 1-4, fixed per entity in every panel
G2, G3, SUPPLY, JUNCTION = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"

F_SYS = 60.0  # [Hz]


def read(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    data = np.array(rows[1:], dtype=float)
    return {name: data[:, k] for k, name in enumerate(rows[0])}


def label_ends(ax, series):
    """Label each line at its end in the right margin, spread to avoid overlap."""
    low, high = ax.get_ylim()
    gap = 0.11 * (high - low)
    ends = sorted((values[-1], name, color) for name, values, color in series)
    placed = []
    for value, name, color in ends:
        y = max(value, placed[-1] + gap) if placed else value
        placed.append(y)
        ax.annotate(name, xy=(1.0, value), xycoords=("axes fraction", "data"), xytext=(1.015, y),
                    textcoords=("axes fraction", "data"), va="center", fontsize=8.5, color=INK,
                    arrowprops={"arrowstyle": "-", "color": color, "lw": 1.0, "relpos": (0.0, 0.5)})


def panel(ax, t, series, ylabel, legend_loc):
    for name, values, color in series:
        ax.plot(t, values, color=color, lw=1.6, label=name)
    ax.set_ylabel(ylabel, color=INK_SECONDARY)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SECONDARY, labelsize=8.5)
    ax.legend(loc=legend_loc, fontsize=8, frameon=False, labelcolor=INK, ncol=2)
    label_ends(ax, series)


def main(csv_path):
    d = read(csv_path)
    t = d["t"]

    pressure = [
        ("Supply (101)", d["GasBus_Gas 101_p"] / 1.0e6, SUPPLY),
        ("Junction (102)", d["GasBus_Gas 102_p"] / 1.0e6, JUNCTION),
        ("G2 gas bus (103)", d["GasBus_Gas 103_p"] / 1.0e6, G2),
        ("G3 gas bus (104)", d["GasBus_Gas 104_p"] / 1.0e6, G3),
    ]
    fuel = [
        ("G2", d["GasBranch_GP_102_103_q2"], G2),
        ("G3", d["GasBranch_GP_102_104_q2"], G3),
    ]
    power = [
        ("G2", d["Genrou_G2_p"], G2),
        ("G3", d["Genrou_G3_p"], G3),
    ]
    speed = [
        ("G2", (d["Genrou_G2_speed"] - 1.0) * F_SYS * 1.0e3, G2),
        ("G3", (d["Genrou_G3_speed"] - 1.0) * F_SYS * 1.0e3, G3),
    ]

    fig, axes = plt.subplots(4, 1, figsize=(8.0, 9.5), sharex=True, facecolor=SURFACE)
    panel(axes[0], t, pressure, "Gas pressure [MPa]", "center right")
    panel(axes[1], t, fuel, "Fuel mass flow [kg/s]", "center right")
    panel(axes[2], t, power, "Electrical power [p.u.]", "center right")
    panel(axes[3], t, speed, "Speed deviation [mHz]", "lower right")

    axes[-1].set_xlabel("Time [s]", color=INK_SECONDARY)
    axes[-1].set_xlim(t[0], t[-1])
    fig.suptitle("SixBusGas: gas-fired machines in steady state", color=INK, fontsize=11)
    fig.subplots_adjust(left=0.1, right=0.8, top=0.95, bottom=0.06, hspace=0.18)

    out = pathlib.Path(__file__).resolve().parent / "SixBusGas.png"
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "SixBusGas.csv")
