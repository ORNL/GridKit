#!/usr/bin/env python3
"""Render the actual fourteen-bus case and its fitted coupled line response."""

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/gridkit-coupledgrid-matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import numpy as np

INK = "#203340"
BLUE = "#0072B2"
GREEN = "#258575"
ORANGE = "#D55E00"
plt.rcParams.update({"font.size": 10, "svg.fonttype": "none", "savefig.dpi": 180})


def save(fig, root, name):
    root.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg", "pdf"):
        fig.savefig(root / f"{name}.{extension}", bbox_inches="tight")
    plt.close(fig)


def draw(case, root):
    positions = {1: (0, 0), 2: (3, 0), 3: (6, 0), 4: (9, 0),
                 5: (9, -4), 6: (6, -4), 7: (3, -4), 8: (0, -4),
                 9: (6, 3), 10: (12, 0), 11: (6, -7), 12: (3, -7),
                 13: (3, 3), 14: (15, 0)}
    bus_number = lambda name: int(name.rsplit("_", 1)[-1])
    devices = case["devices"]
    buses = [d for d in devices if d["class"] == "Bus"]
    if {bus_number(d["id"]) for d in buses} != set(positions):
        raise ValueError("One-line layout requires the specified fourteen bus IDs")
    fig, ax = plt.subplots(figsize=(17, 11))
    ax.set_aspect("equal")
    wire = lambda x, y, color=INK, lw=1.7: ax.plot(x, y, color=color, lw=lw, zorder=1)
    for device in devices:
        if device["class"] not in ("LineLumped", "Switch"):
            continue
        a, b = (bus_number(device["inputs"][key]) for key in ("bus1", "bus2"))
        x1, y1 = positions[a]
        x2, y2 = positions[b]
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        if device["class"] == "Switch":
            opened = device["params"]["open"]
            wire([x1, mx - .2], [y1, my], ORANGE)
            wire([mx + .2, x2], [my, y2], ORANGE)
            wire([mx - .2, mx + .15], [my, my + (.28 if opened else 0)], ORANGE, 2)
            for x in (mx - .2, mx + .2):
                ax.add_patch(Circle((x, my), .035, color=ORANGE, zorder=4))
            ax.text(mx, my + .6, "Load-step switch", ha="center", color=ORANGE, weight="bold")
            ax.text(mx, my - .57, "Initial state: " + ("open" if opened else "closed"), ha="center", fontsize=9)
        else:
            wire([x1, x2], [y1, y2], "#647687")
            length = device["params"]["dx"]
            label = f"{length:g} m"
            if x1 == x2:
                ax.text(mx + .14, my, label, va="center", fontsize=10,
                        bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
            else:
                ax.text(mx, my + .19, label, ha="center", fontsize=10,
                        bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    for bus in buses:
        n = bus_number(bus["id"])
        x, y = positions[n]
        wire([x - .32, x + .32], [y, y], INK, 5)
        ax.text(x - .52, y - .37, f"B{n}", ha="center", weight="bold", color=INK, fontsize=12)

    for device in devices:
        if device["class"] == "LoadZ":
            n = bus_number(device["inputs"]["bus"])
            x, y = positions[n]
            # Leave line connections visually separate from load branches.
            dx, dy = 0., (.95 if n in (4, 6) else -.95)
            wire([x, x + dx], [y, y], ORANGE)
            ax.annotate("", xy=(x + dx, y + dy), xytext=(x + dx, y),
                        arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.7, mutation_scale=14))
            ax.text(x + dx, y + dy + (.25 if dy > 0 else -.25), f"Load {n}", ha="center", color=ORANGE, fontsize=10)
        elif device["class"] == "Machine":
            n = bus_number(device["inputs"]["bus"])
            x, y = positions[n]
            sx, sy = x, y + (-1.25 if n == 3 else 1.25)
            wire([x, sx], [y, sy + .27 if n == 3 else sy - .27], BLUE)
            ax.add_patch(Circle((sx, sy), .27, facecolor="white", edgecolor=BLUE, lw=1.9, zorder=4))
            ax.text(sx, sy, "G", ha="center", va="center", weight="bold", color=BLUE, fontsize=13, zorder=5)
            offset = -.62 if n == 3 else .62
            ax.text(sx, sy + offset, f"Machine {n}", ha="center", color=BLUE, weight="bold")
            ax.text(sx, sy + offset + (-.28 if n == 3 else .28), "TGOV1 + IEEET1", ha="center", color=BLUE, fontsize=9)
        elif device["class"] == "DependentVoltageSource":
            n = bus_number(device["inputs"]["bus"])
            x, y = positions[n]
            sign = 1 if n in (9, 13) else -1
            sy = y + sign * 1.
            wire([x, x], [y, sy - sign * .3], GREEN)
            ax.add_patch(Rectangle((x - .5, sy - .3), 1., .6, edgecolor=GREEN, facecolor="white", lw=1.9))
            ax.plot([x - .5, x + .5], [sy - .3, sy + .3], color=GREEN, lw=1)
            ax.text(x - .24, sy + .12, "=", ha="center", color=GREEN, fontsize=14)
            ax.text(x + .22, sy - .17, "~", ha="center", color=GREEN, fontsize=14)
            ax.text(x, sy + sign * .68, f"IBR {n}", ha="center", color=GREEN, weight="bold")
            ax.text(x, sy + sign * .96, "PWM → Converter → DVS", ha="center", color=GREEN, fontsize=9)
        elif device["class"] == "VoltageSource":
            n = bus_number(device["inputs"]["bus"])
            x, y = positions[n]
            sx, sy = x - 1.4, y + 1.4
            wire([x, sx, sx], [y, y, sy - .3], BLUE)
            ax.add_patch(Circle((sx, sy), .3, facecolor="white", edgecolor=BLUE, lw=1.9))
            ax.text(sx, sy, "~", ha="center", va="center", color=BLUE, fontsize=20)
            ax.text(sx, sy + .62, "Grid source", ha="center", color=BLUE, weight="bold")
            ax.text(sx, sy + .9, "Frequency dependent Y(s)", ha="center", color=BLUE, fontsize=9)
    ax.text(10.1, -3.5, "COUPLED LINE MODEL", fontsize=11, color=INK, weight="bold")
    ax.text(10.1, -3.95, "13 three-phase lumped π sections", color=INK)
    ax.text(10.1, -4.4, "Full 3 × 3 Z(s) and Y(s) matrices", color=INK)
    ax.text(10.1, -4.85, "Untransposed ACSR conductors", color=INK)
    ax.text(10.1, -5.3, "Skin effect + lossy earth return", color=INK)
    ax.text(10.1, -6.15, "SOURCE / LOAD MODELS", fontsize=11, color=INK, weight="bold")
    ax.text(10.1, -6.6, "2 synchronous machines + exciters", color=INK)
    ax.text(10.1, -7.05, "3 open-loop IBRs · 900 Hz PWM", color=INK)
    ax.text(10.1, -7.5, "7 R–L loads + 1 switched resistor", color=INK)
    ax.set_xlim(-3, 16.3)
    ax.set_ylim(-9.5, 5.8)
    ax.axis("off")
    fig.suptitle("14-bus coupled EMT network · 13.8 kV · 60 Hz", fontsize=18, color=INK, y=.98)
    fig.text(.5, .028, "One line represents three conductors. Labels give physical length. Grid-source, converter-filter, line and seven load models have frequency dependence.",
             ha="center", fontsize=10, color="#526372")
    fig.tight_layout(rect=(0, .04, 1, .96))
    save(fig, root, "one_line")


def line_response(path, root):
    with path.open(newline="") as stream:
        names = next(csv.reader(stream))
    values = np.loadtxt(path, delimiter=",", skiprows=1)
    data = dict(zip(names, values.T))
    frequency = data["freq_hz"]
    fig, axs = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    for j, (row, col, label) in enumerate(((0, 0, "AA"), (1, 1, "BB"), (2, 2, "CC"),
                                          (0, 1, "AB"), (0, 2, "AC"), (1, 2, "BC"))):
        column = 0 if row == col else 1
        color = ("#0072B2", "#D55E00", "#8E5BB7")[j % 3]
        for part, axrow, scale in (("re", 0, 1000.), ("im", 1, 1e6 / (2 * np.pi * frequency))):
            raw = data[f"raw_Z{row}{col}_{part}"] * scale
            fitted = data[f"fit_Z{row}{col}_{part}"] * scale
            ax = axs[axrow, column]
            ax.semilogx(frequency, fitted, color=color, label=label, lw=1.5)
            ax.semilogx(frequency[::40], raw[::40], "o", color=color, ms=4, fillstyle="none")
    for column, title in enumerate(("Self terms", "Mutual terms")):
        axs[0, column].set_title(title)
        axs[0, column].set_ylabel("Resistance [Ω/km]")
        axs[1, column].set_ylabel("Inductance X / (2πf) [mH/km]")
        axs[1, column].set_xlabel("Frequency [Hz]")
        axs[0, column].legend(ncol=3, frameon=False)
    for ax in axs.flat:
        ax.grid(True, alpha=.2, which="both")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Physical coupled line impedance · OpenLine data and rational fit", fontsize=16)
    fig.text(.5, .018, "Solid: rational model. Open markers: physical calculation. Conductors: A(−1.2,10.5), B(0,11), C(1.2,10.5) m; earth 100 Ω·m.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .035, 1, .96))
    save(fig, root, "line_frequency_response")

    fig, axs = plt.subplots(1, 3, figsize=(14, 6), gridspec_kw={"width_ratios": [.75, 1.2, 1.2]})
    ax = axs[0]
    ax.axhline(0., color=INK, lw=2)
    for x, height, name, color in ((-1.2, 10.5, "A", BLUE), (0., 11., "B", ORANGE), (1.2, 10.5, "C", "#8E5BB7")):
        ax.plot([x, x], [0, height], ls=":", color="#a8b0b8", lw=.9)
        ax.scatter([x], [height], s=100, color=color, zorder=3)
        ax.text(x, height + .55, name, ha="center", weight="bold", color=color)
    ax.set_xlim(-2., 2.)
    ax.set_ylim(-.4, 12.5)
    ax.set_xticks([-1.2, 0., 1.2])
    ax.set_xlabel("Horizontal position [m]")
    ax.set_ylabel("Height above earth [m]")
    ax.set_title("Untransposed geometry")
    ax.text(0., 5., "ACSR Linnet · 20 °C\nRadius 9.144 mm\nEarth 100 Ω·m", ha="center", fontsize=9,
            bbox=dict(facecolor="white", edgecolor="none", alpha=.95))
    ax.spines[["top", "right"]].set_visible(False)
    index = int(np.argmin(abs(frequency - 60.)))
    impedance = np.array([[complex(data[f"raw_Z{r}{c}_re"][index], data[f"raw_Z{r}{c}_im"][index])
                           for c in range(3)] for r in range(3)]) * 1000
    capacitance = np.array([[data[f"raw_Y{r}{c}_im"][index] for c in range(3)] for r in range(3)])
    capacitance *= 1e12 / (2 * np.pi * frequency[index])
    for ax, matrix, title, unit in ((axs[1], impedance, f"Series Z at {frequency[index]:.2f} Hz", "Ω/km (R + jX)"),
                                     (axs[2], capacitance, "Shunt capacitance C", "nF/km")):
        color_values = np.abs(matrix) if np.iscomplexobj(matrix) else matrix
        im = ax.imshow(color_values, cmap="Blues" if np.iscomplexobj(matrix) else "coolwarm", alpha=.8)
        for r in range(3):
            for c in range(3):
                value = matrix[r, c]
                label = f"{value.real:.3f}\n+j{value.imag:.3f}" if np.iscomplexobj(matrix) else f"{value:.3f}"
                ax.text(c, r, label, ha="center", va="center", color=INK, fontsize=10,
                        bbox=dict(facecolor="white", alpha=.78, edgecolor="none", pad=2.5))
        ax.set_xticks(range(3), list("ABC"))
        ax.set_yticks(range(3), list("ABC"))
        ax.set_title(title)
        ax.set_xlabel(unit)
        ax.set_ylabel("Phase")
    fig.suptitle("Physical line geometry and full phase coupling", fontsize=16)
    fig.text(.5, .025, "Matrices are per unit length; every line uses its stated physical length. Conductors are shown enlarged; no transposition, shield or neutral.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .06, 1, .93))
    save(fig, root, "line_geometry_and_coupling")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--line-data", type=Path)
    args = parser.parse_args()
    draw(json.loads(args.case.read_text()), args.output)
    if args.line_data:
        line_response(args.line_data, args.output)
    print(f"Rendered network figures to {args.output}")


if __name__ == "__main__":
    main()
