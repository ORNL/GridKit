"""Map the Texas7kGas pressure drop after the pipeline rupture.

Five maps of the whole gas network show the pressure drop since the rupture, on
a log scale, with each pipe colored by the mean drop of its end buses. The last
panel follows the drop at buses increasingly far from the rupture.

Usage: python3 plot.py Texas7kGas.csv

The figure is written next to this script.
"""

import csv
import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
CASE = HERE.parents[3] / "cases" / "PhasorDynamics" / "Texas7kGas" / "Texas7kGas.case.json"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3de"
NETWORK = "#d6d4cc"  # Pipes that dropped less than the color-scale floor
# One-hue orange ramp, light to dark: larger drop is darker
DROP = LinearSegmentedColormap.from_list("drop", ["#fbd3b6", "#f5a77a", "#eb6834", "#b84a17", "#6e2a0a", "#3a1404"])
DROP_FLOOR, DROP_CEILING = 1.0, 1000.0  # [kPa]

SNAPSHOTS = [1.0, 5.0, 10.0, 20.0, 30.0]  # [min] after the rupture
DISTANCES = [0.0, 50.0, 100.0, 200.0]  # [km] from the rupture
EARTH_RADIUS = 6371.0  # [km]


def read(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    data = np.array(rows[1:], dtype=float)
    return rows[0], data


def distance_km(lat0, lon0, lat, lon):
    x = np.radians(lon - lon0) * np.cos(np.radians(0.5 * (lat + lat0)))
    y = np.radians(lat - lat0)
    return EARTH_RADIUS * np.hypot(x, y)


def main(csv_path):
    header, data = read(csv_path)
    t = data[:, 0]
    case = json.load(open(CASE))
    event = json.load(open(HERE / "Texas7kGas.solver.json"))["events"][0]

    buses = case["buses"]
    number = np.array([b["number"] for b in buses])
    lat = np.array([b["extension"]["latitude"] for b in buses])
    lon = np.array([b["extension"]["longitude"] for b in buses])
    pressure = data[:, [header.index(f"GasBus_{b['name']}_p") for b in buses]]
    drop = (pressure[0] - pressure) / 1.0e3  # [kPa]

    index = {n: k for k, n in enumerate(number)}
    pipes = [(index[d["ports"]["bus1"]], index[d["ports"]["bus2"]], d["params"].get("D", 0.6))
             for d in case["devices"] if d["class"] in ("GasBranch", "GasCompressor")]
    segments = np.array([[(lon[a], lat[a]), (lon[b], lat[b])] for a, b, _ in pipes])
    widths = np.array([0.5 + 1.6 * diameter / 1.2 for _, _, diameter in pipes])
    ends = np.array([(a, b) for a, b, _ in pipes])

    rupture = index[event["bus"]]
    far = distance_km(lat[rupture], lon[rupture], lat, lon)
    norm = LogNorm(DROP_FLOOR, DROP_CEILING)

    fig, axes = plt.subplots(2, 3, figsize=(16.0, 10.0), facecolor=SURFACE, layout="constrained")
    maps = axes.flat[: len(SNAPSHOTS)]
    for ax, minutes in zip(maps, SNAPSHOTS):
        k = int(np.argmin(np.abs(t - (event["time"] + 60.0 * minutes))))
        pipe_drop = 0.5 * (drop[k, ends[:, 0]] + drop[k, ends[:, 1]])
        quiet = pipe_drop < DROP_FLOOR
        ax.add_collection(LineCollection(segments[quiet], colors=NETWORK, linewidths=widths[quiet], capstyle="round"))
        order = np.argsort(pipe_drop[~quiet])
        ax.add_collection(LineCollection(segments[~quiet][order], colors=DROP(norm(pipe_drop[~quiet][order])),
                                         linewidths=widths[~quiet][order] + 0.6, capstyle="round"))
        ax.scatter(lon[rupture], lat[rupture], marker="x", s=70, color=INK, linewidths=1.8, zorder=3)
        reached = np.sum(drop[k] >= DROP_FLOOR)
        ax.set_title(f"{minutes:.0f} min after the rupture", color=INK, fontsize=11)
        ax.text(0.03, 0.03, f"{reached} of {len(buses)} buses down by 1 kPa or more\nlargest drop {drop[k].max():.0f} kPa",
                transform=ax.transAxes, fontsize=8.5, color=INK_SECONDARY, va="bottom")
        ax.set_xlim(lon.min() - 0.3, lon.max() + 0.3)
        ax.set_ylim(lat.min() - 0.3, lat.max() + 0.3)
        ax.set_aspect(1.0 / np.cos(np.radians(lat.mean())))
        ax.set_facecolor(SURFACE)
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ax.spines.values():
            side.set_color(GRID)

    # Drop over time at the bus with the largest final drop near each distance
    ax = axes.flat[-1]
    minutes = (t - event["time"]) / 60.0
    after = minutes >= 0.0
    shades = [DROP(norm(v)) for v in (900.0, 200.0, 40.0, 8.0)]
    for target, color in zip(DISTANCES, shades):
        ring = np.flatnonzero(np.abs(far - target) <= max(5.0, 0.1 * target))
        bus = ring[np.argmax(drop[-1, ring])]
        ax.plot(minutes[after], drop[after, bus], color=color, lw=2.0)
        name = "At the rupture" if target == 0.0 else f"{target:.0f} km away"
        ax.annotate(f"{name}  {drop[-1, bus]:.0f} kPa", xy=(minutes[-1], drop[-1, bus]), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=8.5, color=INK)
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlim(0.0, minutes[-1])
    ax.set_xlabel("Time after the rupture [min]", color=INK_SECONDARY)
    ax.set_ylabel("Pressure drop [kPa]", color=INK_SECONDARY)
    ax.set_title("Drop spreading outward", color=INK, fontsize=11)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=INK_SECONDARY, labelsize=8.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=DROP)
    cbar = fig.colorbar(sm, ax=list(axes.flat), location="bottom", shrink=0.45, aspect=50)
    cbar.set_label("Pressure drop since the rupture [kPa], log scale;  gray: under 1 kPa;  x: rupture", color=INK_SECONDARY)
    cbar.ax.tick_params(colors=INK_SECONDARY, labelsize=8.5)
    cbar.outline.set_edgecolor(GRID)
    fig.suptitle(f"Texas7kGas: {event['w']:.0f} kg/s rupture at junction {event['bus']}, 30 minutes at a fixed 5 s step",
                 color=INK, fontsize=13, y=1.03)

    out = HERE / "Texas7kGas.png"
    fig.savefig(out, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Texas7kGas.csv")
