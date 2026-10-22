#!/usr/bin/env python3
"""Plot the system response and protection relays of a DynamicSimulation run."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e3e2df"
CONTEXT = "#a3a29d"
COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
FIELDS = {
    "Vm": ("|V|", "pu", "Bus voltage magnitude"),
    "omega": ("ω", "pu", "Generator speed deviation"),
    "speed": ("ω", "pu", "Generator speed"),
}
RELAYS = ("DefiniteTimeRelay", "InverseTimeRelay")
MEASURED = ("m", "im")
COMMIT_LEVEL = 0.5
TRIP_LEVEL = 0.75


def load(path):
    """Read column names, times, and values from a CSV."""
    with path.open() as stream:
        columns = stream.readline().strip().split(",")[1:]
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return columns, data[:, 0], data[:, 1:]


def split(column):
    """Split a monitor column into model class, device name, and variable."""
    model, rest = column.split("_", 1)
    name, variable = rest.rsplit("_", 1)
    return model, name, variable


def fault(solver, case):
    """Fault interval and description lines; a fault without fault_off is permanent."""
    on = next((event for event in solver["events"] if event["type"] == "fault_on"), None)
    if on is None:
        return None
    off = next((event for event in solver["events"] if event["type"] == "fault_off"), None)
    bus = next(bus for bus in case["buses"] if bus["number"] == on["bus"])
    interval = f"from {on['time']:g} s, permanent"
    if off is not None:
        interval = f"{on['time']:g}–{off['time']:g} s  ({(off['time'] - on['time']) * 1e3:g} ms)"
    return (
        on["time"],
        None if off is None else off["time"],
        (
            f"bus {bus['number']}  {bus['name'].strip()}  ({bus['params']['kv']:g} kV)",
            interval,
            f"Z = {on['R']:g} + j{on['X']:g} pu",
        ),
    )


def style(axis, title, ylabel, fault_data, end):
    axis.set_facecolor(SURFACE)
    if fault_data is not None:
        fault_on, fault_off, _ = fault_data
        axis.axvspan(fault_on, end if fault_off is None else fault_off, color=INK, alpha=0.05, lw=0)
    axis.set_title(title, loc="left", color=INK, fontsize=12, fontweight="semibold", pad=8)
    axis.grid(True, color=GRID, lw=0.6)
    axis.set_axisbelow(True)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set_color(GRID)
    axis.tick_params(colors=MUTED, labelsize=9)
    axis.set_ylabel(ylabel, color=MUTED, fontsize=9)


def header(figure, title, subtitle, fault_data):
    figure.text(0.055, 0.978, title, color=INK, fontsize=15, fontweight="semibold", va="top")
    figure.text(0.055, 0.978 - 0.031 * 9.5 / figure.get_figheight(), subtitle, color=MUTED, fontsize=9.5, va="top")
    if fault_data is None:
        return
    figure.text(0.985, 0.978, "Fault", color=INK, fontsize=10, fontweight="semibold", ha="right", va="top")
    for index, line in enumerate(fault_data[2]):
        y = 0.978 - (0.027 + index * 0.0225) * 9.5 / figure.get_figheight()
        figure.text(0.985, y, line, color=MUTED, fontsize=9.5, ha="right", va="top")


def save(figure, png):
    figure.tight_layout(rect=(0, 0, 1, 1 - 0.085 * 9.5 / figure.get_figheight()))
    figure.savefig(png, facecolor=SURFACE)
    plt.close(figure)
    print(f"Saved {png}")


def plot_system(columns, time, values, png, case_name, fault_data):
    """One panel per monitored system field, every series in one thin color."""
    fields = [field for field in FIELDS if any(split(column)[2] == field for column in columns)]
    if not fields:
        return
    figure, axes = plt.subplots(len(fields), 1, figsize=(11, 0.6 + 3.0 * len(fields)), sharex=True, dpi=150, squeeze=False)
    figure.patch.set_facecolor(SURFACE)
    series = 0
    for axis, field in zip(axes[:, 0], fields):
        symbol, unit, description = FIELDS[field]
        indices = [index for index, column in enumerate(columns) if split(column)[2] == field]
        series += len(indices)
        style(axis, f"{description}  ·  {len(indices)} series", f"{symbol}  ({unit})", fault_data, time[-1])
        if len(indices) > len(COLORS):
            axis.plot(time, values[:, indices], color=COLORS[0], lw=0.6, alpha=0.35, solid_joinstyle="round")
            continue
        for index, color in zip(indices, COLORS):
            model, name, _ = split(columns[index])
            axis.plot(time, values[:, index], color=color, lw=1.4, solid_joinstyle="round", label=f"{model} {name}")
        legend = axis.legend(loc="lower right", fontsize=9, frameon=True, labelcolor=INK)
        legend.get_frame().set(facecolor=SURFACE, edgecolor=GRID, lw=0.8)
    axes[-1, 0].set_xlim(time[0], time[-1])
    axes[-1, 0].set_xlabel("time (s)", color=MUTED, fontsize=9)
    header(figure, f"{case_name} · system response", f"{series} series · {len(time)} samples", fault_data)
    save(figure, png)


def plot_protection(columns, time, values, png, case, fault_data):
    """Measured current, latch, and trip command of every relay; tripped relays in color."""
    relays = {}
    pickups = {
        device["id"]: device["params"].get("Ipickup", device["params"].get("Mpickup"))
        for device in case["devices"]
        if device["class"] in RELAYS
    }
    for index, column in enumerate(columns):
        model, name, variable = split(column)
        if model in RELAYS:
            relays.setdefault(name, {})["measured" if variable in MEASURED else variable] = values[:, index]
    if not relays:
        return

    trip_times = {}
    for name, relay in relays.items():
        if relay["trip"].max() > COMMIT_LEVEL:
            trip_times[name] = time[np.argmax(relay["trip"] > COMMIT_LEVEL)]
    tripped = sorted(trip_times, key=trip_times.get)[: len(COLORS)]
    others = [name for name in relays if name not in tripped]

    start = time[0]
    stop = time[-1]
    if fault_data is not None:
        start = max(time[0], fault_data[0] - 0.1)
        stop = min(time[-1], np.ceil(10 * (max([fault_data[0], *trip_times.values()]) + 0.4)) / 10)

    figure, axes = plt.subplots(3, 1, figsize=(11, 9.5), sharex=True, dpi=150)
    figure.patch.set_facecolor(SURFACE)
    panels = (
        ("measured", "Measured magnitude", "|m|  (pu)"),
        ("x", "Relay latch", "x"),
        ("trip", "Trip command", "s"),
    )
    for axis, (key, title, ylabel) in zip(axes, panels):
        style(axis, title, ylabel, fault_data, time[-1])
        for name in others:
            axis.plot(time, relays[name][key], color=CONTEXT, lw=0.6, alpha=0.6, solid_joinstyle="round")
        for name, color in zip(tripped, COLORS):
            axis.plot(time, relays[name][key], color=color, lw=1.6, solid_joinstyle="round")

    for name, color in zip(tripped, COLORS):
        if pickups.get(name) is not None:
            axes[0].axhline(pickups[name], color=color, lw=0.9, ls=(0, (4, 3)))
    for level, label in ((COMMIT_LEVEL, "lockout commit"), (TRIP_LEVEL, "trip level")):
        axes[1].axhline(level, color=MUTED, lw=0.8, ls=(0, (2, 3)))
        axes[1].text(start, level, f" {label}", color=MUTED, fontsize=8.5, va="bottom", ha="left")
    for axis in axes[1:]:
        axis.set_ylim(-0.05, 1.08)

    handles = [
        Line2D([], [], color=color, lw=1.6, label=f"{name}   trip {trip_times[name]:.3f} s")
        for name, color in zip(tripped, COLORS)
    ]
    if others:
        handles.append(Line2D([], [], color=CONTEXT, lw=1.0, label=f"not tripped ({len(others)})"))
    if any(pickups.get(name) is not None for name in tripped):
        handles.append(Line2D([], [], color=MUTED, lw=0.9, ls=(0, (4, 3)), label="pickup"))
    legend = axes[0].legend(handles=handles, loc="upper right", fontsize=9, frameon=True, labelcolor=INK)
    legend.get_frame().set(facecolor=SURFACE, edgecolor=GRID, lw=0.8)

    axes[-1].set_xlim(start, stop)
    axes[-1].set_xlabel("time (s)", color=MUTED, fontsize=9)
    case_name = case.get("header", {}).get("case_name", png.stem)
    subtitle = f"{len(relays)} relays · {len(trip_times)} tripped · {start:.1f}–{stop:.1f} s of {time[-1]:g} s"
    header(figure, f"{case_name} · protection", subtitle, fault_data)
    save(figure, png)


def main():
    """Plot the output CSV of a solver JSON.

    After running the simulation:

        python plot.py path/to/ThreeBusRelay.solver.json --dir figures

    Leave out --dir to save the plots beside the output CSV.
    """
    parser = argparse.ArgumentParser(description="Plot a DynamicSimulation run.")
    parser.add_argument("solver", type=Path, help="*.solver.json from the completed run's build directory")
    parser.add_argument("--dir", type=Path, help="PNG folder (default: beside the output CSV)")
    args = parser.parse_args()
    solver_path = args.solver.resolve()
    solver_dir = solver_path.parent
    try:
        solver = json.loads(solver_path.read_text())
        case = json.loads((solver_dir / solver["system_model_file"]).read_text())
        case_name = case.get("header", {}).get("case_name", solver_path.stem)
        fault_data = fault(solver, case)
        output = solver_dir / solver["output_file"]
        folder = (args.dir or output.parent).resolve()
        folder.mkdir(parents=True, exist_ok=True)
        columns, time, values = load(output)
        plot_system(columns, time, values, folder / f"{output.stem}.system.png", case_name, fault_data)
        plot_protection(columns, time, values, folder / f"{output.stem}.protection.png", case, fault_data)
    except KeyError as error:
        parser.error(f"Missing JSON field: {error}")
    except (OSError, ValueError, StopIteration) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
