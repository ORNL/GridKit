"""Figures of an EMT fault study: python3 plot.py <Case> <output.csv>, run from the example folder."""
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[2]

BACKGROUND = "#FCFCFA"
INK = "#0B0B0B"
MUTED = "#555555"
GRID = "#E5E5E3"
SPINE = "#DCDCDA"
SHADE = "#EFEFED"
FAMILY = "#6F8493"
HIGHLIGHT = "#11639A"
EXTREME = "#A6492B"
PHASES = {"A": "#2B78D9", "B": "#EC6A3A", "C": "#1FA97A"}

DPI = 200
WIDTH = 9.12
LEFT = 0.0948
RIGHT = 0.970

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "figure.facecolor": BACKGROUND,
    "axes.facecolor": BACKGROUND,
    "savefig.facecolor": BACKGROUND,
    "axes.edgecolor": SPINE,
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "xtick.color": "#444444",
    "ytick.color": "#444444",
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.labelcolor": MUTED,
    "axes.labelsize": 9,
})


def place_name(name):
    """Leading words of a bus name, without the voltage and circuit numbers."""
    words = []
    for word in name.split():
        word = re.sub(r"\d+(\.\d+)?$", "", word)
        if not word:
            break
        words.append(word.capitalize())
    return " ".join(words)


def number(value):
    return f"{value:g}"


class Study:
    def __init__(self, case, output):
        solver = json.loads(next(Path(".").glob("*.solver.json")).read_text())
        folder = ROOT / "cases" / "EMT" / case
        self.case = case
        self.data = json.loads((folder / f"{case}.case.json").read_text())
        self.buses = {str(bus["number"]): bus for bus in self.data["buses"]}
        self.devices = {device["id"]: device for device in self.data["devices"]}
        self.fault = (solver["events"][0]["time"], solver["events"][1]["time"])
        self.tmax = solver["tmax"]
        self.fault_bus = str(self.devices["fault_switch"]["terminals"]["bus1"])

        keep = ("t", "Bus_", "Machine_", "Tgov1_", "GastPti_", "Hygov_", "SexsPti_", "Ieeet1_", "Esdc1a_", "Switch_", "Pmu_")
        self.frame = pd.read_csv(output, usecols=lambda column: column.startswith(keep))
        self.t = self.frame["t"].to_numpy()

    def columns(self, pattern):
        return [column for column in self.frame.columns if re.fullmatch(pattern, column)]

    def bus_voltage(self, bus):
        return self.frame[[f"Bus_{bus}_v[{n}]" for n in (1, 2, 3)]].to_numpy()

    def measured_buses(self):
        """RMS voltage column of one Pmu on each bus that has one."""
        columns = {}
        for column in self.frame.columns:
            if match := re.fullmatch(r"Pmu_(\w+)_vrms", column):
                columns.setdefault(str(self.devices[match.group(1)]["inputs"]["v"]), column)
        return columns

    def window(self, lower, upper):
        return (self.t >= lower) & (self.t <= upper)


def figure(title, tag, panels):
    """Figure with stacked panels in the layout of the reference figures."""
    height = {2: 5.605, 3: 7.125}[panels]
    fig = plt.figure(figsize=(WIDTH, height), dpi=DPI)
    if panels == 2:
        boxes = [(0.543, 0.2935), (0.098, 0.2935)]
        title_y = 0.944
    else:
        boxes = [(0.676, 0.193), (0.381, 0.193), (0.0856, 0.193)]
        title_y = 0.946
    axes = [fig.add_axes([LEFT, bottom, RIGHT - LEFT, size]) for bottom, size in boxes]
    fig.text(LEFT, title_y, title, fontsize=16, fontweight="bold", color=INK, ha="left", va="center")
    fig.text(RIGHT, title_y, tag, fontsize=9, color=MUTED, ha="right", va="center")
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(length=3, width=0.8)
        ax.set_axisbelow(True)
    axes[-1].set_xlabel("Time (s)")
    return fig, axes


def heading(ax, left, right=None):
    ax.text(0, 1.045, left, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, ha="left", va="bottom")
    if right:
        ax.text(1, 1.045, right, transform=ax.transAxes, fontsize=9, color=HIGHLIGHT, ha="right", va="bottom")


def time_axis(ax, study, lower, upper, ticks):
    ax.axvspan(*study.fault, color=SHADE, zorder=0)
    ax.set_xlim(lower, upper)
    ax.set_xticks(ticks)
    ax.set_xticklabels([number(tick) for tick in ticks])


def fault_ticks(study, lower, upper):
    start, end = study.fault
    return [lower] + list(np.round(np.arange(start, end + 1e-9, 0.025), 3)) + [upper]


def record_ticks(study):
    start, end = study.fault
    step = 2.0
    if study.tmax <= 2:
        step = 0.2
    ticks = [tick for tick in np.arange(0, study.tmax + 1e-9, step) if abs(tick - start) > 0.4 * step and abs(tick - end) > 0.4 * step]
    marks = [start, end]
    if end - start < 0.4 * step:
        marks = [start]
    return sorted(ticks + marks)


def family(ax, t, values, highlight):
    """Every column in the family colour, and column `highlight` on top."""
    ax.plot(t, values, color=FAMILY, linewidth=0.7, alpha=0.4)
    ax.plot(t, values[:, highlight], color=HIGHLIGHT, linewidth=1.5)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))


def thin(study, lower, upper):
    """Every sample near the fault and every fourth elsewhere."""
    near = study.window(study.fault[0] - 0.1, study.fault[1] + 0.2)
    index = np.arange(len(study.t))
    return (study.window(lower, upper)) & (near | (index % 4 == 0))


def voltage(study, window):
    columns = study.measured_buses()
    buses = list(columns)
    magnitude = np.empty((len(study.t), len(buses)))
    for k, bus in enumerate(buses):
        magnitude[:, k] = study.frame[columns[bus]].to_numpy() / (1.0e3 * study.buses[bus]["params"]["kv"])

    fault = buses.index(study.fault_bus)
    label = f"Bus {study.fault_bus} · {place_name(study.buses[study.fault_bus]['name'])}"
    row, column = np.unravel_index(np.argmax(magnitude), magnitude.shape)

    fig, axes = figure("Bus voltage magnitude", f"{study.case}  /  EMT", 2)
    spans = [(0.0, study.tmax, record_ticks(study), "full record"), (*window, fault_ticks(study, *window), "fault and clearing")]
    for ax, (lower, upper, ticks, name) in zip(axes, spans):
        mask = thin(study, lower, upper)
        family(ax, study.t[mask], magnitude[mask], fault)
        time_axis(ax, study, lower, upper, ticks)
        ax.set_ylabel("Voltage (pu)")
        heading(ax, f"{len(buses)} buses · {name}", label)
        ax.plot(study.t[row], magnitude[row, column], "o", color=EXTREME, markersize=4, clip_on=False, zorder=5)
        offset, side = -9, "right"
        if study.t[row] - lower < 0.3 * (upper - lower):
            offset, side = 9, "left"
        ax.annotate(f"Bus {buses[column]} · {magnitude[row, column]:.2f} pu", (study.t[row], magnitude[row, column]),
                    xytext=(offset, -9), textcoords="offset points", ha=side, va="center", fontsize=8.5, color=EXTREME)
    return fig


def machine_bus(column):
    return re.match(r"[A-Za-z0-9]+_(\w+?)_", column).group(1)


def machines(study, window):
    omega = study.columns(r"Machine_.*_omega")
    torque = [column.replace("_omega", "_te") for column in omega]
    mask = study.window(*window)
    frequency = study.data["params"]["freq_base"] * study.frame[omega].to_numpy()[mask]
    te = study.frame[torque].to_numpy()[mask]
    t = study.t[mask]

    highlight = int(np.argmax(np.abs(frequency - frequency[0]).max(axis=0)))
    label = f"Machine at bus {machine_bus(omega[highlight])}"

    fig, axes = figure("Synchronous machines", f"{study.case}  /  EMT", 2)
    panels = [(frequency, "Hz", "Rotor electrical frequency"), (te, "pu", "Electromagnetic torque")]
    for ax, (values, unit, name) in zip(axes, panels):
        family(ax, t, values, highlight)
        time_axis(ax, study, *window, fault_ticks(study, *window))
        ax.set_ylabel(unit)
        heading(ax, f"{name} · {len(omega)} machines", label)
    return fig, machine_bus(omega[highlight])


def kinds(columns, names):
    """Count and name of the model classes among `columns`."""
    classes = sorted({column.split("_")[0] for column in columns})
    if len(classes) == 1:
        return f"{len(columns)} {names[classes[0]]}"
    return f"{len(columns)} {names['']}"


def controllers(study, window):
    power = study.columns(r"(Tgov1|GastPti|Hygov)_.*_(PM|PMECH)")
    field = study.columns(r"(SexsPti|Ieeet1|Esdc1a)_.*_EFD")
    mask = study.window(*window)
    t = study.t[mask]
    before = study.window(0.0, study.fault[0])
    change = study.frame[power].to_numpy()
    change = (change - change[before][-1])[mask]
    efd = study.frame[field].to_numpy()[mask]

    highlight = int(np.argmax(np.abs(change).max(axis=0)))
    unit = power[highlight].split("_", 1)[1].rsplit("_", 2)[0]
    exciter = next((k for k, column in enumerate(field) if column.split("_", 1)[1].startswith(unit + "_")), 0)
    label = f"Unit at bus {machine_bus(power[highlight])}"

    governors = kinds(power, {"Tgov1": "TGOV1 governors", "GastPti": "GAST-PTI governors", "Hygov": "HYGOV governors", "": "governors"})
    exciters = kinds(field, {"SexsPti": "SEXS-PTI exciters", "Ieeet1": "IEEET1 exciters", "Esdc1a": "ESDC1A exciters", "": "exciters"})

    fig, axes = figure("Governors and exciters", f"{study.case}  /  EMT", 2)
    panels = [(change, highlight, "ΔP (pu)", f"Change in mechanical power from pre-fault · {governors}"),
              (efd, exciter, "EFD (pu)", f"Field voltage · {exciters}")]
    for ax, (values, column, unit_label, name) in zip(axes, panels):
        family(ax, t, values, column)
        time_axis(ax, study, *window, fault_ticks(study, *window))
        ax.set_ylabel(unit_label)
        heading(ax, name, label)
    return fig


def phase_labels(ax, t, values, after):
    """Letter of each phase above its first crest after time `after`."""
    span = values.max() - values.min()
    for k, (letter, colour) in enumerate(PHASES.items()):
        y = values[:, k]
        crests = np.flatnonzero((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (t[1:-1] > after) & (y[1:-1] > 0.15 * np.abs(values).max())) + 1
        if len(crests) == 0:
            continue
        crest = crests[0]
        ax.text(t[crest], y[crest] + 0.06 * span, letter, color=colour, fontsize=8.5, fontweight="bold", ha="center", va="bottom",
                bbox={"boxstyle": "round,pad=0.1", "facecolor": BACKGROUND, "edgecolor": "none", "alpha": 0.8})


def waveforms(study, window, terminal_bus):
    mask = study.window(*window)
    t = study.t[mask]
    fault_bus = study.buses[study.fault_bus]
    terminal = study.buses[terminal_bus]
    current = study.frame[[f"Switch_fault_switch_i12[{n}]" for n in (1, 2, 3)]].to_numpy()[mask]

    name = place_name(fault_bus["name"])
    fig, axes = figure(f"{name} three-phase fault", f"{study.case}  /  EMT", 3)
    panels = [
        (study.bus_voltage(study.fault_bus)[mask] / 1.0e3, "kV", f"Fault-bus voltage · bus {study.fault_bus}, {number(fault_bus['params']['kv'])} kV", window[0] + 0.004),
        (study.bus_voltage(terminal_bus)[mask] / 1.0e3, "kV", f"Generator-terminal voltage · bus {terminal_bus}, {number(terminal['params']['kv'])} kV", window[0] + 0.004),
        (current / 1.0e3, "kA", f"Fault-breaker current · bus {study.fault_bus}", study.fault[0]),
    ]
    for ax, (values, unit, title, after) in zip(axes, panels):
        for k, colour in enumerate(PHASES.values()):
            ax.plot(t, values[:, k], color=colour, linewidth=1.1)
        phase_labels(ax, t, values, after)
        time_axis(ax, study, *window, fault_ticks(study, *window))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.margins(y=0.14)
        ax.set_ylabel(unit)
        heading(ax, title)
    return fig


def main():
    case, output = sys.argv[1], sys.argv[2]
    study = Study(case, output)
    window = (study.fault[0] - 0.04, study.fault[1] + 0.03)

    figures = {"voltage": voltage(study, window)}
    figures["machines"], terminal_bus = machines(study, window)
    figures["controllers"] = controllers(study, window)
    figures["waveforms"] = waveforms(study, window, terminal_bus)
    for name, fig in figures.items():
        fig.savefig(f"{case}.{name}.png", dpi=DPI)
        plt.close(fig)
        print(f"{case}.{name}.png")


if __name__ == "__main__":
    main()
