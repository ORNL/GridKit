"""Figures of an EMT fault study in figures/: python3 ../plot.py <Case> <output.csv>, run from the study folder."""
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[3]

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
        self.reference = None
        files = sorted(Path("reference").glob(f"{case}.*.ref.csv"))
        if files:
            self.reference = pd.concat([pd.read_csv(path, index_col="t") for path in files], axis=1).sort_index()

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

    def overlay(self, ax, column, lower, upper, scale=1.0, color=EXTREME):
        if self.reference is not None and column in self.reference:
            values = self.reference[column].dropna().loc[lower - 1e-9:upper + 1e-9]
            ax.plot(values.index, values * scale, color=color, linewidth=1.4, linestyle=(0, (3, 2)), zorder=4)


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

    fig, axes = figure("Bus voltage magnitude", f"{study.case}  /  EMT", 3 if study.reference is not None else 2)
    spans = [(0.0, study.tmax, record_ticks(study), "full record"), (*window, fault_ticks(study, *window), "fault and clearing")]
    for ax, (lower, upper, ticks, name) in zip(axes, spans):
        mask = thin(study, lower, upper)
        family(ax, study.t[mask], magnitude[mask], fault)
        study.overlay(ax, columns[study.fault_bus], lower, upper, 1.0 / (1e3 * study.buses[study.fault_bus]["params"]["kv"]))
        time_axis(ax, study, lower, upper, ticks)
        ax.set_ylabel("Voltage (pu)")
        heading(ax, f"{len(buses)} buses · {name}", label)
        ax.plot(study.t[row], magnitude[row, column], "o", color=EXTREME, markersize=4, clip_on=False, zorder=5)
        offset, side = -9, "right"
        if study.t[row] - lower < 0.3 * (upper - lower):
            offset, side = 9, "left"
        ax.annotate(f"Bus {buses[column]} · {magnitude[row, column]:.2f} pu", (study.t[row], magnitude[row, column]),
                    xytext=(offset, -9), textcoords="offset points", ha=side, va="center", fontsize=8.5, color=EXTREME)
    if study.reference is not None:
        ax = axes[2]
        start = study.fault[0] - 2.0 / study.data["params"]["freq_base"]
        take = study.window(start, study.fault[0])
        for phase, color in enumerate(PHASES.values(), 1):
            key = f"Bus_{study.fault_bus}_v[{phase}]"
            ax.plot(study.t[take], study.frame.loc[take, key] / 1e3, color=color, linewidth=1.6)
            study.overlay(ax, key, start, study.fault[0], 1e-3, np.array(to_rgb(color)) * 0.55)
        phase_labels(ax, study.t[take], study.bus_voltage(study.fault_bus)[take] / 1e3, start)
        ax.margins(y=0.15)
        ax.set(xlim=(start, study.fault[0]), ylabel="Phase voltage (kV)")
        heading(ax, "Two cycles before the fault · phases A, B, C", label)
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

    fig, axes = figure("Synchronous machines", f"{study.case}  /  EMT", 3 if study.reference is not None else 2)
    panels = [(frequency, "Hz", "Rotor electrical frequency", omega[highlight], study.data["params"]["freq_base"]),
              (te, "pu", "Electromagnetic torque", torque[highlight], 1)]
    for ax, (values, unit, name, column, scale) in zip(axes, panels):
        family(ax, t, values, highlight)
        study.overlay(ax, column, *window, scale)
        time_axis(ax, study, *window, fault_ticks(study, *window))
        ax.set_ylabel(unit)
        heading(ax, f"{name} · {len(omega)} machines", label)
    if study.reference is not None:
        ax = axes[2]
        scale = study.data["params"]["freq_base"]
        ax.plot(study.t, study.frame[omega[highlight]] * scale, color=HIGHLIGHT, linewidth=1.5)
        study.overlay(ax, omega[highlight], 0, study.tmax, scale)
        time_axis(ax, study, 0, study.tmax, record_ticks(study))
        ax.set_ylabel("Hz")
        heading(ax, "Rotor electrical frequency · full record", label)
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
    ticks = fault_ticks(study, *window)
    if study.reference is not None:
        window = (0, study.tmax)
        ticks = record_ticks(study)
        keep = thin(study, *window)
        t = study.t[keep]
        panels = [(study.frame.loc[keep, power].to_numpy(), highlight, "P (pu)", f"Mechanical power · {governors}"),
                  (study.frame.loc[keep, field].to_numpy(), exciter, "EFD (pu)", f"Field voltage · {exciters}")]
    for ax, key, (values, column, unit_label, name) in zip(axes, [power[highlight], field[exciter]], panels):
        family(ax, t, values, column)
        study.overlay(ax, key, *window)
        time_axis(ax, study, *window, ticks)
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


def phases(ax, study, signal, lower, upper, scale, unit, title):
    """Phases A, B and C of `signal` from `lower` to `upper`."""
    mask = study.window(lower, upper)
    t = study.t[mask]
    values = study.frame[[f"{signal}[{n}]" for n in (1, 2, 3)]].to_numpy()[mask] * scale
    for k, colour in enumerate(PHASES.values()):
        ax.plot(t, values[:, k], color=colour, linewidth=1.1)
        study.overlay(ax, f"{signal}[{k + 1}]", lower, upper, scale, np.array(to_rgb(colour)) * 0.55)
    phase_labels(ax, t, values, lower)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.margins(y=0.14)
    ax.set_ylabel(unit)
    heading(ax, title)


def detail(time):
    """Window and ticks around a switching time."""
    lower, upper = time - 0.002, time + 0.010
    return lower, upper, list(np.round(np.arange(lower, upper + 1e-9, 0.002), 3))


def waveforms(study, window, terminal_bus):
    """Phase voltages of the fault bus and of a generator terminal."""
    fault_bus = study.buses[study.fault_bus]
    terminal = study.buses[terminal_bus]
    fig, axes = figure(f"{place_name(fault_bus['name'])} three-phase fault", f"{study.case}  /  EMT", 2)
    panels = [(study.fault_bus, f"Fault-bus voltage · bus {study.fault_bus}, {number(fault_bus['params']['kv'])} kV"),
              (terminal_bus, f"Generator-terminal voltage · bus {terminal_bus}, {number(terminal['params']['kv'])} kV")]
    for ax, (bus, title) in zip(axes, panels):
        phases(ax, study, f"Bus_{bus}_v", *window, 1e-3, "kV", title)
        time_axis(ax, study, *window, fault_ticks(study, *window))
    return fig


def transients(study):
    """Phase voltages of the fault bus as the fault is applied and as it clears."""
    fault_bus = study.buses[study.fault_bus]
    fig, axes = figure(f"{place_name(fault_bus['name'])} fault transients", f"{study.case}  /  EMT", 2)
    for ax, time, name in zip(axes, study.fault, ("fault application", "clearing")):
        lower, upper, ticks = detail(time)
        phases(ax, study, f"Bus_{study.fault_bus}_v", lower, upper, 1e-3, "kV", f"Fault-bus voltage at {name} · bus {study.fault_bus}")
        time_axis(ax, study, lower, upper, ticks)
    return fig


def currents(study, window):
    """Phase currents of the fault breaker, each pole opening at its current zero."""
    fault_bus = study.buses[study.fault_bus]
    fig, axes = figure(f"{place_name(fault_bus['name'])} fault current", f"{study.case}  /  EMT", 2)
    lower, upper, ticks = detail(study.fault[1])
    panels = [(*window, fault_ticks(study, *window), f"Fault-breaker current · bus {study.fault_bus}"),
              (lower, upper, ticks, f"Fault-breaker current at clearing · bus {study.fault_bus}")]
    for ax, (start, end, marks, title) in zip(axes, panels):
        phases(ax, study, "Switch_fault_switch_i12", start, end, 1e-3, "kA", title)
        time_axis(ax, study, start, end, marks)
    return fig


def main():
    case, output = sys.argv[1], sys.argv[2]
    study = Study(case, output)
    window = (study.fault[0] - 0.04, study.fault[1] + 0.03)

    figures = {"voltage": voltage(study, window)}
    figures["machines"], terminal_bus = machines(study, window)
    figures["controllers"] = controllers(study, window)
    figures["waveforms"] = waveforms(study, window, terminal_bus)
    figures["transients"] = transients(study)
    figures["currents"] = currents(study, window)
    Path("figures").mkdir(exist_ok=True)
    for name, fig in figures.items():
        if study.reference is not None:
            colors = (HIGHLIGHT, EXTREME)
            if name in ("waveforms", "transients", "currents"):
                colors = (INK, INK)
            handles = [Line2D([], [], color=colors[0], label="GridKit"),
                       Line2D([], [], color=colors[1], linestyle=(0, (3, 2)), label="PSCAD")]
            fig.legend(handles=handles, loc="center right", bbox_to_anchor=(0.855, fig.texts[0].get_position()[1]),
                       frameon=False, ncol=2, fontsize=8, handlelength=2.4)
        fig.savefig(f"figures/{case}.{name}.png", dpi=DPI)
        plt.close(fig)
        print(f"figures/{case}.{name}.png")


if __name__ == "__main__":
    main()
