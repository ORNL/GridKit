#!/usr/bin/env python3
"""Plot the WECC240Relay line trip next to this script.

Usage: python3 plot.py [run_dir]
  run_dir holds WECC240Relay.case.json, WECC240Relay.solver.json, and the
  WECC240Relay.csv output; it defaults to the current directory.
"""

import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
OUT = Path(__file__).resolve().parent / "WECC240Relay.png"

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT2 = "#52514e"
GRID = "#e6e5e1"
MUTED = "#c9c8c3"
SLOTS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")

F0 = 60.0
FREQ_FILTER = 0.05  # s, first-order filter on d(theta)/dt
ZOOM = (0.97, 1.25)
WIDE = (0.5, 6.0)

case = json.loads((RUN / "WECC240Relay.case.json").read_text())
solver = json.loads((RUN / "WECC240Relay.solver.json").read_text())
with (RUN / "WECC240Relay.csv").open() as stream:
    header = stream.readline().strip().split(",")
data = np.loadtxt(RUN / "WECC240Relay.csv", delimiter=",", skiprows=1)
keep = np.append(np.diff(data[:, 0]) > 0, True)  # drop the duplicate sample at the event restart
data = data[keep]
t = data[:, 0]
column = {name: i for i, name in enumerate(header)}

# Bus names repeat in WECC240, so bus columns are matched by case order.
buses = [b for b in case["buses"] if b["mon"]]
name = {b["number"]: b["name"] for b in case["buses"]}
bus_vm = dict(zip((b["number"] for b in buses), (i for i, h in enumerate(header) if h.startswith("Bus_") and h.endswith("_Vm"))))
bus_va = dict(zip((b["number"] for b in buses), (i for i, h in enumerate(header) if h.startswith("Bus_") and h.endswith("_Va"))))

devices = case["devices"]
relays = [d for d in devices if d["class"] == "OvercurrentRelay"]
faults = [d for d in devices if d["class"] == "BusFault"]
fault_event = next(e for e in solver["events"] if e["type"] == "fault_on")
fault_bus = faults[fault_event["element_id"]]["ports"]["bus"]
fault_time = fault_event["time"]
primary_ttrip = min(r["params"]["Ttrip"] for r in relays)


def series(label):
    return data[:, column[label]]


def at(values, time):
    return values[np.argmin(np.abs(t - time))]


def relay_location(relay):
    bus, remote = (int(n) for n in relay["id"].split("_")[1:3])
    role = "primary"
    if relay["params"]["Ttrip"] > primary_ttrip:
        role = "backup"
    toward = f" (to {name[remote]})"
    if remote == fault_bus:
        toward = ""
    return bus, f"{name[bus]}{toward}, {role}"


def bus_frequency(bus):
    theta = np.unwrap(data[:, bus_va[bus]])
    raw = F0 + np.gradient(theta, t) / (2.0 * np.pi)
    out = np.empty_like(raw)
    out[0] = raw[0]
    for k in range(1, len(raw)):
        dt = t[k] - t[k - 1]
        out[k] = out[k - 1] + dt / (FREQ_FILTER + dt) * (raw[k] - out[k - 1])
    return out


primary = relays[0]["id"]
trip_signal = series(f"OvercurrentRelay_{primary}_trip")
k = np.argmax((t >= fault_time) & (trip_signal >= 0.5))
trip_time = t[k - 1] + (0.5 - trip_signal[k - 1]) / (trip_signal[k] - trip_signal[k - 1]) * (t[k] - t[k - 1])
breaker = next(d for d in devices if d["class"] == "BranchBreakers" and fault_bus in d["ports"].values())
events = [(fault_time, "fault"), (trip_time, "trip"), (trip_time + breaker["params"]["Tbrk"], "open")]
locations = [relay_location(r) for r in relays]

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.edgecolor": GRID,
        "axes.labelcolor": TEXT2,
        "xtick.color": TEXT2,
        "ytick.color": TEXT2,
        "axes.titlecolor": TEXT,
        "axes.titleweight": "bold",
        "axes.titlesize": 11,
        "axes.titlelocation": "left",
        "axes.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "legend.labelcolor": TEXT2,
        "lines.linewidth": 1.6,
    }
)

fig, axes = plt.subplots(3, 2, figsize=(13, 12), constrained_layout=True)
ax_ct, ax_x, ax_v, ax_f, ax_i, ax_w = axes.ravel()


def mark_events(ax, label=False):
    for time, text in events:
        ax.axvline(time, color=TEXT2, linestyle=":", linewidth=1.0)
        if label:
            ax.annotate(text, (time, 0.98), xycoords=("data", "axes fraction"), xytext=(3, 0),
                        textcoords="offset points", rotation=90, ha="left", va="top", color=TEXT2, fontsize=9)


def note(ax, text, xy):
    ax.annotate(text, xy, xytext=(-3, 3), textcoords="offset points", ha="right", va="bottom", color=TEXT2, fontsize=9)


# Relay CT currents normalized by pickup; each location keeps its color in every panel
for slot, (relay, (_, label)) in enumerate(zip(relays, locations)):
    params = relay["params"]
    ax_ct.plot(t, series(f"OvercurrentRelay_{relay['id']}_im") / params["Ipickup"], color=SLOTS[slot],
               label=f"{label} {params['Ipickup']:g} p.u. / {params['Ttrip']:g} s")
ax_ct.axhline(1.0, color=TEXT2, linestyle="--", linewidth=1.0)
note(ax_ct, "pickup", (ZOOM[1], 1.0))
mark_events(ax_ct, label=True)
ax_ct.set(xlim=ZOOM, ylim=(0, 3.4), title="Relay CT current / pickup", ylabel="|I| / Ipickup (-)")
ax_ct.legend(loc="upper right", bbox_to_anchor=(1.0, 0.86), fontsize=8.5, handlelength=1.5)

# Relay latches; equal settings give coincident curves, so every second relay is dashed
for slot, relay in enumerate(relays):
    ax_x.plot(t, series(f"OvercurrentRelay_{relay['id']}_x"), color=SLOTS[slot], linestyle=["-", (0, (4, 3))][slot % 2])
for level, text in [(0.75, "trip level x = 3/4"), (0.5, "commit x = 1/2")]:
    ax_x.axhline(level, color=TEXT2, linestyle="--", linewidth=1.0)
    note(ax_x, text, (ZOOM[1], level))
backups = [r for r in relays if r["params"]["Ttrip"] > primary_ttrip]
peak = max(series(f"OvercurrentRelay_{r['id']}_x").max() for r in backups)
ax_x.annotate(f"backups peak at x = {peak:.2f}, then reset", (events[2][0] + 0.005, peak), xytext=(8, 10),
              textcoords="offset points", color=TEXT2, fontsize=9,
              arrowprops={"arrowstyle": "-", "color": TEXT2, "linewidth": 0.8})
ax_x.annotate("primaries lock out", (1.15, 1.0), xytext=(0, -14), textcoords="offset points", color=TEXT2, fontsize=9)
mark_events(ax_x)
ax_x.set(xlim=ZOOM, ylim=(-0.03, 1.05), title="Relay lockout latches", ylabel="x (-)")

# Bus voltages: every bus in gray, relay locations and the fault bus highlighted
for bus in bus_vm:
    ax_v.plot(t, data[:, bus_vm[bus]], color=MUTED, linewidth=0.6, zorder=1)
ax_v.plot(t, data[:, bus_vm[fault_bus]], color=TEXT2, linestyle="--", linewidth=1.2,
          label=f"{name[fault_bus]} ({fault_bus}), fault", zorder=2)
for slot, (bus, _) in enumerate(locations):
    ax_v.plot(t, data[:, bus_vm[bus]], color=SLOTS[slot], label=f"{name[bus]} ({bus})", zorder=3)
ax_v.plot([], [], color=MUTED, linewidth=0.6, label=f"other {len(bus_vm) - len(locations) - 1} buses")
mark_events(ax_v)
ax_v.set(xlim=WIDE, ylim=(0, None), title="Bus voltage magnitudes", ylabel="Vm (p.u.)")
ax_v.legend(loc="lower right", fontsize=9)

# Bus frequencies from the filtered angle derivative; the dead fault bus has no frequency
for bus in bus_va:
    if bus != fault_bus:
        ax_f.plot(t, bus_frequency(bus), color=MUTED, linewidth=0.6, zorder=1)
for slot, (bus, _) in enumerate(locations):
    ax_f.plot(t, bus_frequency(bus), color=SLOTS[slot], label=f"{name[bus]} ({bus})", zorder=3)
ax_f.plot([], [], color=MUTED, linewidth=0.6, label="other buses")
mark_events(ax_f)
ax_f.set(xlim=WIDE, title=f"Bus frequencies (dθ/dt, {FREQ_FILTER * 1e3:g} ms filter)", ylabel="frequency (Hz)")
ax_f.legend(loc="upper right", fontsize=9, ncol=2)

# Branch-current change from pre-fault: the cleared line and the monitored lines
pre = fault_time - 0.01
cleared = series(f"OvercurrentRelay_{primary}_im")
ax_i.plot(t, cleared - at(cleared, pre), color=TEXT2, linestyle="--", linewidth=1.2,
          label=f"faulted line (cleared): {at(cleared, pre):.1f} → 0")
lines = [d for d in devices if d["class"] == "Branch" and d.get("mon")]
for slot, line in enumerate(lines):
    current = series(f"Branch_{line['id']}_im1")
    ax_i.plot(t, current - at(current, pre), color=SLOTS[4 + slot],
              label=f"{name[line['ports']['bus1']]}-{name[line['ports']['bus2']]}: {at(current, pre):.1f} → {current[-1]:.1f}")
ax_i.axhline(0.0, color=TEXT2, linewidth=0.8)
ax_i.annotate("cleared-line fault current off scale", (events[2][0] + 0.005, 6.5), xytext=(4, -4),
              textcoords="offset points", ha="left", va="top", color=TEXT2, fontsize=9)
mark_events(ax_i)
ax_i.set(xlim=WIDE, ylim=(-6.5, 6.5), title="Branch-current change from pre-fault",
         ylabel="ΔI, bus-1 end (p.u.)", xlabel="time (s)")
ax_i.legend(loc="upper right", bbox_to_anchor=(1.0, 0.92), fontsize=9)

# Generator speeds in Hz; the units at the JOHN DAY plant are highlighted
speeds = {d["id"]: d["ports"]["bus"] for d in devices if d["class"] == "Genrou"}
plant = next(d["ports"]["bus2"] for d in devices if "Branch" in d["class"] and d["ports"]["bus1"] == locations[0][0]
             and d["ports"]["bus2"] in speeds.values())
for unit, bus in speeds.items():
    highlight = bus == plant
    ax_w.plot(t, F0 * (1.0 + series(f"Genrou_{unit}_omega")), color=[MUTED, SLOTS[0]][highlight],
              linewidth=[0.6, 1.6][highlight], zorder=[1, 3][highlight])
units = sum(bus == plant for bus in speeds.values())
ax_w.plot([], [], color=SLOTS[0], label=f"{name[plant]} plant ({plant}), {units} coincident units")
ax_w.plot([], [], color=MUTED, linewidth=0.6, label=f"other {len(speeds) - units} GENROU units")
mark_events(ax_w)
ax_w.set(xlim=WIDE, title="Generator speeds", ylabel="speed (Hz)", xlabel="time (s)")
ax_w.legend(loc="upper right", fontsize=9)

fig.suptitle(f"WECC240: fault at {name[fault_bus]} at {fault_time:g} s, cleared by relays at both line ends",
             color=TEXT, fontweight="bold", fontsize=13, x=0.01, ha="left")
fig.savefig(OUT, dpi=150)
print(OUT)
