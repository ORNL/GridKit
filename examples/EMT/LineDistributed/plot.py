"""Plot the EMT demo. Figures always go beside this source script, in plots/."""

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/gridkit-line-distributed-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
COLORS = ["#0072B2", "#D55E00", "#009E73"]
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": .18, "svg.fonttype": "none",
                     "svg.hashsalt": "line-distributed"})


def save(fig, name, destination):
    fig.savefig(destination / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    fig.savefig(destination / f"{name}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def events(axes, data):
    for ax in np.asarray(axes).flat:
        for event in data["events"]:
            if event["on"] < data["solver"]["tf"]:
                ax.axvspan(event["on"]*1000, min(event["off"], data["solver"]["tf"])*1000,
                           color="0.5", alpha=.08, zorder=0)
        ax.set_xlabel("Time (ms)")
        ax.set_xlim(data["solver"]["t0"]*1000, data["solver"]["tf"]*1000)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, nargs="?", default=HERE / "results")
    args = parser.parse_args()
    destination = HERE / "plots"
    if "build" in destination.parts:
        parser.error("Run the plot.py in examples/EMT/LineDistributed; figures belong in the source example directory")
    destination.mkdir(exist_ok=True)
    runs = []
    for path in sorted(args.results.iterdir()):
        if not (path / "statistics.json").exists():
            continue
        data = json.loads((path / "input.json").read_text())
        stats = json.loads((path / "statistics.json").read_text())
        samples = np.genfromtxt(path / "samples.csv", delimiter=",", names=True)
        steps = np.genfromtxt(path / "steps.csv", delimiter=",", names=True)
        t = samples["t"]*1000
        fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), layout="constrained")
        for end, ax in enumerate(axes[0], 1):
            for color, phase in zip(COLORS, "abc"):
                ax.plot(t, samples[f"v{end}_{phase}"]/1000, color=color, lw=1.0, label=phase.upper())
            ax.set_title(f"({'a' if end == 1 else 'b'}) {'Sending' if end == 1 else 'Receiving'} voltage")
            ax.set_ylabel("Voltage (kV)")
            ax.legend(ncol=3, loc="upper right", frameon=False)
        ax = axes[1, 0]
        ax.plot(t, samples["i1_a"], color=COLORS[0], lw=1, label="Sending end")
        ax.plot(t, samples["i2_a"], color=COLORS[1], lw=1, label="Receiving end")
        ax.set_title("(c) Phase-A line current injections")
        ax.set_ylabel("Current into bus (A)")
        ax.legend(frameon=False)
        ax = axes[1, 1]
        if steps.size:
            ax.semilogy(steps["t"]*1000, steps["h_over_tau"], color=COLORS[0], lw=.8)
            ax.axhline(1, color="0.3", ls="--", lw=.8)
        ax.set_title("(d) Accepted steps / shortest positive delay")
        ax.set_ylabel(r"$h/\tau_{min}$ (dimensionless)")
        events(axes, data)
        fig.suptitle(f"{path.name.replace('_', ' ').title()}  |  {stats['unknowns']} unknowns, {stats['steps']:,} accepted steps")
        save(fig, path.name, destination)
        runs.append((path, data, samples))

    if not runs:
        raise RuntimeError("No completed simulations to plot")
    rows = (len(runs) + 2)//3
    fig, axes = plt.subplots(rows, 3, figsize=(13, 2.5*rows), layout="constrained", squeeze=False)
    for ax, (path, data, samples) in zip(axes.flat, runs):
        for color, phase in zip(COLORS, "abc"):
            ax.plot(samples["t"]*1000, samples[f"v2_{phase}"]/1000, color=color, lw=.8, label=phase.upper())
        events([ax], data)
        ax.set_title(path.name.replace("_", " ").title(), fontsize=10)
        ax.set_ylabel("Receiving voltage (kV)")
    for ax in list(axes.flat)[len(runs):]:
        ax.set_visible(False)
    fig.suptitle("LineDistributed: three-phase source, distributed line, impedance load")
    save(fig, "overview", destination)

    path = args.results / "energization"
    if (path / "samples.csv").exists():
        data = json.loads((path / "input.json").read_text())
        samples = np.genfromtxt(path / "samples.csv", delimiter=",", names=True)
        tau = data["line"]["H"]["delay"]
        on = data["events"][0]["on"]
        fig, axes = plt.subplots(2, 1, figsize=(10, 6), layout="constrained", sharex=True)
        axes[0].plot(samples["t"]*1000, samples["v1_a"]/1000, color=COLORS[0], label="Sending voltage")
        axes[0].plot(samples["t"]*1000, samples["v2_a"]/1000, color=COLORS[1], label="Receiving voltage")
        reference = path / "reference.csv"
        if reference.exists():
            ref = np.genfromtxt(reference, delimiter=",", names=True)
            axes[0].plot(ref["t"]*1000, ref["v2_a"]/1000, color="0.15", ls="--", lw=.8, label="Exact characteristic reference")
        axes[0].set_ylabel("Phase-A voltage (kV)")
        axes[0].set_title("(a) Energization and load reflections")
        axes[1].plot(samples["t"]*1000, samples["ref1_1"], color=COLORS[0], label="Launched at sending end")
        axes[1].plot(samples["t"]*1000, samples["inc2_1"], color=COLORS[1], label="Arriving at receiving end")
        axes[1].set_ylabel("Wave current (A)")
        axes[1].set_xlabel("Time (ms)")
        axes[1].set_title("(b) Actual wave delay, including subsequent reflections")
        for ax in axes:
            ax.set_xlim((on-.0001)*1000, (on+7*tau)*1000)
            for arrival in (on+tau, on+3*tau, on+5*tau):
                ax.axvline(arrival*1000, color="0.5", ls=":", lw=.8)
            ax.legend(frameon=False, fontsize=9)
        visible = (samples["t"] >= on-.0001) & (samples["t"] <= on+7*tau)
        for ax, columns, scale in ((axes[0], ("v1_a", "v2_a"), 1000),
                                   (axes[1], ("ref1_1", "inc2_1"), 1)):
            low = min(samples[c][visible].min()/scale for c in columns)
            high = max(samples[c][visible].max()/scale for c in columns)
            margin = .08*max(high-low, 1)
            ax.set_ylim(low-margin, high+margin)
        save(fig, "travelling-waves", destination)
    print(f"Saved PNG and editable SVG figures in {destination}")


if __name__ == "__main__":
    main()
