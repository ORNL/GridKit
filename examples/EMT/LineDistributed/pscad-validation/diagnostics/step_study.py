"""Which conditions let GridKit's adaptive step exceed the line delay in energized operation?

Variants of the 69 kV wood-pole line (length, conductor losses, soil resistivity)
are fitted with the same pipeline as the validation lines and simulated in the
``sag`` scenario with two terminations: the R-L load of the validation suite and
a purely resistive load of the same resistance. The accepted step over the
shortest modal delay is then summarized in the pre-event steady state (0 to 20
ms) and after the recovery (60 to 100 ms).

Usage, with the GridWorkbench environment:

    python diagnostics/step_study.py fit      # geometry variants, fits, run files
    python diagnostics/step_study.py report   # after running runs/*.json with line_distributed

The GridKit runs happen in the container:

    docker exec gridkit-ld bash -lc 'cd /work/examples/EMT/LineDistributed/pscad-validation/diagnostics \
      && for f in runs/*.json; do /home/gridkit/build/examples/EMT/LineDistributed/line_distributed $f results; done'
"""

import copy
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import make_lines  # noqa: E402
from scenarios import scenarios  # noqa: E402

LINES = HERE / "lines"
RUNS = HERE / "runs"
RESULTS = HERE / "results"
BASE = json.loads((ROOT / "lines" / "69kv-wood-pole.line.json").read_text())
#: name -> (length_m, conductor conductivity factor, earth conductivity S/m)
VARIANTS = {
    "l2": (2e3, 1.0, 0.01),
    "l4": (4e3, 1.0, 0.01),
    "l8": (8e3, 1.0, 0.01),
    "l16": (16e3, 1.0, 0.01),
    "l50": (50e3, 1.0, 0.01),
    "l100": (100e3, 1.0, 0.01),
    "l200": (200e3, 1.0, 0.01),
    "l16_lossy4": (16e3, 0.25, 0.01),
    "l100_lossy4": (100e3, 0.25, 0.01),
    "l16_soil1000": (16e3, 1.0, 0.001),
    "l100_soil1000": (100e3, 1.0, 0.001),
}
#: Terminations. Loads: the suite's series R-L (open at high frequency) or the same
#: resistance alone. Sources: the suite's X/R = 10 branch (inductive at high
#: frequency) or a resistor of the same short-circuit impedance magnitude.
LOADS = ("rl", "r")
SOURCES = ("xr10", "r")
WINDOWS = {"steady 0-20 ms": (0.0, 0.02), "recovered 60-100 ms": (0.06, 0.1)}
STYLE = {("xr10", "rl"): ("#D55E00", "-"), ("xr10", "r"): ("#0072B2", "-"), ("r", "rl"): ("#D55E00", ":"), ("r", "r"): ("#0072B2", ":")}


def geometry(name, length, loss_factor, earth):
    doc = copy.deepcopy(BASE)
    doc["name"] = f"{BASE['name']} variant {name}"
    doc["include"] = []
    catalog = json.loads((ROOT / "lines" / "catalog.json").read_text())
    conductors = {}
    for key in {c["type"] for c in doc["conductors"]}:
        entry = copy.deepcopy(catalog["conductors"][key])
        entry["conductivity"] *= loss_factor
        conductors[key] = entry
    doc["catalog"] = {"conductors": conductors, "towers": {doc["tower"]: catalog["towers"][doc["tower"]]}}
    doc["path"]["length"] = length
    doc["earth"]["conductivity"] = earth
    return doc


def fit():
    LINES.mkdir(exist_ok=True)
    RUNS.mkdir(exist_ok=True)
    for name, (length, loss, earth) in VARIANTS.items():
        (LINES / f"{name}.line.json").write_text(json.dumps(geometry(name, length, loss, earth), indent=2) + "\n")
    make_lines.LINES = LINES
    make_lines.main(list(VARIANTS))
    for name in VARIANTS:
        sequence = json.loads((LINES / name / "sequence.json").read_text())
        runs = scenarios("69kv-wood-pole", sequence, "fd")
        for source in SOURCES:
            for load in LOADS:
                spec = copy.deepcopy(runs["sag"])
                spec.pop("description", None)
                spec["name"] = f"{name}__{source}__{load}"
                spec["line"] = f"../lines/{name}/fd.fit.json"
                if load == "r":
                    spec["load"]["inductance"] = [0.0, 0.0, 0.0]
                if source == "r":
                    r, l = spec["source"]["resistance"][0], spec["source"]["inductance"][0]
                    magnitude = float(np.hypot(r, 2 * np.pi * 60.0 * l))
                    spec["source"]["resistance"] = [magnitude] * 3
                    spec["source"]["inductance"] = [0.0] * 3
                (RUNS / f"{spec['name']}.json").write_text(json.dumps(spec, indent=2) + "\n")
    print("fits and run files written; now run the GridKit cases in the container")


def analyse(path):
    steps = np.genfromtxt(path / "steps.csv", delimiter=",", names=True)
    data = json.loads((path / "input.json").read_text())
    tau = min(d for d in data["line"]["H"]["delays"] if d > 0)
    t, h = steps["t"], steps["h"]
    result = {"tau_s": tau, "accepted_steps": int(len(h)), "max_h_over_tau": float(h.max() / tau)}
    for label, (lo, hi) in WINDOWS.items():
        start = np.maximum(t - h, lo)
        end = np.minimum(t, hi)
        duration = np.clip(end - start, 0, None)
        keep = duration > 0
        hw, dw = h[keep], duration[keep]
        result[label] = {"median_h_over_tau": float(np.median(hw) / tau), "max_h_over_tau": float(hw.max() / tau),
                         "time_fraction_above_tau": float(dw[hw > tau].sum() / dw.sum()), "steps": int(keep.sum())}
    return result, t, h, tau


def report():
    summary = {}
    names = [n for n in VARIANTS if any((RESULTS / f"{n}__{s}__{l}" / "steps.csv").exists() for s in SOURCES for l in LOADS)]
    fig, axes = plt.subplots(len(names), 1, figsize=(13, 2.8 * len(names)), layout="constrained", squeeze=False)
    for ax, name in zip(axes[:, 0], names):
        length, loss, earth = VARIANTS[name]
        tau = None
        for source in SOURCES:
            for load in LOADS:
                path = RESULTS / f"{name}__{source}__{load}"
                if not (path / "steps.csv").exists():
                    continue
                stats, t, h, tau = analyse(path)
                summary[f"{name}__{source}__{load}"] = stats
                color, style = STYLE[(source, load)]
                ax.semilogy(t * 1e3, h / tau, color=color, ls=style, lw=0.7,
                            label=f"source {'R-L' if source == 'xr10' else 'R'}, load {'R-L' if load == 'rl' else 'R'}: median h/tau "
                                  f"{stats['steady 0-20 ms']['median_h_over_tau']:.2f} / {stats['recovered 60-100 ms']['median_h_over_tau']:.2f}, "
                                  f"{stats['accepted_steps']} steps")
        ax.axhline(1, color="black", ls="--", lw=1)
        ax.axvspan(20, 50, color="0.6", alpha=0.12)
        ax.set(ylabel=r"$h/\tau$", ylim=(1e-4, 1e3))
        ax.set_title(f"{name}: {length / 1e3:g} km, conductor conductivity x{loss:g}, soil {1 / earth:g} ohm m, tau = {tau * 1e6:.0f} us", fontsize=10)
        ax.grid(color="0.9", lw=0.5)
        ax.legend(frameon=False, fontsize=7, loc="lower right", ncol=2)
    axes[-1, 0].set_xlabel("Study time (ms)")
    fig.suptitle("69 kV wood-pole variants, FD line, sag scenario: accepted step over the shortest modal delay "
                 "(medians: 0-20 ms steady state / 60-100 ms after recovery)")
    fig.savefig(HERE / "step_study.png", dpi=130)
    (HERE / "step_study.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"{'variant':16s} {'source':6s} {'load':4s} {'tau us':>7s} {'steps':>7s} {'median h/tau 0-20':>18s} {'max 0-20':>9s} {'%t>tau 0-20':>12s} {'median 60-100':>14s} {'max 60-100':>11s} {'%t>tau 60-100':>14s}")
    for key, s in summary.items():
        name, source, load = key.split("__")
        a, b = s["steady 0-20 ms"], s["recovered 60-100 ms"]
        print(f"{name:16s} {source:6s} {load:4s} {s['tau_s'] * 1e6:7.1f} {s['accepted_steps']:7d} {a['median_h_over_tau']:18.3f} {a['max_h_over_tau']:9.2f} {100 * a['time_fraction_above_tau']:12.1f} "
              f"{b['median_h_over_tau']:14.3f} {b['max_h_over_tau']:11.2f} {100 * b['time_fraction_above_tau']:14.1f}")


if __name__ == "__main__":
    {"fit": fit, "report": report}[sys.argv[1]]()
