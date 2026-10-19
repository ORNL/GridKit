"""Plot GridKit's accepted step over the line propagation delay against time.

For every ``results/gridkit/<line>__<model>/<scenario>/steps.csv`` the accepted step
h is divided by the shortest modal delay (solid) and the longest (dotted) and drawn
against time with the h/tau = 1 line. One figure per line and model goes to
``plots/steps/<line>_<model>.png``; ``plots/steps/summary.json`` lists, per run, the
maximum h/tau, the share of steps and of elapsed time with h above the shortest and
the longest delay, and the intervals where h exceeded the shortest delay.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "lib"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from scenarios import CODES  # noqa: E402

GRIDKIT = HERE / "results" / "gridkit"
OUT = HERE / "plots" / "steps"
SCENARIOS = ("energize", "pulse", "fault", "sag")


def intervals(t, h, above):
    """Contiguous (start, end) intervals of accepted steps with ``above`` true."""
    result = []
    start = None
    for k in range(len(t)):
        if above[k] and start is None:
            start = t[k] - h[k]
        if not above[k] and start is not None:
            result.append((float(start), float(t[k - 1])))
            start = None
    if start is not None:
        result.append((float(start), float(t[-1])))
    return result


def analyse(path):
    steps = np.genfromtxt(path / "steps.csv", delimiter=",", names=True)
    data = json.loads((path / "input.json").read_text())
    delays = np.array(data["line"]["H"]["delays"])
    delays = delays[delays > 0]
    tau_min, tau_max = float(delays.min()), float(delays.max())
    t, h = steps["t"], steps["h"]
    duration = np.minimum(t, data["solver"]["tf"]) - np.maximum(t - h, data["solver"]["t0"])
    keep = duration > 0
    t, h, duration = t[keep], h[keep], duration[keep]
    above_min, above_max = h > tau_min, h > tau_max
    return {
        "tau_min_s": tau_min, "tau_max_s": tau_max, "accepted_steps": int(len(h)),
        "max_h_s": float(h.max()), "median_h_s": float(np.median(h)),
        "max_h_over_tau_min": float(h.max() / tau_min), "max_h_over_tau_max": float(h.max() / tau_max),
        "median_h_over_tau_min": float(np.median(h) / tau_min),
        "steps_above_tau_min": int(above_min.sum()), "steps_above_tau_max": int(above_max.sum()),
        "step_fraction_above_tau_min": float(above_min.mean()),
        "time_fraction_above_tau_min": float(duration[above_min].sum() / duration.sum()),
        "time_fraction_above_tau_max": float(duration[above_max].sum() / duration.sum()),
        "intervals_above_tau_min_s": intervals(t, h, above_min),
        "events": [(e["type"], e["on"], e["off"]) for e in data["events"]],
    }, t, h


def main(argv):
    models = argv or ["fd"]
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {}
    for line in CODES:
        for model in models:
            folder = GRIDKIT / f"{line}__{model}"
            runs = [s for s in SCENARIOS if (folder / s / "steps.csv").exists()]
            if not runs:
                continue
            fig, axes = plt.subplots(len(runs), 1, figsize=(13, 2.9 * len(runs)), layout="constrained", squeeze=False)
            for ax, scenario in zip(axes[:, 0], runs):
                stats, t, h = analyse(folder / scenario)
                summary[f"{line}/{model}/{scenario}"] = stats
                ax.semilogy(t * 1e3, h / stats["tau_min_s"], color="#0072B2", lw=0.7, label=r"$h/\tau_{min}$")
                if stats["tau_max_s"] > 1.001 * stats["tau_min_s"]:
                    ax.semilogy(t * 1e3, h / stats["tau_max_s"], color="#D55E00", lw=0.6, alpha=0.7, ls=":", label=r"$h/\tau_{max}$")
                ax.axhline(1.0, color="black", ls="--", lw=1.0)
                for kind, on, off in stats["events"]:
                    ax.axvspan(on * 1e3, min(off, t[-1]) * 1e3, color="0.6", alpha=0.12)
                for start, end in stats["intervals_above_tau_min_s"]:
                    ax.axvspan(start * 1e3, end * 1e3, color="#009E73", alpha=0.25)
                ax.set_ylabel(r"$h/\tau$")
                ax.set_title(f"{scenario}: {stats['accepted_steps']} steps, max h/tau(min) = {stats['max_h_over_tau_min']:.2f}, "
                             f"median {stats['median_h_over_tau_min']:.3f}, time with h > tau(min): {stats['time_fraction_above_tau_min']:.1%}",
                             fontsize=10)
                ax.set_ylim(1e-5, 1e2)
                ax.grid(color="0.9", lw=0.5)
                ax.legend(frameon=False, fontsize=8, loc="upper right")
            axes[-1, 0].set_xlabel("Study time (ms)")
            fig.suptitle(f"{line} / {model}: accepted GridKit step over the modal propagation delays "
                         f"(tau min {summary[f'{line}/{model}/{runs[0]}']['tau_min_s'] * 1e6:.1f} us, max {summary[f'{line}/{model}/{runs[0]}']['tau_max_s'] * 1e6:.1f} us); "
                         "shaded grey: events, shaded green: h above tau(min)")
            fig.savefig(OUT / f"{line}_{model}.png", dpi=130)
            plt.close(fig)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"{'run':40s} {'steps':>8s} {'max h/tau_min':>14s} {'median':>8s} {'%time>tau_min':>14s} {'%time>tau_max':>14s}  intervals above tau_min")
    for key, s in summary.items():
        spans = ", ".join(f"{a * 1e3:.2f}-{b * 1e3:.2f} ms" for a, b in s["intervals_above_tau_min_s"][:6])
        more = f" (+{len(s['intervals_above_tau_min_s']) - 6} more)" if len(s["intervals_above_tau_min_s"]) > 6 else ""
        print(f"{key:40s} {s['accepted_steps']:8d} {s['max_h_over_tau_min']:14.3f} {s['median_h_over_tau_min']:8.4f} "
              f"{100 * s['time_fraction_above_tau_min']:14.1f} {100 * s['time_fraction_above_tau_max']:14.1f}  {spans}{more}")


if __name__ == "__main__":
    main(sys.argv[1:])
