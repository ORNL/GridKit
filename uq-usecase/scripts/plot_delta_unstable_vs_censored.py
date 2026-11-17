"""Rotor-angle deviation from COI for one unstable and one censored replay.

Inputs are replay folders holding `case.json`, `solver.json` and `replay.csv`
(DynamicSimulation output). The fault on/off times come from `solver.json`.

python plot_delta_unstable_vs_censored.py <unstable_dir> <censored_dir> --out <png>
"""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedFormatter, FixedLocator

PI = np.pi


def load_replay(run_dir):
    run_dir = Path(run_dir)
    case = json.loads((run_dir / "case.json").read_text())
    solver = json.loads((run_dir / "solver.json").read_text())
    gens = [d for d in case["devices"] if d["class"] == "Genrou"]
    ids = [g["id"] for g in gens]
    w = np.array([g["params"]["H"] * g["params"]["mva"] for g in gens])
    w = w / w.sum()
    df = pd.read_csv(run_dir / "replay.csv").drop_duplicates(subset="t", keep="last")
    df = df.sort_values("t")
    t = df["t"].to_numpy()
    delta = df[[f"Genrou_{i}_delta" for i in ids]].to_numpy()
    dev = delta - (delta @ w)[:, None]
    events = {e["type"]: e["time"] for e in solver["events"]}
    return t, dev, events["fault_on"], events["fault_off"]


def fmt_pi(v, _):
    k = round(v / (PI / 2))
    if abs(v - k * PI / 2) > 1e-6:
        return ""
    sign = "-" if k < 0 else ""
    m = abs(k)
    if m == 0:
        return "0"
    if m == 2:
        return f"{sign}π"
    if m == 4:
        return f"{sign}2π"
    if m == 1:
        return f"{sign}π/2"
    return f"{sign}{m}π/2"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("unstable_dir")
    p.add_argument("censored_dir")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    ticks = [k * PI / 2 for k in range(-4, 5)]
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    for ax, run_dir in zip(axes, [args.unstable_dir, args.censored_dir]):
        t, dev, t_on, t_off = load_replay(run_dir)
        ax.plot(t, dev, color="0.75", lw=0.7)
        mx = np.abs(dev).max(axis=1)
        ax.plot(t, mx, color="tab:red", lw=1.8,
                label=r"max$_i$ $|\delta_i-\delta_{COI}|$")
        ax.axhline(PI, color="k", ls="--", lw=1)
        ax.axhline(-PI, color="k", ls="--", lw=1)
        ax.axvspan(t_on, t_off, color="orange", alpha=0.15, label="fault on")
        ax.set_ylim(-2.5 * PI, 2.5 * PI)
        ax.yaxis.set_major_locator(FixedLocator(ticks))
        ax.yaxis.set_major_formatter(FixedFormatter([fmt_pi(v, None) for v in ticks]))
        ax.set_ylabel("rad")
        ax.set_xlim(0.9, 3.0)
        ax.tick_params(labelsize=12)
        ax.legend(loc="upper left", fontsize=11)

        if (mx > PI).any():
            cross = t[np.argmax(mx > PI)]
            ax.axvline(cross, color="tab:red", ls=":", lw=1)
            ax.text(cross + 0.05, -2.0 * PI, "pole slip: passes π and keeps growing",
                    color="tab:red", fontsize=11)
        else:
            ax.text(1.6, -2.0 * PI, "stays below π: recovers", color="tab:green", fontsize=11)

    axes[1].set_xlabel("time (s)", fontsize=12)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print("saved", args.out)


if __name__ == "__main__":
    main()
