#!/usr/bin/env python3
"""Plot every element of the propagation fitting targets Gin and Gout.

Reproduces the UniversalLineModel construction exactly from a sweep that
monitors Tv, Ti, H, and Tau: per mode the delay is the minimum of the
monitored tau trace, the backwound propagation is
hmps_m = H_m exp(+j omega tau_m), and the factors are

    Gin[row, col]  = hmps_row * Tv[col, row]
    Gout[row, col] = conj(Ti[row, col])

Each factor gets one figure with magnitude, unwrapped phase, real, and
imaginary panels for all K x K elements, written to output/<line>/. These
are the raw fitter inputs, so the rationality of the data can be judged
directly.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gallery import line_model_file

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output"
BUILD = HERE.parents[2] / "build"
DEFAULT_SWEEP_APP = BUILD / "application" / "EMT" / "FrequencyResponse" / "FrequencyResponse"

FREQUENCY = {"start": 10.0, "stop": 1.0e8, "points": 401, "scale": "log"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--line", default="345kv-horizontal")
    parser.add_argument("--sweep-app", type=Path, default=DEFAULT_SWEEP_APP)
    return parser.parse_args()


def run_sweep(name: str, line_dir: Path, sweep_app: Path) -> Path:
    # Paths inside the spec resolve relative to the spec file in
    # output/<line>/.
    spec = {
        "model": str(line_model_file(name)),
        "frequency": FREQUENCY,
        "output": {
            "file": "factors.csv",
            "variables": ["Tv", "Ti", "H", "Tau"],
        },
    }
    spec_path = line_dir / "factors.solver.json"
    spec_path.write_text(json.dumps(spec, indent=4) + "\n")
    subprocess.run([str(sweep_app), str(spec_path)],
                   check=True, cwd=HERE, stdout=subprocess.DEVNULL)
    return line_dir / "factors.csv"


def read_factors(path: Path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        names = reader.fieldnames or []
    omega = np.array([float(r["omega"]) for r in rows])
    m = len(rows)

    k = 0
    while f"Overhead_H_real_{k}" in names:
        k += 1
    if k == 0:
        raise ValueError(f"{path} lacks Overhead_H columns")

    def matrix(prefix: str) -> np.ndarray:
        out = np.zeros((m, k, k), complex)
        for s, r in enumerate(rows):
            for i in range(k):
                for j in range(k):
                    out[s, i, j] = complex(
                        float(r[f"{prefix}_real_{i}_{j}"]),
                        float(r[f"{prefix}_imag_{i}_{j}"]))
        return out

    h = np.array([[complex(float(r[f"Overhead_H_real_{q}"]),
                           float(r[f"Overhead_H_imag_{q}"]))
                   for q in range(k)] for r in rows])
    tau = np.array([[float(r[f"Overhead_Tau_{q}"]) for q in range(k)]
                    for r in rows])
    return omega, k, matrix("Overhead_Tv"), matrix("Overhead_Ti"), h, tau


def build_targets(omega, k, tv, ti, h, tau):
    delays = tau.min(axis=0)
    hmps = h * np.exp(1j * omega[:, None] * delays[None, :])

    gin = np.zeros((len(omega), k, k), complex)
    gout = np.zeros((len(omega), k, k), complex)
    for row in range(k):
        for col in range(k):
            gin[:, row, col] = hmps[:, row] * tv[:, col, row]
            gout[:, row, col] = np.conj(ti[:, row, col])
    return gin, gout, delays


def factor_figure(name: str, omega, factor, title: str, path: Path) -> None:
    k = factor.shape[1]
    fig, axes = plt.subplots(4, 1, figsize=(10, 16), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.0, 0.85, k * k))

    for i in range(k):
        for j in range(k):
            entry = factor[:, i, j]
            color = colors[i * k + j]
            label = f"({i},{j})"
            axes[0].plot(omega, np.abs(entry), color=color, lw=1.1,
                         label=label)
            axes[1].plot(omega, np.degrees(np.unwrap(np.angle(entry))),
                         color=color, lw=1.1)
            axes[2].plot(omega, entry.real, color=color, lw=1.1)
            axes[3].plot(omega, entry.imag, color=color, lw=1.1)

    axes[0].set_yscale("log")
    axes[0].set_ylabel("|entry|")
    axes[1].set_ylabel("unwrapped phase [deg]")
    axes[2].set_yscale("symlog", linthresh=1e-8)
    axes[2].set_ylabel("Re")
    axes[3].set_yscale("symlog", linthresh=1e-8)
    axes[3].set_ylabel("Im")
    for ax in axes:
        ax.set_xscale("log")
        ax.grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=8, ncol=3)
    axes[-1].set_xlabel("omega [rad/s]")
    fig.suptitle(f"{name}: {title}", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"wrote {path}")


def main() -> None:
    args = parse_args()
    line_dir = OUTPUT_DIR / args.line
    line_dir.mkdir(parents=True, exist_ok=True)

    factors_csv = run_sweep(args.line, line_dir, args.sweep_app)
    omega, k, tv, ti, h, tau = read_factors(factors_csv)
    gin, gout, delays = build_targets(omega, k, tv, ti, h, tau)
    print(f"{args.line}: {len(omega)} samples, {k} modes, "
          f"delays [us] {1e6 * delays}")

    factor_figure(args.line, omega, gin,
                  "Gin = diag(Hmps) Tv^T (fitter input)",
                  line_dir / "gin.png")
    factor_figure(args.line, omega, gout,
                  "Gout = conj(Ti) (fitter input)",
                  line_dir / "gout.png")


if __name__ == "__main__":
    main()
