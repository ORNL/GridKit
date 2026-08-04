#!/usr/bin/env python3
"""Compare the true propagation function against its fitted approximation.

Runs the UniversalLineModel application for one gallery line (or all of
them with --all), assembles the model-free reference

    P_ref(j omega) = conj(Ti) diag(H) Tv^T

from the sweep monitors, and composes the fitted approximation

    P_fit(j omega) = Gout_fit diag(exp(-j omega tau)) Gin_fit

from propagation.model.json. The composition is gauge invariant, so the
comparison shows real fitting error only. Two figures land in
output/<line>/: element magnitudes and element phases, each with the
actual response on the left and the approximation on the right.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gallery import ULM_ARGS, ULM_LINE_ARGS, line_model_file
from plot_factors import read_factors
from sweep_vs_fit import eval_model

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output"
BUILD = HERE.parents[2] / "build"
DEFAULT_ULM_APP = BUILD / "application" / "EMT" / "UniversalLineModel" / "UniversalLineModel"


def gallery_lines() -> list[str]:
    return sorted(spec.name.replace(".solver.json", "")
                  for spec in HERE.glob("*.solver.json"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--line", default="345kv-horizontal")
    parser.add_argument("--all", action="store_true",
                        help="run every line with a gallery solver spec")
    parser.add_argument("--ulm-app", type=Path, default=DEFAULT_ULM_APP)
    return parser.parse_args()


def run_ulm(name: str, ulm_app: Path) -> Path:
    line_dir = OUTPUT_DIR / name
    command = [
        str(ulm_app),
        str(line_model_file(name)),
        "-o",
        str(line_dir),
        *ULM_ARGS,
        *ULM_LINE_ARGS.get(name, []),
    ]
    result = subprocess.run(command, capture_output=True, text=True, cwd=HERE)
    for line in result.stdout.splitlines():
        if "VectorFitting:" in line:
            print(line.strip())
    # Exit codes: 0 all targets met and passive, 2 a target missed,
    # 3 targets met but nonpassive. Anything else is a hard failure.
    if result.returncode not in (0, 2, 3):
        raise RuntimeError(
            f"UniversalLineModel failed ({result.returncode}):\n"
            f"{result.stdout}\n{result.stderr}")
    return line_dir


def compose(omega, k, gout, tau, gin):
    """P[m] = Gout[m] diag(exp(-j omega[m] tau)) Gin[m]."""
    out = np.zeros((len(omega), k, k), complex)
    for m in range(len(omega)):
        delay = np.exp(-1j * omega[m] * tau)
        out[m] = gout[m] @ np.diag(delay) @ gin[m]
    return out


def side_by_side(name, omega, reference, fitted, transform, ylabel, title,
                 path, log_y):
    k = reference.shape[1]
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharex=True, sharey=True)
    colors = plt.cm.viridis(np.linspace(0.0, 0.85, k * k))

    for i in range(k):
        for j in range(k):
            color = colors[i * k + j]
            label = f"({i},{j})" if k <= 3 else None
            axes[0].plot(omega, transform(reference[:, i, j]),
                         color=color, lw=1.1, label=label)
            axes[1].plot(omega, transform(fitted[:, i, j]),
                         color=color, lw=1.1)

    for ax, panel in ((axes[0], "actual"), (axes[1], "approximation")):
        ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")
        ax.set_title(f"{title} — {panel}", fontsize=11)
        ax.set_xlabel("omega [rad/s]")
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel(ylabel)
    if k <= 3:
        axes[0].legend(fontsize=8, ncol=3)

    error = np.sqrt(np.mean(np.abs(fitted - reference) ** 2)
                    / max(np.mean(np.abs(reference) ** 2), 1e-30))
    axes[1].annotate(f"composed rel rms {error:.2e}",
                     xy=(0.03, 0.03), xycoords="axes fraction", fontsize=9)

    fig.suptitle(f"{name}: propagation function P = Gout diag(H) Gin",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"wrote {path}")


def run_line(name: str, ulm_app: Path) -> None:
    print(f"=== {name} ===")
    line_dir = run_ulm(name, ulm_app)
    omega, k, tv, ti, h, _ = read_factors(line_dir / "response.csv")

    propagation = json.loads((line_dir / "propagation.model.json").read_text())
    tau = np.array(propagation["delays"]["tau"])
    gin_fit = eval_model(propagation["input"], omega).reshape(len(omega), k, k)
    gout_fit = eval_model(propagation["output"], omega).reshape(len(omega), k, k)

    reference = np.zeros((len(omega), k, k), complex)
    for m in range(len(omega)):
        reference[m] = np.conj(ti[m]) @ np.diag(h[m]) @ tv[m].T

    fitted = compose(omega, k, gout_fit, tau, gin_fit)

    side_by_side(name, omega, reference, fitted,
                 np.abs, "|P entry|", "element magnitudes",
                 line_dir / "propagation_mag.png", log_y=True)
    side_by_side(name, omega, reference, fitted,
                 lambda entry: np.degrees(np.unwrap(np.angle(entry))),
                 "unwrapped phase [deg]", "element phases",
                 line_dir / "propagation_phase.png", log_y=False)


def main() -> None:
    args = parse_args()
    OUTPUT_DIR.mkdir(exist_ok=True)
    for name in gallery_lines() if args.all else [args.line]:
        run_line(name, args.ulm_app)


if __name__ == "__main__":
    main()
