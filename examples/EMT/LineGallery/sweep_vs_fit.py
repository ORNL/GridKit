#!/usr/bin/env python3
"""Fit the monitored sweep quantities directly and plot data beside fit.

For one catalog line the script sweeps Alpha, Beta, R, L, G, C with the
FrequencyResponse application, hands gamma = alpha + j beta, Z = R + j w L,
and Y = G + j w C straight to the standalone VectorFitting application, and
renders the monitored curves (left column) next to the fitted rational
model (right column). Everything generated lands in output/<line>/. This
isolates the estimator and the sweep from the propagation factor
construction: if these fits are clean, any remaining fitting trouble lives
in the quantities the UniversalLineModel builds, not here.
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
DEFAULT_FIT_APP = BUILD / "application" / "Fitting" / "VectorFitting"

FREQUENCY = {"start": 10.0, "stop": 1.0e8, "points": 401, "scale": "log"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--line", default="345kv-horizontal")
    parser.add_argument("--sweep-app", type=Path, default=DEFAULT_SWEEP_APP)
    parser.add_argument("--fit-app", type=Path, default=DEFAULT_FIT_APP)
    return parser.parse_args()


def run_sweep(name: str, line_dir: Path, sweep_app: Path) -> Path:
    # Paths inside the spec resolve relative to the spec file, which
    # lives in output/<line>/ alongside every other generated artifact.
    spec = {
        "model": str(line_model_file(name)),
        "frequency": FREQUENCY,
        "output": {
            "file": "rlgc.csv",
            "variables": ["Alpha", "Beta", "R", "L", "G", "C"],
        },
    }
    spec_path = line_dir / "rlgc.solver.json"
    spec_path.write_text(json.dumps(spec, indent=4) + "\n")
    subprocess.run([str(sweep_app), str(spec_path)],
                   check=True, cwd=HERE, stdout=subprocess.DEVNULL)
    return line_dir / "rlgc.csv"


def read_sweep(path: Path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        names = reader.fieldnames or []
    omega = np.array([float(r["omega"]) for r in rows])
    m = len(rows)

    n = 0
    while f"Overhead_R_{n}_{n}" in names:
        n += 1
    k = 0
    while f"Overhead_Alpha_{k}" in names:
        k += 1
    if n == 0 or k == 0:
        raise ValueError(f"{path} lacks R or Alpha columns")

    def matrix(prefix: str) -> np.ndarray:
        out = np.zeros((m, n, n))
        for s, r in enumerate(rows):
            for i in range(n):
                for j in range(n):
                    out[s, i, j] = float(r[f"{prefix}_{i}_{j}"])
        return out

    alpha = np.array([[float(r[f"Overhead_Alpha_{q}"]) for q in range(k)] for r in rows])
    beta = np.array([[float(r[f"Overhead_Beta_{q}"]) for q in range(k)] for r in rows])
    return omega, n, k, matrix("Overhead_R"), matrix("Overhead_L"), \
        matrix("Overhead_G"), matrix("Overhead_C"), alpha, beta


def write_samples(path: Path, omega: np.ndarray, values: np.ndarray) -> None:
    """values: (samples, channels) complex, channels row-major."""
    channels = values.shape[1]
    with path.open("w") as stream:
        stream.write("omega," + ",".join(
            f"re_{c},im_{c}" for c in range(channels)) + "\n")
        for s in range(len(omega)):
            row = [repr(omega[s])]
            for c in range(channels):
                row.append(repr(values[s, c].real))
                row.append(repr(values[s, c].imag))
            stream.write(",".join(row) + "\n")


def run_fit(fit_app: Path, line_dir: Path, name: str,
            rows: int, cols: int) -> dict:
    spec = {
        "samples": f"{name}.samples.csv",
        "rows": rows,
        "cols": cols,
        "fit": {
            "min_poles": 2,
            "max_poles": 30,
            "target_rel_rms": 1.0e-3,
            "terms": "linear",
            "weighting": "inverse_magnitude",
        },
        "output": {"model": f"{name}.model.json"},
    }
    spec_path = line_dir / f"{name}.fit.json"
    spec_path.write_text(json.dumps(spec, indent=4) + "\n")
    result = subprocess.run([str(fit_app), str(spec_path)],
                            capture_output=True, text=True, timeout=600)
    stats = [l for l in result.stdout.splitlines() if "VectorFitting:" in l]
    print(f"{name}: " + (stats[-1].strip() if stats else f"exit {result.returncode}"))
    if result.returncode not in (0, 1, 2):
        raise RuntimeError(f"{name} fit failed:\n{result.stdout}\n{result.stderr}")
    return json.loads((line_dir / f"{name}.model.json").read_text())


def eval_model(model: dict, omega: np.ndarray) -> np.ndarray:
    poles = np.array([complex(p[0], p[1]) for p in model["poles"]])
    channels = model["rows"] * model["cols"]

    def flat_real(entry) -> np.ndarray:
        out: list[float] = []

        def walk(node):
            if isinstance(node, (int, float)):
                out.append(float(node))
            else:
                for child in node:
                    walk(child)

        walk(entry)
        return np.array(out).reshape(channels)

    def flat_complex(entry) -> np.ndarray:
        out: list[complex] = []

        def walk(node):
            if (isinstance(node, list) and len(node) == 2
                    and all(isinstance(x, (int, float)) for x in node)):
                out.append(complex(node[0], node[1]))
            else:
                for child in node:
                    walk(child)

        walk(entry)
        return np.array(out).reshape(channels)

    d = flat_real(model["D"]) if "D" in model else np.zeros(channels)
    e = flat_real(model["E"]) if "E" in model else np.zeros(channels)
    residues = np.array([flat_complex(r) for r in model["residues"]])

    s = 1j * omega
    out = np.zeros((len(omega), channels), complex)
    for m in range(len(omega)):
        out[m] = d + s[m] * e + residues.T @ (1.0 / (s[m] - poles))
    return out


def rel_rms(data: np.ndarray, fit: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.abs(fit - data) ** 2))
                 / max(np.sqrt(np.mean(np.abs(data) ** 2)), 1e-30))


def panel_pair(axes_row, omega, data_curves, fit_curves, title, ylabel,
               log_y=True, labels=None) -> None:
    """Left: monitored. Right: fit (dashed) over data (light gray)."""
    left, right = axes_row
    n = len(data_curves)
    colors = plt.cm.viridis(np.linspace(0.0, 0.85, n))
    for i, curve in enumerate(data_curves):
        left.plot(omega, curve, color=colors[i], lw=1.2,
                  label=labels[i] if labels else None)
    for i in range(n):
        right.plot(omega, data_curves[i], color="0.8", lw=2.2)
        right.plot(omega, fit_curves[i], color=colors[i], lw=1.0, ls="--")
    for ax, name in ((left, f"{title} — monitored"),
                     (right, f"{title} — VF model")):
        ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")
        ax.set_title(name, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.grid(True, which="both", alpha=0.25)
    if labels:
        left.legend(fontsize=7, ncol=2)


def main() -> None:
    args = parse_args()
    line_dir = OUTPUT_DIR / args.line
    line_dir.mkdir(parents=True, exist_ok=True)

    sweep_csv = run_sweep(args.line, line_dir, args.sweep_app)
    omega, n, k, R, L, G, C, alpha, beta = read_sweep(sweep_csv)
    print(f"{args.line}: {len(omega)} samples, {n}x{n} conductor matrices, {k} modes")

    gamma = alpha + 1j * beta
    Z = (R + 1j * omega[:, None, None] * L).reshape(len(omega), n * n)
    Y = (G + 1j * omega[:, None, None] * C).reshape(len(omega), n * n)

    write_samples(line_dir / "gamma.samples.csv", omega, gamma)
    write_samples(line_dir / "z.samples.csv", omega, Z)
    write_samples(line_dir / "y.samples.csv", omega, Y)

    gamma_fit = eval_model(
        run_fit(args.fit_app, line_dir, "gamma", k, 1), omega)
    z_fit = eval_model(run_fit(args.fit_app, line_dir, "z", n, n), omega)
    y_fit = eval_model(run_fit(args.fit_app, line_dir, "y", n, n), omega)

    R_fit = z_fit.reshape(len(omega), n, n).real
    L_fit = z_fit.reshape(len(omega), n, n).imag / omega[:, None, None]
    G_fit = y_fit.reshape(len(omega), n, n).real
    C_fit = y_fit.reshape(len(omega), n, n).imag / omega[:, None, None]

    err = {
        "gamma": rel_rms(gamma, gamma_fit),
        "Z": rel_rms(Z, z_fit),
        "Y": rel_rms(Y, y_fit),
    }
    print("rel rms:", {q: f"{v:.3e}" for q, v in err.items()})

    fig, axes = plt.subplots(6, 2, figsize=(11, 21), sharex=True)
    mode_labels = [f"mode {q}" for q in range(k)]

    panel_pair(axes[0], omega, [alpha[:, q] for q in range(k)],
               [gamma_fit[:, q].real for q in range(k)],
               "alpha (Re gamma)", "Np/m", log_y=True, labels=mode_labels)
    panel_pair(axes[1], omega, [beta[:, q] for q in range(k)],
               [gamma_fit[:, q].imag for q in range(k)],
               "beta (Im gamma)", "rad/m", log_y=True, labels=mode_labels)

    def entries(mat: np.ndarray):
        curves = [mat[:, i, i] for i in range(n)] + [mat[:, 0, 1]]
        labels = [f"({i},{i})" if i < 3 else None for i in range(n)] + ["(0,1)"]
        return curves, labels

    for row, (name, unit, data_m, fit_m, logy) in enumerate([
        ("R", "ohm/m", R, R_fit, True),
        ("L", "H/m", L, L_fit, False),
        ("G", "S/m", G, G_fit, False),
        ("C", "F/m", C, C_fit, False),
    ], start=2):
        curves, labels = entries(data_m)
        fits, _ = entries(fit_m)
        panel_pair(axes[row], omega, curves, fits, name, unit,
                   log_y=logy, labels=labels)

    annotations = [
        (axes[0, 1], f"gamma rel rms {err['gamma']:.2e}"),
        (axes[2, 1], f"Z rel rms {err['Z']:.2e}"),
        (axes[4, 1], f"Y rel rms {err['Y']:.2e}"),
    ]
    for ax, text in annotations:
        ax.annotate(text, xy=(0.03, 0.05), xycoords="axes fraction", fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel("omega [rad/s]")
    fig.suptitle(f"{args.line}: monitored sweep vs VectorFitting rational model",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    plot_file = line_dir / "sweep_vs_fit.png"
    fig.savefig(plot_file, dpi=130)
    print(f"wrote {plot_file}")


if __name__ == "__main__":
    main()
