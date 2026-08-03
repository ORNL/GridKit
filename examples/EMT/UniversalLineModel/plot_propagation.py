#!/usr/bin/env python3
"""Validate the fitted propagation factors by composition.

The factor split Gout, delays, Gin is gauge dependent, but the composed
phase-domain propagation matrix

    P(j omega) = Gout(j omega) diag(exp(-j omega tau)) Gin(j omega)

is not: per mode the eigenvector phase cancels between the input and
output factors. This script assembles the model-free reference
P_ref = conj(Ti) diag(H) Tv^T from the sweep monitors and compares it
against the composition of the three fitted pieces, so any real fitting
or convention error shows up as composed error rather than as gauge
noise.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE / "response.csv"
DEFAULT_MODEL = HERE / "propagation.model.json"
DEFAULT_OUTPUT = HERE / "propagation_check.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def read_sweep(csv_path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{csv_path} must contain a header")
        rows = list(reader)

    k = 0
    while f"Overhead_Yc_real_{k}_{k}" in reader.fieldnames:
        k += 1
    if k == 0:
        raise ValueError(f"{csv_path} contains no Overhead_Yc columns")

    omega = np.array([float(r["omega"]) for r in rows])
    m = len(rows)

    def matrix(prefix: str) -> np.ndarray:
        out = np.zeros((m, k, k), dtype=np.complex128)
        for sample, r in enumerate(rows):
            for i in range(k):
                for j in range(k):
                    out[sample, i, j] = complex(
                        float(r[f"{prefix}_real_{i}_{j}"]),
                        float(r[f"{prefix}_imag_{i}_{j}"]),
                    )
        return out

    h = np.array(
        [
            [
                complex(
                    float(r[f"Overhead_H_real_{i}"]),
                    float(r[f"Overhead_H_imag_{i}"]),
                )
                for i in range(k)
            ]
            for r in rows
        ]
    )
    return omega, {"Tv": matrix("Overhead_Tv"), "Ti": matrix("Overhead_Ti"), "H": h}


def load_factor(model: dict) -> dict:
    rows, cols = model["rows"], model["cols"]
    poles = np.array([complex(re, im) for re, im in model["poles"]])
    residues = np.array(
        [
            [complex(re, im) for row in matrix for re, im in row]
            for matrix in model["residues"]
        ],
        dtype=np.complex128,
    ).reshape(len(poles), rows, cols)
    d = np.array(model["D"]) if "D" in model else np.zeros((rows, cols))
    return {"poles": poles, "residues": residues, "d": d}


def evaluate_factor(factor: dict, omega: np.ndarray) -> np.ndarray:
    s = 1j * omega[:, None]
    weights = 1.0 / (s - factor["poles"][None, :])
    values = np.einsum("mq,qrc->mrc", weights, factor["residues"])
    return values + factor["d"][None, :, :]


def main() -> None:
    args = parse_args()

    omega, sweep = read_sweep(args.csv)
    spec = json.loads(args.model.read_text())
    delays = np.array(spec["delays"])
    gin = evaluate_factor(load_factor(spec["input"]), omega)
    gout = evaluate_factor(load_factor(spec["output"]), omega)

    # Model-free reference: gauge-invariant dyad assembly of the raw sweep.
    reference = np.einsum(
        "mik,mk,mjk->mij", np.conj(sweep["Ti"]), sweep["H"], sweep["Tv"]
    )
    # Composition of the fitted factors with the extracted delays.
    phase = np.exp(-1j * omega[:, None] * delays[None, :])
    composed = np.einsum("mik,mk,mkj->mij", gout, phase, gin)

    frequency = omega / (2.0 * math.pi)
    k = reference.shape[1]
    labels = [f"{i}{j}" for i in range(k) for j in range(k)]

    scale = np.sqrt(np.mean(np.abs(reference.reshape(len(omega), -1)) ** 2, axis=1))
    error = (
        np.abs((composed - reference).reshape(len(omega), -1)).max(axis=1) / scale
    )
    rel_rms = float(
        np.sqrt(
            np.mean(np.abs(composed - reference) ** 2)
            / np.mean(np.abs(reference) ** 2)
        )
    )

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8))
    for ax, values, title in (
        (axes[0], reference, "Sweep dyad assembly"),
        (axes[1], composed, "Composed fitted factors"),
    ):
        flat = values.reshape(len(omega), -1)
        for channel, label in enumerate(labels):
            ax.plot(frequency, np.abs(flat[:, channel]), linewidth=1.2, label=label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(title)
        ax.set_xlabel("frequency [Hz]")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(title="P", fontsize=7, ncol=3)
    axes[1].sharey(axes[0])
    axes[0].set_ylabel("|P| [-]")

    axes[2].plot(frequency, error, linewidth=1.2, color="tab:red")
    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].set_title("Pointwise composition error")
    axes[2].set_xlabel("frequency [Hz]")
    axes[2].set_ylabel("max |error| / rms reference")
    axes[2].grid(True, which="both", alpha=0.3)

    fig.suptitle(
        "Phase-domain propagation: composed fitted factors vs sweep, "
        f"relative RMS {rel_rms:.2e}"
    )
    fig.tight_layout()
    fig.savefig(args.output, dpi=180)
    plt.close(fig)

    print(f"wrote {args.output} (composed relative RMS {rel_rms:.2e})")


if __name__ == "__main__":
    main()
