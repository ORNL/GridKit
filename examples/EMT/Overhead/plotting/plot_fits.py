#!/usr/bin/env python3
"""Plot the sampled responses against their fitted rational models.

Reads the sampled-response CSVs and the rational-model JSON files
produced by fit.py and writes one source-versus-fit figure per target:
magnitude overlays side by side, with the pole count and the relative
RMS error reported in the title and the pointwise relative error in a
third panel.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE.parent / "output"


@dataclass(frozen=True)
class Dataset:
    name: str
    csv_file: Path
    model_file: Path
    plot_file: Path
    legend_title: str
    ylabel: str
    source_title: str
    fit_title: str


DATASETS = [
    Dataset(
        name="Yc",
        csv_file=OUTPUT_DIR / "yc.csv",
        model_file=OUTPUT_DIR / "yc.model.json",
        plot_file=OUTPUT_DIR / "yc_fit.png",
        legend_title="Yc",
        ylabel="|Yc| [S]",
        source_title="Yc samples",
        fit_title="Rational approximation",
    ),
    Dataset(
        name="Hmin",
        csv_file=OUTPUT_DIR / "hmin.csv",
        model_file=OUTPUT_DIR / "hmin.model.json",
        plot_file=OUTPUT_DIR / "hmin_fit.png",
        legend_title="Hmin",
        ylabel="|Hmin| [-]",
        source_title="Hmin samples",
        fit_title="Rational approximation",
    ),
]


def read_samples_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a sampled-response CSV: omega, then one re/im pair per channel."""
    with path.open(newline="") as stream:
        reader = csv.reader(stream)
        next(reader)
        table = np.array([[float(value) for value in row] for row in reader])

    omega = table[:, 0]
    values = table[:, 1::2] + 1j * table[:, 2::2]
    return omega, values


def load_model(path: Path) -> dict:
    """Load a rational model JSON with matrix-shaped D, E, and residues."""
    model = json.loads(path.read_text())
    rows, cols = model["rows"], model["cols"]
    channels = rows * cols

    model["_poles"] = np.array(
        [complex(re, im) for re, im in model["poles"]], dtype=np.complex128
    )
    model["_residues"] = np.array(
        [
            [complex(re, im) for row in matrix for re, im in row]
            for matrix in model["residues"]
        ],
        dtype=np.complex128,
    ).reshape(len(model["_poles"]), channels)
    for key, slot in (("D", "_d"), ("E", "_e")):
        flat = np.zeros(channels)
        if key in model:
            flat = np.array(model[key], dtype=np.float64).reshape(channels)
        model[slot] = flat
    return model


def evaluate_model(model: dict, omega: np.ndarray) -> np.ndarray:
    s = 1j * omega[:, None]
    weights = 1.0 / (s - model["_poles"][None, :])
    values = weights @ model["_residues"]
    values += model["_d"][None, :] + s * model["_e"][None, :]
    return values


def channel_labels(model: dict) -> list[str]:
    if model["cols"] == 1:
        return [f"{row}" for row in range(model["rows"])]
    return [
        f"{row}{col}"
        for row in range(model["rows"])
        for col in range(model["cols"])
    ]


def relative_rms(source: np.ndarray, fitted: np.ndarray) -> float:
    return float(
        np.sqrt(np.mean(np.abs(fitted - source) ** 2) / np.mean(np.abs(source) ** 2))
    )


def plot_panel(ax, dataset, frequency, labels, values, title) -> None:
    for channel, label in enumerate(labels):
        ax.plot(
            frequency,
            np.abs(values[:, channel]),
            linewidth=1.2,
            label=label,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.set_xlabel("frequency [Hz]")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(title=dataset.legend_title, fontsize=7, ncol=3)


def plot_error_panel(ax, dataset, frequency, source, fitted) -> None:
    scale = np.sqrt(np.mean(np.abs(source) ** 2, axis=1))
    error = np.abs(fitted - source).max(axis=1) / scale
    ax.plot(frequency, error, linewidth=1.2, color="tab:red")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Pointwise error")
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel("max |error| / rms sample")
    ax.grid(True, which="both", alpha=0.3)


def write_fit_plot(dataset: Dataset) -> None:
    omega, source = read_samples_csv(dataset.csv_file)
    model = load_model(dataset.model_file)
    fitted = evaluate_model(model, omega)
    frequency = omega / (2.0 * math.pi)
    labels = channel_labels(model)
    rel_rms = relative_rms(source, fitted)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8))
    plot_panel(axes[0], dataset, frequency, labels, source, dataset.source_title)
    plot_panel(axes[1], dataset, frequency, labels, fitted, dataset.fit_title)
    axes[1].sharey(axes[0])
    plot_error_panel(axes[2], dataset, frequency, source, fitted)
    axes[0].set_ylabel(dataset.ylabel)

    fig.suptitle(
        f"{dataset.name} fit: {len(model['_poles'])} poles, "
        f"relative RMS {rel_rms:.2e}"
    )
    fig.tight_layout()
    fig.savefig(dataset.plot_file, dpi=180)
    plt.close(fig)

    print(f"wrote {dataset.plot_file} ({len(model['_poles'])} poles, rel rms {rel_rms:.2e})")


def main() -> None:
    for dataset in DATASETS:
        write_fit_plot(dataset)


if __name__ == "__main__":
    main()
