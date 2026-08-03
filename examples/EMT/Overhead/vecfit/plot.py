#!/usr/bin/env python3
"""Plot the sampled responses against their fitted rational models.

Reads the sampled-response CSVs and the rational-model JSON files
produced by fit.py and writes one source-versus-fit magnitude figure per
target.
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
    model = json.loads(path.read_text())
    channels = model["rows"] * model["cols"]
    model["_poles"] = np.array(
        [complex(re, im) for re, im in model["poles"]], dtype=np.complex128
    )
    model["_residues"] = np.array(
        [[complex(re, im) for re, im in per_channel] for per_channel in model["residues"]],
        dtype=np.complex128,
    ).reshape(len(model["_poles"]), channels)
    model["_d"] = np.array(model.get("D", np.zeros(channels)))
    model["_e"] = np.array(model.get("E", np.zeros(channels)))
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


def write_fit_plot(dataset: Dataset) -> None:
    omega, source = read_samples_csv(dataset.csv_file)
    model = load_model(dataset.model_file)
    fitted = evaluate_model(model, omega)
    frequency = omega / (2.0 * math.pi)
    labels = channel_labels(model)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), sharex=True, sharey=True)
    plot_panel(axes[0], dataset, frequency, labels, source, dataset.source_title)
    plot_panel(axes[1], dataset, frequency, labels, fitted, dataset.fit_title)
    axes[0].set_ylabel(dataset.ylabel)

    fig.suptitle(f"{dataset.name} fit")
    fig.tight_layout()
    fig.savefig(dataset.plot_file, dpi=180)
    plt.close(fig)

    print(f"wrote {dataset.plot_file}")


def main() -> None:
    for dataset in DATASETS:
        write_fit_plot(dataset)


if __name__ == "__main__":
    main()
