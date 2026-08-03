from __future__ import annotations

import csv
import json
import math
import re
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
    column_pattern: re.Pattern[str]
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
        column_pattern=re.compile(r"^re_Yc_(\d+)_(\d+)$"),
        legend_title="Yc",
        ylabel="|Yc| [S]",
        source_title="Yc CSV",
        fit_title="Rational approximation",
    ),
    Dataset(
        name="Hmin",
        csv_file=OUTPUT_DIR / "hmin.csv",
        model_file=OUTPUT_DIR / "hmin.model.json",
        plot_file=OUTPUT_DIR / "hmin_fit.png",
        column_pattern=re.compile(r"^re_Hmin_(\d+)_(\d+)$"),
        legend_title="Hmin",
        ylabel="|Hmin| [-]",
        source_title="Hmin CSV",
        fit_title="Rational approximation",
    ),
]


def complex_pair(value: list[float]) -> complex:
    return complex(value[0], value[1])


def complex_matrix(rows: list[list[list[float]]]) -> np.ndarray:
    return np.array(
        [[complex_pair(value) for value in row] for row in rows],
        dtype=np.complex128,
    )


def channel_label(entry: tuple[int, int], dataset: Dataset) -> str:
    return f"{entry[0]}{entry[1]}"


def matrix_entries(fieldnames: list[str], dataset: Dataset) -> list[tuple[int, int]]:
    entries = []
    for name in fieldnames:
        match = dataset.column_pattern.match(name)
        if match is None:
            continue
        row, col = match.groups()
        entry = (int(row), int(col))
        imag_name = dataset.column_pattern.sub(
            f"im_{dataset.name}_{entry[0]}_{entry[1]}", name
        )
        if imag_name not in fieldnames:
            raise ValueError(f"missing imaginary column for {dataset.name}_{entry[0]}_{entry[1]}")
        entries.append(entry)

    if not entries:
        raise ValueError(f"{dataset.csv_file} contains no {dataset.name} matrix columns")

    return sorted(entries)


def read_vecfit_csv(dataset: Dataset) -> tuple[list[float], list[tuple[int, int]], list[list[complex]]]:
    with dataset.csv_file.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{dataset.csv_file} must contain a header")
        entries = matrix_entries(reader.fieldnames, dataset)

        frequency = []
        values = []
        for row in reader:
            frequency.append(float(row["freq_Hz"]))
            values.append(
                [
                    complex(
                        float(row[f"re_{dataset.name}_{entry[0]}_{entry[1]}"]),
                        float(row[f"im_{dataset.name}_{entry[0]}_{entry[1]}"]),
                    )
                    for entry in entries
                ]
            )

    return frequency, entries, values


def evaluate_model(dataset: Dataset, frequency: list[float]) -> list[list[complex]]:
    model = json.loads(dataset.model_file.read_text())
    rows, cols = model["shape"]
    poles = np.array([complex_pair(value) for value in model["poles"]], dtype=np.complex128)
    direct = np.array([complex_pair(value) for value in model["d"]], dtype=np.complex128).reshape(rows, cols)
    proportional = np.array([complex_pair(value) for value in model["e"]], dtype=np.complex128).reshape(rows, cols)
    residues = None
    c_matrix = None
    b_matrix = None

    if "residues" in model:
        residues = np.array(
            [[complex_pair(value) for value in row] for row in model["residues"]],
            dtype=np.complex128,
        ).reshape(len(poles), rows, cols)
    elif "C" in model and "B" in model:
        c_matrix = complex_matrix(model["C"])
        b_matrix = complex_matrix(model["B"])
        expected_c_shape = (rows, len(poles))
        expected_b_shape = (len(poles), cols)
        if c_matrix.shape != expected_c_shape:
            raise ValueError(
                f"{dataset.model_file} C shape {c_matrix.shape} "
                f"does not match {expected_c_shape}"
            )
        if b_matrix.shape != expected_b_shape:
            raise ValueError(
                f"{dataset.model_file} B shape {b_matrix.shape} "
                f"does not match {expected_b_shape}"
            )
    else:
        raise ValueError(f"{dataset.model_file} must contain either residues or C/B state-space factors")

    fitted = []
    for freq_hz in frequency:
        s = complex(0.0, 2.0 * math.pi * freq_hz)
        value = direct + s * proportional
        if residues is not None:
            for pole, residue in zip(poles, residues):
                value += residue / (s - pole)
        else:
            weights = 1.0 / (s - poles)
            value += (c_matrix * weights[np.newaxis, :]) @ b_matrix
        fitted.append(list(value.reshape(-1)))

    return fitted


def plot_panel(
    ax,
    dataset: Dataset,
    frequency: list[float],
    entries: list[tuple[int, int]],
    values: list[list[complex]],
    title: str,
) -> None:
    for channel, entry in enumerate(entries):
        magnitude = [abs(row[channel]) for row in values]
        ax.plot(frequency, magnitude, linewidth=1.2, label=channel_label(entry, dataset))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.set_xlabel("frequency [Hz]")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(title=dataset.legend_title, fontsize=7, ncol=3)


def write_fit_plot(dataset: Dataset) -> None:
    frequency, entries, source = read_vecfit_csv(dataset)
    fitted = evaluate_model(dataset, frequency)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), sharex=True, sharey=True)
    plot_panel(axes[0], dataset, frequency, entries, source, dataset.source_title)
    plot_panel(axes[1], dataset, frequency, entries, fitted, dataset.fit_title)
    axes[0].set_ylabel(dataset.ylabel)

    fig.suptitle(f"{dataset.name} fit")
    fig.tight_layout()
    fig.savefig(dataset.plot_file, dpi=180)
    plt.close(fig)

    print(f"wrote {dataset.plot_file}")


def write_fit_plots() -> None:
    for dataset in DATASETS:
        write_fit_plot(dataset)


write_fit_plots()
