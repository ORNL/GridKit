#!/usr/bin/env python3
"""Plot Overhead R, L, G, and C entries versus frequency."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT_DIR = ROOT / "output"
DEFAULT_CSV = OUTPUT_DIR / "overhead.response.csv"
DEFAULT_PNG = OUTPUT_DIR / "Overhead_Romega.png"
PARAMETERS = (
    ("R", "R(w) [ohm/m]"),
    ("L", "L(w) [H/m]"),
    ("G", "G(w) [S/m]"),
    ("C", "C(w) [F/m]"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--output", type=Path, default=DEFAULT_PNG)
    return parser.parse_args()


def coordinate_column(fieldnames: list[str], csv_path: Path) -> str:
    if "omega" in fieldnames:
        return "omega"
    if "t" in fieldnames:
        return "t"
    raise ValueError(f"CSV must contain an omega column: {csv_path}")


def read_columns(csv_path: Path) -> tuple[list[float], dict[str, dict[str, list[float]]]]:
    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"CSV must contain a header: {csv_path}")
        omega_column = coordinate_column(reader.fieldnames, csv_path)

        columns: dict[str, list[str]] = {}
        for parameter, _ in PARAMETERS:
            prefix = f"Overhead_{parameter}_"
            columns[parameter] = [
                name for name in reader.fieldnames if name.startswith(prefix)
            ]
            if not columns[parameter]:
                raise ValueError(f"No {parameter} columns found in {csv_path}")

        frequency: list[float] = []
        series = {
            parameter: {name: [] for name in parameter_columns}
            for parameter, parameter_columns in columns.items()
        }
        for row in reader:
            omega = float(row[omega_column])
            if omega <= 0.0:
                raise ValueError("omega values must be positive for log-scale plotting")
            frequency.append(omega / (2.0 * 3.14159265358979323846))
            for parameter, parameter_columns in columns.items():
                for name in parameter_columns:
                    series[parameter][name].append(float(row[name]))

    return frequency, series


def main() -> None:
    args = parse_args()
    frequency, series = read_columns(args.csv)

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0), sharex=True)
    for ax, (parameter, ylabel) in zip(axes.ravel(), PARAMETERS):
        prefix = f"Overhead_{parameter}_"
        for name, values in series[parameter].items():
            ax.plot(
                frequency,
                values,
                linewidth=1.1,
                label=name.removeprefix(prefix),
            )

        ax.set_xscale("log")
        values = [value for column in series[parameter].values() for value in column]
        if all(value > 0.0 for value in values):
            ax.set_yscale("log")
        elif any(value != 0.0 for value in values):
            ax.set_yscale("symlog", linthresh=1.0e-15)
        else:
            ax.set_ylim(-1.0e-15, 1.0e-15)
        ax.set_title(parameter)
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7, ncol=2)

    for ax in axes[-1, :]:
        ax.set_xlabel("frequency [Hz]")

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)


if __name__ == "__main__":
    main()
