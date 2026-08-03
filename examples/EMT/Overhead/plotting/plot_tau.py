#!/usr/bin/env python3
"""Plot modal Overhead phase delays versus frequency."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESPONSE = ROOT / "response"
DEFAULT_CSV = RESPONSE / "overhead.response.csv"
DEFAULT_OUTPUT = RESPONSE / "Overhead_Tau.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def coordinate_column(fieldnames: list[str], csv_path: Path) -> str:
    if "omega" in fieldnames:
        return "omega"
    if "t" in fieldnames:
        return "t"
    raise ValueError(f"CSV must contain an omega column: {csv_path}")


def modal_indices(fieldnames: list[str], csv_path: Path) -> list[int]:
    prefix = "Overhead_Tau_"
    modes = sorted(
        int(name.removeprefix(prefix))
        for name in fieldnames
        if name.startswith(prefix)
    )
    if not modes:
        raise ValueError(f"No Tau columns found in {csv_path}")
    return modes


def read_columns(csv_path: Path) -> tuple[list[float], dict[int, list[float]]]:
    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"CSV must contain a header: {csv_path}")

        omega_column = coordinate_column(reader.fieldnames, csv_path)
        modes = modal_indices(reader.fieldnames, csv_path)

        frequency: list[float] = []
        series = {mode: [] for mode in modes}
        for row in reader:
            omega = float(row[omega_column])
            if omega <= 0.0:
                raise ValueError("omega values must be positive for log-scale plotting")
            frequency.append(omega / (2.0 * math.pi))
            for mode in modes:
                series[mode].append(float(row[f"Overhead_Tau_{mode}"]))

    return frequency, series


def main() -> None:
    args = parse_args()
    frequency, series = read_columns(args.csv)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for mode, values in series.items():
        ax.plot(frequency, values, linewidth=1.4, label=f"mode {mode}")

    ax.set_xscale("log")
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel("phase delay [s]")
    ax.set_title("Modal phase delay")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
