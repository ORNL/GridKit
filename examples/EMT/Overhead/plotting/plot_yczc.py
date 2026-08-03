#!/usr/bin/env python3
"""Plot Overhead Yc and Zc magnitude and phase versus frequency."""

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
DEFAULT_PNG = RESPONSE / "Overhead_YcZc.png"
PARAMETERS = (
    ("Yc", "|Yc(f)| [S]", "angle(Yc(f)) [rad]"),
    ("Zc", "|Zc(f)| [ohm]", "angle(Zc(f)) [rad]"),
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


def read_complex_columns(
    csv_path: Path,
) -> tuple[list[float], dict[str, dict[str, dict[str, list[float]]]]]:
    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"CSV must contain a header: {csv_path}")
        omega_column = coordinate_column(reader.fieldnames, csv_path)

        columns: dict[str, list[str]] = {}
        for parameter, _, _ in PARAMETERS:
            real_prefix = f"Overhead_{parameter}_real_"
            imag_prefix = f"Overhead_{parameter}_imag_"
            real_suffixes = {
                name.removeprefix(real_prefix)
                for name in reader.fieldnames
                if name.startswith(real_prefix)
            }
            imag_suffixes = {
                name.removeprefix(imag_prefix)
                for name in reader.fieldnames
                if name.startswith(imag_prefix)
            }
            suffixes = sorted(real_suffixes & imag_suffixes)
            if not suffixes:
                raise ValueError(f"No complex {parameter} columns found in {csv_path}")
            if real_suffixes != imag_suffixes:
                raise ValueError(f"Mismatched real/imag {parameter} columns in {csv_path}")
            columns[parameter] = suffixes

        frequency: list[float] = []
        series = {
            parameter: {
                suffix: {"magnitude": [], "phase": []} for suffix in suffixes
            }
            for parameter, suffixes in columns.items()
        }

        for row in reader:
            omega = float(row[omega_column])
            if omega <= 0.0:
                raise ValueError("omega values must be positive for log-scale plotting")
            frequency.append(omega / (2.0 * math.pi))

            for parameter, suffixes in columns.items():
                for suffix in suffixes:
                    real = float(row[f"Overhead_{parameter}_real_{suffix}"])
                    imag = float(row[f"Overhead_{parameter}_imag_{suffix}"])
                    series[parameter][suffix]["magnitude"].append(math.hypot(real, imag))
                    series[parameter][suffix]["phase"].append(math.atan2(imag, real))

    return frequency, series


def main() -> None:
    args = parse_args()
    frequency, series = read_complex_columns(args.csv)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), sharex=True)
    for row, (parameter, mag_ylabel, phase_ylabel) in enumerate(PARAMETERS):
        mag_ax = axes[row, 0]
        phs_ax = axes[row, 1]
        for suffix, values in series[parameter].items():
            mag_ax.plot(frequency, values["magnitude"], linewidth=1.1, label=suffix)
            phs_ax.plot(frequency, values["phase"], linewidth=1.1, label=suffix)

        mag_ax.set_xscale("log")
        phs_ax.set_xscale("log")
        mag_ax.set_title(f"|{parameter}|")
        phs_ax.set_title(f"angle({parameter})")
        mag_ax.set_ylabel(mag_ylabel)
        phs_ax.set_ylabel(phase_ylabel)
        mag_ax.grid(True, which="both", alpha=0.3)
        phs_ax.grid(True, which="both", alpha=0.3)
        mag_ax.legend(fontsize=7, ncol=2)
        phs_ax.legend(fontsize=7, ncol=2)

    for ax in axes[-1, :]:
        ax.set_xlabel("frequency [Hz]")

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)


if __name__ == "__main__":
    main()
