#!/usr/bin/env python3
"""Plot modal Gamma attenuation, phase, and angle versus frequency."""

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
DEFAULT_OUTPUT = RESPONSE / "Overhead_Gamma_modes.png"


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
    alpha_prefix = "Overhead_Alpha_"
    beta_prefix = "Overhead_Beta_"
    alpha_modes = {
        int(name.removeprefix(alpha_prefix))
        for name in fieldnames
        if name.startswith(alpha_prefix)
    }
    beta_modes = {
        int(name.removeprefix(beta_prefix))
        for name in fieldnames
        if name.startswith(beta_prefix)
    }
    if not alpha_modes:
        raise ValueError(f"No Alpha columns found in {csv_path}")
    if alpha_modes != beta_modes:
        raise ValueError(f"Mismatched Alpha/Beta columns in {csv_path}")
    return sorted(alpha_modes)


def read_modal_columns(
    csv_path: Path,
) -> tuple[list[float], dict[int, dict[str, list[float]]]]:
    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"CSV must contain a header: {csv_path}")
        omega_column = coordinate_column(reader.fieldnames, csv_path)
        modes = modal_indices(reader.fieldnames, csv_path)

        frequency: list[float] = []
        series = {mode: {"alpha": [], "beta": [], "angle": []} for mode in modes}
        for row in reader:
            omega = float(row[omega_column])
            if omega <= 0.0:
                raise ValueError("omega values must be positive for log-scale plotting")
            frequency.append(omega / (2.0 * math.pi))

            for mode in modes:
                alpha = float(row[f"Overhead_Alpha_{mode}"])
                beta = float(row[f"Overhead_Beta_{mode}"])
                series[mode]["alpha"].append(alpha)
                series[mode]["beta"].append(beta)
                series[mode]["angle"].append(math.atan2(beta, alpha))

    return frequency, series


def require_positive(values: list[float], label: str) -> None:
    if any(value <= 0.0 for value in values):
        raise ValueError(f"{label} values must be positive for log-scale plotting")


def positive_limits(values: list[float], label: str) -> tuple[float, float]:
    require_positive(values, label)
    lower = min(values)
    upper = max(values)
    return 0.8 * lower, 1.25 * upper


def save_gamma_plot(
    frequency: list[float],
    series: dict[int, dict[str, list[float]]],
    output: Path,
) -> None:
    modes = list(series)
    if len(modes) != 3:
        raise ValueError("Expected exactly three modes for the 2x2 Gamma plot")

    require_positive(frequency, "frequency")
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.5), sharex=True)
    mode_axes = axes.ravel()[:3]
    angle_ax = axes.ravel()[3]

    for ax, mode in zip(mode_axes, modes):
        values = series[mode]
        ax.plot(frequency, values["alpha"], color="red", linewidth=1.4, label="alpha")
        ax.plot(frequency, values["beta"], color="blue", linewidth=1.4, label="beta")
        ax.set_title(f"Mode {mode}")
        ax.set_ylabel("propagation constant [1/m]")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(*positive_limits(values["alpha"] + values["beta"], f"mode {mode}"))
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()

    for mode, values in series.items():
        angle_ax.plot(frequency, values["angle"], linewidth=1.4, label=f"mode {mode}")

    angle_ax.set_title("Modal angle(gamma)")
    angle_ax.set_ylabel("angle(gamma) [rad]")
    angle_ax.set_xscale("log")
    angle_ax.grid(True, which="both", alpha=0.3)
    angle_ax.legend()

    for ax in axes[-1, :]:
        ax.set_xlabel("frequency [Hz]")

    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    frequency, series = read_modal_columns(args.csv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    save_gamma_plot(frequency, series, args.output)


if __name__ == "__main__":
    main()
