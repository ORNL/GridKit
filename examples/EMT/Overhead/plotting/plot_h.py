#!/usr/bin/env python3
"""Plot modal MPS H magnitude and phase versus frequency."""

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
OUTPUT_DIR = ROOT / "output"
DEFAULT_CSV = OUTPUT_DIR / "overhead.response.csv"
DEFAULT_OUTPUT = OUTPUT_DIR / "Overhead_H.png"


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
    real_prefix = "Overhead_H_real_"
    imag_prefix = "Overhead_H_imag_"
    tau_prefix = "Overhead_Tau_"
    real_modes = {
        int(name.removeprefix(real_prefix))
        for name in fieldnames
        if name.startswith(real_prefix)
    }
    imag_modes = {
        int(name.removeprefix(imag_prefix))
        for name in fieldnames
        if name.startswith(imag_prefix)
    }
    tau_modes = {
        int(name.removeprefix(tau_prefix))
        for name in fieldnames
        if name.startswith(tau_prefix)
    }
    if not real_modes:
        raise ValueError(f"No H columns found in {csv_path}")
    if real_modes != imag_modes:
        raise ValueError(f"Mismatched H real/imag columns in {csv_path}")
    if real_modes != tau_modes:
        raise ValueError(f"Mismatched H/Tau modal columns in {csv_path}")
    return sorted(real_modes)


def unwrap_phase(values: list[float]) -> list[float]:
    if not values:
        return []

    unwrapped = [values[0]]
    offset = 0.0
    previous = values[0]
    for value in values[1:]:
        delta = value - previous
        if delta > math.pi:
            offset -= 2.0 * math.pi
        elif delta < -math.pi:
            offset += 2.0 * math.pi
        unwrapped.append(value + offset)
        previous = value

    return unwrapped


def read_h_columns(
    csv_path: Path,
) -> tuple[list[float], dict[int, dict[str, list[float]]]]:
    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"CSV must contain a header: {csv_path}")
        omega_column = coordinate_column(reader.fieldnames, csv_path)
        modes = modal_indices(reader.fieldnames, csv_path)

        omega_values: list[float] = []
        h_values = {mode: [] for mode in modes}
        tau_values = {mode: [] for mode in modes}
        for row in reader:
            omega = float(row[omega_column])
            if omega <= 0.0:
                raise ValueError("omega values must be positive for log-scale plotting")
            omega_values.append(omega)

            for mode in modes:
                real = float(row[f"Overhead_H_real_{mode}"])
                imag = float(row[f"Overhead_H_imag_{mode}"])
                tau = float(row[f"Overhead_Tau_{mode}"])
                h_values[mode].append((real, imag))
                tau_values[mode].append(tau)

        frequency = [omega / (2.0 * math.pi) for omega in omega_values]
        tau_min = {mode: min(tau_values[mode]) for mode in modes}
        series = {mode: {"magnitude": [], "phase": []} for mode in modes}
        for row_index, omega in enumerate(omega_values):
            for mode in modes:
                real, imag = h_values[mode][row_index]
                angle = omega * tau_min[mode]
                cos_angle = math.cos(angle)
                sin_angle = math.sin(angle)
                mps_real = real * cos_angle - imag * sin_angle
                mps_imag = real * sin_angle + imag * cos_angle
                series[mode]["magnitude"].append(math.hypot(mps_real, mps_imag))
                series[mode]["phase"].append(math.atan2(mps_imag, mps_real))

        for values in series.values():
            values["phase"] = unwrap_phase(values["phase"])

    return frequency, series


def main() -> None:
    args = parse_args()
    frequency, series = read_h_columns(args.csv)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.5), sharex=True)
    mag_ax, phase_ax = axes

    for mode, values in series.items():
        mag_ax.plot(frequency, values["magnitude"], linewidth=1.4, label=f"mode {mode}")
        phase_ax.plot(frequency, values["phase"], linewidth=1.4, label=f"mode {mode}")

    mag_ax.set_xscale("log")
    mag_ax.set_xlabel("frequency [Hz]")
    mag_ax.set_ylabel("|H_mps| [-]")
    mag_ax.set_title("Modal MPS propagation function")
    mag_ax.grid(True, which="both", alpha=0.3)
    mag_ax.legend()

    phase_ax.set_xscale("log")
    phase_ax.set_xlabel("frequency [Hz]")
    phase_ax.set_ylabel("unwrapped angle(H_mps) [rad]")
    phase_ax.grid(True, which="both", alpha=0.3)
    phase_ax.legend()

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
