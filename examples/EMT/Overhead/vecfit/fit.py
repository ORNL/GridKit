#!/usr/bin/env python3
"""Prepare vector-fitting inputs from the Overhead sweep and fit them.

Reads the monitor CSV written by the FrequencyResponse application,
extracts the characteristic admittance and the minimum-phase-shifted
propagation function, writes the sampled-response CSVs and solver
specifications consumed by the VectorFitting application, and runs the
application on both targets.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXAMPLE_DIR = HERE.parent
INPUT_CSV = EXAMPLE_DIR / "response" / "overhead.response.csv"
OUTPUT_DIR = EXAMPLE_DIR / "output"
DEFAULT_APP = HERE.parents[3] / "build" / "application" / "Fitting" / "VectorFitting"

YC_POLES = 10
HMIN_POLES = 20


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--app",
        type=Path,
        default=DEFAULT_APP,
        help="Path to the VectorFitting application binary",
    )
    return parser.parse_args()


def read_monitor_csv() -> tuple[list[str], list[dict[str, str]]]:
    with INPUT_CSV.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{INPUT_CSV} must contain a header")
        if "omega" not in reader.fieldnames:
            raise ValueError(f"{INPUT_CSV} must contain an omega column")
        return reader.fieldnames, list(reader)


def conductor_count(fieldnames: list[str]) -> int:
    count = 0
    while f"Overhead_Yc_real_{count}_{count}" in fieldnames:
        count += 1
    if count == 0:
        raise ValueError(f"{INPUT_CSV} contains no Overhead_Yc monitor columns")
    return count


def write_samples_csv(
    path: Path, rows: list[dict[str, str]], channels: list[tuple[str, str]]
) -> None:
    """Write the VectorFitting sampled-response CSV: omega, re, im per channel."""
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        header = ["omega"]
        for real_name, _ in channels:
            label = real_name.removeprefix("Overhead_").replace("real_", "")
            header.extend([f"re_{label}", f"im_{label}"])
        writer.writerow(header)

        for row in rows:
            record = [row["omega"]]
            for real_name, imag_name in channels:
                record.extend([row[real_name], row[imag_name]])
            writer.writerow(record)


def write_hmin_csv(
    path: Path, rows: list[dict[str, str]], modes: int
) -> list[float]:
    """Shift each propagation mode to minimum phase and write its samples."""
    tau_min = [
        min(float(row[f"Overhead_Tau_{mode}"]) for row in rows)
        for mode in range(modes)
    ]

    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        header = ["omega"]
        for mode in range(modes):
            header.extend([f"re_Hmin_{mode}", f"im_Hmin_{mode}"])
        writer.writerow(header)

        for row in rows:
            omega = float(row["omega"])
            record = [row["omega"]]
            for mode in range(modes):
                real = float(row[f"Overhead_H_real_{mode}"])
                imag = float(row[f"Overhead_H_imag_{mode}"])
                angle = omega * tau_min[mode]
                record.extend(
                    [
                        f"{real * math.cos(angle) - imag * math.sin(angle):.17e}",
                        f"{real * math.sin(angle) + imag * math.cos(angle):.17e}",
                    ]
                )
            writer.writerow(record)

    return tau_min


def fit_model(
    app: Path,
    samples: Path,
    model: Path,
    rows: int,
    cols: int,
    poles: int,
    terms: str,
) -> None:
    solver = {
        "samples": samples.name,
        "rows": rows,
        "cols": cols,
        "fit": {
            "poles": poles,
            "terms": terms,
            "weighting": "inverse_magnitude",
        },
        "output": {"model": model.name},
    }
    solver_file = model.with_suffix("").with_suffix(".solver.json")
    solver_file.write_text(json.dumps(solver, indent=2) + "\n")

    subprocess.run(
        [str(app), str(solver_file)],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    print(f"wrote {rows}x{cols} {poles}-pole model: {model}")


def main() -> None:
    args = parse_args()
    if not args.app.exists():
        raise RuntimeError(f"VectorFitting binary not found: {args.app}")

    fieldnames, rows = read_monitor_csv()
    k = conductor_count(fieldnames)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    yc_channels = [
        (f"Overhead_Yc_real_{i}_{j}", f"Overhead_Yc_imag_{i}_{j}")
        for i in range(k)
        for j in range(k)
    ]
    write_samples_csv(OUTPUT_DIR / "yc.csv", rows, yc_channels)
    tau_min = write_hmin_csv(OUTPUT_DIR / "hmin.csv", rows, k)
    (OUTPUT_DIR / "delay.json").write_text(json.dumps(tau_min, indent=2) + "\n")

    print(f"wrote {k}x{k} Yc samples and {k}x1 Hmin samples to {OUTPUT_DIR}")
    print(f"wrote modal tau_min delays: {OUTPUT_DIR / 'delay.json'}")

    fit_model(
        args.app, OUTPUT_DIR / "yc.csv", OUTPUT_DIR / "yc.model.json",
        k, k, YC_POLES, "constant",
    )
    fit_model(
        args.app, OUTPUT_DIR / "hmin.csv", OUTPUT_DIR / "hmin.model.json",
        k, 1, HMIN_POLES, "none",
    )


if __name__ == "__main__":
    main()
