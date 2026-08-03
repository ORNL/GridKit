#!/usr/bin/env python3
"""Prepare vector-fitting inputs from the Overhead sweep and fit them.

Reads the monitor CSV written by the FrequencyResponse application,
extracts the characteristic admittance and the minimum-phase-shifted
propagation function, writes the sampled-response CSVs consumed by the
committed solver specifications, and runs the VectorFitting application
on both.
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
INPUT_CSV = EXAMPLE_DIR / "output" / "overhead.response.csv"
OUTPUT_DIR = EXAMPLE_DIR / "output"
SOLVER_SPECS = [EXAMPLE_DIR / "yc.solver.json", EXAMPLE_DIR / "hmin.solver.json"]
DEFAULT_APP = HERE.parents[3] / "build" / "application" / "Fitting" / "VectorFitting"


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


def write_yc_csv(path: Path, rows: list[dict[str, str]], k: int) -> None:
    """Write Yc samples: omega, then one re/im pair per entry.

    The monitored characteristic admittance solves the sandwich equation
    Yc Z Yc = Y and is symmetric to sweep tolerance, so the samples pass
    through unchanged.
    """
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        header = ["omega"]
        for i in range(k):
            for j in range(k):
                header.extend([f"re_Yc_{i}_{j}", f"im_Yc_{i}_{j}"])
        writer.writerow(header)

        for row in rows:
            record = [row["omega"]]
            for i in range(k):
                for j in range(k):
                    record.extend(
                        [
                            row[f"Overhead_Yc_real_{i}_{j}"],
                            row[f"Overhead_Yc_imag_{i}_{j}"],
                        ]
                    )
            writer.writerow(record)


def write_hmin_csv(path: Path, rows: list[dict[str, str]], modes: int) -> list[float]:
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


def main() -> None:
    args = parse_args()
    if not args.app.exists():
        raise RuntimeError(f"VectorFitting binary not found: {args.app}")

    fieldnames, rows = read_monitor_csv()
    k = conductor_count(fieldnames)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_yc_csv(OUTPUT_DIR / "yc.csv", rows, k)
    tau_min = write_hmin_csv(OUTPUT_DIR / "hmin.csv", rows, k)
    (OUTPUT_DIR / "delay.json").write_text(json.dumps(tau_min, indent=2) + "\n")

    print(f"wrote {k}x{k} Yc samples and {k}x1 Hmin samples to {OUTPUT_DIR}")
    print(f"wrote modal tau_min delays: {OUTPUT_DIR / 'delay.json'}")

    for spec in SOLVER_SPECS:
        subprocess.run(
            [str(args.app), str(spec)], check=True, stdout=subprocess.DEVNULL
        )
        model = json.loads(spec.read_text())["output"]["model"]
        print(f"fitted {spec.name}: {EXAMPLE_DIR / model}")


if __name__ == "__main__":
    main()
