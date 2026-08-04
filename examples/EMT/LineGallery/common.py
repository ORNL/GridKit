#!/usr/bin/env python3
"""Shared paths, sweep configuration, and application invocation for the
line gallery scripts.

Every generated artifact for a line lives in output/<line>/; the catalog
line descriptions live in ../Lines and gallery-only variants live beside
the solver specs in this directory.
"""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output"
BUILD = HERE.parents[2] / "build"

SWEEP_APP = BUILD / "application" / "EMT" / "FrequencyResponse" / "FrequencyResponse"
FIT_APP = BUILD / "application" / "Fitting" / "VectorFitting"
ULM_APP = BUILD / "application" / "EMT" / "UniversalLineModel" / "UniversalLineModel"

FREQUENCY = {"start": 10.0, "stop": 1.0e8, "points": 401, "scale": "log"}

# The propagation error target reflects what factor-based fitting supports
# for untransposed lines with nearly degenerate aerial modes.
ULM_ARGS = ["--h-target", "0.05"]
ULM_LINE_ARGS = {"500kv-double-circuit": ["--yc-max-poles", "16"]}


def line_model_file(name: str) -> Path:
    """Catalog lines live in ../Lines; gallery-only variants live here."""
    local = HERE / f"{name}.line.json"
    return local if local.exists() else HERE.parent / "Lines" / f"{name}.line.json"


def gallery_lines() -> list[str]:
    """Every line with a committed gallery solver specification."""
    return sorted(spec.name.replace(".solver.json", "")
                  for spec in HERE.glob("*.solver.json"))


def line_output_dir(name: str) -> Path:
    line_dir = OUTPUT_DIR / name
    line_dir.mkdir(parents=True, exist_ok=True)
    return line_dir


def read_modal_response(path: Path):
    """Monitored transforms and modal propagation from a sweep CSV.

    Returns omega, the conductor count, Tv, Ti, the modal propagation H,
    and the modal delay traces tau.
    """
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        names = reader.fieldnames or []
    omega = np.array([float(r["omega"]) for r in rows])
    m = len(rows)

    k = 0
    while f"Overhead_H_real_{k}" in names:
        k += 1
    if k == 0:
        raise ValueError(f"{path} lacks Overhead_H columns")

    def matrix(prefix: str) -> np.ndarray:
        out = np.zeros((m, k, k), complex)
        for s, r in enumerate(rows):
            for i in range(k):
                for j in range(k):
                    out[s, i, j] = complex(
                        float(r[f"{prefix}_real_{i}_{j}"]),
                        float(r[f"{prefix}_imag_{i}_{j}"]))
        return out

    h = np.array([[complex(float(r[f"Overhead_H_real_{q}"]),
                           float(r[f"Overhead_H_imag_{q}"]))
                   for q in range(k)] for r in rows])
    tau = np.array([[float(r[f"Overhead_Tau_{q}"]) for q in range(k)]
                    for r in rows])
    return omega, k, matrix("Overhead_Tv"), matrix("Overhead_Ti"), h, tau


def run_ulm(name: str,
            ulm_app: Path = ULM_APP) -> tuple[Path, subprocess.CompletedProcess]:
    """Fit one line into output/<line>/, echoing the fit report lines.

    Exit codes 0 (targets met, passive), 2 (target missed), and
    3 (nonpassive) all leave usable artifacts; anything else raises.
    """
    line_dir = line_output_dir(name)
    command = [
        str(ulm_app),
        str(line_model_file(name)),
        "-o",
        str(line_dir),
        *ULM_ARGS,
        *ULM_LINE_ARGS.get(name, []),
    ]
    result = subprocess.run(command, capture_output=True, text=True, cwd=HERE)
    for line in result.stdout.splitlines():
        if "VectorFitting:" in line or line.startswith("Propagation:"):
            print(line.strip())
    if result.returncode not in (0, 2, 3) or "failed" in result.stdout + result.stderr:
        raise RuntimeError(
            f"UniversalLineModel failed for {name} "
            f"(exit {result.returncode}):\n{result.stdout}\n{result.stderr}")
    return line_dir, result
