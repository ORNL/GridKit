#!/usr/bin/env python3
"""Render the DOE runtime table from saved GridKit timings."""

import argparse
import sys
from math import isfinite
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "HICSS-60"))
from common.simulation import read_json

HERE = Path(__file__).resolve().parent
CASES = ("Hawaii", "IEEE39", "ACTIVSg200", "WECC240", "ACTIVSg2000", "ACTIVSg10k", "ACTIVSg25k", "ACTIVSg70k")
# Fixed-step times: the five paper cases from the HICSS-60 benchmark table,
# ACTIVSg10k, ACTIVSg25k, and ACTIVSg70k supplied separately; not remeasured here.
FIXED_STEP_SECONDS = {
    "Hawaii": 2.172,
    "IEEE39": 6.969,
    "ACTIVSg200": 4.032,
    "WECC240": 4.375,
    "ACTIVSg2000": 14.703,
    "ACTIVSg10k": 62.735,
    "ACTIVSg25k": 109.5,
    "ACTIVSg70k": 385.562,
}
MISSING = "-"


def seconds(value):
    return f"{value:.3f} [s]"


def row(name, measurements, validation):
    gridkit, fixed_step, speedup, rmse, rel = MISSING, MISSING, MISSING, MISSING, MISSING
    if name in FIXED_STEP_SECONDS:
        fixed_step = seconds(FIXED_STEP_SECONDS[name])
    if name in measurements:
        trials = measurements[name]["cpu_seconds"]
        if not trials or any(not isfinite(t) or t <= 0 for t in trials):
            raise ValueError(f"{name}: expected positive CPU times")
        elapsed = median(trials)
        gridkit = seconds(elapsed)
        if name in FIXED_STEP_SECONDS:
            speedup = f"{FIXED_STEP_SECONDS[name] / elapsed:.1f}×"
    if name in validation:
        rmse, rel = (f"{value:.2e}" for value in validation[name])
    return f"| {name} | {gridkit} | {fixed_step} | {speedup} | {rmse} | {rel} |"


def timing_note(measurements):
    runs = {len(result["cpu_seconds"]) for result in measurements.values()}
    if len(runs) != 1:
        raise ValueError("cases have different numbers of timed runs")
    count = runs.pop()
    if count == 1:
        return "GridKit: CPU time of one run with monitoring off."
    return f"GridKit: median CPU time of {count} runs with monitoring off."


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=HERE / "data/measurements.json")
    parser.add_argument("--validation", type=Path, default=HERE / "data/validation.json")
    args = parser.parse_args()
    measurements = read_json(args.data)
    validation = read_json(args.validation)
    lines = [
        "# DOE",
        "",
        r"| Case | GridKit | Fixed-Step | Speedup | $\epsilon_{\mathrm{RMSE}}^{\text{abs}}$ | $\epsilon_{\infty}^{\text{rel}}$ |",
        "|:--|--:|--:|--:|--:|--:|",
        *(row(name, measurements, validation) for name in CASES),
        "",
        f"{timing_note(measurements)} Study settings are in [solvers](solvers/).",
        "",
        "Errors: generator speed against PowerWorld from [Validation](../Validation/).",
    ]
    (HERE / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
