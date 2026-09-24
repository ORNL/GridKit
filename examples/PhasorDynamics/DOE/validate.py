#!/usr/bin/env python3
"""Collect the Validation/ generator speed errors for the DOE table."""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "HICSS-60"))
from common import signals
from common.simulation import ROOT, read_json, write_json
from table import CASES, HERE

VALIDATION = Path("examples/PhasorDynamics/Validation")


def main():
    binary = ROOT / "build/application/PhasorDynamics/DynamicSimulation"
    errors = {}
    for name in (name for name in CASES if (ROOT / VALIDATION / name).is_dir()):
        source, work = ROOT / VALIDATION / name, ROOT / "build" / VALIDATION / name
        config = read_json(source / f"{name}.solver.json")
        subprocess.run([str(binary), f"{name}.solver.json"], cwd=work, capture_output=True, check=True)
        columns, output = signals.load_csv(work / config["output_file"])
        errors[name], _, _ = signals.compare(columns, output, source / config["reference_file"], "_omega")
        print(f"{name}: {errors[name]}", flush=True)
    write_json(HERE / "data/validation.json", errors)


if __name__ == "__main__":
    main()
