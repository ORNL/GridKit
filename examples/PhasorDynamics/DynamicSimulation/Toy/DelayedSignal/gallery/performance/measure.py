"""Refresh the application timings with the current DynamicSimulation build."""
import json
from pathlib import Path
import re
import statistics
import subprocess

HERE = Path(__file__).resolve().parent
REPO = next(path for path in HERE.parents if (path / "GridKit/Solver").is_dir())
results = json.loads((HERE / "results.json").read_text())
for case in results["cases"]:
    times = []
    for _ in range(3):
        run = subprocess.run(
            [str(REPO / "build/application/PhasorDynamics/DynamicSimulation"), "solver.json"],
            cwd=HERE / case["case"], text=True, capture_output=True, check=True
        )
        times.append(float(re.search(r"Complete in ([\d.e+-]+)", run.stdout + run.stderr)[1]))
    case["after_seconds"] = times
    case["after_median"] = statistics.median(times)
    case["speedup"] = case["before_median"] / case["after_median"]
    case["cleanup_speedup"] = case["cleanup_before_median"] / case["after_median"]
    print(case["case"], times)
(HERE / "results.json").write_text(json.dumps(results, indent=2) + "\n")
