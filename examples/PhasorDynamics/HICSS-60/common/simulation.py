"""Shared simulation helpers."""

import json
import math
import os
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
TRIALS = 7
GENERATORS = {"genrou", "gensal", "genclassical"}
SOLVER_DEFAULTS = {"dt_fixed": 0.0, "max_order": 5, "consistent_ic_type": "ya_ydp"}


def check_settings(recorded, expected):
    for key in ("tmax", "mu", "rel_tol", "abs_tol", "dt_monitor", "max_steps", "events", *SOLVER_DEFAULTS):
        if recorded.get(key, SOLVER_DEFAULTS.get(key)) != expected.get(key, SOLVER_DEFAULTS.get(key)):
            raise ValueError(f"Saved run does not match solver setting {key}")


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_suffix(path.suffix + ".tmp")
    pending.write_text(json.dumps(value, indent=4, allow_nan=False) + "\n")
    pending.replace(path)


def disable_monitors(case):
    case.pop("monitors", None)
    for element in [*case["buses"], *case["devices"]]:
        element["mon"] = []


def monitor_signals(case, generator_signals=("omega",)):
    disable_monitors(case)
    for bus in case["buses"]:
        bus["name"] = str(bus["number"])
        bus["mon"] = ["Vm"]
    for device in case["devices"]:
        if device["class"].lower() in GENERATORS:
            device["mon"] = list(generator_signals)


def simulation_result(log):
    counts = re.search(r"Variables: (\d+) differential, (\d+) algebraic", log)
    timing = re.search(r"Complete in ([\d.eE+-]+) seconds", log)
    if not counts or not timing:
        raise ValueError("Missing variable counts or CPU runtime")
    seconds = float(timing[1])
    if not math.isfinite(seconds) or seconds <= 0:
        raise ValueError("CPU runtime must be positive and finite")
    return tuple(map(int, counts.groups())), seconds


def simulate(binary, work, config, label):
    work.mkdir(parents=True, exist_ok=True)
    source = work / f"{label}.solver.json"
    write_json(source, config)
    log_path = work / f"{label}.log"
    with log_path.open("w") as log:
        subprocess.run([str(binary), source.name], cwd=work, stdout=log,
                       stderr=subprocess.STDOUT, check=True, timeout=300)
    return simulation_result(log_path.read_text())


def pin_cpu(cpu=None):
    allowed = os.sched_getaffinity(0)
    cpu = min(allowed) if cpu is None else cpu
    if cpu not in allowed:
        raise ValueError(f"CPU {cpu} is outside the allowed affinity")
    os.sched_setaffinity(0, {cpu})
    return cpu


def measure(binary, work, config, label, *, reuse=False):
    """One warm-up and seven CPU timings, or read those same recorded trials."""
    if not reuse:
        simulate(binary, work, config, f"{label}-warmup")
    trials = []
    for trial in range(1, TRIALS + 1):
        name = f"{label}-{trial}"
        if reuse:
            check_settings(read_json(work / f"{name}.solver.json"), config)
            result = simulation_result((work / f"{name}.log").read_text())
        else:
            result = simulate(binary, work, config, name)
        trials.append(result)
    if any(counts != trials[0][0] for counts, _ in trials):
        raise ValueError("Variable counts changed between timing trials")
    return {"counts": trials[0][0], "cpu_seconds": [seconds for _, seconds in trials]}
