#!/usr/bin/env python3
"""Run sequential, interleaved mu trials and tighter-tolerance comparisons."""

import argparse
import hashlib
import json
import platform
import re
import shutil
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_environment(binary):
    result = {}
    if shutil.which('ldd'):
        linked = subprocess.run(['ldd', str(binary)], text=True, capture_output=True)
        if linked.returncode == 0:
            libraries = re.findall(r'=> (/\S+)', linked.stdout)
            result['shared_libraries_sha256'] = {path: digest(Path(path)) for path in libraries}
    cpuinfo = Path('/proc/cpuinfo')
    if cpuinfo.exists():
        match = re.search(r'^model name\s*:\s*(.*)', cpuinfo.read_text(), re.MULTILINE)
        if match:
            result['cpu_model'] = match[1]
    cache_path = next((parent / 'CMakeCache.txt' for parent in binary.parents
                       if (parent / 'CMakeCache.txt').exists()), None)
    if cache_path is not None:
        cache = dict(re.findall(r'^(\w+):[^=]*=(.*)$', cache_path.read_text(), re.MULTILINE))
        result['build_type'] = cache.get('CMAKE_BUILD_TYPE')
        result['enzyme_enabled'] = cache.get('GridKit_ENABLE_ENZYME')
        sundials_dir = cache.get('SUNDIALS_DIR')
        if sundials_dir:
            config = next((parent / 'include/sundials/sundials_config.h'
                           for parent in Path(sundials_dir).parents
                           if (parent / 'include/sundials/sundials_config.h').exists()), None)
            if config is not None:
                version = re.search(r'#define SUNDIALS_VERSION "([^"]+)"', config.read_text())
                if version:
                    result['sundials_version'] = version[1]
        compiler = cache.get('CMAKE_CXX_COMPILER')
        if compiler and shutil.which(compiler):
            version = subprocess.run([compiler, '--version'], text=True, capture_output=True, check=True)
            result['compiler'] = version.stdout.splitlines()[0]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build/emt-rational/application/EMT/EMTDynamicSimulation")
    parser.add_argument("--output", type=Path, default=ROOT / "build/emt-coupled-results")
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--names", nargs="+", choices=("low", "middle", "high"), default=["low", "middle", "high"])
    parser.add_argument("--skip-tight", action="store_true")
    args = parser.parse_args()
    binary, output = args.binary.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.trials < 1:
        parser.error("--trials must be positive")
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True)
    manifest = {"frequency_hz": 60., "switching_frequency_hz": 900., "voltage_ll_v": 13800.,
                "remote_bus": "bus_10", "events": [{"time": .5, "label": "Load closes"},
                                                     {"time": .8, "label": "Load opens"}],
                "runs": [], "provenance": {"binary": str(binary), "binary_sha256": digest(binary),
                   "case_sha256": digest(HERE / "CoupledGrid.case.json"),
                   "state_sha256": digest(HERE / "CoupledGrid.state.json"),
                   "platform": platform.platform(),
                   "git_revision": revision.stdout.strip() if revision.returncode == 0 else None,
                   "timing": "Sequential rotated order. Complete in is application CPU seconds including consistent initialization, integration and monitoring; wall includes model construction. Tight runs excluded from timing medians.",
                   **build_environment(binary)}}
    records = {}
    for name in args.names:
        study = json.loads((HERE / f"{name}.solver.json").read_text())
        record = {"name": name, "label": name, "mu": study["mu"], "csv": f"{name}/monitor.csv",
                  "complete_seconds": [], "wall_seconds": [], "trial_stats": [],
                  "rel_tol": study["rel_tol"], "abs_tol": study["abs_tol"], "dt_monitor": study["dt_monitor"]}
        manifest["runs"].append(record)
        records[name] = record

    def save():
        (output / "summary.json").write_text(json.dumps(manifest, indent=2) + "\n")

    def simulate(name, trial, tight=False):
        record = records[name]
        folder = output / name
        folder.mkdir(exist_ok=True)
        tag = "tight" if tight else f"trial_{trial + 1}"
        study = json.loads((HERE / f"{name}.solver.json").read_text())
        csv_name = "tight.csv" if tight else "monitor.csv" if trial == 0 else "repeat.csv"
        study.update(system_model_file=str(HERE / "CoupledGrid.case.json"),
                     state_file=str(HERE / "CoupledGrid.state.json"), output_file=str(folder / csv_name))
        if tight:
            study["rel_tol"] *= .1
            study["abs_tol"] *= .1
        study_path = folder / f"{tag}.solver.json"
        study_path.write_text(json.dumps(study, indent=2) + "\n")
        print(f"Starting {name}, {tag}, mu={record['mu']:.6f}", flush=True)
        start = time.perf_counter()
        with (folder / f"{tag}.log").open("w") as log:
            subprocess.run([str(binary), str(study_path)], cwd=folder, stdout=log, stderr=subprocess.STDOUT, check=True)
        elapsed = time.perf_counter() - start
        text = (folder / f"{tag}.log").read_text()
        match = re.search(r"Complete in ([\d.eE+-]+) seconds", text)
        stats_match = re.search(r"IDA statistics: (.*)", text)
        if match is None or stats_match is None or "IDA ERROR" in text:
            raise RuntimeError(f"Incomplete simulation; inspect {folder / (tag + '.log')}")
        seconds = float(match[1])
        stats = {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", stats_match[1])}
        if tight:
            record.update(tight_csv=f"{name}/tight.csv", tight_complete_seconds=seconds,
                          tight_wall_seconds=elapsed, tight_stats=stats,
                          tight_rel_tol=study["rel_tol"], tight_abs_tol=study["abs_tol"])
        else:
            record["complete_seconds"].append(seconds)
            record["wall_seconds"].append(elapsed)
            record["trial_stats"].append(stats)
            record["stats"] = record["trial_stats"][0]
            if trial:
                if digest(folder / csv_name) != digest(folder / "monitor.csv"):
                    raise RuntimeError(f"Nondeterministic repeated output: {name}, {tag}")
                (folder / csv_name).unlink()
        size = re.search(r"DAE variables: (\d+), Jacobian nonzeros: (\d+)", text)
        if size:
            record.update(dae_variables=int(size[1]), jacobian_nonzeros=int(size[2]))
        save()
        print(f"Finished {name}, {tag}: CPU={seconds:.6f}s wall={elapsed:.6f}s steps={stats['steps']}", flush=True)

    for trial in range(args.trials):
        order = args.names[trial % len(args.names):] + args.names[:trial % len(args.names)]
        for name in order:
            simulate(name, trial)
    if not args.skip_tight:
        for name in args.names:
            simulate(name, 0, tight=True)


if __name__ == "__main__":
    main()
