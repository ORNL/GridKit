#!/usr/bin/env python3
"""Measure monitor-free median CPU times and plot runtime sensitivity to mu."""
import argparse
import csv
import subprocess
import sys
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import plot as paper
from common.simulation import (ROOT, TRIALS, check_settings, disable_monitors,
                               measure, pin_cpu, read_json, write_json)
import numpy as np

HERE = Path(__file__).resolve().parent
MUS = sorted({float(f"{mu:g}") for mu in np.geomspace(10, 10000, 40)} | {60.0, 240.0})
FIELDNAMES = ("case", "mu", "cpu", *[f"trial_{i}_seconds" for i in range(1, TRIALS + 1)],
              "cpu_seconds", "status", "notes")
PLOT_CASES = tuple(paper.CASE_STYLES)


def read_results(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def run_sweep(args):
    binary = args.binary.resolve(strict=True)
    cpu = pin_cpu(args.cpu)
    print(f"CPU affinity: {cpu}", flush=True)
    if args.data.exists() and not args.resume:
        raise FileExistsError(f"{args.data}: choose a new output or use --resume")
    rows = read_results(args.data) if args.resume and args.data.exists() else []
    collect_runs(rows)  # Reject duplicate keys before resuming.
    completed = {(row["case"], float(row["mu"])) for row in rows
                 if row["status"] == "ok"}
    args.data.parent.mkdir(parents=True, exist_ok=True)
    for name in args.cases:
        path = HERE / "solvers" / f"{name}.solver.json"
        study = read_json(path)
        benchmark = HERE.parents[1] / "benchmark/solvers" / path.name
        reference = read_json(benchmark)
        check_settings(study, reference)
        case_path = (path.parent / study["system_model_file"]).resolve()
        reference_path = (benchmark.parent / reference["system_model_file"]).resolve()
        if case_path != reference_path:
            raise ValueError(f"{name}: runtime sweep and benchmark use different cases")
        case = read_json(case_path)
        disable_monitors(case)
        work = (args.workdir / name).resolve()
        work.mkdir(parents=True, exist_ok=True)
        if args.resume and (work / "case.json").exists() and read_json(work / "case.json") != case:
            raise ValueError(f"{name}: recorded case differs from the runtime-sweep case")
        write_json(work / "case.json", case)
        for mu in args.mus:
            if (name, mu) in completed:
                measure(binary, work, dict(study, mu=mu), f"mu-{mu:g}", reuse=True)
                continue
            rows = [row for row in rows if (row["case"], float(row["mu"])) != (name, mu)]
            run_study = dict(study, system_model_file="case.json", mu=mu)
            times = []
            failure = ""
            try:
                times = measure(binary, work, run_study, f"mu-{mu:g}")["cpu_seconds"]
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired, ValueError) as error:
                failure = f"{type(error).__name__}: mu-{mu:g}"
            row = dict(case=name, mu=f"{mu:g}", cpu=cpu, status="failed" if failure else "ok",
                       notes=failure, cpu_seconds="nan" if failure else f"{median(times):.9g}")
            for trial in range(TRIALS):
                row[f"trial_{trial+1}_seconds"] = f"{times[trial]:.9g}" if trial < len(times) else "nan"
            rows.append(row)
            pending = args.data.with_suffix(".csv.tmp")
            with pending.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
                writer.writeheader()
                writer.writerows(sorted(rows, key=lambda row: (PLOT_CASES.index(row["case"]), float(row["mu"]))))
            pending.replace(args.data)
            print(f"{name} mu={mu:g}: {failure or row['cpu_seconds'] + ' s CPU'}", flush=True)
    return rows


def collect_runs(rows):
    runs = {}
    for row in rows:
        key = (row["case"], float(row["mu"]))
        value = (float(row["cpu_seconds"]), row.get("status", "ok"))
        if key in runs:
            raise ValueError(f"Duplicate run data for {key}")
        runs[key] = value
    return runs


def plot_runtime_results(rows, path):
    runs = collect_runs(rows)
    case_names = [name for name in PLOT_CASES
                  if any(case_name == name for case_name, _mu in runs)]
    mus = sorted({mu for case_name, mu in runs if case_name in case_names})

    figure, axis = paper.subplots(legend=True, figsize=(7.0, 3.5))
    for case_name in case_names:
        display_name, color = paper.CASE_STYLES[case_name]
        points = [
            (mu, runs[(case_name, mu)][0] if (case_name, mu) in runs
             and runs[(case_name, mu)][1] == "ok" else float("nan"))
            for mu in mus
        ]
        isolated = [i for i, (_, value) in enumerate(points)
                    if np.isfinite(value)
                    and (i == 0 or not np.isfinite(points[i - 1][1]))
                    and (i == len(points) - 1 or not np.isfinite(points[i + 1][1]))]
        axis.plot(
            [mu for mu, _runtime in points],
            [runtime for _mu, runtime in points],
            color=color, label=display_name,
        )
        axis.scatter([points[i][0] for i in isolated], [points[i][1] for i in isolated],
                     color=color, s=9)

    paper.set_mu_axis(axis, mus)
    axis.set_yscale("log")
    paper.labels(figure, x=r"$\mu$ – Smoothing scale", y="Median CPU runtime [sec]")
    paper.legend(figure)
    paper.save_figure(figure, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path,
                        default=ROOT / "build/paper/application/PhasorDynamics/DynamicSimulation")
    parser.add_argument("--workdir", type=Path, default=ROOT / "build/HICSS-60/mu-sweep/runtime-sensitivity")
    parser.add_argument("--mus", nargs="+", type=float, default=MUS)
    parser.add_argument("--cases", nargs="+", choices=PLOT_CASES, default=PLOT_CASES)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--run", action="store_true", help="Collect new data before plotting")
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).parent / "data/mu_sweep_runtime.csv")
    parser.add_argument("--figure", type=Path,
                        default=Path(__file__).parent / "figures/mu_sweep_runtime.png")
    args = parser.parse_args()
    if len(set(args.mus)) != len(args.mus) or any(not np.isfinite(mu) or mu <= 0 for mu in args.mus):
        parser.error("mu values must be unique, positive, and finite")
    if args.resume and not args.run:
        parser.error("--resume requires --run")
    rows = read_results(args.data) if not args.run else run_sweep(args)
    plot_runtime_results(rows, args.figure)
    return int(any(row["status"] != "ok" for row in rows))


if __name__ == "__main__":
    raise SystemExit(main())
