#!/usr/bin/env python3
"""Run the validation cases at several CommonMath smoothing scales."""

import argparse
import csv
import json
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common import plot as paper
from common import signals
from common.signals import column_key, load_csv as load_output
from common.simulation import ROOT, check_settings, monitor_signals, read_json
import numpy as np


CASES = (
    "ACTIVSg200",
    "ACTIVSg2000",
    "Hawaii",
    "IEEE39",
    "WECC240",
)
MUS = sorted({float(f"{mu:g}") for mu in np.geomspace(10, 10000, 40)} | {60.0, 240.0})
PLOT_CASES = tuple(paper.CASE_STYLES)
SIGNALS = (
    ("omega", r"$\Delta\omega$", "_omega"),
    ("p", r"$P$", "_p"),
    ("q", r"$Q$", "_q"),
    ("vmag", r"$|V|$", "_Vm"),
)
FAILURE_COLOR = "#a51c30"

FIELDNAMES = (
    "case",
    "mu",
    "signal",
    "series",
    "samples",
    "relative_max_error",
    "relative_l2_error",
    "absolute_max_error",
    "absolute_l2_error",
    "max_error_time",
    "wall_seconds",
    "status",
    "notes",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path,
                        default=ROOT / "build/paper/application/PhasorDynamics/DynamicSimulation")
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    parser.add_argument("--mus", nargs="+", type=float, default=MUS)
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).parent / "data/mu_sweep_errors.csv")
    parser.add_argument("--figure", type=Path,
                        default=Path(__file__).parent / "figures/mu_sweep_errors.png")
    parser.add_argument("--workdir", type=Path,
                        default=ROOT / "build/HICSS-60/mu-sweep/error-sensitivity")
    parser.add_argument("--run", action="store_true", help="Collect new data before plotting")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def select_output(columns, data, suffix):
    indices = [index for index, column in enumerate(columns)
               if index > 0 and column.endswith(suffix)]
    return [column_key(columns[index]) for index in indices], data[:, indices]


def load_reference(path):
    columns, data = signals.load_csv(path)
    return columns[1:], data[:, 0], data[:, 1:]


def error_record(case_name, mu, signal, output_time, output_values,
                 reference_time, reference_values, threshold, run_seconds):
    signals.check_grid(output_time, output_values, reference_time, reference_values)

    absolute_error = np.abs(output_values - reference_values)
    error_by_time = np.max(absolute_error, axis=1)
    reference_by_time = np.max(np.abs(reference_values), axis=1)
    absolute_max = float(np.max(error_by_time))
    absolute_l2 = float(np.linalg.norm(error_by_time))
    reference_max = float(np.max(reference_by_time))
    reference_l2 = float(np.linalg.norm(reference_by_time))
    relative_max = absolute_max / reference_max if reference_max > threshold else absolute_max
    relative_l2 = absolute_l2 / reference_l2 if reference_l2 > threshold else absolute_l2
    max_row = int(np.argmax(error_by_time))

    return {
        "case": case_name,
        "mu": f"{mu:g}",
        "signal": signal,
        "series": output_values.shape[1],
        "samples": output_values.shape[0],
        "relative_max_error": f"{relative_max:.12e}",
        "relative_l2_error": f"{relative_l2:.12e}",
        "absolute_max_error": f"{absolute_max:.12e}",
        "absolute_l2_error": f"{absolute_l2:.12e}",
        "max_error_time": f"{output_time[max_row]:.12e}",
        "wall_seconds": f"{run_seconds:.6f}",
        "status": "ok",
        "notes": "",
    }


def failure_record(case_name, mu, signal, run_seconds, returncode, notes):
    return {
        "case": case_name,
        "mu": f"{mu:g}",
        "signal": signal,
        "series": "",
        "samples": "",
        "relative_max_error": "nan",
        "relative_l2_error": "nan",
        "absolute_max_error": "nan",
        "absolute_l2_error": "nan",
        "max_error_time": "nan",
        "wall_seconds": f"{run_seconds:.6f}",
        "status": f"solver_failed_{returncode}",
        "notes": notes,
    }


def failure_note(path):
    lines = [re.sub(r"\x1b\[[0-9;]*m", "", line).strip()
             for line in path.read_text().splitlines() if line.strip()]
    for pattern in ("At t =", "Function ", "linesearch"):
        found = [line for line in lines if pattern in line]
        if found:
            return found[-1]
    return lines[-1] if lines else "No diagnostic output"


def sort_key(row):
    signals = [signal for signal, _label, _suffix in SIGNALS]
    return (CASES.index(row["case"]), float(row["mu"]), signals.index(row["signal"]))


def write_results(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_suffix(".csv.tmp")
    with pending.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(rows, key=sort_key))
    pending.replace(path)


def run_sweep(args):
    experiment_root = Path(__file__).resolve().parent
    binary = args.binary.resolve()
    if not binary.is_file():
        raise FileNotFoundError(binary)

    if args.data.exists() and not args.resume:
        raise FileExistsError(f"{args.data}: choose a new output or use --resume")
    rows = read_results(args.data) if args.resume and args.data.exists() else []
    for row in rows:
        row.setdefault("status", "ok")
        row.setdefault("notes", "")
    groups = {}
    for row in rows:
        groups.setdefault((row["case"], float(row["mu"])), []).append(row)
    completed = {key for key, records in groups.items()
                 if len(records) == len(SIGNALS)
                 and {r["signal"] for r in records} == {s[0] for s in SIGNALS}
                 and all(r["status"] == "ok" for r in records)}
    work_root = args.workdir.resolve()
    for case_name in args.cases:
        source_study = experiment_root / "solvers" / f"{case_name}.solver.json"
        study = json.loads(source_study.read_text())
        source_case = (source_study.parent / study["system_model_file"]).resolve()
        references = (source_study.parent / study["reference_file"]).resolve().parent
        case = json.loads(source_case.read_text())
        threshold = float(study.get("abs_err_threshold", np.finfo(float).eps))
        monitor_signals(case, ("p", "q", "omega"))

        work_dir = work_root / case_name
        work_dir.mkdir(parents=True, exist_ok=True)
        case_path = work_dir / source_case.name
        case_path.write_text(json.dumps(case, separators=(",", ":")))

        for mu in args.mus:
            if (case_name, mu) in completed:
                check_settings(read_json(work_dir / f"{case_name}.mu-{mu:g}.solver.json"),
                               dict(study, mu=mu))
                print(f"Skipping completed {case_name} at mu={mu:g}", flush=True)
                continue
            rows = [row for row in rows
                    if (row["case"], float(row["mu"])) != (case_name, mu)]
            output_path = work_dir / f"{case_name}.mu-{mu:g}.csv"
            study_path = work_dir / f"{case_name}.mu-{mu:g}.solver.json"
            log_path = work_dir / f"{case_name}.mu-{mu:g}.log"
            run_study = dict(study)
            run_study["system_model_file"] = case_path.name
            run_study["output_file"] = output_path.name
            run_study["mu"] = mu
            run_study.pop("reference_file", None)
            run_study.pop("error_tolerance", None)
            run_study.pop("error_type", None)
            run_study.pop("abs_err_threshold", None)
            study_path.write_text(json.dumps(run_study, indent=2) + "\n")

            print(f"Running {case_name} at mu={mu:g}...", flush=True)
            start = time.perf_counter()
            try:
                with log_path.open("w") as log:
                    result = subprocess.run(
                        (str(binary), study_path.name), cwd=work_dir, stdout=log,
                        stderr=subprocess.STDOUT, check=False, timeout=300)
                returncode = result.returncode
            except subprocess.TimeoutExpired:
                returncode = "timeout"
            run_seconds = time.perf_counter() - start
            if returncode != 0:
                notes = "Exceeded 300 s" if returncode == "timeout" else failure_note(log_path)
                print(
                    f"  solver failed with status {returncode}: {notes}",
                    flush=True,
                )
                rows.extend(
                    failure_record(case_name, mu, signal, run_seconds,
                                   returncode, notes)
                    for signal, _label, _suffix in SIGNALS
                )
                write_results(rows, args.data)
                if output_path.exists():
                    output_path.unlink()
                continue

            output_columns, output_data = load_output(output_path)
            output_time = output_data[:, 0]
            for signal, _label, suffix in SIGNALS:
                labels, output_values = select_output(output_columns, output_data, suffix)
                reference_labels, reference_time, reference_values = load_reference(
                    references / f"{case_name}.{signal}.ref.csv")
                output_values = signals.align(labels, output_values, reference_labels)
                record = error_record(
                    case_name, mu, signal, output_time, output_values,
                    reference_time, reference_values, threshold, run_seconds)
                rows.append(record)
                print(
                    f"  {signal:5s}: max={float(record['relative_max_error']):.6e}, "
                    f"L2={float(record['relative_l2_error']):.6e}",
                    flush=True,
                )
            write_results(rows, args.data)
            output_path.unlink()
    return rows


def read_results(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    keys = [(r["case"], float(r["mu"]), r["signal"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate case/mu/signal measurements")
    return rows


def plot_results(rows, path):
    lookup = {(row["case"], float(row["mu"]), row["signal"]): row for row in rows}
    if len(lookup) != len(rows):
        raise ValueError("Duplicate case/mu/signal measurements")
    for case, mu in {(r["case"], float(r["mu"])) for r in rows}:
        if any((case, mu, signal) not in lookup for signal, _, _ in SIGNALS):
            raise ValueError(f"{case} mu={mu}: incomplete signal measurements")
    case_names = [name for name in PLOT_CASES
                  if any(row["case"] == name for row in rows)]
    mus = sorted({float(row["mu"]) for row in rows if row["case"] in case_names})

    figure, axes = paper.subplots(2, 2, sharex=True, sharey=True, legend=True)
    handles = []
    for panel, (axis, (signal, label, _suffix)) in enumerate(zip(axes.ravel(), SIGNALS)):
        failed_mus = set()
        for case_name in case_names:
            display_name, color = paper.CASE_STYLES[case_name]
            found = [(mu, lookup.get((case_name, mu, signal))) for mu in mus]
            present = [(mu, row) for mu, row in found if row is not None]
            failed_mus.update(mu for mu, row in present
                              if row.get("status", "ok") != "ok")
            line, = axis.plot(
                [mu for mu, _row in present],
                [float(row["relative_max_error"]) for _mu, row in present],
                color=color, label=display_name,
                marker="o" if len(mus) <= 8 else None, markersize=3.0,
            )
            if panel == 0:
                handles.append(line)
        paper.set_mu_axis(axis, mus)
        axis.set_yscale("log")
        paper.panel_label(axis, panel, label)
        for mu in failed_mus:
            axis.plot(mu, 0.02, marker="x", color=FAILURE_COLOR, markersize=5,
                      markeredgewidth=1.0, transform=axis.get_xaxis_transform(),
                      clip_on=False)

    paper.legend(figure, handles)
    paper.labels(figure, x=r"$\mu$ – Smoothing scale", y="Max relative error")
    paper.save_figure(figure, path)


def main():
    args = parse_args()
    if len(set(args.mus)) != len(args.mus) or any(not np.isfinite(mu) or mu <= 0 for mu in args.mus):
        raise ValueError("mu values must be unique, positive, and finite")
    if args.resume and not args.run:
        raise ValueError("--resume requires --run")
    rows = read_results(args.data) if not args.run else run_sweep(args)
    plot_results(rows, args.figure)
    return int(any(row.get("status", "ok") != "ok" for row in rows))


if __name__ == "__main__":
    raise SystemExit(main())
