#!/usr/bin/env python3
"""Run the validation cases at several CommonMath smoothing scales."""

import argparse
import csv
import json
import subprocess
import tempfile
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CASES = (
    "ACTIVSg200",
    "ACTIVSg2000",
    "Hawaii",
    "IEEE39",
    "WECC240",
)
PLOT_CASES = ("Hawaii", "IEEE39", "ACTIVSg200", "ACTIVSg2000", "WECC240")
SIGNALS = (
    ("omega", r"$\Delta\omega$", "_omega"),
    ("p", r"$P$", "_p"),
    ("q", r"$Q$", "_q"),
    ("vmag", r"$|V|$", "_Vm"),
)
GENERATOR_CLASSES = {"genrou", "gensal", "genclassical"}
CASE_STYLES = {
    "Hawaii": ("Hawaii", "#3b6ea5"),
    "IEEE39": ("NE", "#2f9184"),
    "ACTIVSg200": ("Illinois", "#8663a6"),
    "ACTIVSg2000": ("Texas", "#cc7330"),
    "WECC240": ("WECC", "#b84a47"),
}
GRID_COLOR = "#dde2e7"
FAILURE_COLOR = "#a51c30"
RC_PARAMS = {
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "Times New Roman", "Nimbus Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 15.0,
    "axes.labelsize": 15.0,
    "axes.edgecolor": "#222222",
    "axes.linewidth": 1.0,
    "xtick.labelsize": 14.0,
    "ytick.labelsize": 14.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 5.0,
    "ytick.major.size": 5.0,
    "xtick.minor.size": 2.8,
    "ytick.minor.size": 2.8,
    "legend.fontsize": 15.0,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}
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
    "run_seconds",
    "status",
    "notes",
)


def parse_args():
    root = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path,
                        default=root / "build/application/PhasorDynamics/DynamicSimulation")
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    parser.add_argument("--mus", nargs="+", type=float, default=(60.0, 120.0, 240.0))
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("mu_sweep_errors.csv"))
    parser.add_argument("--figure", type=Path,
                        default=Path(__file__).with_name("mu_sweep_errors.png"))
    parser.add_argument("--runtime-figure", type=Path,
                        default=Path(__file__).with_name("mu_sweep_runtime.png"))
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def prepare_case(case):
    for bus in case["buses"]:
        bus["mon"] = ["Vm"]
    for device in case["devices"]:
        device["mon"] = (["p", "q", "omega"]
                         if device["class"].lower() in GENERATOR_CLASSES else [])


def column_key(column):
    parts = column.split("_")
    body = parts[1:-1]
    if body and body[-1].lower() == parts[0].lower():
        body.pop()
    return " ".join(body)


def load_output(path):
    with path.open() as stream:
        columns = stream.readline().strip().split(",")
    return columns, np.loadtxt(path, delimiter=",", skiprows=1)


def select_output(columns, data, suffix):
    indices = [index for index, column in enumerate(columns)
               if index > 0 and column.endswith(suffix)]
    return [column_key(columns[index]) for index in indices], data[:, indices]


def load_reference(path):
    with path.open() as stream:
        labels = stream.readline().strip().split(",")[1:]
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return labels, data[:, 0], data[:, 1:]


def align_output(signal, labels, values, reference_labels, bus_numbers):
    if signal == "vmag":
        if reference_labels != bus_numbers:
            raise ValueError("Bus-voltage reference order does not match the case bus order")
        if len(labels) != len(reference_labels):
            raise ValueError("Bus-voltage output and reference column counts differ")
        return values

    if labels == reference_labels:
        return values
    if sorted(labels) != sorted(reference_labels):
        raise ValueError(f"{signal} output and reference labels differ")
    return values[:, [labels.index(label) for label in reference_labels]]


def error_record(case_name, mu, signal, output_time, output_values,
                 reference_time, reference_values, threshold, run_seconds):
    if output_values.shape != reference_values.shape:
        raise ValueError(
            f"{case_name} {signal}: shape mismatch {output_values.shape} != {reference_values.shape}")
    if not np.allclose(output_time, reference_time, rtol=0.0, atol=1.0e-6):
        raise ValueError(f"{case_name} {signal}: output and reference time grids differ")

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
        "run_seconds": f"{run_seconds:.6f}",
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
        "run_seconds": f"{run_seconds:.6f}",
        "status": f"solver_failed_{returncode}",
        "notes": notes,
    }


def sort_key(row):
    signals = [signal for signal, _label, _suffix in SIGNALS]
    return (CASES.index(row["case"]), float(row["mu"]), signals.index(row["signal"]))


def write_results(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(sorted(rows, key=sort_key))


def run_sweep(args):
    root = Path(__file__).resolve().parents[3]
    validation_root = Path(__file__).resolve().parent
    binary = args.binary.resolve()
    if not binary.is_file():
        raise FileNotFoundError(binary)

    rows = read_results(args.output) if args.resume and args.output.exists() else []
    for row in rows:
        row.setdefault("status", "ok")
        row.setdefault("notes", "")
    completed = {(row["case"], float(row["mu"])) for row in rows}
    with tempfile.TemporaryDirectory(prefix="gridkit-mu-sweep-") as temporary:
        work_root = Path(temporary)
        for case_name in args.cases:
            source_case = root / "cases/PhasorDynamics" / case_name / f"{case_name}.case.json"
            source_study = validation_root / case_name / f"{case_name}.solver.json"
            references = validation_root / case_name / "reference"
            case = json.loads(source_case.read_text())
            study = json.loads(source_study.read_text())
            threshold = float(study.get("abs_err_threshold", np.finfo(float).eps))
            bus_numbers = [str(bus["number"]) for bus in case["buses"]]
            prepare_case(case)

            work_dir = work_root / case_name
            work_dir.mkdir()
            case_path = work_dir / source_case.name
            case_path.write_text(json.dumps(case, separators=(",", ":")))

            for mu in args.mus:
                if (case_name, mu) in completed:
                    print(f"Skipping completed {case_name} at mu={mu:g}", flush=True)
                    continue
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
                with log_path.open("w") as log:
                    result = subprocess.run(
                        (str(binary), study_path.name),
                        cwd=work_dir,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=False,
                    )
                run_seconds = time.perf_counter() - start
                if result.returncode != 0:
                    log_lines = log_path.read_text().splitlines()
                    notes = next(
                        (line.strip() for line in reversed(log_lines) if "At t =" in line),
                        log_lines[-1].strip() if log_lines else "No diagnostic output",
                    )
                    print(
                        f"  solver failed with status {result.returncode}: {notes}",
                        flush=True,
                    )
                    rows.extend(
                        failure_record(case_name, mu, signal, run_seconds,
                                       result.returncode, notes)
                        for signal, _label, _suffix in SIGNALS
                    )
                    write_results(rows, args.output)
                    if output_path.exists():
                        output_path.unlink()
                    continue

                output_columns, output_data = load_output(output_path)
                output_time = output_data[:, 0]
                for signal, _label, suffix in SIGNALS:
                    labels, output_values = select_output(output_columns, output_data, suffix)
                    reference_labels, reference_time, reference_values = load_reference(
                        references / f"{case_name}.{signal}.ref.csv")
                    output_values = align_output(
                        signal, labels, output_values, reference_labels, bus_numbers)
                    record = error_record(
                        case_name, mu, signal, output_time, output_values,
                        reference_time, reference_values, threshold, run_seconds)
                    rows.append(record)
                    print(
                        f"  {signal:5s}: max={float(record['relative_max_error']):.6e}, "
                        f"L2={float(record['relative_l2_error']):.6e}",
                        flush=True,
                    )
                write_results(rows, args.output)
                output_path.unlink()
    return rows


def read_results(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def set_mu_axis(axis, mus):
    """Label every sampled mu when they are few; fall back to decades when dense."""
    from matplotlib.ticker import (FixedFormatter, FixedLocator, LogFormatterMathtext,
                                   LogLocator, NullFormatter, NullLocator)

    axis.set_xscale("log")
    axis.set_xlim(min(mus), max(mus))
    if len(mus) > 8:
        axis.xaxis.set_major_locator(LogLocator(base=10.0))
        axis.xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
        axis.xaxis.set_minor_locator(LogLocator(base=10.0, subs=tuple(range(2, 10))))
        axis.xaxis.set_minor_formatter(NullFormatter())
        return
    axis.xaxis.set_major_locator(FixedLocator(mus))
    axis.xaxis.set_major_formatter(FixedFormatter([f"{mu:g}" for mu in mus]))
    axis.xaxis.set_minor_locator(NullLocator())


def style_axis(axis):
    axis.grid(True, which="major", color=GRID_COLOR, linewidth=0.6)
    axis.set_axisbelow(True)
    axis.tick_params(which="both", top=True, right=True)


def save_figure(figure, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=600)
    outputs = [path]
    if path.suffix.lower() == ".png":
        pdf = path.with_suffix(".pdf")
        figure.savefig(pdf)
        outputs.append(pdf)
    plt.close(figure)
    print("Wrote " + " and ".join(str(output) for output in outputs))


def plot_results(rows, path):
    lookup = {(row["case"], float(row["mu"]), row["signal"]): row for row in rows}
    case_names = [name for name in PLOT_CASES
                  if any(row["case"] == name for row in rows)]
    mus = sorted({float(row["mu"]) for row in rows if row["case"] in case_names})
    plt.rcParams.update(RC_PARAMS)

    figure, axes = plt.subplots(2, 2, figsize=(7.0, 7.7), sharex=True, sharey=True)
    handles = []
    for panel, (axis, (signal, label, _suffix)) in enumerate(zip(axes.ravel(), SIGNALS)):
        failed_mus = set()
        for case_name in case_names:
            display_name, color = CASE_STYLES[case_name]
            found = [(mu, lookup.get((case_name, mu, signal))) for mu in mus]
            present = [(mu, row) for mu, row in found if row is not None]
            failed_mus.update(mu for mu, row in present
                              if row.get("status", "ok") != "ok")
            line, = axis.plot(
                [mu for mu, _row in present],
                [float(row["relative_max_error"]) for _mu, row in present],
                color=color, linewidth=2.2, solid_capstyle="round", label=display_name,
            )
            if panel == 0:
                handles.append(line)
        set_mu_axis(axis, mus)
        axis.set_yscale("log")
        style_axis(axis)
        letter = chr(ord("a") + panel)
        axis.text(0.965, 0.94, rf"$\mathbf{{({letter})}}$ {label}",
                  transform=axis.transAxes, ha="right", va="top", fontsize=15.0)
        for mu in failed_mus:
            axis.plot(mu, 0.02, marker="x", color=FAILURE_COLOR, markersize=5,
                      markeredgewidth=1.0, transform=axis.get_xaxis_transform(),
                      clip_on=False)

    figure.legend(handles=handles, loc="upper center", ncol=len(case_names),
                  frameon=False, bbox_to_anchor=(0.5, 0.985),
                  handlelength=1.5, columnspacing=1.3)
    for axis in axes[-1, :]:
        axis.set_xlabel(r"$\mu$ - Smoothing scale")
    figure.supylabel("Max relative error", x=0.025)
    # The mu axis ends on real data, so the end labels sit at the panel edges
    # and need room between the columns.
    figure.subplots_adjust(top=0.90, bottom=0.105, left=0.125, right=0.975,
                           hspace=0.18, wspace=0.16)
    save_figure(figure, path)


def collect_runs(rows):
    runs = {}
    for row in rows:
        key = (row["case"], float(row["mu"]))
        value = (float(row["run_seconds"]), row.get("status", "ok"))
        if key in runs and runs[key] != value:
            raise ValueError(f"Inconsistent run data for {key}")
        runs[key] = value
    return runs


def plot_runtime_results(rows, path):
    runs = collect_runs(rows)
    case_names = [name for name in PLOT_CASES
                  if any(case_name == name for case_name, _mu in runs)]
    mus = sorted({mu for case_name, mu in runs if case_name in case_names})
    plt.rcParams.update(RC_PARAMS)

    figure, axis = plt.subplots(figsize=(7.0, 3.5))
    for case_name in case_names:
        display_name, color = CASE_STYLES[case_name]
        points = [
            (mu, runs[(case_name, mu)][0])
            for mu in mus
            if (case_name, mu) in runs and runs[(case_name, mu)][1] == "ok"
        ]
        axis.plot(
            [mu for mu, _runtime in points],
            [runtime for _mu, runtime in points],
            color=color, linewidth=2.2, solid_capstyle="round", label=display_name,
        )

    set_mu_axis(axis, mus)
    axis.set_yscale("log")
    axis.set_xlabel(r"$\mu$ - Smoothing scale")
    axis.set_ylabel("Elapsed runtime [sec]")
    style_axis(axis)

    figure.legend(loc="upper center", ncol=len(case_names), frameon=False,
                  bbox_to_anchor=(0.5, 0.985), handlelength=1.5, columnspacing=1.3)
    figure.subplots_adjust(top=0.78, bottom=0.20, left=0.13, right=0.985)
    save_figure(figure, path)


def main():
    args = parse_args()
    rows = read_results(args.output) if args.plot_only else run_sweep(args)
    plot_results(rows, args.figure)
    plot_runtime_results(rows, args.runtime_figure)


if __name__ == "__main__":
    main()
