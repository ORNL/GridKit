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
    "ACTIVSg500",
    "ACTIVSg2000",
    "ACTIVSg10k",
    "Hawaii",
    "IEEE39",
    "WECC240",
)
SIGNALS = (
    ("omega", r"$\Delta\omega$", "_omega"),
    ("p", r"$P$", "_p"),
    ("q", r"$Q$", "_q"),
    ("vmag", r"$|V|$", "_Vm"),
)
GENERATOR_CLASSES = {"genrou", "gensal", "genclassical"}
COLORS = ("#173665", "#079b88", "#d45b00", "#7851a9")
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


def write_results(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


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


def plot_results(rows, path):
    mus = sorted({float(row["mu"]) for row in rows})
    case_names = [name for name in CASES if any(row["case"] == name for row in rows)]
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Nimbus Roman", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "axes.edgecolor": "#333333",
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    figure, axes = plt.subplots(4, 2, figsize=(7.2, 9.2), sharex=True, sharey=True)
    axes = axes.ravel()
    handles = []
    for panel, (axis, case_name) in enumerate(zip(axes, case_names)):
        failed_mus = set()
        for color, (signal, label, _suffix) in zip(COLORS, SIGNALS):
            values = []
            for mu in mus:
                row = next(row for row in rows
                           if row["case"] == case_name
                           and float(row["mu"]) == mu
                           and row["signal"] == signal)
                values.append(100.0 * float(row["relative_max_error"]))
                if row.get("status", "ok") != "ok":
                    failed_mus.add(mu)
            line, = axis.plot(
                mus, values, color=color, marker="o", markersize=4.5,
                markeredgecolor="white", markeredgewidth=0.6, linewidth=1.8,
                label=label,
            )
            if panel == 0:
                handles.append(line)
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_xticks(mus, [f"{mu:g}" for mu in mus])
        axis.grid(True, color="#dddddd", linewidth=0.55)
        axis.set_axisbelow(True)
        axis.text(0.04, 0.94, f"({chr(ord('a') + panel)})", transform=axis.transAxes,
                  va="top", fontweight="bold")
        axis.text(0.96, 0.08, case_name, transform=axis.transAxes,
                  ha="right", va="bottom", fontsize=11, fontweight="bold")
        for mu in failed_mus:
            axis.plot(mu, 0.06, marker="x", color="#a51c30", markersize=6,
                      markeredgewidth=1.2, transform=axis.get_xaxis_transform(),
                      clip_on=False)
        if failed_mus:
            axis.text(0.04, 0.08, r"$\times$ solver failure", color="#a51c30",
                      transform=axis.transAxes, va="bottom", fontsize=7.5)

    for axis in axes[len(case_names):]:
        axis.axis("off")
    figure.legend(handles=handles, loc="upper center", ncol=len(SIGNALS),
                  frameon=False, bbox_to_anchor=(0.5, 0.995),
                  handlelength=2.4, columnspacing=1.8)
    figure.supxlabel(r"CommonMath smoothing scale, $\mu$", y=0.025)
    figure.supylabel(r"Maximum relative error [\%]", x=0.025)
    figure.subplots_adjust(top=0.94, bottom=0.07, left=0.12, right=0.98,
                           hspace=0.14, wspace=0.13)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=600, bbox_inches="tight", pad_inches=0.03)
    plt.close(figure)
    print(f"Wrote {path}")


def main():
    args = parse_args()
    rows = read_results(args.output) if args.plot_only else run_sweep(args)
    plot_results(rows, args.figure)


if __name__ == "__main__":
    main()
