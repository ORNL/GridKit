#!/usr/bin/env python3
"""Collect statistics and traces once per contingency for Figures 4 and 5."""

import argparse
import csv
import math
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import plot as paper
from common import trace as ida
from common.simulation import ROOT, check_settings, disable_monitors, read_json, write_json


HERE = Path(__file__).resolve().parent
CASES = tuple(paper.CASE_STYLES)

COUNTERS = (
    "accepted_steps", "jacobian_evals", "residual_evals", "linear_solver_setups",
    "error_test_failures", "nonlinear_iterations", "nonlinear_convergence_failures",
)
FIELDS = ("case", "bus", "status", *COUNTERS, "diagnostic")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path,
                        default=ROOT / "build/paper/application/PhasorDynamics/ContingencyAnalysis")
    parser.add_argument("--cases", nargs="+", choices=CASES, default=CASES)
    parser.add_argument("--buses", nargs="+", type=int, help="Restrict to these bus numbers")
    parser.add_argument("--reactance", type=float, default=0.5, help="Fault reactance in p.u.")
    parser.add_argument("--workdir", type=Path, default=ROOT / "build/HICSS-60/contingency")
    parser.add_argument("--data", type=Path,
                        default=HERE / "data/contingency_stats.csv")
    parser.add_argument("--timeout", type=float, default=120.0,
                        help="Maximum seconds per bus run (default: 120)")
    parser.add_argument("--resume", action="store_true", help="Reuse verified runs with matching inputs")
    parser.add_argument("--retry-failed", action="store_true", help="Retry failed buses when resuming")
    args = parser.parse_args()
    if args.retry_failed and not args.resume:
        parser.error("--retry-failed requires --resume")
    if len(set(args.cases)) != len(args.cases):
        parser.error("--cases must not contain duplicates")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--timeout must be positive and finite")
    if not math.isfinite(args.reactance) or args.reactance <= 0:
        parser.error("--reactance must be positive and finite")
    return args


def collect_record(stats_path, bus, returncode, endpoints):
    document = read_json(stats_path)
    if document.get("schema_version") != 1 or len(document.get("records", [])) != 1:
        raise ValueError("Expected exactly one schema-version-1 fault record")
    record = document["records"][0]
    if record["bus"] != bus or record["fault_id"] != 0:
        raise ValueError("Statistics identify a different bus or fault")
    if record["status"] != "ok":
        return {"status": "failed", "diagnostic": record.get("diagnostic", "")}
    if returncode != 0:
        raise ValueError("Successful record with a failing application exit status")
    stats = record["stats"]
    segments = record["segments"]
    if [item["end_time"] for item in segments] != endpoints:
        raise ValueError("Missing or unexpected event segments")
    if [item["start_time"] for item in segments] != [0.0, *endpoints[:-1]]:
        raise ValueError("Event segments are not contiguous")
    for key in COUNTERS:
        values = [stats[key], *(segment[key] for segment in segments)]
        if any(type(value) is not int or value < 0 for value in values):
            raise ValueError(f"Invalid {key} counter")
        if stats[key] != sum(segment[key] for segment in segments):
            raise ValueError(f"Segment sum disagrees with total {key}")
    if stats["accepted_steps"] == 0:
        raise ValueError("Successful run has no accepted steps")
    return {"status": "ok", **stats, "diagnostic": ""}


def run_sweep(args, solvers=HERE / "solvers"):
    binary = args.binary.resolve(strict=True)
    resume = args.resume
    if args.data.exists() and not resume:
        raise FileExistsError("Results already exist; choose a new --data or use --resume")

    inputs = {}
    for name in args.cases:
        path = solvers / f"{name}.solver.json"
        study = read_json(path)
        case_path = (path.parent / study["system_model_file"]).resolve(strict=True)
        case = read_json(case_path)
        buses = [bus["number"] for bus in case["buses"]]
        if len(set(buses)) != len(buses):
            raise ValueError(f"Duplicate buses in {name}")
        if args.buses:
            missing = set(args.buses) - set(buses)
            if missing:
                raise ValueError(f"{name} has no buses {sorted(missing)}")
            buses = [bus for bus in buses if bus in args.buses]
        inputs[name] = {"solver": study, "case_file": str(case_path.relative_to(ROOT)),
                        "buses": buses}

    rows = []
    if resume and args.data.exists():
        with args.data.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        keys = [(row["case"], int(row["bus"])) for row in rows]
        if len(keys) != len(set(keys)):
            raise ValueError("Duplicate case/bus measurements")
    existing = {(row["case"], int(row["bus"])): row for row in rows}
    args.data.parent.mkdir(parents=True, exist_ok=True)

    def save_row(row):
        existing[(row["case"], int(row["bus"]))] = row
        pending = args.data.with_suffix(".csv.tmp")
        with pending.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(existing.values())
        pending.replace(args.data)

    for name, spec in inputs.items():
        case = read_json(ROOT / spec["case_file"])
        disable_monitors(case)
        case["devices"] = [device for device in case["devices"]
                           if device["class"].lower() != "busfault"]
        fault = {"class": "BusFault", "id": "paper_fault", "ports": {},
                 "params": {"state0": False, "R": 0.0, "X": args.reactance}}
        case["devices"].append(fault)
        work = (args.workdir / name).resolve()
        work.mkdir(parents=True, exist_ok=True)
        if resume and (work / "case.json").exists():
            saved_case = read_json(work / "case.json")
            saved_case["devices"][-1]["ports"] = {}
            if saved_case != case:
                raise ValueError(f"{name}: saved case or fault settings differ")
        endpoints = [event["time"] for event in spec["solver"]["events"]]
        endpoints.append(spec["solver"]["tmax"])
        for index, bus in enumerate(spec["buses"], 1):
            key = (name, bus)
            old = existing.get(key)
            if old and (old["status"] == "ok" or not args.retry_failed):
                check_settings(read_json(work / f"bus-{bus}.solver.json"), spec["solver"])
                if old["status"] == "ok":
                    record = collect_record(work / f"bus-{bus}.stats.json", bus,
                                            0, endpoints)
                    if record["status"] != "ok" or any(str(record[k]) != str(old[k]) for k in COUNTERS):
                        raise ValueError(f"{key}: saved statistics changed")
                    if "solver_trace_file" in spec["solver"]:
                        validate_trace(work / f"bus-{bus}.trace.csv", record, endpoints)
                print(f"{name}: bus {bus}: reused", flush=True)
                continue
            fault["ports"]["bus"] = bus
            write_json(work / "case.json", case)
            study = dict(spec["solver"], system_model_file="case.json",
                         contingency_stats_file=f"bus-{bus}.stats.json")
            if "solver_trace_file" in study:
                study["solver_trace_file"] = f"bus-{bus}.trace.csv"
                (work / study["solver_trace_file"]).unlink(missing_ok=True)
            study_path = work / f"bus-{bus}.solver.json"
            stats_path = work / study["contingency_stats_file"]
            log_path = work / f"bus-{bus}.log"
            stats_path.unlink(missing_ok=True)
            write_json(study_path, study)
            row = {"case": name, "bus": bus}
            try:
                with log_path.open("w") as log:
                    result = subprocess.run(
                        [str(binary), study_path.name], cwd=work, stdout=log,
                        stderr=subprocess.STDOUT, timeout=args.timeout, check=False)
                row.update(collect_record(stats_path, bus, result.returncode, endpoints))
                if row["status"] == "ok" and "solver_trace_file" in study:
                    validate_trace(work / study["solver_trace_file"], row, endpoints)
            except subprocess.TimeoutExpired:
                row.update(status="timeout", diagnostic=f"Exceeded {args.timeout:g} s")
            except (OSError, ValueError, KeyError, TypeError) as error:
                row.update(status="invalid_stats", diagnostic=str(error))
            if row["status"] != "ok":
                lines = log_path.read_text().splitlines() if log_path.exists() else []
                detail = next((line.strip() for line in reversed(lines)
                               if "failed with flag" in line or "At t =" in line), "")
                if detail:
                    row["diagnostic"] = "; ".join(filter(None, (row["diagnostic"], detail)))
                row["diagnostic"] = " ".join(
                    re.sub(r"\x1b\[[0-9;]*m", "", row["diagnostic"]).split())
            save_row(row)
            if index % 25 == 0 or index == len(spec["buses"]) or row["status"] != "ok":
                print(f"{name}: {index}/{len(spec['buses'])}, bus {bus}: {row['status']}", flush=True)
    return list(existing.values())


def validate_trace(path, record, endpoints):
    ends = [0.0, *endpoints]
    trace = ida.load(path, ends)
    _, counts = ida.cumulative_counts(trace, ends)
    if counts[-1].tolist() != [int(record[key]) for key in ida.COUNTERS]:
        raise ValueError(f"{path}: trace totals disagree with contingency statistics")


def main():
    rows = run_sweep(parse_args())
    return int(any(row["status"] != "ok" for row in rows))


if __name__ == "__main__":
    raise SystemExit(main())
