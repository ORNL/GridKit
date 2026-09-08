#!/usr/bin/env python3
"""Generate and run adaptive EMT line comparisons, retaining all study artifacts."""

import argparse
import copy
import hashlib
import json
import math
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def diag(value):
    return [[value if i == j else 0.0 for j in range(3)] for i in range(3)]


def fit(value):
    return {"rows": 3, "cols": 3, "D": diag(value)}


def history():
    return {
        "omega": 0.0,
        "i_ref1": [0.0] * 3,
        "i_ref2": [0.0] * 3,
        "d_i_ref1": [0.0] * 3,
        "d_i_ref2": [0.0] * 3,
    }


def bus(name):
    return {"class": "Bus", "id": name, "mon": ["va", "vb", "vc"]}


def distributed(name, a, b, yc, h, capped=False):
    h = copy.deepcopy(h)
    h["limit_step"] = capped
    return {
        "class": "LineDistributed",
        "id": name,
        "params": {"N": 3, "K": 3, "conductors": [1, 2, 3]},
        "inputs": {"bus1": a, "bus2": b},
        "submodels": {"Yc": yc, "H": h},
        "mon": ["i_inc1a", "i_inc2a", "i_ref1a", "i_ref2a"],
    }


def pi_line(name, a, b, length, rlgc):
    return {
        "class": "LineLumped",
        "id": name,
        "params": {
            "N": 3,
            "K": 3,
            "conductors": [1, 2, 3],
            "dx": length,
            **{key + "p": value for key, value in rlgc.items()},
        },
        "inputs": {"bus1": a, "bus2": b},
        "mon": ["i12a"],
    }


def bergeron(kind, sine=False, short=False):
    tau, zc = (20e-6 if short else 300e-6), 300.0
    dev = [bus("b0"), bus("b1")]
    signals = []
    if sine:
        dev.append(
            {
                "class": "VoltageSource",
                "id": "source",
                "inputs": {"bus": "b0"},
                "params": {
                    "E": [1000 / math.sqrt(2)] * 3,
                    "phi": [0, -2 * math.pi / 3, 2 * math.pi / 3],
                    "omega": 2 * math.pi * 60,
                    "Rs": diag(100),
                },
                "mon": ["ia", "ib", "ic"],
            }
        )
        events = []
    else:
        signals = [{"id": "command", "value": 0.0}]
        dev.append(
            {
                "class": "DependentVoltageSource",
                "id": "source",
                "params": {"N": 3},
                "inputs": {
                    "bus": "b0",
                    "ea": "command",
                    "eb": "command",
                    "ec": "command",
                },
                "submodels": {"Y": fit(0.01)},
                "mon": ["ia", "ib", "ic"],
            }
        )
        events = [
            {
                "time": 0.001,
                "type": "signal_step",
                "signal_id": "command",
                "value": 1000.0,
            }
        ]
    dev.append(
        {
            "class": "LoadZ",
            "id": "load",
            "inputs": {"bus": "b1"},
            "params": {"R": diag(600)},
            "mon": ["ia", "ib", "ic"],
        }
    )
    state = {}
    if kind.startswith("distributed"):
        dev.append(
            distributed(
                "line",
                "b0",
                "b1",
                fit(1 / zc),
                {"K": 3, "modes": [{"tau": tau, "H": fit(1)}]},
                kind.endswith("capped"),
            )
        )
        state["history"] = {"line": history()}
    else:
        n = int(kind.split("_")[1])
        rlgc = {"R": diag(0), "L": diag(zc * tau), "G": diag(0), "C": diag(tau / zc)}
        for i in range(1, n):
            dev.append(bus(f"p{i}"))
        nodes = ["b0"] + [f"p{i}" for i in range(1, n)] + ["b1"]
        for i in range(n):
            dev.append(pi_line(f"line{i}", nodes[i], nodes[i + 1], 1 / n, rlgc))
    case = {
        "header": {
            "case_name": "Lossless line comparison",
            "case_description": "100 ohm source, 600 ohm termination",
            "case_comments": "Analytic Bergeron and pi sections",
        },
        "signals": signals,
        "devices": dev,
    }
    meta = {
        "family": "bergeron",
        "kind": kind,
        "tau": tau,
        "sine": sine,
        "short": short,
    }
    return case, state, events, meta


def frequency_line(kind, fits):
    case, state, events, meta = bergeron(kind)
    length = 60000
    yc = json.loads((fits / "yc.json").read_text())
    h = json.loads((fits / f"h{length}.json").read_text())
    rlgc = json.loads((fits / "summary.json").read_text())["pi_per_m"]
    for device in case["devices"]:
        if device["class"] == "LineDistributed":
            device["submodels"] = {
                "Yc": yc,
                "H": {**h, "limit_step": kind.endswith("capped")},
            }
        elif device["class"] == "LineLumped":
            n = int(kind.split("_")[1])
            device["params"] = {
                "N": 3,
                "K": 3,
                "conductors": [1, 2, 3],
                "dx": length / n,
                **{key + "p": value for key, value in rlgc.items()},
            }
    case["header"]["case_name"] = "Frequency-dependent zero-sequence line energization"
    meta.update(family="frequency", length_m=length, tau=length / 299792458)
    return case, state, events, meta


def network(count, kind, hybrid, fits):
    source = json.loads((ROOT / "cases/EMT/IBR/TenBus.case.json").read_text())
    template = lambda kind: copy.deepcopy(
        next(d for d in source["devices"] if d["class"] == kind)
    )
    main = count - 2
    dev = [bus(f"b{i}") for i in range(count)]
    signals = []
    state = {"buses": {}, "devices": {}, "history": {}}
    vpeak = 13800 * math.sqrt(2 / 3)
    initial_voltage = dict(
        zip(
            ["va", "vb", "vc"],
            [0.0, -vpeak * math.sqrt(3) / 2, vpeak * math.sqrt(3) / 2],
        )
    )
    for i in range(main):
        state["buses"][f"b{i}"] = initial_voltage
    generator_nodes = [0, main // 2] if count <= 10 else [0, 6, 12]
    for i in generator_nodes:
        machine, gov = template("Machine"), template("Tgov1")
        machine.update(
            id=f"g{i}",
            inputs={"bus": f"b{i}", "pm": f"pm{i}"},
            outputs={"speed": f"w{i}"},
        )
        machine["mon"] = ["omega", "theta", "ia", "ib", "ic", "p", "q"]
        gov.update(id=f"gov{i}", inputs={"speed": f"w{i}"}, outputs={"pmech": f"pm{i}"})
        dev += [machine, gov]
        signals += [{"id": f"w{i}"}, {"id": f"pm{i}"}]
        # The machine initializer reconciles field and mechanical inputs.
        power = main * 250000 / len(generator_nodes)
        peak_current = power / (3 * vpeak / 2)
        state["devices"][f"g{i}"] = dict(
            zip(
                ["ia", "ib", "ic"],
                [0, -peak_current * math.sqrt(3) / 2, peak_current * math.sqrt(3) / 2],
            )
        )
    for i in range(main):
        dev.append(
            {
                "class": "LoadZ",
                "id": f"load{i}",
                "params": {"R": diag(13800**2 / 250000)},
                "inputs": {"bus": f"b{i}"},
            }
        )
    edges = [(i, (i + 1) % main) for i in range(main)]
    edges += [(i, (i + main // 2) % main) for i in range(main // 2) if i % 2 == 0]
    yc = json.loads((fits / "yc.json").read_text())
    coefficients = json.loads((fits / "summary.json").read_text())
    lengths = [5000, 20000, 40000, 60000]
    edge_data = []
    for k, (a, b) in enumerate(edges):
        length = lengths[k % len(lengths)]
        name = f"line{k}"
        h = json.loads((fits / f"h{length}.json").read_text())
        if kind.startswith("distributed"):
            dev.append(
                distributed(name, f"b{a}", f"b{b}", yc, h, kind.endswith("capped"))
            )
            state["history"][name] = history()
        else:
            dev.append(
                pi_line(name, f"b{a}", f"b{b}", length, coefficients["pi_per_m"])
            )
        edge_data.append(
            {
                "from": a,
                "to": b,
                "length_m": length,
                "delays_s": [m["tau"] for m in h["modes"]],
            }
        )
    # Separate switched shunts keep the advertised total bus count exact.
    for name, target, resistance in [
        ("fault", main, 10.0),
        ("load_step", main + 1, 13800**2 / 500000),
    ]:
        dev.append(
            {
                "class": "LoadZ",
                "id": name + "_load",
                "inputs": {"bus": f"b{target}"},
                "params": {"R": diag(resistance)},
            }
        )
        dev.append(
            {
                "class": "Switch",
                "id": name,
                "inputs": {"bus1": f"b{main-1}", "bus2": f"b{target}"},
                "params": {"open": True},
                "mon": ["open", "i12a"],
            }
        )
    converter_nodes = [main - 2] if count <= 10 else [3, 9, 15]
    if hybrid:
        for i in converter_nodes:
            signals += [{"id": f"dc{i}", "value": 28918.846170570516}]
            signals += [
                {"id": f"{prefix}{phase}{i}"}
                for prefix in ["s", "e", "if"]
                for phase in "abc"
            ]
            pwm, converter, filt = (
                template("PWM"),
                template("Converter"),
                template("DependentVoltageSource"),
            )
            pwm.update(id=f"pwm{i}", outputs={"s": [f"s{p}{i}" for p in "abc"]})
            converter.update(
                id=f"converter{i}",
                inputs={
                    "s": [f"s{p}{i}" for p in "abc"],
                    "vdc": f"dc{i}",
                    "i": [f"if{p}{i}" for p in "abc"],
                },
                outputs={"e": [f"e{p}{i}" for p in "abc"]},
            )
            filt.update(
                id=f"filter{i}",
                inputs={"bus": f"b{i}", **{f"e{p}": f"e{p}{i}" for p in "abc"}},
                outputs={f"i{p}": f"if{p}{i}" for p in "abc"},
            )
            dev += [pwm, converter, filt]
    case = {
        "header": {
            "case_name": f"{count}-bus frequency-dependent EMT network",
            "case_description": "Governed machines and switched shunts"
            + (" with open-loop PWM converters" if hybrid else ""),
            "case_comments": "Synthetic network; cyclically transposed GridWorkbench 345 kV geometry used at 13.8 kV for an EMT demonstration",
        },
        "signals": signals,
        "devices": dev,
    }
    meta = {
        "family": "network",
        "buses": count,
        "kind": kind,
        "hybrid": hybrid,
        "edges": edge_data,
        "generators": generator_nodes,
        "converters": converter_nodes if hybrid else [],
    }
    return case, state, meta


def run_case(args, name, case, state, events, meta, horizon, spacing):
    folder = args.results / name
    folder.mkdir(parents=True, exist_ok=True)
    meta = {
        **meta,
        "events": events,
        "horizon_s": horizon,
        "monitor_interval_s": spacing,
    }
    write(folder / "case.json", case)
    write(folder / "state.json", state)
    write(folder / "model.json", meta)
    records = []
    for trial in range(args.trials):
        if (args.results / "STOP").exists():
            print(
                "Stopping at a completed-trial boundary (results/STOP exists).",
                flush=True,
            )
            return
        out = folder / f"trial{trial+1}"
        out.mkdir(exist_ok=True)
        solver = {
            "system_model_file": "../case.json",
            "state_file": "../state.json",
            "tmax": horizon,
            "dt_monitor": spacing,
            "rel_tol": args.rtol,
            "abs_tol": args.atol,
            "max_steps": 1000000,
            "mu": 50000 if meta.get("hybrid") else 240,
            "events": events,
            "output_file": str((out / "response.csv").resolve()),
            "step_output_file": str((out / "steps.csv").resolve()),
        }
        write(out / "solver.json", solver)
        command = [
            "/usr/bin/time",
            "-f",
            "%e,%U,%M",
            "-o",
            str(out / "resources.txt"),
            str(args.exe),
            str(out / "solver.json"),
        ]
        start = time.perf_counter()
        with (out / "run.log").open("w") as log:
            proc = subprocess.run(
                command, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT
            )
        wall = time.perf_counter() - start
        text = (out / "run.log").read_text()
        record = {
            "returncode": proc.returncode,
            "wall_s": wall,
            "trial": trial + 1,
            "rtol": args.rtol,
            "atol": args.atol,
        }
        for key, pattern in [
            ("cpu_s", r"Complete in ([\d.e+-]+)"),
            ("variables", r"DAE variables: (\d+)"),
            ("nnz", r"Jacobian nonzeros: (\d+)"),
        ]:
            match = re.search(pattern, text)
            if match:
                record[key] = float(match[1])
        for key, value in re.findall(
            r"(steps|residual_evals|linear_setups|error_test_fails|nonlinear_iters|nonlinear_convergence_fails)=(\d+)",
            text,
        ):
            record[key] = int(value)
        resources = (out / "resources.txt").read_text().splitlines()[-1].split(",")
        if len(resources) == 3:
            record["peak_rss_kb"] = int(resources[2])
        record["executable_sha256"] = hashlib.sha256(args.exe.read_bytes()).hexdigest()
        records.append(record)
        write(folder / "statistics.json", records)
        print(name, record, flush=True)
        if proc.returncode:
            print(text[-2000:], flush=True)
            break


def main(args):
    args.results = args.results.resolve()
    if (args.results / "STOP").exists():
        print("Study batch stopped: remove results/STOP to resume.", flush=True)
        return
    args.exe = args.exe.resolve()
    args.fits = args.fits.resolve()
    if args.models is None:
        args.models = (
            ["distributed", "pi_1"]
            if args.family == "network"
            else ["distributed", "distributed_capped", "pi_1", "pi_10"]
        )
    if args.family == "network" and "pi_10" in args.models:
        raise ValueError(
            "Network studies use one pi section per edge to preserve the bus count"
        )
    if args.family == "frequency" and (args.sine or args.short):
        raise ValueError(
            "The frequency-dependent line study uses a 60 km line and a step source"
        )
    if args.trials <= 0 or args.rtol <= 0 or args.atol <= 0 or args.horizon <= 0:
        raise ValueError("Trials, tolerances, and horizon must be positive")
    if args.family in ["bergeron", "frequency"]:
        for kind in args.models:
            case, state, events, meta = (
                frequency_line(kind, args.fits)
                if args.family == "frequency"
                else bergeron(kind, args.sine, args.short)
            )
            name = f'{args.family}_{"sine" if args.sine else "step"}_{"short_" if args.short else ""}{kind}'
            run_case(
                args,
                name,
                case,
                state,
                events,
                meta,
                0.12 if args.sine else 0.008,
                10e-6 if args.sine else 2e-6,
            )
    else:
        for count in args.buses:
            for kind in args.models:
                case, state, meta = network(count, kind, args.hybrid, args.fits)
                meta["event"] = args.event
                events = []
                if args.event == "load":
                    events = [
                        {
                            "time": args.horizon * 0.5,
                            "type": "switch",
                            "element_id": "load_step",
                            "open": False,
                        },
                        {
                            "time": args.horizon * 0.75,
                            "type": "switch",
                            "element_id": "load_step",
                            "open": True,
                        },
                    ]
                elif args.event in ["fault", "fault_a"]:
                    events = [
                        {
                            "time": args.horizon * 0.5,
                            "type": "switch",
                            "element_id": "fault",
                            "open": False,
                        },
                        {
                            "time": args.horizon * 0.5 + 0.02,
                            "type": "switch",
                            "element_id": "fault",
                            "open": True,
                        },
                    ]
                if args.event == "fault_a":
                    fault = next(d for d in case["devices"] if d["id"] == "fault_load")
                    fault["params"]["R"] = [
                        [10.0, 0.0, 0.0],
                        [0.0, 1e9, 0.0],
                        [0.0, 0.0, 1e9],
                    ]
                    for event in events:
                        event["time"] += 1 / 240
                name = f'network{count}_{"hybrid" if args.hybrid else "machine"}_{args.event}_{kind}'
                run_case(args, name, case, state, events, meta, args.horizon, 20e-6)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--exe", type=Path, default=ROOT / "build/application/EMT/EMTDynamicSimulation"
    )
    p.add_argument(
        "--results", type=Path, default=ROOT / "local/emt-distributed/studies"
    )
    p.add_argument(
        "--fits", type=Path, default=ROOT / "local/emt-distributed/fits-final"
    )
    p.add_argument(
        "--family", choices=["bergeron", "frequency", "network"], default="bergeron"
    )
    p.add_argument(
        "--models",
        nargs="+",
        choices=["distributed", "distributed_capped", "pi_1", "pi_10"],
    )
    p.add_argument(
        "--buses", type=int, nargs="+", choices=[7, 8, 9, 10, 20], default=[8, 20]
    )
    p.add_argument("--sine", action="store_true")
    p.add_argument("--short", action="store_true")
    p.add_argument("--hybrid", action="store_true")
    p.add_argument(
        "--event", choices=["baseline", "load", "fault", "fault_a"], default="load"
    )
    p.add_argument("--horizon", type=float, default=0.4)
    p.add_argument("--rtol", type=float, default=1e-6)
    p.add_argument("--atol", type=float, default=1e-7)
    p.add_argument("--trials", type=int, default=3)
    main(p.parse_args())
