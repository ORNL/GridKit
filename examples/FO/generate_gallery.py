#!/usr/bin/env python3
"""Generate the 100-study WECC240 forced-oscillation response gallery."""

import argparse
import copy
import csv
import json
import math
import os
import re
import subprocess
import tempfile
import time
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "gridkit-fo-matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


SURFACE = "#fcfcfb"
INK = "#171717"
MUTED = "#65625e"
GRID = "#dedbd5"
FLEET = "#777777"
LOCAL = "#d1495b"
SYNCHRONOUS = "#3567a8"
RENEWABLE = "#25886f"

DT_MONITOR = 0.05
TMAX = 10.0
TON = 1.0

SUFFIXES = {
    "GastPti": "_gast_pti",
    "Genrou": "_genrou",
    "Hygov": "_hygov",
    "Ieeest": "_ieeest",
    "SexsPti": "_sexs_pti",
    "Tgov1": "_tgov1",
    "Regca": "_regca",
    "Reecb": "_reecb",
    "Repca": "_repca",
}

BASE_AMPLITUDE = {
    ("GastPti", "speed"): 5.0e-4,
    ("Genrou", "efd"): 2.0e-2,
    ("Genrou", "pmech"): 4.0e-2,
    ("Hygov", "speed"): 4.0e-4,
    ("Ieeest", "input"): 1.0e-3,
    ("SexsPti", "vs"): 7.5e-4,
    ("Tgov1", "speed"): 5.0e-4,
    ("Regca", "ipcmd"): 1.5e-2,
    ("Regca", "iqcmd"): 1.5e-2,
    ("Reecb", "pe"): 1.0e-2,
    ("Reecb", "pref"): 1.0e-2,
    ("Reecb", "qext"): 1.0e-2,
    ("Reecb", "qgen"): 1.0e-2,
    ("Repca", "ii"): 1.0e-2,
    ("Repca", "ir"): 1.0e-2,
    ("Repca", "p"): 1.0e-2,
    ("Repca", "q"): 1.0e-2,
}

PROFILES = (
    {"name": "stationary sine", "waveform": 0, "f": 0.35, "Kf": 0.0, "Kd": 0.0, "Phi": 0.0},
    {"name": "chirped sine", "waveform": 0, "f": 0.20, "Kf": 0.12, "Kd": 0.0, "Phi": 0.2},
    {"name": "decaying sine", "waveform": 0, "f": 0.80, "Kf": 0.0, "Kd": 0.18, "Phi": 0.4},
    {"name": "chirped decaying sine", "waveform": 0, "f": 0.25, "Kf": 0.09, "Kd": 0.08, "Phi": 0.6},
    {"name": "stationary square", "waveform": 1, "f": 0.35, "Kf": 0.0, "Kd": 0.0, "Phi": 0.0},
    {"name": "chirped square", "waveform": 1, "f": 0.20, "Kf": 0.10, "Kd": 0.0, "Phi": 0.2},
    {"name": "stationary triangle", "waveform": 2, "f": 0.55, "Kf": 0.0, "Kd": 0.0, "Phi": 0.0},
    {"name": "decaying triangle", "waveform": 2, "f": 0.90, "Kf": 0.0, "Kd": 0.12, "Phi": 0.3},
    {"name": "stationary sawtooth", "waveform": 3, "f": 0.40, "Kf": 0.0, "Kd": 0.0, "Phi": 0.0},
    {"name": "chirped sawtooth", "waveform": 3, "f": 0.18, "Kf": 0.12, "Kd": 0.0, "Phi": 0.1},
)

SYNCHRONOUS_TARGETS = (
    "GastPti:1131_G_gast_pti.speed",
    "GastPti:2233_EG_gast_pti.speed",
    "GastPti:3133_NG_gast_pti.speed",
    "GastPti:3234_NG_gast_pti.speed",
    "GastPti:3836_DG_gast_pti.speed",
    "GastPti:4039_G_gast_pti.speed",
    "GastPti:5032_G_gast_pti.speed",
    "GastPti:6433_G_gast_pti.speed",
    "Genrou:1131_C_genrou.efd",
    "Genrou:1331_G_genrou.efd",
    "Genrou:2438_RG_genrou.efd",
    "Genrou:2637_H_genrou.efd",
    "Genrou:3135_MG_genrou.efd",
    "Genrou:3333_NG_genrou.efd",
    "Genrou:3631_NB_genrou.efd",
    "Genrou:3831_NN_genrou.efd",
    "Genrou:3933_CG_genrou.efd",
    "Genrou:4031_G_genrou.efd",
    "Genrou:4132_G_genrou.efd",
    "Genrou:4231_C_genrou.efd",
    "Genrou:6132_B_genrou.efd",
    "Genrou:6231_C_genrou.efd",
    "Genrou:6335_H_genrou.efd",
    "Genrou:8033_H_genrou.efd",
    "Genrou:1032_C_genrou.pmech",
    "Genrou:1232_C_genrou.pmech",
    "Genrou:2130_E_genrou.pmech",
    "Genrou:2634_C_genrou.pmech",
    "Genrou:2638_H_genrou.pmech",
    "Genrou:3333_CG_genrou.pmech",
    "Genrou:3433_NG_genrou.pmech",
    "Genrou:3731_NH_genrou.pmech",
    "Genrou:3835_NG_genrou.pmech",
    "Genrou:3933_NH_genrou.pmech",
    "Genrou:4035_C_genrou.pmech",
    "Genrou:4132_N_genrou.pmech",
    "Genrou:4232_G_genrou.pmech",
    "Genrou:6132_G_genrou.pmech",
    "Genrou:6333_C_genrou.pmech",
    "Genrou:7032_C_genrou.pmech",
    "Hygov:2130_H_hygov.speed",
    "Hygov:3531_NH_hygov.speed",
    "Hygov:4039_H_hygov.speed",
    "Hygov:5031_H_hygov.speed",
    "Hygov:7032_H_hygov.speed",
    "Ieeest:1333_G_ieeest.input",
    "Ieeest:2438_WG_ieeest.input",
    "Ieeest:3931_NH_ieeest.input",
    "Ieeest:6235_H_ieeest.input",
    "Ieeest:7031_G_ieeest.input",
    "SexsPti:2030_G_sexs_pti.vs",
    "SexsPti:2630_G_sexs_pti.vs",
    "SexsPti:4131_H_sexs_pti.vs",
    "SexsPti:6533_C_sexs_pti.vs",
    "SexsPti:8034_G_sexs_pti.vs",
    "Tgov1:1034_C_tgov1.speed",
    "Tgov1:1431_N_tgov1.speed",
    "Tgov1:3531_CE_tgov1.speed",
    "Tgov1:3931_NB_tgov1.speed",
    "Tgov1:4131_B_tgov1.speed",
    "Tgov1:5032_R_tgov1.speed",
    "Tgov1:6335_C_tgov1.speed",
    "Tgov1:7031_P_tgov1.speed",
)

RENEWABLE_TARGETS = (
    "Regca:1032_S_regca.ipcmd",
    "Regca:1034_W_regca.iqcmd",
    "Reecb:1333_S_reecb.pref",
    "Reecb:1431_S_reecb.qext",
    "Reecb:2130_S_reecb.pe",
    "Reecb:2332_S_reecb.qgen",
    "Repca:2431_S_repca.ir",
    "Repca:2434_S_repca.ii",
    "Repca:2438_S_repca.p",
    "Repca:2438_SW_repca.q",
    "Regca:2439_S_regca.ipcmd",
    "Regca:2533_S_regca.iqcmd",
    "Reecb:2631_S_reecb.pref",
    "Reecb:3234_NW_reecb.qext",
    "Reecb:3433_S_reecb.pe",
    "Reecb:3835_S_reecb.qgen",
    "Repca:3932_S_repca.ir",
    "Repca:3933_NW_repca.ii",
    "Repca:3933_S_repca.p",
    "Repca:4031_S_repca.q",
    "Regca:4031_W_regca.ipcmd",
    "Regca:4035_W_regca.iqcmd",
    "Reecb:4039_W_reecb.pref",
    "Reecb:4131_W_reecb.qext",
    "Reecb:4132_W_reecb.pe",
    "Reecb:4232_W_reecb.qgen",
    "Repca:5032_S_repca.ir",
    "Repca:5032_W_repca.ii",
    "Repca:6132_S_repca.p",
    "Repca:6132_W_repca.q",
    "Regca:6235_W_regca.ipcmd",
    "Regca:6333_W_regca.iqcmd",
    "Reecb:6433_W_reecb.pref",
    "Reecb:6533_S_reecb.qext",
    "Reecb:6533_W_reecb.pe",
    "Reecb:7031_W_reecb.qgen",
    "Repca:7032_S_repca.ir",
)

PLOT_NAMES = (
    "01-voltage-magnitude.png",
    "02-bus-angle.png",
    "03-active-power.png",
    "04-reactive-power.png",
)


def parse_target(target):
    component_class, rest = target.split(":", 1)
    component_id, port = rest.rsplit(".", 1)
    suffix = SUFFIXES[component_class]
    if not component_id.endswith(suffix):
        raise ValueError(f"Unexpected component id {component_id}")
    unit = component_id[: -len(suffix)]
    bus = int(unit.split("_", 1)[0])
    return component_class, component_id, port, unit, bus


def study_definitions(case):
    devices = {(d["class"], d["id"]): d for d in case["devices"]}
    targets = tuple(("synchronous", t) for t in SYNCHRONOUS_TARGETS)
    targets += tuple(("renewable", t) for t in RENEWABLE_TARGETS)
    if len(targets) != 100:
        raise ValueError(f"Expected 100 targets, found {len(targets)}")

    studies = []
    units = set()
    for number, (domain, target) in enumerate(targets, start=1):
        component_class, component_id, port, unit, bus = parse_target(target)
        device = devices.get((component_class, component_id))
        if device is None or port not in device.get("ports", {}):
            raise ValueError(f"Target is not connected in WECC240: {target}")
        if unit in units:
            raise ValueError(f"Duplicate physical unit target: {unit}")
        units.add(unit)

        profile = dict(PROFILES[(number - 1) % len(PROFILES)])
        profile["Phi"] += 0.13 * (((number - 1) // len(PROFILES)) % 5)
        amplitude_scale = 0.85 + 0.05 * ((number * 3) % 7)
        amplitude = BASE_AMPLITUDE[(component_class, port)] * amplitude_scale
        slug = re.sub(r"[^a-z0-9]+", "_", f"{component_class}_{unit}_{port}".lower()).strip("_")

        studies.append(
            {
                "number": number,
                "id": f"fo_{number:03d}",
                "domain": domain,
                "target": target,
                "component_class": component_class,
                "component_id": component_id,
                "port": port,
                "unit": unit,
                "bus": bus,
                "signal_id": device["ports"][port],
                "profile": profile["name"],
                "params": {
                    "A": amplitude,
                    "f": profile["f"],
                    "Kf": profile["Kf"],
                    "Phi": profile["Phi"],
                    "Ton": TON,
                    "Toff": -1.0,
                    "Tr": 0.15,
                    "Tf": 0.0,
                    "Kd": profile["Kd"],
                    "waveform": profile["waveform"],
                },
                "relative_directory": f"{domain}/{number:03d}_{slug}",
            }
        )
    return studies


def instrument_case(case, csv_path):
    instrumented = copy.deepcopy(case)
    instrumented["monitors"] = [{"file_name": str(csv_path), "format": "csv"}]
    for bus in instrumented["buses"]:
        bus["name"] = f'{bus["number"]}_{bus["name"].strip()}'
        bus["mon"] = ["Vm", "Va"]
    for device in instrumented["devices"]:
        device.pop("mon", None)
        if device["class"] in {"Genrou", "Regca"}:
            device["mon"] = ["p", "q"]
    return instrumented


def load_output(csv_path):
    with csv_path.open(newline="") as stream:
        header = next(csv.reader(stream))
    values = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    if values.ndim == 1:
        values = values.reshape(1, -1)
    if values.shape[1] != len(header):
        raise ValueError("Monitor header and data column counts differ")
    return header, values


def column_indices(header, prefix, suffix):
    return [i for i, name in enumerate(header) if name.startswith(prefix) and name.endswith(suffix)]


def prepare_response(header, values, study):
    t = values[:, 0]
    before = np.flatnonzero(t < TON)
    if not len(before):
        raise ValueError("No pre-activation monitor sample")
    baseline = before[-1]

    vm_indices = column_indices(header, "Bus_", "_Vm")
    va_indices = column_indices(header, "Bus_", "_Va")
    gp_indices = column_indices(header, "Genrou_", "_p")
    gq_indices = column_indices(header, "Genrou_", "_q")
    rp_indices = column_indices(header, "Regca_", "_p")
    rq_indices = column_indices(header, "Regca_", "_q")

    expected = (243, 243, 103, 103, 37, 37)
    actual = tuple(map(len, (vm_indices, va_indices, gp_indices, gq_indices, rp_indices, rq_indices)))
    if actual != expected:
        raise ValueError(f"Unexpected response-channel counts: {actual}")

    source_name = f'ForcedOscillation_{study["id"]}_output'
    try:
        source_index = header.index(source_name)
    except ValueError as error:
        raise ValueError(f"Missing source monitor {source_name}") from error

    vm = values[:, vm_indices]
    va = np.unwrap(values[:, va_indices], axis=0)
    gp = values[:, gp_indices]
    gq = values[:, gq_indices]
    rp = values[:, rp_indices]
    rq = values[:, rq_indices]

    local_bus_prefix = f'Bus_{study["bus"]}_'
    local_vm = next(i for i, column in enumerate(vm_indices) if header[column].startswith(local_bus_prefix))
    local_va = next(i for i, column in enumerate(va_indices) if header[column].startswith(local_bus_prefix))

    if study["domain"] == "synchronous":
        local_prefix = f'Genrou_{study["unit"]}_genrou_'
        local_p = ("synchronous", next(i for i, column in enumerate(gp_indices) if header[column].startswith(local_prefix)))
        local_q = ("synchronous", next(i for i, column in enumerate(gq_indices) if header[column].startswith(local_prefix)))
    else:
        local_prefix = f'Regca_{study["unit"]}_regca_'
        local_p = ("renewable", next(i for i, column in enumerate(rp_indices) if header[column].startswith(local_prefix)))
        local_q = ("renewable", next(i for i, column in enumerate(rq_indices) if header[column].startswith(local_prefix)))

    return {
        "t": t,
        "source": values[:, source_index],
        "vm": vm - vm[baseline],
        "va": np.rad2deg(va - va[baseline]),
        "gp": gp - gp[baseline],
        "gq": gq - gq[baseline],
        "rp": rp - rp[baseline],
        "rq": rq - rq[baseline],
        "local_vm": local_vm,
        "local_va": local_va,
        "local_p": local_p,
        "local_q": local_q,
        "baseline_time": t[baseline],
    }


def style_axis(axis):
    axis.set_facecolor(SURFACE)
    axis.grid(True, color=GRID, linewidth=0.6)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    axis.spines[["left", "bottom"]].set_color(GRID)
    axis.tick_params(colors=MUTED, labelsize=8)


def finish_limits(axis, arrays):
    limit = max(float(np.nanmax(np.abs(array))) for array in arrays)
    if not math.isfinite(limit) or limit < 1.0e-12:
        limit = 1.0e-12
    axis.set_ylim(-1.08 * limit, 1.08 * limit)
    axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    axis.yaxis.get_offset_text().set(color=MUTED, fontsize=8)
    return limit


def figure_axes(study, response):
    fig = plt.figure(figsize=(9.5, 5.8), dpi=120, facecolor=SURFACE)
    grid = fig.add_gridspec(2, 1, height_ratios=(1.0, 3.4), hspace=0.16)
    source_axis = fig.add_subplot(grid[0])
    response_axis = fig.add_subplot(grid[1], sharex=source_axis)

    source_axis.axvspan(TON, TMAX, color="#ece8df", alpha=0.65, linewidth=0)
    source_axis.axvline(TON, color=MUTED, linewidth=0.8, linestyle="--")
    source_axis.plot(response["t"], response["source"], color=INK, linewidth=1.2)
    source_axis.axhline(0.0, color=GRID, linewidth=0.7)
    source_axis.set_ylabel("forcing", color=MUTED, fontsize=8)
    source_axis.tick_params(labelbottom=False)
    style_axis(source_axis)

    parameters = study["params"]
    source_axis.set_title(
        f'{study["profile"]}  ·  A={parameters["A"]:.4g}  ·  '
        f'f₀={parameters["f"]:.3g} Hz  ·  Kf={parameters["Kf"]:.3g} Hz/s  ·  '
        f'Kd={parameters["Kd"]:.3g} 1/s',
        loc="left",
        color=MUTED,
        fontsize=8.5,
        pad=5,
    )

    fig.suptitle(
        f'{study["number"]:03d} · {study["domain"]} · {study["target"]}',
        x=0.075,
        y=0.988,
        ha="left",
        color=INK,
        fontsize=12,
        fontweight="semibold",
    )
    return fig, response_axis


def save_bus_plot(path, study, response, field, local_field, ylabel):
    fig, axis = figure_axes(study, response)
    values = response[field]
    axis.plot(response["t"], values, color=FLEET, linewidth=0.45, alpha=0.12)
    axis.plot(response["t"], values[:, response[local_field]], color=LOCAL, linewidth=1.5)
    axis.axvline(TON, color=MUTED, linewidth=0.8, linestyle="--")
    axis.axhline(0.0, color=GRID, linewidth=0.7)
    axis.set_ylabel(ylabel, color=MUTED, fontsize=9)
    axis.set_xlabel("time (s)", color=MUTED, fontsize=9)
    axis.set_xlim(0.0, TMAX)
    style_axis(axis)
    finish_limits(axis, (values,))
    axis.legend(
        handles=(
            Line2D([], [], color=FLEET, linewidth=1.0, label="all 243 buses"),
            Line2D([], [], color=LOCAL, linewidth=1.5, label=f'bus {study["bus"]}'),
        ),
        loc="upper right",
        frameon=False,
        fontsize=8,
    )
    fig.text(
        0.075,
        0.012,
        f'response relative to t={response["baseline_time"]:.2f} s',
        color=MUTED,
        fontsize=7.5,
    )
    fig.savefig(path, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def save_power_plot(path, study, response, sync_field, renewable_field, local_field, ylabel):
    fig, axis = figure_axes(study, response)
    sync = response[sync_field]
    renewable = response[renewable_field]
    axis.plot(response["t"], sync, color=SYNCHRONOUS, linewidth=0.45, alpha=0.11)
    axis.plot(response["t"], renewable, color=RENEWABLE, linewidth=0.5, alpha=0.15)
    domain, local_index = response[local_field]
    local = sync[:, local_index] if domain == "synchronous" else renewable[:, local_index]
    axis.plot(response["t"], local, color=LOCAL, linewidth=1.5)
    axis.axvline(TON, color=MUTED, linewidth=0.8, linestyle="--")
    axis.axhline(0.0, color=GRID, linewidth=0.7)
    axis.set_ylabel(ylabel, color=MUTED, fontsize=9)
    axis.set_xlabel("time (s)", color=MUTED, fontsize=9)
    axis.set_xlim(0.0, TMAX)
    style_axis(axis)
    finish_limits(axis, (sync, renewable))
    axis.legend(
        handles=(
            Line2D([], [], color=SYNCHRONOUS, linewidth=1.0, label="103 synchronous machines"),
            Line2D([], [], color=RENEWABLE, linewidth=1.0, label="37 renewable converters"),
            Line2D([], [], color=LOCAL, linewidth=1.5, label="associated resource"),
        ),
        loc="upper right",
        frameon=False,
        fontsize=8,
    )
    fig.text(
        0.075,
        0.012,
        f'response relative to t={response["baseline_time"]:.2f} s',
        color=MUTED,
        fontsize=7.5,
    )
    fig.savefig(path, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def save_plots(output_directory, study, response):
    output_directory.mkdir(parents=True, exist_ok=True)
    save_bus_plot(
        output_directory / PLOT_NAMES[0],
        study,
        response,
        "vm",
        "local_vm",
        "Δ|V| (pu)",
    )
    save_bus_plot(
        output_directory / PLOT_NAMES[1],
        study,
        response,
        "va",
        "local_va",
        "Δ∠V (deg)",
    )
    save_power_plot(
        output_directory / PLOT_NAMES[2],
        study,
        response,
        "gp",
        "rp",
        "local_p",
        "ΔP (pu, system base)",
    )
    save_power_plot(
        output_directory / PLOT_NAMES[3],
        study,
        response,
        "gq",
        "rq",
        "local_q",
        "ΔQ (pu, system base)",
    )


def run_study(binary, case_path, csv_path, solver_path, study):
    solver = {
        "system_model_file": str(case_path),
        "dt_monitor": DT_MONITOR,
        "tmax": TMAX,
        "max_steps": 10000,
        "events": [],
        "forced_oscillations": [
            {
                "id": study["id"],
                "target": study["target"],
                "mode": "add",
                "params": study["params"],
                "mon": ["output"],
            }
        ],
    }
    solver_path.write_text(json.dumps(solver, indent=2) + "\n")
    csv_path.unlink(missing_ok=True)

    started = time.perf_counter()
    result = subprocess.run(
        (str(binary), str(solver_path)),
        cwd=solver_path.parent,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    elapsed = time.perf_counter() - started
    if result.returncode != 0 or not csv_path.exists():
        diagnostics = (result.stdout + "\n" + result.stderr).strip().splitlines()
        raise RuntimeError(
            f'Study {study["number"]:03d} failed (exit {result.returncode}):\n'
            + "\n".join(diagnostics[-30:])
        )
    return elapsed


def extrema(response):
    return {
        "peak_input": float(np.max(np.abs(response["source"]))),
        "max_abs_dvm": float(np.max(np.abs(response["vm"]))),
        "max_abs_dva_deg": float(np.max(np.abs(response["va"]))),
        "max_abs_dp": float(max(np.max(np.abs(response["gp"])), np.max(np.abs(response["rp"])))),
        "max_abs_dq": float(max(np.max(np.abs(response["gq"])), np.max(np.abs(response["rq"])))),
    }


def write_manifest(path, rows):
    fields = (
        "study",
        "domain",
        "target",
        "signal_id",
        "bus",
        "profile",
        "waveform",
        "A",
        "f",
        "Kf",
        "Kd",
        "Phi",
        "Ton",
        "Toff",
        "Tr",
        "samples",
        "simulation_seconds",
        "peak_input",
        "max_abs_dvm",
        "max_abs_dva_deg",
        "max_abs_dp",
        "max_abs_dq",
        "directory",
    )
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    repository = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary",
        type=Path,
        default=repository / "build/application/PhasorDynamics/DynamicSimulation",
    )
    parser.add_argument("--output-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--limit", type=int, help="generate only the first N studies")
    args = parser.parse_args()

    binary = args.binary.resolve()
    if not binary.is_file():
        raise SystemExit(f"DynamicSimulation not found: {binary}")

    case_path = repository / "cases/PhasorDynamics/WECC240/WECC240.case.json"
    base_case = json.loads(case_path.read_text())
    studies = study_definitions(base_case)
    if args.limit is not None:
        if args.limit < 1 or args.limit > len(studies):
            raise SystemExit("--limit must be between 1 and 100")
        studies = studies[: args.limit]

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    rows = []

    with tempfile.TemporaryDirectory(prefix="gridkit-fo-gallery-") as temporary:
        temporary = Path(temporary)
        csv_path = temporary / "response.csv"
        instrumented_path = temporary / "WECC240.fo-monitor.case.json"
        solver_path = temporary / "study.json"
        instrumented = instrument_case(base_case, csv_path)
        instrumented_path.write_text(json.dumps(instrumented))

        for study in studies:
            elapsed = run_study(binary, instrumented_path, csv_path, solver_path, study)
            header, values = load_output(csv_path)
            csv_path.unlink()
            response = prepare_response(header, values, study)
            save_plots(output_root / study["relative_directory"], study, response)
            metrics = extrema(response)

            parameters = study["params"]
            rows.append(
                {
                    "study": f'{study["number"]:03d}',
                    "domain": study["domain"],
                    "target": study["target"],
                    "signal_id": study["signal_id"],
                    "bus": study["bus"],
                    "profile": study["profile"],
                    "waveform": parameters["waveform"],
                    "A": f'{parameters["A"]:.8g}',
                    "f": f'{parameters["f"]:.8g}',
                    "Kf": f'{parameters["Kf"]:.8g}',
                    "Kd": f'{parameters["Kd"]:.8g}',
                    "Phi": f'{parameters["Phi"]:.8g}',
                    "Ton": parameters["Ton"],
                    "Toff": parameters["Toff"],
                    "Tr": parameters["Tr"],
                    "samples": len(response["t"]),
                    "simulation_seconds": f"{elapsed:.6f}",
                    "peak_input": f'{metrics["peak_input"]:.8e}',
                    "max_abs_dvm": f'{metrics["max_abs_dvm"]:.8e}',
                    "max_abs_dva_deg": f'{metrics["max_abs_dva_deg"]:.8e}',
                    "max_abs_dp": f'{metrics["max_abs_dp"]:.8e}',
                    "max_abs_dq": f'{metrics["max_abs_dq"]:.8e}',
                    "directory": study["relative_directory"],
                }
            )
            write_manifest(output_root / "manifest.csv", rows)
            print(
                f'[{study["number"]:03d}/100] {study["target"]} · '
                f'{study["profile"]} · {elapsed:.3f} s'
            )

    expected_plots = 4 * len(studies)
    actual_plots = sum(1 for path in output_root.rglob("*.png") if path.is_file())
    if actual_plots != expected_plots:
        raise RuntimeError(f"Expected {expected_plots} plots, found {actual_plots}")
    print(f"Generated {actual_plots} plots for {len(studies)} studies in {output_root}")


if __name__ == "__main__":
    main()
