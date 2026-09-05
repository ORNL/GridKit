#!/usr/bin/env python3
"""Write a measured CoupledGrid report from summary.json and metrics.json."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import statistics

import numpy as np

from pwm_analysis import attenuation

HERE = Path(__file__).resolve().parent


def number(value, digits=5):
    return "—" if value is None else f"{value:.{digits}g}"


def table(headers, rows):
    return "\n".join([
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
        *("| " + " | ".join(str(value) for value in row) + " |" for row in rows),
    ])


def impedance(model, frequency):
    s = 2j * np.pi * frequency
    result = np.array(model["D"], complex) + s * np.array(model["E"])
    for pole, residue in zip(model["poles"], model["residues"]):
        matrix = np.array(residue)
        result += (matrix[..., 0] + 1j * matrix[..., 1]) / (s - complex(*pole))
    return result


def sequence_impedance(matrix):
    a = np.exp(2j * np.pi / 3)
    transform = np.array([[1, 1, 1], [1, a, a*a], [1, a*a, a]]) / 3
    return transform @ matrix @ np.linalg.inv(transform)


def rms_reactor(window, name):
    currents = window["ibrs"][name]["phase_current_rms_a"]
    return math.sqrt(statistics.mean(value * value for value in currents))


def window_label(name, window):
    label = {"before_event_1": "Before close", "before_event_2": "Load connected",
             "final": "Final"}.get(name, name)
    return f"{label} [{window['start_s']:.3f}, {window['end_s']:.3f}) s"


def statistic(record, name):
    values = [entry[name] for entry in record.get("trial_stats", []) if name in entry]
    if not values:
        value = record.get("stats", {}).get(name)
        return "—" if value is None else f"{value:,}"
    return f"{statistics.median(values):,.0f}"


def make_report(summary, metrics, case, fit, line_model):
    records = {float(row["mu"]): row for row in summary["runs"]}
    runs = sorted(metrics["runs"], key=lambda row: row["mu"])
    if not runs or any(float(run["mu"]) not in records for run in runs):
        raise ValueError("Metrics and manifest have incompatible runs")
    for run in runs:
        record = records[float(run["mu"])]
        if not record["complete_seconds"] or record["complete_seconds"] != run["complete_seconds"]:
            raise ValueError("Metrics timing data differs from the manifest; rerun analyze.py")
    frequency = summary["frequency_hz"]
    carrier = summary["switching_frequency_hz"]
    remote_bus = summary.get("remote_bus", "bus_10")
    baseline_mu = metrics["baseline_mu"]
    baseline = next(run for run in runs if run["mu"] == baseline_mu)
    baseline_cpu = statistics.median(records[baseline_mu]["complete_seconds"])
    devices = case["devices"]
    lines = [device for device in devices if device["class"] == "LineLumped"]
    machines = sorted(device["id"] for device in devices if device["class"] == "Machine")
    reactors = sorted((device["id"] for device in devices if device["class"] == "DependentVoltageSource"),
                      key=lambda name: int(name.rsplit("_", 1)[-1]))
    physical = fit["physical_parameters"]
    source = fit["source"]
    expected = {4 * frequency, 4 * math.sqrt(frequency * carrier), 4 * carrier}
    present = {float(run["mu"]) for run in runs}
    full_sweep = len(present) == len(expected) and all(
        any(math.isclose(value, target, rel_tol=1e-12) for value in present) for target in expected)
    replicated = all(len(records[run["mu"]]["complete_seconds"]) >= 3 for run in runs)
    refined = all("tolerance_refinement" in run for run in runs)
    complete = full_sweep and replicated and refined
    text = ["# CoupledGrid: measured EMT mu comparison", "",
            ("Three mu values, at least three primary trials per value, and same-mu tolerance refinements are available."
             if complete else
             "**Partial result set.** Tables below describe only the available completed runs; missing trials or tolerance refinements are not inferred."), "",
            f"Synthetic {sum(device['class'] == 'Bus' for device in devices)}-bus, "
            f"{summary['voltage_ll_v']/1000:g} kV, {frequency:g} Hz network with {len(lines)} coupled "
            f"overhead-line pi sections, {len(machines)} synchronous machines with governors and IEEET1 exciters, "
            f"and {len(reactors)} PWM–Converter–DVS sources. Carrier frequency: {carrier:g} Hz. "
            f"The recorded duration is {number(runs[0]['final_time_s'])} s.", "",
            "[Rendered one-line](one_line.svg) · [Plot gallery](index.html) · "
            "[Raw run manifest](summary.json) · [Numerical metrics](metrics.json)", "",
            "## Runtime and solver work", ""]
    rows = []
    for run in runs:
        record = records[run["mu"]]
        cpu = record["complete_seconds"]
        wall = record.get("wall_seconds", [])
        rows.append([number(run["mu"], 8), len(cpu), number(statistics.median(cpu), 7),
                     f"{number(min(cpu), 6)}–{number(max(cpu), 6)}",
                     number(statistics.median(wall), 7) if wall else "—",
                     statistic(record, "steps"), statistic(record, "residual_evals"),
                     number(statistics.median(cpu) / baseline_cpu, 4)])
    text += [table(["Mu", "Primary trials", "CPU median (s)", "CPU range (s)",
                    "Wall median (s)", "IDA steps", "Residual evaluations", "CPU / high-mu"], rows), "",
             "Application `Complete in` measures process CPU time spent in consistent initialization, "
             "integration, event handling, and monitoring I/O. Wall time is measured separately around the "
             "whole application and also includes model construction. IDA counters include all event segments "
             "and initialization work; the table uses their median across available trials. Tight runs are "
             "excluded from the primary runtime medians. Solver work and waveform accuracy are separate measurements.", ""]
    rows = []
    for run in runs:
        record = records[run["mu"]]
        calls = record["stats"]["residual_evals"]
        rows.append([number(run["mu"], 8), statistic(record, "linear_setups"),
                     statistic(record, "error_test_fails"), statistic(record, "nonlinear_convergence_fails"),
                     number(1e6 * statistics.median(record["complete_seconds"]) / calls, 6)])
    text += [table(["Mu", "Linear setups", "Error-test failures", "Nonlinear convergence failures",
                    "Total CPU / residual evaluation (us)"], rows), "",
             "The last column amortizes total application CPU time over residual calls; it is not a profile "
             "of the residual kernel alone. Lower mu widens PWM's pulse-history window, increasing the "
             "number of pulse evaluations per call. Error tests and linear setups also contribute to "
             "runtime, so smoother waveforms do not imply a faster simulation.", ""]
    rows = []
    for run in runs:
        record = records[run["mu"]]
        rows.append([number(run["mu"], 8), number(record.get("rel_tol")), number(record.get("abs_tol")),
                     number(record.get("tight_rel_tol")), number(record.get("tight_abs_tol")),
                     number(run["monitor_dt_s"] * 1e6), record.get("dae_variables", "—"),
                     record.get("jacobian_nonzeros", "—")])
    text += [table(["Mu", "Primary rel tol", "Primary abs tol", "Tight rel tol", "Tight abs tol",
                    "Monitor spacing (us)", "DAE variables", "Initial Jacobian nnz"], rows), "",
             "[Runtime and solver-work plot](plots/runtime_comparison.svg)", "",
             "## What changing mu does", "",
             "Mu changes the applied PWM voltage and the global CommonMath smoothing in machine saturation "
             "and controller limits. It changes the model equations; it is not an integration tolerance. "
             "All runs keep the physical case and ideal DC-source voltage fixed. The highest mu is a "
             "comparison baseline, not a hard-switching limit or an accuracy oracle.", "",
             "For a smoothed pulse train, the harmonic amplitude multiplier is "
             "`A(f,mu) = x/sinh(x)`, with `x = 2*pi^2*f/mu`. The isolated sigmoid edge's "
             "10–90% rise time is `2*ln(9)/mu`. The independent pulse-edge Fourier calculation "
             "also includes PWM sampling and alignment.", ""]
    rows = []
    for run in runs:
        converters = run["converter_validation"]
        name = sorted(converters)[0]
        validation = converters[name]
        harmonics = {h["harmonic"]: h for h in validation["selected_harmonics"]}
        rows.append([number(run["mu"], 8), number(attenuation(frequency, run["mu"]), 7),
                     number(2000 * math.log(9) / run["mu"], 6),
                     number(validation["dc_voltage_v"], 8),
                     *[number(harmonics[h]["predicted_peak_v"], 7) for h in (1, 13, 15, 17)]])
    text += [table(["Mu", "60 Hz gain", "Edge rise (ms)", "DC voltage (V)",
                    "60 Hz peak (V)", "780 Hz peak (V)", "900 Hz peak (V)", "1020 Hz peak (V)"], rows), "",
             "Peak predictions refer to one converter phase after common-mode removal. The converters have "
             "identical PWM settings; the 900 Hz triplen carrier cancels from their phase voltages. "
             "The 780 Hz and 1020 Hz sidebands remain. These are open-loop ideal-DC sources behind RL reactors, "
             "without PLLs, current regulation, DC-link dynamics, or protection. Large reactive absorption "
             "at low mu is behavior of this ideal-source model, not credible protected IBR operation.", "",
             "[Overlaid switching waveforms](plots/switching_waveforms.svg) · "
             "[Harmonic spectra](plots/switching_spectrum.svg)", "",
             "## Measured operating windows", "",
             f"`V` is {remote_bus}'s three-phase mean-square line-to-line RMS voltage on the "
             f"{summary['voltage_ll_v']/1000:g} kV base. `P` is the mean of instantaneous terminal abc power; "
             "`Q60` uses RMS 60 Hz phasors. Positive power means injection into the network. Reactor current "
             "is `sqrt(mean(Ia_rms^2, Ib_rms^2, Ic_rms^2))`; values include all resolved waveform content. "
             "These finite windows do not establish equilibrium.", ""]
    rows, machine_rows = [], []
    for run in runs:
        for name, window in run["windows"].items():
            label = window_label(name, window)
            rows.append([number(run["mu"], 8), label,
                         number(window["bus_voltage"][remote_bus]["mean_square_line_to_line_rms_pu"], 7),
                         number(window["total_ibr_p_mw"], 7), number(window["total_ibr_fundamental_q_mvar"], 7),
                         " / ".join(number(rms_reactor(window, reactor), 6) for reactor in reactors)])
            machine_rows.append([number(run["mu"], 8), label,
                                 *[number(window["machines"][machine]["frequency_mean_hz"], 9) for machine in machines],
                                 *[number(window["machines"][machine]["efd_mean_pu"], 7) for machine in machines]])
    text += [table(["Mu", "Window", "V (pu)", "Total IBR P (MW)", "Total IBR Q60 (Mvar)",
                    "Reactor RMS A: " + " / ".join(reactors)], rows), "",
             "Machine frequency and field voltage below are window means. `Efd` uses each machine's own per-unit base.", "",
             table(["Mu", "Window", *[f"{machine} f (Hz)" for machine in machines],
                    *[f"{machine} Efd (pu)" for machine in machines]], machine_rows), ""]
    if len(runs) > 1:
        low, high = runs[0]["windows"]["final"], baseline["windows"]["final"]
        dv = 100 * (low["bus_voltage"][remote_bus]["mean_square_line_to_line_rms_pu"]
                    - high["bus_voltage"][remote_bus]["mean_square_line_to_line_rms_pu"])
        text += [f"In the final window, mu={number(runs[0]['mu'], 8)} minus mu={number(baseline_mu, 8)} "
                 f"changes {remote_bus} voltage by {number(dv, 6)} percentage points of rated voltage, "
                 f"total IBR active injection by {number(low['total_ibr_p_mw'] - high['total_ibr_p_mw'], 6)} MW, "
                 f"and total IBR fundamental reactive injection by "
                 f"{number(low['total_ibr_fundamental_q_mvar'] - high['total_ibr_fundamental_q_mvar'], 6)} Mvar.", ""]
    text += ["[Response overlays](plots/response_comparison.svg) · "
             "[Event waveforms](plots/event_waveforms.svg) · "
             "[Differences from the high-mu run](plots/response_differences.svg)", "",
             "## Tolerance sensitivity versus mu effects", "",
             "Each entry below is the largest absolute waveform difference over the indicated monitor group. "
             "The mu effect compares different smoothed models; the tolerance check compares the same mu "
             "at primary and tighter integration tolerances. Their maxima can occur at different nodes or times. "
             "Tolerance sensitivity is evidence about numerical convergence, not an exact error bound.", ""]
    comparison_rows = []
    groups = [("bus_voltage", "Bus phase voltage"), ("ibr_current", "IBR reactor current"),
              ("machine_frequency", "Machine frequency"), ("machine_field_voltage", "Machine Efd")]
    for run in runs:
        for region, label in (("after_0p1s", "t >= 0.1 s"), ("final_window", "Final window")):
            for group, quantity in groups:
                effect = run["difference_from_high_mu"][region][group]["worst_absolute"]
                refined_value = run.get("tolerance_refinement", {}).get(region, {}).get(group)
                error = None if refined_value is None else refined_value["worst_absolute"]["max_abs"]
                ratio = None if error is None or error == 0 or run["mu"] == baseline_mu else effect["max_abs"] / error
                comparison_rows.append([number(run["mu"], 8), label, f"{quantity} ({effect['unit']})",
                                        number(effect["max_abs"], 6), number(error, 6), number(ratio, 5)])
    text += [table(["Mu", "Region", "Quantity", "Max difference from high mu", "Max same-mu tolerance change",
                    "Mu effect / tolerance change"], comparison_rows), ""]
    if not refined:
        text += ["Same-mu tighter-tolerance results are missing for at least one run; no numerical convergence conclusion is made for those rows.", ""]
    text += ["## Independent checks", ""]
    validation_rows = []
    for run in runs:
        checks = list(run["converter_validation"].values())
        validation_rows.append([number(run["mu"], 8), number(run["kcl"]["max_abs_a"], 6),
                                number(max(c["all_harmonics_1_to_49_max_abs_error_v"] for c in checks), 6),
                                number(max(c["all_harmonics_1_to_49_max_gate_abs_error"] for c in checks), 6),
                                number(max(c["common_mode_sum_max_abs_v"] for c in checks), 6)])
    text += [table(["Mu", "Max sampled KCL residual (A)", "Max converter harmonic error (V peak)",
                    "Max gate harmonic error", "Max abs(va+vb+vc) (V)"], validation_rows), "",
             "KCL sums every monitored source, load, machine, line-series, line-shunt, and switch terminal "
             "injection at every recorded time. PWM validation compares independently integrated ideal pulse "
             "edges with measured final-window Fourier coefficients for all harmonics 1–49. It validates "
             "the applied PWM/bridge waveform, not the whole network solution.", "",
             "## Physical line coupling and fit", "",
             f"OpenLine supplied untransposed {physical['conductor']} parameters at "
             f"{physical['conductor_temperature_C']:g} C and earth resistivity "
             f"{physical['earth_resistivity_ohm_m']:g} ohm m. Schelkunoff internal conductor impedance, "
             "Carson earth return, external inductance, and Maxwell capacitance retain all nine phase-matrix entries. "
             f"The {len(lines)} pi sections span {min(line['params']['dx'] for line in lines):g}–"
             f"{max(line['params']['dx'] for line in lines):g} m.", ""]
    z = impedance(line_model["Zp"], frequency) * 1000
    full = sequence_impedance(z)
    diagonal = sequence_impedance(np.diag(np.diag(z)))
    text += [f"At {frequency:g} Hz, phase-A self impedance is {number(z[0, 0].real, 6)} + "
             f"j{number(z[0, 0].imag, 6)} ohm/km; A–B mutual impedance is "
             f"{number(z[0, 1].real, 6)} + j{number(z[0, 1].imag, 6)} ohm/km. "
             "Setting mutual terms to zero produces the following analytical change; no separate diagonal-only simulation is implied.", "",
             table(["Impedance matrix", "R0 (ohm/km)", "X0 (ohm/km)", "R1 (ohm/km)", "X1 (ohm/km)"],
                   [[label, *[number(value, 7) for value in (seq[0, 0].real, seq[0, 0].imag,
                                                           seq[1, 1].real, seq[1, 1].imag)]]
                    for label, seq in (("Full phase coupling", full), ("Self terms only", diagonal))]), "",
             f"Removing mutual terms changes the positive-sequence diagonal reactance by "
             f"{number(100*(diagonal[1, 1].imag/full[1, 1].imag-1), 6)}% and the zero-sequence "
             f"diagonal reactance by {number(100*(diagonal[0, 0].imag/full[0, 0].imag-1), 6)}%. "
             "These are projections of an untransposed matrix; off-diagonal sequence coupling also remains in the full model.", "",
             f"The passive rational impedance has {len(line_model['Zp']['poles'])} negative real poles. "
             f"Maximum held-out entrywise error is {number(100*fit['max_held_out_entry_relative_error'], 6)}%; "
             f"RMS entrywise error is {number(100*fit['rms_held_out_entry_relative_error'], 6)}%, over "
             f"{fit['frequency_band_Hz'][0]:g}–{fit['frequency_band_Hz'][1]:g} Hz. "
             "Passivity follows from positive-semidefinite Foster coefficient matrices. "
             f"The exact full capacitance representation has maximum relative error "
             f"{number(fit['shunt_max_entry_relative_error'], 5)}; an independent electrostatic "
             f"potential-matrix inversion agrees to {number(fit['independent_capacitance_max_entry_relative_error'], 5)} relative.", ""]
    pi_errors = fit.get("pi_600m_relative_nodal_admittance_Frobenius_error_by_Hz", {})
    if pi_errors:
        text += ["For the longest section, comparison with the exact uniform-line telegrapher transfer matrix "
                 "using identical per-meter parameters isolates pi-section discretization error:", "",
                 table(["Frequency (Hz)", "600 m relative terminal-admittance matrix error (%)"],
                       [[key, number(100*value, 6)] for key, value in sorted(pi_errors.items(), key=lambda item: float(item[0]))]), ""]
    text += ["[Conductor geometry and coupled matrices](line_geometry_and_coupling.svg) · "
             "[Raw versus fitted self/mutual line responses](line_frequency_response.svg)", "",
             "## Limits and provenance", "",
             "The common initial bus voltages and machine dispatch are a fundamental-frequency estimate; "
             "line and filter memories start from model defaults. This launches artificial energization "
             "transients. The small physical shunt capacitances create modes above the 30 kHz line-fit band "
             "and the monitor Nyquist frequency. IDA resolves those modes internally, but recorded startup "
             "waveforms cannot validate them. Ideal switch edges also excite frequencies beyond the validated "
             "band. Later-window comparisons are the primary operating evidence; sub-millisecond peaks "
             "carry no physical high-frequency accuracy claim. Propagation and delay models are not used.", "",
             f"Line calculator: local `/home/lukel/openline`, "
             f"[OpenLine revision `{source['revision']}`]({source['repository']}/tree/{source['revision']}). "
             "The requested GridWorkbench checkout could not be located; OpenLine was used as the available "
             "home-repository parameter calculator. Home repositories and their branches were not modified.", "",
             f"Raw line sweep SHA-256: `{source['sweep_sha256']}`.", ""]
    provenance = summary.get("provenance", {})
    for key, label in (("git_revision", "GridKit base revision"), ("binary_sha256", "Binary SHA-256"),
                       ("case_sha256", "Case SHA-256"), ("state_sha256", "Initial-state SHA-256"),
                       ("platform", "Platform"), ("cpu_model", "CPU"),
                       ("build_type", "Build type"), ("compiler", "Compiler"),
                       ("sundials_version", "SUNDIALS version"), ("enzyme_enabled", "Enzyme enabled")):
        if provenance.get(key):
            text.append(f"- {label}: `{provenance[key]}`")
    if provenance.get('shared_libraries_sha256'):
        text += ['', f"The manifest also records SHA-256 hashes for {len(provenance['shared_libraries_sha256'])} loaded shared libraries."]
    text += ["", "Timing and validation values above were read from the linked manifest and metrics; "
             "raw monitor CSVs, solver configurations, and logs remain in each run directory.", ""]
    return "\n".join(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--case", type=Path, default=HERE / "CoupledGrid.case.json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.results.resolve()
    summary = json.loads((root / "summary.json").read_text())
    metrics = json.loads((root / "metrics.json").read_text())
    case_text = args.case.read_bytes()
    recorded_hash = summary.get("provenance", {}).get("case_sha256")
    if recorded_hash and hashlib.sha256(case_text).hexdigest() != recorded_hash:
        raise ValueError("Current case differs from the simulated case hash")
    case = json.loads(case_text)
    fit = json.loads((HERE / "lines/fit_validation.json").read_text())
    line_model = json.loads((HERE / "lines/model.json").read_text())
    if any(device["submodels"] != line_model for device in case["devices"] if device["class"] == "LineLumped"):
        raise ValueError("Case line coefficients differ from the validated physical model")
    output = args.output or root / "report.md"
    text = make_report(summary, metrics, case, fit, line_model)
    test_log = root / 'verification/emt_ctest.log'
    if test_log.exists():
        test_summary = re.search(r'^\d+% tests passed, .*$', test_log.read_text(), re.MULTILINE)
        if test_summary:
            text += '\n## Repository checks\n\n' + test_summary[0] + '. '
            text += '[Captured EMT CTest output](verification/emt_ctest.log).\n'
    output.write_text(text)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
