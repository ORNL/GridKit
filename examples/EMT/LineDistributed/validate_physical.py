"""Physical-line validation, independent unfitted telegrapher reference and studies.

Usage: python3 validate_physical.py [--executable path/to/line_distributed]
All generated reports are written under the source example's results/ directory.
"""

import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from line_constants import HERE, terminal_admittance, wave_delays
from validate import errors, forcing, frequency_solution, load, passivity, terminal_voltages

RESULTS = HERE / "results/physical"


def length_of(data):
    return data["line"]["H"]["delays"][0]/wave_delays(1)[0]


def response(data, s, excitation, raw):
    """Voltages and source/load currents from frequency-domain terminal KCL."""
    s = np.atleast_1d(s)
    e = np.broadcast_to(excitation, (len(s), 3))
    zs = np.array(data["source"]["resistance"]) + s[:, None]*data["source"].get("inductance", [0]*3)
    zl = np.array(data["load"]["resistance"]) + s[:, None]*data["load"]["inductance"]
    if raw:
        a = terminal_admittance(s, length_of(data))
        rhs = np.zeros((len(s), 6), complex)
        rhs[:, :3] = e/zs
        for k in range(3):
            a[:, k, k] += 1/zs[:, k]
            a[:, 3+k, 3+k] += 1/zl[:, k]
        voltage = np.linalg.solve(a, rhs)
    else:
        voltage = frequency_solution(data, s, e)[:, :6]
    return np.column_stack((voltage, (e-voltage[:, :3])/zs, voltage[:, 3:]/zl))


def steady(data, t, raw):
    if data["initialization"] == "zero":
        return np.zeros((len(t), 12))
    source = data["source"]
    amplitude = np.sqrt(2)*np.asarray(source["rms"])*np.exp(1j*np.asarray(source["phase"]))
    jw = 2j*np.pi*source["frequency"]
    return (np.exp(jw*t[:, None])*response(data, jw, amplitude, raw)[0]).real


def fourier_reference(data, t, dt, raw, padding=4):
    """Source perturbation through a fixed linear network, with zero padding."""
    assert all(e["type"] in ("voltage", "source") for e in data["events"])
    count = 1 << int(np.ceil(np.log2(padding*data["solver"]["tf"]/dt)))
    grid = np.arange(count)*dt
    baseline = copy.deepcopy(data)
    baseline["events"] = []
    forcing_spectrum = np.fft.rfft(forcing(data, grid)[0]-forcing(baseline, grid)[0], axis=0)
    s = 2j*np.pi*np.fft.rfftfreq(count, dt)
    spectrum = np.empty((len(s), 12), complex)
    for start in range(0, len(s), 2048):
        stop = start+2048
        spectrum[start:stop] = response(data, s[start:stop], forcing_spectrum[start:stop], raw)
    waveform = np.fft.irfft(spectrum, n=count, axis=0)
    return np.column_stack([np.interp(t, grid, waveform[:, k]) for k in range(12)]) + steady(data, t, raw)


def observed(samples):
    return np.column_stack((terminal_voltages(samples),
                           *[samples[f"source_{p}"] for p in "abc"],
                           *[samples[f"load_{p}"] for p in "abc"]))


def step_metrics(path, start, end):
    data, _ = load(path)
    steps = np.genfromtxt(path / "steps.csv", delimiter=",", names=True)
    delays = np.array(data["line"]["H"]["delays"])
    delays = delays[delays > 0]
    # Clip intervals at window boundaries so elapsed-time fractions cover exactly
    # the requested window (a sample-count fraction would overweight tiny steps).
    duration = np.maximum(0, np.minimum(steps["t"], end)-np.maximum(steps["t"]-steps["h"], start))
    keep = duration > 0
    h, duration = steps["h"][keep], duration[keep]
    assert abs(duration.sum()-(end-start)) < 1e-10
    above = h > max(delays)
    streak = longest = 0
    for width, accepted in zip(duration, above):
        streak = streak+width if accepted else 0
        longest = max(longest, streak)
    return {"window_s": [start, end], "steps": len(h), "delays_s": delays.tolist(),
            "median_h_over_tau_min": float(np.median(h/min(delays))),
            "median_h_over_tau_max": float(np.median(h/max(delays))),
            "max_h_over_tau_max": float(max(h/max(delays))),
            "time_fraction_above_all_delays": float(duration[above].sum()/duration.sum()),
            "step_fraction_above_all_delays": float(above.mean()),
            "longest_interval_above_all_delays_s": float(longest)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path)
    args = parser.parse_args()
    report_path = RESULTS / "validation.json"
    report_path.unlink(missing_ok=True)
    suite = json.loads((HERE / "physical.json").read_text())
    report = {}
    for run in suite["runs"]:
        name = run["name"]
        path = RESULTS / name
        stats = json.loads((path / "statistics.json").read_text())
        data, samples = load(path)
        actual, t = observed(samples), samples["t"]
        assert np.isfinite(actual).all()
        scale = max(np.max(np.abs(forcing(data, t)[0])), 1)
        entry = {"passivity": passivity(data), "statistics": stats}
        assert entry["passivity"]["violations_beyond_roundoff"] == 0
        assert entry["passivity"]["reciprocity_error_S"] < 1e-8
        assert stats["history_accepted_steps"] == stats["steps"]
        kcl = np.column_stack([samples[f"source_{p}"]+samples[f"i1_{p}"] for p in "abc"]
                              + [samples[f"i2_{p}"]-samples[f"load_{p}"]-samples[f"fault_{p}"] for p in "abc"])
        entry["kcl_max_A"] = float(abs(kcl).max())
        assert entry["kcl_max_A"] < .02
        if all(e["type"] != "fault" for e in data["events"]):
            dt = .5e-6 if data["initialization"] == "zero" else 2e-6
            fitted = fourier_reference(data, t, dt, False)
            raw = fourier_reference(data, t, dt, True)
            finer = fourier_reference(data, t, dt/2, True, padding=8)
            entry["reference_refinement"] = errors(raw[:, :6], finer[:, :6], scale)
            entry["numerical_error"] = errors(actual[:, :6], fitted[:, :6], scale)
            entry["physical_model_error"] = errors(actual[:, :6], finer[:, :6], scale)
            entry["current_error_A"] = float(np.max(abs(actual[:, 6:]-finer[:, 6:])))
            assert entry["current_error_A"] < .005*max(np.max(abs(finer[:, 6:])), 1), (name, entry)
            assert entry["reference_refinement"]["relative_peak"] < 2e-4, (name, entry)
            assert entry["numerical_error"]["relative_peak"] < 2e-4, (name, entry)
            assert entry["physical_model_error"]["relative_peak"] < .005, (name, entry)
            np.savetxt(path / "reference.csv", np.column_stack((t, finer)), delimiter=",", comments="",
                       header="t,v1_a,v1_b,v1_c,v2_a,v2_b,v2_c,source_a,source_b,source_c,load_a,load_b,load_c")
            if data["initialization"] == "zero":
                arrival = min(e["on"] for e in data["events"])+min(data["line"]["H"]["delays"])
                entry["prearrival_max_V"] = float(abs(actual[t < arrival-1e-12, 3:6]).max())
                assert entry["prearrival_max_V"] < .01
        else:
            entry["reference"] = "Time-varying fault: refinement/cap checks, no independent transient reference"
        if data["solver"]["tf"] >= .2:
            mask = t >= .1
            entry["settled_voltage_error"] = errors(actual[mask, :6], steady(data, t[mask], False)[:, :6], scale)
            entry["settled_steps"] = step_metrics(path, .1, .2)
            assert entry["settled_voltage_error"]["relative_peak"] < 2e-4, (name, entry)
        if name == "physical_short_recovery":
            assert entry["settled_steps"]["time_fraction_above_all_delays"] > .9, entry
        report[name] = entry
        print(name, json.dumps(entry), flush=True)

    if args.executable:
        executable = args.executable.resolve()
        checks = RESULTS / "checks"
        checks.mkdir(exist_ok=True)
        report["executable_sha256"] = hashlib.sha256(executable.read_bytes()).hexdigest()
        for name in ("physical_pulse", "physical_fault", "physical_short_recovery"):
            base, original = load(RESULTS / name)
            comparisons = {}
            for variant in ("refined", "capped", "compressed"):
                spec = copy.deepcopy(base)
                spec["name"] = name+"_"+variant
                if variant == "refined":
                    spec["solver"]["rel_tol"] *= .1
                    spec["solver"]["history_abs_tol"] *= .1
                elif variant == "capped":
                    spec["solver"]["h_max"] = min(base["line"]["H"]["delays"])/4
                else:
                    spec["solver"]["compress"] = True
                input_path = checks / (spec["name"]+".json")
                input_path.write_text(json.dumps(spec)+"\n")
                subprocess.run([str(executable), str(input_path), str(checks)], check=True)
                _, trial = load(checks / spec["name"])
                interpolated = np.column_stack([np.interp(original["t"], trial["t"], trial[f"v{e}_{p}"])
                                                for e in (1, 2) for p in "abc"])
                scale = max(np.max(abs(forcing(base, original["t"])[0])), 1)
                comparisons[variant] = errors(terminal_voltages(original), interpolated, scale)
                comparisons[variant]["statistics"] = json.loads((checks / spec["name"] / "statistics.json").read_text())
                assert comparisons[variant]["relative_peak"] < 2e-4, (name, variant, comparisons)
            report[name]["comparisons"] = comparisons
        # Timings use identical monitors and include initialization-excluded CSV
        # runtime, matching the executable's existing statistics contract.
        times = {"uncapped": [], "capped": []}
        for label in times:
            base, _ = load(RESULTS / "physical_short_recovery")
            base["name"] = "timing_"+label
            if label == "capped":
                base["solver"]["h_max"] = min(base["line"]["H"]["delays"])/4
            path = checks / (base["name"]+".json")
            path.write_text(json.dumps(base)+"\n")
            for _ in range(3):
                subprocess.run([str(executable), str(path), str(checks)], check=True)
                times[label].append(json.loads((checks / base["name"] / "statistics.json").read_text())["seconds_including_csv"])
        report["timing"] = {"trials_s": times, "median_s": {k: float(np.median(v)) for k, v in times.items()}}
    report_path.write_text(json.dumps(report, indent=2)+"\n")
    print("Validated:", report_path)


if __name__ == "__main__":
    main()
