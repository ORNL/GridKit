"""Independent traveling-wave and frequency-domain checks; no C++ state equations.

Usage: python3 validate.py results [--executable path/to/line_distributed]
The optional executable also runs tolerance/compression and timing comparisons.
"""

import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np


def matrix(values):
    return np.array([[complex(*x) if isinstance(x, list) else complex(x)
                      for x in row] for row in values])


def fit(data, s):
    s = np.atleast_1d(s)
    result = np.broadcast_to(matrix(data.get("D", np.zeros((data["rows"], data["cols"])))),
                             (len(s), data["rows"], data["cols"])).copy()
    for pole, residue in zip(data.get("poles", []), data.get("residues", [])):
        pole = complex(*pole) if isinstance(pole, list) else pole
        result += matrix(residue)[None] / (s - pole)[:, None, None]
    return result


def operators(data, s):
    line = data["line"]
    yc = fit(line["Yc"], s)
    h = line["H"]
    if "fit" in h:
        propagation = fit(h["fit"], s) * np.exp(-np.atleast_1d(s) * h["delay"])[:, None, None]
    else:
        propagation = ((fit(h["output"], s) * np.exp(-np.outer(np.atleast_1d(s), h["delays"]))[:, None, :])
                       @ fit(h["input"], s))
    return yc, propagation


def frequency_solution(data, s, excitation):
    """Solve terminal KCL and traveling-wave equations directly in frequency.

    No rational memory states or delay interpolation enter this calculation.
    Keeping wave variables avoids a singular nodal stamp at lossless DC.
    """
    s = np.atleast_1d(s)
    n, k = data["line"]["N"], data["line"]["K"]
    yc, h = operators(data, s)
    p = np.zeros((n, k))
    p[np.asarray(data["line"]["conductors"]) - 1, np.arange(k)] = 1
    source_g = np.zeros((len(s), n, n), complex)
    for phase in range(n):
        source_g[:, phase, phase] = 1/(data["source"]["resistance"][phase]
            + s*data["source"].get("inductance", [0]*n)[phase])
    load_g = 1 / (np.asarray(data["load"]["resistance"])[None]
                  + s[:, None] * np.asarray(data["load"]["inductance"])[None])
    size = 2 * n + 2 * k
    a = np.zeros((len(s), size, size), dtype=complex)
    rhs = np.zeros((len(s), size), dtype=complex)
    v1, v2, r1, r2 = slice(0, n), slice(n, 2*n), slice(2*n, 2*n+k), slice(2*n+k, size)
    yphase = p @ yc @ p.T
    a[:, v1, v1] = -source_g - yphase
    a[:, v1, r2] = p @ h
    a[:, v2, v2] = -yphase
    for phase in range(n):
        a[:, n+phase, n+phase] -= load_g[:, phase]
    a[:, v2, r1] = p @ h
    a[:, r1, v1] = 2 * yc @ p.T
    a[:, r1, r1] = -np.eye(k)
    a[:, r1, r2] = -h
    a[:, r2, v2] = 2 * yc @ p.T
    a[:, r2, r2] = -np.eye(k)
    a[:, r2, r1] = -h
    rhs[:, v1] = -np.einsum("fij,fj->fi", source_g, np.broadcast_to(excitation, (len(s), n)))
    return np.linalg.solve(a, rhs)


def phasors(data):
    source = data["source"]
    e = np.sqrt(2) * np.asarray(source["rms"]) * np.exp(1j * np.asarray(source["phase"]))
    return frequency_solution(data, 2j * np.pi * source["frequency"], e)[0]


def ramp(t, start, duration):
    x = np.clip((t - start) / duration, 0, 1)
    return x**3 * (10 + x * (-15 + 6*x))


def forcing(data, t):
    t = np.atleast_1d(t)
    amplitude = np.full((len(t), 3), float(data["initialization"] == "sinusoidal"))
    fault = np.zeros_like(amplitude)
    pulse = np.zeros_like(amplitude)
    for event in data["events"]:
        envelope = ramp(t, event["on"], event["rise"]) - ramp(t, event["off"], event["fall"])
        array = {"source": amplitude, "fault": fault, "voltage": pulse}[event["type"]]
        array[:, np.asarray(event["phases"]) - 1] += event["value"] * envelope[:, None]
    source = data["source"]
    voltage = amplitude * np.sqrt(2) * np.asarray(source["rms"]) * np.cos(
        2*np.pi*source["frequency"] * t[:, None] + np.asarray(source["phase"]))
    return voltage + pulse, fault


def bergeron_reference(data, t):
    """Exact characteristic recursion for the uncoupled constant-Yc R circuit.

    Evaluates each source/event at its exact retarded time, without a time grid.
    """
    tau = data["line"]["H"]["delay"]
    g = np.diag(matrix(data["line"]["Yc"]["D"])).real
    rs = np.asarray(data["source"]["resistance"])
    gl = 1 / np.asarray(data["load"]["resistance"])
    initial = phasors(data) if data["initialization"] == "sinusoidal" else np.zeros(12, complex)
    omega = 2*np.pi*data["source"]["frequency"]
    r1 = np.zeros((len(t), 3))
    r2 = r1.copy()
    depth = int(np.ceil((max(t) - data["solver"]["t0"]) / tau)) + 1
    for i in range(depth, -1, -1):
        time = t - i*tau
        e, fault = forcing(data, time)
        v1 = (e/rs + r2) / (1/rs + g)
        v2 = r1 / (gl + fault + g)
        r1, r2 = 2*g*v1 - r2, 2*g*v2 - r1
        before = time < data["solver"]["t0"]
        rotation = np.exp(1j*omega*time[before, None])
        r1[before] = (rotation * initial[6:9]).real
        r2[before] = (rotation * initial[9:12]).real
    return np.column_stack((v1, v2))


def spectral_reference(data, t, dt):
    """Independent inverse-Fourier response to a compact source disturbance.

    Zero-padding separates the response from periodic copies. Two resolutions
    are compared by the caller to check this numerical reference itself.
    """
    assert all(e["type"] == "source" for e in data["events"])
    assert data["initialization"] == "sinusoidal" and data["solver"]["t0"] == 0
    count = 1 << int(np.ceil(np.log2(4 * data["solver"]["tf"] / dt)))
    grid = np.arange(count)*dt
    e, _ = forcing(data, grid)
    baseline = copy.deepcopy(data)
    baseline["events"] = []
    e0, _ = forcing(baseline, grid)
    spectrum = np.fft.rfft(e - e0, axis=0)
    s = 2j*np.pi*np.fft.rfftfreq(count, dt)
    result = np.empty((len(s), 6), complex)
    for start in range(0, len(s), 512):
        stop = start + 512
        result[start:stop] = frequency_solution(data, s[start:stop], spectrum[start:stop])[:, :6]
    response = np.fft.irfft(result, n=count, axis=0)
    transient = np.column_stack([np.interp(t, grid, response[:, i]) for i in range(6)])
    return transient + (np.exp(2j*np.pi*data["source"]["frequency"]*t[:, None]) * phasors(data)[:6]).real


def terminal_voltages(samples):
    return np.column_stack([samples[f"v{end}_{phase}"] for end in (1, 2) for phase in "abc"])


def load(path):
    return json.loads((path / "input.json").read_text()), np.genfromtxt(path / "samples.csv", delimiter=",", names=True)


def errors(actual, expected, scale):
    difference = actual - expected
    return {"max_V": float(np.max(np.abs(difference))),
            "rms_V": float(np.sqrt(np.mean(difference**2))),
            "relative_peak": float(np.max(np.abs(difference)) / scale)}


def passivity(data):
    s = 2j*np.pi*np.geomspace(.01, 1e7, 2000)
    yc, h = operators(data, s)
    k = data["line"]["K"]
    waves = np.zeros((len(s), 2*k, 2*k), complex)
    y = waves.copy()
    waves[:, :k, :k] = waves[:, k:, k:] = np.eye(k)
    waves[:, :k, k:] = waves[:, k:, :k] = h
    y[:, :k, :k] = y[:, k:, k:] = yc
    # Lossless nodal admittance has poles at resonances. Do not interpret
    # cancellation in an ill-conditioned conversion as a passivity violation.
    condition = np.linalg.cond(waves)
    keep = condition < 1e9
    condition, waves, y = condition[keep], waves[keep], y[keep]
    admittance = 2*np.linalg.solve(waves, y) - y
    hermitian = (admittance + admittance.conj().transpose(0, 2, 1))/2
    minimum = np.linalg.eigvalsh(hermitian)[:, 0]
    roundoff = 64*np.finfo(float).eps*condition*np.maximum(1, np.max(np.abs(admittance), axis=(1, 2)))
    return {"sampled_min_eigenvalue_S": float(minimum.min()),
            "violations_beyond_roundoff": int(np.count_nonzero(minimum < -roundoff-1e-10)),
            "max_roundoff_allowance_S": float(roundoff.max()),
            "excluded_nodal_poles": int(np.count_nonzero(~keep)),
            "reciprocity_error_S": float(np.max(np.abs(admittance - admittance.transpose(0, 2, 1)))),
            "frequency_Hz": [.01, 1e7], "points": 2000}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--executable", type=Path)
    args = parser.parse_args()
    (args.results / "validation.json").unlink(missing_ok=True)
    report = {}
    for path in sorted(args.results.iterdir()):
        if not (path / "statistics.json").exists():
            continue
        data, samples = load(path)
        actual = terminal_voltages(samples)
        scale = np.sqrt(2)*max(data["source"]["rms"])
        assert np.isfinite(actual).all(), path.name
        entry = {"passivity": passivity(data)}
        stats = json.loads((path / "statistics.json").read_text())
        assert stats["history_accepted_steps"] == stats["steps"]
        assert entry["passivity"]["violations_beyond_roundoff"] == 0
        assert entry["passivity"]["reciprocity_error_S"] < 1e-8
        source_kcl = np.column_stack([(samples[f"e_{p}"] - samples[f"v1_{p}"])/data["source"]["resistance"][i]
                                     + samples[f"i1_{p}"] for i, p in enumerate("abc")])
        load_kcl = np.column_stack([samples[f"i2_{p}"] - samples[f"load_{p}"] - samples[f"fault_{p}"] for p in "abc"])
        entry["kcl_max_A"] = float(max(np.abs(source_kcl).max(), np.abs(load_kcl).max()))
        assert entry["kcl_max_A"] < .02, (path.name, entry)
        reference = None
        if "delay" in data["line"]["H"] and not any(data["load"]["inductance"]):
            reference = bergeron_reference(data, samples["t"])
            entry["reference"] = "Exact retarded-time Bergeron recursion"
        elif not data["events"]:
            reference = (np.exp(2j*np.pi*data["source"]["frequency"] * samples["t"][:, None]) * phasors(data)[:6]).real
            entry["reference"] = "Direct frequency-domain circuit solution"
        elif all(e["type"] == "source" for e in data["events"]):
            reference = spectral_reference(data, samples["t"], 2e-6)
            finer = spectral_reference(data, samples["t"], 1e-6)
            entry["reference_refinement"] = errors(reference, finer, scale)
            assert entry["reference_refinement"]["relative_peak"] < 1e-5
            reference = finer
            entry["reference"] = "Inverse Fourier transform of terminal transfer"
        if reference is not None:
            entry["voltage_error"] = errors(actual, reference, scale)
            assert entry["voltage_error"]["relative_peak"] < 2e-4, (path.name, entry)
            np.savetxt(path / "reference.csv", np.column_stack((samples["t"], reference)), delimiter=",",
                       header="t,v1_a,v1_b,v1_c,v2_a,v2_b,v2_c", comments="")
        if path.name in ("energization", "matched"):
            first = data["events"][0]["on"] + data["line"]["H"]["delay"]
            entry["prearrival_max_V"] = float(np.max(np.abs(actual[samples["t"] < first-1e-12, 3:])))
            assert entry["prearrival_max_V"] < .01
        if path.name == "short_delay_sag":
            steps = np.genfromtxt(path / "steps.csv", delimiter=",", names=True)
            entry["max_h_over_tau"] = float(steps["h_over_tau"].max())
            assert entry["max_h_over_tau"] > 1
        report[path.name] = entry

    if args.executable:
        executable = args.executable.resolve()
        report["executable_sha256"] = hashlib.sha256(executable.read_bytes()).hexdigest()
        checks = args.results / "checks"
        checks.mkdir(exist_ok=True)
        for name in ("energization", "frequency_dependent_sag", "coupled_fault", "short_delay_sag"):
            base, original = load(args.results / name)
            comparisons = {}
            for suffix in ("refined", "compressed", "capped"):
                spec = copy.deepcopy(base)
                spec["name"] = f"{name}_{suffix}"
                if suffix == "refined":
                    spec["solver"]["rel_tol"] *= .1
                    spec["solver"]["history_abs_tol"] *= .1
                elif suffix == "compressed":
                    spec["solver"]["compress"] = True
                else:
                    h = spec["line"]["H"]
                    tau = h["delay"] if "delay" in h else min(h["delays"])
                    spec["solver"]["h_max"] = tau/4
                input_path = checks / f"{spec['name']}.json"
                input_path.write_text(json.dumps(spec, indent=2) + "\n")
                subprocess.run([str(executable), str(input_path), str(checks)], check=True)
                _, trial = load(checks / spec["name"])
                interpolated = np.column_stack([np.interp(original["t"], trial["t"], trial[f"v{e}_{p}"])
                                                for e in (1, 2) for p in "abc"])
                comparison = errors(terminal_voltages(original), interpolated, np.sqrt(2)*max(base["source"]["rms"]))
                assert comparison["relative_peak"] < 2e-4, (name, suffix, comparison)
                comparisons[suffix] = comparison
            report[name]["comparisons"] = comparisons
        spec, _ = load(args.results / "frequency_dependent_sag")
        spec["name"] = "timing"
        input_path = checks / "timing.json"
        input_path.write_text(json.dumps(spec, indent=2) + "\n")
        timings = []
        for _ in range(3):
            subprocess.run([str(executable), str(input_path), str(checks)], check=True)
            timings.append(json.loads((checks / "timing/statistics.json").read_text())["seconds_including_csv"])
        report["timing"] = {"seconds_including_csv": timings, "median": float(np.median(timings))}
    (args.results / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
