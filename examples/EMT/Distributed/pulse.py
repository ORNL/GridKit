#!/usr/bin/env python3
"""Run and validate the 60 km frequency-dependent square-pulse comparison."""

import argparse
import csv
import hashlib
import json
import platform
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import scipy
from scipy.optimize import nnls

from run import ROOT, frequency_line, run_case, write

LENGTH = 60000.0
START = 0.001
WIDTH = 100e-6
AMPLITUDE = 1000.0
HORIZON = 0.003
SPACING = 0.5e-6
SECTIONS = (1, 5, 20, 100)
MODELS = ("distributed",) + tuple(f"pi_{n}" for n in SECTIONS)
GEOMETRY = ROOT / "cases/EMT/Distributed/overhead.line.json"


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scalar_response(fit, s, zero=True):
    """Evaluate one eigenvalue of a transposed phase-matrix fit."""

    def modal(matrix):
        return np.sum(matrix[0]) if zero else matrix[0, 0] - matrix[0, 1]

    value = np.full(s.shape, modal(np.array(fit["D"])), complex)
    value += s * modal(np.array(fit.get("E", np.zeros((3, 3)))))
    for pole, residue in zip(fit.get("poles", []), fit.get("residues", [])):
        r = np.array(residue)
        value += modal(r[..., 0] + 1j * r[..., 1]) / (s - complex(*pole))
    return value


def fit_lumped(fits, output, gridworkbench=None):
    """Fit Z'=R+sL+sum r*s/(s+p), with nonnegative modal coefficients.

    Each term is a passive parallel RL branch in series with R and L.
    Positive semidefinite modal projectors preserve matrix positive realness.
    Y'=sC is exact for this geometry's frequency-independent shunt data.
    """
    raw = np.load(fits / "parameters.npz")
    w, z, y = raw["omega"], raw["Z"], raw["Y"]
    s = 1j * w
    # Avoid unnecessarily fast poles far above the 10 MHz data band; those
    # worsen DAE conditioning without improving the pulse-band fit.
    p = 2 * np.pi * np.geomspace(0.001, 3e7, 32)
    basis = np.column_stack((np.ones(len(s)), s, s[:, None] / (s[:, None] + p)))
    p0 = np.ones((3, 3)) / 3
    zp = {
        "rows": 3,
        "cols": 3,
        "D": np.zeros((3, 3)),
        "E": np.zeros((3, 3)),
        "poles": [],
        "residues": [],
    }
    metrics = {}
    for name, projector, values in (
        ("aerial", np.eye(3) - p0, z[:, 0, 0] - z[:, 0, 1]),
        ("zero", p0, z[:, 0, 0] + 2 * z[:, 0, 1]),
    ):
        a, b = basis / abs(values[:, None]), values / abs(values)
        a, b = np.vstack((a.real, a.imag)), np.r_[b.real, b.imag]
        scale = np.linalg.norm(a, axis=0)
        coefficients = nnls(a / scale, b, maxiter=2000)[0] / scale
        resistance, inductance, *weights = coefficients
        if resistance <= 0 or inductance <= 0:
            raise ValueError("The series fit requires positive modal R and L")
        zp["D"] += (resistance + sum(weights)) * projector
        zp["E"] += inductance * projector
        for rate, weight in zip(p, weights):
            if weight > 0:
                zp["poles"].append([-float(rate), 0.0])
                residue = -weight * rate * projector
                zp["residues"].append(
                    np.stack((residue, np.zeros((3, 3))), axis=-1).tolist()
                )
        error = abs(basis @ coefficients - values) / abs(values)
        metrics[name] = {
            "relative_max": float(max(error)),
            "relative_rms": float(np.sqrt(np.mean(error**2))),
            "R_ohm_per_m": float(resistance),
            "L_H_per_m": float(inductance),
            "positive_RL_branches": int(np.count_nonzero(weights)),
        }
        if max(error) > 2e-4:
            raise ValueError(f"{name} series fit exceeds 0.02% error")
    for key in ("D", "E"):
        zp[key] = zp[key].tolist()
    capacitance = (y[0] / s[0]).real
    if not np.allclose(y, s[:, None, None] * capacitance, rtol=1e-10, atol=1e-20):
        raise ValueError(
            "This example requires frequency-independent shunt capacitance"
        )
    if min(np.linalg.eigvalsh(capacitance)) <= 0:
        raise ValueError("Shunt capacitance must be positive definite")
    yp = {
        "rows": 3,
        "cols": 3,
        "D": np.zeros((3, 3)).tolist(),
        "E": capacitance.tolist(),
    }
    if gridworkbench is not None:
        sys.path.insert(0, str(gridworkbench.resolve()))
        from gridworkbench.emt.parameters import load, sweep

        # Recompute the geometry between the original fitting frequencies.
        held_w = np.sqrt(w[:-1] * w[1:])
        fresh = sweep(load(GEOMETRY), held_w)
        physical = fresh.R + 1j * held_w[:, None, None] * fresh.L
        held_z = (
            sum(np.roll(np.roll(physical, k, axis=1), k, axis=2) for k in range(3)) / 3
        )
        for name, zero in (("aerial", False), ("zero", True)):
            target = (
                held_z[:, 0, 0] + 2 * held_z[:, 0, 1]
                if zero
                else held_z[:, 0, 0] - held_z[:, 0, 1]
            )
            error = abs(scalar_response(zp, 1j * held_w, zero) - target) / abs(target)
            metrics[name].update(
                holdout_relative_max=float(max(error)), holdout_samples=len(held_w)
            )
            if max(error) > 2e-4:
                raise ValueError(
                    f"{name} series fit exceeds 0.02% on the independent geometry grid"
                )
    write(output / "zp.json", zp)
    write(output / "yp.json", yp)
    write(output / "fit_quality.json", metrics)
    return metrics


def make_case(kind, fits, coefficients):
    case, state, events, meta = frequency_line(kind, fits)
    case["header"] = {
        "case_name": "60 km frequency-dependent square pulse",
        "case_description": "1 kV, 100 us equal-phase source pulse; 100 ohm source, 600 ohm load",
        "case_comments": "Illustrative GridWorkbench geometry, cyclically transposed; frequency-dependent pi sections",
    }
    for device in case["devices"]:
        if device["class"] == "LineLumped":
            device["params"] = {
                "N": 3,
                "K": 3,
                "conductors": [1, 2, 3],
                "dx": LENGTH / int(kind.split("_")[1]),
            }
            device["submodels"] = {
                "Zp": read(coefficients / "zp.json"),
                "Yp": read(coefficients / "yp.json"),
            }
        elif device["class"] == "LineDistributed":
            device["mon"] = [
                f"{signal}{phase}" for signal in ("i_ref1", "i_inc2") for phase in "abc"
            ]
        elif device["class"] == "Bus" and device["id"] not in ("b0", "b1"):
            device.pop("mon", None)
    events.append(
        {
            "time": START + WIDTH,
            "type": "signal_step",
            "signal_id": "command",
            "value": 0.0,
        }
    )
    meta.update(
        pulse_width_s=WIDTH,
        amplitude_V=AMPLITUDE,
        sections=0 if kind == "distributed" else int(kind.split("_")[1]),
    )
    return case, state, events, meta


def simulate(args):
    for precision, rtol, atol in (("nominal", 1e-7, 1e-9), ("fine", 1e-8, 1e-10)):
        options = SimpleNamespace(
            results=args.output / precision,
            exe=args.exe,
            trials=1,
            rtol=rtol,
            atol=atol,
        )
        for kind in args.models:
            case, state, events, meta = make_case(
                kind, args.fits, args.output / "coefficients"
            )
            run_case(options, kind, case, state, events, meta, HORIZON, SPACING)
            if read(options.results / kind / "statistics.json")[-1]["returncode"]:
                raise RuntimeError(
                    f"Simulation failed: {precision}/{kind}/trial1/run.log"
                )


def distributed_transfer(s, yc, h):
    """Independent scalar two-port solution with 100/600 ohm terminations."""
    zc = 1 / yc
    denominator = (1 + 100 / 600) * (1 + h * h) + (zc / 600 + 100 / zc) * (1 - h * h)
    receiving = 2 * h / denominator
    sending = ((1 + h * h) + zc / 600 * (1 - h * h)) / denominator
    reflected = (yc + 1 / 600) * 2 / denominator
    return {
        "sending_V": sending,
        "receiving_V": receiving,
        "i_ref1_A": reflected,
        "i_inc2_A": h * reflected,
    }


def reference(fits, coefficients, spacing, kinds):
    """Invert circuit transfer functions on a Laplace contour, without EMT.

    A long transform window suppresses periodic wraparound. Repeating at half
    the transform spacing checks convergence independently of IDA/history.
    """
    samples = 2 ** int(np.ceil(np.log2(0.016 / spacing)))
    sigma = 24 / (samples * spacing)
    s = sigma + 2j * np.pi * np.fft.rfftfreq(samples, spacing)
    excitation = AMPLITUDE * np.exp(-START * s) * (-np.expm1(-WIDTH * s)) / s
    time = np.arange(int(round(HORIZON / spacing)) + 1) * spacing

    def inverse(transfer, feedthrough=0.0, relaxation=0.0):
        # Integrate known high-frequency terms analytically. This removes
        # Fourier ringing from source jumps and slope changes without
        # smoothing the applied pulse or the simulated waveforms.
        remainder = transfer - feedthrough
        if relaxation:
            remainder = remainder - relaxation / (s + relaxation)
        values = np.fft.irfft(excitation * remainder, n=samples)[: len(time)]
        values = values / spacing * np.exp(sigma * time)
        for event, sign in ((START, 1), (START + WIDTH, -1)):
            lag = np.maximum(time - event, 0)
            values += sign * AMPLITUDE * feedthrough * (time >= event)
            if relaxation:
                values += sign * AMPLITUDE * (-np.expm1(-relaxation * lag))
        return values

    yc_fit = read(fits / "yc.json")
    yc = scalar_response(yc_fit, s)
    h = sum(
        scalar_response(mode["H"], s) * np.exp(-s * mode["tau"])
        for mode in read(fits / "h60000.json")["modes"]
    )
    result = {"time_s": time}
    yc_infinity = sum(yc_fit["D"][0])
    feedthrough = {
        "sending_V": 1 / (1 + 100 * yc_infinity),
        "i_ref1_A": 2 * yc_infinity / (1 + 100 * yc_infinity),
    }
    for name, value in distributed_transfer(s, yc, h).items():
        result[f"distributed_{name}"] = inverse(value, feedthrough.get(name, 0.0))
    z = scalar_response(read(coefficients / "zp.json"), s)
    y = scalar_response(read(coefficients / "yp.json"), s)
    continuum = distributed_transfer(
        s, np.sqrt(y / z), np.exp(-LENGTH * np.sqrt(z * y))
    )
    result["continuum_receiving_V"] = inverse(continuum["receiving_V"])
    for kind in kinds:
        if kind == "distributed":
            continue
        n = int(kind.split("_")[1])
        series, shunt = z * LENGTH / n, y * LENGTH / (2 * n)
        admittance = np.full(s.shape, 1 / 600, complex)
        gain = np.ones(s.shape, complex)
        # Work from the load toward the source. This avoids overflow in
        # ABCD matrix powers for many electrically long sections at high f.
        for _ in range(n):
            right = admittance + shunt
            factor = 1 / (1 + series * right)
            gain *= factor
            admittance = shunt + right * factor
        sending = 1 / (1 + 100 * admittance)
        rate = 1 / (
            100 * sum(read(coefficients / "yp.json")["E"][0]) * LENGTH / (2 * n)
        )
        result[f"{kind}_sending_V"] = inverse(sending, relaxation=rate)
        result[f"{kind}_receiving_V"] = inverse(sending * gain)
    return result


def waveform(folder):
    values = np.genfromtxt(folder / "trial1/response.csv", delimiter=",", names=True)
    if not all(np.all(np.isfinite(values[key])) for key in values.dtype.names):
        raise ValueError(f"Nonfinite waveform: {folder}")
    if np.any(np.diff(values["t"]) < 0) or abs(values["t"][-1] - HORIZON) > 1e-10:
        raise ValueError(f"Incomplete or nonmonotonic waveform: {folder}")
    steps = np.genfromtxt(folder / "trial1/steps.csv", delimiter=",", names=True)
    stats = read(folder / "statistics.json")[-1]
    if len(steps) != stats["steps"]:
        raise ValueError(f"Accepted-step count mismatch: {folder}")
    return values, stats


def rms(values):
    return float(np.sqrt(np.mean(values**2)))


def validate(args):
    references = [
        reference(args.fits, args.output / "coefficients", spacing, MODELS)
        for spacing in (0.1e-6, 0.05e-6)
    ]
    time = np.arange(int(round(HORIZON / SPACING)) + 1) * SPACING
    keep = time >= START
    # Event rows have left/right values, whereas a transform converges to
    # their midpoint. Exclude a 1 us neighborhood only for numerical checks.
    away = keep & (abs(time - START) > 1e-6) & (abs(time - START - WIDTH) > 1e-6)
    oracle = {
        key: np.interp(time, references[-1]["time_s"], value)
        for key, value in references[-1].items()
        if key != "time_s"
    }
    convergence = {
        key: rms(
            (np.interp(time, references[0]["time_s"], references[0][key]) - value)[away]
        )
        for key, value in oracle.items()
    }
    ref = oracle["distributed_receiving_V"]
    summary = {
        "reference_refinement_rms": convergence,
        "continuum_fit_difference_rms_V": rms(
            (oracle["continuum_receiving_V"] - ref)[keep]
        ),
        "metric_window_s": [START, HORIZON],
        "models": {},
    }
    failures = []
    curves = {"time_s": time, **oracle}
    for kind in MODELS:
        nominal, stats = waveform(args.output / "nominal" / kind)
        fine, fine_stats = waveform(args.output / "fine" / kind)
        signals = {"sending_V": "Bus_b0_va", "receiving_V": "Bus_b1_va"}
        if kind == "distributed":
            signals.update(
                i_ref1_A="LineDistributed_line_i_ref1a",
                i_inc2_A="LineDistributed_line_i_inc2a",
            )
        metrics = {
            "sections": 0 if kind == "distributed" else int(kind.split("_")[1]),
            "nominal_statistics": stats,
            "fine_statistics": fine_stats,
            "numerical_checks": {},
        }
        for name, column in signals.items():
            v = np.interp(time, nominal["t"], nominal[column])
            vf = np.interp(time, fine["t"], fine[column])
            target = oracle[f"{kind}_{name}"]
            curves[f"emt_{kind}_{name}"] = vf
            metrics["numerical_checks"][name] = {
                "nominal_reference_rms": rms((v - target)[away]),
                "fine_reference_rms": rms((vf - target)[away]),
                "fine_reference_max": float(max(abs(vf - target)[away])),
                "tolerance_refinement_rms": rms((v - vf)[away]),
                "at_source_events_max": float(max(abs(vf - target)[keep & ~away])),
            }
            limit = 1e-3 if name.endswith("_V") else 1e-5
            check = metrics["numerical_checks"][name]
            if any(
                check[key] > limit
                for key in ("fine_reference_rms", "tolerance_refinement_rms")
            ):
                failures.append(f"{kind}/{name}: numerical RMS exceeds {limit}")
            if convergence[f"{kind}_{name}"] > limit:
                failures.append(
                    f"{kind}/{name}: inverse-transform refinement exceeds {limit}"
                )
        receiving = curves[f"emt_{kind}_receiving_V"]
        peak = np.flatnonzero(keep)[np.argmax(receiving[keep])]
        metrics.update(
            peak_V=float(receiving[peak]),
            peak_time_after_launch_us=float((time[peak] - START) * 1e6),
            waveform_error_rms_V=rms((receiving - ref)[keep]),
            waveform_error_percent=100 * rms((receiving - ref)[keep]) / rms(ref[keep]),
        )
        metrics["phase_symmetry_max_V"] = float(
            max(
                np.max(abs(fine[f"Bus_{bus}_va"] - fine[f"Bus_{bus}_v{phase}"]))
                for bus in ("b0", "b1")
                for phase in "bc"
            )
        )
        if metrics["phase_symmetry_max_V"] > 1e-4:
            failures.append(f"{kind}: phase symmetry error exceeds 0.1 mV")
        summary["models"][kind] = metrics
    errors = [summary["models"][f"pi_{n}"]["waveform_error_rms_V"] for n in SECTIONS]
    if np.any(np.diff(errors) >= 0):
        failures.append(
            "Receiving-waveform RMS error did not decrease with spatial refinement"
        )
    summary["checks"] = {
        "passed": not failures,
        "failures": failures,
        "voltage_numerical_rms_limit_V": 1e-3,
        "current_numerical_rms_limit_A": 1e-5,
        "phase_symmetry_limit_V": 1e-4,
    }
    write(args.output / "validation.json", summary)
    if failures:
        raise ValueError("; ".join(failures))
    np.savez(args.output / "waveforms.npz", **curves)
    with (args.output / "waveforms.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(curves)
        writer.writerows(zip(*curves.values()))
    return summary


def main(args):
    args.output = args.output.resolve()
    args.fits = args.fits.resolve()
    args.exe = args.exe.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.stage in ("all", "run"):
        quality = fit_lumped(
            args.fits, args.output / "coefficients", args.gridworkbench
        )
        print("Series fit quality:", quality, flush=True)
        files = [
            GEOMETRY,
            GEOMETRY.parent / "north-american.catalog.json",
            Path(__file__),
            Path(__file__).with_name("pulse_plot.py"),
            Path(__file__).with_name("run.py"),
            args.exe,
            *(
                args.fits / name
                for name in ("parameters.npz", "yc.json", "h60000.json", "summary.json")
            ),
        ]
        if args.gridworkbench is not None:
            files.extend(
                sorted(
                    (args.gridworkbench / "gridworkbench/emt/parameters").glob("*.py")
                )
            )
        write(
            args.output / "provenance.json",
            {
                "files_sha256": {str(p): digest(p) for p in files},
                "geometry": str(GEOMETRY),
                "source_fits": str(args.fits),
                "python_version": platform.python_version(),
                "numpy_version": np.__version__,
                "scipy_version": scipy.__version__,
                "fit_band_hz": read(args.fits / "summary.json")["fit_band_hz"],
            },
        )
        simulate(args)
    if args.stage in ("all", "plot", "diagrams"):
        if tuple(args.models) != MODELS:
            raise ValueError("The pulse report requires all five models")
        from pulse_plot import diagrams, plots

        if args.stage == "diagrams":
            diagrams(args.output, read(args.output / "validation.json"))
        else:
            plots(args.output, validate(args))
        print(f"Report: {args.output / 'index.html'}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fits", type=Path, default=ROOT / "local/emt-distributed/fits-final"
    )
    parser.add_argument(
        "--output", type=Path, default=ROOT / "local/emt-distributed/pulse"
    )
    parser.add_argument(
        "--exe", type=Path, default=ROOT / "build/application/EMT/EMTDynamicSimulation"
    )
    parser.add_argument("--stage", choices=("all", "run", "plot", "diagrams"), default="all")
    parser.add_argument(
        "--gridworkbench",
        type=Path,
        help="Also check series fits against fresh geometry samples on an independent frequency grid",
    )
    parser.add_argument("--models", nargs="+", choices=MODELS, default=MODELS)
    main(parser.parse_args())
