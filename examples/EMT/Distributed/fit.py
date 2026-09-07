#!/usr/bin/env python3
"""Fit GridWorkbench overhead-line operators for the distributed-line studies."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def response(fit, omega):
    value = np.broadcast_to(np.array(fit["D"], complex), (len(omega), 3, 3)).copy()
    for pole, residue in zip(fit.get("poles", []), fit.get("residues", [])):
        r = np.array(residue)
        r = r[..., 0] + 1j * r[..., 1]
        value += r / (1j * omega[:, None, None] - complex(*pole))
    return value


def fit_scalar(omega, values, projector, poles, constant):
    """Fit a real causal scalar function, then lift it to a phase projector."""
    import skrf
    from skrf.vectorFitting import VectorFitting

    # The Network object carries generic samples; no S/Y/Z conversion is used.
    network = skrf.Network(
        frequency=skrf.Frequency.from_f(omega / (2 * np.pi), unit="hz"),
        s=values[:, None, None],
        z0=1,
    )
    model = VectorFitting(network)
    model.vector_fit(
        n_poles_real=poles,
        n_poles_cmplx=0,
        init_pole_spacing="log",
        fit_constant=constant,
        fit_proportional=False,
    )
    fit = {
        "rows": 3,
        "cols": 3,
        "D": (model.constant_coeff[0] * projector).tolist(),
        "poles": [],
        "residues": [],
    }
    for pole, residue in zip(model.poles, model.residues[0]):
        if pole.real >= 0:
            raise ValueError("Unstable vector fit")
        for p, r in (
            [(pole, residue), (pole.conjugate(), residue.conjugate())]
            if pole.imag
            else [(pole, residue)]
        ):
            fit["poles"].append([float(p.real), float(p.imag)])
            matrix = r * projector
            fit["residues"].append(
                np.stack((matrix.real, matrix.imag), axis=-1).tolist()
            )
    error = np.linalg.norm(
        response(fit, omega) - values[:, None, None] * projector, axis=(1, 2)
    )
    return fit, {
        "absolute_max": float(np.max(error)),
        "absolute_rms": float(np.sqrt(np.mean(error**2))),
    }


def generate(args):
    sys.path.insert(0, str(args.gridworkbench))
    from gridworkbench.emt.parameters import load, sweep

    args.output.mkdir(parents=True, exist_ok=True)
    omega = 2 * np.pi * np.unique(np.r_[np.geomspace(0.01, 1e7, 701), 60])
    line = load(args.geometry)
    data = sweep(line, omega)

    # Cyclic transposition of the physical GridWorkbench matrices gives a
    # reciprocal, uniform transposed line with constant real modal projectors.
    def transpose(value):
        return sum(np.roll(np.roll(value, k, axis=1), k, axis=2) for k in range(3)) / 3

    Z = transpose(data.R + 1j * omega[:, None, None] * data.L)
    Y = transpose(data.G + 1j * omega[:, None, None] * data.C)
    P0 = np.ones((3, 3)) / 3
    projectors = [np.eye(3) - P0, P0]
    zm = np.column_stack((Z[:, 0, 0] - Z[:, 0, 1], Z[:, 0, 0] + 2 * Z[:, 0, 1]))
    ym = np.column_stack((Y[:, 0, 0] - Y[:, 0, 1], Y[:, 0, 0] + 2 * Y[:, 0, 1]))
    lam = np.sqrt(zm * ym)
    characteristic = sum(
        np.sqrt(ym[:, m] / zm[:, m])[:, None, None] * projectors[m] for m in range(2)
    )
    modal_yc = [
        fit_scalar(omega, np.sqrt(ym[:, m] / zm[:, m]), projectors[m], 18, True)
        for m in range(2)
    ]
    yc = {
        "rows": 3,
        "cols": 3,
        "D": sum(np.array(f[0]["D"]) for f in modal_yc).tolist(),
        "poles": sum((f[0]["poles"] for f in modal_yc), []),
        "residues": sum((f[0]["residues"] for f in modal_yc), []),
    }
    metrics = [f[1] for f in modal_yc]
    (args.output / "yc.json").write_text(json.dumps(yc, indent=2) + "\n")
    np.savez(
        args.output / "parameters.npz",
        omega=omega,
        Z=Z,
        Y=Y,
        lam=lam,
        yc=characteristic,
    )
    summary = {
        "source": "abirchfield/GridWorkbench, cyclically transposed 345kv-horizontal geometry",
        "geometry": str(args.geometry),
        "fit_band_hz": [0.01, 1e7],
        "fitter": "scikit-rf real-coefficient vector fitting",
        "yc": metrics,
        "lengths": {},
    }
    # Remove only the vacuum flight time; retain the dispersive modal lag in H.
    tau_per_m = np.full(2, 1 / 299792458.0)
    check_omega = 2 * np.pi * np.r_[0, np.geomspace(1e-5, 1e10, 10001), 60]
    check_yc = response(yc, check_omega)
    for length in args.lengths:
        modes, errors = [], []
        for m in range(2):
            tau = float(length * tau_per_m[m])
            values = np.exp(-length * lam[:, m] + 1j * omega * tau)
            fit, error = fit_scalar(omega, values, projectors[m], args.poles, False)
            modes.append({"tau": tau, "H": fit})
            errors.append(error)
        propagation = {"K": 3, "modes": modes}
        check_h = sum(
            response(mode["H"], check_omega)
            * np.exp(-1j * check_omega * mode["tau"])[:, None, None]
            for mode in modes
        )

        def minimum_conductance(scale):
            result = float("inf")
            for projector in projectors:
                rank = np.trace(projector)
                y = np.einsum("fij,ji->f", check_yc, projector) / rank
                h = np.einsum("fij,ji->f", check_h, projector) / rank * scale
                result = min(
                    result,
                    np.min((y * (1 - h) / (1 + h)).real),
                    np.min((y * (1 + h) / (1 - h)).real),
                )
            return float(result)

        before = minimum_conductance(1)
        scale = 1.0
        if before < 0:
            lo, hi = 0.0, 1.0
            for _ in range(40):
                mid = (lo + hi) / 2
                if minimum_conductance(mid) >= 1e-12:
                    lo = mid
                else:
                    hi = mid
            scale = lo
            if scale < 0.999:
                raise ValueError(f"Fit needs excessive attenuation: {scale}")
            for mode in modes:
                for key in ["D", "residues"]:
                    mode["H"][key] = (np.array(mode["H"][key]) * scale).tolist()
        prediction = sum(
            response(mode["H"], omega)
            * np.exp(-1j * omega * mode["tau"])[:, None, None]
            for mode in modes
        )
        target = sum(
            np.exp(-length * lam[:, m])[:, None, None] * projectors[m] for m in range(2)
        )
        relative = np.linalg.norm(prediction - target, axis=(1, 2)) / np.linalg.norm(
            target, axis=(1, 2)
        )
        metrics_h = {
            "delays_s": [m["tau"] for m in modes],
            "modal_fit_errors": errors,
            "relative_max": float(np.max(relative)),
            "relative_rms": float(np.sqrt(np.mean(relative**2))),
            "attenuation_scale": scale,
            "sampled_min_terminal_conductance_before_S": before,
            "sampled_min_terminal_conductance_after_S": minimum_conductance(scale),
            "sampled_passivity_band_hz": [0, 1e10],
            "passivity_samples": len(check_omega),
        }
        print(length, metrics_h, flush=True)
        (args.output / f"h{length:g}.json").write_text(
            json.dumps(propagation, indent=2) + "\n"
        )
        summary["lengths"][str(length)] = metrics_h
    at60 = int(np.argmin(abs(omega - 2 * np.pi * 60)))
    summary["pi_per_m"] = {
        "R": Z[at60].real.tolist(),
        "L": (Z[at60].imag / omega[at60]).tolist(),
        "G": Y[at60].real.tolist(),
        "C": (Y[at60].imag / omega[at60]).tolist(),
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gridworkbench", type=Path, required=True)
    p.add_argument("--poles", type=int, default=32)
    p.add_argument(
        "--geometry",
        type=Path,
        default=Path(__file__).resolve().parents[3]
        / "cases/EMT/Distributed/overhead.line.json",
    )
    p.add_argument("--output", type=Path, default=Path("local/emt-distributed/fits"))
    p.add_argument(
        "--lengths", type=float, nargs="+", default=[5000, 20000, 40000, 60000]
    )
    generate(p.parse_args())
