#!/usr/bin/env python3
"""Fit passive, coupled EMT line operators to the OpenLine parameter sweep."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
OPENLINE_REVISION = "56c405add9f082fa0816338978f25a8242c2cc15"
FREQUENCY_SCALE = 2 * np.pi * 30000.0


def recompute(openline, output):
    """Build a small driver in /tmp, leaving the supplied repository untouched."""
    revision = subprocess.check_output(
        ["git", "-C", str(openline), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != OPENLINE_REVISION:
        raise ValueError(f"Expected OpenLine {OPENLINE_REVISION}, found {revision}")
    changes = subprocess.check_output(
        ["git", "-C", str(openline), "status", "--porcelain", "--untracked-files=no"],
        text=True,
    )
    if changes:
        raise ValueError("OpenLine has tracked local changes; use a clean source checkout")
    with tempfile.TemporaryDirectory(prefix="gridkit-line-driver-") as directory:
        manifest = Path(directory) / "Cargo.toml"
        manifest.write_text(
            '[package]\nname="gridkit-line-driver"\nversion="0.1.0"\nedition="2021"\n'
            '[dependencies]\n'
            f'openline={{path={json.dumps(str(openline / "crates/openline"))}}}\n'
            f'openline-compute={{path={json.dumps(str(openline / "crates/openline-compute"))}}}\n'
            '[[bin]]\nname="gridkit-line-sweep"\n'
            f'path={json.dumps(str(HERE / "openline_sweep.rs"))}\n'
        )
        shutil.copyfile(HERE / "Cargo.lock", Path(directory) / "Cargo.lock")
        env = os.environ.copy()
        env.setdefault("CARGO_HOME", "/tmp/gridkit-cargo")
        env.setdefault("CARGO_TARGET_DIR", "/tmp/gridkit-line-target")
        subprocess.run(
            ["cargo", "run", "--locked", "--manifest-path", str(manifest), "--bin",
             "gridkit-line-sweep", "--", str(output)],
            env=env, check=True,
        )


def read_sweep(path):
    with path.open() as stream:
        records = list(csv.DictReader(stream))
    frequencies = np.array(sorted({float(row["frequency_hz"]) for row in records}))
    indices = {frequency: index for index, frequency in enumerate(frequencies)}
    values = {name: np.zeros((len(frequencies), 3, 3), complex) for name in ("Z", "Y")}
    for row in records:
        # CsvWriter exports ohm/km and S/km; LineLumped dx is supplied in meters.
        value = complex(float(row["real"]), float(row["imag"])) / 1000.0
        values[row["quantity"]][indices[float(row["frequency_hz"])],
                                int(row["row"]), int(row["col"])] = value
    return frequencies, values


def basis(frequencies, rates):
    s = 2j * np.pi * np.asarray(frequencies)
    return np.column_stack((np.ones(len(s)), s / FREQUENCY_SCALE,
                            s[:, None] / (s[:, None] + rates)))


def fit_series(frequencies, data):
    """Shared stable poles, real symmetric coefficients, relative entry weighting."""
    training = np.arange(0, len(frequencies), 2)
    weights = 1.0 / np.abs(data)

    def solve_coefficients(log_rates):
        matrix = basis(frequencies, np.exp(log_rates))
        coefficients = np.empty((len(log_rates) + 2, 3, 3))
        for row in range(3):
            for col in range(row, 3):
                weighted = matrix[training] * weights[training, row, col, None]
                target = data[training, row, col] * weights[training, row, col]
                answer = np.linalg.lstsq(
                    np.concatenate((weighted.real, weighted.imag)),
                    np.concatenate((target.real, target.imag)), rcond=None,
                )[0]
                coefficients[:, row, col] = coefficients[:, col, row] = answer
        prediction = np.einsum("fk,kij->fij", matrix, coefficients)
        return coefficients, prediction

    def residual(log_rates):
        _, prediction = solve_coefficients(log_rates)
        error = ((prediction - data) * weights)[training]
        return np.concatenate((error.real.ravel(), error.imag.ravel()))

    solution = least_squares(
        residual, np.log(2 * np.pi * np.geomspace(0.5, 1e5, 8)),
        max_nfev=140, ftol=1e-9, xtol=1e-9, gtol=1e-9,
    )
    if not solution.success:
        raise RuntimeError(f"Pole optimization did not converge: {solution.message}")
    rates = np.exp(solution.x)
    coefficients, _ = solve_coefficients(solution.x)
    before_projection = coefficients.copy()
    # Z(s)=R0+sLinf+sum A_k*s/(s+a_k) is positive real when a_k>0 and
    # R0,Linf,A_k are positive semidefinite. Project small negative eigenvalues
    # from the unconstrained fit, then validate the resulting approximation.
    for index, coefficient in enumerate(coefficients):
        eigenvalues, vectors = np.linalg.eigh(coefficient)
        coefficients[index] = (vectors * np.maximum(eigenvalues, 0.0)) @ vectors.T
    prediction = np.einsum("fk,kij->fij", basis(frequencies, rates), coefficients)
    relative = np.abs((prediction - data) / data)
    if relative[1::2].max() > 0.003:
        raise ValueError("Held-out maximum entrywise relative error exceeds 0.3%")
    metrics = {
        "training_points": len(training),
        "held_out_points": len(frequencies[1::2]),
        "max_held_out_entry_relative_error": float(relative[1::2].max()),
        "rms_held_out_entry_relative_error": float(np.sqrt(np.mean(relative[1::2] ** 2))),
        "max_held_out_matrix_relative_error": float(np.max(
            np.linalg.norm((prediction - data)[1::2], axis=(1, 2))
            / np.linalg.norm(data[1::2], axis=(1, 2)))),
        "pole_rates_per_s": rates.tolist(),
        "pole_fit_function_evaluations": solution.nfev,
        "passivity": {
            "realization": "R0 + s*Linf + sum A_k*s/(s+a_k), all matrices symmetric PSD",
            "R0_eigenvalues_ohm_per_m": np.linalg.eigvalsh(coefficients[0]).tolist(),
            "Linf_eigenvalues_H_per_m": np.linalg.eigvalsh(coefficients[1] / FREQUENCY_SCALE).tolist(),
            "A_eigenvalues_ohm_per_m": np.linalg.eigvalsh(coefficients[2:]).tolist(),
            "max_projection_entry_change_ohm_per_m": float(np.max(np.abs(
                coefficients[2:] - before_projection[2:]))),
        },
    }
    return rates, coefficients, prediction, metrics


def complex_array(array):
    return np.stack((np.real(array), np.imag(array)), axis=-1).tolist()


def vectorfit_model(direct, proportional, rates=(), residues=()):
    return {"rows": 3, "cols": 3, "D": np.asarray(direct).tolist(),
            "E": np.asarray(proportional).tolist(),
            "poles": [[-float(rate), 0.0] for rate in rates],
            "residues": complex_array(np.asarray(residues)) if len(residues) else []}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recompute", action="store_true", help="Run OpenLine before fitting")
    parser.add_argument("--openline", type=Path, default=Path.home() / "openline")
    args = parser.parse_args()
    source = HERE / "openline.csv"
    if args.recompute:
        recompute(args.openline.resolve(), source)
    frequencies, raw = read_sweep(source)
    rates, coefficients, fitted_z, metrics = fit_series(frequencies, raw["Z"])
    capacitance = (raw["Y"] / (2j * np.pi * frequencies[:, None, None])).real.mean(axis=0)
    fitted_y = 2j * np.pi * frequencies[:, None, None] * capacitance
    shunt_error = float(np.max(np.abs((fitted_y - raw["Y"]) / raw["Y"])))
    if shunt_error > 1e-10 or np.linalg.eigvalsh(capacitance).min() <= 0:
        raise ValueError("Expected exact positive-definite Maxwell capacitance matrix")
    # Independent electrostatic check: invert the image-method potential matrix.
    coordinates = np.array([[-1.2, 10.5], [0.0, 11.0], [1.2, 10.5]])
    potential = np.empty((3, 3))
    for row in range(3):
        for col in range(3):
            x1, h1 = coordinates[row]
            x2, h2 = coordinates[col]
            distance = 0.009144 if row == col else np.hypot(x1 - x2, h1 - h2)
            image_distance = np.hypot(x1 - x2, h1 + h2)
            potential[row, col] = np.log(image_distance / distance) / (2 * np.pi * 8.8541878128e-12)
    electrostatic_error = float(np.max(np.abs((np.linalg.inv(potential) - capacitance) / capacitance)))
    if electrostatic_error > 1e-8:
        raise ValueError("OpenLine capacitance disagrees with independent image-method calculation")
    model = {
        "Zp": vectorfit_model(
            coefficients[0] + coefficients[2:].sum(axis=0),
            coefficients[1] / FREQUENCY_SCALE, rates,
            -rates[:, None, None] * coefficients[2:]),
        "Yp": vectorfit_model(np.zeros((3, 3)), capacitance),
    }
    (HERE / "model.json").write_text(json.dumps(model, indent=2) + "\n")
    # Numerically cross-check the exported standard VectorFit coefficients.
    exported = np.asarray(model["Zp"]["D"]) + (
        2j * np.pi * frequencies[:, None, None] * np.asarray(model["Zp"]["E"]))
    for rate, residue in zip(rates, -rates[:, None, None] * coefficients[2:]):
        exported += residue / (2j * np.pi * frequencies[:, None, None] + rate)
    if not np.allclose(exported, fitted_z, rtol=1e-10, atol=1e-13):
        raise ValueError("Exported pole-residue realization differs from fitted model")
    scan_frequencies = np.geomspace(1e-4, 1e9, 20001)
    scan_z = np.einsum("fk,kij->fij", basis(scan_frequencies, rates), coefficients)
    minimum_resistance = float(np.linalg.eigvalsh(scan_z.real).min())
    if minimum_resistance <= 0:
        raise ValueError("Exported impedance fails numerical passivity check")
    # Compare a 600 m pi section with exact uniform-line telegrapher equations,
    # keeping the same fitted per-meter Z and Y in both calculations.
    pi_errors = {}
    for frequency in (60, 900, 2700, 9000, 30000):
        series = np.einsum("fk,kij->fij", basis([frequency], rates), coefficients)[0]
        shunt = 2j * np.pi * frequency * capacitance
        transfer = expm(np.block([[np.zeros((3, 3)), series],
                                  [shunt, np.zeros((3, 3))]]) * 600)
        aa, bb = transfer[:3, :3], transfer[:3, 3:]
        cc, dd = transfer[3:, :3], transfer[3:, 3:]
        inverse = np.linalg.inv(bb)
        exact = np.block([[dd @ inverse, cc - dd @ inverse @ aa],
                          [-inverse, inverse @ aa]])
        branch = np.linalg.inv(series * 600)
        lumped = np.block([[branch + shunt * 300, -branch],
                           [-branch, branch + shunt * 300]])
        pi_errors[str(frequency)] = float(np.linalg.norm(lumped - exact) / np.linalg.norm(exact))
    metrics.update({
        "source": {"repository": "https://github.com/lukelowry/openline",
                   "revision": OPENLINE_REVISION,
                   "sweep_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                   "driver_sha256": hashlib.sha256((HERE / "openline_sweep.rs").read_bytes()).hexdigest()},
        "frequency_band_Hz": [float(frequencies[0]), float(frequencies[-1])],
        "physical_parameters": {
            "system_voltage_kV": 13.8,
            "conductor": "ACSR_Linnet, 336.4 kcmil, 26 aluminum / 7 steel strands",
            "conductor_diameter_m": 0.018288,
            "core_diameter_m": 0.00673608,
            "catalog_GMR_m": 0.007415784,
            "catalog_DC_resistance_ohm_per_m_at_20C": 0.0001660105,
            "conductor_temperature_C": 20,
            "earth_resistivity_ohm_m": 100,
            "phase_coordinates_m": coordinates.tolist(),
            "transposition": "Untransposed",
            "formulations": ["SchelkunoffInternal", "PerfectEarthExternal",
                             "CarsonEarthReturn", "MaxwellPerfectEarth"],
        },
        "shunt_max_entry_relative_error": shunt_error,
        "independent_capacitance_max_entry_relative_error": electrostatic_error,
        "capacitance_F_per_m": capacitance.tolist(),
        "sampled_passivity_band_Hz": [1e-4, 1e9],
        "sampled_passivity_points": len(scan_frequencies),
        "minimum_sampled_real_Z_eigenvalue_ohm_per_m": minimum_resistance,
        "pi_600m_relative_nodal_admittance_Frobenius_error_by_Hz": pi_errors,
        "definition": "Entry relative error = abs(fit-OpenLine)/abs(OpenLine); odd-index points held out",
    })
    (HERE / "fit_validation.json").write_text(json.dumps(metrics, indent=2) + "\n")
    # Wide-form input for plots; the final case uses model.json directly.
    columns = ["freq_hz"]
    arrays = [frequencies]
    for name, values in (("raw_Z", raw["Z"]), ("fit_Z", fitted_z),
                         ("raw_Y", raw["Y"]), ("fit_Y", fitted_y)):
        for row in range(3):
            for col in range(3):
                for part, operation in (("re", np.real), ("im", np.imag)):
                    columns.append(f"{name}{row}{col}_{part}")
                    arrays.append(operation(values[:, row, col]))
    np.savetxt(HERE / "response.csv", np.column_stack(arrays), delimiter=",",
               header=",".join(columns), comments="", fmt="%.15e")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
