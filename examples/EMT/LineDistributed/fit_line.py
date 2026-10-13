"""Regenerate physical line samples and pole/residue JSON (offline only).

Requires numpy, scipy and scikit-rf==1.13.0. Uses scikit-rf's established vector
fitting implementation; the Network.s array is only a container for normalized
scalar transfer responses, NOT a physical scattering matrix. Passivity is
checked on the assembled line, never inferred from that container.
"""

import json
import hashlib
import os
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/tmp/gridkit-line-distributed-matplotlib")
import numpy as np
import skrf
import scipy
from skrf.vectorFitting import VectorFitting

from line_constants import HERE, T, constants, modal, parameters, wave_delays
from validate import fit, passivity

DATA = HERE / "data/overhead"


def encode(a):
    a = np.asarray(a)
    if np.iscomplexobj(a) and np.any(a.imag != 0):
        return np.stack((a.real, a.imag), axis=-1).tolist()
    return a.real.tolist()


def export(vf, scales, left, right):
    def transform(v):
        # The first two fitted channels are the common and differential modes.
        return left @ np.diag([v[0]*scales[0], v[3]*scales[1], v[3]*scales[1]]) @ right
    poles, residues = [], []
    for i, p in enumerate(vf.poles):
        residue = transform(vf.residues[:, i])
        if p.imag == 0:
            poles.append(float(p.real))
            residues.append(encode(residue.real))
        else:
            poles.extend(([float(p.real), float(p.imag)], [float(p.real), float(-p.imag)]))
            residues.extend((encode(residue), encode(residue.conj())))
    return {"rows": 3, "cols": 3, "D": encode(transform(vf.constant_coeff)),
            "poles": poles, "residues": residues}


def vector_fit(frequency, values, order, left, right):
    scales = np.max(np.abs(values[:, :2]), axis=0)
    response = np.zeros((len(frequency), 2, 2), complex)
    response[:, 0, 0] = values[:, 0]/scales[0]
    response[:, 1, 1] = values[:, 1]/scales[1]
    network = skrf.Network(frequency=skrf.Frequency.from_f(frequency, unit="hz"), s=response)
    vf = VectorFitting(network)
    vf.max_iterations = 100
    with warnings.catch_warnings():
        # This container is not an S-parameter model; its passivity is irrelevant.
        warnings.filterwarnings("ignore", message="The fitted network is passive.*", category=UserWarning)
        vf.vector_fit(n_poles_real=4, n_poles_cmplx=(order-4)//2, init_pole_spacing="log",
                      fit_constant=True, fit_proportional=False, enforce_dc=False)
    # Refit residues with a genuine DC equality: eliminate D from every row,
    # F(s)-F(0) = sum r_k [1/(s-p_k) + 1/p_k]. Column scaling avoids disparate
    # pole magnitudes. Relative weighting preserves the low-frequency response.
    s = 2j*np.pi*frequency
    basis = []
    for p in vf.poles:
        a, b = abs(p)/(s-p), abs(p)/(s-p.conjugate())
        basis.extend([a] if p.imag == 0 else [a+b, 1j*(a-b)])
    basis = np.array(basis).T
    dc_basis = basis[0].real.copy()
    basis -= dc_basis
    for mode, channel in enumerate((0, 3)):
        target = values[:, mode]/scales[mode]
        weight = 1/np.maximum(abs(target), 1e-3)
        a = basis[1:]*weight[1:, None]
        b = (target[1:]-target[0])*weight[1:]
        norm = np.linalg.norm(a, axis=0)
        a /= norm
        x = np.linalg.lstsq(np.vstack((a.real, a.imag)), np.r_[b.real, b.imag], rcond=None)[0]/norm
        vf.constant_coeff[channel] = target[0].real - dc_basis@x
        index = 0
        for i, p in enumerate(vf.poles):
            vf.residues[channel, i] = abs(p)*(x[index] if p.imag == 0 else x[index]+1j*x[index+1])
            index += 1 if p.imag == 0 else 2
    return export(vf, scales, left, right)


def error(actual, expected, floor):
    return float(np.max(np.abs(actual-expected)/np.maximum(np.abs(expected), floor)))


def main():
    report_path = HERE / "results/fit-report.json"
    report_path.parent.mkdir(exist_ok=True)
    report_path.unlink(missing_ok=True)
    d = parameters()
    # Low-frequency guard points and DC are needed to recover series resistance.
    f = np.unique(np.r_[0, np.geomspace(1e-4, 1e6, 701), 60])
    check_f = np.unique(np.r_[0, np.geomspace(1.013e-4, 1e6, 1801), 60])
    s, check_s = 2j*np.pi*f, 2j*np.pi*check_f
    z, y = constants(s)
    columns = np.column_stack((f, *[v for a in (z, y) for m in (0, 1) for v in (a[:, m].real, a[:, m].imag)]))
    np.savetxt(DATA / "zy.csv", columns, delimiter=",", fmt="%.16e", comments="",
               header="f_Hz,z0_real_ohm_m,z0_imag_ohm_m,z1_real_ohm_m,z1_imag_ohm_m,y0_real_S_m,y0_imag_S_m,y1_real_S_m,y1_imag_S_m")
    yc = np.sqrt(y/z)
    expected_yc, _ = modal(check_s, 0)
    identity = {"rows": 3, "cols": 3, "D": T.T.tolist(), "poles": [], "residues": []}
    report = {"fitter": "scikit-rf " + skrf.__version__, "numpy": np.__version__, "scipy": scipy.__version__, "fit_band_Hz": [1e-4, 1e6],
              "validation_points": len(check_f), "physical_assumptions": d["assumptions"]}
    for order in (8, 12, 16, 20, 24, 28, 32, 36, 40, 48):
        yc_fit = vector_fit(f, yc, order, T, T.T)
        actual_yc = np.diagonal(T.T @ fit(yc_fit, check_s) @ T, axis1=1, axis2=2)
        yc_error = error(actual_yc, expected_yc, 1e-6)
        print("Yc candidate", order, yc_error, flush=True)
        if yc_error < .0002:
            break
    else:
        raise RuntimeError(f"Yc fit did not meet error budget: {yc_error}")
    print("Yc", len(yc_fit["poles"]), yc_error, flush=True)
    (DATA / "yc.fit.json").write_text(json.dumps(yc_fit, indent=2)+"\n")
    report["Yc"] = {"poles": len(yc_fit["poles"]), "max_scaled_error": yc_error}
    for name, length in d["lengths_m"].items():
        tau = wave_delays(length)
        _, propagation = modal(s, length)
        attenuation = propagation*np.exp(s[:, None]*tau[None])
        _, expected = modal(check_s, length)
        for order in (8, 12, 16, 20, 24, 28, 32, 36, 40, 48):
            output = vector_fit(f, attenuation, order, T, np.eye(3))
            h = {"input": identity, "delays": tau.tolist(), "output": output}
            line = {"N": 3, "K": 3, "conductors": [1, 2, 3], "Yc": yc_fit, "H": h}
            actual = np.diagonal(T.T @ fit(output, check_s), axis1=1, axis2=2)*np.exp(-check_s[:, None]*tau[None])
            h_error = error(actual, expected, 1e-3)
            passive = passivity({"line": line})
            print(name, len(output["poles"]), h_error, passive["violations_beyond_roundoff"], flush=True)
            if h_error < .002 and passive["violations_beyond_roundoff"] == 0:
                break
        else:
            raise RuntimeError(f"{name} fit did not meet error/passivity budget")
        # Scalar modal terminal eigenadmittances avoid the cancellation of a
        # six-terminal nodal conversion near DC and weak-loss resonances.
        guard_s = 2j*np.pi*np.r_[0, np.geomspace(1e-6, 1e8, 6001)]
        guard_y = np.diagonal(T.T@fit(yc_fit, guard_s)@T, axis1=1, axis2=2)
        guard_h = np.diagonal(T.T@fit(output, guard_s), axis1=1, axis2=2)*np.exp(-guard_s[:, None]*tau)
        eigen = np.minimum((guard_y*(1+guard_h)/(1-guard_h)).real,
                           (guard_y*(1-guard_h)/(1+guard_h)).real)
        assert eigen.min() >= -1e-12, (name, "modal terminal passivity", eigen.min())
        raw_y0, raw_h0 = modal([0j], length)
        expected_dc = np.r_[raw_y0*(1+raw_h0)/(1-raw_h0), raw_y0*(1-raw_h0)/(1+raw_h0)]
        actual_dc = np.r_[guard_y[:1]*(1+guard_h[:1])/(1-guard_h[:1]),
                          guard_y[:1]*(1-guard_h[:1])/(1+guard_h[:1])]
        dc_error = error(actual_dc, expected_dc, 1e-15)
        assert dc_error < 1e-7, (name, "DC two-port", dc_error)
        line["Yc"] = "yc.fit.json"
        (DATA / f"{name}.fit.json").write_text(json.dumps(line, indent=2)+"\n")
        report[name] = {"length_m": length, "delays_s": tau.tolist(), "poles": len(output["poles"]),
                        "max_scaled_error": h_error, "passivity": passive,
                        "modal_passivity_guard": {"min_terminal_eigenadmittance_real_S": float(eigen.min()),
                                                   "frequency_Hz": [0, 1e8], "points": len(guard_s)},
                        "dc_terminal_relative_error": dc_error}
    f_ground = np.geomspace(.01, 1e6, 61)
    approximate, _ = constants(2j*np.pi*f_ground)
    exact, _ = constants(2j*np.pi*f_ground, carson=True)
    report["ground_return_check"] = {"method": "Independent Carson integral quadrature", "points": len(f_ground),
                                      "max_relative_modal_Z_error": error(approximate, exact, 1e-15)}
    report["input_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in sorted(DATA.glob("*")) if p.suffix in (".json", ".csv")}
    report_path.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
