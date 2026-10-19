"""Line data for both simulators from the geometry documents in lines/.

For every ``lines/<name>.line.json`` this writes ``lines/<name>/``:

- ``constants.csv``: transposed modal z and y per metre over frequency (GridWorkbench).
- ``sequence.json``: 60 Hz sequence R, L, C, surge impedances and lossless delays.
- ``tower.json``: the PSCAD universal-tower description of the same geometry.
- ``fd.fit.json``: GridKit LineDistributed input, vector-fitted Yc and propagation.
- ``bergeron.fit.json``: GridKit pole-free (lossless Bergeron) input from 60 Hz L and C.
- ``fit-report.json`` and ``fits.png``: fitting orders, held-out errors, passivity.

Run with the GridWorkbench environment: ``python make_lines.py [name ...]``.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parent
sys.path.insert(0, str(HERE / "lib"))
sys.path.insert(0, str(EXAMPLE))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from linedata import T, compute, tower  # noqa: E402
from validate import fit as evaluate, passivity  # noqa: E402
from vectorfit import scaled_error, vector_fit  # noqa: E402

LINES = HERE / "lines"
FREQUENCY = np.unique(np.r_[0.0, np.geomspace(1e-3, 1e6, 601), 60.0])  # 0 Hz pins the fits at DC
CHECK = np.unique(np.r_[np.geomspace(1.07e-3, 0.97e6, 1201), 60.0])
YC_ORDERS = (6, 8, 10, 12, 14, 16, 20, 24)
H_ORDERS = (8, 12, 16, 20, 24, 28, 32, 36, 40, 48, 56)
YC_TARGET = 2e-4   # held-out relative error of the modal characteristic admittances
H_TARGET = 2e-3    # held-out error of the modal propagation functions, magnitude floor 1e-3
COLORS = ("#0072B2", "#D55E00")


def encode(a):
    a = np.asarray(a)
    if np.iscomplexobj(a) and np.any(a.imag != 0):
        return np.stack((a.real, a.imag), axis=-1).tolist()
    return a.real.tolist()


def export(fit, left, right):
    """A GridKit VectorFit of the modal fit (ground, aerial) placed on the diagonal
    (ground, aerial, aerial) between the factors ``left`` and ``right``."""
    def transform(v):
        return left @ np.diag([v[0], v[1], v[1]]) @ right

    poles, residues = [], []
    q = 0
    while q < len(fit.poles):
        p = fit.poles[q]
        residue = transform(fit.residues[:, q])
        if p.imag == 0:
            poles.append(float(p.real))
            residues.append(encode(residue.real))
            q += 1
        else:
            poles.extend(([float(p.real), float(p.imag)], [float(p.real), float(-p.imag)]))
            residues.extend((encode(residue), encode(residue.conj())))
            q += 2
    return {"rows": 3, "cols": 3, "D": encode(transform(fit.d)), "poles": poles, "residues": residues}


def modal_of(matrix_response):
    """Modal diagonal (ground, aerial) of evaluated 3x3 phase responses T diag T^T."""
    return np.diagonal(T.T @ matrix_response @ T, axis1=1, axis2=2)[:, :2]


def modal_output(matrix_response):
    """Modal diagonal (ground, aerial) of an evaluated output factor T diag."""
    return np.diagonal(T.T @ matrix_response, axis1=1, axis2=2)[:, :2]


def fit_admittance(data):
    yc = data.characteristic_admittance().T            # (2, frequencies)
    expected = data.characteristic_admittance(2j * np.pi * CHECK).T
    best = None
    for order in YC_ORDERS:
        fit = vector_fit(FREQUENCY, yc, order, iterations=30)
        exported = export(fit, T, T.T)
        actual = modal_of(evaluate(exported, 2j * np.pi * CHECK)).T
        error = scaled_error(actual, expected, 1e-6)
        print(f"  Yc order {order:2d}: held-out relative error {error:.2e}", flush=True)
        if best is None or error < best[1]:
            best = (exported, error, order)
        if error < YC_TARGET:
            break
    return best


def fit_propagation(data, yc_fit):
    delays = data.wave_delays
    attenuation = data.attenuation().T                 # (2, frequencies)
    expected = data.propagation(2j * np.pi * CHECK).T
    identity = {"rows": 3, "cols": 3, "D": T.T.tolist(), "poles": [], "residues": []}
    best = None
    for order in H_ORDERS:
        fit = vector_fit(FREQUENCY, attenuation, order, iterations=30)
        output = export(fit, T, np.eye(3))
        line = {"N": 3, "K": 3, "conductors": [1, 2, 3], "Yc": yc_fit,
                "H": {"input": identity, "delays": [delays[0], delays[1], delays[1]], "output": output}}
        actual = modal_output(evaluate(output, 2j * np.pi * CHECK)).T * np.exp(-2j * np.pi * CHECK[None, :] * delays[:, None])
        error = scaled_error(actual, expected, 1e-3)
        passive = passivity({"line": line})
        print(f"  H  order {order:2d}: held-out scaled error {error:.2e}, passivity violations {passive['violations_beyond_roundoff']}", flush=True)
        candidate = (line, error, order, passive)
        if best is None or (error < best[1] and passive["violations_beyond_roundoff"] <= best[3]["violations_beyond_roundoff"]):
            best = candidate
        if error < H_TARGET and passive["violations_beyond_roundoff"] == 0:
            break
    return best


def bergeron(data, sequence):
    l = np.array([sequence["L0"], sequence["L1"]])
    c = np.array([sequence["C0"], sequence["C1"]])
    surge = np.sqrt(l / c)
    delays = data.length * np.sqrt(l * c)
    yc = {"rows": 3, "cols": 3, "D": (T @ np.diag([1 / surge[0], 1 / surge[1], 1 / surge[1]]) @ T.T).tolist(),
          "poles": [], "residues": []}
    return {
        "N": 3, "K": 3, "conductors": [1, 2, 3], "Yc": yc,
        "H": {"input": {"rows": 3, "cols": 3, "D": T.T.tolist(), "poles": [], "residues": []},
              "delays": [float(delays[0]), float(delays[1]), float(delays[1])],
              "output": {"rows": 3, "cols": 3, "D": T.tolist(), "poles": [], "residues": []}},
    }, surge, delays


def figure(data, fd, report, destination):
    s = 2j * np.pi * CHECK
    yc_raw = data.characteristic_admittance(s)
    h_raw = data.propagation(s)
    yc_fit = modal_of(evaluate(fd["Yc"], s))
    delays = np.array(fd["H"]["delays"])
    h_fit = modal_output(evaluate(fd["H"]["output"], s)) * np.exp(-s[:, None] * delays[None, :2])
    ac = data.frequency > 0
    z, y, f_ac = data.z[ac], data.y[ac], data.frequency[ac]
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), layout="constrained")
    for m, color, name in zip((0, 1), COLORS, ("Ground mode", "Aerial mode")):
        axes[0, 0].loglog(CHECK, 1e3 * abs(yc_raw[:, m]), color=color, label=name + " (GridWorkbench)")
        axes[0, 0].loglog(CHECK[::25], 1e3 * abs(yc_fit[::25, m]), ".", color=color, ms=4, label=name + " fit")
        axes[0, 1].semilogx(CHECK, abs(h_raw[:, m]), color=color, label=name + " (GridWorkbench)")
        axes[0, 1].semilogx(CHECK[::25], abs(h_fit[::25, m]), ".", color=color, ms=4, label=name + " fit")
        axes[0, 2].semilogx(CHECK, 100 * abs(yc_fit[:, m] - yc_raw[:, m]) / abs(yc_raw[:, m]), color=color, label=name + " Yc")
        axes[0, 2].semilogx(CHECK, 100 * abs(h_fit[:, m] - h_raw[:, m]) / np.maximum(abs(h_raw[:, m]), 1e-3), "--", color=color, label=name + " H")
        axes[1, 0].loglog(f_ac, 1e3 * z[:, m].real, color=color, label=name)
        axes[1, 1].semilogx(f_ac, 1e6 * z[:, m].imag / (2 * np.pi * f_ac), color=color, label=name)
        phase = np.unwrap(np.angle(h_raw[:, m] * np.exp(s * delays[m])))
        group = delays[m] - np.gradient(phase, 2 * np.pi * CHECK)
        axes[1, 2].semilogx(CHECK, 1e6 * group, color=color, label=name + " group delay")
        axes[1, 2].axhline(1e6 * delays[m], color=color, ls=":", lw=1)
    axes[0, 0].set(title="Characteristic admittance", ylabel="|Yc| (mS)")
    axes[0, 1].set(title=f"Propagation over {data.length / 1e3:g} km", ylabel="|H|", ylim=(0, 1.05))
    axes[0, 2].set(title=f"Held-out fitting error: Yc {report['Yc']['poles']} poles, H {report['H']['poles']} poles",
                   ylabel="Error (%)", yscale="log")
    axes[1, 0].set(title="Modal resistance (skin and ground return)", ylabel="ohm/km")
    axes[1, 1].set(title="Modal inductance", ylabel="mH/km")
    axes[1, 2].set(title="Dispersion: group delay and wavefront delay", ylabel="us")
    axes[1, 2].set_ylim(0.9e6 * delays.min(), 1e6 * max(delays.max() * 1.6, delays.min() * 1.1))
    for ax in axes.flat:
        ax.set_xlabel("Frequency (Hz)")
        ax.grid(color="0.9", lw=0.5)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle(f"{data.document['name']}: {data.length / 1e3:g} km, ideally transposed, "
                 f"{1 / data.document['earth']['conductivity']:g} ohm m soil")
    fig.savefig(destination, dpi=150)
    plt.close(fig)


def main(names):
    for path in sorted(LINES.glob("*.line.json")):
        name = path.name.removesuffix(".line.json")
        if names and name not in names:
            continue
        print(f"== {name}", flush=True)
        data = compute(path, FREQUENCY)
        out = LINES / name
        out.mkdir(exist_ok=True)
        columns = np.column_stack((data.frequency, data.z[:, 0].real, data.z[:, 0].imag, data.z[:, 1].real, data.z[:, 1].imag,
                                   data.y[:, 0].real, data.y[:, 0].imag, data.y[:, 1].real, data.y[:, 1].imag))
        np.savetxt(out / "constants.csv", columns, delimiter=",", fmt="%.16e", comments="",
                   header="f_Hz,z0_real_ohm_m,z0_imag_ohm_m,z1_real_ohm_m,z1_imag_ohm_m,y0_real_S_m,y0_imag_S_m,y1_real_S_m,y1_imag_S_m")
        sequence = data.sequence(60.0)
        fd_delays = data.wave_delays
        berg, surge, berg_delays = bergeron(data, sequence)
        sequence.update({
            "length_m": data.length,
            "voltage_kv": data.document["voltage_kv"],
            "short_circuit_mva": data.document["short_circuit_mva"],
            "surge_impedance_ohm": {"zero": float(surge[0]), "positive": float(surge[1])},
            "bergeron_delays_s": {"zero": float(berg_delays[0]), "positive": float(berg_delays[1])},
            "wavefront_delays_s": {"ground": float(fd_delays[0]), "aerial": float(fd_delays[1])},
            "characteristic_admittance_60hz_S": {
                "ground": [float(v) for v in (data.characteristic_admittance(2j * np.pi * np.array([60.0]))[0, 0].real,
                                              data.characteristic_admittance(2j * np.pi * np.array([60.0]))[0, 0].imag)],
                "aerial": [float(v) for v in (data.characteristic_admittance(2j * np.pi * np.array([60.0]))[0, 1].real,
                                              data.characteristic_admittance(2j * np.pi * np.array([60.0]))[0, 1].imag)]},
        })
        (out / "sequence.json").write_text(json.dumps(sequence, indent=2) + "\n")
        (out / "tower.json").write_text(json.dumps(tower(data), indent=2) + "\n")
        (out / "bergeron.fit.json").write_text(json.dumps(berg, indent=2) + "\n")
        yc_fit, yc_error, yc_order = fit_admittance(data)
        line, h_error, h_order, passive = fit_propagation(data, yc_fit)
        (out / "fd.fit.json").write_text(json.dumps(line, indent=2) + "\n")
        report = {
            "line": data.document["name"], "length_m": data.length,
            "fit_band_Hz": [float(FREQUENCY[0]), float(FREQUENCY[-1])], "fit_points": len(FREQUENCY),
            "held_out_points": len(CHECK),
            "wavefront_delays_s": fd_delays.tolist(),
            "Yc": {"poles": len(yc_fit["poles"]), "max_relative_error": yc_error, "target": YC_TARGET},
            "H": {"poles": len(line["H"]["output"]["poles"]), "max_scaled_error": h_error, "target": H_TARGET,
                  "magnitude_floor": 1e-3},
            "passivity_scan": passive,
            "sequence_60Hz": sequence,
        }
        (out / "fit-report.json").write_text(json.dumps(report, indent=2) + "\n")
        figure(data, line, report, out / "fits.png")
        print(f"  delays: ground {fd_delays[0] * 1e6:.2f} us, aerial {fd_delays[1] * 1e6:.2f} us; "
              f"Zc {surge[0]:.1f}/{surge[1]:.1f} ohm; Yc {yc_order} poles ({yc_error:.1e}), H {h_order} poles ({h_error:.1e})",
              flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
