"""Compare PSCAD's line constants with GridWorkbench's for every frequency-dependent line.

PSCAD's Line Constants Program writes the phase-domain series impedance and
shunt admittance it computed from the tower (``<line>_zm/zp.out`` and
``<line>_ym/yp.out``: magnitudes and arguments per element, SI per metre) and
the modal characteristic admittance it fitted (``<line>_ycmp.out``). This
script reduces both tools' balanced matrices to the ground and aerial modes
and reports the relative differences of modal R, L and C over frequency, so
the part of the frequency-dependent mismatch that comes from the two
line-constants programs is known separately from the time-domain solvers.

Outputs ``lines/<line>/constants-comparison.png`` and ``.json``.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "lib"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from scenarios import CODES, project_name  # noqa: E402

LINES = HERE / "lines"
PSCAD = HERE / "results" / "pscad"
COLORS = ("#0072B2", "#D55E00")
REPORT_HZ = (60.0, 1e3, 1e4, 1e5)


def read_matrix(path):
    data = np.loadtxt(path, skiprows=1)
    return data[:, 1], data[:, 2:11].reshape(-1, 3, 3)


def pscad_modal(folder, name):
    f, zmag = read_matrix(folder / f"{name}_zm.out")
    _, zarg = read_matrix(folder / f"{name}_zp.out")
    _, ymag = read_matrix(folder / f"{name}_ym.out")
    _, yarg = read_matrix(folder / f"{name}_yp.out")
    z = zmag * np.exp(1j * np.radians(zarg))
    y = ymag * np.exp(1j * np.radians(yarg))
    # Balanced matrices: ground mode = self + 2 mutual, aerial = self - mutual.
    self_z, mutual_z = z[:, 0, 0], (z[:, 0, 1] + z[:, 0, 2] + z[:, 1, 2]) / 3
    self_y, mutual_y = y[:, 0, 0], (y[:, 0, 1] + y[:, 0, 2] + y[:, 1, 2]) / 3
    return f, np.column_stack((self_z + 2 * mutual_z, self_z - mutual_z)), np.column_stack((self_y + 2 * mutual_y, self_y - mutual_y))


def gridworkbench_modal(line, f):
    table = np.genfromtxt(LINES / line / "constants.csv", delimiter=",", names=True)
    ac = table["f_Hz"] > 0
    source = table["f_Hz"][ac]
    def interp(name):
        return np.interp(np.log10(f), np.log10(source), table[name][ac])
    z = np.column_stack((interp("z0_real_ohm_m") + 1j * interp("z0_imag_ohm_m"), interp("z1_real_ohm_m") + 1j * interp("z1_imag_ohm_m")))
    y = np.column_stack((interp("y0_real_S_m") + 1j * interp("y0_imag_S_m"), interp("y1_real_S_m") + 1j * interp("y1_imag_S_m")))
    return z, y


def main(argv):
    for line, code in CODES.items():
        if argv and line not in argv:
            continue
        folder = PSCAD / project_name(line, "fd", "fault")
        name = f"{code}_fd"
        if not (folder / f"{name}_zm.out").exists():
            print("missing PSCAD constants for", line)
            continue
        f, z_ps, y_ps = pscad_modal(folder, name)
        z_gw, y_gw = gridworkbench_modal(line, f)
        omega = 2 * np.pi * f
        r_ps, l_ps, c_ps = z_ps.real, z_ps.imag / omega[:, None], y_ps.imag / omega[:, None]
        r_gw, l_gw, c_gw = z_gw.real, z_gw.imag / omega[:, None], y_gw.imag / omega[:, None]
        report = {"line": line, "frequency_band_Hz": [float(f[0]), float(f[-1])], "modes": ["ground", "aerial"], "at_frequency": {}}
        for hz in REPORT_HZ:
            k = int(np.argmin(np.abs(f - hz)))
            report["at_frequency"][f"{hz:g}"] = {
                "R_relative_difference": ((r_ps[k] - r_gw[k]) / r_gw[k]).tolist(),
                "L_relative_difference": ((l_ps[k] - l_gw[k]) / l_gw[k]).tolist(),
                "C_relative_difference": ((c_ps[k] - c_gw[k]) / c_gw[k]).tolist(),
                "Z_magnitude_relative_difference": (np.abs(z_ps[k] - z_gw[k]) / np.abs(z_gw[k])).tolist(),
                "R_pscad_ohm_per_km": (1e3 * r_ps[k]).tolist(), "R_gridworkbench_ohm_per_km": (1e3 * r_gw[k]).tolist(),
                "L_pscad_mH_per_km": (1e6 * l_ps[k]).tolist(), "L_gridworkbench_mH_per_km": (1e6 * l_gw[k]).tolist(),
                "C_pscad_nF_per_km": (1e12 * c_ps[k]).tolist(), "C_gridworkbench_nF_per_km": (1e12 * c_gw[k]).tolist(),
            }
        report["max_Z_magnitude_relative_difference_over_band"] = np.max(np.abs(z_ps - z_gw) / np.abs(z_gw), axis=0).tolist()
        report["max_C_relative_difference_over_band"] = np.max(np.abs(c_ps - c_gw) / c_gw, axis=0).tolist()
        (LINES / line / "constants-comparison.json").write_text(json.dumps(report, indent=2) + "\n")
        fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), layout="constrained")
        for m, color, label in zip((0, 1), COLORS, ("Ground mode", "Aerial mode")):
            axes[0, 0].loglog(f, 1e3 * r_gw[:, m], color=color, label=label + " GridWorkbench")
            axes[0, 0].loglog(f[::25], 1e3 * r_ps[::25, m], ".", color=color, ms=4, label=label + " PSCAD")
            axes[0, 1].semilogx(f, 1e6 * l_gw[:, m], color=color, label=label + " GridWorkbench")
            axes[0, 1].semilogx(f[::25], 1e6 * l_ps[::25, m], ".", color=color, ms=4, label=label + " PSCAD")
            axes[0, 2].semilogx(f, 1e12 * c_gw[:, m], color=color, label=label + " GridWorkbench")
            axes[0, 2].semilogx(f[::25], 1e12 * c_ps[::25, m], ".", color=color, ms=4, label=label + " PSCAD")
            axes[1, 0].semilogx(f, 100 * (r_ps[:, m] - r_gw[:, m]) / r_gw[:, m], color=color, label=label)
            axes[1, 1].semilogx(f, 100 * (l_ps[:, m] - l_gw[:, m]) / l_gw[:, m], color=color, label=label)
            axes[1, 2].semilogx(f, 100 * np.abs(z_ps[:, m] - z_gw[:, m]) / np.abs(z_gw[:, m]), color=color, label=label + " |Z|")
            axes[1, 2].semilogx(f, 100 * np.abs(c_ps[:, m] - c_gw[:, m]) / c_gw[:, m], "--", color=color, label=label + " C")
        axes[0, 0].set(title="Modal resistance", ylabel="ohm/km")
        axes[0, 1].set(title="Modal inductance", ylabel="mH/km")
        axes[0, 2].set(title="Modal capacitance", ylabel="nF/km")
        axes[1, 0].set(title="PSCAD - GridWorkbench, resistance", ylabel="%")
        axes[1, 1].set(title="PSCAD - GridWorkbench, inductance", ylabel="%")
        axes[1, 2].set(title="PSCAD - GridWorkbench, |Z| and C", ylabel="%", yscale="log")
        for ax in axes.flat:
            ax.set_xlabel("Frequency (Hz)")
            ax.grid(color="0.9", lw=0.5)
            ax.legend(frameon=False, fontsize=8)
        fig.suptitle(f"{line}: line constants from the same tower geometry, PSCAD LCP against GridWorkbench emt.parameters")
        fig.savefig(LINES / line / "constants-comparison.png", dpi=140)
        plt.close(fig)
        at60 = report["at_frequency"]["60"]
        print(f"{line:18s} 60 Hz: dR {np.array(at60['R_relative_difference']) * 100} %  dL {np.array(at60['L_relative_difference']) * 100} %  "
              f"dC {np.array(at60['C_relative_difference']) * 100} %   max |dZ|/|Z| over band {np.array(report['max_Z_magnitude_relative_difference_over_band']) * 100} %")


if __name__ == "__main__":
    main(sys.argv[1:])
