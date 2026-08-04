#!/usr/bin/env python3
"""Time-domain form of every fitted function in this directory.

Each *.model.json holds one or more pole-residue models

    F(s) = E s + D + sum_q R_q / (s - p_q)

whose time-domain form is

    f(t) = E delta'(t) + D delta(t) + sum_q R_q exp(p_q t),   t >= 0.

The two singular terms live only at t = 0 and cannot be drawn, so they
are reported in the corner of each figure and what is plotted is the
regular part, entry by entry: a linear time axis on top, and log time
against log magnitude below, where each pole's decay rate reads off as
a straight segment.

propagation.model.json carries one fitted matrix and one transport
delay per mode, so it yields a figure per unwound modal kernel
Hmin_m(t) and one for the propagation function itself,

    H(t) = sum_m Hmin_m(t - tau_m),

the modal kernels standing off the origin by their own delays.

Recursive convolution is what an EMT step does with these kernels, so
the same recursion

    z_q[n] = alpha_q z_q[n-1] + beta_q x[n] + gamma_q x[n-1]
    alpha_q = exp(p_q dt),  I0 = (alpha_q - 1) / p_q
    I1      = (dt alpha_q - I0) / p_q
    beta_q  = I0 - I1 / dt,  gamma_q = I1 / dt

is driven here with a numerical impulse and compared against the
analytic kernel, so the plotted functions and the recursion that will
consume them are checked against each other.

Coefficients came from the 345 kV horizontal line of the line gallery.
Yc and the propagation function are 3 by 3 after reduction; Z, Y, and
gamma are fit on the raw 8 conductor description.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import lfilter

HERE = Path(__file__).resolve().parent

# Resolve the fastest pole, run until the kernel has fallen to a
# thousandth of its peak, and draw that window with this many points.
SAMPLES_PER_PERIOD = 10.0
ENVELOPE_FLOOR = 1.0e-3
BULK_FLOOR = 1.0e-2
PLOT_SAMPLES = 20_000
LOG_SAMPLES = 800
CHECK_SAMPLES = 4_000

UNITS = {
    "z": "ohm/m",
    "y": "S/m",
    "gamma": "1/m",
    "yc": "S",
    "hmin": "-",
    "h": "-",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", nargs="*", type=Path,
                        default=sorted(HERE.glob("*.model.json")))
    return parser.parse_args()


def load_model(spec: dict) -> dict:
    """Pole-residue model as arrays: poles (Q), residues (Q,K,C), D, E."""
    rows, cols = spec["rows"], spec["cols"]
    poles = np.array([complex(re, im) for re, im in spec["poles"]])
    residues = np.array(
        [[complex(re, im) for row in matrix for re, im in row]
         for matrix in spec["residues"]],
        dtype=complex,
    ).reshape(len(poles), rows, cols)

    def constant(key: str) -> np.ndarray:
        return (np.array(spec[key], dtype=float).reshape(rows, cols)
                if key in spec else np.zeros((rows, cols)))

    return {"poles": poles, "residues": residues, "D": constant("D"),
            "E": constant("E"), "rows": rows, "cols": cols}


def functions(path: Path) -> list[tuple[str, list[tuple[dict, float]]]]:
    """The named functions a model file defines, each a list of delayed
    pole-residue components (model, tau)."""
    spec = json.loads(path.read_text())
    if "modes" in spec:
        modes = [(load_model(m["Hmin"]), m["delay"]["tau"])
                 for m in spec["modes"]]
        named = [(f"Hmin{q}", [(model, 0.0)])
                 for q, (model, _) in enumerate(modes)]
        return named + [("H", modes)]
    return [(path.name.replace(".model.json", ""), [(load_model(spec), 0.0)])]


def kernel(model: dict, t: np.ndarray) -> np.ndarray:
    """Regular part sum_q R_q exp(p_q t), shape (N, rows, cols)."""
    values = np.einsum("qt,qij->tij", np.exp(np.outer(model["poles"], t)),
                       model["residues"])
    residual = np.abs(values.imag).max() / max(np.abs(values.real).max(), 1e-300)
    if residual > 1.0e-10:
        raise ValueError(f"kernel is not real, imaginary part {residual:.2e}")
    return values.real


def evaluate(components: list[tuple[dict, float]], t: np.ndarray) -> np.ndarray:
    """Sum of the delayed components on an absolute time grid."""
    model = components[0][0]
    total = np.zeros((len(t), model["rows"], model["cols"]))
    for model, tau in components:
        alive = t >= tau
        total[alive] += kernel(model, t[alive] - tau)
    return total


def memory(components: list[tuple[dict, float]],
           floor: float = ENVELOPE_FLOOR) -> float:
    """Time past the first arrival after which every component has
    fallen below that fraction of its own peak."""
    arrival = min(tau for _, tau in components)
    extent = 0.0
    for model, tau in components:
        poles = model["poles"]
        probe = np.geomspace(1.0e-3 / np.abs(poles).max(),
                             1.0e3 / np.abs(poles.real).min(), 2000)
        envelope = np.abs(kernel(model, probe)).max(axis=(1, 2))
        last = probe[np.flatnonzero(envelope > floor * envelope.max())[-1]]
        extent = max(extent, tau - arrival + float(last))
    return extent


def fastest_dt(components: list[tuple[dict, float]]) -> float:
    return 1.0 / (SAMPLES_PER_PERIOD
                  * max(np.abs(m["poles"]).max() for m, _ in components))


def recursion_check(model: dict, dt: float) -> float:
    """Recursive convolution of a numerical impulse against the kernel.

    The impulse is the unit-area sample x[0] = 1/dt, so from the second
    sample on the recursion has to track sum_q R_q exp(p_q t) to the
    order of the step size. This is the independent path: the figures
    come from the closed form, this number says the recursion an EMT
    step would run agrees with it. A transport delay shifts the kernel
    without touching the recursion, so each component checks alone.
    """
    poles = model["poles"]
    alpha = np.exp(poles * dt)
    i0 = (alpha - 1.0) / poles
    i1 = (dt * alpha - i0) / poles
    beta = i0 - i1 / dt
    gamma = i1 / dt

    t = dt * np.arange(CHECK_SAMPLES)
    x = np.zeros((CHECK_SAMPLES, model["cols"]), dtype=complex)
    x[0, :] = 1.0 / dt

    response = np.zeros((CHECK_SAMPLES, model["rows"]), dtype=complex)
    for q in range(len(poles)):
        state = lfilter([beta[q], gamma[q]], [1.0, -alpha[q]], x, axis=0)
        response += state @ model["residues"][q].T

    exact = kernel(model, t).sum(axis=2)
    return float(np.abs(response.real[1:] - exact[1:]).max()
                 / np.abs(exact[1:]).max())


def figure(name: str, components: list[tuple[dict, float]], unit: str,
           path: Path) -> None:
    """Linear time on top, log time against log magnitude below."""
    model = components[0][0]
    arrival = min(tau for _, tau in components)
    delays = sorted(tau for _, tau in components)
    unit = f"({unit})/s" if "/" in unit else f"{unit}/s"

    # The linear axis shows the bulk of the function, the log axis its
    # whole memory; delayed kernels keep their dead time on the linear
    # axis.
    linear = np.linspace(0.0, memory(components, BULK_FLOOR), PLOT_SAMPLES)
    logarithmic = np.geomspace(fastest_dt(components), memory(components),
                               LOG_SAMPLES)
    entries = model["rows"] * model["cols"]
    colors = plt.cm.viridis(np.linspace(0.0, 0.85, entries))

    lead = np.linspace(0.0, arrival, 200, endpoint=False)
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
    for grid, ax, transform in ((linear, axes[0], lambda v: v),
                                (logarithmic, axes[1], np.abs)):
        values = evaluate(components, arrival + grid).reshape(len(grid), -1)
        # The linear axis carries absolute time, the log axis time since
        # the first arrival, so delayed kernels stay readable on both.
        stamp = 1.0e6 * (grid + arrival if ax is axes[0] else grid)
        if ax is axes[0] and arrival > 0.0:
            stamp = np.concatenate([1.0e6 * lead, stamp])
            values = np.vstack([np.zeros((len(lead), entries)), values])
        for entry in range(entries):
            label = (f"({entry // model['cols']},{entry % model['cols']})"
                     if entries <= 9 else None)
            ax.plot(stamp, transform(values[:, entry]),
                    color=colors[entry], lw=1.0, label=label)
        ax.grid(True, which="both", alpha=0.25)
        ax.set_ylabel(f"f(t) [{unit}]")

    axes[0].set_xlabel("time [us]")
    axes[0].set_title("regular part on a linear time axis", fontsize=11)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("time [us]" if arrival == 0.0
                       else "time since first arrival [us]")
    axes[1].set_ylabel(f"|f(t)| [{unit}]")
    axes[1].set_title("magnitude against log time, one segment per pole",
                      fontsize=11)
    if entries <= 9:
        axes[0].legend(fontsize=8, ncol=min(entries, 3), loc="lower right")

    peak_d = max(np.abs(m["D"]).max() for m, _ in components)
    peak_e = max(np.abs(m["E"]).max() for m, _ in components)
    singular = [rf"$\max|D|$ = {peak_d:.4g} on the $\delta(t)$ term"]
    if peak_e > 0.0:
        singular.append(rf"$\max|E|$ = {peak_e:.4g} on the $\delta'(t)$ term")
    if arrival > 0.0:
        listed = ", ".join(f"{1.0e6 * tau:.2f}" for tau in delays)
        singular.append(rf"delays $\tau_m$ = {listed} us")
    axes[0].annotate("\n".join(singular), xy=(0.985, 0.95),
                     xycoords="axes fraction", ha="right", va="top", fontsize=9,
                     bbox={"boxstyle": "round", "fc": "white", "ec": "0.7"})

    poles = sum(len(m["poles"]) for m, _ in components)
    parts = ("" if len(components) == 1
             else f"{len(components)} delayed modes, ")
    fig.suptitle(f"{name}(t): {model['rows']}x{model['cols']} kernel, "
                 f"{parts}{poles} poles", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}")


def main() -> None:
    args = parse_args()
    if not args.models:
        raise RuntimeError(f"no *.model.json found in {HERE}")

    for path in args.models:
        for name, components in functions(path):
            agreement = max(recursion_check(model, fastest_dt([(model, tau)]))
                            for model, tau in components)
            print(f"{name}: {components[0][0]['rows']}x"
                  f"{components[0][0]['cols']}, "
                  f"{sum(len(m['poles']) for m, _ in components)} poles, "
                  f"memory {1.0e6 * memory(components):.2f} us, recursive "
                  f"convolution agrees to {agreement:.2e}")
            figure(name, components,
                   UNITS.get(name.lower().rstrip("0123456789"), "-"),
                   HERE / f"{name.lower()}_time.png")


if __name__ == "__main__":
    main()
