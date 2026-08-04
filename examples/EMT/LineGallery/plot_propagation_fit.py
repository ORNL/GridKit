#!/usr/bin/env python3
"""Compare the true propagation function against its fitted approximation.

Assembles the model-free reference

    P_ref(j omega) = conj(Ti) diag(H) Tv^T

from the UniversalLineModel sweep monitors and composes the fitted
approximation from propagation.model.json according to its domain:

    modal: P_fit = Gout_fit diag(exp(-j omega tau_m)) Gin_fit
    phase: P_fit = H_fit exp(-j omega tau)

The composition is gauge invariant, so the comparison shows real fitting
error only. Two figures land in output/<line>/: element magnitudes and
element phases, each with the actual response on the left and the
approximation on the right.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import OUTPUT_DIR, ULM_APP, gallery_lines, run_ulm
from plot_factors import read_factors
from sweep_vs_fit import eval_model


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--line", default="345kv-horizontal")
    parser.add_argument("--all", action="store_true",
                        help="run every line with a gallery solver spec")
    parser.add_argument("--h-domain", default="phase",
                        choices=("modal", "phase"))
    parser.add_argument("--ulm-app", type=Path, default=ULM_APP)
    return parser.parse_args()


def compose_fit(propagation: dict, omega: np.ndarray, k: int) -> np.ndarray:
    """The fitted propagation matrix per sample, by artifact domain."""
    domain = propagation["domain"]
    if domain == "modal":
        tau = np.array(propagation["delays"]["tau"])
        gin = eval_model(propagation["input"], omega).reshape(len(omega), k, k)
        gout = eval_model(propagation["output"], omega).reshape(len(omega), k, k)
        out = np.zeros((len(omega), k, k), complex)
        for m in range(len(omega)):
            out[m] = gout[m] @ np.diag(np.exp(-1j * omega[m] * tau)) @ gin[m]
        return out
    if domain == "phase":
        tau = propagation["delay"]["tau"]
        h = eval_model(propagation["H"], omega).reshape(len(omega), k, k)
        return h * np.exp(-1j * omega * tau)[:, None, None]
    raise ValueError(f"unknown propagation domain: {domain}")


def side_by_side(name, omega, reference, fitted, transform, ylabel, title,
                 path, log_y, suptitle="propagation function P = Gout diag(H) Gin"):
    k = reference.shape[1]
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharex=True, sharey=True)
    colors = plt.cm.viridis(np.linspace(0.0, 0.85, k * k))

    for i in range(k):
        for j in range(k):
            color = colors[i * k + j]
            label = f"({i},{j})" if k <= 3 else None
            axes[0].plot(omega, transform(reference[:, i, j]),
                         color=color, lw=1.1, label=label)
            axes[1].plot(omega, transform(fitted[:, i, j]),
                         color=color, lw=1.1)

    for ax, panel in ((axes[0], "actual"), (axes[1], "approximation")):
        ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")
        ax.set_title(f"{title} — {panel}", fontsize=11)
        ax.set_xlabel("omega [rad/s]")
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel(ylabel)
    if k <= 3:
        axes[0].legend(fontsize=8, ncol=3)

    error = np.sqrt(np.mean(np.abs(fitted - reference) ** 2)
                    / max(np.mean(np.abs(reference) ** 2), 1e-30))
    axes[1].annotate(f"composed rel rms {error:.2e}",
                     xy=(0.03, 0.03), xycoords="axes fraction", fontsize=9)

    fig.suptitle(f"{name}: {suptitle}", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"wrote {path}")


def render(name: str) -> float:
    """Figures from the artifacts already in output/<line>/.

    Returns the composed relative rms error of the approximation.
    """
    line_dir = OUTPUT_DIR / name
    omega, k, tv, ti, h, _ = read_factors(line_dir / "response.csv")
    propagation = json.loads((line_dir / "propagation.model.json").read_text())

    reference = np.zeros((len(omega), k, k), complex)
    for m in range(len(omega)):
        reference[m] = np.conj(ti[m]) @ np.diag(h[m]) @ tv[m].T

    fitted = compose_fit(propagation, omega, k)

    domain = propagation["domain"]
    side_by_side(name, omega, reference, fitted,
                 np.abs, "|P entry|", f"element magnitudes ({domain})",
                 line_dir / "propagation_mag.png", log_y=True)
    side_by_side(name, omega, reference, fitted,
                 lambda entry: np.degrees(np.unwrap(np.angle(entry))),
                 "unwrapped phase [deg]", f"element phases ({domain})",
                 line_dir / "propagation_phase.png", log_y=False)

    if domain == "phase":
        # The fitter's own view: the time-shifted target
        # Hmps = P exp(+j omega tau) against the rational approximation,
        # before the delay is reapplied.
        tau = propagation["delay"]["tau"]
        target = reference * np.exp(1j * omega * tau)[:, None, None]
        approx = eval_model(propagation["H"], omega).reshape(len(omega), k, k)
        suptitle = "time-shifted target Hmps = P exp(+s tau) vs rational fit"
        side_by_side(name, omega, target, approx,
                     np.abs, "|Hmps entry|", "element magnitudes",
                     line_dir / "hmps_mag.png", log_y=True,
                     suptitle=suptitle)
        side_by_side(name, omega, target, approx,
                     lambda entry: np.degrees(np.unwrap(np.angle(entry))),
                     "unwrapped phase [deg]", "element phases",
                     line_dir / "hmps_phase.png", log_y=False,
                     suptitle=suptitle)

    return float(np.sqrt(np.mean(np.abs(fitted - reference) ** 2)
                         / max(np.mean(np.abs(reference) ** 2), 1e-30)))


def run(name: str, h_domain: str = "phase",
        ulm_app: Path = ULM_APP) -> float:
    print(f"=== {name} ===")
    run_ulm(name, h_domain=h_domain, ulm_app=ulm_app)
    return render(name)


def main() -> None:
    args = parse_args()
    for name in gallery_lines() if args.all else [args.line]:
        run(name, args.h_domain, args.ulm_app)


if __name__ == "__main__":
    main()
