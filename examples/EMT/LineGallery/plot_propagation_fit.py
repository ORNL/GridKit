#!/usr/bin/env python3
"""Compare the true propagation function against its fitted approximation.

Assembles the model-free reference

    H_ref(j omega) = conj(Ti) diag(h) Tv^T

from the UniversalLineModel sweep monitors, where h holds the modal
propagation, and composes the fitted approximation from
propagation.model.json as

    H_fit = Hmin_fit exp(-j omega tau)

with the one shared delay the fit removed. Four figures land in
output/<line>/: the propagation element magnitudes and phases, actual
beside approximation, and the same pair for the unwound matrix Hmin the
fitter actually sees.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import OUTPUT_DIR, ULM_APP, gallery_lines, read_modal_response, run_ulm
from sweep_vs_fit import eval_model


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--line", default="345kv-horizontal")
    parser.add_argument("--all", action="store_true",
                        help="run every line with a gallery solver spec")
    parser.add_argument("--ulm-app", type=Path, default=ULM_APP)
    return parser.parse_args()


def side_by_side(name, omega, reference, fitted, transform, ylabel, title,
                 path, log_y, suptitle):
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
    axes[1].annotate(f"rel rms {error:.2e}",
                     xy=(0.03, 0.03), xycoords="axes fraction", fontsize=9)

    fig.suptitle(f"{name}: {suptitle}", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"wrote {path}")


def render(name: str) -> float:
    """Figures from the artifacts already in output/<line>/.

    Returns the relative rms error of the fitted propagation function.
    """
    line_dir = OUTPUT_DIR / name
    omega, k, tv, ti, h, _ = read_modal_response(line_dir / "response.csv")
    propagation = json.loads((line_dir / "propagation.model.json").read_text())
    tau = propagation["delay"]["tau"]

    reference = np.zeros((len(omega), k, k), complex)
    for m in range(len(omega)):
        reference[m] = np.conj(ti[m]) @ np.diag(h[m]) @ tv[m].T

    hmin = eval_model(propagation["Hmin"], omega).reshape(len(omega), k, k)
    fitted = hmin * np.exp(-1j * omega * tau)[:, None, None]

    suptitle = "propagation function H = conj(Ti) diag(h) Tv^T"
    side_by_side(name, omega, reference, fitted,
                 np.abs, "|H entry|", "element magnitudes",
                 line_dir / "propagation_mag.png", log_y=True,
                 suptitle=suptitle)
    side_by_side(name, omega, reference, fitted,
                 lambda entry: np.degrees(np.unwrap(np.angle(entry))),
                 "unwrapped phase [deg]", "element phases",
                 line_dir / "propagation_phase.png", log_y=False,
                 suptitle=suptitle)

    # The fitter's own view: the unwound target against its rational
    # approximation, before the delay is reapplied.
    target = reference * np.exp(1j * omega * tau)[:, None, None]
    suptitle = "unwound target Hmin = H exp(+s tau) vs rational fit"
    side_by_side(name, omega, target, hmin,
                 np.abs, "|Hmin entry|", "element magnitudes",
                 line_dir / "hmin_mag.png", log_y=True, suptitle=suptitle)
    side_by_side(name, omega, target, hmin,
                 lambda entry: np.degrees(np.unwrap(np.angle(entry))),
                 "unwrapped phase [deg]", "element phases",
                 line_dir / "hmin_phase.png", log_y=False, suptitle=suptitle)

    return float(np.sqrt(np.mean(np.abs(fitted - reference) ** 2)
                         / max(np.mean(np.abs(reference) ** 2), 1e-30)))


def run(name: str, ulm_app: Path = ULM_APP) -> float:
    print(f"=== {name} ===")
    run_ulm(name, ulm_app=ulm_app)
    return render(name)


def main() -> None:
    args = parse_args()
    for name in gallery_lines() if args.all else [args.line]:
        run(name, args.ulm_app)


if __name__ == "__main__":
    main()
