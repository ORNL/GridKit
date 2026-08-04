#!/usr/bin/env python3
"""Compare the true propagation function against its fitted approximation.

Assembles the model-free reference

    H_ref(j omega) = conj(Ti) diag(h) Tv^T

from the UniversalLineModel sweep monitors, where h holds the modal
propagation, and composes the fitted approximation from
propagation.model.json as the modal sum

    H_fit = sum_m Hmin_m exp(-j omega tau_m)

with the per-mode delays the fit removed. Figures land in
output/<line>/: the propagation element magnitudes and phases, actual
beside approximation, and the same pair per mode for the unwound
rank-one dyad each Hmin_m fit actually sees.
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
    modes = propagation["modes"]

    reference = np.zeros((len(omega), k, k), complex)
    for m in range(len(omega)):
        reference[m] = np.conj(ti[m]) @ np.diag(h[m]) @ tv[m].T

    fitted_modes = [
        eval_model(mode["Hmin"], omega).reshape(len(omega), k, k)
        for mode in modes
    ]
    delays = [mode["delay"]["tau"] for mode in modes]
    fitted = sum(hmin * np.exp(-1j * omega * tau)[:, None, None]
                 for hmin, tau in zip(fitted_modes, delays))

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

    # The fitter's own view, mode by mode: each unwound rank-one dyad
    # against its rational approximation, before the delay and the
    # other modes are reapplied.
    for stale in line_dir.glob("hmin_*.png"):
        stale.unlink()
    for q, (hmin, tau) in enumerate(zip(fitted_modes, delays)):
        dyad = np.einsum("mi,m,mj->mij", np.conj(ti[:, :, q]), h[:, q],
                         tv[:, :, q])
        target = dyad * np.exp(1j * omega * tau)[:, None, None]
        suptitle = (f"unwound mode {q} dyad Hmin{q} = "
                    "conj(ti_m) h_m exp(+s tau_m) tv_m^T vs rational fit")
        side_by_side(name, omega, target, hmin,
                     np.abs, f"|Hmin{q} entry|", "element magnitudes",
                     line_dir / f"hmin{q}_mag.png", log_y=True,
                     suptitle=suptitle)
        side_by_side(name, omega, target, hmin,
                     lambda entry: np.degrees(np.unwrap(np.angle(entry))),
                     "unwrapped phase [deg]", "element phases",
                     line_dir / f"hmin{q}_phase.png", log_y=False,
                     suptitle=suptitle)

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
