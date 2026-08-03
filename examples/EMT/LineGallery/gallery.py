#!/usr/bin/env python3
"""Sweep and fit every catalog line, then plot responses and fit quality.

For each line in ../Lines the gallery runs the FrequencyResponse
application on the committed solver specification, renders a pre-fitting
response overview (Yc, Zc, gamma as alpha and beta, tau, H), runs the
UniversalLineModel application to fit the line, renders the Yc fit
accuracy figure, and collects pole counts and error statistics into
output/stats.json and output/stats.md.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output"
BUILD = HERE.parents[2] / "build"
DEFAULT_SWEEP_APP = BUILD / "application" / "EMT" / "FrequencyResponse" / "FrequencyResponse"
DEFAULT_ULM_APP = BUILD / "application" / "EMT" / "UniversalLineModel" / "UniversalLineModel"

LINES = [
    "69kv-wood-pole",
    "138kv-delta",
    "345kv-horizontal",
    "500kv-double-circuit",
    "765kv-horizontal",
]

# The propagation factor target reflects what factor-based fitting
# supports for untransposed lines with nearly degenerate aerial modes.
ULM_ARGS = ["--h-target", "0.05"]
ULM_LINE_ARGS = {"500kv-double-circuit": ["--yc-max-poles", "16"]}

FIT_LINE = re.compile(
    r"^(?P<label>\w+): VectorFitting: .*rel rms (?P<rms>[0-9.e+-]+), "
    r"worst channel (?P<worst>[0-9.e+-]+), order (?P<order>\d+)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep-app", type=Path, default=DEFAULT_SWEEP_APP)
    parser.add_argument("--ulm-app", type=Path, default=DEFAULT_ULM_APP)
    parser.add_argument("--lines", nargs="*", default=LINES)
    return parser.parse_args()


def read_response(csv_path: Path) -> tuple[np.ndarray, int, dict[str, np.ndarray]]:
    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"{csv_path} must contain a header")
        rows = list(reader)

    k = 0
    while f"Overhead_Yc_real_{k}_{k}" in reader.fieldnames:
        k += 1
    if k == 0:
        raise ValueError(f"{csv_path} contains no Overhead_Yc columns")

    omega = np.array([float(r["omega"]) for r in rows])
    m = len(rows)

    def matrix(prefix: str) -> np.ndarray:
        out = np.zeros((m, k * k), dtype=np.complex128)
        for sample, r in enumerate(rows):
            for i in range(k):
                for j in range(k):
                    out[sample, i * k + j] = complex(
                        float(r[f"{prefix}_real_{i}_{j}"]),
                        float(r[f"{prefix}_imag_{i}_{j}"]),
                    )
        return out

    def modal(prefix: str, complex_: bool) -> np.ndarray:
        if complex_:
            return np.array(
                [
                    [
                        complex(
                            float(r[f"{prefix}_real_{i}"]),
                            float(r[f"{prefix}_imag_{i}"]),
                        )
                        for i in range(k)
                    ]
                    for r in rows
                ]
            )
        return np.array(
            [[float(r[f"{prefix}_{i}"]) for i in range(k)] for r in rows]
        )

    data = {
        "Yc": matrix("Overhead_Yc"),
        "Zc": matrix("Overhead_Zc"),
        "Alpha": modal("Overhead_Alpha", False),
        "Beta": modal("Overhead_Beta", False),
        "Tau": modal("Overhead_Tau", False),
        "H": modal("Overhead_H", True),
    }
    return omega, k, data


def matrix_panel(ax, frequency, values, k, title, ylabel) -> None:
    for channel in range(values.shape[1]):
        label = f"{channel // k}{channel % k}" if k <= 3 else None
        ax.plot(frequency, np.abs(values[:, channel]), linewidth=1.0, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(title if k <= 3 else f"{title} ({values.shape[1]} entries)")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.3)
    if k <= 3:
        ax.legend(fontsize=6, ncol=3)


def modal_panel(ax, frequency, values, title, ylabel, log_y=True) -> None:
    for mode in range(values.shape[1]):
        ax.plot(frequency, values[:, mode], linewidth=1.2, label=f"{mode}")
    ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(title="mode", fontsize=6, ncol=2)


def write_responses_figure(name: str, omega, k, data, plot_file: Path) -> None:
    frequency = omega / (2.0 * math.pi)

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.0), sharex=True)
    matrix_panel(axes[0, 0], frequency, data["Yc"], k, "Characteristic admittance", "|Yc| [S]")
    matrix_panel(axes[0, 1], frequency, data["Zc"], k, "Characteristic impedance", "|Zc| [ohm]")
    modal_panel(axes[0, 2], frequency, data["Alpha"], "Gamma: attenuation", "alpha [Np/m]")
    modal_panel(axes[1, 0], frequency, data["Beta"], "Gamma: phase", "beta [rad/m]")
    modal_panel(
        axes[1, 1], frequency, 1.0e6 * data["Tau"], "Modal delay", "tau [us]", log_y=False
    )
    modal_panel(axes[1, 2], frequency, np.abs(data["H"]), "Propagation magnitude", "|H| [-]")
    for ax in axes[1]:
        ax.set_xlabel("frequency [Hz]")

    fig.suptitle(f"{name}: frequency response ({k} conductors after reduction)")
    fig.tight_layout()
    fig.savefig(plot_file, dpi=180)
    plt.close(fig)
    print(f"wrote {plot_file}")


def load_factor(model: dict) -> dict:
    rows, cols = model["rows"], model["cols"]
    poles = np.array([complex(re, im) for re, im in model["poles"]])
    residues = np.array(
        [
            [complex(re, im) for row in matrix for re, im in row]
            for matrix in model["residues"]
        ],
        dtype=np.complex128,
    ).reshape(len(poles), rows * cols)
    d = np.zeros(rows * cols)
    if "D" in model:
        d = np.array(model["D"], dtype=np.float64).reshape(rows * cols)
    return {"poles": poles, "residues": residues, "d": d}


def evaluate_factor(factor: dict, omega: np.ndarray) -> np.ndarray:
    s = 1j * omega[:, None]
    return (1.0 / (s - factor["poles"][None, :])) @ factor["residues"] + factor["d"][None, :]


def write_ycfit_figure(name: str, omega, k, yc_samples, model_file: Path, plot_file: Path) -> float:
    model = json.loads(model_file.read_text())
    factor = load_factor(model)
    fitted = evaluate_factor(factor, omega)
    source = yc_samples
    frequency = omega / (2.0 * math.pi)
    rel_rms = float(
        np.sqrt(np.mean(np.abs(fitted - source) ** 2) / np.mean(np.abs(source) ** 2))
    )

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8))
    matrix_panel(axes[0], frequency, source, k, "Yc samples", "|Yc| [S]")
    matrix_panel(axes[1], frequency, fitted, k, "Rational approximation", "|Yc| [S]")
    axes[1].sharey(axes[0])

    scale = np.sqrt(np.mean(np.abs(source) ** 2, axis=1))
    error = np.abs(fitted - source).max(axis=1) / scale
    axes[2].plot(frequency, error, linewidth=1.2, color="tab:red")
    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].set_title("Pointwise error")
    axes[2].set_ylabel("max |error| / rms sample")
    axes[2].grid(True, which="both", alpha=0.3)
    for ax in axes:
        ax.set_xlabel("frequency [Hz]")

    fig.suptitle(
        f"{name}: Yc fit, {len(factor['poles'])} poles, relative RMS {rel_rms:.2e}"
    )
    fig.tight_layout()
    fig.savefig(plot_file, dpi=180)
    plt.close(fig)
    print(f"wrote {plot_file}")
    return rel_rms


def run_line(name: str, args: argparse.Namespace) -> dict:
    print(f"=== {name} ===")
    subprocess.run(
        [str(args.sweep_app), f"{name}.solver.json"],
        check=True,
        cwd=HERE,
        stdout=subprocess.DEVNULL,
    )
    omega, k, data = read_response(OUTPUT_DIR / f"{name}.response.csv")
    write_responses_figure(name, omega, k, data, OUTPUT_DIR / f"{name}.responses.png")

    model_dir = OUTPUT_DIR / name
    command = [
        str(args.ulm_app),
        str(HERE.parent / "Lines" / f"{name}.line.json"),
        "-o",
        str(model_dir),
        *ULM_ARGS,
        *ULM_LINE_ARGS.get(name, []),
    ]
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    if result.returncode < 0 or "failed" in result.stdout + result.stderr:
        raise RuntimeError(f"UniversalLineModel failed for {name}:\n{result.stdout}")

    fits = {}
    for line in result.stdout.splitlines():
        match = FIT_LINE.match(line.strip())
        if match:
            fits[match.group("label")] = {
                "poles": int(match.group("order")),
                "rel_rms": float(match.group("rms")),
                "worst_channel_rel_rms": float(match.group("worst")),
            }
    passive = "Yc fit is passive" in result.stdout

    yc_rel_rms = write_ycfit_figure(
        name, omega, k, data["Yc"], model_dir / "yc.model.json",
        OUTPUT_DIR / f"{name}.ycfit.png",
    )

    propagation = json.loads((model_dir / "propagation.model.json").read_text())
    return {
        "line": name,
        "conductors": k,
        "targets_met": result.returncode == 0,
        "yc_passive": passive,
        "yc_rel_rms_check": yc_rel_rms,
        "delays_us": [1.0e6 * tau for tau in propagation["delays"]],
        "fits": fits,
    }


def write_stats(records: list[dict]) -> None:
    (OUTPUT_DIR / "stats.json").write_text(json.dumps(records, indent=2) + "\n")

    lines = [
        "# Line gallery fit statistics",
        "",
        "| Line | K | Yc poles | Yc rel RMS | Gin poles | Gin rel RMS "
        "| Gout poles | Gout rel RMS | Targets met | Yc passive |",
        "| ---- | - | -------- | ---------- | --------- | ----------- "
        "| ---------- | ------------ | ----------- | ---------- |",
    ]
    for r in records:
        fits = r["fits"]

        def cell(label: str, key: str) -> str:
            entry = fits.get(label)
            if entry is None:
                return "-"
            if key == "poles":
                return str(entry["poles"])
            return f"{entry['rel_rms']:.2e}"

        lines.append(
            f"| {r['line']} | {r['conductors']} "
            f"| {cell('Yc', 'poles')} | {cell('Yc', 'rms')} "
            f"| {cell('Gin', 'poles')} | {cell('Gin', 'rms')} "
            f"| {cell('Gout', 'poles')} | {cell('Gout', 'rms')} "
            f"| {'yes' if r['targets_met'] else 'no'} "
            f"| {'yes' if r['yc_passive'] else 'no'} |"
        )
    lines += [
        "",
        "Modal delays [us]: "
        + "; ".join(
            f"{r['line']}: " + ", ".join(f"{tau:.2f}" for tau in r["delays_us"])
            for r in records
        ),
        "",
    ]
    (OUTPUT_DIR / "stats.md").write_text("\n".join(lines))
    print(f"wrote {OUTPUT_DIR / 'stats.json'}")
    print(f"wrote {OUTPUT_DIR / 'stats.md'}")


def main() -> None:
    args = parse_args()
    for app in (args.sweep_app, args.ulm_app):
        if not app.exists():
            raise RuntimeError(f"application binary not found: {app}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    records = [run_line(name, args) for name in args.lines]
    write_stats(records)


if __name__ == "__main__":
    main()
