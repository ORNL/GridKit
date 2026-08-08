#!/usr/bin/env python3
"""Plot every validated monitor in a PhasorDynamics study against its reference.

The study file is the single source of truth: each entry in its `monitors` list
that carries a `reference_file` becomes one figure of three panels -- GridKit,
PowerWorld, and their difference -- all sharing a y-scale so the difference is
read against the magnitude it is a difference of. Fault details are taken from
the study and case files so the annotation cannot drift out of step with what
was simulated.

Columns are matched by machine or bus rather than by position, falling back to
position with a warning if the keys do not line up. Rows are matched
positionally, the way compareCSV does it.

Usage:
    plot_validation.py <study.json> [--rundir DIR] [--outdir DIR]

Figures are written as <output-file-stem>.png, e.g. wecc.omega.png.
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#e3e2df"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")

# Axis label and description per monitored variable; anything unlisted falls
# back to the variable name in per-unit.
FIELDS = {
    "omega": ("Δω  (pu)", "generator speed deviation"),
    "speed": ("ω  (pu)", "generator speed"),
    "delta": ("δ  (rad)", "generator rotor angle"),
    "p": ("P  (pu)", "generator real power"),
    "q": ("Q  (pu)", "generator reactive power"),
    "Vm": ("|V|  (pu)", "bus voltage magnitude"),
    "Va": ("∠V  (rad)", "bus voltage angle"),
}


def sci(value, digits=3):
    """Format a value as mantissa x 10^exp for mathtext."""
    exp = int(np.floor(np.log10(abs(value)))) if value else 0
    return f"{value / 10.0 ** exp:.{digits}f} \\times 10^{{{exp}}}"


def column_key(column):
    """Reduce a GridKit column to the bare machine or bus a reference names.

    `Genrou_1032_C_genrou_omega` -> `1032 C`,  `Bus_1001_Vm` -> `1001`.
    """
    parts = column.split("_")
    if len(parts) < 3:
        return column.strip()
    device_class, body = parts[0], parts[1:-1]
    # Component ids repeat the class as a suffix; bus ids do not
    if body and body[-1].lower() == device_class.lower():
        body = body[:-1]
    return " ".join(body) if body else column.strip()


def load(path, key=column_key):
    with open(path) as f:
        labels = [key(c) for c in f.readline().strip().split(",")[1:]]
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return labels, data[:, 0], data[:, 1:]


def fault_description(study, case, study_dir):
    """Where the fault is, how long it lasts, and its impedance."""
    on = next(e for e in study["events"] if e["type"] == "fault_on")
    off = next(e for e in study["events"] if e["type"] == "fault_off")

    faults = [d for d in case["devices"] if d["class"] == "BusFault"]
    fault = faults[on["element_id"]]
    bus = next(b for b in case["buses"] if b["number"] == fault["ports"]["bus"])

    r, x = fault["params"]["R"], fault["params"]["X"]
    return on["time"], off["time"], [
        f"bus {bus['number']}  {bus['name'].strip()}  ({bus['params']['kv']:g} kV)",
        f"{on['time']:g}–{off['time']:g} s  ({(off['time'] - on['time']) * 1e3:g} ms)",
        f"Z = {r:g} + j{x:g} pu",
    ]


def metrics(err, ref, t):
    """eps_rmse = sqrt(dt / (N T)) ||X - Xhat||_2 ;  eps_inf = ||X - Xhat||_inf

    Each is also normalised by itself applied to the reference. The sqrt(dt/NT)
    cancels in the RMSE ratio, leaving a Frobenius ratio; the L-inf ratio is
    what DynamicSimulation reports as "Total max" under error_type "relative".
    dt comes from the median positive step, since event times appear twice.
    """
    steps = np.diff(t)
    dt = np.median(steps[steps > 0])
    duration = t[-1] - t[0]
    n_series = err.shape[1]

    return {
        "rmse": np.sqrt(dt / (n_series * duration)) * np.linalg.norm(err),
        "inf": np.abs(err).max(),
        "rel_rmse": np.linalg.norm(err) / np.linalg.norm(ref),
        "rel_inf": np.abs(err).max() / np.abs(ref).max(),
        "dt": dt,
        "duration": duration,
        "n_series": n_series,
    }


def plot_field(out_csv, ref_csv, png, title, axis_label, fault, quiet=False):
    gk_labels, gk_t, gk = load(out_csv)
    pw_labels, pw_t, pw = load(ref_csv)

    if gk_labels != pw_labels:
        if sorted(gk_labels) == sorted(pw_labels):
            gk = gk[:, [gk_labels.index(m) for m in pw_labels]]
        else:
            print(f"  note: column keys differ; matching {png.name} by position")

    n = min(len(gk_t), len(pw_t))
    if len(gk_t) != len(pw_t):
        print(f"  note: row counts differ ({len(gk_t)} vs {len(pw_t)}); using first {n}")
    err = gk[:n] - pw[:n]
    m = metrics(err, pw[:n], pw_t[:n])

    fault_on, fault_off, fault_lines = fault

    fig, axes = plt.subplots(3, 1, figsize=(11, 9.5), sharex=True, dpi=150)
    fig.patch.set_facecolor(SURFACE)

    panels = (
        (gk_t, gk, SERIES[0], "GridKit"),
        (pw_t, pw, SERIES[1], "PowerWorld"),
        (pw_t[:n], err, SERIES[2], "Difference"),
    )

    for ax, (t, y, color, panel_title) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        ax.axvspan(fault_on, fault_off, color=INK, alpha=0.05, lw=0)
        ax.plot(t, y, color=color, lw=0.6, alpha=0.35, solid_joinstyle="round")
        ax.set_title(panel_title, loc="left", color=INK, fontsize=12,
                     fontweight="semibold", pad=8)
        ax.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=INK_MUTED, labelsize=9)
        ax.set_ylabel(axis_label, color=INK_MUTED, fontsize=9)

        # The power-of-ten rides on the axis as an offset, not in the label
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
        ax.yaxis.get_offset_text().set(color=INK_MUTED, fontsize=9)

    # The signal panels take the data's own range; the difference is meant to
    # read as flat against it, so it keeps the identical span. Where the signal
    # already straddles zero all three share limits outright. Where it does not
    # (voltage, real power) the difference panel is recentred on zero, which
    # holds units-per-pixel equal without dragging an unused decade of empty
    # axis into the signal panels.
    lo = min(gk.min(), pw.min())
    hi = max(gk.max(), pw.max())
    pad = 0.08 * (hi - lo)
    lo, hi = lo - pad, hi + pad

    axes[0].set_ylim(lo, hi)
    axes[1].set_ylim(lo, hi)
    if lo < 0.0 < hi:
        axes[2].set_ylim(lo, hi)
    else:
        axes[2].set_ylim(-0.5 * (hi - lo), 0.5 * (hi - lo))

    axes[2].set_xlim(gk_t[0], gk_t[-1])
    axes[2].set_xlabel("time (s)", color=INK_MUTED, fontsize=9)

    X_LABEL, X_RMSE, X_INF = 0.578, 0.800, 0.985
    Y_HEAD, Y_ABS, Y_REL = 0.945, 0.845, 0.755

    axes[2].add_patch(Rectangle(
        (0.558, 0.700), 0.440, 0.290, transform=axes[2].transAxes,
        facecolor=SURFACE, edgecolor=GRID, lw=0.8, zorder=3,
    ))

    cells = (
        (Y_HEAD, "", "$\\epsilon_{\\mathrm{RMSE}}$", "$\\epsilon_{\\infty}$", INK_MUTED),
        (Y_ABS, "Absolute", f"${sci(m['rmse'])}$ pu", f"${sci(m['inf'])}$ pu", INK),
        (Y_REL, "Relative", f"{m['rel_rmse'] * 100:.3f} %", f"{m['rel_inf'] * 100:.3f} %", INK),
    )
    for y, row_label, rmse_cell, inf_cell, color in cells:
        axes[2].text(X_LABEL, y, row_label, transform=axes[2].transAxes, ha="left",
                     va="center", color=INK_MUTED, fontsize=10, zorder=4)
        for x, cell in ((X_RMSE, rmse_cell), (X_INF, inf_cell)):
            axes[2].text(x, y, cell, transform=axes[2].transAxes, ha="right",
                         va="center", color=color, fontsize=10, zorder=4)

    # Selective direct label: where the worst disagreement happens. The zero
    # line sits anywhere from mid-panel to the floor depending on the field, so
    # lead away from whichever edge the point is near.
    r, c = np.unravel_index(np.abs(err).argmax(), err.shape)
    y0, y1 = axes[2].get_ylim()
    x0, x1 = axes[2].get_xlim()
    below = (err[r, c] - y0) / (y1 - y0) < 0.5
    left = (pw_t[r] - x0) / (x1 - x0) > 0.6
    axes[2].annotate(
        f"max at {pw_labels[c]},  t = {pw_t[r]:.3f} s",
        xy=(pw_t[r], err[r, c]),
        xytext=(-18 if left else 18, 46 if below else -46),
        textcoords="offset points",
        color=INK, fontsize=9,
        ha="right" if left else "left",
        va="bottom" if below else "top",
        bbox=dict(facecolor=SURFACE, edgecolor="none", alpha=0.85, pad=2),
        arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=0.8),
    )

    fig.text(0.055, 0.978, title, color=INK, fontsize=15,
             fontweight="semibold", va="top")
    fig.text(0.055, 0.947, f"{m['n_series']} series · {n} samples",
             color=INK_MUTED, fontsize=9.5, va="top")

    fig.text(0.985, 0.978, "Fault", color=INK, fontsize=10,
             fontweight="semibold", ha="right", va="top")
    for i, line in enumerate(fault_lines):
        fig.text(0.985, 0.951 - i * 0.0225, line, color=INK_MUTED, fontsize=9.5,
                 ha="right", va="top")

    fig.tight_layout(rect=(0, 0, 1, 0.915))
    fig.savefig(png, facecolor=SURFACE)
    plt.close(fig)

    if not quiet:
        print(f"wrote {png}")
        print(f"  eps_rmse : {m['rmse']:.6e} pu   ({m['rel_rmse'] * 100:.3f} %)")
        print(f"  eps_inf  : {m['inf']:.6e} pu   ({m['rel_inf'] * 100:.3f} %)"
              f"  at {pw_labels[c]}, t = {pw_t[r]:.3f} s")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("study", type=Path, help="study (solver) JSON file")
    parser.add_argument("--rundir", type=Path,
                        help="directory holding the monitor output files "
                             "(default: the study file's directory)")
    parser.add_argument("--outdir", type=Path,
                        help="directory to write figures to "
                             "(default: the study file's directory)")
    args = parser.parse_args()

    study_dir = args.study.resolve().parent
    rundir = args.rundir or study_dir
    outdir = args.outdir or study_dir
    outdir.mkdir(parents=True, exist_ok=True)

    study = json.loads(args.study.read_text())
    case = json.loads((study_dir / study["system_model_file"]).read_text())
    fault = fault_description(study, case, study_dir)
    case_name = case.get("header", {}).get("case_name", args.study.stem)

    validated = [m for m in study.get("monitors", []) if m.get("reference_file")]
    if not validated:
        sys.exit(f"{args.study}: no monitors carry a reference_file")

    for mon in validated:
        out_csv = rundir / mon["file_name"]
        ref_csv = study_dir / mon["reference_file"]
        if not out_csv.exists():
            print(f"skipping {mon['file_name']}: not found in {rundir}")
            continue

        # `Genrou_*_omega` -> `omega`;  the variable names the field
        include = mon.get("include", "")
        patterns = [include] if isinstance(include, str) else include
        var = patterns[0].split("_")[-1] if patterns else Path(out_csv).stem
        axis_label, description = FIELDS.get(var, (f"{var}  (pu)", var))

        plot_field(out_csv, ref_csv, outdir / (out_csv.stem + ".png"),
                   f"{case_name} · {description}", axis_label, fault)


if __name__ == "__main__":
    main()
