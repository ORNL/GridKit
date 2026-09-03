#!/usr/bin/env python3
"""Plot IDA solver statistics as a function of simulation time.

Reads the ida_steps.json accepted-step history emitted by DynamicSimulation
and renders, per validation case, the cumulative and per-step counter
trajectories along with the step size. Cases with only aggregate statistics
(ida_stats.json) appear in the summary table but have no time-series figure.
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STATS_DIR = Path(sys.argv[1])
OUT_DIR = Path(sys.argv[2])
OUT_DIR.mkdir(parents=True, exist_ok=True)

CASES = [
    "IEEE39",
    "ACTIVSg200",
    "ACTIVSg500",
    "WECC240",
    "Hawaii",
    "ACTIVSg2000",
    "ACTIVSg10k",
]

# counter_delta group -> key, and the label used in plots
COUNTERS = [
    ("nonlinear_solver", "iterations", "Nonlinear iterations"),
    ("linear_solver", "jacobian_evals", "Jacobian evaluations"),
    ("integrator", "residual_evals", "Residual evaluations"),
    ("integrator", "linear_solver_setups", "Linear solver setups"),
    ("integrator", "error_test_failures", "Error test failures"),
    ("nonlinear_solver", "convergence_failures", "Nonlinear conv. failures"),
]

COLORS = plt.get_cmap("tab10")


def load_steps(case):
    path = STATS_DIR / f"{case}.ida_steps.json"
    if not path.exists():
        return None
    with open(path) as f:
        doc = json.load(f)

    t, dt, order, seg_id = [], [], [], []
    deltas = {label: [] for _, _, label in COUNTERS}
    seg_bounds = []

    # The last accepted step of a segment normally lands past the segment end:
    # IDA_ONE_STEP takes a full step and the remaining output times are then
    # produced by interpolation. Its work belongs to the segment, but plotting
    # it at its raw end time would run the time axis backwards when the next
    # segment resumes at the boundary. Clamp those samples to the segment end.
    for seg in doc["segments"]:
        seg_start, seg_end = seg["start_time"], seg["end_time"]
        seg_bounds.append((seg_start, seg_end))
        for s in seg["steps"]:
            t.append(min(s["step_end_time"], seg_end))
            dt.append(s["last_step"])
            order.append(s["current_order"])
            seg_id.append(seg["segment_index"])
            cd = s["counter_delta"]
            for grp, key, label in COUNTERS:
                deltas[label].append(cd.get(grp, {}).get(key, 0))

    return {
        "t": np.asarray(t),
        "dt": np.asarray(dt),
        "order": np.asarray(order),
        "seg": np.asarray(seg_id),
        "deltas": {k: np.asarray(v) for k, v in deltas.items()},
        "seg_bounds": seg_bounds,
        "doc": doc,
    }


def load_summary(case):
    path = STATS_DIR / f"{case}.ida_stats.json"
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def event_times(d):
    """Interior segment boundaries are the discrete events (fault on/off)."""
    return [b[0] for b in d["seg_bounds"][1:]]


def aggregate_row(case, s):
    """Row for a case that has aggregate stats but no per-step history."""
    i, n, l = s["integrator"], s["nonlinear_solver"], s["linear_solver"]
    return {
        "case": case,
        "accepted_steps": i["steps"],
        "segments": s["segment_count"],
        "min_dt": float("nan"),
        "max_dt": float("nan"),
        "Nonlinear iterations": n["iterations"],
        "Jacobian evaluations": l["jacobian_evals"],
        "Residual evaluations": i["residual_evals"],
        "Linear solver setups": i["linear_solver_setups"],
        "Error test failures": i["error_test_failures"],
        "Nonlinear conv. failures": n["convergence_failures"],
        "summary": s,
    }


def plot_case(case, d):
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    evs = event_times(d)

    # --- cumulative counters -------------------------------------------
    ax = axes[0]
    for i, (_, _, label) in enumerate(COUNTERS):
        cum = np.cumsum(d["deltas"][label])
        if cum[-1] == 0:
            continue
        style = "--" if label == "Residual evaluations" else "-"
        # a log axis cannot render the leading zeros, so mask them out
        y = np.where(cum > 0, cum.astype(float), np.nan)
        ax.plot(d["t"], y, label=f"{label} ({cum[-1]:,})",
                color=COLORS(i), lw=1.4, ls=style)
    ax.set_yscale("log")
    ax.set_ylabel("Cumulative count")
    ax.set_title(f"{case} — IDA solver statistics vs. simulation time")
    ax.legend(fontsize=8, loc="lower right", ncol=2)
    ax.grid(alpha=0.3, which="both")

    # --- per-step counters ---------------------------------------------
    ax = axes[1]
    for i, (_, _, label) in enumerate(COUNTERS[:3]):
        style = "--" if label == "Residual evaluations" else "-"
        ax.plot(d["t"], d["deltas"][label], label=label,
                color=COLORS(i), lw=0.8, alpha=0.85, ls=style)
    ax.set_ylabel("Per accepted step")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    # --- step size ------------------------------------------------------
    ax = axes[2]
    ax.semilogy(d["t"], d["dt"], color="k", lw=1.0)
    ax.set_ylabel("Step size (s)")
    ax.set_xlabel("Simulation time (s)")
    ax.grid(alpha=0.3)

    for a in axes:
        for ev in evs:
            a.axvline(ev, color="tab:gray", ls="--", lw=0.9, alpha=0.8)

    fig.tight_layout()
    out = OUT_DIR / case / "figures" / f"{case}.ida_stats.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_overview(data):
    """One figure comparing cumulative work across all cases."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    panels = [
        ("Nonlinear iterations", axes[0][0]),
        ("Jacobian evaluations", axes[0][1]),
        ("Residual evaluations", axes[1][0]),
        ("Error test failures", axes[1][1]),
    ]
    for label, ax in panels:
        for i, (case, d) in enumerate(data.items()):
            cum = np.cumsum(d["deltas"][label])
            y = np.where(cum > 0, cum.astype(float), np.nan)
            ax.plot(d["t"], y, label=case, color=COLORS(i), lw=1.3)
        ax.set_title(label)
        ax.set_xlabel("Simulation time (s)")
        ax.set_ylabel("Cumulative count")
        ax.set_yscale("log")
        ax.grid(alpha=0.3, which="both")
    axes[0][0].legend(fontsize=8, ncol=2)
    fig.suptitle("IDA cumulative solver work vs. simulation time — all validated cases")
    fig.tight_layout()
    out = OUT_DIR / "IdaStatistics" / "all_cases.ida_stats.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main():
    data = {}
    rows = []
    for case in CASES:
        d = load_steps(case)
        if d is None:
            s = load_summary(case)
            if s is None:
                print(f"skip {case}: no statistics")
                continue
            print(f"{case}: aggregate only (no step history)")
            rows.append(aggregate_row(case, s))
            continue
        data[case] = d
        out = plot_case(case, d)
        print(f"wrote {out}")

        s = load_summary(case)
        tot = {label: int(np.sum(d["deltas"][label])) for _, _, label in COUNTERS}
        rows.append(
            {
                "case": case,
                "accepted_steps": int(len(d["t"])),
                "segments": len(d["seg_bounds"]),
                "min_dt": float(np.min(d["dt"])),
                "max_dt": float(np.max(d["dt"])),
                **tot,
                "summary": s,
            }
        )

    if data:
        print(f"wrote {plot_overview(data)}")

    with open(OUT_DIR / "IdaStatistics" / "ida_stats_table.json", "w") as f:
        json.dump(rows, f, indent=2)

    # console table
    hdr = ["case", "steps", "nonlin it", "jac", "resid", "lsetup", "errfail", "min dt", "max dt"]
    print()
    print(" | ".join(f"{h:>12}" for h in hdr))
    for r in rows:
        print(" | ".join(f"{v:>12}" for v in [
            r["case"], r["accepted_steps"],
            r["Nonlinear iterations"], r["Jacobian evaluations"],
            r["Residual evaluations"], r["Linear solver setups"],
            r["Error test failures"],
            f"{r['min_dt']:.2e}", f"{r['max_dt']:.2e}"]))


if __name__ == "__main__":
    main()
