#!/usr/bin/env python3
"""Compare accepted IDA adaptive step sizes from multiple trace CSV files."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from plot_common import (
    INK_MUTED,
    apply_style,
    case_name,
    configured_output,
    event_times,
    load_study,
    plt,
)


SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "studies", nargs="+", type=Path, help="solver JSON files used for the runs"
    )
    parser.add_argument("--output", type=Path, help="output PNG path")
    args = parser.parse_args()

    traces = []
    for path in args.studies:
        study, case, study_path = load_study(path)
        trace_path = configured_output(study, study_path, "step_trace_file")
        frame = pd.read_csv(trace_path, usecols=["segment", "t", "h"])
        frame = frame.apply(pd.to_numeric, errors="coerce").dropna()
        frame = frame[frame["h"] > 0.0]
        if frame.empty:
            raise ValueError(f"no positive IDA steps in {trace_path}")
        traces.append((case_name(case, study), study, case, frame))

    output = args.output or Path("ACTIVSg_adaptive_step.png")
    first_study = traces[0][1]
    events = event_times(first_study)
    nominal_frequency = float(traces[0][2].get("params", {}).get("freq_base", 60.0))

    apply_style()
    figure, axis = plt.subplots(figsize=(7.0, 3.6))
    for index, (label, _, _, frame) in enumerate(traces):
        first_segment = frame["segment"].min()
        for segment, rows in frame.groupby("segment", sort=True):
            axis.step(
                rows["t"],
                rows["h"],
                where="post",
                color=SERIES[index % len(SERIES)],
                linewidth=0.9,
                alpha=0.9,
                solid_joinstyle="round",
                label=label if segment == first_segment else None,
            )

    reference = 1.0 / (4.0 * nominal_frequency)
    axis.axhline(
        reference,
        color=INK_MUTED,
        linestyle=(0, (4, 3)),
        linewidth=0.7,
    )
    axis.annotate(
        r"$1/(4f_0)$",
        xy=(1.01, reference),
        xycoords=("axes fraction", "data"),
        va="center",
        color=INK_MUTED,
        annotation_clip=False,
    )
    fault_on = [time for time, kind in events if kind.lower() == "fault_on"]
    fault_off = [time for time, kind in events if kind.lower() == "fault_off"]
    if fault_on and fault_off:
        axis.axvspan(
            min(fault_on), max(fault_off), color=INK_MUTED, alpha=0.12, linewidth=0
        )
    else:
        for event_time, _ in events:
            axis.axvline(event_time, color=INK_MUTED, linestyle=":", linewidth=0.8)

    axis.set_yscale("log")
    axis.set_xlim(0.0, max(float(study["tmax"]) for _, study, _, _ in traces))
    axis.set_ylim(1.0e-3, 1.0e0)
    axis.set_xlabel(r"$t$ $-$ Time [s]")
    axis.set_ylabel(r"$h$ $-$ Time step [s]")
    axis.grid(True, which="major")
    axis.set_axisbelow(True)
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.14),
        ncol=min(len(traces), 4),
        frameon=False,
        handlelength=1.5,
        columnspacing=1.6,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output)
    plt.close(figure)
    for label, _, _, frame in traces:
        print(
            f"{label:<12} {len(frame):>6,} steps, "
            f"h={frame['h'].min():.6g} to {frame['h'].max():.6g} s"
        )
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
