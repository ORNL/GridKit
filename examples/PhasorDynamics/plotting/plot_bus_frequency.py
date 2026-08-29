#!/usr/bin/env python3
"""Plot bus frequency from a DynamicSimulation monitor CSV."""

from __future__ import annotations

import argparse
from pathlib import Path

from plot_common import (
    BLUE_RAMP,
    GRID,
    INK_MUTED,
    apply_style,
    case_name,
    configured_output,
    default_plot_path,
    event_times,
    load_study,
    matplotlib,
    plt,
)

import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection


def load_angles(csv_path: Path) -> tuple[np.ndarray, np.ndarray]:
    raw_columns = list(pd.read_csv(csv_path, nrows=0).columns)
    time_column = next(
        (column for column in raw_columns if column.strip() in {"t", "time"}),
        None,
    )
    angle_columns = [column for column in raw_columns if column.strip().endswith("_Va")]
    if time_column is None:
        raise ValueError(f"missing t or time column in {csv_path}")
    if not angle_columns:
        raise ValueError(f"no bus voltage-angle columns ending in _Va in {csv_path}")

    selected = [time_column, *angle_columns]
    frame = pd.read_csv(csv_path, usecols=selected, dtype=str)[selected]
    time_values = pd.to_numeric(frame[time_column], errors="coerce")
    valid = time_values.notna()
    time = time_values[valid].to_numpy(dtype=float)
    angle = frame.loc[valid, angle_columns].to_numpy(dtype=float)

    order = np.argsort(time, kind="stable")
    time = time[order]
    angle = angle[order]
    keep = np.concatenate((time[1:] != time[:-1], np.array([True])))
    return time[keep], angle[keep]


def physical_samples(
    time: np.ndarray, events: list[tuple[float, str]]
) -> np.ndarray:
    keep = np.ones_like(time, dtype=bool)
    if len(time) < 2:
        return keep
    reach = 1.5 * float(np.median(np.diff(time)))
    for event_time, _ in events:
        keep &= np.abs(time - event_time) > reach
    return keep


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path, help="solver JSON used for the run")
    parser.add_argument("--input", type=Path, help="override monitor CSV path")
    parser.add_argument("--output", type=Path, help="output PNG path")
    parser.add_argument("--ylim", type=float, nargs=2, metavar=("LOW", "HIGH"))
    args = parser.parse_args()

    study, case, study_path = load_study(args.study)
    csv_path = args.input or configured_output(study, study_path, "output_file")
    output = args.output or default_plot_path(study_path, study, "bus_frequency")
    nominal_frequency = float(case.get("params", {}).get("freq_base", 60.0))
    events = event_times(study)

    time, angle = load_angles(csv_path)
    frequency = nominal_frequency + np.gradient(
        np.unwrap(angle, axis=0), time, axis=0
    ) / (2.0 * np.pi)
    physical = physical_samples(time, events)

    peak = np.abs(frequency[physical] - nominal_frequency).max(axis=0)
    peak_min = float(peak.min())
    peak_max = float(peak.max())
    if peak_max <= peak_min:
        peak_max = peak_min + np.finfo(float).eps
    color_map = matplotlib.colors.LinearSegmentedColormap.from_list(
        "gridkit_blue", BLUE_RAMP
    )
    color_norm = matplotlib.colors.Normalize(vmin=peak_min, vmax=peak_max)
    draw_order = np.argsort(peak)
    segments = [
        np.column_stack([time, frequency[:, index]]) for index in draw_order
    ]

    apply_style()
    figure, axis = plt.subplots(figsize=(7.0, 3.8))
    axis.add_collection(
        LineCollection(
            segments,
            linewidths=0.35,
            colors=color_map(color_norm(peak[draw_order])),
            alpha=0.55,
            capstyle="round",
            zorder=2,
        )
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

    axis.axhline(
        nominal_frequency,
        color=INK_MUTED,
        linestyle=(0, (4, 3)),
        linewidth=0.7,
        zorder=3,
    )
    axis.set_xlim(float(time[0]), float(time[-1]))
    if args.ylim:
        axis.set_ylim(*args.ylim)
    else:
        band = float(
            np.percentile(np.abs(frequency[physical] - nominal_frequency), 99.9)
        )
        band = max(band, 1.0e-6)
        axis.set_ylim(nominal_frequency - 1.12 * band, nominal_frequency + 1.12 * band)
    axis.ticklabel_format(axis="y", style="plain", useOffset=False)
    axis.set_xlabel(r"$t$ $-$ Time [s]")
    axis.set_ylabel(r"$f$ $-$ Bus frequency [Hz]")
    axis.set_title(f"{case_name(case, study)}   ·   {frequency.shape[1]:,} buses", pad=8)
    axis.grid(True)
    axis.set_axisbelow(True)

    colorbar = figure.colorbar(
        plt.cm.ScalarMappable(norm=color_norm, cmap=color_map),
        ax=axis,
        fraction=0.032,
        pad=0.015,
    )
    colorbar.set_label(r"Peak $|\Delta f|$ [Hz]", labelpad=6, fontsize=9)
    colorbar.ax.tick_params(labelsize=8, length=0, colors=INK_MUTED)
    colorbar.outline.set_visible(False)

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output)
    plt.close(figure)
    print(f"buses           {frequency.shape[1]:,}")
    print(f"samples         {len(time):,} ({time[0]:g} to {time[-1]:g} s)")
    print(
        f"physical range  {frequency[physical].min():.4f} to "
        f"{frequency[physical].max():.4f} Hz"
    )
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
