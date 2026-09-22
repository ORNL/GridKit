"""Shared publication style and layouts for the HICSS-60 figures."""

from pathlib import Path
from shutil import copy2

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import (FixedFormatter, FixedLocator, LogFormatterMathtext,
                               LogLocator, MultipleLocator, NullFormatter, NullLocator)


CASE_STYLES = {
    "Hawaii": ("Hawaii", "#3b6ea5"),
    "IEEE39": ("NE", "#2f9184"),
    "ACTIVSg200": ("Illinois", "#8663a6"),
    "ACTIVSg2000": ("Texas", "#cc7330"),
    "WECC240": ("WECC", "#b84a47"),
}
RC_PARAMS = {
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "Times New Roman", "Nimbus Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 15.0,
    "axes.labelsize": 15.0,
    "axes.edgecolor": "#222222",
    "axes.linewidth": 1.0,
    "xtick.labelsize": 14.0,
    "ytick.labelsize": 14.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 5.0,
    "ytick.major.size": 5.0,
    "xtick.minor.size": 2.8,
    "ytick.minor.size": 2.8,
    "lines.linewidth": 1.6,
    "lines.solid_capstyle": "round",
    "legend.fontsize": "medium",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}


def subplots(nrows=1, ncols=1, *, sharex=False, sharey=False, legend=False,
             figsize=None):
    """Use compact gaps, with optional dimensions for a fixed paper footprint."""
    width = figsize[0] if figsize else 7.0
    scale = width / 7.0
    plt.rcParams.update({key: value * scale if isinstance(value, float) else value
                         for key, value in RC_PARAMS.items()})
    left, right, bottom = (value * scale for value in (0.85, 0.50, 0.70))
    if figsize:
        right = 0.14 * scale
    top = (0.58 if legend else 0.18) * scale
    gap_x = 0.22 * scale
    gap_y = (0.18 if sharex else 0.42) * scale
    panel_width = (width - left - right - (ncols - 1) * gap_x) / ncols
    panel_height = panel_width * (1.0 if ncols > 1 else 0.55)
    height = bottom + top + nrows * panel_height + (nrows - 1) * gap_y
    if figsize:
        height = figsize[1]
        panel_height = (height - bottom - top - (nrows - 1) * gap_y) / nrows
    figure, axes = plt.subplots(nrows, ncols, figsize=(width, height),
                                sharex=sharex, sharey=sharey)
    figure.subplots_adjust(left=left / width, right=1 - right / width,
                           bottom=bottom / height, top=1 - top / height,
                           wspace=gap_x / panel_width, hspace=gap_y / panel_height)
    for axis in figure.axes:
        axis.set_box_aspect(panel_height / panel_width)
        axis.grid(True, which="major", color="#dde2e7", linewidth=0.6 * scale)
        axis.set_axisbelow(True)
        axis.tick_params(which="both", top=True, right=True)
    return figure, axes


def labels(figure, *, x=None, y=None):
    bounds = figure.subplotpars
    if x:
        figure.supxlabel(x, x=(bounds.left + bounds.right) / 2,
                          y=0.04 / figure.get_figheight())
    if y:
        figure.supylabel(y, x=0.02 / figure.get_figwidth(),
                          y=(bounds.bottom + bounds.top) / 2)


def legend(figure, handles=None):
    if handles is None:
        handles = figure.axes[0].get_legend_handles_labels()[0]
    center = (figure.subplotpars.left + figure.subplotpars.right) / 2
    top = figure.subplotpars.top + 0.14 * figure.get_figwidth() / 7.0 / figure.get_figheight()
    figure.legend(handles=handles, loc="lower center", ncol=len(handles),
                  frameon=False, bbox_to_anchor=(center, top), borderaxespad=0,
                  handlelength=1.2, handletextpad=0.5, columnspacing=1.2)


def panel_label(axis, index, text, *, bottom=False):
    axis.text(0.965, 0.04 if bottom else 0.96, rf"$\mathbf{{({chr(97 + index)})}}$ {text}",
              transform=axis.transAxes, ha="right", va="bottom" if bottom else "top",
              fontsize=plt.rcParams["xtick.labelsize"])


def step_axis(figure, axis, duration, events, frequency):
    event_markers(axis, events)
    reference = 1.0 / (4.0 * frequency)
    axis.axhline(reference, color="#777777", linestyle="--", linewidth=1.0)
    axis.annotate(r"$\frac{1}{4f}$", (1.01, reference),
                  xycoords=("axes fraction", "data"), va="center", annotation_clip=False)
    axis.set_yscale("log")
    axis.set_xlim(0, duration)
    axis.set_ylim(1e-3, 1e0)
    axis.xaxis.set_major_locator(MultipleLocator(5))
    axis.xaxis.set_minor_locator(MultipleLocator(2.5))
    labels(figure, x=r"$t$ – Time [sec]", y=r"$h$ – Time step [sec]")


def event_markers(axis, times):
    for time in times:
        axis.axvline(time, color="#777777", linestyle=":",
                     linewidth=0.7 * plt.rcParams["axes.linewidth"], zorder=0)


def set_mu_axis(axis, mus):
    """Use sampled values for previews and decades for dense sweeps."""
    axis.set_xscale("log")
    lower, upper = min(mus), max(mus)
    if lower == upper:
        lower, upper = lower / 1.2, upper * 1.2
    axis.set_xlim(lower, upper)
    if len(mus) > 8:
        axis.xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
        axis.xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
        axis.xaxis.set_minor_locator(LogLocator(base=10.0, subs=tuple(range(2, 10))))
        axis.xaxis.set_minor_formatter(NullFormatter())
    else:
        axis.xaxis.set_major_locator(FixedLocator(mus))
        axis.xaxis.set_major_formatter(FixedFormatter([f"{mu:g}" for mu in mus]))
        axis.xaxis.set_minor_locator(NullLocator())
    # Keep endpoint labels inside the panel when columns sit close together.
    for value, label in zip(axis.get_xticks(), axis.get_xticklabels()):
        if value == min(mus):
            label.set_horizontalalignment("left")
        elif value == max(mus):
            label.set_horizontalalignment("right")


def save_figure(figure, path, *, tight=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    bounds = "tight" if tight else None
    figure.savefig(path, dpi=600, bbox_inches=bounds, pad_inches=0.03)
    if path.suffix.lower() == ".png":
        figure.savefig(path.with_suffix(".pdf"), bbox_inches=bounds, pad_inches=0.03)
    plt.close(figure)
    output = Path(__file__).resolve().parents[1] / "output"
    output.mkdir(exist_ok=True)
    for source in ([path, path.with_suffix(".pdf")] if path.suffix.lower() == ".png" else [path]):
        target = output / source.name
        if source.resolve() != target.resolve():
            copy2(source, target)
    print(f"Wrote {path}")
