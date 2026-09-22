#!/usr/bin/env python3
"""Generate Figure 6: accepted IDA step sizes for Texas at five tolerances."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import plot as paper
from common import trace as ida
from common.simulation import ROOT, read_json
from matplotlib.colors import ListedColormap
import numpy as np


HERE = Path(__file__).resolve().parent


def read_studies():
    studies = [(path, read_json(path)) for path in (HERE / "solvers").glob("*.solver.json")]
    studies.sort(key=lambda item: item[1]["rel_tol"], reverse=True)
    if len(studies) != 5:
        raise ValueError("Expected five tolerance-pair solver files")
    common = None
    for path, study in studies:
        if not np.isclose(study["rel_tol"], 100 * study["abs_tol"], rtol=1e-12, atol=0):
            raise ValueError(f"{path}: rel_tol must be 100 times abs_tol")
        settings = {key: value for key, value in study.items()
                    if key not in {"rel_tol", "abs_tol", "solver_trace_file"}}
        if common is not None and settings != common:
            raise ValueError("Tolerance studies must otherwise share the same settings")
        common = settings
    return studies


def plot(studies):
    source, study = studies[0]
    case = read_json(source.parent / study["system_model_file"])
    ends = ida.boundaries(study)
    rtols = [study["rel_tol"] for _, study in studies]
    colors = [color for _, color in paper.CASE_STYLES.values()]
    figure, axis = paper.subplots(legend=True)
    for (source, study), color in zip(studies, colors):
        trace = ida.load(source.parent / study["solver_trace_file"], ends)
        for times, sizes in ida.step_segments(trace, ends):
            axis.step(times, sizes, where="pre", color=color)
    paper.step_axis(figure, axis, study["tmax"], ends[1:-1],
                    float(case.get("params", {}).get("freq_base", 60.0)))
    bounds = axis.get_position()
    strip = figure.add_axes([bounds.x0, bounds.y1 + 0.04 / figure.get_figheight(),
                             bounds.width, 0.10 / figure.get_figheight()])
    strip.imshow([np.arange(len(rtols))], cmap=ListedColormap(colors),
                 aspect="auto", interpolation="nearest")
    strip.set_yticks([])
    strip.set_xticks(range(len(rtols)), [rf"$10^{{{int(np.log10(rtol))}}}$" for rtol in rtols])
    strip.tick_params(axis="x", top=True, bottom=False, labeltop=True,
                      labelbottom=False, length=0, pad=3)
    strip.set_xlabel(r"$r_{\mathrm{tol}}$", labelpad=3)
    strip.xaxis.set_label_position("top")
    for spine in strip.spines.values():
        spine.set_linewidth(0.6)
    paper.save_figure(figure, HERE / "figures/ida_tol.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path,
                        default=ROOT / "build/paper/application/PhasorDynamics/DynamicSimulation")
    parser.add_argument("--workdir", type=Path, default=ROOT / "build/HICSS-60/ida-tol")
    parser.add_argument("--run", action="store_true", help="Collect new data before plotting")
    args = parser.parse_args()
    studies = read_studies()
    if args.run:
        binary = args.binary.resolve(strict=True)
        for source, study in studies:
            elapsed = ida.collect(source, binary,
                                  args.workdir / source.name.removesuffix(".solver.json"))
            print(f"Texas rel_tol={study['rel_tol']:g}, abs_tol={study['abs_tol']:g}: "
                  f"{elapsed:.6f} s CPU", flush=True)
    plot(studies)


if __name__ == "__main__":
    main()
