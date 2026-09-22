#!/usr/bin/env python3
"""Render the paper table from saved benchmark and validation measurements."""

import argparse
import sys
from math import isfinite
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.simulation import TRIALS, read_json

HERE = Path(__file__).resolve().parent
# PowerWorld times supplied with the paper; not remeasured here.
CASES = (
    ("IEEE39", "IEEE N.E.", 6.969),
    ("ACTIVSg200", "Synth. Illinois", 4.032),
    ("Hawaii", "Synth. Hawaii", 2.172),
    ("WECC240", "WECC", 4.375),
    ("ACTIVSg2000", "Synth. Texas", 14.703),
)


def scientific(value):
    mantissa, exponent = f"{value:.2e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def write_tables(results):
    markdown = [
        r"| | Variables | | Runtime (s) | | | $\epsilon_{\mathrm{RMSE}}^{\text{abs}}$ | | $\epsilon_{\infty}^{\text{rel}}$ | |",
        "|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|",
        r"| Case | Diff. | Alg. | PowerWorld | GridKit | Speedup | $\omega$ | $\lvert V\rvert$ | $\omega$ | $\lvert V\rvert$ |",
    ]
    latex = [
        r"\begin{table*}[!t]",
        r"\centering",
        r"\caption{Runtimes, speedup, and validation error of GridKit against the PowerWorld reference solution.}",
        r"\label{tab:verification-errors}",
        r"\begin{tabular}{l rr rrr rr rr}",
        r"\toprule",
        r" & \multicolumn{2}{c}{Variables} & \multicolumn{2}{c}{Runtime (s)}",
        r" & & \multicolumn{2}{c}{$\epsilon_{\mathrm{RMSE}}^\text{abs}$}",
        r" & \multicolumn{2}{c}{$\epsilon_{\infty}^\text{rel}$} \\",
        r"\cmidrule(lr){2-3} \cmidrule(lr){4-5} \cmidrule(lr){7-8} \cmidrule(lr){9-10}",
        r"Case & Diff. & Alg. & PowerWorld & GridKit & Speedup & $\omega$ & $|V|$ & $\omega$ & $|V|$ \\",
        r"\midrule",
    ]
    for result in results:
        elapsed = median(result["gridkit_seconds"])
        runtime = f"{elapsed:.3f}"
        speedup = f"{result['powerworld_seconds'] / elapsed:.1f}" + r"$\times$"
        omega, voltage = result["errors"]["omega"], result["errors"]["vmag"]
        cells = [result["case"], str(result["differential"]), str(result["algebraic"]),
                 f"{result['powerworld_seconds']:.3f}", runtime, speedup,
                 *(scientific(value) for value in (omega[0], voltage[0], omega[1], voltage[1]))]
        markdown.append("| " + " | ".join(cells).replace(r"$\times$", "x") + " |")
        latex.append(" & ".join(cells) + r" \\")
    latex.extend([r"\bottomrule", r"\end{tabular}", r"\end{table*}"])
    readme = HERE / "README.md"
    introduction = readme.read_text().split("\n## Results\n", 1)[0].rstrip()
    readme.write_text(introduction + "\n\n## Results\n\n" + "\n".join(markdown) + "\n")
    output = HERE.parent / "output"
    output.mkdir(exist_ok=True)
    (output / "latex.txt").write_text("\n".join(latex) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, default=HERE / "data/measurements.json")
    parser.add_argument("--validation", type=Path, default=HERE.parent / "validation/data/measurements.json")
    args = parser.parse_args()
    benchmark, validation = read_json(args.benchmark), read_json(args.validation)
    results = []
    for name, label, seconds in CASES:
        timed, checked = benchmark[name], validation[name]
        if timed["counts"] != checked["counts"]:
            raise ValueError(f"{name}: benchmark/validation variable counts differ")
        if checked["status"] != "ok":
            raise ValueError(f"{name}: validation did not pass")
        if len(timed["cpu_seconds"]) != TRIALS or any(not isfinite(t) or t <= 0 for t in timed["cpu_seconds"]):
            raise ValueError(f"{name}: expected {TRIALS} positive CPU times")
        results.append(dict(case=label, differential=checked["counts"][0],
                            algebraic=checked["counts"][1], powerworld_seconds=seconds,
                            gridkit_seconds=timed["cpu_seconds"],
                            errors=checked["errors"]))
    write_tables(results)


if __name__ == "__main__":
    main()
