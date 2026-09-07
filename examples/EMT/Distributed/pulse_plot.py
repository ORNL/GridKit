"""TeX diagrams, pulse responses, and a local gallery for the pulse study."""

import csv
import hashlib
import html
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from zipfile import ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pulse import AMPLITUDE, GEOMETRY, LENGTH, MODELS, ROOT, SECTIONS, START, WIDTH

DISPLAY_SECTIONS = (1, 5, 20)
DISPLAY_MODELS = ("distributed",) + tuple(f"pi_{n}" for n in DISPLAY_SECTIONS)

COLORS = {
    "distributed": "#202b35",
    "pi_1": "#d55e00",
    "pi_5": "#a64b89",
    "pi_20": "#0072b2",
}
LABELS = {
    "distributed": "Distributed",
    **{f"pi_{n}": f"{n} π sections" for n in SECTIONS},
}
LABELS["pi_1"] = "1 π section"


def save(fig, output, name):
    for extension in ("svg", "pdf", "png"):
        fig.savefig(output / f"{name}.{extension}", dpi=200, facecolor="white")
    plt.close(fig)


def tex_figure(output, template, name):
    """Render editable TeX with the shared README style to PDF, SVG, and PNG."""
    source = Path(__file__).resolve().parent
    style = source / "diagram-style.tex"
    if not style.exists():
        style = ROOT / "docs/Figures/EMT/diagram-style.tex"
    shutil.copyfile(style, output / "diagram-style.tex")
    shutil.copyfile(source / template, output / f"{name}.tex")
    with tempfile.TemporaryDirectory(prefix=f"{name}-", dir=output) as temporary:
        temporary = Path(temporary).resolve()
        commands = [[
            "pdflatex", "-interaction=nonstopmode", "-halt-on-error",
            f"-output-directory={temporary}", f"{name}.tex",
        ]] * 2
        commands += [
            ["dvisvgm", "--pdf", "--no-fonts", f"--output={name}.svg",
             str(temporary / f"{name}.pdf")],
            ["pdftoppm", "-png", "-r", "600", "-singlefile",
             str(temporary / f"{name}.pdf"), name],
        ]
        for command in commands:
            result = subprocess.run(command, cwd=output, capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(
                    f"{' '.join(command)} failed:\n{result.stdout}\n{result.stderr}"
                )
        shutil.copyfile(temporary / f"{name}.pdf", output / f"{name}.pdf")


def geometry(output):
    doc = json.loads(GEOMETRY.read_text())
    catalog = json.loads((GEOMETRY.parent / doc["include"][0]).read_text())
    attachments = catalog["towers"][doc["tower"]]["attachments"]
    span = doc["path"]["span"]
    values, macros = [], {}
    x = np.linspace(0, span, 400)
    for wire in doc["conductors"]:
        at = attachments[wire["at"]]
        material = catalog["conductors"][wire["type"]]
        a = wire["tension"] / material["weight"]
        eta = span / (2 * a)
        mean = at["h"] - a * (np.cosh(eta) - np.sinh(eta) / eta)
        middle = at["h"] - a * (np.cosh(eta) - 1)
        values.append({
            "attachment": wire["at"], "phase": wire["phase"], "x_m": at["x"],
            "attachment_height_m": at["h"], "mean_height_m": mean,
            "midspan_height_m": middle, "sag_m": at["h"] - middle,
        })
        group = "shield" if wire["phase"] == "g" else "phase"
        if f"{group}Mean" in macros:
            continue
        macros.update({
            f"{group}Mean": f"{mean:.12g}", f"{group}Mid": f"{middle:.12g}",
            f"{group}MeanLabel": f"{mean:.2f}",
            f"{group}SagLabel": f"{at['h'] - middle:.2f}",
        })
        height = middle + a * (np.cosh((x - span / 2) / a) - 1)
        np.savetxt(output / f"{group}-span.dat", np.c_[x, height], fmt="%.10g")
    (output / "pulse-geometry-data.tex").write_text(
        "".join(r"\newcommand{\%s}{%s}" % item + "\n" for item in macros.items())
    )
    (output / "geometry.json").write_text(json.dumps(values, indent=2) + "\n")
    tex_figure(output, "pulse_geometry.tex", "01_geometry")


def circuits(output):
    macros = {
        "studyLength": f"{LENGTH / 1000:g}",
        "circuitRows": ",".join(
            f"{n}/{LENGTH / n / 1000:g}/{row}"
            for row, n in enumerate(DISPLAY_SECTIONS)
        ),
    }
    (output / "pulse-circuit-data.tex").write_text(
        "".join(r"\newcommand{\%s}{%s}" % item + "\n" for item in macros.items())
    )
    tex_figure(output, "pulse_circuits.tex", "02_circuits")


def responses(output, data):
    style = {
        "font.family": "STIXGeneral", "font.size": 12,
        "mathtext.fontset": "stix", "axes.labelsize": 12,
        "axes.titlesize": 13, "axes.linewidth": 0.7,
        "svg.fonttype": "none", "pdf.fonttype": 42,
    }
    with plt.rc_context(style):
        fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.3), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.085, right=0.98, bottom=0.10, top=0.87, hspace=0.26)
        t = (data["time_s"] - START) * 1e6
        source, = axes[0].plot(
            [0, 0, WIDTH * 1e6, WIDTH * 1e6, 1000],
            [0, AMPLITUDE, AMPLITUDE, 0, 0],
            color="#777777", lw=1.1, ls=(0, (5, 3)),
            label=r"Source voltage $e(t)$", zorder=2,
        )
        handles = {}
        for ax, name, title in zip(
            axes, ("sending_V", "receiving_V"),
            ("(a) Sending end", "(b) Receiving end"),
        ):
            for kind in (*DISPLAY_MODELS[1:], "distributed"):
                line, = ax.plot(
                    t, data[f"emt_{kind}_{name}"], color=COLORS[kind],
                    lw=2.0 if kind == "distributed" else 1.5,
                    label="Distributed" if kind == "distributed" else rf"$n={kind[3:]}$",
                    zorder=4 if kind == "distributed" else 3,
                )
                handles[kind] = line
            ax.spines[["top", "right"]].set_visible(False)
            ax.spines[["left", "bottom"]].set_color("#555555")
            ax.grid(axis="y", color="#e5e5e5", linewidth=0.55)
            ax.set_axisbelow(True)
            ax.axhline(0, color="#777777", lw=0.65, zorder=1)
            ax.set_title(title, loc="left", pad=10)
            ax.set(xlim=(0, 1000), ylim=(-150, 1050), ylabel="Voltage [V]")
            ax.set_xticks(np.arange(0, 1001, 200))
            ax.set_yticks(np.arange(0, 1001, 250))
            ax.tick_params(direction="out", length=3, width=0.7, pad=5)
        delay = LENGTH / 299792458 * 1e6
        axes[1].axvline(delay, color="#777777", ls=":", lw=1.0, zorder=1)
        axes[1].text(
            delay + 12, 0.93, rf"$\tau = {delay:.2f}\,\mu\mathrm{{s}}$",
            transform=axes[1].get_xaxis_transform(), color="#555555", fontsize=11,
        )
        axes[1].set_xlabel("Time after pulse launch [µs]", labelpad=8)
        fig.text(0.085, 0.96, "60 km line · 100 µs pulse", fontsize=15, va="center")
        legend = [handles[kind] for kind in DISPLAY_MODELS] + [source]
        fig.legend(
            legend, [line.get_label() for line in legend], ncol=5,
            loc="center right", bbox_to_anchor=(0.98, 0.96), borderaxespad=0,
            frameon=False, fontsize=10.5, handlelength=1.7,
            columnspacing=1.0, handletextpad=0.6,
        )
        save(fig, output, "03_pulse_response")


def report(output, summary):
    entries = [
        (
            "01_geometry",
            "Line geometry",
            "Illustrative twin-Drake geometry. The cross-section shows the sag-averaged wire heights used to calculate the parameters; the span shows the physical catenaries and their sag below the attachments. Cross-section distances share one scale; the span and bundle detail use separate horizontal scales. The electrical study uses a uniform, cyclically averaged 60 km line.",
        ),
        (
            "02_circuits",
            "Line discretization",
            "The n = 1, 5, and 20 circuits share one horizontal scale: Δx = 60/n km. Series symbols denote the fitted frequency-domain impedance: R(ω) = Δx Re Z′₀(jω), L(ω) = Δx Im Z′₀(jω)/ω. Each end shunt is C/2 and each interior shunt is C, with C = C′₀Δx. Bottom nodes denote the common reference potential. Equal phase excitation selects this zero-sequence equivalent; the source and load are omitted.",
        ),
        (
            "03_pulse_response",
            "Sending and receiving pulses",
            "The dashed trace is the prescribed source voltage e(t): a 1 kV, 100 µs square pulse applied through 100 Ω per phase. The source current is i(t) = [e(t) − vₛ(t)] / 100 Ω. Solid traces show the computed sending and receiving voltages for the distributed model and n = 1, 5, and 20 lumped sections. Equal phase excitation makes the three phase voltages coincide. The dotted marker denotes explicit transport delay.",
        ),
    ]
    rows, records = [], []
    for kind in MODELS:
        m = summary["models"][kind]
        row = [
            LABELS[kind],
            f"{m['peak_V']:.3f}",
            f"{m['peak_time_after_launch_us']:.2f}",
            f"{m['waveform_error_rms_V']:.6g}",
            f"{m['waveform_error_percent']:.4g}",
        ]
        if kind in DISPLAY_MODELS:
            rows.append("<tr>" + "".join(f"<td>{value}</td>" for value in row) + "</tr>")
        records.append(row)
    with (output / "metrics.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "model",
                "peak_V",
                "peak_time_after_launch_us",
                "rms_error_V",
                "rms_error_percent",
            ]
        )
        writer.writerows(records)
    max_error = max(
        m["numerical_checks"]["receiving_V"]["fine_reference_rms"]
        for m in summary["models"].values()
    )
    max_refine = max(
        m["numerical_checks"]["receiving_V"]["tolerance_refinement_rms"]
        for m in summary["models"].values()
    )
    max_inverse = max(
        v
        for k, v in summary["reference_refinement_rms"].items()
        if k.endswith("receiving_V")
    )
    quality = json.loads((output / "coefficients/fit_quality.json").read_text())
    max_fit = 100 * max(m["relative_max"] for m in quality.values())
    evidence = (
        f"All five circuits completed at both tolerance levels (10 simulations). "
        f"Maximum fine-run receiving-voltage RMS error against each circuit’s independent Laplace reference: {max_error:.4g} V. "
        f"Maximum nominal-to-fine RMS difference: {max_refine:.4g} V. "
        f"Halving inverse-transform spacing from 0.1 to 0.05 µs changed receiving references by at most {max_inverse:.4g} V RMS. "
        f"The series-impedance fit has at most {max_fit:.4g}% relative error at the geometry sample frequencies (0.01 Hz–10 MHz); shunt capacitance is represented exactly."
    )
    if all("holdout_relative_max" in m for m in quality.values()):
        held = 100 * max(m["holdout_relative_max"] for m in quality.values())
        evidence += f" Fresh geometry calculations at {quality['zero']['holdout_samples']} intervening frequencies give a maximum series-fit error of {held:.4g}%."
    limitations = (
        "An ideal square pulse has unlimited bandwidth; the geometry fit covers 0.01 Hz–10 MHz. "
        "Numerical reference checks exclude 1 µs around the two source discontinuities; their errors are recorded separately. "
        "Spatial RMS errors include the entire 2 ms response window, including the wavefronts. "
        "The new series fits use positive RL terms and positive modal inductances; the inherited distributed fit has sampled passivity checks, not an all-frequency proof. "
        "Peak times are measured on the 0.5 µs output grid. These illustrative conductor and earth data are not manufacturer or field measurements."
    )

    def source_link(name):
        if name == "03_pulse_response":
            return ""
        return f' · <a href="{name}.tex">TeX source</a>'

    figures = "\n".join(
        f'<section><h2>{i}. {title}</h2><img src="{name}.svg" alt="{html.escape(caption)}"><p>{caption}</p><p class="links"><a href="{name}.svg">SVG</a> · <a href="{name}.pdf">PDF</a> · <a href="{name}.png">PNG</a>{source_link(name)}</p></section>'
        for i, (name, title, caption) in enumerate(entries, 1)
    )
    page = f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>60 km line · square-pulse study</title><style>
body{{font:16px/1.6 system-ui,sans-serif;color:#202b35;background:#f5f7f8;margin:0}}main{{max-width:1100px;margin:auto;padding:40px 24px}}
h1{{font-size:32px;line-height:1.2;margin-bottom:12px}}h2{{font-size:22px}}section{{background:white;border:1px solid #dde3e7;border-radius:8px;padding:24px;margin:24px 0}}
img{{width:100%;height:auto}}a{{color:#006c9e}}p{{max-width:95ch}}.meta,.links{{color:#64717b;font-size:14px}}table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}
th,td{{padding:9px 12px;text-align:right;border-bottom:1px solid #e2e7ea}}th:first-child,td:first-child{{text-align:left}}.scroll{{overflow-x:auto}}
</style><main><p class="meta">GRIDKIT · EMT EXAMPLE</p><h1>A short pulse through a 60 km overhead line</h1>
<p>Distributed propagation and progressively finer frequency-dependent π circuits, derived from the same overhead geometry.</p>
<p class="meta">1 kV source pulse · 100 µs duration · zero-sequence excitation · 100 Ω source · 600 Ω load · zero initial history</p>
{figures}<section><h2>Response and convergence</h2><div class="scroll"><table><thead><tr><th>Model</th><th>Peak [V]</th><th>Peak time [µs]</th><th>RMS error [V]</th><th>RMS error [%]</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
<p class="meta">Peak time is relative to pulse launch. RMS error compares receiving voltage with the independent distributed reference over 0–2 ms after launch; percent error is normalized by that reference’s RMS. The continuum coefficient-fit difference is {summary["continuum_fit_difference_rms_V"]:.4g} V RMS.</p></section>
<section><h2>Validation</h2><p>{evidence}</p><p>{limitations}</p>
<p>Adaptive IDA with sparse KLU; nominal tolerances 10⁻⁷ / 10⁻⁹ and fine tolerances 10⁻⁸ / 10⁻¹⁰. Monitor spacing is not a fixed solver step.</p>
<p class="links"><a href="waveforms.csv">Waveforms CSV</a> · <a href="metrics.csv">Metrics CSV</a> · <a href="validation.json">Validation details</a> · <a href="coefficients/fit_quality.json">Fit quality</a> · <a href="provenance.json">Simulation input hashes</a> · <a href="figure_provenance.json">Figure source hashes</a> · <a href="diagram_sources.zip">Editable diagrams (TeX + data)</a></p></section></main></html>"""
    (output / "index.html").write_text(page)
    markdown = (
        "# 60 km square-pulse study\n\n[Open the three-figure gallery](index.html). [Editable diagrams: TeX and data](diagram_sources.zip).\n\n"
    )
    for name, title, caption in entries:
        markdown += f"## {title}\n\n![{title}]({name}.png)\n\n{caption}\n\n"
    markdown += "## Response and convergence\n\n| Model | Peak [V] | Peak time [µs] | RMS error [V] | RMS error [%] |\n| --- | ---: | ---: | ---: | ---: |\n"
    markdown += "".join(
        "| " + " | ".join(row) + " |\n"
        for kind, row in zip(MODELS, records) if kind in DISPLAY_MODELS
    )
    markdown += f"\nPeak time is relative to launch. RMS errors use the complete 2 ms response window. The continuum coefficient-fit difference is {summary['continuum_fit_difference_rms_V']:.4g} V RMS.\n\n"
    markdown += "## Validation\n\n" + evidence + "\n\n" + limitations + "\n"
    (output / "REPORT.md").write_text(markdown)
    diagram_sources = [
        "01_geometry.tex", "02_circuits.tex", "diagram-style.tex",
        "pulse-geometry-data.tex", "pulse-circuit-data.tex",
        "phase-span.dat", "shield-span.dat",
    ]
    with ZipFile(output / "diagram_sources.zip", "w") as bundle:
        for name in diagram_sources:
            bundle.write(output / name, name)
    sources = [
        Path(__file__).resolve(), GEOMETRY,
        GEOMETRY.parent / "north-american.catalog.json",
        *(output / name for name in diagram_sources),
        *(output / name for name in (
            "waveforms.npz",
            "validation.json", "coefficients/fit_quality.json", "coefficients/yp.json",
        )),
    ]
    (output / "figure_provenance.json").write_text(json.dumps({
        "files_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "authoritative_diagram_format": "TeX",
        "shared_style": "docs/Figures/EMT/diagram-style.tex",
    }, indent=2) + "\n")


def diagrams(output, summary):
    """Redraw the diagrams and report using retained validation and waveforms."""
    geometry(output)
    circuits(output)
    report(output, summary)


def plots(output, summary):
    data = np.load(output / "waveforms.npz")
    geometry(output)
    circuits(output)
    responses(output, data)
    report(output, summary)
