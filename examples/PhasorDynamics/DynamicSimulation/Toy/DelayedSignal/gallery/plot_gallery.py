"""Redraw the recorded delay experiments as focused, PNG-only figures."""
import hashlib
import html
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/gridkit-delay-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, NullLocator
import numpy as np

from generate import read_data

HERE = Path(__file__).resolve().parent
FIGURES = HERE / "figures"

# Shared styling for every PNG. Each ordinary view contains a single experiment.
BLUE = "#0072B2"
ORANGE = "#D55E00"
GRAY = "#747B82"
INK = "#25313B"
DPI = 180
STYLE = {
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.labelsize": 10, "axes.titlesize": 10,
    "axes.labelcolor": INK, "text.color": INK,
    "axes.edgecolor": "#A6ADB4", "axes.linewidth": 0.65,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": GRAY, "ytick.color": GRAY,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "legend.frameon": False, "legend.fontsize": 9,
    "grid.color": "#E5E8EB", "grid.linewidth": 0.6,
    "axes.axisbelow": True, "axes.formatter.useoffset": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
}

# Original figure names remain the default view of each family.
FAMILIES = [
    ("ramp", "01_smooth_ramp", "Smooth rising transition", "ramp_0.25",
     "The delayed signal follows the same transition at a later time."),
    ("width", "02_pulse_widths", "Gaussian pulse widths", "width_0.1",
     "Choose a pulse width to see how the accepted steps follow its two edges."),
    ("features", "03_multiple_timescales", "Broad and narrow features", "features_0.15",
     "A broad feature and a narrow pulse share the same source history."),
    ("packet", "04_oscillatory_packet", "Oscillatory packet", "packet_0.1",
     "The oscillatory packet retains its shape as it propagates."),
    ("pulse", "05_adaptive_pulse", "Adaptive pulse propagation", "pulse_0.1",
     "Steps contract near the source pulse and its delayed arrival."),
    ("long", "06_long_history", "Long propagation delay", "long_8",
     "The input and its delayed arrival are shown on the same time axis."),
    ("edges", "08_sharp_edges", "Sharp smooth edges", "sharp_rise",
     "These uncapped runs take steps longer than the delay, then refine near the edges. Edge widths are 10–90% times."),
    ("sine", "09_bipolar_sine", "Bipolar sine waves", "sine_5",
     "Each view shows up to three cycles. At 20 Hz, the 50 ms delay is one full period, so the curves coincide."),
    ("chirp", "10_chirp", "Increasing-frequency waveform", "chirp_0.05",
     "The source frequency rises from 1 to 11 Hz. Each delay is shown separately."),
    ("corners", "11_derivative_corners", "Continuous signals with corners", "triangle_adaptive",
     "Signal values remain continuous. The variants compare corner shapes and source-stencil restarts."),
    ("narrow", "12_narrow_pulses", "Narrow pulses at different offsets", "narrow_2.013",
     "A 1 ms Gaussian pulse is resolved with a maximum step of 0.5 ms."),
]


def duration(value):
    return f"{value * 1000:g} ms" if value < 1 else f"{value:g} s"


def label(run):
    if run["group"] == "sine":
        return f"{run['frequency']:g} Hz"
    if "label" in run:
        return run["label"]
    if run["group"] == "width":
        return f"Width {duration(run['width'])}"
    return f"τ = {duration(run['tau'])}"


def subtitle(run):
    parts = [f"τ = {duration(run['tau'])}"]
    if run["shape"] == "sine":
        parts.append(f"{run['frequency']:g} Hz")
        if np.isclose(run["frequency"] * run["tau"], 1):
            parts.append("Delay equals one period")
    elif run["shape"] == "chirp":
        parts.append(f"{run['frequency']:g} Hz + {run['chirp_rate']:g} Hz/s")
    elif run["shape"] in ("rise", "fall", "asymmetric"):
        parts.append(f"10–90% edge width {duration(run['width'])}" if run["shape"] != "asymmetric"
                     else f"Rise {duration(run['width'])} · fall {duration(run['fall_width'])}")
    elif run["shape"] == "pulse":
        parts.append(f"Gaussian width {duration(run['width'])}")
    if "label" in run and run["group"] != "sine":
        parts.append(run["label"])
    return "  ·  ".join(parts)


def window(run):
    if "zoom" in run:
        return tuple(run["zoom"])
    if run["shape"] == "sine":
        return 0, min(run["tf"], 3 / run["frequency"])
    return 0, run["tf"]


def time_axis(interval):
    """Use milliseconds near short features, with an explicit time origin."""
    if interval[1] - interval[0] < 0.25:
        origin = round(sum(interval) / 2)
        label = f"Time from {origin:g} s (ms)" if origin else "Time (ms)"
        return lambda t: 1000 * (t - origin), label
    return lambda t: t, "Time (s)"


def canvas(title, description, *, equal=False, history=False):
    height = 6.9 if history else 5.3
    ratios = [2.5, 1, 1] if history else ([1, 1] if equal else [2.5, 1])
    fig, axes = plt.subplots(
        len(ratios), 1, sharex=True, figsize=(8, height),
        gridspec_kw={"height_ratios": ratios},
    )
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.64 / height,
                        top=1 - 1.11 / height, hspace=0.30 if history else 0.24)
    fig.text(0.12, 1 - 0.244 / height, title, fontsize=14, weight="semibold", va="top")
    fig.text(0.12, 1 - 0.52 / height, description, fontsize=9, color=GRAY, va="top")
    for index, ax in enumerate(axes):
        ax.grid(axis="y", which="major")
        ax.text(-0.105, 1, f"({chr(97 + index)})", transform=ax.transAxes,
                va="top", color=GRAY, fontsize=9)
        ax.yaxis.set_major_locator(MaxNLocator(4))
    axes[-1].xaxis.set_major_locator(MaxNLocator(5))
    return fig, axes


def legend(fig, ax):
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.11, 1 - 1 / fig.get_figheight()),
               ncol=len(labels), handlelength=2.7, columnspacing=2.4)


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=DPI)
    plt.close(fig)


def plot_ratio(ax, run, interval, transform, color=BLUE, linestyle="-"):
    steps = read_data(run["name"], "steps")
    ax.step(transform(np.r_[0, steps["t"]]),
            np.r_[steps["h_over_tau"][0], steps["h_over_tau"]],
            where="pre", color=color, lw=1.1, ls=linestyle)
    visible = steps[(steps["t"] > interval[0]) & (steps["t"] - steps["h"] < interval[1])]
    upper = max(visible["h_over_tau"])
    if run.get("h_max", 0):
        cap = run["h_max"] / run["tau"]
        ax.axhline(cap, color=GRAY, ls=":", lw=0.9)
        ax.annotate(r"$h = h_{\max}$", (transform(interval[1]), cap), xytext=(-4, 4),
                    textcoords="offset points", ha="right", va="bottom", color=GRAY, fontsize=8,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
        ax.text(1, 1.08, f"$h_{{\\max}}$ = {duration(run['h_max'])} · {run['h_max_reason']}",
                transform=ax.transAxes, ha="right", fontsize=8, color=GRAY)
        upper = max(upper, cap)
    return min(visible["h_over_tau"]), upper


def finish_ratio(ax, limits, interval, transform):
    ax.set_yscale("log")
    lower, upper = min(limits[0], 1) * 0.65, max(limits[1], 1) * 1.6
    exponents = np.arange(np.ceil(np.log10(lower)), np.floor(np.log10(upper)) + 1)
    if len(exponents) > 4:
        exponents = np.unique(np.r_[0, exponents[np.linspace(0, len(exponents) - 1, 3, dtype=int)]])
    ax.set_yticks(10.0 ** exponents)
    ax.set_ylim(lower, upper)
    ax.yaxis.set_minor_locator(NullLocator())
    ax.axhline(1, color=GRAY, ls=(0, (4, 4)), lw=0.8)
    ax.annotate("h = τ", (transform(interval[1]), 1), xytext=(-4, -4),
                textcoords="offset points", ha="right", va="top", color=GRAY, fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    ax.set_ylabel(r"$h/\tau$")
    ax.set_xlim(*transform(np.asarray(interval)))


def plot_storage(ax, steps, interval, transform, column="segments", **style):
    times = np.r_[0, steps["t"]]
    counts = np.r_[0, steps[column]]
    ax.step(transform(times), counts, where="post", lw=1.3, **style)
    first = max(0, np.searchsorted(times, interval[0], side="right") - 1)
    last = np.searchsorted(times, interval[1], side="right")
    return max(counts[first:last])


def finish_storage(ax, maximum):
    ax.set_ylabel("Retained\nsegments")
    ax.set_ylim(0, max(1, maximum * 1.18))
    ax.yaxis.set_major_locator(MaxNLocator(3, integer=True))


def retained_history(ax, run, interval, transform):
    steps = read_data(run["name"], "steps")
    maxima = []
    for column, color, linestyle, label in (
        ("uncompressed_segments", GRAY, "--", "Uncompressed"),
        ("compressed_segments", BLUE, "-", "Compressed"),
    ):
        maxima.append(plot_storage(ax, steps, interval, transform, column,
                                   color=color, ls=linestyle, label=label))
    finish_storage(ax, max(maxima))
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2,
              borderaxespad=0, fontsize=8, handlelength=2.7, columnspacing=2.4)
    ax.text(1, 1.08, "Same accepted knots", transform=ax.transAxes,
            ha="right", fontsize=8, color=GRAY)


def signal_figure(run, name, title):
    interval = window(run)
    transform, xlabel = time_axis(interval)
    samples = read_data(run["name"])
    visible = samples[(samples["t"] >= interval[0]) & (samples["t"] <= interval[1])]
    fig, (waveform, ratio, storage) = canvas(title, subtitle(run), history=True)
    waveform.plot(transform(visible["t"]), visible["u"], color=GRAY, ls=(0, (5, 3)),
                  lw=1.3, label="Input")
    waveform.plot(transform(visible["t"]), visible["z"], color=BLUE, lw=1.7,
                  label="Delayed output")
    waveform.set_ylabel("Signal (a.u.)")
    waveform.margins(y=0.15)
    bounds = plot_ratio(ratio, run, interval, transform)
    finish_ratio(ratio, bounds, interval, transform)
    if not run.get("h_max", 0):
        ratio.text(1, 1.08, "No step cap", transform=ratio.transAxes,
                   ha="right", fontsize=8, color=GRAY)
    retained_history(storage, run, interval, transform)
    storage.set_xlabel(xlabel)
    legend(fig, waveform)
    save(fig, name)


def compression_figure(runs):
    interval = (0, runs[0]["tf"])
    transform, xlabel = time_axis(interval)
    fig, (storage, ratio) = canvas("History compression", "Same input and delay · τ = 2 s", equal=True)
    bounds = []
    for run, color, linestyle in zip(runs, (GRAY, BLUE), ("--", "-")):
        steps = read_data(run["name"], "steps")
        storage.step(steps["t"], steps["segments"], where="post", color=color, ls=linestyle,
                     lw=1.5, label="Compressed" if run["compress"] else "Uncompressed")
        bounds.append(plot_ratio(ratio, run, interval, transform, color, linestyle))
    storage.set_ylabel("Retained segments")
    storage.set_ylim(bottom=0)
    finish_ratio(ratio, (min(b[0] for b in bounds), max(b[1] for b in bounds)), interval, transform)
    ratio.set_xlabel(xlabel)
    legend(fig, storage)
    save(fig, "07_compression")


def resolution_figure(runs):
    interval = tuple(runs[0]["zoom"])
    transform, xlabel = time_axis(interval)
    fig, (waveform, ratio, storage) = canvas("Resolving a narrow source pulse", "Gaussian width 1 ms · τ = 8 s · delayed arrival", history=True)
    bounds = []
    maxima = []
    for run, color, linestyle in zip(runs, (ORANGE, BLUE), ("--", "-")):
        samples = read_data(run["name"])
        visible = samples[(samples["t"] >= interval[0]) & (samples["t"] <= interval[1])]
        waveform.plot(transform(visible["t"]), visible["z"], color=color, lw=1.7, ls=linestyle,
                      label="Uncapped — pulse missed" if run.get("expected_missed")
                      else f"Maximum step {duration(run['h_max'])}")
        bounds.append(plot_ratio(ratio, run, interval, transform, color, linestyle))
        maxima.append(plot_storage(storage, read_data(run["name"], "steps"), interval,
                                   transform, color=color, ls=linestyle))
    # A few reference markers verify the shape without another overlapping line.
    reference = visible[::max(1, len(visible) // 11)]
    waveform.plot(transform(reference["t"]), reference["exact_z"], "o", ms=4,
                  mfc="white", mec=INK, mew=0.8, label="Analytic reference")
    waveform.set_ylabel("Signal (a.u.)")
    waveform.margins(y=0.15)
    finish_ratio(ratio, (min(b[0] for b in bounds), max(b[1] for b in bounds)), interval, transform)
    finish_storage(storage, max(maxima))
    storage.text(1, 1.08, "Actual solver history · compression off", transform=storage.transAxes,
                 ha="right", fontsize=8, color=GRAY)
    storage.set_xlabel(xlabel)
    legend(fig, waveform)
    save(fig, "13_source_resolution")


def history_figure():
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.3))
    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.16, top=0.71, wspace=0.29)
    fig.text(0.10, 0.95, "History reconstruction", fontsize=14, weight="semibold", va="top")
    fig.text(0.10, 0.88, "Accepted source values and the polynomials connecting them", fontsize=9, color=GRAY)
    for index, (name, interval, title) in enumerate((
        ("sharp_rise", (1.999, 2.001), "Smooth edge"),
        ("triangle_restart", (0.96, 1.04), "Derivative corner"),
    )):
        ax = axes[index]
        transform, xlabel = time_axis(interval)
        segments = read_data(name, "segments")
        selected = segments[(segments["end"] > interval[0]) & (segments["begin"] < interval[1])]
        for i, segment in enumerate(selected):
            x = np.linspace(0, 1, 20)
            t = segment["begin"] + x * (segment["end"] - segment["begin"])
            y = ((1 - x) ** 3 * segment["c0"] + 3 * x * (1 - x) ** 2 * segment["c1"]
                 + 3 * x ** 2 * (1 - x) * segment["c2"] + x ** 3 * segment["c3"])
            ax.plot(transform(t), y, color=BLUE, lw=1.6, label="History polynomial" if i == 0 else None)
        ax.scatter(transform(selected["begin"]), selected["c0"], s=18, color="white",
                   edgecolor=BLUE, linewidth=1, label="Accepted value", zorder=3)
        ax.set(xlim=transform(np.asarray(interval)), xlabel=xlabel)
        # Set the range from samples within the displayed window, not full segments.
        samples = read_data(name)
        visible = samples[(samples["t"] >= interval[0]) & (samples["t"] <= interval[1])]
        low, high = min(visible["exact_u"]), max(visible["exact_u"])
        ax.set_ylim(low - 0.1 * (high - low), high + 0.12 * (high - low))
        ax.set_title(f"({chr(97 + index)}) {title}", loc="left", pad=12)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.grid(axis="y")
    axes[0].set_ylabel("Signal (a.u.)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.09, 0.77), ncol=2, handlelength=2.7)
    save(fig, "14_history_detail")


def write_index(entries, runs):
    parts = ['''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Delayed signal gallery · GridKit</title><style>
:root{font-family:system-ui,sans-serif;color:#25313b;background:white;line-height:1.6}
body{max-width:960px;margin:48px auto;padding:0 24px}h1{font-size:30px;line-height:1.2;margin:0 0 12px}
h2{font-size:20px;font-weight:600;margin:0}p{margin:10px 0 20px;color:#626b74}a{color:#0072b2}
header{padding-bottom:24px}nav{display:flex;flex-wrap:wrap;gap:6px 18px;font-size:13px;margin:20px 0}
figure{margin:0;padding:32px 0;border-top:1px solid #e5e8eb}.number{color:#89929a;font-weight:400;margin-right:12px}
.variants{display:flex;flex-wrap:wrap;gap:8px;margin:18px 0 8px}
.variants a{color:#53616c;border:1px solid #d5dce1;padding:4px 12px;border-radius:4px;text-decoration:none;font-size:13px}
.variants a[aria-current=true]{color:#0072b2;border-color:#0072b2;background:#f0f7fb}
.variants a:focus-visible{outline:2px solid #0072b2;outline-offset:3px}
img{width:100%;height:auto;display:block}figcaption{font-size:14px;color:#626b74;margin:4px 0 12px}
.download{font-size:13px}details{border-top:1px solid #e5e8eb;margin-top:28px;padding:18px 0;font-size:13px}
summary{cursor:pointer;font-size:15px}table{border-collapse:collapse;width:100%;margin:18px 0}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid #e5e8eb}th{font-weight:600}
.scroll{overflow-x:auto}@media(max-width:600px){body{padding:0 12px;margin-top:28px}h1{font-size:25px}}
</style></head><body><header><h1>Delayed signal gallery</h1>
<p>Three aligned panels show propagation, the accepted step relative to the delay, h/τ, and retained history. The dashed line at 1 marks h = τ. History panels compare compressed and uncompressed storage on the same accepted source knots.</p>
<p>Steps are uncapped by default. Only the narrow Gaussian source-resolution cases use h_max; the PNGs state its value and purpose and mark it with a dotted line.</p>
<p>Choose a variant beneath each heading. Every view is a standalone PNG.</p><nav>''']
    for entry in entries:
        parts.append(f'<a href="#{entry["name"]}">{html.escape(entry["title"])}</a>')
    parts.append('</nav></header><main>')
    for entry in entries:
        variants = entry["variants"]
        default = variants[0]
        number = entry["name"].split('_')[0]
        parts.append(f'<figure id="{entry["name"]}"><h2><span class="number">{number}</span>{html.escape(entry["title"])}</h2>')
        if len(variants) > 1:
            parts.append('<div class="variants" aria-label="Signal variants">')
            for i, variant in enumerate(variants):
                parts.append(f'<a href="figures/{variant["file"]}.png" data-variant '
                             f'aria-current="{str(i == 0).lower()}" data-label="{html.escape(variant["label"])}">'
                             f'{html.escape(variant["label"])}</a>')
            parts.append('</div>')
        parts.append(f'<img src="figures/{default["file"]}.png" loading="lazy" '
                     f'alt="{html.escape(entry["title"] + ": " + default["label"])}">'
                     f'<figcaption>{html.escape(entry["caption"])}</figcaption>'
                     f'<a class="download" href="figures/{default["file"]}.png">Open PNG</a></figure>')
    parts.append('''</main><details><summary>Numerical results and methods</summary>
<p>All 37 recorded experiments remain available. Error checks, step rejections, native solver retries, and storage statistics are in the data below. The uncapped pulse deliberately demonstrates undersampling; every other run passes its sampled analytic-error check.</p>
<p><a href="summary.csv">CSV summary</a> · <a href="summary.json">JSON summary</a> · <a href="data.zip">Recorded data</a> · <a href="README.md">Methods</a></p>
<div class="scroll"><table><tr><th>Run</th><th>Sampled max error</th><th>History rejections</th><th>Peak segments</th></tr>''')
    for run in runs:
        parts.append(f'<tr><td>{html.escape(run["name"])}</td><td>{run["max_abs_z_error"]:.3g}</td>'
                     f'<td>{run["rejected"]}</td><td>{run["peak_segments"]}</td></tr>')
    parts.append('</table></div></details><details><summary>Application timings</summary>')
    parts.append('<p>Median of three app-reported CPU times with identical inputs.</p>'
                 '<table><tr><th>Case</th><th>Before cleanup</th><th>Current</th></tr>')
    for case in json.loads((HERE / "performance/results.json").read_text())["cases"]:
        parts.append(f'<tr><td>{html.escape(case["case"])}</td><td>{case["cleanup_before_median"]:.4f} s</td>'
                     f'<td>{case["after_median"]:.4f} s</td></tr>')
    parts.append('''</table></details><script>
document.querySelectorAll('[data-variant]').forEach(link => {
  link.addEventListener('click', event => {
    if (event.ctrlKey || event.metaKey || event.shiftKey || event.altKey) return;
    event.preventDefault();
    const figure = link.closest('figure');
    figure.querySelector('img').src = link.href;
    figure.querySelector('img').alt = figure.querySelector('h2').textContent + ': ' + link.dataset.label;
    figure.querySelector('.download').href = link.href;
    figure.querySelectorAll('[data-variant]').forEach(item => item.setAttribute('aria-current', String(item === link)));
  });
});
</script></body></html>''')
    (HERE / "index.html").write_text('\n'.join(parts) + '\n')


def main():
    FIGURES.mkdir(exist_ok=True)
    runs = json.loads((HERE / "summary.json").read_text())
    entries = []
    with plt.rc_context(STYLE):
        for group, name, title, default, caption in FAMILIES:
            selected = [run for run in runs if run["group"] == group]
            selected.sort(key=lambda run: run["name"] != default)
            variants = []
            for i, run in enumerate(selected):
                filename = name if i == 0 else f"{name}__{run['name'].replace('.', 'p')}"
                signal_figure(run, filename, title)
                variants.append({"label": label(run), "file": filename})
            entries.append({"name": name, "title": title, "caption": caption, "variants": variants})
        compression_figure([run for run in runs if run["group"] == "compression"])
        resolution_figure([run for run in runs if run["group"] == "resolution"])
        history_figure()
    for name, title, caption in (
        ("07_compression", "History compression", "Compression reduces retained segments. The second panel shows its effect on the accepted steps; both runs use the same source and delay."),
        ("13_source_resolution", "Resolving a narrow source pulse", "The uncapped solve misses the pulse. A 0.5 ms maximum step resolves it. Reference markers show the analytic delayed waveform. The third panel shows actual retained segments for each solve, both with compression off."),
        ("14_history_detail", "History reconstruction", "Each circle is an accepted source value. The connecting polynomials preserve the continuous value through a marked derivative corner. These are the original segments, before compression."),
    ):
        entries.append({"name": name, "title": title, "caption": caption,
                        "variants": [{"label": title, "file": name}]})
    entries.sort(key=lambda entry: entry["name"])
    write_index(entries, runs)
    # Rendering has its own provenance; keep the original simulation record.
    path = HERE / "provenance.json"
    provenance = json.loads(path.read_text())
    provenance["rendering"] = {
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "data_sha256": hashlib.sha256((HERE / "data.zip").read_bytes()).hexdigest(),
        "matplotlib": matplotlib.__version__, "format": "png", "dpi": DPI,
        "files": [f"figures/{v['file']}.png" for entry in entries for v in entry["variants"]],
    }
    path.write_text(json.dumps(provenance, indent=2) + '\n')
    print(f"Wrote {len(provenance['rendering']['files'])} PNGs in {len(entries)} gallery sections")


if __name__ == "__main__":
    main()
