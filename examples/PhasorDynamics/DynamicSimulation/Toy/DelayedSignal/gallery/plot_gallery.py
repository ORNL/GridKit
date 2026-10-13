"""Draw PNG figures from recorded GridKit output; no simulations run here."""
import html
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/gridkit-delay-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from generate import read_data

HERE = Path(__file__).resolve().parent
FIGURES = HERE / "figures"
FIGURES.mkdir(exist_ok=True)
RUNS = json.loads((HERE / "summary.json").read_text())
COLORS = ["#176b9d", "#d26432", "#25836c"]
PAGES = []
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.16, "grid.linewidth": 0.7,
    "legend.frameon": False, "savefig.facecolor": "white",
    "axes.formatter.useoffset": False,
})


def group(name):
    return [run for run in RUNS if run["group"] == name]


def delay_label(tau):
    return f"τ = {tau * 1000:g} ms" if tau < 1 else f"τ = {tau:g} s"


def layout(title, subtitle, rows=3, columns=1):
    figure, axes = plt.subplots(
        rows, columns, figsize=(11.4, 8.5), squeeze=False, sharex="col"
    )
    figure.suptitle(title, x=0.085, y=0.975, ha="left", fontsize=18, weight="bold")
    figure.text(0.085, 0.934, subtitle, color="#495461", fontsize=11)
    return figure, axes


def waveforms(axis, runs):
    source = read_data(runs[0]["name"])
    axis.plot(source["t"], source["u"], "--", color="#303944", lw=1.8,
              label="Input u(t)", zorder=5)
    for color, run in zip(COLORS, runs):
        values = read_data(run["name"])
        axis.plot(values["t"], values["z"], color=color, lw=1.55,
                  label=delay_label(run["tau"]))
    axis.set_ylabel("Signal value")
    axis.legend(loc="lower right", bbox_to_anchor=(1, 1.02), ncol=4, fontsize=9)


def step_ratios(axis, runs):
    for color, run in zip(COLORS, runs):
        steps = read_data(run["name"], "steps")
        # Each recorded h belongs to the interval ENDING at its recorded t.
        axis.step(steps["t"], steps["h_over_tau"], where="pre", color=color, lw=1.1)
        attempts = read_data(run["name"], "attempts")
        rejected = attempts[attempts["rejected"] > 0]
        axis.scatter(rejected["t"], rejected["h"] / run["tau"],
                     marker="x", s=25, color=color, zorder=4)
    axis.axhline(1, color="#59616b", ls=":", lw=1)
    axis.set(yscale="log", ylabel="Accepted h / τ")
    axis.set_title("Step size relative to delay · dotted line: h = τ · ×: rejected trial", loc="left", fontsize=10)


def storage(axis, runs):
    for color, run in zip(COLORS, runs):
        steps = read_data(run["name"], "steps")
        axis.step(steps["t"], steps["segments"], where="post", color=color, lw=1.4)
    axis.set(ylabel="Stored history segments")
    axis.set_title("Accepted source history retained for future delayed reads", loc="left", fontsize=10)


def save(figure, filename, title, caption):
    figure.tight_layout(rect=(0.015, 0.07, 0.99, 0.98), h_pad=1.7)
    figure.text(0.085, 0.025, "GridKit output · adaptive IDAS · rtol 10⁻⁶, atol 10⁻⁸ · h_min 1 µs · analytic prehistory",
                fontsize=9, color="#586473")
    figure.savefig(FIGURES / f"{filename}.png", dpi=170)
    plt.close(figure)
    PAGES.append((filename, title, caption))


def adaptive_figure(group_name, filename, title, subtitle, caption, memory_first=False):
    runs = group(group_name)
    figure, axes = layout(title, subtitle)
    waveforms(axes[0, 0], runs)
    if memory_first:
        storage(axes[1, 0], runs)
        step_ratios(axes[2, 0], runs)
    else:
        step_ratios(axes[1, 0], runs)
        if group_name == "pulse":
            for color, run in zip(COLORS, runs):
                steps = read_data(run["name"], "steps")
                axes[2, 0].step(steps["t"], steps["h"], where="pre", color=color, lw=1.1)
            axes[2, 0].set(yscale="log", ylabel="Accepted step h (s)")
            axes[2, 0].set_title("The same accepted steps in seconds", loc="left", fontsize=10)
        else:
            storage(axes[2, 0], runs)
    axes[2, 0].set_xlabel("Time (s)")
    for axis in axes[:, 0]:
        axis.set_xlim(0, runs[0]["tf"])
    measured = "; ".join(
        f"{delay_label(run['tau'])}: max h/τ {run['max_h_tau']:.3g}, "
        f"max sampled error {run['max_abs_z_error']:.2g}"
        for run in runs
    )
    save(figure, filename, title, caption + "\n\n" + measured)


adaptive_figure(
    "ramp", "01_smooth_ramp", "A smooth transition, delayed by different amounts",
    "Ramp width 150 ms · delays from 10 ms to 1.5 s",
    "Input u(t) = 1 + 0.3[1 + tanh((t − 2)/0.15)]. The output settles to the same "
    "value at a later time. The solver refines around both the source transition "
    "and its delayed appearance. Compression is off; pruning remains active."
)

figure, axes = layout("The same delay, three different pulse widths",
                      "Delay fixed at 250 ms · only the input width changes", columns=2)
for row, run in enumerate(group("width")):
    values = read_data(run["name"])
    color = COLORS[row]
    axes[row, 0].plot(values["t"], values["u"], "--", color="#303944", lw=1.7, label="Input u(t)")
    axes[row, 0].plot(values["t"], values["z"], color=color, lw=1.5, label="Delayed output")
    axes[row, 0].set(ylabel="Signal value", title=f"Gaussian width w = {run['width'] * 1000:g} ms")
    axes[row, 0].legend(loc="upper right", fontsize=9)
    steps = read_data(run["name"], "steps")
    axes[row, 1].step(steps["t"], steps["h_over_tau"], where="pre", color=color, lw=1.1)
    axes[row, 1].axhline(1, color="#59616b", ls=":", lw=1)
    axes[row, 1].set(yscale="log", ylabel="Accepted h / τ",
                     title=f"{run['accepted']} steps · {run['rejected']} history rejections")
    for axis in axes[row]:
        axis.set_xlim(0, run["tf"])
for axis in axes[-1]:
    axis.set_xlabel("Time (s)")
save(figure, "02_pulse_widths", "The same delay, three different pulse widths",
     "Input u(t) = 1 + 0.6 exp[−((t − 2)/w)²]. Widths are 25, 100, and 400 ms; "
     "w is the Gaussian scale parameter, not its full width at half maximum. "
     "All runs have the same 250 ms delay. The required step size changes with "
     "the signal feature width even though the delay stays fixed.")

adaptive_figure(
    "features", "03_multiple_timescales", "Broad and narrow features in one signal",
    "A broad pulse followed by a sharp, smooth pulse · three delay durations",
    "Input u(t) = 1 + 0.25 exp[−((t − 2)/0.3)²] + 0.6 exp[−((t − 4)/0.04)²]. "
    "The narrow feature requires more temporal resolution. Stored history preserves "
    "it until each delayed output has passed it. Compression is off."
)
adaptive_figure(
    "packet", "04_oscillatory_packet", "An oscillatory packet keeps its shape through a delay",
    "A smooth envelope around oscillations · delays from 5 ms to 800 ms",
    "Input u(t) = 1 + 0.5 exp[−((t − 3)/0.6)²] sin(20(t − 3)). The oscillation "
    "period is about 314 ms. A delay shifts both the envelope and the oscillations. "
    "Compression is off. Steps grow in quiet parts and shrink to resolve the oscillations."
)
adaptive_figure(
    "pulse", "05_adaptive_pulse", "Adaptive steps around a smooth pulse",
    "Gaussian width 120 ms · delays of 1 ms, 100 ms, and 500 ms",
    "Input u(t) = 1 + 0.6 exp[−((t − 2)/0.12)²]. The dashed trace is the actual "
    "GridKit input. Accepted h/τ can exceed one; "
    "accuracy still requires refinement around the pulse. Compression is off."
)
adaptive_figure(
    "long", "06_long_history", "A short feature can require a long history",
    "The same 120 ms pulse · delays of 100 ms, 2 s, and 8 s",
    "The source pulse occurs near 2 s, while the 8 s tap reproduces it near 10 s. "
    "Its history must survive the intervening quiet period. Falling segment counts "
    "show pruning after the feature has left the required delay window. "
    "Compression is off in every run.", memory_first=True
)

runs = group("compression")
figure, axes = layout("History compression around a delayed pulse",
                      "Gaussian width 120 ms · delay 2 s · compression is the only changed option")
source = read_data(runs[0]["name"])
axes[0, 0].plot(source["t"], source["u"], "--", color="#303944", lw=1.8, label="Input u(t)")
for color, run, style in zip(COLORS, runs, ["-", "--"]):
    values = read_data(run["name"])
    steps = read_data(run["name"], "steps")
    label = "Compression on" if run["compress"] else "Compression off"
    axes[0, 0].plot(values["t"], values["z"], style, color=color, lw=1.7, label=label)
    axes[1, 0].step(steps["t"], steps["segments"], where="post", color=color, lw=1.4,
                     label=f"{label}: peak {run['peak_segments']} segments")
    axes[2, 0].plot(values["t"], abs(values["z"] - values["exact_z"]) * 1e6,
                    style, color=color, lw=1.15, label=label)
axes[0, 0].set_ylabel("Signal value")
axes[0, 0].legend(loc="lower right", bbox_to_anchor=(1, 1.02), ncol=3, fontsize=9)
axes[1, 0].set_ylabel("Stored history segments")
axes[1, 0].set_title("Both runs prune history outside the required delay window", loc="left", fontsize=10)
axes[1, 0].legend(loc="upper right", fontsize=9)
axes[2, 0].set(xlabel="Time (s)", ylabel="Absolute delay error (10⁻⁶)")
axes[2, 0].set_title("Compare each output with the analytic input evaluated at t − τ", loc="left", fontsize=10)
for axis in axes[:, 0]:
    axis.set_xlim(0, runs[0]["tf"])
change = 100 * (1 - runs[1]["peak_segments"] / runs[0]["peak_segments"])
save(figure, "07_compression", "History compression around a delayed pulse",
     f"Input u(t) = 1 + 0.6 exp[−((t − 2)/0.12)²]. Peak segment count decreases by "
     f"{change:.1f}%. Accepted steps: {runs[0]['accepted']} off, {runs[1]['accepted']} on. "
     "Compression may change the adaptive step sequence. Segment count is not total "
     "process memory, and a compression-error bound is not a bound on the difference "
     "between complete numerical trajectories.")

parts = ["""<!doctype html><html lang="en"><meta charset="utf-8">
<title>GridKit signal delay gallery</title><style>
body{font:16px/1.6 system-ui,sans-serif;color:#24313e;max-width:1180px;margin:36px auto;padding:0 24px;background:#f5f7fa}
h1{font-size:32px}h2{font-size:23px}a{color:#176b9d}nav{display:flex;gap:14px;flex-wrap:wrap}
figure{margin:30px 0;padding:22px;background:white;border:1px solid #dce2e9;border-radius:8px}
img{width:100%;height:auto}figcaption{color:#465362}
table{border-collapse:collapse}td,th{padding:7px 22px 7px 0;text-align:left;border-bottom:1px solid #dce2e9}
</style><h1>Smooth signals, delayed in GridKit</h1>
<p>Seven PNG figures from 20 runs of the actual GridKit delay implementation.
Every waveform comparison shows the input. Step sizes are accepted internal solver
steps; they are distinct from the monitor interval.</p>
<p><a href="summary.csv">Numerical summary</a> · <a href="README.md">Methods and reproduction</a>
· <a href="performance/results.json">Application timing comparison</a></p><nav>"""]
parts.extend(f'<a href="#{name}">{index + 1}. {html.escape(title)}</a>'
             for index, (name, title, _) in enumerate(PAGES))
parts.append("</nav>")
performance = json.loads((HERE / "performance/results.json").read_text())
parts.append('<h2>Application timings</h2><p>Median of three app-reported CPU times, using identical case settings. '
             'Native retries preserve multistep history; all crossed delay intervals receive an accuracy check.</p>'
             '<table><tr><th>Case</th><th>Original</th><th>Before this change</th><th>Native check</th><th>Time change</th></tr>')
for case in performance["cases"]:
    change = 100 * (case["after_median"] / case["previous_median"] - 1)
    parts.append(f'<tr><td>{html.escape(case["case"])}</td><td>{case["before_median"]:.3f} s</td>'
                 f'<td>{case["previous_median"]:.3f} s</td><td>{case["after_median"]:.3f} s</td><td>{change:+.1f}%</td></tr>')
parts.append('</table>')
for name, title, caption in PAGES:
    parts.append(f'<figure id="{name}"><h2>{html.escape(title)}</h2>'
                 f'<a href="figures/{name}.png"><img src="figures/{name}.png" alt="{html.escape(title)}"></a>'
                 f'<figcaption>{html.escape(caption).replace(chr(10), "<br>")}</figcaption></figure>')
parts.append("</html>")
(HERE / "index.html").write_text("\n".join(parts))
print(f"Wrote {len(PAGES)} PNG figures and index.html")
