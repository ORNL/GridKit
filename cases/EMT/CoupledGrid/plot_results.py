#!/usr/bin/env python3
"""Render measured EMT mu comparisons from a run manifest and monitor CSVs."""

import argparse
import csv
import html
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/gridkit-coupledgrid-matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLORS = ("#0072B2", "#D55E00", "#8E5BB7")
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .2,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 180, "svg.fonttype": "none"})


def read_csv(path):
    with path.open(newline="") as stream:
        header = next(csv.reader(stream))
    values = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if not np.isfinite(values).all() or np.any(np.diff(values[:, 0]) < 0):
        raise ValueError(f"Nonfinite output or unordered time samples: {path}")
    # Event output may repeat a timestamp; retain the post-event row.
    values = values[np.r_[np.diff(values[:, 0]) > 0, True]]
    return values[:, 0], {key: values[:, i] for i, key in enumerate(header)}


def phases(data, device, stem):
    prefix = f"{device['class']}_{device['id']}_"
    return np.column_stack([data[prefix + stem + phase] for phase in "abc"])


def monitor(data, device, name):
    return data[f"{device['class']}_{device['id']}_{name}"]


def cycle_mean(t, value, frequency):
    """Trailing trapezoidal mean over one period; omit incomplete windows."""
    integral = np.r_[0., np.cumsum(np.diff(t) * (value[:-1] + value[1:]) / 2)]
    mean = frequency * (integral - np.interp(t - 1 / frequency, t, integral))
    mean[t < t[0] + 1 / frequency] = np.nan
    return mean


def powers(voltage, current):
    """Three-phase P and instantaneous alpha/beta Q; positive injection."""
    va = (2 * voltage[:, 0] - voltage[:, 1] - voltage[:, 2]) / 3
    vb = (voltage[:, 1] - voltage[:, 2]) / np.sqrt(3)
    ia = (2 * current[:, 0] - current[:, 1] - current[:, 2]) / 3
    ib = (current[:, 1] - current[:, 2]) / np.sqrt(3)
    return np.sum(voltage * current, axis=1), 1.5 * (vb * ia - va * ib)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    root = args.results.resolve()
    manifest = json.loads((args.manifest or root / "summary.json").read_text())
    case = json.loads(args.case.read_text())
    devices = case["devices"]
    buses = {d["id"]: d for d in devices if d["class"] == "Bus"}
    first = lambda cls: next(d for d in devices if d["class"] == cls)
    machine, ibr, pwm, converter, line, load = [first(cls) for cls in
        ("Machine", "DependentVoltageSource", "PWM", "Converter", "LineLumped", "LoadZ")]
    m_bus, i_bus = (buses[d["inputs"]["bus"]] for d in (machine, ibr))
    ibr_name = "IBR " + i_bus["id"].rsplit("_", 1)[-1]
    machine_name = "Machine " + m_bus["id"].rsplit("_", 1)[-1]
    remote = buses[manifest.get("remote_bus", "bus_10")]
    frequency = manifest.get("frequency_hz", 60.)
    carrier = manifest.get("switching_frequency_hz", 900.)
    base_voltage = manifest.get("voltage_ll_v", 13800.)
    events = manifest.get("events", [])
    event_time = events[0]["time"] if events else 0.
    startup_end = manifest.get("startup_end_s", .2)
    plots = root / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    runs = sorted(manifest["runs"], key=lambda run: run["mu"])
    gallery = []
    measured = []
    for color, run in zip(COLORS, runs):
        t, data = read_csv(root / run["csv"])
        run.update(t=t, data=data, color=color, label=rf"$\mu={run['mu']:,.6g}\;\mathrm{{s}}^{{-1}}$")
        run["voltage"] = {bus: phases(data, device, "v") for bus, device in buses.items()}
        run["rms"] = {bus: np.sqrt(np.maximum(0., cycle_mean(t,
            np.mean((v - np.roll(v, -1, axis=1)) ** 2, axis=1), frequency))) / base_voltage
            for bus, v in run["voltage"].items()}
        current = phases(data, ibr, "i")
        run["ibr_current"] = current
        run["ibr_pq"] = [cycle_mean(t, x, frequency) / 1e6
                         for x in powers(run["voltage"][i_bus["id"]], current)]

    def save(fig, name, title, description, time_axis=True, window=None):
        if time_axis:
            for ax in fig.axes:
                ax.axvspan(0., startup_end, color="#8796a6", alpha=.1, lw=0.)
                for event in events:
                    ax.axvline(event["time"], color="#58636d", ls=":", lw=.8)
                if ax.get_subplotspec().is_last_row():
                    ax.set_xlabel("Time [s]")
                ax.ticklabel_format(axis="y", style="plain", useOffset=False)
                if window:
                    ax.set_xlim(*window)
        handles, labels = next((ax.get_legend_handles_labels() for ax in fig.axes
                                if ax.get_legend_handles_labels()[0]), ([], []))
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5, .955),
                   ncol=3, frameon=False)
        fig.suptitle(title, fontsize=15, y=.995)
        fig.tight_layout(rect=(0, .025, 1, .925))
        fig.text(.5, .012, description, fontsize=9, ha="center", color="#46515d")
        for extension in ("png", "svg", "pdf"):
            fig.savefig(plots / f"{name}.{extension}")
        plt.close(fig)
        gallery.append((name, title, description))

    fig, axs = plt.subplots(4, 2, figsize=(14, 12), sharex=True)
    for run in runs:
        t, data = run["t"], run["data"]
        series = [run["rms"][d["id"]] for d in (m_bus, i_bus, remote)]
        series += [monitor(data, machine, "omega") * frequency,
                   monitor(data, machine, "efd"), run["ibr_pq"][0], run["ibr_pq"][1],
                   cycle_mean(t, monitor(data, machine, "p"), frequency) / 1e6]
        labels = [f"{d['id']}: line–line RMS [pu]" for d in (m_bus, i_bus, remote)]
        labels += [f"{machine_name}: electrical frequency [Hz]", f"{machine_name} field voltage [pu]",
                   f"{ibr_name} injected P [MW]", f"{ibr_name} mean Clarke Q [Mvar]", f"{machine_name} injected P [MW]"]
        for ax, value, label in zip(axs.flat, series, labels):
            ax.plot(t, value, color=run["color"], label=run["label"], lw=1.15)
            ax.set_ylabel(label)
    save(fig, "response_comparison", "Coupled EMT grid · response to the same event",
         "Same circuit and inputs. Gray: initial startup interval (0–0.2 s), not a convergence claim. Dotted: events. RMS/P/Q: trailing 60 Hz cycle.")

    fig, axs = plt.subplots(3, 1, figsize=(14, 9), sharex=True)
    reference = runs[-1]
    reference_series = [reference["rms"][remote["id"]],
                        monitor(reference["data"], machine, "omega") * frequency,
                        reference["ibr_pq"][0]]
    for run in runs:
        series = [run["rms"][remote["id"]], monitor(run["data"], machine, "omega") * frequency,
                  run["ibr_pq"][0]]
        for ax, value, ref, label in zip(axs, series, reference_series,
                [f"{remote['id']}: ΔRMS voltage [pu]", f"{machine_name} Δfrequency [Hz]", f"{ibr_name} ΔP [MW]"]):
            difference = value - np.interp(run["t"], reference["t"], ref)
            ax.plot(run["t"], difference, color=run["color"], label=run["label"], lw=1.2)
            ax.set_ylabel(label)
    save(fig, "response_differences", f"Response differences relative to μ = {reference['mu']:g} s⁻¹",
         "Highest μ is a comparison baseline. Differences include changed PWM excitation, not only numerical error. Gray: initial startup interval.")

    t_end = min(run["t"][-1] for run in runs)
    event_window = (max(0., event_time - 2 / frequency), min(t_end, event_time + 3 / frequency))
    fig, axs = plt.subplots(3, 2, figsize=(14, 10), sharex=True)
    for run in runs:
        t, data = run["t"], run["data"]
        series = [run["voltage"][d["id"]][:, 0] / 1000 for d in (m_bus, i_bus, remote)]
        series += [run["ibr_current"][:, 0], monitor(data, line, "i12a"), monitor(data, load, "ia")]
        labels = [f"{d['id']}: phase-a voltage [kV]" for d in (m_bus, i_bus, remote)]
        labels += [f"{ibr_name} phase-a injected current [A]", f"{line['id']}: phase-a series current [A]",
                   f"{load['id']}: phase-a current [A]"]
        mask = (t >= event_window[0]) & (t <= event_window[1])
        for ax, value, label in zip(axs.flat, series, labels):
            ax.plot(t[mask], value[mask], color=run["color"], label=run["label"], lw=1.)
            ax.set_ylabel(label)
    save(fig, "event_waveforms", "Event waveforms · all μ values overlaid",
         "Raw samples. Positive source/load current injects into its bus; line current flows from its first bus to its second.", window=event_window)

    zoom_end = max(5 / carrier, event_time - 1 / frequency) if events else t_end - 1 / frequency
    zoom_end = min(zoom_end, t_end)
    zoom = (max(0., zoom_end - 5 / carrier), zoom_end)
    fig, axs = plt.subplots(4, 1, figsize=(14, 10), sharex=True)
    for run in runs:
        t, data = run["t"], run["data"]
        mask = (t >= zoom[0]) & (t <= zoom[1])
        series = [monitor(data, pwm, "sa"), monitor(data, converter, "voa") / 1000,
                   run["ibr_current"][:, 0], run["voltage"][i_bus["id"]][:, 0] / 1000]
        for ax, value, label in zip(axs, series, ["PWM $s_a$ [1]",
            "Bridge $v_{oa}$ [kV]", f"{ibr_name} $i_a$ [A]", f"{i_bus['id']} $v_a$ [kV]"]):
            ax.plot(t[mask], value[mask], color=run["color"], label=run["label"], lw=1.3)
            ax.set_ylabel(label)
    save(fig, "switching_waveforms", f"{ibr_name} switching detail · five {carrier:g} Hz carrier periods",
         "Identical DC voltage and modulation. μ changes the smooth PWM approximation and therefore the applied bridge waveform.", window=zoom)

    fig, axs = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    last_event = max((event["time"] for event in events), default=0.)
    for run in runs:
        t, data = run["t"], run["data"]
        dt = float(np.median(np.diff(t)))
        if np.max(np.abs(np.diff(t) - dt)) > dt * 1e-5:
            raise ValueError("Spectrum requires a uniform monitor interval")
        start = max(last_event + 1 / frequency, t_end - 6 / frequency)
        n = int(round((t_end - start) * frequency))
        mask = (t >= t_end - n / frequency) & (t < t_end - dt / 4)
        if n < 1 or mask.sum() < 16:
            raise ValueError("At least one post-event cycle is needed for the spectrum")
        signals = [monitor(data, pwm, "sa")[mask], monitor(data, converter, "voa")[mask], run["ibr_current"][mask, 0]]
        for ax, signal, label in zip(axs, signals, ["PWM peak [1]", "Bridge peak [V]", f"{ibr_name} current peak [A]"]):
            window = np.hanning(len(signal))
            amplitude = 2 * np.abs(np.fft.rfft(signal * window)) / window.sum()
            amplitude[0] *= .5
            if len(signal) % 2 == 0:
                amplitude[-1] *= .5
            hz = np.fft.rfftfreq(len(signal), dt)
            ax.semilogy(hz, np.maximum(amplitude, 1e-12), color=run["color"], label=run["label"], lw=1.1)
            ax.set_ylabel(label)
            ax.set_xlabel("Frequency [Hz]")
            ax.set_xlim(0, min(5 * carrier, .5 / dt))
            ax.axvline(carrier, color="#58636d", ls=":", lw=.8)
    save(fig, "switching_spectrum", f"{ibr_name} · PWM, bridge and current spectra",
         "Qualitative final-window spectra: Hann window with coherent-gain correction. Dotted: carrier. Exact harmonic amplitudes are in metrics.json.", time_axis=False)

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    for j, run in enumerate(runs):
        samples = np.atleast_1d(run["complete_seconds"]).astype(float)
        axs[0].bar(j, np.median(samples), color=run["color"], width=.6, alpha=.85)
        axs[0].scatter(j + np.linspace(-.12, .12, len(samples)), samples, color="#23313a", s=24, zorder=4)
        axs[0].text(j, max(samples), f"{np.median(samples):.3g} s", ha="center", va="bottom", fontsize=10)
        stats = run.get("stats", {})
        counts = [stats.get(key, np.nan) for key in ("steps", "residual_evals", "nonlinear_iters")]
        axs[1].bar(np.arange(3) + (j - 1) * .24, counts, width=.24, color=run["color"], label=run["label"])
        measured.append({"mu": run["mu"], "median_complete_seconds": float(np.median(samples)),
                         "all_complete_seconds": samples.tolist(), "samples": len(run["t"]),
                         "monitor_dt_s": float(np.median(np.diff(run["t"]))),
                         "stats": stats})
    axs[0].set_xticks(range(len(runs)), [f"{run['mu']:,.3f}" for run in runs])
    axs[0].set_xlabel("μ [s⁻¹]")
    axs[0].set_ylabel("App-reported Complete in [s]")
    axs[0].set_ylim(top=axs[0].get_ylim()[1] * 1.14)
    axs[1].set_xticks(range(3), ["IDA steps", "Residual evaluations", "Nonlinear iterations"])
    axs[1].set_ylabel("Solver work [count]")
    save(fig, "runtime_comparison", "Runtime and solver work",
         "Bars show median timings; dots show individual trials. Counts describe the retained waveform run.", time_axis=False)

    (plots / "plot_metadata.json").write_text(json.dumps(measured, indent=2) + "\n")
    cards = "".join(f'<figure><a href="{name}.svg"><img src="{name}.png" alt="{html.escape(title)}"></a>'
                    f'<figcaption><b>{html.escape(title)}</b><p>{html.escape(description)}</p>'
                    f'<a href="{name}.pdf">PDF</a> · <a href="{name}.svg">SVG</a></figcaption></figure>'
                    for name, title, description in gallery)
    (plots / "index.html").write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1"><title>Coupled EMT grid: μ comparison</title>'
        '<style>body{font:16px system-ui;color:#203340;background:#f5f7fa;margin:2rem auto;max-width:1500px;padding:0 1rem}'
        'figure{background:white;border:1px solid #dce2e8;margin:1.5rem 0;padding:1rem}img{width:100%}'
        'a{color:#0072b2}figcaption p{margin:.5rem 0}</style><h1>Coupled EMT grid: μ comparison</h1>'
        '<p>Figures use actual monitor CSVs. Each color denotes the same μ throughout. '
        'The circuit, source commands and solver tolerances are held fixed.</p>' + cards + '</html>')
    run_rows = "".join(f'<tr><td><span style="color:{run["color"]}">●</span> {run["mu"]:,.6g}</td>'
                       f'<td>{np.median(run["complete_seconds"]):.4g}</td>'
                       f'<td>{np.median(run.get("wall_seconds", [float("nan")])):.4g}</td>'
                       f'<td>{run.get("stats", {}).get("steps", "—")}</td>'
                       f'<td><a href="{html.escape(run["csv"])}">Monitor CSV</a></td></tr>' for run in runs)
    root_cards = "".join(f'<figure id="{name}"><a href="plots/{name}.svg"><img src="plots/{name}.png" '
                        f'alt="{html.escape(title)}"></a><figcaption><b>{html.escape(title)}</b>'
                        f'<p>{html.escape(description)}</p><a href="plots/{name}.pdf">PDF</a> · '
                        f'<a href="plots/{name}.svg">SVG</a></figcaption></figure>' for name, title, description in gallery)
    event_text = "; ".join(f'{event["time"]:g} s: {html.escape(event.get("label", event.get("type", "event")))}'
                           for event in events)
    (root / "index.html").write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1"><title>14-bus coupled EMT study</title>'
        '<style>body{font:16px system-ui;color:#203340;background:#f5f7fa;margin:2rem auto;max-width:1500px;padding:0 1rem;line-height:1.5}'
        'figure{background:white;border:1px solid #dce2e8;margin:1.5rem 0;padding:1rem}img{width:100%}'
        'a{color:#0072b2}figcaption p{margin:.5rem 0}table{border-collapse:collapse;background:white;width:100%}'
        'td,th{border:1px solid #dce2e8;padding:.65rem;text-align:left}nav{padding:1rem;background:#eaf0f5}'
        '.pair{display:grid;grid-template-columns:1fr 1fr;gap:1rem}.pair figure{margin:.5rem 0}'
        '@media(max-width:800px){.pair{display:block}table{font-size:.85rem}}</style>'
        '<h1>14-bus coupled EMT study</h1><p>13.8 kV · 60 Hz · 13 physical coupled line sections · '
        '2 synchronous machines with TGOV1 and IEEET1 · 3 open-loop PWM converter sources.</p>'
        '<nav><a href="report.md">Measured results and interpretation</a> · <a href="metrics.json">Measured metrics</a> · '
        '<a href="summary.json">Run metadata</a> · '
        '<a href="#response_comparison">Response overlays</a> · <a href="#switching_waveforms">Switching detail</a> · '
        '<a href="#runtime_comparison">Runtime</a></nav>'
        f'<p>The three runs retain identical network parameters, DC voltage, modulation and initial state. Events: {event_text}. '
        'The gray 0–0.2 s interval marks initial startup; it does not assert subsequent equilibrium. '
        'Changing μ changes the PWM excitation. The highest μ is a comparison baseline, not an exact reference.</p>'
        '<table><tr><th>μ [s⁻¹]</th><th>Median Complete in [s]</th><th>Median wall time [s]</th><th>IDA steps</th><th>Waveform data</th></tr>'
        + run_rows + '</table><figure><a href="one_line.svg"><img src="one_line.png" alt="14-bus network one-line"></a>'
        '<figcaption><b>Actual case topology</b> · <a href="one_line.pdf">PDF</a> · <a href="one_line.svg">SVG</a></figcaption></figure>'
        '<div class="pair"><figure><a href="line_geometry_and_coupling.svg"><img src="line_geometry_and_coupling.png" alt="Physical geometry and full phase matrices"></a>'
        '<figcaption>Physical conductor geometry and full phase coupling · <a href="line_geometry_and_coupling.pdf">PDF</a></figcaption></figure>'
        '<figure><a href="line_frequency_response.svg"><img src="line_frequency_response.png" alt="Physical and fitted self and mutual line impedance"></a>'
        '<figcaption>Physical self/mutual line response and rational fit · <a href="line_frequency_response.pdf">PDF</a></figcaption></figure></div>'
        + root_cards + '</html>')
    print(f"Rendered {len(gallery)} comparison figures to {plots}")


if __name__ == "__main__":
    main()
