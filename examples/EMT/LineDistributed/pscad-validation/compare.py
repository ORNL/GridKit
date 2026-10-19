"""Compare the GridKit and PSCAD recordings of every scenario; write plots and the report.

For each line, model and scenario the PSCAD channels are interpolated onto
GridKit's output times and compared on the study window [0, tf]:

- peak error: max |GridKit - PSCAD| over the three phases, divided by the peak
  |PSCAD| of the same quantity (terminal voltages v1/v2, source and receiving
  currents), in percent;
- RMS error: RMS difference over RMS of the PSCAD signal, in percent;
- the same metrics between PSCAD at its nominal step and at half the step, where a
  refined run exists, which bounds PSCAD's own discretization error;
- GridKit step statistics: accepted steps, median step over the shortest modal
  delay, and the fraction of elapsed time with steps above every delay.

Outputs: ``plots/<line>/<model>_<scenario>.png``, ``plots/summary.png``,
``validation.json`` and ``validation.md``. Run with the GridWorkbench environment.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "lib"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from scenarios import CODES, MODELS, PSCAD_STEPS, project_name  # noqa: E402

GRIDKIT = HERE / "results" / "gridkit"
PSCAD = HERE / "results" / "pscad"
PLOTS = HERE / "plots"
#: GridKit column prefix -> PSCAD column prefix, and the quantity's label.
QUANTITIES = {"v1": ("v1", "Sending-end voltage"), "v2": ("v2", "Receiving-end voltage"),
              "source": ("is", "Source current"), "i2": ("i2", "Receiving-end current")}
COLORS = ("#0072B2", "#D55E00", "#009E73")
SCENARIO_ORDER = ("energize", "pulse", "fault", "sag")
#: Where each scenario's detail window sits: (start, width) in seconds, None for the full run.
ZOOM = {"energize": (0.9e-3, 1.6e-3), "pulse": (0.9e-3, 1.4e-3), "fault": (19.5e-3, 2.5e-3), "sag": (19.5e-3, 2.5e-3)}


def load_gridkit(line, model, scenario):
    path = GRIDKIT / f"{line}__{model}" / scenario
    if not (path / "statistics.json").exists():
        return None
    samples = np.genfromtxt(path / "samples.csv", delimiter=",", names=True)
    stats = json.loads((path / "statistics.json").read_text())
    data = json.loads((path / "input.json").read_text())
    steps = np.genfromtxt(path / "steps.csv", delimiter=",", names=True)
    return {"samples": samples, "statistics": stats, "input": data, "steps": steps}


def load_pscad(line, model, scenario, refined=False):
    name = project_name(line, model, scenario, refined)
    path = PSCAD / name / f"{name}.csv"
    if not path.exists():
        return None
    return np.genfromtxt(path, delimiter=",", names=True)


def on_grid(pscad, t):
    """PSCAD channels interpolated onto the times ``t``."""
    return {name: np.interp(t, pscad["t"], pscad[name]) for name in pscad.dtype.names if name != "t"}


def errors(reference, actual, names_reference, names_actual, mask):
    """Peak- and RMS-normalized errors of three phases against ``reference``."""
    ref = np.column_stack([reference[n][mask] for n in names_reference])
    act = np.column_stack([actual[n][mask] for n in names_actual])
    difference = act - ref
    peak = np.max(np.abs(ref))
    rms = np.sqrt(np.mean(ref**2))
    index = np.unravel_index(np.argmax(np.abs(difference)), difference.shape)
    return {"peak_percent": float(100 * np.max(np.abs(difference)) / peak) if peak > 0 else 0.0,
            "rms_percent": float(100 * np.sqrt(np.mean(difference**2)) / rms) if rms > 0 else 0.0,
            "reference_peak": float(peak), "max_abs_difference": float(np.max(np.abs(difference))),
            "worst_phase": "abc"[index[1]]}


def step_metrics(data, steps):
    delays = np.array(data["line"]["H"]["delays"])
    delays = delays[delays > 0]
    h, t = steps["h"], steps["t"]
    duration = np.minimum(t, data["solver"]["tf"]) - np.maximum(t - h, data["solver"]["t0"])
    keep = duration > 0
    h, duration = h[keep], duration[keep]
    above = h > delays.max()
    return {"accepted_steps": int(len(steps)), "min_delay_s": float(delays.min()), "max_delay_s": float(delays.max()),
            "median_h_over_min_delay": float(np.median(h / delays.min())),
            "max_h_over_max_delay": float(np.max(h / delays.max())),
            "time_fraction_above_all_delays": float(duration[above].sum() / duration.sum()),
            "max_history_segments": int(np.max(steps["segments"]))}


def compare_one(line, model, scenario):
    gk = load_gridkit(line, model, scenario)
    ps = load_pscad(line, model, scenario)
    if gk is None or ps is None:
        return None
    samples = gk["samples"]
    t = samples["t"]
    tf = gk["input"]["solver"]["tf"]
    grid = on_grid(ps, t)
    window = (t >= 0) & (t <= tf)
    entry = {"line": line, "model": model, "scenario": scenario, "tf_s": tf,
             "pscad_step_s": PSCAD_STEPS[scenario][0],
             "gridkit": {"seconds": gk["statistics"]["seconds_including_csv"], "unknowns": gk["statistics"]["unknowns"],
                         **step_metrics(gk["input"], gk["steps"]), "history_rejections": gk["statistics"]["history_rejections"],
                         "error_test_failures": gk["statistics"]["error_test_failures_including_history"]},
             "forcing_check_percent": errors(grid, samples, ["e_a", "e_b", "e_c"], ["e_a", "e_b", "e_c"], window)["peak_percent"],
             "errors": {}}
    for key, (ps_key, _) in QUANTITIES.items():
        entry["errors"][key] = errors(grid, samples, [f"{ps_key}_{p}" for p in "abc"], [f"{key}_{p}" for p in "abc"], window)
    events = gk["input"]["events"]
    if gk["input"]["initialization"] == "sinusoidal" and events:
        before = (t >= 0) & (t < min(e["on"] for e in events))
        entry["steady_state_errors"] = {
            key: errors(grid, samples, [f"{ps_key}_{p}" for p in "abc"], [f"{key}_{p}" for p in "abc"], before)
            for key, (ps_key, _) in QUANTITIES.items()}
    fine = load_pscad(line, model, scenario, refined=True)
    if fine is not None:
        fine_grid = on_grid(fine, t)
        entry["pscad_refinement"] = {
            key: {"coarse_vs_fine_peak_percent": errors(fine_grid, grid, [f"{ps_key}_{p}" for p in "abc"], [f"{ps_key}_{p}" for p in "abc"], window)["peak_percent"],
                  "gridkit_vs_fine_peak_percent": errors(fine_grid, samples, [f"{ps_key}_{p}" for p in "abc"], [f"{key}_{p}" for p in "abc"], window)["peak_percent"]}
            for key, (ps_key, _) in QUANTITIES.items()}
    figure(line, model, scenario, gk, grid, entry)
    return entry


def figure(line, model, scenario, gk, grid, entry):
    samples = gk["samples"]
    t = samples["t"]
    tf = gk["input"]["solver"]["tf"]
    unit_v, scale_v = ("kV", 1e-3)
    peak_i = entry["errors"]["i2"]["reference_peak"]
    unit_i, scale_i = ("kA", 1e-3) if peak_i >= 2000 else ("A", 1.0)
    fig, axes = plt.subplots(3, 2, figsize=(14, 9.5), layout="constrained")
    start, width = ZOOM[scenario]
    windows = ((0, tf), (start, min(start + width, tf)))
    for column, (lo, hi) in enumerate(windows):
        m = (t >= lo) & (t <= hi)
        for p, color in zip("abc", COLORS):
            axes[0, column].plot(t[m] * 1e3, samples[f"v2_{p}"][m] * scale_v, color=color, lw=1.0, label=f"GridKit {p.upper()}")
            axes[0, column].plot(t[m] * 1e3, grid[f"v2_{p}"][m] * scale_v, "--", color="black", lw=0.7, label="PSCAD" if p == "a" else None)
            axes[1, column].plot(t[m] * 1e3, samples[f"i2_{p}"][m] * scale_i, color=color, lw=1.0, label=f"GridKit {p.upper()}")
            axes[1, column].plot(t[m] * 1e3, grid[f"i2_{p}"][m] * scale_i, "--", color="black", lw=0.7, label="PSCAD" if p == "a" else None)
        for (key, (ps_key, label)), color in zip(QUANTITIES.items(), [COLORS[0], COLORS[1], COLORS[2], "#CC79A7"]):
            ref_peak = entry["errors"][key]["reference_peak"]
            difference = np.max(np.abs(np.column_stack([samples[f"{key}_{p}"][m] - grid[f"{ps_key}_{p}"][m] for p in "abc"])), axis=1)
            axes[2, column].plot(t[m] * 1e3, 100 * difference / ref_peak, color=color, lw=0.8, label=label)
        axes[0, column].set_ylabel(f"v2 ({unit_v})")
        axes[1, column].set_ylabel(f"i2 ({unit_i})")
        axes[2, column].set(ylabel="|GridKit - PSCAD| / peak (%)", xlabel="Study time (ms)", yscale="log")
        axes[2, column].set_ylim(1e-5, 100)
        for ax in axes[:, column]:
            ax.grid(color="0.9", lw=0.5)
            ax.set_xlim(lo * 1e3, hi * 1e3)
    axes[0, 0].set_title("Receiving-end voltage, full run")
    axes[0, 1].set_title("Detail")
    axes[1, 0].set_title("Receiving-end current")
    axes[2, 0].set_title("Normalized difference")
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=4)
    axes[2, 0].legend(frameon=False, fontsize=8, ncol=2)
    e = entry["errors"]
    fig.suptitle(f"{line} / {MODELS[model].split(' (')[0]} / {scenario}  -  peak errors: v1 {e['v1']['peak_percent']:.3f}%, "
                 f"v2 {e['v2']['peak_percent']:.3f}%, source current {e['source']['peak_percent']:.3f}%, i2 {e['i2']['peak_percent']:.3f}%"
                 f"  |  GridKit {entry['gridkit']['accepted_steps']} steps, median h/tau(min) {entry['gridkit']['median_h_over_min_delay']:.2f}; "
                 f"PSCAD {int(round(tf / PSCAD_STEPS[scenario][0]))} steps")
    destination = PLOTS / line
    destination.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination / f"{model}_{scenario}.png", dpi=130)
    plt.close(fig)


def summary_figure(entries):
    lines = list(CODES)
    rows = [(line, model) for line in lines for model in REPORTED_MODELS]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), layout="constrained")
    for ax, (key, title) in zip(axes, (("v2", "Receiving-end voltage"), ("i2", "Receiving-end current"))):
        matrix = np.full((len(rows), len(SCENARIO_ORDER)), np.nan)
        for entry in entries:
            matrix[rows.index((entry["line"], entry["model"])), SCENARIO_ORDER.index(entry["scenario"])] = entry["errors"][key]["peak_percent"]
        image = ax.imshow(np.log10(np.maximum(matrix, 1e-4)), cmap="viridis_r", vmin=-3, vmax=1, aspect="auto")
        for i in range(len(rows)):
            for j in range(len(SCENARIO_ORDER)):
                if not np.isnan(matrix[i, j]):
                    ax.text(j, i, f"{matrix[i, j]:.2g}", ha="center", va="center", fontsize=8,
                            color="white" if matrix[i, j] > 0.3 else "black")
        ax.set_xticks(range(len(SCENARIO_ORDER)), SCENARIO_ORDER)
        ax.set_yticks(range(len(rows)), [f"{line} / {model}" for line, model in rows])
        ax.set_title(f"{title}: peak |GridKit - PSCAD| / peak |PSCAD| (%)")
    fig.colorbar(image, ax=axes, label="log10 of peak error (%)", shrink=0.8)
    PLOTS.mkdir(exist_ok=True)
    fig.savefig(PLOTS / "summary.png", dpi=140)
    plt.close(fig)


def markdown(entries):
    out = ["# GridKit LineDistributed against PSCAD: validation summary", "",
           "Peak errors are max |GridKit - PSCAD| over the three phases and the study window, normalized by the peak "
           "|PSCAD| of the same quantity. RMS errors are normalized by the PSCAD RMS. GridKit steps are the accepted "
           "adaptive steps; h/tau is the accepted step over the shortest modal delay (median) and 'time above all "
           "delays' is the fraction of elapsed time with steps longer than every modal delay. PSCAD steps are tf over the fixed step.", "",
           "| Line | Model | Scenario | v1 peak % | v2 peak % | source i peak % | i2 peak % | v2 RMS % | i2 RMS % | PSCAD h/2 effect on v2 % | GridKit steps | median h/tau | time above all delays | GridKit s | PSCAD steps |",
           "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for entry in entries:
        e, g = entry["errors"], entry["gridkit"]
        refinement = entry.get("pscad_refinement", {}).get("v2", {}).get("coarse_vs_fine_peak_percent")
        out.append("| {line} | {model} | {scenario} | {v1:.3f} | {v2:.3f} | {src:.3f} | {i2:.3f} | {v2r:.3f} | {i2r:.3f} | {ref} | {steps} | {h:.2f} | {above:.0%} | {sec:.1f} | {psteps} |".format(
            line=entry["line"], model=entry["model"], scenario=entry["scenario"],
            v1=e["v1"]["peak_percent"], v2=e["v2"]["peak_percent"], src=e["source"]["peak_percent"], i2=e["i2"]["peak_percent"],
            v2r=e["v2"]["rms_percent"], i2r=e["i2"]["rms_percent"],
            ref=f"{refinement:.3f}" if refinement is not None else "-",
            steps=g["accepted_steps"], h=g["median_h_over_min_delay"], above=g["time_fraction_above_all_delays"],
            sec=g["seconds"], psteps=int(round(entry["tf_s"] / entry["pscad_step_s"]))))
    steady = [entry for entry in entries if "steady_state_errors" in entry]
    if steady:
        out += ["", "## Pre-event steady state (sinusoidal cases, 0 to 20 ms)", "",
                "| Line | Model | Scenario | v2 peak % | i2 peak % | source i peak % |", "|---|---|---|---:|---:|---:|"]
        for entry in steady:
            s = entry["steady_state_errors"]
            out.append(f"| {entry['line']} | {entry['model']} | {entry['scenario']} | {s['v2']['peak_percent']:.3f} | {s['i2']['peak_percent']:.3f} | {s['source']['peak_percent']:.3f} |")
    return "\n".join(out) + "\n"


#: The line models the report covers: the frequency-dependent model is the general case.
REPORTED_MODELS = ("fd",)


def main(argv):
    entries = []
    for line in CODES:
        for model in REPORTED_MODELS:
            for scenario in SCENARIO_ORDER:
                if argv and not any(a in f"{line} {model} {scenario}" for a in argv):
                    continue
                entry = compare_one(line, model, scenario)
                if entry is None:
                    print(f"missing results for {line} {model} {scenario}", flush=True)
                    continue
                e = entry["errors"]
                print(f"{line:18s} {model} {scenario:9s} v1 {e['v1']['peak_percent']:7.3f}%  v2 {e['v2']['peak_percent']:7.3f}%  "
                      f"is {e['source']['peak_percent']:7.3f}%  i2 {e['i2']['peak_percent']:7.3f}%  forcing {entry['forcing_check_percent']:.2e}%  "
                      f"steps {entry['gridkit']['accepted_steps']:7d}  h/tau {entry['gridkit']['median_h_over_min_delay']:.2f}", flush=True)
                entries.append(entry)
    if entries:
        summary_figure(entries)
        (HERE / "validation.json").write_text(json.dumps(entries, indent=2) + "\n")
        (HERE / "validation.md").write_text(markdown(entries))
        print("Wrote", HERE / "validation.md")


if __name__ == "__main__":
    main(sys.argv[1:])
