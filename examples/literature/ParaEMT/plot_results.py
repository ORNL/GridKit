"""Plot ParaEMT-only results and report differences when halving its time step."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
COLORS = ("#2874a6", "#d35400", "#229954")


def compare(directory):
    coarse, fine = [pd.read_csv(directory / step / "reference.csv.gz")
                    for step in ("dt50us", "dt25us")]
    assert list(coarse.columns) == list(fine.columns)
    assert np.array_equal(coarse.time_s, fine.time_s)
    assert np.isfinite(coarse).all().all() and np.isfinite(fine).all().all()
    assert len(coarse) == 6001 and coarse.time_s.iloc[-1] == 3.0
    metrics = {}
    for group, pattern in {
            "phase_voltage_pu": r"_v[abc]_pu$", "voltage_magnitude_pu": r"_vm_pu$",
            "rotor_speed_pu": "_speed_pu$", "exciter_voltage_pu": "_efd_pu$",
            "mechanical_power_pu": "_pm_pu$", "stabilizer_output_pu": "_pss_vs_pu$"}.items():
        columns = coarse.filter(regex=pattern).columns
        difference = coarse[columns] - fine[columns]
        metrics[group] = {}
        for window, mask in {
                "all_samples": coarse.time_s >= 0,
                "before_event": coarse.time_s < 1,
                "event_through_20ms": (coarse.time_s >= 1) & (coarse.time_s <= 1.02),
                "after_20ms": coarse.time_s > 1.02}.items():
            values = difference.loc[mask].to_numpy()
            row, col = np.unravel_index(np.abs(values).argmax(), values.shape)
            metrics[group][window] = {
                "max_absolute_difference": float(np.abs(values).max()),
                "rms_difference": float(np.sqrt(np.mean(values**2))),
                "max_at_time_s": float(coarse.loc[mask, "time_s"].iloc[row]),
                "max_at_column": columns[col],
            }
    metrics["sampled_bus1_peak_magnitude_pu"] = {
        "dt50us": float(coarse.bus1_vm_pu.max()), "dt25us": float(fine.bus1_vm_pu.max())}
    metrics["interpretation"] = "ParaEMT time-step comparison only; no GridKit or independently verified truth data"
    return coarse, fine, metrics


def style(axes):
    for ax in np.asarray(axes).flat:
        ax.axvline(1, color="0.5", linestyle=":", linewidth=1)
        ax.set_xlabel("Time (s)")
        ax.grid(alpha=0.22)
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)


def main():
    plots = ROOT / "plots"
    plots.mkdir(exist_ok=True)
    all_metrics = {}
    for event in ("governor_step", "trip"):
        coarse, fine, metrics = compare(ROOT / "results" / event)
        all_metrics[event] = metrics
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
        if event == "governor_step":
            title = "ParaEMT 9-bus: generator 1 governor reference −0.02 pu at 1 s"
            generators = (1, 2, 3)
        else:
            title = "ParaEMT 9-bus: generator 1 trip at 1 s — exploratory reference"
            # The upstream kernel continues updating disconnected machine states.
            generators = (2, 3)
        for ax, variable, ylabel in zip(axes.flat,
                ("speed", "vm", "pm", "efd"),
                ("Rotor speed (pu)", "Bus magnitude (pu)",
                 "Mechanical power (pu, machine base)", "Field voltage (pu, exciter base)")):
            for i, gen in enumerate(generators):
                key = f"bus{gen}_vm_pu" if variable == "vm" else f"gen{gen}_{variable}_pu"
                label = f"Bus {gen}" if variable == "vm" else f"Gen {gen}"
                ax.plot(fine.time_s, fine[key], color=COLORS[i], label=label, linewidth=1.5)
                ax.plot(coarse.time_s, coarse[key], color=COLORS[i], linestyle="--", linewidth=1)
            ax.set_ylabel(ylabel)
            ax.legend(fontsize=9)
        style(axes)
        fig.suptitle(title + "\nSolid: 25 µs; dashed: 50 µs. ParaEMT only; GridKit equivalence blocked.", fontsize=12)
        fig.savefig(plots / f"{event}_response.png", dpi=170)
        plt.close(fig)

    coarse, fine, _ = compare(ROOT / "results" / "trip")
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), layout="constrained")
    mask = (fine.time_s >= 0.99) & (fine.time_s <= 1.03)
    for frame, label, linestyle in ((coarse, "50 µs", "--"), (fine, "25 µs", "-")):
        axes[0].plot(frame.time_s[mask], frame.bus1_vm_pu[mask], linestyle, marker=".", label=label)
    axes[0].set_ylabel("Bus 1 magnitude (pu)")
    axes[0].legend()
    for i, phase in enumerate("abc"):
        key = f"bus4_v{phase}_pu"
        axes[1].plot(fine.time_s[mask], fine[key][mask], color=COLORS[i], label=f"Phase {phase}")
        axes[1].plot(coarse.time_s[mask], coarse[key][mask], color=COLORS[i], linestyle="--")
    axes[1].set_ylabel("Bus 4 phase voltage (pu)")
    axes[1].legend()
    style(axes)
    fig.suptitle("ParaEMT trip: the sampled voltage spike grows when the step is halved\n"
                 "Output spacing is 0.5 ms; sub-sample peaks are not resolved. GridKit not run.", fontsize=12)
    fig.savefig(plots / "trip_waveforms.png", dpi=170)
    plt.close(fig)
    (ROOT / "results/refinement.json").write_text(json.dumps(all_metrics, indent=2) + "\n")
    print(json.dumps({event: {group: data["all_samples"]["max_absolute_difference"]
                            for group, data in metrics.items() if "all_samples" in data}
                      for event, metrics in all_metrics.items()}, indent=2))


if __name__ == "__main__":
    main()
