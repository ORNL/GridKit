"""Physical line figures; always written to this source example's plots/."""

import json

import numpy as np

from line_constants import HERE, T, constants, modal
from plot import COLORS, plt, save
from validate import fit, load, operators
from validate_physical import RESULTS, length_of, step_metrics

DESTINATION = HERE / "plots"


def finish(fig, name):
    save(fig, name, DESTINATION)


def fit_figures():
    f = np.geomspace(.01, 1e6, 1600)
    s = 2j*np.pi*f
    data, _ = load(RESULTS / "physical_pulse")
    yraw, hraw = modal(s, length_of(data))
    yc, h = operators(data, s)
    yc = np.diagonal(T.T@yc@T, axis1=1, axis2=2)
    h = np.diagonal(T.T@h@T, axis1=1, axis2=2)
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), layout="constrained")
    for m, color, name in zip((0, 1), COLORS, ("Common mode", "Differential mode")):
        for ax, original, fitted, scale in ((axes[0, 0], yraw[:, m], yc[:, m], 1000),
                                            (axes[0, 1], hraw[:, m], h[:, m], 1)):
            ax.loglog(f, abs(original)*scale, color=color, label=name+" physical")
            ax.loglog(f[::20], abs(fitted[::20])*scale, ".", color=color, ms=3, label=name+" fitted")
        axes[1, 0].semilogx(f, 100*abs(yc[:, m]-yraw[:, m])/abs(yraw[:, m]), color=color, label=name)
        axes[1, 1].semilogx(f, 100*abs(h[:, m]-hraw[:, m])/np.maximum(abs(hraw[:, m]), 1e-3), color=color, label=name)
    axes[0, 0].set(title="Characteristic admittance", ylabel=r"$|Y_c|$ (mS)")
    axes[0, 1].set(title="Propagation through 80 km", ylabel=r"$|H|$", ylim=(1e-4, 1.1))
    axes[1, 0].set(title="Held-out complex admittance error", ylabel="Relative error (%)")
    axes[1, 1].set(title="Held-out complex propagation error", ylabel="Scaled error (%) — magnitude floor 0.001")
    for ax in axes.flat:
        ax.set_xlabel("Frequency (Hz)")
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Physical overhead line: raw line constants and actual vector fits")
    finish(fig, "physical_fits")

    z, y = constants(s)
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), layout="constrained")
    for m, color, name in zip((0, 1), COLORS, ("Common mode", "Differential mode")):
        axes[0, 0].loglog(f, 1000*z[:, m].real, color=color, label=name)
        axes[0, 1].semilogx(f, 1e6*z[:, m].imag/(2*np.pi*f), color=color, label=name)
        axes[1, 0].semilogx(f, 1e12*y[:, m].imag/(2*np.pi*f), color=color, label=name)
        # Unwrap raw phase before differentiation; this is band-dependent group
        # delay, distinct from the explicit high-frequency wavefront delay.
        delay = data["line"]["H"]["delays"][m] - np.gradient(np.unwrap(np.angle(hraw[:, m]*np.exp(s*data["line"]["H"]["delays"][m]))), 2*np.pi*f)
        axes[1, 1].semilogx(f, delay*1e6, color=color, label=name+" group delay")
        axes[1, 1].axhline(data["line"]["H"]["delays"][m]*1e6, color=color, ls=":", lw=1)
    for ax, title, units in zip(axes.flat,
            ("Resistance including skin and ground return", "Inductance including ground return", "Capacitance", "80 km delay: dispersion and wavefront"),
            ("Resistance (ohm/km)", "Inductance (mH/km)", "Capacitance (nF/km)", "Delay (microseconds)")):
        ax.set(title=title, ylabel=units, xlabel="Frequency (Hz)")
        ax.legend(frameon=False, fontsize=8)
    axes[1, 1].set_ylim(200, 650)
    fig.suptitle("Documented 22 kV design: solid aluminium, ideal transposition, 100 ohm m soil")
    finish(fig, "physical_line_constants")


def pulse_figures():
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), layout="constrained")
    for column, name in enumerate(("physical_pulse", "physical_reflections")):
        data, samples = load(RESULTS / name)
        reference = np.genfromtxt(RESULTS / name / "reference.csv", delimiter=",", names=True)
        t = samples["t"]*1e3
        ax = axes[0, column]
        ax.plot(t, samples["v1_a"]/1000, color=COLORS[0], label="Sending A")
        ax.plot(t, samples["v2_a"]/1000, color=COLORS[1], label="Receiving A")
        ax.plot(reference["t"]*1e3, reference["v2_a"]/1000, "k--", lw=.8, label="Unfitted telegrapher reference")
        ax.set_title(f"Load {data['load']['resistance'][0]:g} ohm/phase")
        ax.set_ylabel("Voltage (kV)")
        ax = axes[1, column]
        for phase, color in zip("abc", COLORS):
            ax.plot(t, samples[f"v2_{phase}"]/1000, color=color, label="Receiving "+phase.upper())
        ax.set_ylabel("Receiving voltage (kV)")
        ax.set_title("Cross-phase response (B = C by symmetry)")
        for ax in axes[:, column]:
            for k in (1, 3, 5):
                ax.axvline((.001+k*min(data["line"]["H"]["delays"]))*1000, color=".6", ls=":", lw=.7)
            ax.set(xlim=(.95, 3), xlabel="Time (ms)")
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle("80 km line: phase-A pulse, 5 us rise, 50 us FWHM — finite source impedance")
    finish(fig, "physical_travelling_waves")

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), layout="constrained", sharex=True)
    for name, m, color, label in (("physical_common_mode", 0, COLORS[0], "Common-mode pulse"),
                                  ("physical_differential_mode", 1, COLORS[1], "Differential-mode pulse")):
        data, samples = load(RESULTS / name)
        reference = np.genfromtxt(RESULTS / name / "reference.csv", delimiter=",", names=True)
        v1 = np.column_stack([samples[f"v1_{p}"] for p in "abc"])@T
        v2 = np.column_stack([samples[f"v2_{p}"] for p in "abc"])@T
        ref = np.column_stack([reference[f"v2_{p}"] for p in "abc"])@T
        scale = np.max(abs(v1[:, m]))
        axes[0].plot(samples["t"]*1000, v1[:, m]/scale, color=color, label=label)
        axes[1].plot(samples["t"]*1000, v2[:, m]/scale, color=color, label=label)
        axes[1].plot(reference["t"]*1000, ref[:, m]/scale, "k--", lw=.7,
                     label="Unfitted reference" if m == 0 else None)
    for ax in axes:
        ax.set(xlim=(.95, 2.2), ylabel="Modal voltage / launched peak", xlabel="Time (ms)")
        ax.legend(frameon=False)
    axes[0].set_title("Launched modal voltage")
    axes[1].set_title("Received response — ground-return mode broadens and attenuates")
    fig.suptitle("Same physical line and termination; different modal excitations")
    finish(fig, "physical_modes")


def circuit_figures():
    for name in ("physical_fault", "physical_short_recovery"):
        data, samples = load(RESULTS / name)
        steps = np.genfromtxt(RESULTS / name / "steps.csv", delimiter=",", names=True)
        stats = step_metrics(RESULTS / name, .1, .2)
        fig, axes = plt.subplots(3, 1, figsize=(11, 8), layout="constrained", sharex=True)
        for p, color in zip("abc", COLORS):
            axes[0].plot(samples["t"]*1000, samples[f"v2_{p}"]/1000, color=color, lw=.8, label=p.upper())
            axes[1].plot(samples["t"]*1000, samples[f"source_{p}"], color=color, lw=.8, label=p.upper())
        maximum = max(data["line"]["H"]["delays"])
        minimum = min(data["line"]["H"]["delays"])
        axes[2].semilogy(steps["t"]*1000, steps["h"]/maximum, color=COLORS[0], lw=.8, label=r"$h/\tau_{max}$")
        axes[2].semilogy(steps["t"]*1000, steps["h"]/minimum, color=COLORS[1], lw=.6, alpha=.6, label=r"$h/\tau_{min}$")
        axes[2].axhline(1, color=".2", ls="--", lw=1)
        for ax in axes:
            for event in data["events"]:
                ax.axvspan(event["on"]*1000, event["off"]*1000, color=".6", alpha=.12)
            ax.axvspan(100, 200, color=COLORS[2], alpha=.07)
            ax.set_xlim(0, 200)
            ax.legend(frameon=False, ncol=3)
        axes[0].set_ylabel("Receiving voltage (kV)")
        axes[1].set_ylabel("Source current (A)")
        axes[2].set(ylabel="Accepted step / delay", xlabel="Time (ms)")
        axes[2].set_title(f"100–200 ms: {stats['time_fraction_above_all_delays']:.1%} of elapsed time above every delay; "
                          f"median h/tau(max) = {stats['median_h_over_tau_max']:.2f}")
        title = "80 km: phase-A 10 ohm terminal fault and clearing" if name == "physical_fault" else "1 km: phase-A voltage sag and recovered energized operation"
        load_label = "400 ohm + 50 mH" if name == "physical_fault" else "400 ohm"
        fig.suptitle(title+f"\n22 kV source, 250 MVA short-circuit strength, X/R=10; {load_label} load/phase")
        finish(fig, name)

    path = RESULTS / "physical_short_recovery"
    data, samples = load(path)
    steps = np.genfromtxt(path / "steps.csv", delimiter=",", names=True)
    reference = np.genfromtxt(path / "reference.csv", delimiter=",", names=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 6.5), layout="constrained")
    for col, window in enumerate(((19.98, 20.25), (100, 150))):
        axes[0, col].plot(samples["t"]*1000, samples["v2_a"]/1000, color=COLORS[0], label="Simulated receiving A")
        axes[0, col].plot(reference["t"]*1000, reference["v2_a"]/1000, "k--", lw=.7, label="Unfitted reference")
        axes[1, col].semilogy(steps["t"]*1000, steps["h"]/max(data["line"]["H"]["delays"]), color=COLORS[0])
        axes[1, col].axhline(1, color=".3", ls="--")
        for ax in axes[:, col]:
            ax.set(xlim=window, xlabel="Time (ms)")
        axes[0, col].legend(frameon=False, fontsize=8)
        visible = (samples["t"]*1000 >= window[0]) & (samples["t"]*1000 <= window[1])
        low, high = np.min(samples["v2_a"][visible])/1000, np.max(samples["v2_a"][visible])/1000
        margin = .1*max(high-low, 1)
        axes[0, col].set_ylim(low-margin, high+margin)
        axes[1, col].set_ylim(.005, 100)
        axes[0, col].set_ylabel("Receiving voltage (kV)")
        axes[1, col].set_ylabel(r"$h/\tau_{max}$")
    metrics = step_metrics(path, .1, .2)
    axes[1, 1].text(.03, .78, f"100–200 ms: {metrics['time_fraction_above_all_delays']:.0%} of elapsed time above all delays\n"
                    f"Median h/tau(max) = {metrics['median_h_over_tau_max']:.2f}",
                    transform=axes[1, 1].transAxes, va="top", fontsize=10)
    axes[0, 0].set_title("Resolved disturbance")
    axes[0, 1].set_title("Three recovered 60 Hz cycles")
    fig.suptitle("Adaptive delay integration: transient resolution and steps exceeding all modal delays")
    finish(fig, "physical_adaptivity")


def main():
    if "build" in HERE.parts:
        raise RuntimeError("Run from the source example; plots must not go in build/")
    DESTINATION.mkdir(exist_ok=True)
    fit_figures()
    pulse_figures()
    circuit_figures()
    print("Saved physical PNG/SVG figures in", DESTINATION)


if __name__ == "__main__":
    main()
