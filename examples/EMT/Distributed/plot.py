#!/usr/bin/env python3
"""Plot retained EMT studies and collect solver, delay, waveform, and fit statistics."""

import argparse
import csv as csv_module
import html
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fit import response

COLORS = {
    "distributed": "#176b91",
    "distributed_capped": "#479b70",
    "pi_1": "#d06030",
    "pi_10": "#9471ba",
}
LABELS = {
    "distributed": "Distributed",
    "distributed_capped": "Distributed, delay cap",
    "pi_1": "One π section",
    "pi_10": "Ten π sections",
}


def csv(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def step_stats(steps, minimum, maximum, start=0):
    mask = steps["time"] > start
    h = steps["step"][mask]
    advanced = np.minimum(h, steps["time"][mask] - start)
    return {
        "steps_observed": len(h),
        "minimum_s": float(np.min(h)),
        "mean_s": float(np.mean(h)),
        "time_percent_above_shortest_delay": float(
            100 * np.sum(advanced[h > minimum * (1 + 1e-8)]) / np.sum(advanced)
        ),
        "time_percent_above_longest_delay": float(
            100 * np.sum(advanced[h > maximum * (1 + 1e-8)]) / np.sum(advanced)
        ),
        "median_s": float(np.median(h)),
        "p95_s": float(np.quantile(h, 0.95)),
        "maximum_s": float(np.max(h)),
        "maximum_over_shortest_delay": float(np.max(h) / minimum),
        "percent_above_shortest_delay": float(100 * np.mean(h > minimum * (1 + 1e-8))),
        "percent_above_longest_delay": float(100 * np.mean(h > maximum * (1 + 1e-8))),
    }


def collect(folder):
    meta = json.loads((folder / "model.json").read_text())
    records = json.loads((folder / "statistics.json").read_text())
    successful = [r for r in records if r["returncode"] == 0]
    if not successful:
        return None
    trial = successful[0]["trial"]
    solver = json.loads((folder / f"trial{trial}/solver.json").read_text())
    meta.setdefault("events", solver["events"])
    meta.setdefault("horizon_s", solver["tmax"])
    values = csv(folder / f"trial{trial}/response.csv")
    steps = csv(folder / f"trial{trial}/steps.csv")
    summary = {
        "name": folder.name,
        "source_folder": str(folder.resolve()),
        **meta,
        "successful_trials": len(successful),
        "solver": {
            key: float(np.median([r[key] for r in successful]))
            for key in successful[0]
            if key not in ["returncode", "trial", "executable_sha256"]
        },
        "cpu_range_s": [
            min(r["cpu_s"] for r in successful),
            max(r["cpu_s"] for r in successful),
        ],
        "executable_sha256": successful[0]["executable_sha256"],
    }
    tau = (
        [meta["tau"]]
        if meta["family"] in ["bergeron", "frequency"]
        else [t for e in meta["edges"] for t in e["delays_s"]]
    )
    summary["delay"] = {
        "min_s": min(tau),
        "median_s": float(np.median(tau)),
        "max_s": max(tau),
    }
    summary["solver"]["atol"] = solver["abs_tol"]
    summary["output_samples"] = len(values)
    summary["accepted_step_count_matches_solver"] = len(steps) == successful[0]["steps"]
    summary["bdf_order_counts"] = {
        str(int(k)): int(v)
        for k, v in zip(*np.unique(steps["order"], return_counts=True))
    }
    summary["internal_steps"] = step_stats(steps, min(tau), max(tau))
    summary["late_internal_steps"] = step_stats(
        steps, min(tau), max(tau), values["t"][-1] * 0.8
    )
    voltage = np.column_stack(
        [
            values[key]
            for key in values.dtype.names
            if key.startswith("Bus_b") and key.endswith(("_va", "_vb", "_vc"))
        ]
    )
    summary["maximum_bus_phase_voltage_V"] = float(np.max(abs(voltage)))
    for key in values.dtype.names:
        if key.endswith("_omega"):
            summary.setdefault("machine_speed", {})[key] = {
                "min_pu": float(np.min(values[key])),
                "max_pu": float(np.max(values[key])),
            }
    return summary, values, steps


def bergeron(records, output):
    for sine, short in [(False, False), (True, True)]:
        cases = [
            r
            for r in records
            if r[0]["family"] == "bergeron"
            and r[0]["sine"] == sine
            and r[0]["short"] == short
        ]
        if not cases:
            continue
        fig, axes = plt.subplots(3, 1, figsize=(11, 9))
        for meta, d, steps in cases:
            kind = meta["kind"]
            color = COLORS[kind]
            label = LABELS[kind]
            tau = meta["tau"]
            t = d["t"]
            voltage = d["Bus_b1_va"]
            if sine:
                w = 2 * np.pi * 60
                phasor = 1000 / (
                    (1 + 100 / 600) * np.cos(w * tau)
                    + 1j * (300 / 600 + 100 / 300) * np.sin(w * tau)
                )
                exact = np.real(phasor * np.exp(1j * w * t))
                mask = t > 0.06
                error = voltage - exact
                meta["analytic_steady_error"] = {
                    "rms_V": float(np.sqrt(np.mean(error[mask] ** 2))),
                    "maximum_V": float(np.max(abs(error[mask]))),
                    "window_s": [0.06, float(t[-1])],
                    "receiving_peak_V": float(abs(phasor)),
                }
                axes[0].plot(t * 1e3, voltage, color=color, label=label, lw=1)
                axes[1].plot(
                    t[mask] * 1e3, error[mask] * 1e3, color=color, label=label, lw=0.9
                )
            else:
                exact = sum(
                    1000 * (-1 / 6) ** n * (t > 0.001 + (2 * n + 1) * tau)
                    for n in range(15)
                )
                distance = abs(
                    (t - 0.001 - tau) / (2 * tau)
                    - np.round((t - 0.001 - tau) / (2 * tau))
                )
                mask = distance > 1e-5
                meta["analytic_step_error_away_from_fronts"] = {
                    "rms_V": float(
                        np.sqrt(np.mean((voltage[mask] - exact[mask]) ** 2))
                    ),
                    "maximum_V": float(np.max(abs(voltage[mask] - exact[mask]))),
                }
                axes[0].plot(t * 1e3, voltage, color=color, label=label, lw=1)
                axes[1].plot(t * 1e3, d["Bus_b0_va"], color=color, label=label, lw=1)
            axes[2].semilogy(
                steps["time"] * 1e3,
                steps["step"] / tau,
                color=color,
                label=label,
                lw=0.7,
                alpha=0.85,
            )
        if sine:
            axes[0].plot(
                t * 1e3, exact, "k--", lw=0.8, label="Analytic steady sinusoid"
            )
            axes[0].set_xlim(60, 100)
            axes[1].set_ylabel("Steady voltage error [mV]")
            title = "Lossless 20 µs line: 60 Hz source energization"
        else:
            axes[0].step(
                t * 1e3,
                exact,
                "k--",
                where="post",
                lw=0.8,
                label="Analytic reflections",
            )
            axes[0].set_xlim(0.8, 3.4)
            axes[1].set_xlim(0.8, 3.4)
            axes[1].set_ylabel("Sending voltage [V]")
            title = "Lossless 300 µs line: source step and repeated reflections"
        axes[0].set_ylabel("Receiving voltage [V]")
        axes[0].set_title(title)
        axes[0].legend(ncol=2)
        axes[2].axhline(1, color="k", ls="--", lw=0.8)
        axes[2].set_ylabel("Accepted step / delay")
        for ax in axes:
            ax.set_xlabel("Time [ms]")
            ax.grid(alpha=0.2)
        save(fig, output / f'bergeron_{"sine" if sine else "step"}.png')


def network(records, output):
    groups = {}
    for r in records:
        if r[0]["family"] == "network":
            key = (
                r[0]["buses"],
                r[0]["hybrid"],
                r[0].get("event", r[0]["name"].split("_")[2]),
            )
            groups.setdefault(key, []).append(r)
    for (buses, hybrid, event), cases in groups.items():
        fig, axes = plt.subplots(4, 1, figsize=(12, 11))
        target = buses - 3
        phase = "a" if event == "fault_a" else "b"
        reference = None
        for meta, d, steps in cases:
            kind = meta["kind"]
            color = COLORS[kind]
            label = LABELS[kind]
            t = d["t"]
            v = d[f"Bus_b{target}_v{phase}"]
            event_time = meta["events"][0]["time"] if meta["events"] else 0.0
            mask = (t >= event_time - 0.0003) & (t <= event_time + 0.002)
            rms = np.sqrt(sum(d[f"Bus_b{target}_v{phase}"] ** 2 for phase in "abc") / 3)
            axes[0].plot(t * 1e3, rms / 1e3, color=color, label=label, lw=0.7)
            axes[1].plot(t[mask] * 1e3, v[mask] / 1e3, color=color, label=label, lw=0.9)
            for k, g in enumerate(meta["generators"]):
                axes[2].plot(
                    t * 1e3,
                    d[f"Machine_g{g}_omega"],
                    color=color,
                    ls=["-", "--", ":"][k],
                    label=label + f", g{g}",
                    lw=0.9,
                )
            axes[3].semilogy(
                steps["time"] * 1e3,
                steps["step"] * 1e6,
                color=color,
                label=label,
                lw=0.6,
                alpha=0.8,
            )
            if kind == "distributed":
                reference = (t, v)
        if reference:
            rt, rv = reference
            for meta, d, _ in cases:
                if meta["kind"] == "pi_1":
                    difference = d[f"Bus_b{target}_v{phase}"] - np.interp(
                        d["t"], rt, rv
                    )
                    mask = d["t"] >= max(0, event_time - 0.01)
                    meta["difference_from_distributed"] = {
                        "voltage_phase": phase,
                        "event_window_voltage_rms_V": float(
                            np.sqrt(np.mean(difference[mask] ** 2))
                        ),
                        "event_window_voltage_max_V": float(
                            np.max(abs(difference[mask]))
                        ),
                        "window_s": [max(0, event_time - 0.01), float(d["t"][-1])],
                    }
        for tau in sorted(set(t for e in cases[0][0]["edges"] for t in e["delays_s"])):
            axes[3].axhline(tau * 1e6, color="gray", ls=":", lw=0.7)
        axes[0].set_title(
            f"{buses}-bus network: {event} switching"
            + (
                " with open-loop PWM converters"
                if hybrid
                else " with governed machines"
            )
        )
        axes[0].set_ylabel(f"Bus {target} phase norm [kV]")
        axes[1].set_ylabel(f"Phase {phase} event detail [kV]")
        axes[2].set_ylabel("Machine speed [pu]")
        axes[3].set_ylabel("Accepted step [µs]")
        axes[0].legend()
        axes[2].legend(ncol=2, fontsize=8)
        for ax in axes:
            ax.set_xlabel("Time [ms]")
            ax.grid(alpha=0.2)
        save(
            fig,
            output / f'network{buses}_{"hybrid" if hybrid else "machine"}_{event}.png',
        )
        if hybrid:
            fig, axes = plt.subplots(2, 1, figsize=(11, 6))
            for meta, d, _ in cases:
                color = COLORS[meta["kind"]]
                label = LABELS[meta["kind"]]
                g = meta["converters"][0]
                t = d["t"]
                mask = (t > meta["horizon_s"] * 0.4) & (t < meta["horizon_s"] * 0.45)
                axes[0].plot(
                    t[mask] * 1e3,
                    d[f"Converter_converter{g}_ea"][mask] / 1e3,
                    color=color,
                    label=label,
                    lw=0.8,
                )
                axes[1].plot(
                    t[mask] * 1e3,
                    d[f"DependentVoltageSource_filter{g}_ia"][mask],
                    color=color,
                    label=label,
                    lw=0.8,
                )
            axes[0].set_title(
                f"{buses}-bus open-loop converter, 900 Hz PWM ({event} case)"
            )
            for ax, label in zip(
                axes, ["Converter phase a [kV]", "Filter phase a [A]"]
            ):
                ax.set_ylabel(label)
                ax.set_xlabel("Time [ms]")
                ax.grid(alpha=0.2)
            axes[0].legend()
            save(fig, output / f"network{buses}_converter_{event}.png")


def fits(folder, output):
    data = np.load(folder / "parameters.npz")
    w = data["omega"]
    yc = json.loads((folder / "yc.json").read_text())
    pred = response(yc, w)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes[0, 0].loglog(
        w / (2 * np.pi), abs(data["yc"][:, 0, 0]), "k", label="GridWorkbench"
    )
    axes[0, 0].loglog(w / (2 * np.pi), abs(pred[:, 0, 0]), "--", label="Rational fit")
    axes[1, 0].loglog(
        w / (2 * np.pi),
        np.linalg.norm(pred - data["yc"], axis=(1, 2))
        / np.linalg.norm(data["yc"], axis=(1, 2)),
    )
    axes[0, 0].set_ylabel("|Yc(1,1)| [S]")
    axes[1, 0].set_ylabel("Yc relative matrix error")
    p0 = np.ones((3, 3)) / 3
    ps = [np.eye(3) - p0, p0]
    for length in [5000, 20000, 40000, 60000]:
        model = json.loads((folder / f"h{length}.json").read_text())
        p = sum(
            response(m["H"], w) * np.exp(-1j * w * m["tau"])[:, None, None]
            for m in model["modes"]
        )
        v = sum(
            np.exp(-length * data["lam"][:, m])[:, None, None] * ps[m] for m in range(2)
        )
        axes[0, 1].loglog(w / (2 * np.pi), abs(p[:, 0, 0]), label=f"{length/1000:g} km")
        axes[1, 1].loglog(
            w / (2 * np.pi),
            np.linalg.norm(p - v, axis=(1, 2)) / np.linalg.norm(v, axis=(1, 2)),
            label=f"{length/1000:g} km",
        )
    axes[0, 1].set_ylabel("|H(1,1)|")
    axes[1, 1].set_ylabel("H relative matrix error")
    for ax in axes.flat:
        ax.set_xlabel("Frequency [Hz]")
        ax.grid(alpha=0.2)
    axes[0, 0].legend()
    axes[0, 1].legend()
    fig.suptitle("Transposed GridWorkbench overhead geometry: frequency-dependent fits")
    save(fig, output / "line_fits.png")


def topology(records, output):
    selected = {}
    for meta, _, _ in records:
        if meta["family"] == "network":
            selected[meta["buses"]] = meta
    for count, meta in selected.items():
        main = count - 2
        angle = 2 * np.pi * np.arange(main) / main
        position = np.column_stack((np.cos(angle), np.sin(angle)))
        position = np.vstack((position, [[1.75, -0.4], [1.75, 0.4]]))
        colors = {5000: "#478fa8", 20000: "#64a56a", 40000: "#cf9241", 60000: "#b96783"}
        fig, ax = plt.subplots(figsize=(10, 7))
        for edge in meta["edges"]:
            a, b = edge["from"], edge["to"]
            ax.plot(*position[[a, b]].T, color=colors[edge["length_m"]], lw=2, zorder=1)
        for b in [main, main + 1]:
            ax.plot(*position[[main - 1, b]].T, "k--", lw=1, zorder=1)
        for b in range(count):
            marker = (
                "s"
                if b in meta["generators"]
                else ("^" if b in meta["converters"] else "o")
            )
            ax.scatter(
                *position[b],
                s=500,
                marker=marker,
                facecolor="white",
                edgecolor="#253b50",
                zorder=2,
            )
            ax.text(*position[b], str(b), ha="center", va="center", zorder=3)
        for length, color in colors.items():
            ax.plot(
                [],
                [],
                color=color,
                lw=2,
                label=f"{length/1000:g} km, {length/299792458*1e6:.2f} µs",
            )
        ax.plot([], [], "ks", mfc="white", label="Governed machine")
        ax.plot([], [], "k^", mfc="white", label="Converter in hybrid variant")
        ax.text(1.75, -0.62, "Fault shunt", ha="center")
        ax.text(1.75, 0.62, "Switched load", ha="center")
        ax.set_title(
            f'{count} buses, {len(meta["edges"])} distributed lines; transport delays shown'
        )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.legend(loc="lower left", bbox_to_anchor=(0, -0.2), ncol=2)
        save(fig, output / f"network{count}_topology.png")


def refinement(records, folder, output):
    cases = list(records)
    for level in ["coarse", "fine"]:
        for case in sorted((folder / level).glob("bergeron_sine_short_distributed*")):
            if (case / "statistics.json").exists():
                record = collect(case)
                if record:
                    cases.append(record)
    rows = []
    for meta, values, steps in cases:
        if not (
            meta["family"] == "bergeron"
            and meta["sine"]
            and meta["short"]
            and meta["kind"].startswith("distributed")
        ):
            continue
        w = 2 * np.pi * 60
        phasor = 1000 / (
            (1 + 100 / 600) * np.cos(w * meta["tau"])
            + 1j * (300 / 600 + 100 / 300) * np.sin(w * meta["tau"])
        )
        mask = values["t"] > 0.06
        error = values["Bus_b1_va"][mask] - np.real(
            phasor * np.exp(1j * w * values["t"][mask])
        )
        rows.append(
            {
                "rtol": meta["solver"]["rtol"],
                "atol": meta["solver"]["atol"],
                "kind": meta["kind"],
                "steps": len(steps),
                "max_h_over_tau": float(np.max(steps["step"]) / meta["tau"]),
                "rms_V": float(np.sqrt(np.mean(error**2))),
                "max_V": float(np.max(abs(error))),
                "source_folder": meta["source_folder"],
            }
        )
    if len({row["rtol"] for row in rows}) < 2:
        return
    (output / "refinement.json").write_text(json.dumps(rows, indent=2) + "\n")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for kind in ["distributed", "distributed_capped"]:
        values = sorted([r for r in rows if r["kind"] == kind], key=lambda r: r["rtol"])
        tolerance = [r["rtol"] for r in values]
        axes[0].loglog(
            tolerance,
            [r["rms_V"] for r in values],
            "o-",
            color=COLORS[kind],
            label=LABELS[kind],
        )
        axes[1].semilogx(
            tolerance,
            [r["max_h_over_tau"] for r in values],
            "o-",
            color=COLORS[kind],
            label=LABELS[kind],
        )
    axes[0].set_ylabel("Analytic steady voltage RMS error [V]")
    axes[1].set_ylabel("Largest accepted step / delay")
    for ax in axes:
        ax.set_xlabel("Relative tolerance")
        ax.grid(alpha=0.2)
    axes[0].legend()
    fig.suptitle("60 Hz, 20 µs lossless line: tolerance refinement")
    save(fig, output / "delay_refinement.png")


def frequency_response(records, fits, output):
    cases = [r for r in records if r[0]["family"] == "frequency"]
    if not cases:
        return
    # Bromwich inversion on a uniformly sampled vertical Laplace contour.
    # The exponential window suppresses circular wraparound by exp(-24).
    samples = 2**18
    spacing = 2e-7
    duration = samples * spacing
    sigma = 24 / duration
    omega = 2 * np.pi * np.fft.rfftfreq(samples, spacing)
    s = sigma + 1j * omega

    def zero_sequence(fit):
        value = np.full(s.shape, np.sum(fit["D"][0]), complex)
        for pole, residue in zip(fit["poles"], fit["residues"]):
            row = np.array(residue)[0]
            value += np.sum(row[:, 0] + 1j * row[:, 1]) / (s - complex(*pole))
        return value

    yc = zero_sequence(json.loads((fits / "yc.json").read_text()))
    modes = json.loads((fits / "h60000.json").read_text())["modes"]
    h = sum(zero_sequence(m["H"]) * np.exp(-s * m["tau"]) for m in modes)
    zc = 1 / yc
    denominator = (1 + 100 / 600) * (1 + h * h) + (zc / 600 + 100 / zc) * (1 - h * h)
    transfer = 2 * h / denominator
    time = np.arange(samples) * spacing
    reference = (
        np.fft.irfft(1000 * transfer * np.exp(-s * 0.001) / s, n=samples)
        / spacing
        * np.exp(sigma * time)
    )
    keep = time <= 0.008
    np.savetxt(
        output / "frequency_reference.csv",
        np.c_[time[keep], reference[keep]],
        delimiter=",",
        header="time,receiving_voltage",
        comments="",
    )
    fig, axes = plt.subplots(3, 1, figsize=(11, 9))
    for meta, d, steps in cases:
        kind = meta["kind"]
        t = d["t"]
        expected = np.interp(t, time, reference)
        error = d["Bus_b1_va"] - expected
        meta["frequency_domain_error"] = {
            "rms_V": float(np.sqrt(np.mean(error**2))),
            "maximum_V": float(np.max(abs(error))),
            "reference": "inverse Laplace transform of the fitted terminal transfer",
            "inverse_transform_spacing_s": spacing,
            "exponential_window": 24,
        }
        axes[0].plot(
            t * 1e3, d["Bus_b1_va"], color=COLORS[kind], label=LABELS[kind], lw=1
        )
        if kind.startswith("distributed"):
            axes[1].plot(t * 1e3, error, color=COLORS[kind], label=LABELS[kind], lw=0.8)
        axes[2].semilogy(
            steps["time"] * 1e3,
            steps["step"] / meta["tau"],
            color=COLORS[kind],
            label=LABELS[kind],
            lw=0.8,
        )
    axes[0].plot(
        time[keep] * 1e3,
        reference[keep],
        "k--",
        label="Frequency-domain reference",
        lw=0.8,
    )
    axes[0].set_xlim(0.9, 2.5)
    axes[0].set_ylabel("Receiving voltage [V]")
    axes[0].set_title("60 km frequency-dependent line: zero-sequence energization")
    axes[1].set_ylabel("Distributed minus reference [V]")
    axes[2].set_ylabel("Accepted step / delay")
    axes[2].axhline(1, color="k", ls="--", lw=0.8)
    for ax in axes:
        ax.set_xlabel("Time [ms]")
        ax.grid(alpha=0.2)
    axes[0].legend(ncol=2)
    save(fig, output / "frequency_step.png")


def main(args):
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    folders = [folder for root in args.results for folder in root.iterdir()]
    for folder in sorted(folders):
        if folder.is_dir() and (folder / "statistics.json").exists():
            record = collect(folder)
            if record:
                records.append(record)
    topology(records, args.output)
    refinement(records, args.refinement, args.output)
    bergeron(records, args.output)
    network(records, args.output)
    fits(args.fits, args.output)
    frequency_response(records, args.fits, args.output)
    summaries = [r[0] for r in records]
    (args.output / "statistics.json").write_text(json.dumps(summaries, indent=2) + "\n")
    lines = [
        "| Study | Horizon [ms] | Trials | CPU median [s] | DAE rows | Steps | Rejected error tests | h median / max [µs] | h > shortest delay [%] |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in summaries:
        s = r["solver"]
        h = r["internal_steps"]
        lines.append(
            f"| {r['name']} | {r['horizon_s']*1000:g} | {r['successful_trials']} | {s['cpu_s']:.4g} | {s['variables']:.0f} | {s['steps']:.0f} | {s['error_test_fails']:.0f} | {h['median_s']*1e6:.3g} / {h['maximum_s']*1e6:.3g} | {h['percent_above_shortest_delay']:.3g} |"
        )
    (args.output / "statistics.md").write_text("\n".join(lines) + "\n")
    networks = [r for r in summaries if r["family"] == "network"]
    if networks:
        fig, axes = plt.subplots(1, 2, figsize=(13, max(5, len(networks) * 0.3)))
        labels = [r["name"].replace("network", "").replace("_", " ") for r in networks]
        for ax, key, title in zip(
            axes,
            ["cpu_s", "steps"],
            ["Median process CPU [s]", "Accepted internal steps"],
        ):
            ax.barh(
                labels,
                [r["solver"][key] for r in networks],
                color=[COLORS[r["kind"]] for r in networks],
            )
            ax.set_xscale("log")
            ax.set_xlabel(title)
            ax.grid(axis="x", alpha=0.2)
        save(fig, args.output / "network_cost.png")

    def flatten(value, prefix=""):
        result = {}
        for key, item in value.items():
            name = prefix + key
            if isinstance(item, dict):
                result.update(flatten(item, name + "."))
            else:
                result[name] = json.dumps(item) if isinstance(item, list) else item
        return result

    flat = [flatten(r) for r in summaries]
    columns = list(dict.fromkeys(key for row in flat for key in row))
    with (args.output / "statistics.csv").open("w", newline="") as stream:
        writer = csv_module.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(flat)
    cards = []
    for path in sorted(args.output.glob("*.png")):
        label = html.escape(path.stem.replace("_", " "))
        name = html.escape(path.name)
        cards.append(
            f'<article><h2>{label}</h2><a href="{name}"><img src="{name}" loading="lazy" alt="{label}"></a></article>'
        )
    page = (
        """<!doctype html><html lang="en"><meta charset="utf-8"><title>Distributed EMT studies</title>
<style>body{font:16px system-ui;margin:2rem auto;max-width:1500px;padding:0 1rem;color:#213244;background:#f5f7fa}a{color:#176b91}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(500px,1fr));gap:1.4rem}article{background:white;padding:1rem;border:1px solid #dce2e8;border-radius:8px}h2{font-size:1.1rem}img{max-width:100%;height:auto}input{padding:.7rem;margin:1rem 0;width:20rem;max-width:90%}@media(max-width:550px){main{display:block}article{margin-bottom:1rem}}</style>
<h1>Distributed EMT line studies</h1><p>Adaptive stepping throughout. Select a plot to open its full-size image.</p>
<p><a href="../REPORT.md">Report and architecture assessment</a> · <a href="statistics.csv">All statistics (CSV)</a> · <a href="statistics.json">All statistics (JSON)</a></p>
<label>Filter plots <input id="filter" type="search" placeholder="20, fault, converter, Bergeron…"></label><main>"""
        + "".join(cards)
        + """</main>
<script>document.getElementById('filter').addEventListener('input',function(){let q=this.value.toLowerCase();document.querySelectorAll('article').forEach(a=>a.hidden=!a.textContent.toLowerCase().includes(q));});</script></html>"""
    )
    (args.output / "index.html").write_text(page)
    print(f"Wrote {len(records)} study summaries and plots to {args.output}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[3]
    p.add_argument(
        "--results",
        type=Path,
        nargs="+",
        default=[root / "local/emt-distributed/final-studies"],
    )
    p.add_argument(
        "--fits", type=Path, default=root / "local/emt-distributed/fits-final"
    )
    p.add_argument(
        "--refinement", type=Path, default=root / "local/emt-distributed/refinement"
    )
    p.add_argument("--output", type=Path, default=root / "local/emt-distributed/plots")
    main(p.parse_args())
