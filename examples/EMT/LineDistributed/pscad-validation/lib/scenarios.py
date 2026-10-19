"""The validation scenarios, expressed as GridKit LineDistributed run objects.

Every scenario is the example circuit: a 60 Hz three-phase source behind a
short-circuit impedance (X/R = 10), the line, and a grounded-wye load, with
the example's finite C2 transitions for every event. The same objects drive
the PSCAD forcing script, so the two simulators integrate identical inputs.
"""

import math

CODES = {"69kv-wood-pole": "l69", "138kv-delta": "l138", "345kv-h-frame": "l345", "765kv-horizontal": "l765"}
MODELS = {"fd": "Frequency-dependent (GridKit vector fits / PSCAD FD phase model)",
          "bg": "Lossless Bergeron (constant 60 Hz L and C in both tools)"}
#: Fixed PSCAD steps per scenario: EMTDC integration step and recording interval, seconds.
PSCAD_STEPS = {"energize": (1e-6, 2e-6), "fault": (2e-6, 1e-5), "pulse": (5e-7, 1e-6), "sag": (2e-6, 1e-5)}
#: PSCAD start-up before study time zero for sinusoidally initialized scenarios.
STARTUP = 1.0
PHASES = [0.0, -2 * math.pi / 3, 2 * math.pi / 3]


def round_to(value, digits=1):
    return float(f"{value:.{digits}g}")


def source(sequence):
    kv, mva = sequence["voltage_kv"], sequence["short_circuit_mva"]
    z = (kv * 1e3) ** 2 / (mva * 1e6)
    r = z / math.sqrt(101.0)
    return {"rms": [kv * 1e3 / math.sqrt(3)] * 3, "phase": PHASES, "frequency": 60.0,
            "resistance": [r] * 3, "inductance": [10 * r / (2 * math.pi * 60.0)] * 3}


def rl_load(sequence, model):
    """2 Zc per phase: at 0.95 power factor for the frequency-dependent tier, resistive
    for the lossless tier. A lossless line between an inductive source and an
    inductive load reflects every high-frequency component fully at both ends, so
    nothing damps numerical high-frequency content and the adaptive step collapses;
    a resistive load absorbs two thirds of each reflection."""
    r = 2 * sequence["surge_impedance_ohm"]["positive"]
    inductance = r * math.tan(math.acos(0.95)) / (2 * math.pi * 60.0) if model == "fd" else 0.0
    return {"resistance": [r] * 3, "inductance": [inductance] * 3}


def solver(sequence, tf, dt):
    v_peak = math.sqrt(2) * sequence["voltage_kv"] * 1e3 / math.sqrt(3)
    scale = v_peak / sequence["surge_impedance_ohm"]["positive"]
    # History compression is on: without it, fronts recirculating between two
    # near-reflecting terminations (open line, inductive source) spawn a new
    # history segment per accepted step every round trip, and the step count grows
    # without bound. Merged segments spend at most a quarter of the tolerance.
    return {"t0": 0.0, "tf": tf, "dt": dt, "rel_tol": 1e-6, "history_abs_tol": round_to(1e-6 * scale), "h_max": 0.0,
            "compress": True}


def delay(sequence):
    return max(sequence["bergeron_delays_s"]["zero"], sequence["bergeron_delays_s"]["positive"],
               sequence["wavefront_delays_s"]["ground"], sequence["wavefront_delays_s"]["aerial"])


def scenarios(name, sequence, model="fd"):
    """The four runs of one line and line model, keyed by scenario name."""
    zc = sequence["surge_impedance_ohm"]["positive"]
    v_peak = math.sqrt(2) * sequence["voltage_kv"] * 1e3 / math.sqrt(3)
    tau = delay(sequence)
    three_phase = name in ("138kv-delta", "765kv-horizontal")
    nominal = source(sequence)
    load = rl_load(sequence, model)
    load_text = "R-L load (2 Zc, pf 0.95)" if model == "fd" else "resistive load (2 Zc)"
    return {
        "energize": {
            "name": "energize",
            "description": "Energization of the open-ended line from zero state: the source amplitude rises from 0 to 1 "
                           "over 50 us at 1 ms, launching a steep front that reflects off the open receiving end.",
            "source": nominal, "load": {"resistance": [1e5] * 3, "inductance": [0.0] * 3},
            "initialization": "zero",
            "events": [{"type": "source", "value": 1.0, "on": 1e-3, "off": 1.0, "rise": 5e-5, "fall": 5e-5, "phases": [1, 2, 3]}],
            "solver": solver(sequence, 1e-3 + max(0.03, 12 * tau), 2e-6),
        },
        "fault": {
            "name": "fault",
            "description": f"Sinusoidal steady state into the {load_text}, then a "
                           f"{'three-phase' if three_phase else 'phase-A'} {zc / 50:.1f} ohm receiving-end fault applied at 20 ms "
                           "(50 us transition) and removed over 2 ms at 60 ms.",
            "source": nominal, "load": load,
            "initialization": "sinusoidal",
            "events": [{"type": "fault", "value": 50 / zc, "on": 0.02, "off": 0.06, "rise": 5e-5, "fall": 2e-3,
                        "phases": [1, 2, 3] if three_phase else [1]}],
            "solver": solver(sequence, 0.12, 1e-5),
        },
        "pulse": {
            "name": "pulse",
            "description": f"Surge test from zero state: a phase-A pulse of {0.1 * v_peak / 1e3:.1f} kV (5 us rise, 50 us FWHM) "
                           f"from a 30 ohm + 20 uH source into a {1.2 * zc:.0f} ohm termination.",
            "source": {"rms": [0.0] * 3, "phase": PHASES, "frequency": 60.0, "resistance": [30.0] * 3, "inductance": [2e-5] * 3},
            "load": {"resistance": [1.2 * zc] * 3, "inductance": [0.0] * 3},
            "initialization": "zero",
            "events": [{"type": "voltage", "value": 0.1 * v_peak, "on": 1e-3, "off": 1.05e-3, "rise": 5e-6, "fall": 5e-6, "phases": [1]}],
            "solver": solver(sequence, 1e-3 + max(2.5e-3, 5 * tau), 1e-6),
        },
        "sag": {
            "name": "sag",
            "description": f"Sinusoidal steady state into the {load_text}; the phase-A source amplitude sags by 40% from 20 to "
                           "50 ms (50 us transitions). Smooth operation lets the adaptive solver step over the line delay.",
            "source": nominal, "load": load,
            "initialization": "sinusoidal",
            "events": [{"type": "source", "value": -0.4, "on": 0.02, "off": 0.05, "rise": 5e-5, "fall": 5e-5, "phases": [1]}],
            "solver": solver(sequence, 0.1, 1e-5),
        },
    }


def project_name(line, model, scenario, refined=False):
    return f"V_{CODES[line]}_{model}_{scenario}" + ("_h" if refined else "")
