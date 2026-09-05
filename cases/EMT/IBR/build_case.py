#!/usr/bin/env python3
"""Rebuild the synthetic case and its fundamental-frequency initial estimate."""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
VLL, FREQUENCY, MU = 13800.0, 60.0, 240.0
OMEGA = 2 * np.pi * FREQUENCY
M, FC = 0.8, 900.0


def diagonal(value):
    return (np.eye(3) * value).tolist()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=4) + "\n")


def main():
    # Exact fundamental coefficient of the ideal pulse train, followed by the
    # Fourier transform of the logistic smoothing kernel. This estimates only
    # the operating point; all time-domain results come from GridKit/IDA.
    k = np.arange(round(FC / FREQUENCY))
    duty = (1 + M * np.sin(OMEGA * (k + 0.5) / FC)) / 2
    on = (k + 0.5 * (1 - duty)) / FC
    off = (k + 0.5 + 0.5 * duty) / FC
    coefficient = np.sum(np.exp(-1j * OMEGA * on) - np.exp(-1j * OMEGA * off)) / (1j * OMEGA / FREQUENCY)
    attenuation = (np.pi * OMEGA / MU) / np.sinh(np.pi * OMEGA / MU)
    unit_peak = 2 * coefficient * attenuation
    dc = float(1.02 * VLL * np.sqrt(2 / 3) / abs(unit_peak))
    source_voltage = dc * unit_peak / np.sqrt(2)
    machine_voltage = VLL / np.sqrt(3) * np.exp(1j * (np.angle(source_voltage) - 0.025))

    # SI series R [ohm], L [H]; dx=1 makes the per-length matrices equal totals.
    lines = [(1, 7, .12, .0015), (2, 7, .15, .002), (3, 8, .12, .0015),
             (4, 7, .10, .001), (5, 8, .10, .001), (6, 8, .14, .0015),
             (1, 2, .30, .003), (2, 3, .30, .003), (7, 8, .50, .004)]
    loads = {4: 380.88, 5: 380.88, 6: 380.88, 7: 31.74, 8: 38.088, 9: 95.22, 10: 2.0}
    source_r, source_l = 1.0, .016
    # Buses 7 and 8 are initially joined by the closed ideal tie switch.
    group = [0, 1, 2, 3, 4, 5, 6, 6]
    admittance = np.zeros((7, 7), complex)
    rhs = np.zeros(7, complex)
    for a, b, r, l in lines:
        a, b = group[a - 1], group[b - 1]
        if a == b:
            continue
        y = 1 / (r + 1j * OMEGA * l)
        admittance[a, a] += y
        admittance[b, b] += y
        admittance[a, b] -= y
        admittance[b, a] -= y
        if a < 3 and b < 3:
            admittance[a, a] += 0.5j * OMEGA * 1e-6
            admittance[b, b] += 0.5j * OMEGA * 1e-6
    for bus, resistance in loads.items():
        if bus <= 8:
            admittance[group[bus - 1], group[bus - 1]] += 1 / resistance
    for bus in (4, 5, 6):
        a = group[bus - 1]
        y = 1 / (source_r + 1j * OMEGA * source_l)
        admittance[a, a] += y
        rhs[a] += y * source_voltage
    v = np.full(7, machine_voltage, complex)
    v[3:] = np.linalg.solve(admittance[3:, 3:], rhs[3:] - admittance[3:, :3] @ v[:3])
    machine_power = 3 * v[:3] * np.conj((admittance @ v - rhs)[:3])
    voltages = [v[g] for g in group] + [0j, 0j]

    case = {"header": {"case_name": "Synthetic EMT 10-bus hybrid grid",
                       "case_description": "Three governed synchronous machines, three PWM converters, seven resistive loads and three event switches",
                       "case_comments": "13.8 kV, 60 Hz synthetic demonstration; open-loop ideal DC sources; shared PWM mu=240 1/s. See README.md for smoothing compensation and limitations."},
            "signals": [], "devices": []}
    devices = case["devices"]
    for bus in range(1, 11):
        devices.append({"class": "Bus", "id": f"bus_{bus}", "mon": ["va", "vb", "vc"]})
    machine_params = dict(N=3, S=10e6, V=VLL, f=FREQUENCY, F=0.0, Rs=.003,
                          Ll=.15, Lmd=1.66, Lmq=1.61, L0=.15, Rfd=.0006,
                          Llfd=.165, R1d=.0284, Ll1d=.1713, R1q=.0062,
                          Ll1q=.7252, R2q=.0237, Ll2q=.125, S10=.1, S12=.5)
    for bus, inertia in zip((1, 2, 3), (3.7, 4.5, 3.0)):
        case["signals"] += [{"id": f"speed_{bus}"}, {"id": f"pmech_{bus}"}]
        devices.append({"class": "Machine", "id": f"machine_{bus}",
                        "params": dict(machine_params, H=inertia),
                        "inputs": {"bus": f"bus_{bus}", "pm": f"pmech_{bus}"},
                        "outputs": {"speed": f"speed_{bus}"},
                        "mon": ["theta", "omega", "te", "ifd", "efd", "ks", "psi_at", "ia", "ib", "ic", "p", "q", "id", "iq"]})
        devices.append({"class": "Tgov1", "id": f"governor_{bus}",
                        "params": dict(Trate=10.0, R=.05, T1=.2, T2=.5, T3=2.0, Pvmax=1.0, Pvmin=0.0, Dt=0.0),
                        "inputs": {"speed": f"speed_{bus}"}, "outputs": {"pmech": f"pmech_{bus}"}})
    for bus in (4, 5, 6):
        gates = [f"s{p}_{bus}" for p in "abc"]
        emf = [f"e{p}_{bus}" for p in "abc"]
        case["signals"] += [{"id": f"dc_{bus}", "value": dc}] + [{"id": s} for s in gates + emf]
        devices += [
            {"class": "PWM", "id": f"pwm_{bus}", "params": dict(M=M, fm=FREQUENCY, fc=FC, alignment=.5), "outputs": {"s": gates}, "mon": ["s"]},
            {"class": "Converter", "id": f"converter_{bus}", "inputs": {"s": gates, "vdc": f"dc_{bus}"}, "outputs": {"vo": emf}, "mon": ["vo"]},
            {"class": "DependentVoltageSource", "id": f"filter_{bus}",
             "params": {"Rs": diagonal(source_r), "Ls": diagonal(source_l)},
             "inputs": dict(bus=f"bus_{bus}", **dict(zip(("ea", "eb", "ec"), emf))),
             "mon": ["ea", "eb", "ec", "ia", "ib", "ic"]}]
    for a, b, r, l in lines:
        devices.append({"class": "LineLumped", "id": f"line_{a}_{b}",
                        "params": {"N": 3, "K": 3, "conductors": [1, 2, 3], "dx": 1.0,
                                   "Rp": diagonal(r), "Lp": diagonal(l), "Gp": diagonal(0.0), "Cp": diagonal(1e-6 if a <= 3 and b <= 3 else 0.0)},
                        "inputs": {"bus1": f"bus_{a}", "bus2": f"bus_{b}"},
                        "mon": ["i12a", "i12b", "i12c", "i_sh1a", "i_sh1b", "i_sh1c", "i_sh2a", "i_sh2b", "i_sh2c"]})
    for bus, resistance in loads.items():
        devices.append({"class": "LoadZ", "id": f"load_{bus}", "params": {"R": diagonal(resistance)},
                        "inputs": {"bus": f"bus_{bus}"}, "mon": ["ia", "ib", "ic"]})
    for name, a, b, is_open in [("tie", 7, 8, False), ("load_step", 8, 9, True), ("fault", 7, 10, True)]:
        devices.append({"class": "Switch", "id": name, "params": {"open": is_open},
                        "inputs": {"bus1": f"bus_{a}", "bus2": f"bus_{b}"}, "mon": ["open", "i12a", "i12b", "i12c"]})
    state = {"header": {"version": 1, "time": 0.0, "description": "Fundamental phasor estimate; line/filter currents are initialized by the models, so startup transients remain."},
             "buses": {}, "devices": {}}
    for bus, voltage in enumerate(voltages, 1):
        state["buses"][f"bus_{bus}"] = dict(zip(("va", "vb", "vc"), (np.sqrt(2) * np.real(voltage * np.exp(1j * np.array([0, -2*np.pi/3, 2*np.pi/3])))).tolist()))
    for bus, power in enumerate(machine_power, 1):
        state["devices"][f"machine_{bus}"] = {"p": float(power.real), "q": float(power.imag)}
    write_json(HERE / "TenBus.case.json", case)
    write_json(HERE / "TenBus.state.json", state)
    print(f"Smoothing fundamental gain: {attenuation:.8f}; effective DC: {dc:.3f} V")
    print("Machine initial P+jQ [MVA]:", machine_power / 1e6)


if __name__ == "__main__":
    main()
