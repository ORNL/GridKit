#!/usr/bin/env python3
"""Build a 14-bus hybrid EMT case from OpenLine's coupled line fit."""

import json
from pathlib import Path

import numpy as np
from scipy.optimize import root

from pwm_analysis import pwm_peak_coefficients

HERE = Path(__file__).resolve().parent
VLL, F, FC, M = 13800.0, 60.0, 900.0, .8
W = 2 * np.pi * F
PHASES = np.array([0., -2 * np.pi / 3, 2 * np.pi / 3])
LINES = [(1, 2, 350), (2, 3, 500), (3, 4, 450), (4, 5, 600),
         (5, 6, 500), (6, 7, 450), (7, 8, 400), (8, 1, 550),
         (3, 9, 300), (4, 10, 400), (6, 11, 350), (7, 12, 450),
         (2, 13, 300)]
MACHINES = {3: (5e6, 3.7, 1.2e6), 7: (4e6, 4.2, 1e6)}
LOADS = {2: .6e6, 4: .7e6, 5: .65e6, 6: .55e6, 8: .5e6,
         10: .6e6, 12: .5e6, 14: .8e6}
MU = [4 * F, float(np.sqrt(4 * F * 4 * FC)), 4 * FC]


def write_json(path, value):
    path.write_text(json.dumps(value, indent=4) + "\n")


def diagonal(value):
    return (np.eye(3) * value).tolist()


def rl_admittance(r, l):
    """Exact three-phase RL admittance in the VectorFit representation."""
    residue = np.eye(3) / l
    return dict(rows=3, cols=3, D=diagonal(0.), E=diagonal(0.),
                poles=[[-r / l, 0.]],
                residues=[np.stack((residue, np.zeros((3, 3))), axis=-1).tolist()])


def response(model, s):
    y = np.array(model["D"], complex) + s * np.array(model["E"])
    for pole, residue in zip(model["poles"], model["residues"]):
        r = np.array(residue)
        y += (r[..., 0] + 1j * r[..., 1]) / (s - complex(*pole))
    return y


def main():
    line_model = json.loads((HERE / "lines/model.json").read_text())
    z, ysh = (response(line_model[key], 1j * W) for key in ("Zp", "Yp"))
    # DC voltage and every physical parameter stay fixed in the mu sweep.
    gate_peak = pwm_peak_coefficients(MU[-1], max_harmonic=1)[0]["gate_peak"]
    dc = 1.05 * VLL * np.sqrt(2 / 3) / gate_peak
    e_ibr = 1.05 * VLL / np.sqrt(3) * np.exp(1j * (-np.pi / 2 + PHASES))
    grid_angle = -np.pi / 2 - .07
    e_grid = VLL / np.sqrt(3) * np.exp(1j * (grid_angle + PHASES))
    grid_r, grid_l = .12, .001
    filter_r = .02 * VLL ** 2 / .75e6
    filter_l = .15 * VLL ** 2 / .75e6 / W

    # Unbalanced linear loads, specified by nominal power at rated voltage.
    # Phase permutations give unequal mutual-coupling excitation throughout.
    load_z = {}
    for bus, p in LOADS.items():
        phase_p = p / 3 * np.roll([1.06, .98, .96], bus % 3)
        q = phase_p * np.tan(np.arccos(.97)) if bus != 14 else 0.
        load_z[bus] = (VLL ** 2 / 3) / (phase_p - 1j * q)

    # Independent full three-phase fundamental network estimate. Machine
    # terminals impose balanced voltage; solve two angles for dispatch.
    admittance = np.zeros((39, 39), complex)
    rhs = np.zeros(39, complex)
    sl = lambda bus: slice(3 * (bus - 1), 3 * bus)
    for a, b, length in LINES:
        series = np.linalg.inv(z * length)
        admittance[sl(a), sl(a)] += series + ysh * length / 2
        admittance[sl(b), sl(b)] += series + ysh * length / 2
        admittance[sl(a), sl(b)] -= series
        admittance[sl(b), sl(a)] -= series
    for bus, impedance in load_z.items():
        if bus != 14:
            admittance[sl(bus), sl(bus)] += np.diag(1 / impedance)
    for bus, emf, impedance in [(1, e_grid, grid_r + 1j * W * grid_l)] + [
            (bus, e_ibr, filter_r + 1j * W * filter_l) for bus in (9, 11, 13)]:
        admittance[sl(bus), sl(bus)] += np.eye(3) / impedance
        rhs[sl(bus)] += emf / impedance
    fixed = np.array([i for bus in MACHINES for i in range(sl(bus).start, sl(bus).stop)])
    free = np.array([i for i in range(39) if i not in fixed])

    def network(angles):
        v = np.zeros(39, complex)
        v[fixed] = np.concatenate([VLL / np.sqrt(3) * np.exp(1j * (angle + PHASES))
                                   for angle in angles])
        v[free] = np.linalg.solve(admittance[np.ix_(free, free)],
                                  rhs[free] - admittance[np.ix_(free, fixed)] @ v[fixed])
        current = admittance @ v - rhs
        powers = np.array([np.sum(v[sl(bus)] * np.conj(current[sl(bus)])) for bus in MACHINES])
        return v, powers

    desired = np.array([values[2] for values in MACHINES.values()])
    solution = root(lambda angles: (network(angles)[1].real - desired) / 1e6,
                    [grid_angle, grid_angle])
    if not solution.success or np.max(np.abs(network(solution.x)[1].real - desired)) > 1.:
        raise RuntimeError(f"Fundamental network initialization failed: {solution.message}")
    v, power = network(solution.x)

    case = {"header": {"case_name": "CoupledGrid: 14-bus hybrid EMT network",
                       "case_description": "Two governed and excited machines, three PWM converter sources, and untransposed frequency dependent overhead lines",
                       "case_comments": "Synthetic 13.8 kV network using physical OpenLine conductor geometry. Fixed DC voltage across all mu studies. See README.md."},
            "signals": [], "devices": []}
    devices = case["devices"]
    for bus in range(1, 15):
        devices.append(dict(id=f"bus_{bus}", **{"class": "Bus"}, mon=["va", "vb", "vc"]))
    devices.append({"class": "VoltageSource", "id": "source_1",
                    "params": {"E": [VLL / np.sqrt(3)] * 3,
                               "phi": (grid_angle + PHASES).tolist(), "omega": W},
                    "submodels": {"Y": rl_admittance(grid_r, grid_l)},
                    "inputs": {"bus": "bus_1"}, "mon": ["ea", "eb", "ec", "ia", "ib", "ic"]})
    machine_params = dict(N=3, V=VLL, f=F, F=0.0, Rs=.003, Ll=.15,
                          Lmd=1.66, Lmq=1.61, L0=.15, Rfd=.0006, Llfd=.165,
                          R1d=.0284, Ll1d=.1713, R1q=.0062, Ll1q=.7252,
                          R2q=.0237, Ll2q=.125, S10=.1, S12=.5)
    for bus, (rating, inertia, _) in MACHINES.items():
        case["signals"] += [{"id": f"{name}_{bus}"} for name in ("speed", "pmech", "efd")]
        devices += [
            {"class": "Machine", "id": f"machine_{bus}",
             "params": dict(machine_params, S=rating, H=inertia),
             "inputs": {"bus": f"bus_{bus}", "pm": f"pmech_{bus}", "efd": f"efd_{bus}"},
             "outputs": {"speed": f"speed_{bus}"},
             "mon": ["theta", "omega", "te", "ifd", "efd", "ks", "psi_at", "ia", "ib", "ic", "p", "q"]},
            {"class": "Tgov1", "id": f"governor_{bus}",
             "params": dict(Trate=rating / 1e6, R=.05, T1=.2, T2=.5, T3=2., Pvmax=1., Pvmin=0., Dt=0.),
             "inputs": {"speed": f"speed_{bus}"}, "outputs": {"pmech": f"pmech_{bus}"}},
            {"class": "Ieeet1", "id": f"exciter_{bus}",
             "params": dict(V=VLL, Tr=.02, Ka=20., Ta=.05, Vrmax=5., Vrmin=-5.,
                            Ke=1., Te=.5, Kf=.05, Tf=1., E1=2.8, Se1=.04, E2=3.73, Se2=.33, Ispdlim=0.),
             "inputs": {"bus": f"bus_{bus}", "speed": f"speed_{bus}"},
             "outputs": {"efd": f"efd_{bus}"}, "mon": ["efd", "vts", "vr"]}]
    for bus in (9, 11, 13):
        gates, emf = ([f"{name}{p}_{bus}" for p in "abc"] for name in ("s", "e"))
        case["signals"] += [{"id": f"dc_{bus}", "value": dc}] + [{"id": s} for s in gates + emf]
        devices += [
            {"class": "PWM", "id": f"pwm_{bus}", "params": dict(M=M, fm=F, fc=FC, alignment=.5),
             "outputs": {"s": gates}, "mon": ["s"]},
            {"class": "Converter", "id": f"converter_{bus}", "inputs": {"s": gates, "vdc": f"dc_{bus}"},
             "outputs": {"vo": emf}, "mon": ["vo"]},
            {"class": "DependentVoltageSource", "id": f"filter_{bus}",
             "submodels": {"Y": rl_admittance(filter_r, filter_l)},
             "inputs": dict(bus=f"bus_{bus}", **dict(zip(("ea", "eb", "ec"), emf))),
             "mon": ["ea", "eb", "ec", "ia", "ib", "ic"]}]
    for a, b, length in LINES:
        devices.append({"class": "LineLumped", "id": f"line_{a}_{b}",
                        "params": {"N": 3, "K": 3, "conductors": [1, 2, 3], "dx": float(length)},
                        "submodels": line_model, "inputs": {"bus1": f"bus_{a}", "bus2": f"bus_{b}"},
                        "mon": ["i12a", "i12b", "i12c", "i_sh1a", "i_sh1b", "i_sh1c", "i_sh2a", "i_sh2b", "i_sh2c"]})
    for bus, impedance in load_z.items():
        params = {"R": np.diag(impedance.real).tolist(), "L": np.diag(impedance.imag / W).tolist()}
        devices.append({"class": "LoadZ", "id": f"load_{bus}", "params": params,
                        "inputs": {"bus": f"bus_{bus}"}, "mon": ["ia", "ib", "ic"]})
    devices.append({"class": "Switch", "id": "load_step", "params": {"open": True},
                    "inputs": {"bus1": "bus_10", "bus2": "bus_14"},
                    "mon": ["open", "i12a", "i12b", "i12c"]})
    state = {"header": {"version": 1, "time": 0., "description": "Common high-mu fundamental estimate; internal line/filter states start at model defaults, so energization transients remain."},
             "buses": {}, "devices": {}}
    for bus in range(1, 15):
        voltage = v[sl(bus)] if bus != 14 else np.zeros(3)
        state["buses"][f"bus_{bus}"] = dict(zip(("va", "vb", "vc"), (np.sqrt(2) * voltage.real).tolist()))
    for bus, p in zip(MACHINES, power):
        state["devices"][f"machine_{bus}"] = {"p": float(p.real), "q": float(p.imag)}
    write_json(HERE / "CoupledGrid.case.json", case)
    write_json(HERE / "CoupledGrid.state.json", state)
    for mu, name in zip(MU, ("low", "middle", "high")):
        write_json(HERE / f"{name}.solver.json", {
            "system_model_file": "CoupledGrid.case.json", "state_file": "CoupledGrid.state.json",
            "dt_monitor": 1 / 60000, "tmax": 1.2, "rel_tol": 1e-6, "abs_tol": 1e-7,
            "mu": mu, "max_steps": 2000000, "consistent_ic_type": "ya_ydp",
            "events": [{"time": .5, "type": "switch_close", "element_id": "load_step"},
                       {"time": .8, "type": "switch_open", "element_id": "load_step"}],
            "output_file": f"{name}.csv"})
    print(f"Fixed DC voltage: {dc:.9f} V; filter R={filter_r:.6f} ohm, L={filter_l:.9f} H")
    print("Machine initial P+jQ [MVA]:", power / 1e6)


if __name__ == "__main__":
    main()
