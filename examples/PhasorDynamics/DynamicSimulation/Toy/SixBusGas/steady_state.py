"""Write the SixBusGas case in steady state.

The power flow gives each machine's mechanical power P_m. A machine on
GasInputPm burns k_f * P_m of gas, and the gas-network steady state with those
withdrawals gives the gas-bus pressures. Both are written to the case file,
with the gas law derived from T, G, and the equation of state as in GasTranSim.

Usage: python3 steady_state.py [case_file]
"""

import json
import pathlib
import sys

import numpy as np
from scipy.optimize import fsolve

S_SYS = 100.0e6  # System power base [VA]
F_SYS = 60.0  # System frequency base [Hz]

# Electrical network, per unit on S_SYS
SLACK_BUS = 1
V_SET = {1: 1.00, 2: 1.02, 3: 1.02}
LOADS = {4: (1.0, 0.3), 5: (1.5, 0.5), 6: (1.2, 0.4)}
BRANCHES = [  # id, bus1, bus2, R, X, B
    ("BR_1_4", 1, 4, 0.01, 0.10, 0.02),
    ("BR_1_6", 1, 6, 0.01, 0.10, 0.02),
    ("BR_4_5", 4, 5, 0.01, 0.08, 0.02),
    ("BR_5_6", 5, 6, 0.01, 0.08, 0.02),
    ("BR_2_5", 2, 5, 0.00, 0.06, 0.00),
    ("BR_3_6", 3, 6, 0.00, 0.06, 0.00),
]

# Gas-fired GENROU machines: bus -> (id, active power, MVA base, gas bus)
MACHINES = {
    2: ("G2", 1.6, 200.0, 103),
    3: ("G3", 1.2, 150.0, 104),
}
GENROU = {
    "H": 4.0,
    "D": 0.0,
    "Ra": 0.0,
    "Tdop": 6.0,
    "Tdopp": 0.04,
    "Tqop": 0.8,
    "Tqopp": 0.06,
    "Xd": 1.8,
    "Xdp": 0.3,
    "Xdpp": 0.25,
    "Xq": 1.7,
    "Xqp": 0.55,
    "Xqpp": 0.25,
    "Xl": 0.15,
    "S10": 0.0,
    "S12": 0.0,
}
KF = 0.06  # Specific fuel consumption [kg/MJ]

# Gas network, SI units
SUPPLY_BUS = 101
SUPPLY_PRESSURE = 5.0e6  # [Pa]
GAS_BUSES = [101, 102, 103, 104]
GAS_BRANCHES = [  # id, bus1, bus2, length [m], diameter [m]
    ("GP_101_102", 101, 102, 5000.0, 0.30),
    ("GP_102_103", 102, 103, 2000.0, 0.25),
    ("GP_102_104", 102, 104, 3000.0, 0.25),
]
FRICTION = 0.01  # Darcy friction factor [-]
TEMPERATURE = 288.706  # [K]
SPECIFIC_GRAVITY = 0.6  # [-]
EOS = "simple_cnga"  # ideal, simple_cnga, or full_cnga
DX_MAX = 250.0  # Maximum segment length [m]


def gas_law(temperature, specific_gravity, eos):
    """Return the sound speed a [m/s] and density coefficients b1 [-], b2 [1/Pa]."""
    a = np.sqrt(8.314 * temperature / (0.02896 * specific_gravity))
    if eos == "ideal":
        return a, 1.0, 0.0
    if eos == "simple_cnga":
        return a, 1.00300865, 2.96848838e-8
    c = 344400.0 * 10.0 ** (1.785 * specific_gravity) / (6894.75729 * (1.8 * temperature) ** 3.825)
    return a, 1.0 + c * 101350.0, c


SOUND_SPEED, B1, B2 = gas_law(TEMPERATURE, SPECIFIC_GRAVITY, EOS)


def solve_power_flow():
    """Return complex bus voltages and machine power injections."""
    buses = sorted({b for _, b1, b2, *_ in BRANCHES for b in (b1, b2)})
    index = {bus: k for k, bus in enumerate(buses)}

    Y = np.zeros((len(buses), len(buses)), dtype=complex)
    for _, b1, b2, r, x, b in BRANCHES:
        i, j = index[b1], index[b2]
        y = 1.0 / complex(r, x)
        Y[i, i] += y + 0.5j * b
        Y[j, j] += y + 0.5j * b
        Y[i, j] -= y
        Y[j, i] -= y

    p_set = np.zeros(len(buses))
    q_set = np.zeros(len(buses))
    for bus, (_, p, _, _) in MACHINES.items():
        p_set[index[bus]] += p
    for bus, (p, q) in LOADS.items():
        p_set[index[bus]] -= p
        q_set[index[bus]] -= q

    magnitude = np.ones(len(buses))
    for bus, v in V_SET.items():
        magnitude[index[bus]] = v
    angle = np.zeros(len(buses))
    with_angle = [index[b] for b in buses if b != SLACK_BUS]
    with_magnitude = [index[b] for b in buses if b not in V_SET]

    def mismatch(x):
        angle[with_angle] = x[: len(with_angle)]
        magnitude[with_magnitude] = x[len(with_angle) :]
        v = magnitude * np.exp(1j * angle)
        s = v * np.conj(Y @ v)
        return np.concatenate([s.real[with_angle] - p_set[with_angle], s.imag[with_magnitude] - q_set[with_magnitude]])

    x0 = np.concatenate([angle[with_angle], magnitude[with_magnitude]])
    x = fsolve(mismatch, x0, xtol=1.0e-12)
    assert np.max(np.abs(mismatch(x))) < 1.0e-10, "power flow did not converge"

    v = magnitude * np.exp(1j * angle)
    s = v * np.conj(Y @ v)
    voltage = {bus: v[index[bus]] for bus in buses}
    injection = {bus: s[index[bus]] + complex(*LOADS.get(bus, (0.0, 0.0))) for bus in buses}
    return voltage, injection


def potential(p):
    """Integral of the density over pressure from zero to p."""
    return p**2 * (B1 / 2.0 + B2 * p / 3.0) / SOUND_SPEED**2


def pipe_mass_flow(p1, p2, length, diameter):
    """Steady mass flow from bus 1 to bus 2 [kg/s]."""
    area = np.pi * diameter**2 / 4.0
    drop = potential(p1) - potential(p2)
    return area * np.sign(drop) * np.sqrt(2.0 * diameter * abs(drop) / (FRICTION * length))


def solve_gas_network(withdrawal):
    """Return steady gas-bus pressures for the given withdrawals [kg/s]."""
    unknown = [bus for bus in GAS_BUSES if bus != SUPPLY_BUS]

    def pressures(x):
        p = dict(zip(unknown, x))
        p[SUPPLY_BUS] = SUPPLY_PRESSURE
        return p

    def balance(x):
        p = pressures(x)
        q = {bus: -withdrawal.get(bus, 0.0) for bus in unknown}
        for _, b1, b2, length, diameter in GAS_BRANCHES:
            flow = pipe_mass_flow(p[b1], p[b2], length, diameter)
            if b1 in q:
                q[b1] -= flow
            if b2 in q:
                q[b2] += flow
        return [q[bus] for bus in unknown]

    x = fsolve(balance, [0.99 * SUPPLY_PRESSURE] * len(unknown), xtol=1.0e-12)
    assert np.max(np.abs(balance(x))) < 1.0e-9, "gas network did not converge"
    return pressures(x)


def build_case():
    voltage, injection = solve_power_flow()

    # GENROU mechanical power on system base: terminal power plus armature loss
    withdrawal = {}
    for bus, (_, _, mva, gas_bus) in MACHINES.items():
        current = abs(np.conj(injection[bus] / voltage[bus])) * S_SYS / (mva * 1.0e6)
        pmech = injection[bus].real + GENROU["Ra"] * current**2 * mva * 1.0e6 / S_SYS
        withdrawal[gas_bus] = KF * pmech * S_SYS / 1.0e6
    pressure = solve_gas_network(withdrawal)

    buses = []
    for bus in sorted(voltage):
        buses.append(
            {
                "number": bus,
                "class": "BusInfinite" if bus == SLACK_BUS else "Bus",
                "name": f"Bus {bus}",
                "init": {"Vr": voltage[bus].real, "Vi": voltage[bus].imag},
                "params": {"kv": 230.0},
                "mon": ["Vm"],
            }
        )
    for bus in GAS_BUSES:
        buses.append(
            {
                "number": bus,
                "class": "GasBusInfinite" if bus == SUPPLY_BUS else "GasBus",
                "name": f"Gas {bus}",
                "init": {"p": pressure[bus]},
                "mon": ["p"],
            }
        )

    devices = []
    for name, b1, b2, r, x, b in BRANCHES:
        devices.append({"class": "Branch", "id": name, "ports": {"bus1": b1, "bus2": b2}, "params": {"R": r, "X": x, "G": 0.0, "B": b}})
    for bus, (p, q) in LOADS.items():
        admittance = complex(p, -q) / abs(voltage[bus]) ** 2
        impedance = 1.0 / admittance
        devices.append({"class": "LoadZ", "id": f"L{bus}", "ports": {"bus": bus}, "params": {"R": impedance.real, "X": impedance.imag}})
    for bus, (name, _, mva, gas_bus) in MACHINES.items():
        params = {"p0": injection[bus].real, "q0": injection[bus].imag, **GENROU, "mva": mva, "kf": KF}
        devices.append({"class": "Genrou", "id": name, "ports": {"bus": bus, "GasInputPm": gas_bus}, "params": params, "mon": ["p", "speed", "delta"]})
    for name, b1, b2, length, diameter in GAS_BRANCHES:
        params = {"L": length, "D": diameter, "lambda": FRICTION, "a": SOUND_SPEED, "b1": B1, "b2": B2, "dxmax": DX_MAX}
        devices.append({"class": "GasBranch", "id": name, "ports": {"bus1": b1, "bus2": b2}, "params": params, "mon": ["q1", "q2", "linepack"]})

    return {
        "header": {
            "format_version": 0.2,
            "format_revision": 0,
            "case_name": "SixBusGas",
            "case_description": "Six-bus system whose two GENROU machines burn gas from a four-bus gas network",
            "case_comments": "Steady state written by examples/PhasorDynamics/DynamicSimulation/Toy/SixBusGas/steady_state.py",
        },
        "params": {"freq_base": F_SYS, "va_base": S_SYS},
        "buses": buses,
        "devices": devices,
    }


if __name__ == "__main__":
    default = pathlib.Path(__file__).resolve().parents[5] / "cases" / "PhasorDynamics" / "Toy" / "SixBusGas.case.json"
    path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else default
    path.write_text(json.dumps(build_case(), indent=2) + "\n")
    print(f"Wrote {path}")
