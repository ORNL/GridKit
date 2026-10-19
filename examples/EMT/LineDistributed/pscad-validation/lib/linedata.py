"""Line constants of an ideally transposed overhead line from its geometry.

GridWorkbench's ``emt.parameters`` package computes the full conductor Z and Y
matrices (tubular skin effect, Deri-Semlyen earth return with earth
permittivity, image capacitance), reduces bundles and eliminates grounded
shield wires. This module averages the cyclic transpositions, decomposes the
balanced matrices into the ground mode and the (degenerate) aerial mode with
the constant orthogonal matrix T, and derives the quantities both simulators
consume: 60 Hz sequence data, the modal characteristic admittance and
propagation function, and the lossless wavefront delays.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from gridworkbench.emt.parameters import (
    SHUNT,
    admittance,
    capacitive,
    geometric,
    impedance,
    incidence,
    leakage,
    load,
    reduce,
)
from gridworkbench.emt.parameters.series import SERIES

#: Constant real modal transformation of a balanced three-phase line: the first column
#: is the ground (common) mode, the other two span the degenerate aerial modes.
T = np.array([[1 / np.sqrt(3), np.sqrt(2 / 3), 0],
              [1 / np.sqrt(3), -1 / np.sqrt(6), 1 / np.sqrt(2)],
              [1 / np.sqrt(3), -1 / np.sqrt(6), -1 / np.sqrt(2)]])
SHUNT_CONDUCTANCE = 1e-11  # S/m, PSCAD's tower default, and the LineDistributed convention


def transposed(matrices):
    """Average the three cyclic conductor positions of an ideal transposition."""
    return sum(np.roll(np.roll(matrices, k, axis=-1), k, axis=-2) for k in range(3)) / 3


def modal(matrices):
    """Ground and aerial eigenvalues of balanced 3 by 3 matrices (..., 3, 3)."""
    diagonal = np.diagonal(T.T @ matrices @ T, axis1=-2, axis2=-1)
    return diagonal[..., 0], diagonal[..., 1]


@dataclass
class LineConstants:
    name: str
    document: dict
    line: object                      # gridworkbench.emt.parameters.Line
    frequency: np.ndarray             # Hz
    z: np.ndarray                     # (frequencies, 2) modal series impedance, ohm/m
    y: np.ndarray                     # (frequencies, 2) modal shunt admittance, S/m
    z_phase: np.ndarray               # (frequencies, 3, 3) transposed phase impedance
    y_phase: np.ndarray               # (frequencies, 3, 3) transposed phase admittance
    capacitance: np.ndarray           # (2,) modal capacitance, F/m
    inductance_hf: np.ndarray         # (2,) lossless (geometric) modal inductance, H/m

    @property
    def length(self):
        return float(self.line.length)

    @property
    def wave_delays(self):
        """Lossless modal wavefront delays in seconds: [ground, aerial]."""
        return self.length * np.sqrt(self.inductance_hf * self.capacitance)

    def characteristic_admittance(self, s=None):
        z, y = self._at(s)
        return np.sqrt(y / z)

    def propagation(self, s=None):
        z, y = self._at(s)
        return np.exp(-self.length * np.sqrt(z * y))

    def attenuation(self, s=None):
        """Propagation with the wavefront delay removed: H exp(s tau)."""
        s_values = 2j * np.pi * self.frequency if s is None else np.atleast_1d(s)
        return self.propagation(s) * np.exp(s_values[:, None] * self.wave_delays[None, :])

    def sequence(self, frequency=60.0):
        """Positive and zero sequence R (ohm/m), L (H/m), C (F/m) at one frequency."""
        omega = 2 * np.pi * frequency
        z, y = self._at(2j * np.pi * np.array([frequency]))
        z, y = z[0], y[0]
        return {
            "frequency_hz": frequency,
            "R1": float(z[1].real), "L1": float(z[1].imag / omega), "C1": float(y[1].imag / omega),
            "R0": float(z[0].real), "L0": float(z[0].imag / omega), "C0": float(y[0].imag / omega),
            "G": SHUNT_CONDUCTANCE,
        }

    def _at(self, s):
        if s is None:
            return self.z, self.y
        omega = np.imag(np.atleast_1d(s))
        return constants(self.line, omega)


def dc_resistance(line):
    """Per-phase DC resistance in ohm/m: the bundle's sub-conductors in parallel."""
    phases = {}
    for k, phase in enumerate(line.phase):
        if phase in "abc":
            phases.setdefault(phase, []).append(k)
    values = []
    for members in phases.values():
        conductance = sum(line.sigma[m] * np.pi * (line.r[m] ** 2 - line.q[m] ** 2) for m in members)
        values.append(1 / conductance)
    return float(np.mean(values))


def phase_matrices(line, omega):
    """Transposed 3 by 3 phase z (ohm/m) and y (S/m) at angular frequencies ``omega``.

    At zero frequency the geometric and earth-return terms vanish and the shield
    wires decouple, so z is the DC resistance and y the shunt conductance.
    """
    omega = np.atleast_1d(np.asarray(omega, float))
    z = np.zeros((len(omega), 3, 3), complex)
    y = np.zeros((len(omega), 3, 3), complex)
    zero = omega == 0
    if zero.any():
        z[zero] = dc_resistance(line) * np.eye(3)
        y[zero] = SHUNT_CONDUCTANCE * np.eye(3)
    if (~zero).any():
        e = incidence(line)
        w = omega[~zero]
        z_full = impedance(w, line, SERIES)
        c_full = capacitive(w, line) / (1j * w[:, None, None])
        zr, cr = reduce(z_full, c_full, e)
        # Ideal transposition as PSCAD's Line Constants Program applies it: average the
        # series impedance and the Maxwell potential coefficients P = C^-1, then invert.
        # Averaging C itself (equivalent to averaging Y) is the other common convention;
        # it reproduces PSCAD's capacitances only to about 0.9% on flat bundled towers
        # and moves the aerial-mode velocity 0.4% below the speed of light, because the
        # product L C = mu0 eps0 of the lossless limit survives averaging P but not C.
        c_t = np.linalg.inv(transposed(np.linalg.inv(cr)))
        z[~zero] = transposed(zr)
        y[~zero] = 1j * w[:, None, None] * c_t + SHUNT_CONDUCTANCE * np.eye(3)
    return z, y


def constants(line, omega):
    """Transposed modal z (ohm/m) and y (S/m) at angular frequencies ``omega``."""
    z, y = phase_matrices(line, omega)
    z0, z1 = modal(z)
    y0, y1 = modal(y)
    return np.column_stack((z0, z1)), np.column_stack((y0, y1))


def compute(path, frequency):
    """Line constants of the geometry document at ``path`` over ``frequency`` (Hz);
    a zero frequency is the analytic DC limit."""
    path = Path(path)
    document = json.loads(path.read_text())
    line = load(path)
    omega = 2 * np.pi * np.asarray(frequency, float)
    e = incidence(line)
    z, y = phase_matrices(line, omega)
    z0, z1 = modal(z)
    y0, y1 = modal(y)
    # Lossless limit: geometric inductance only, image capacitance only.
    one = np.array([1.0])
    l_full = geometric(one, line) / 1j
    c_full = capacitive(one, line) / 1j
    l_hf, c_hf = reduce(l_full, c_full, e)
    l0, l1 = modal(transposed(l_hf.real))
    c0, c1 = modal(transposed(c_hf.real))
    return LineConstants(
        name=path.name.removesuffix(".line.json"),
        document=document,
        line=line,
        frequency=np.asarray(frequency, float),
        z=np.column_stack((z0, z1)),
        y=np.column_stack((y0, y1)),
        z_phase=z,
        y_phase=y,
        capacitance=np.array([c0[0], c1[0]]),
        inductance_hf=np.array([l0[0], l1[0]]),
    )


def tower(data: LineConstants):
    """The PSCAD universal-tower description of the geometry: phase bundle centres,
    sub-conductor count and spacing, shield wires, conductor data and earth."""
    line = data.line
    phases = {}
    for k, phase in enumerate(line.phase):
        phases.setdefault(phase, []).append(k)
    conductors = []
    bundle = None
    spacing = 0.0
    material = None
    for phase in "abc":
        members = phases[phase]
        x = float(np.mean(line.x[members]))
        h = float(np.mean(line.h[members]))
        conductors.append({"phase": phase, "x": x, "h": h})
        if bundle is None:
            bundle = len(members)
            if bundle > 1:
                first = members[0]
                spacing = float(min(np.hypot(line.x[first] - line.x[m], line.h[first] - line.h[m]) for m in members[1:]))
        elif len(members) != bundle:
            raise ValueError(f"{data.name}: PSCAD towers take one bundle size for all phases")
        for m in members:
            current = (float(line.r[m]), float(line.q[m]), float(line.sigma[m]), float(line.mu[m]))
            if material is None:
                material = current
            elif current != material:
                raise ValueError(f"{data.name}: PSCAD towers take one phase conductor type")
    radius, inner, sigma, mu = material
    ground = []
    for m in phases.get("g", []):
        ground.append({
            "x": float(line.x[m]), "h": float(line.h[m]), "radius": float(line.r[m]),
            "dc_resistance_ohm_per_km": 1e3 / (float(line.sigma[m]) * np.pi * (float(line.r[m]) ** 2 - float(line.q[m]) ** 2)),
            "relative_permeability": float(line.mu[m] / (4e-7 * np.pi)),
        })
    if len(ground) > 2:
        raise ValueError(f"{data.name}: PSCAD towers take at most two ground wires")
    earth = data.document["earth"]
    return {
        "conductors": conductors,
        "sub_conductors": bundle,
        "bundle_spacing_m": spacing,
        "conductor": {
            "radius_m": radius, "inner_radius_m": inner,
            "dc_resistance_ohm_per_km": 1e3 / (sigma * np.pi * (radius ** 2 - inner ** 2)),
            "relative_permeability": float(mu / (4e-7 * np.pi)),
        },
        "ground_wires": ground,
        "earth_resistivity_ohm_m": 1 / earth["conductivity"],
        "earth_relative_permittivity": earth.get("permittivity", 1.0),
        "shunt_conductance_S_per_m": SHUNT_CONDUCTANCE,
        "length_m": data.length,
        "ideal_transposition": True,
    }
