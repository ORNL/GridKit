"""Physical reference for the documented ideally transposed overhead line.

SI units throughout. No rational fits or time-domain solver states are used.
Deri-Semlyen ground return and exact solid-cylinder internal impedance; see
 data/overhead/line.json for assumptions and primary references.
"""

import json
from pathlib import Path

import numpy as np
from scipy.constants import mu_0 as MU, epsilon_0 as EPS
from scipy.integrate import quad
from scipy.special import ive

HERE = Path(__file__).resolve().parent
T = np.array([[1/np.sqrt(3), np.sqrt(2/3), 0],
              [1/np.sqrt(3), -1/np.sqrt(6), 1/np.sqrt(2)],
              [1/np.sqrt(3), -1/np.sqrt(6), -1/np.sqrt(2)]])


def parameters():
    return json.loads((HERE / "data/overhead/line.json").read_text())


def transpose(a):
    """Average the three cyclic conductor positions of an ideal transposition."""
    return sum(np.roll(np.roll(a, k, axis=-1), k, axis=-2) for k in range(3))/3


def constants(s, data=None, carson=False):
    """Return three modal series impedances [ohm/m], admittances [S/m].

    carson=True uses independent adaptive quadrature of Carson's integral for
    the ground correction (slow; intended for checking selected frequencies).
    """
    d = parameters() if data is None else data
    s = np.atleast_1d(s).astype(complex)
    xy = np.asarray(d["positions_m"])
    dx = xy[:, 0, None] - xy[None, :, 0]
    height = xy[:, 1, None] + xy[None, :, 1]
    distance = np.sqrt(dx**2 + (xy[:, 1, None] - xy[None, :, 1])**2)
    np.fill_diagonal(distance, d["radius_m"])
    image = np.sqrt(dx**2 + height**2)
    logarithm = np.log(image/distance)
    potential = logarithm/(2*np.pi*EPS)
    capacitance = transpose(np.linalg.inv(potential))
    y = s[:, None, None]*capacitance + d["shunt_conductance_S_per_m"]*np.eye(3)
    z = np.zeros_like(y)
    nonzero = s != 0
    w = s[nonzero]
    depth = np.sqrt(d["soil_resistivity_ohm_m"]/(w*MU))
    displaced = np.sqrt(dx**2 + (height + 2*depth[:, None, None])**2)
    z[nonzero] = w[:, None, None]*MU/(2*np.pi)*np.log(displaced/distance)
    if carson:
        for index in np.flatnonzero(nonzero):
            for i in range(3):
                for j in range(i, 3):
                    h, x = height[i, j], dx[i, j]
                    a = h*h*s[index]*MU/d["soil_resistivity_ohm_m"]
                    def integrand(u):
                        return np.exp(-u)*np.cos(x*u/h)/(u + np.sqrt(u*u+a))
                    integral = (quad(lambda u: integrand(u).real, 0, np.inf, epsabs=1e-10)[0]
                                + 1j*quad(lambda u: integrand(u).imag, 0, np.inf, epsabs=1e-10)[0])
                    z[index, i, j] = z[index, j, i] = s[index]*MU/(2*np.pi)*(logarithm[i, j]+2*integral)
    radius, sigma = d["radius_m"], d["conductivity_S_per_m"]
    dc = 1/(np.pi*radius**2*sigma)
    internal = np.full(len(s), dc, complex)
    k = np.sqrt(w*MU*d["relative_permeability"]*sigma)
    # Exponentially scaled Bessel functions avoid overflow at high frequency.
    internal[nonzero] = k/(2*np.pi*radius*sigma)*ive(0, k*radius)/ive(1, k*radius)
    z += internal[:, None, None]*np.eye(3)
    z = transpose(z)
    return np.diagonal(T.T @ z @ T, axis1=1, axis2=2).copy(), np.diagonal(T.T @ y @ T, axis1=1, axis2=2).copy()


def modal(s, length, data=None):
    z, y = constants(s, data)
    gamma = np.sqrt(z*y)
    return np.sqrt(y/z), np.exp(-length*gamma)


def operators(s, length):
    yc, h = modal(s, length)
    return (T[None]*yc[:, None, :]) @ T.T, (T[None]*h[:, None, :]) @ T.T


def terminal_admittance(s, length):
    """Unfitted telegrapher two-port, using coth/csch (independent of waves)."""
    z, y = constants(s)
    q = length*np.sqrt(z*y)
    yc = np.sqrt(y/z)
    a = yc/np.tanh(q)
    # csch in exponentials remains finite even when a mode is highly attenuated.
    b = -2*yc*np.exp(-q)/(-np.expm1(-2*q))
    diagonal = (T[None]*a[:, None, :]) @ T.T
    mutual = (T[None]*b[:, None, :]) @ T.T
    return np.concatenate((np.concatenate((diagonal, mutual), axis=2),
                           np.concatenate((mutual, diagonal), axis=2)), axis=1)


def wave_delays(length):
    """High-frequency modal delays of the homogenized transposed construction."""
    d = parameters()
    xy = np.asarray(d["positions_m"])
    dx = xy[:, 0, None] - xy[None, :, 0]
    direct = np.sqrt(dx**2 + (xy[:, 1, None]-xy[None, :, 1])**2)
    np.fill_diagonal(direct, d["radius_m"])
    log = np.log(np.sqrt(dx**2+(xy[:, 1, None]+xy[None, :, 1])**2)/direct)
    inductance = transpose(MU/(2*np.pi)*log)
    capacitance = transpose(np.linalg.inv(log/(2*np.pi*EPS)))
    return length*np.sqrt(np.diag(T.T@inductance@T)*np.diag(T.T@capacitance@T))
