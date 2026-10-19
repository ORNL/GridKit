"""Vector fitting of smooth scalar frequency responses with common poles.

A compact implementation of Gustavsen and Semlyen's vector fitting (pole
relocation with the sigma function, unstable poles flipped) followed by a
weighted residue fit. The residue fit can pin the response exactly at the
lowest frequency sample (the DC-like constraint the LineDistributed fits use).

Only NumPy is required. Poles are returned with each complex pole immediately
followed by its conjugate, the ordering GridKit's VectorFit parser expects.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class Fit:
    poles: np.ndarray      # complex, conjugates adjacent
    residues: np.ndarray   # (functions, poles) complex, conjugates adjacent
    d: np.ndarray          # (functions,) real constant terms

    def evaluate(self, s):
        s = np.atleast_1d(np.asarray(s, complex))
        basis = 1 / (s[:, None] - self.poles[None, :])
        return basis @ self.residues.T + self.d[None, :]


def initial_poles(frequency, count, real_count=2):
    """Log-spaced starting poles: ``real_count`` real poles, the rest complex pairs."""
    pairs = (count - real_count) // 2
    real_count = count - 2 * pairs
    omega = 2 * np.pi * np.geomspace(frequency[0] if frequency[0] > 0 else frequency[1], frequency[-1], pairs + real_count)
    poles = []
    for w in omega[:real_count]:
        poles.append(complex(-w, 0))
    for w in omega[real_count:]:
        poles.extend([complex(-w / 100, w), complex(-w / 100, -w)])
    return np.array(poles, complex)


def _real_basis(s, poles):
    """Real-coefficient basis columns of the pole set: one per real pole, two per pair."""
    columns, kinds = [], []
    q = 0
    while q < len(poles):
        p = poles[q]
        if p.imag == 0:
            columns.append(1 / (s - p))
            kinds.append(("real", q))
            q += 1
        else:
            a, b = 1 / (s - p), 1 / (s - np.conj(p))
            columns.extend([a + b, 1j * (a - b)])
            kinds.append(("pair", q))
            q += 2
    return np.array(columns).T, kinds


def _coefficients_to_residues(x, kinds, count):
    residues = np.zeros(count, complex)
    index = 0
    for kind, q in kinds:
        if kind == "real":
            residues[q] = x[index]
            index += 1
        else:
            residues[q] = x[index] + 1j * x[index + 1]
            residues[q + 1] = np.conj(residues[q])
            index += 2
    return residues


def _stack(matrix):
    return np.vstack((matrix.real, matrix.imag))


def relocate(s, values, poles, weights, iterations=20):
    """Pole relocation: fit sigma*F and sigma with common poles, take sigma's zeros."""
    values = np.atleast_2d(values)
    functions = values.shape[0]
    for _ in range(iterations):
        basis, kinds = _real_basis(s, poles)
        width = basis.shape[1]
        rows = []
        rhs = []
        for f in range(functions):
            w = weights[f][:, None]
            block = np.zeros((len(s), functions * (width + 1) + width), complex)
            block[:, f * (width + 1):(f + 1) * (width + 1) - 1] = basis
            block[:, (f + 1) * (width + 1) - 1] = 1.0
            block[:, functions * (width + 1):] = -basis * values[f][:, None]
            rows.append(_stack(block * w))
            rhs.append(_stack((values[f] * w[:, 0])[:, None])[:, 0])
        a = np.vstack(rows)
        b = np.concatenate(rhs)
        scale = np.linalg.norm(a, axis=0)
        scale[scale == 0] = 1
        x = np.linalg.lstsq(a / scale, b, rcond=None)[0] / scale
        sigma = x[functions * (width + 1):]
        # Zeros of sigma(s) = 1 + sum c/(s - p): eigenvalues of A - b c^T in real form.
        size = len(poles)
        a_matrix = np.zeros((size, size))
        b_vector = np.zeros(size)
        c_vector = np.zeros(size)
        index = 0
        for kind, q in kinds:
            p = poles[q]
            if kind == "real":
                a_matrix[q, q] = p.real
                b_vector[q] = 1
                c_vector[q] = sigma[index]
                index += 1
            else:
                a_matrix[q, q] = a_matrix[q + 1, q + 1] = p.real
                a_matrix[q, q + 1] = p.imag
                a_matrix[q + 1, q] = -p.imag
                b_vector[q] = 2
                c_vector[q] = sigma[index]
                c_vector[q + 1] = sigma[index + 1]
                index += 2
        new = np.linalg.eigvals(a_matrix - np.outer(b_vector, c_vector))
        new = np.where(new.real > 0, -new.real + 1j * new.imag, new)
        poles = _sort_pairs(new)
    return poles


def _sort_pairs(poles):
    """Real poles first, then conjugate pairs with the positive-imaginary member first."""
    real = sorted((p.real for p in poles if abs(p.imag) <= 1e-12 * max(1.0, abs(p.real))), reverse=True)
    pairs = sorted({(round(p.real, 12), round(abs(p.imag), 12)) for p in poles if abs(p.imag) > 1e-12 * max(1.0, abs(p.real))}, key=lambda item: item[1])
    result = [complex(r, 0) for r in real]
    for re, im in pairs:
        result.extend([complex(re, im), complex(re, -im)])
    return np.array(result, complex)


def residues(s, values, poles, weights, pin_first=True):
    """Weighted residue fit with fixed poles; ``pin_first`` makes the fit exact at s[0]."""
    values = np.atleast_2d(values)
    basis, kinds = _real_basis(s, poles)
    result_residues = np.zeros((values.shape[0], len(poles)), complex)
    d = np.zeros(values.shape[0])
    for f in range(values.shape[0]):
        w = weights[f]
        if pin_first:
            shifted = basis - basis[0].real
            target = values[f] - values[f][0]
            a = _stack(shifted[1:] * w[1:, None])
            b = _stack((target[1:] * w[1:])[:, None])[:, 0]
            scale = np.linalg.norm(a, axis=0)
            scale[scale == 0] = 1
            x = np.linalg.lstsq(a / scale, b, rcond=None)[0] / scale
            d[f] = values[f][0].real - basis[0].real @ x
        else:
            a = _stack(np.column_stack((basis, np.ones(len(s)))) * w[:, None])
            b = _stack((values[f] * w)[:, None])[:, 0]
            scale = np.linalg.norm(a, axis=0)
            scale[scale == 0] = 1
            x = np.linalg.lstsq(a / scale, b, rcond=None)[0] / scale
            d[f] = x[-1]
            x = x[:-1]
        result_residues[f] = _coefficients_to_residues(x, kinds, len(poles))
    return Fit(poles, result_residues, d)


def vector_fit(frequency, values, order, *, weights=None, iterations=20, pin_first=True, real_count=2):
    """Fit ``values`` (functions by frequencies) with ``order`` common poles."""
    values = np.atleast_2d(np.asarray(values, complex))
    s = 2j * np.pi * np.asarray(frequency, float)
    if weights is None:
        weights = [1 / np.maximum(np.abs(v), 1e-3 * np.max(np.abs(v))) for v in values]
    poles = initial_poles(np.asarray(frequency, float), order, real_count)
    poles = relocate(s, values, poles, weights, iterations)
    return residues(s, values, poles, weights, pin_first)


def scaled_error(actual, expected, floor):
    return float(np.max(np.abs(actual - expected) / np.maximum(np.abs(expected), floor)))
