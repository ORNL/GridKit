#!/usr/bin/env python3
"""Predict continuous PWM harmonics from an independent periodic sigmoid sum.

Every pulse replica uses the instantaneous sinusoidal duty. Fourier quadrature
resolves the sigmoid edges independently of the simulation's monitor spacing.
"""

import argparse
import json
import math

import numpy as np


def pwm_peak_coefficients(mu, fc=900., f=60., M=.8, vdc=1.,
                          alignment=.5, max_harmonic=49):
    """Return phase-a peak amplitudes for gate s and converter voltage e.

    Every harmonic is included, even zero or tiny coefficients. Peak amplitudes
    are twice the magnitude of the complex Fourier-series coefficient. DC is
    excluded: the gate mean is 1/2 and the converter mean is zero. Converter
    voltage is Vdc * (sa - (sa + sb + sc)/3), which removes triplen harmonics.
    """
    if not all(math.isfinite(x) for x in (mu, fc, f, M, vdc, alignment)):
        raise ValueError('Require finite PWM parameters')
    if not (mu > 0 and fc > f > 0 and 0 <= M <= 1 and 0 <= alignment <= 1 and vdc >= 0):
        raise ValueError('Invalid PWM frequency, modulation, alignment, or DC voltage')
    if not isinstance(max_harmonic, int) or max_harmonic < 1:
        raise ValueError('Require a positive integer maximum harmonic')
    count = round(fc / f)
    if count % 3 or not math.isclose(fc / f, count, rel_tol=1e-13):
        raise ValueError('Require fc / f to be a positive multiple of three')
    points = max(8192, math.ceil(16 * mu / f), 64 * max_harmonic)
    time = np.arange(points) / (points * f)
    duty = (1 + M * np.sin(2 * np.pi * f * time)) / 2
    on = alignment * (1 - duty)
    off = on + duty
    carrier = np.remainder(time * fc, 1)
    radius = math.ceil(math.log(4 / np.finfo(float).eps) * fc / mu)
    gate = np.zeros(points)
    for k in range(-radius, radius + 1):
        gate += .5 * (np.tanh(.5 * mu / fc * (carrier - on - k))
                      - np.tanh(.5 * mu / fc * (carrier - off - k)))
    peaks = 2 * np.abs(np.fft.rfft(gate)) / points
    coefficients = []
    for harmonic in range(1, max_harmonic + 1):
        gate_peak = float(peaks[harmonic])
        coefficients.append({
            'harmonic': harmonic,
            'frequency_hz': harmonic * f,
            'gate_peak': gate_peak,
            'converter_peak_v': 0. if harmonic % 3 == 0 else vdc * gate_peak,
        })
    return coefficients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vdc', type=float, default=1., help='DC link voltage [V]')
    parser.add_argument('--fc', type=float, default=900., help='Carrier frequency [Hz]')
    parser.add_argument('--f', type=float, default=60., help='Modulation frequency [Hz]')
    parser.add_argument('--M', type=float, default=.8, help='Modulation index')
    parser.add_argument('--alignment', type=float, default=.5)
    parser.add_argument('--mu', type=float, nargs='+', help='Default: 4f, 4sqrt(f fc), 4fc')
    args = parser.parse_args()
    mus = args.mu or [4 * args.f, 4 * math.sqrt(args.f * args.fc), 4 * args.fc]
    result = {
        'parameters': {'fc_hz': args.fc, 'f_hz': args.f, 'M': args.M,
                       'vdc_v': args.vdc, 'alignment': args.alignment},
        'interpretation': (
            'Small mu suppresses switching and approaches instantaneous duty; '
            'large mu resolves switching. Fixed-duty carrier means are preserved. '
            'CommonMath saturation and limiters also depend on mu. '
            'Converter common-mode removal cancels triplen harmonics, including '
            'the 900 Hz carrier for fc/f = 15. The 780 and 1020 Hz sidebands remain. '
            'Measure network power from terminal voltages and currents.'),
        'predictions': [],
    }
    for mu in mus:
        result['predictions'].append({
            'mu': mu,
            'edge_10_90_ms': 2000 * math.log(9) / mu,
            'harmonics': pwm_peak_coefficients(mu, args.fc, args.f, args.M,
                                               args.vdc, args.alignment),
        })
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
