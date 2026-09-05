#!/usr/bin/env python3
"""Predict PWM and converter harmonics by integrating the ideal pulse edges.

Pwm.output replaces each ideal pulse edge with sigmoid(mu * time). Equivalently,
the ideal periodic pulse train is convolved with the logistic density
mu / (4 cosh(mu * time / 2)**2). Its Fourier transform is x / sinh(x), where
x = 2 pi**2 frequency / mu. This changes harmonic amplitudes without phase lag.
The calculation below does not evaluate GridKit's sampled switching waveform.
"""

import argparse
import cmath
import json
import math


def attenuation(frequency, mu):
    """Exact amplitude ratio of a smoothed edge train to ideal switching."""
    if not math.isfinite(mu) or mu <= 0 or not math.isfinite(frequency):
        raise ValueError('Require positive finite mu and finite frequency')
    x = 2 * math.pi**2 * abs(frequency) / mu
    if x == 0:
        return 1.
    # Avoid overflow at heavily attenuated harmonics and cancellation at x=0.
    return 2 * x * math.exp(-x) / -math.expm1(-2 * x)


def pwm_peak_coefficients(mu, fc=900., f=60., M=.8, vdc=1.,
                          alignment=.5, max_harmonic=49):
    """Return phase-a peak amplitudes for gate s and converter voltage vo.

    Every harmonic is included, even zero or tiny coefficients. Peak amplitudes
    are twice the magnitude of the complex Fourier-series coefficient. DC is
    excluded: the gate mean is 1/2 and the converter mean is zero. Converter
    voltage is Vdc * (sa - (sa + sb + sc)/3), which removes triplen harmonics.
    """
    if not all(math.isfinite(x) for x in (fc, f, M, vdc, alignment)):
        raise ValueError('Require finite PWM parameters')
    if not (fc > f > 0 and 0 <= M <= 1 and 0 <= alignment <= 1 and vdc >= 0):
        raise ValueError('Invalid PWM frequency, modulation, alignment, or DC voltage')
    count = round(fc / f)
    if count % 3 or not math.isclose(fc / f, count, rel_tol=1e-13):
        raise ValueError('Require fc / f to be a positive multiple of three')
    period = 1 / f
    edges = []
    for k in range(count):
        duty = (1 + M * math.sin(2 * math.pi * (k + alignment) / count)) / 2
        on = (k + alignment * (1 - duty)) / fc
        off = (k + alignment + (1 - alignment) * duty) / fc
        edges.append((on, off))
    coefficients = []
    for harmonic in range(1, max_harmonic + 1):
        frequency = harmonic * f
        omega = 2 * math.pi * frequency
        coefficient = sum((cmath.exp(-1j * omega * on)
                           - cmath.exp(-1j * omega * off)) / (1j * omega)
                          for on, off in edges) / period
        ideal_gate_peak = 2 * abs(coefficient)
        factor = attenuation(frequency, mu)
        gate_peak = ideal_gate_peak * factor
        coefficients.append({
            'harmonic': harmonic,
            'frequency_hz': frequency,
            'attenuation': factor,
            'ideal_gate_peak': ideal_gate_peak,
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
            'Changing mu changes the applied PWM voltage and the CommonMath '
            'saturation and limiter equations, independently of solver tolerances. '
            'Converter common-mode removal cancels triplen harmonics, including '
            'the 900 Hz carrier for fc/f = 15. The 780 and 1020 Hz sidebands remain. '
            'Network real or reactive power cannot be inferred from voltage '
            'attenuation alone; use simulated terminal voltages and currents.'),
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
