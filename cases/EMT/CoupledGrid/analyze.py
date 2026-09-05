#!/usr/bin/env python3
"""Measure the CoupledGrid mu sweep from its manifest and raw monitor CSVs."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from pwm_analysis import pwm_peak_coefficients

HERE = Path(__file__).resolve().parent


def read_csv(path):
    with path.open(newline='') as stream:
        header = next(csv.reader(stream))
    values = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    if values.shape[1] != len(header) or not np.isfinite(values).all():
        raise ValueError(f'Malformed or nonfinite monitor output: {path}')
    if np.any(np.diff(values[:, 0]) < -1e-12):
        raise ValueError(f'Monitor times run backwards: {path}')
    # Keep the right-hand state when both sides of an event are monitored.
    values = values[np.r_[np.diff(values[:, 0]) > 1e-12, True]]
    t = values[:, 0]
    if len(t) < 3:
        raise ValueError(f'Insufficient monitor samples: {path}')
    dt = float(np.median(np.diff(t)))
    if np.max(np.abs(np.diff(t) - dt)) > dt * 1e-5:
        raise ValueError(f'Expected uniform monitor samples: {path}')
    return t, {key: values[:, i] for i, key in enumerate(header)}


def monitor(data, device, name):
    return data[f"{device['class']}_{device['id']}_{name}"]


def phases(data, device, stem):
    return np.column_stack([monitor(data, device, stem + p) for p in 'abc'])


def phasor(t, value, frequency):
    """RMS Fourier coefficient, using a half-open integer-cycle window."""
    return np.sqrt(2.) * np.mean(value * np.exp(-2j * np.pi * frequency * t)[:, None], axis=0)


def phase_metrics(t, value, frequency, voltage_base=None):
    fundamental = phasor(t, value, frequency)
    a = np.exp(2j * np.pi / 3)
    zero, positive, negative = np.array([[1, 1, 1], [1, a, a*a], [1, a*a, a]]) @ fundamental / 3
    reconstruction = np.sqrt(2.) * np.real(np.exp(2j * np.pi * frequency * t)[:, None] * fundamental)
    residual = value - np.mean(value, axis=0) - reconstruction
    result = {
        'phase_rms': np.sqrt(np.mean(value**2, axis=0)).tolist(),
        'fundamental_phase_rms': np.abs(fundamental).tolist(),
        'positive_sequence_rms': float(abs(positive)),
        'negative_sequence_rms': float(abs(negative)),
        'zero_sequence_rms': float(abs(zero)),
        'negative_to_positive_percent': float(100 * abs(negative / positive)) if abs(positive) > 1e-9 else None,
        'residual_without_dc_or_60hz_rms': np.sqrt(np.mean(residual**2, axis=0)).tolist(),
    }
    if voltage_base is not None:
        line_voltage = value - np.roll(value, -1, axis=1)
        line_rms = np.sqrt(np.mean(line_voltage**2, axis=0))
        result['line_to_line_rms_v'] = line_rms.tolist()
        result['mean_square_line_to_line_rms_pu'] = float(np.sqrt(np.mean(line_rms**2)) / voltage_base)
    return result


def powers(t, voltage, current, frequency):
    """Terminal power, with current positive into the network."""
    vf, cf = phasor(t, voltage, frequency), phasor(t, current, frequency)
    fundamental_power = np.sum(vf * np.conj(cf))
    va = (2 * voltage[:, 0] - voltage[:, 1] - voltage[:, 2]) / 3
    vb = (voltage[:, 1] - voltage[:, 2]) / np.sqrt(3)
    ia = (2 * current[:, 0] - current[:, 1] - current[:, 2]) / 3
    ib = (current[:, 1] - current[:, 2]) / np.sqrt(3)
    return {
        'mean_p_mw': float(np.mean(np.sum(voltage * current, axis=1)) / 1e6),
        'fundamental_p_mw': float(fundamental_power.real / 1e6),
        'fundamental_q_mvar': float(fundamental_power.imag / 1e6),
        'mean_clarke_q_mvar': float(np.mean(1.5 * (vb * ia - va * ib)) / 1e6),
        'phase_current_rms_a': np.sqrt(np.mean(current**2, axis=0)).tolist(),
    }


def kcl(t, data, devices):
    """Sum every monitored terminal injection; includes both line shunts."""
    balances = {d['id']: np.zeros((len(t), 3)) for d in devices if d['class'] == 'Bus'}
    for d in devices:
        cls, inputs = d['class'], d.get('inputs', {})
        if cls in ('VoltageSource', 'DependentVoltageSource', 'Machine', 'LoadZ'):
            balances[inputs['bus']] += phases(data, d, 'i')
        elif cls in ('LineLumped', 'Switch'):
            current = phases(data, d, 'i12')
            balances[inputs['bus1']] -= current
            balances[inputs['bus2']] += current
            if cls == 'LineLumped':
                balances[inputs['bus1']] += phases(data, d, 'i_sh1')
                balances[inputs['bus2']] += phases(data, d, 'i_sh2')
    per_bus = {}
    for bus, balance in balances.items():
        row, phase = np.unravel_index(np.argmax(np.abs(balance)), balance.shape)
        per_bus[bus] = {'max_abs_a': float(abs(balance[row, phase])),
                        'rms_a': float(np.sqrt(np.mean(balance**2))),
                        'worst_time_s': float(t[row]), 'worst_phase': 'abc'[phase]}
    return {'max_abs_a': max(v['max_abs_a'] for v in per_bus.values()), 'buses': per_bus}


def window_mask(t, start, end, frequency):
    dt = float(np.median(np.diff(t)))
    if start < t[0] - dt / 4 or end > t[-1] + dt / 4:
        raise ValueError(f'Window [{start}, {end}) is outside the simulation')
    mask = (t >= start - dt / 4) & (t < end - dt / 4)
    expected = round((end - start) / dt)
    if mask.sum() != expected or not np.isclose(expected * dt * frequency, round(expected * dt * frequency), atol=1e-8):
        raise ValueError('Spectral window must contain uniformly sampled complete fundamental cycles')
    return mask


def window_metrics(t, data, devices, frequency, voltage_base, mask):
    time = t[mask]
    buses = {d['id']: phases(data, d, 'v')[mask] for d in devices if d['class'] == 'Bus'}
    result = {'bus_voltage': {bus: phase_metrics(time, v, frequency, voltage_base) for bus, v in buses.items()},
              'machines': {}, 'ibrs': {}, 'voltage_sources': {}, 'loads': {}}
    groups = {'Machine': 'machines', 'DependentVoltageSource': 'ibrs',
              'VoltageSource': 'voltage_sources', 'LoadZ': 'loads'}
    for d in devices:
        if d['class'] not in groups:
            continue
        metrics = powers(time, buses[d['inputs']['bus']], phases(data, d, 'i')[mask], frequency)
        if d['class'] == 'Machine':
            speed = monitor(data, d, 'omega')[mask] * frequency
            field = monitor(data, d, 'efd')[mask]
            metrics.update(frequency_mean_hz=float(np.mean(speed)), frequency_min_hz=float(np.min(speed)),
                           frequency_max_hz=float(np.max(speed)), efd_mean_pu=float(np.mean(field)),
                           efd_min_pu=float(np.min(field)), efd_max_pu=float(np.max(field)),
                           monitored_mean_p_mw=float(np.mean(monitor(data, d, 'p')[mask])) / 1e6,
                           monitored_mean_q_mvar=float(np.mean(monitor(data, d, 'q')[mask])) / 1e6)
        result[groups[d['class']]][d['id']] = metrics
    result['total_ibr_p_mw'] = sum(d['mean_p_mw'] for d in result['ibrs'].values())
    result['total_ibr_fundamental_q_mvar'] = sum(d['fundamental_q_mvar'] for d in result['ibrs'].values())
    return result


def converter_validation(t, data, case, mu, mask):
    signals = {s['id']: s for s in case['signals']}
    devices = case['devices']
    result = {}
    for converter in (d for d in devices if d['class'] == 'Converter'):
        pwm = next(d for d in devices if d['class'] == 'PWM' and d['outputs']['s'] == converter['inputs']['s'])
        vdc = signals[converter['inputs']['vdc']]['value']
        params = pwm['params']
        predicted = pwm_peak_coefficients(mu, fc=params['fc'], f=params['fm'], M=params['M'],
                                          vdc=vdc, alignment=params.get('alignment', .5))
        voltage = phases(data, converter, 'vo')[mask]
        gate = phases(data, pwm, 's')[mask]
        rows = []
        for h in predicted:
            measured_v = float(np.sqrt(2.) * abs(phasor(t[mask], voltage, h['frequency_hz'])[0]))
            measured_s = float(np.sqrt(2.) * abs(phasor(t[mask], gate, h['frequency_hz'])[0]))
            rows.append(dict(harmonic=h['harmonic'], frequency_hz=h['frequency_hz'],
                             predicted_peak_v=h['converter_peak_v'], measured_peak_v=measured_v,
                             voltage_abs_error_v=abs(measured_v - h['converter_peak_v']),
                             predicted_gate_peak=h['gate_peak'], measured_gate_peak=measured_s,
                             gate_abs_error=abs(measured_s - h['gate_peak'])))
        result[converter['id']] = {
            'dc_voltage_v': vdc, 'all_harmonics_1_to_49_max_abs_error_v': max(h['voltage_abs_error_v'] for h in rows),
            'all_harmonics_1_to_49_max_gate_abs_error': max(h['gate_abs_error'] for h in rows),
            'common_mode_sum_max_abs_v': float(np.max(np.abs(np.sum(voltage, axis=1)))),
            'selected_harmonics': [h for h in rows if h['harmonic'] in (1, 13, 15, 17)],
        }
    return result


def compare(t, data, reference_time, reference, devices, frequency, voltage_base, mask):
    """Compare matching raw output samples with explicit physical scales."""
    if t.shape != reference_time.shape or np.max(np.abs(t - reference_time)) > 1e-10:
        raise ValueError('Comparison runs have different monitor times')
    groups = {'bus_voltage': [], 'machine_frequency': [], 'machine_field_voltage': [],
              'machine_active_power': [], 'ibr_current': []}
    for d in devices:
        prefix = f"{d['class']}_{d['id']}_"
        if d['class'] == 'Bus':
            groups['bus_voltage'] += [(prefix + 'v' + p, 1., voltage_base * np.sqrt(2 / 3), 'V') for p in 'abc']
        elif d['class'] == 'Machine':
            groups['machine_frequency'].append((prefix + 'omega', frequency, 1., 'Hz'))
            groups['machine_field_voltage'].append((prefix + 'efd', 1., 1., 'pu'))
            groups['machine_active_power'].append((prefix + 'p', 1e-6, 1., 'MW'))
        elif d['class'] == 'DependentVoltageSource':
            # Reference RMS gives a useful scale even when a phase crosses zero.
            groups['ibr_current'] += [(prefix + 'i' + p, 1., max(1., float(np.sqrt(np.mean(reference[prefix + 'i' + p][mask]**2)))), 'A') for p in 'abc']
    result = {}
    for name, columns in groups.items():
        entries = []
        for column, factor, scale, unit in columns:
            error = (data[column][mask] - reference[column][mask]) * factor
            worst = int(np.argmax(np.abs(error)))
            entries.append({'column': column, 'max_abs': float(abs(error[worst])),
                            'rms': float(np.sqrt(np.mean(error**2))), 'unit': unit,
                            'normalization_scale': scale, 'normalized_max_abs': float(np.max(np.abs(error)) / scale),
                            'time_of_max_s': float(t[mask][worst])})
        result[name] = {'worst_normalized': max(entries, key=lambda x: x['normalized_max_abs']),
                        'worst_absolute': max(entries, key=lambda x: x['max_abs'])}
    return result


def analyze(case, manifest, root):
    devices = case['devices']
    frequency = manifest.get('frequency_hz', 60.)
    voltage_base = manifest.get('voltage_ll_v', 13800.)
    runs = sorted(manifest['runs'], key=lambda r: r['mu'])
    # The high-mu result is a comparison baseline, never an accuracy oracle.
    baseline_time, baseline = read_csv(root / runs[-1]['csv'])
    result = {
        'definitions': {
            'window': 'Half-open 0.1 s windows, six nominal 60 Hz cycles, no FFT taper.',
            'power': 'Mean P uses all instantaneous abc products; positive means injection into the grid. Fundamental Q is imag(sum(V60 * conj(I60))) with RMS phasors.',
            'sequence': 'RMS 60 Hz abc Fourier coefficients; positive sequence has b lagging a by 120 degrees. Voltage units are V; negative/positive ratio is percent.',
            'residual': 'RMS after subtracting DC and the 60 Hz component; includes harmonics and transient drift, so it is not THD.',
            'comparison': 'Primary minus high-mu describes different smoothed models. Primary minus tight uses the same mu and measures sensitivity to integration tolerances.',
            'timing': 'Median of supplied Complete in CPU seconds; wall seconds are separate. IDA counters include every event segment and consistent-initial-condition work.',
        },
        'baseline_mu': runs[-1]['mu'], 'runs': [],
    }
    for run in runs:
        t, data = (baseline_time, baseline) if run is runs[-1] else read_csv(root / run['csv'])
        windows = {f'before_event_{i+1}': (event['time'] - .1, event['time'])
                   for i, event in enumerate(manifest.get('events', [])) if event['time'] >= t[0] + .1}
        windows['final'] = (float(t[-1] - .1), float(t[-1]))
        masks = {name: window_mask(t, start, end, frequency) for name, (start, end) in windows.items()}
        metrics = {'mu': run['mu'], 'samples': len(t), 'final_time_s': float(t[-1]),
                   'monitor_dt_s': float(np.median(np.diff(t))), 'kcl': kcl(t, data, devices),
                   'median_complete_seconds': float(np.median(run['complete_seconds'])),
                   'complete_seconds': run['complete_seconds'], 'stats': run.get('stats', {}),
                   'windows': {}}
        if 'wall_seconds' in run:
            metrics['median_wall_seconds'] = float(np.median(run['wall_seconds']))
            metrics['wall_seconds'] = run['wall_seconds']
        for name, mask in masks.items():
            metrics['windows'][name] = dict(start_s=windows[name][0], end_s=windows[name][1],
                                            **window_metrics(t, data, devices, frequency, voltage_base, mask))
        metrics['converter_validation'] = converter_validation(t, data, case, run['mu'], masks['final'])
        comparison_masks = {'whole_run': np.ones(t.shape, dtype=bool), 'after_0p1s': t >= .1,
                            'final_window': masks['final']}
        metrics['difference_from_high_mu'] = {name: compare(t, data, baseline_time, baseline, devices,
                                                            frequency, voltage_base, mask)
                                              for name, mask in comparison_masks.items()}
        if run.get('tight_csv'):
            tight_time, tight = read_csv(root / run['tight_csv'])
            metrics['tolerance_refinement'] = {name: compare(t, data, tight_time, tight, devices,
                                                            frequency, voltage_base, mask)
                                              for name, mask in comparison_masks.items()}
        result['runs'].append(metrics)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, default=HERE / 'CoupledGrid.case.json')
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = args.results.resolve()
    manifest = json.loads((args.manifest or root / 'summary.json').read_text())
    result = analyze(json.loads(args.case.read_text()), manifest, root)
    output = args.output or root / 'metrics.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(f'Wrote numerical comparisons to {output}')


if __name__ == '__main__':
    main()
