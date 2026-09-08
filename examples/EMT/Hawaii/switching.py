"""Check all nine Hawaii bridges against periodic logistic pulse edges."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent
CASE = ROOT.parents[2] / 'cases/EMT/Hawaii'


def pulse(time, modulation, frequency, mu):
    duty = (1 + modulation) / 2
    phase = (time * frequency) % 1
    start, stop = (1 - duty) / 2, (1 + duty) / 2
    radius = math.ceil(math.log(4 / 2.220446049250313e-16) * frequency / mu)
    return sum(0.5 * (math.tanh(mu / frequency * (phase - start - k) / 2)
                      - math.tanh(mu / frequency * (phase - stop - k) / 2))
               for k in range(-radius, radius + 1))


def amplitude(time, values, frequency):
    rotated = [value * complex(math.cos(2 * math.pi * frequency * t), -math.sin(2 * math.pi * frequency * t))
               for t, value in zip(time, values)]
    integral = sum((time[k + 1] - time[k]) * (rotated[k + 1] + rotated[k]) / 2 for k in range(len(time) - 1))
    return abs(2 * integral / (time[-1] - time[0]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/switching.json')
    args = parser.parse_args()
    case = json.loads((CASE / 'Hawaii.case.json').read_text())
    report = json.loads((CASE / 'conversion.json').read_text())
    study = json.loads((CASE / 'Hawaii.solver.json').read_text())
    input_hashes = {name: hashlib.sha256((CASE / name).read_bytes()).hexdigest()
                    for name in ('Hawaii.case.json', 'Hawaii.state.json', 'Hawaii.solver.json', 'conversion.json')}
    assert report['choices']['carrier_alignment'] == 0.5, 'Pulse oracle assumes centred carriers'
    for device in case['devices']:
        device.pop('mon', None)
        monitors = {'PWM': ['s'], 'Modulation': ['m'], 'Converter': ['e', 'idc'],
                    'Filter': ['i'], 'DCLink': ['vdc']}
        if device['class'] in monitors:
            device['mon'] = monitors[device['class']]
    with tempfile.TemporaryDirectory(prefix='gridkit-hawaii-switching-') as temporary:
        run = Path(temporary)
        (run / 'case.json').write_text(json.dumps(case))
        study.update(system_model_file=str(run / 'case.json'), state_file=str(CASE / 'Hawaii.state.json'),
                     tmax=1 / 60, dt_monitor=1 / 720000, events=[], output_file='switching.csv', step_output_file='')
        (run / 'study.json').write_text(json.dumps(study))
        result = subprocess.run([str(args.exe.resolve()), str(run / 'study.json')], cwd=run,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
        with (run / 'switching.csv').open(newline='') as stream:
            data = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    assert data and all(math.isfinite(value) for row in data for value in row.values()), 'Nonfinite or empty output'
    time = [row['t'] for row in data]
    assert abs(time[0]) < 1e-12 and abs(time[-1] - study['tmax']) < 1e-12, 'Incomplete switching window'
    assert all(b > a for a, b in zip(time, time[1:])), 'Nonmonotone switching samples'
    log = result.stdout.decode()
    fc, mu = report['choices']['carrier_Hz'], study['mu']
    frequencies = [60, fc - 120, fc - 60, fc, fc + 60, fc + 120, 2 * fc - 60, 2 * fc + 60, 3 * fc - 120, 3 * fc + 120]
    metrics = {'mu': mu, 'carrier_Hz': fc, 'edge_10_90_s': 2 * math.log(9) / mu,
               'monitor_interval_s': study['dt_monitor'], 'window_s': [time[0], time[-1]],
               'description': 'One fundamental cycle including startup. Fourier amplitudes compared with the independent periodic logistic edge sum at each monitored duty.',
               'input_sha256': input_hashes,
               'executable_sha256': hashlib.sha256(args.exe.resolve().read_bytes()).hexdigest(),
               'cpu_time_s': float(re.search(r'Complete in (\S+) seconds', log).group(1)),
               'plants': {}}
    for plant, parameters in report['inverters'].items():
        predicted, measured = [], []
        pulse_error, power_error = 0.0, 0.0
        for row in data:
            edges = [pulse(row['t'], row[f'Modulation_{plant}_modulation_m{p}'], fc, mu) for p in 'abc']
            pulse_error = max(pulse_error, *(abs(a - row[f'PWM_{plant}_pwm_s{p}']) for a, p in zip(edges, 'abc')))
            voltage = [row[f'Converter_{plant}_bridge_e{p}'] for p in 'abc']
            current = [row[f'Filter_{plant}_filter_i{p}'] for p in 'abc']
            vdc = row[f'DCLink_{plant}_dc_vdc']
            power_error = max(power_error, abs(sum(v * i for v, i in zip(voltage, current))
                                              - vdc * row[f'Converter_{plant}_bridge_idc']) / parameters['rating_VA'])
            measured.append(voltage[0])
            predicted.append(vdc * (edges[0] - sum(edges) / 3))
        observed = [amplitude(time, measured, f) for f in frequencies]
        expected = [amplitude(time, predicted, f) for f in frequencies]
        error = max(abs(a - b) for a, b in zip(observed, expected))
        assert pulse_error < 5e-13, (plant, pulse_error)
        assert power_error < 1e-13, (plant, power_error)
        assert error < 1e-8, (plant, error)
        metrics['plants'][plant] = {'pulse_max_error': pulse_error, 'bridge_power_identity_error_pu': power_error,
                                    'frequencies_Hz': frequencies, 'measured_peak_V': observed,
                                    'predicted_peak_V': expected, 'maximum_harmonic_error_V': error}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(metrics, indent=2) + '\n')
    print('All nine switching bridges pass pulse-edge, harmonic, and power-identity checks.')


if __name__ == '__main__':
    main()
