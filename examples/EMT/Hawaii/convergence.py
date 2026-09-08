"""Compare validated Hawaii runs at two solver tolerances (standard library)."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_run(directory):
    metrics = json.loads((directory / 'metrics.json').read_text())
    averaged = directory / 'Hawaii.averaged.csv'
    if hashlib.sha256(averaged.read_bytes()).hexdigest() != metrics['averaged_sha256']:
        raise ValueError('Averaged data do not match the validated run')
    with averaged.open(newline='') as stream:
        data = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    if not data or not all(math.isfinite(value) for row in data for value in row.values()):
        raise ValueError('Empty or nonfinite averaged data')
    return metrics, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--refined', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/convergence.json')
    args = parser.parse_args()
    baseline, a = read_run(args.baseline)
    refined, b = read_run(args.refined)
    for key in ('input_sha256', 'executable_sha256', 'library_sha256'):
        if baseline[key] != refined[key]:
            raise ValueError(f'Runs differ in {key}')
    varying = {'rel_tol', 'abs_tol', 'tmax'}
    settings = [{key: value for key, value in run['study'].items() if key not in varying}
                for run in (baseline, refined)]
    if settings[0] != settings[1]:
        raise ValueError('Only tolerances and duration may differ')
    if any(refined['study'][key] >= baseline['study'][key] for key in ('rel_tol', 'abs_tol')):
        raise ValueError('Refined tolerances must be tighter')
    duration = min(baseline['final_time_s'], refined['final_time_s'])
    if duration < 1.5:
        raise ValueError('Both runs must cover fault inception, clearing, and recovery through 1.5 s')
    end = min(a[-1]['time'], b[-1]['time'])
    a = [row for row in a if row['time'] <= end + 1e-10]
    b = [row for row in b if row['time'] <= end + 1e-10]
    if len(a) != len(b) or a[0].keys() != b[0].keys():
        raise ValueError('Averaged channels or sample counts differ')
    maxima = dict.fromkeys(('vmag', 'omega', 'p', 'q'), 0.0)
    for x, y in zip(a, b):
        if abs(x['time'] - y['time']) > 1e-10:
            raise ValueError('Averaging windows differ')
        for key in x.keys() - {'time'}:
            kind = key.split(':')[0]
            maxima[kind] = max(maxima[kind], abs(x[key] - y[key]))
    result = {
        'description': 'Maximum absolute difference of centred one-cycle channels, including both fault events.',
        'maximum_channel_difference_pu': maxima,
        'duration_s': duration, 'averaged_window_s': [a[0]['time'], a[-1]['time']],
        'baseline_tolerances': [baseline['study'][key] for key in ('rel_tol', 'abs_tol')],
        'refined_tolerances': [refined['study'][key] for key in ('rel_tol', 'abs_tol')],
        'input_sha256': baseline['input_sha256'],
        'baseline_averaged_sha256': baseline['averaged_sha256'],
        'refined_averaged_sha256': refined['averaged_sha256'],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(maxima, indent=2))


if __name__ == '__main__':
    main()
