"""Compare Hawaii line cases and PWM smoothing with identical solver settings."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def benchmark(args):
    cases = {'baseline': args.baseline.resolve(), 'coupled': ROOT / 'cases/EMT/Hawaii'}
    study = json.loads((ROOT / 'cases/EMT/Hawaii/Hawaii.solver.json').read_text())
    study.update(tmax=args.tmax, dt_monitor=args.dt_monitor,
                 output_file='Hawaii.csv', step_output_file='')
    study.update(max_order=args.max_order, scaled_abs_tol=args.scaled_abs_tol)
    study['events'] = [event for event in study['events'] if event['time'] <= args.tmax]
    exe = args.exe.resolve()
    libraries = sorted((exe.parents[2] / 'GridKit').rglob('libgridkit*.so'))
    result = {
        'executable_sha256': digest(exe),
        'library_sha256': {str(p.relative_to(exe.parents[2])): digest(p) for p in libraries},
        'study': study, 'runs': [],
    }
    combinations = [(label, mu) for label in cases for mu in args.mu]
    for trial in range(args.trials):
        order = combinations[trial % len(combinations):] + combinations[:trial % len(combinations)]
        for label, mu in order:
            directory = args.output.resolve() / label / f'mu{mu:g}' / f'trial{trial + 1}'
            directory.mkdir(parents=True, exist_ok=True)
            hashes = {}
            for name in ('Hawaii.case.json', 'Hawaii.state.json'):
                data = (cases[label] / name).read_bytes()
                (directory / name).write_bytes(data)
                hashes[name] = hashlib.sha256(data).hexdigest()
            settings = dict(study, mu=mu)
            (directory / 'study.json').write_text(json.dumps(settings, indent=2) + '\n')
            start = time.monotonic()
            with (directory / 'simulation.log').open('w') as log:
                subprocess.run([str(exe), 'study.json'], cwd=directory, stdout=log,
                               stderr=subprocess.STDOUT, check=True, timeout=args.timeout)
            wall = time.monotonic() - start
            log = (directory / 'simulation.log').read_text()
            counters = re.search(r'IDA statistics: (.*)', log).group(1)
            record = {
                'case': label, 'mu': mu, 'trial': trial + 1,
                'input_sha256': hashes,
                'cpu_seconds': float(re.search(r'Complete in (\S+) seconds', log).group(1)),
                'wall_seconds': wall,
                'ida': {key: int(value) for key, value in re.findall(r'(\w+)=(\d+)', counters)},
            }
            result['runs'].append(record)
            (directory / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
            (args.output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
            print(f"{label} mu={mu:g} trial={trial + 1}: {record['cpu_seconds']:.6g} CPU s, "
                  f"{record['ida']['steps']} steps", flush=True)
    result['medians'] = [
        {'case': label, 'mu': mu,
         'cpu_seconds': statistics.median(r['cpu_seconds'] for r in result['runs']
                                          if r['case'] == label and r['mu'] == mu)}
        for label, mu in combinations
    ]
    (args.output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True, help='Saved case/state directory before regeneration')
    parser.add_argument('--exe', type=Path, default=ROOT / 'build/application/EMT/EMTDynamicSimulation')
    parser.add_argument('--output', type=Path, default=HERE / 'results/line-parameters/benchmark')
    parser.add_argument('--mu', type=float, nargs='+', default=[240., 50000.])
    parser.add_argument('--tmax', type=float, default=0.2)
    parser.add_argument('--dt-monitor', type=float, default=0.)
    parser.add_argument('--max-order', type=int, choices=range(1, 6), default=5)
    parser.add_argument('--scaled-abs-tol', action='store_true', help='Both cases must supply the same nominal ratings')
    parser.add_argument('--trials', type=int, default=3)
    parser.add_argument('--timeout', type=float, default=600.)
    args = parser.parse_args()
    if args.trials < 1:
        parser.error('--trials must be positive')
    benchmark(args)


if __name__ == '__main__':
    main()
