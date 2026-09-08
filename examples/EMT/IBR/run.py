#!/usr/bin/env python3
"""Run and check the six EMT examples; no plotting dependencies required."""
import argparse
import csv
import hashlib
import json
import math
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCENARIOS = sorted(p.name.removesuffix('.solver.json') for p in HERE.glob('*.solver.json'))


def validate_csv(path, final_time):
    with path.open(newline='') as stream:
        rows = csv.reader(stream)
        header = next(rows)
        count, previous, duplicates = 0, -1., 0
        for row in rows:
            if len(row) != len(header):
                raise ValueError(f'{path}: incomplete row {count + 1}')
            values = [float(x) for x in row]
            if not all(math.isfinite(x) for x in values) or values[0] < previous:
                raise ValueError(f'{path}: nonfinite values or decreasing time at row {count + 1}')
            duplicates += values[0] == previous
            previous = values[0]
            count += 1
        if not math.isclose(previous, final_time, abs_tol=1e-12):
            raise ValueError(f'{path}: stopped at {previous}, expected {final_time}')
    return dict(rows=count, columns=len(header), duplicate_event_times=duplicates, final_time=previous)


def validate_nominal_voltage(path):
    """The smooth case starts with a balanced 1.02 pu converter voltage."""
    with path.open(newline='') as stream:
        initial = next(csv.DictReader(stream))
    for bus in (4, 5, 6):
        v = [float(initial[f'Converter_converter_{bus}_e{p}']) for p in 'abc']
        vll = math.sqrt(sum((v[p] - v[(p + 1) % 3])**2 for p in range(3)) / 3)
        if not math.isclose(vll, 1.02 * 13800, rel_tol=1e-6):
            raise ValueError(f'{path}: converter {bus} starts at {vll} V line-line RMS; expected 14076 V')


def run_one(name, exe, results, smoke, check_only=False, mu=None):
    source = HERE / (name + '.solver.json')
    study = json.loads(source.read_text())
    directory = results / name
    directory.mkdir(parents=True, exist_ok=True)
    for key in ('system_model_file', 'state_file'):
        study[key] = str((HERE / study[key]).resolve())
    if mu is not None:
        study['mu'] = mu
        edge_step = 10**math.floor(math.log10(2*math.log(9)/mu/8))
        study['dt_monitor'] = min(study['dt_monitor'], 1e-5, edge_step)
    if smoke:
        study['tmax'] *= .01
        study['dt_monitor'] = .001
        for event in study['events']:
            event['time'] *= .01
    study_path = directory / 'run.solver.json'
    if check_only:
        if json.loads(study_path.read_text()) != study:
            raise ValueError(f'{name}: saved run settings differ; rerun the study')
    else:
        study_path.write_text(json.dumps(study, indent=4) + '\n')
    wall = None
    if not check_only:
        start = time.perf_counter()
        with (directory / 'run.log').open('w') as log:
            subprocess.run([str(exe), str(study_path)], cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True)
        wall = time.perf_counter() - start
    monitors = validate_csv(directory / study['output_file'], study['tmax'])
    states = validate_csv(directory / study['state_output_file'], study['tmax'])
    if study.get('mu', 240.) == 240.:
        validate_nominal_voltage(directory / study['output_file'])
    expected_rows = round(study['tmax'] / study['dt_monitor']) + 1
    for event in study['events']:
        sample = event['time'] / study['dt_monitor']
        expected_rows += 1 if math.isclose(sample, round(sample), rel_tol=0, abs_tol=1e-9) else 2
    if monitors['rows'] != expected_rows or states['rows'] != expected_rows:
        raise ValueError(f'{name}: missing samples: expected {expected_rows}')
    if monitors['rows'] != states['rows'] or monitors['final_time'] != states['final_time']:
        raise ValueError(f'{name}: monitor/state sample mismatch')
    schema = json.loads((directory / (study['state_output_file'] + '.json')).read_text())
    variables = schema['variables']
    count = (states['columns'] - 1) // 2
    if sorted(v['index'] for v in variables) != list(range(count)):
        raise ValueError(f'{name}: incomplete DAE index map')
    # Retain only event samples rather than loading the complete monitor file.
    event_samples = [[] for _ in study['events']]
    if event_samples:
        with (directory / study['output_file']).open(newline='') as stream:
            for row in csv.DictReader(stream):
                sample_time = float(row['t'])
                for event, samples in zip(study['events'], event_samples):
                    if abs(sample_time - event['time']) < 1e-12:
                        samples.append(float(row[f"Switch_{event['element_id']}_open"]))
    for event, samples in zip(study['events'], event_samples):
        command = float(event['open'])
        if len(samples) != 2 or samples[-1] != command or samples[0] == command:
            raise ValueError(f'{name}: event command not recorded correctly: {event}')
    hashes = {key: hashlib.sha256(Path(study[key]).read_bytes()).hexdigest() for key in ('system_model_file', 'state_file')}
    if check_only:
        previous = json.loads((directory / 'validation.json').read_text())
        if previous['input_sha256'] != hashes:
            raise ValueError(f'{name}: case/state changed since the recorded run')
        wall = previous['wall_seconds']
    summary = dict(scenario=name, smoke=smoke, executable=str(exe), wall_seconds=wall,
                   input_sha256=hashes, solver_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                   monitors=monitors, states=states, dae_variables=count,
                   differential_variables=sum(v['differential'] for v in variables), events=study['events'])
    (directory / 'validation.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f"{name}: PASS, {monitors['rows']} samples, {count} DAE variables, wall={wall} s", flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True, help='EMTDynamicSimulation executable')
    parser.add_argument('--results', type=Path, default=HERE / 'results')
    parser.add_argument('--scenario', nargs='+', choices=SCENARIOS)
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--mu', type=float, help='Override smoothing; DC voltage stays fixed')
    parser.add_argument('--check-only', action='store_true', help='Validate existing outputs without rerunning')
    parser.add_argument('--smoke', action='store_true', help='Short integration and event checks for CTest')
    args = parser.parse_args()
    if args.mu is not None and (not math.isfinite(args.mu) or args.mu <= 0):
        parser.error('--mu must be positive and finite')
    scenarios = args.scenario or (['04_FaultClearing'] if args.smoke else SCENARIOS)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(lambda name: run_one(name, args.exe.resolve(), args.results.resolve(), args.smoke, args.check_only, args.mu), scenarios))


if __name__ == '__main__':
    main()
