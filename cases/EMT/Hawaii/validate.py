"""Run the Hawaii fault study and check its physical contracts (standard library)."""

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def conversion_checks(case, report):
    counts = collections.Counter(d['class'] for d in case['devices'])
    expected = {'Machine': 30, 'Tgov1': 30, 'Ieeet1': 30, 'Ieeest': 14,
                'LineLumped': 77, 'Transformer': 12, 'LoadZ': 28, 'Switch': 1,
                'PLL': 9, 'OuterPowerControl': 9, 'InnerCurrentControl': 9,
                'PWM': 9, 'Converter': 9, 'DCLink': 9, 'DependentVoltageSource': 9}
    for kind, count in expected.items():
        require(counts[kind] == count, f'{kind} count: {counts[kind]} != {count}')
    require(not counts['Regfma'] and not counts['REGFMA'], 'Unexpected REGFMA replacement')
    regularized = 0
    largest_time_error = 0.0
    for machine in report['machines'].values():
        for axis in ('d', 'q'):
            data = machine[axis]
            require(abs(data['reconstructed_Xp'] - data['source_Xp']) < 1e-12, 'Transient reactance conversion')
            require(abs(data['reconstructed_Xpp'] - data['effective_Xpp']) < 1e-12, 'Subtransient reactance conversion')
            require(all(t > 0 for t in data['winding_open_times_s']), 'Unstable winding pole')
            largest_time_error = max(largest_time_error, *(abs(a / b - 1) for a, b in
                                     zip(data['winding_open_times_s'], data['source_open_times_s'])))
            regularized += data['effective_Xpp'] != data['source_Xpp']
    require(regularized == 18, 'Damper regularization count')
    require(report['power_flow']['max_KCL_pu'] < 1e-10, 'Initial network current balance')
    return {'regularized_machines': regularized,
            'largest_open_time_relative_change': largest_time_error,
            'initial_max_KCL_pu': report['power_flow']['max_KCL_pu']}


def analyze(csv_path, step_path, study, case, state, report, record, output):
    metrics = {'conversion': conversion_checks(case, report), 'study': study,
               'source_revision': report['source_revision'], 'cpu_time_s': record['cpu_time_s'],
               'study_sha256': record['study_sha256']}
    metrics['input_sha256'] = record['input_sha256']
    metrics['executable_sha256'] = record['executable_sha256']
    metrics['library_sha256'] = record['library_sha256']
    buses = [d for d in case['devices'] if d['class'] == 'Bus' and d['id'] != 'fault_bus']
    machines = list(report['machines'])
    inverters = list(report['inverters'])
    vb = {d['inputs']['bus']: d['params']['V'] for d in case['devices'] if d['class'] == 'Machine'}
    # Buses without machines still have source voltage bases in their initial phasors.
    for bus in buses:
        name = bus['id']
        if name not in vb:
            values = state['buses'][name]
            number = name.removeprefix('bus_')
            vm = math.hypot(*report['power_flow']['bus_phasors_pu'][number])
            vb[name] = math.sqrt(sum(values['v' + p]**2 for p in 'abc')) / vm
    channels = ([f'vmag:{b["id"].removeprefix("bus_")}' for b in buses]
                + [f'{kind}:{m.removesuffix("_genrou").replace("_", " ")}'
                   for kind in ('omega', 'p', 'q') for m in machines])
    sample_period = study['dt_monitor']
    window_count = round(1 / (60 * sample_period))
    require(window_count > 0 and math.isclose(window_count * sample_period, 1 / 60, rel_tol=1e-10),
            'Monitor interval must divide one fundamental period')
    decimation = max(1, round(1 / (240 * sample_period)))
    window = collections.deque()
    sums = [0j] * len(channels)
    averaged = []
    first = None
    last_time = -1.0
    switch_changes = []
    switch_state = None
    imax_ratio = 0.0
    dc_range = [math.inf, -math.inf]
    omega_range = [math.inf, -math.inf]
    pq_initial_error = 0.0
    v1_min = math.inf
    samples_seen = 0
    with csv_path.open(newline='') as stream:
        for raw in csv.DictReader(stream):
            row = {key: float(value) for key, value in raw.items()}
            require(all(math.isfinite(v) for v in row.values()), 'Nonfinite monitor value')
            time = row['t']
            require(time >= last_time - 1e-12, 'Monitor time runs backwards')
            opened = row['Switch_fault_switch_open'] > 0.5
            if switch_state is not None and opened != switch_state:
                switch_changes.append([time, opened])
            switch_state = opened
            if first is None:
                first = row
            for m in machines:
                speed = row[f'Machine_{m}_omega'] - 1
                omega_range = [min(omega_range[0], speed), max(omega_range[1], speed)]
            for plant in inverters:
                data = report['inverters'][plant]
                limit = math.hypot(row[f'InnerCurrentControl_{plant}_inner_ilimd'],
                                   row[f'InnerCurrentControl_{plant}_inner_ilimq']) / data['Imax_A']
                imax_ratio = max(imax_ratio, limit)
                dc = row[f'DCLink_{plant}_dc_vdc'] / data['initial_vdc_V']
                dc_range = [min(dc_range[0], dc), max(dc_range[1], dc)]
            instantaneous_v1 = math.sqrt(sum(row[f'Bus_bus_1_v{p}']**2 for p in 'abc')) / vb['bus_1']
            if 1.02 < time < 1.09:
                v1_min = min(v1_min, instantaneous_v1)
            # Duplicate event rows are checked above but do not add weight to cycle averages.
            if time <= last_time + 1e-12:
                continue
            if samples_seen:
                require(abs(time - last_time - sample_period) < 1e-9, 'Nonuniform monitor samples')
            last_time = time
            samples_seen += 1
            rotation = complex(math.cos(2 * math.pi * 60 * time), -math.sin(2 * math.pi * 60 * time))
            values = []
            for bus in buses:
                name = bus['id']
                a, b, c = (row[f'Bus_{name}_v{p}'] for p in 'abc')
                space_vector = complex(math.sqrt(2 / 3) * (a - (b + c) / 2), (b - c) / math.sqrt(2))
                values.append(space_vector * rotation / vb[name])
            for kind in ('omega', 'p', 'q'):
                values.extend(row[f'Machine_{m}_{kind}'] - 1 if kind == 'omega'
                              else row[f'Machine_{m}_{kind}'] / 1e8 for m in machines)
            window.append((time, values))
            for k, value in enumerate(values):
                sums[k] += value
            if len(window) > window_count:
                _, old = window.popleft()
                for k, value in enumerate(old):
                    sums[k] -= value
            if len(window) == window_count and samples_seen % decimation == 0:
                average = [abs(value / window_count) if k < len(buses) else (value / window_count).real
                           for k, value in enumerate(sums)]
                averaged.append([(window[0][0] + time) / 2] + average)
    require(first is not None, 'Empty monitor output')
    require(abs(last_time - study['tmax']) < 1e-9, 'Study stopped before its final time')
    require(len(switch_changes) == 2, f'Expected two fault events, got {switch_changes}')
    for actual, expected in zip(switch_changes, ((1.0, False), (1.1, True))):
        require(abs(actual[0] - expected[0]) < 1e-10 and actual[1] == expected[1], 'Fault event timing')
    for m, data in report['machines'].items():
        for kind, unit in (('p', 'W'), ('q', 'var')):
            error = abs(first[f'Machine_{m}_{kind}'] - data[f'dispatch_{unit}']) / data['rating_VA']
            pq_initial_error = max(pq_initial_error, error)
    for plant, data in report['inverters'].items():
        vd, vq = (first[f'Park_{plant}_voltage_y{n}'] for n in (1, 2))
        id_, iq = (first[f'Park_{plant}_current_y{n}'] for n in (1, 2))
        power = {'p': vd * id_ + vq * iq, 'q': vq * id_ - vd * iq}
        for kind, unit in (('p', 'W'), ('q', 'var')):
            error = abs(power[kind] - data[f'dispatch_{unit}']) / data['rating_VA']
            pq_initial_error = max(pq_initial_error, error)
    require(pq_initial_error < 1e-10, f'Initial dispatch error: {pq_initial_error}')
    # The full baseline run has a 5.543e-5 relative algebraic interpolation
    # error at 2.000556 s. Pin a 1e-4 allowance for monitored limiter values.
    current_limit_tolerance = 1e-4
    require(imax_ratio <= 1 + current_limit_tolerance, f'Current limiter exceeded: {imax_ratio}')
    require(0.8 < dc_range[0] and dc_range[1] < 1.2, f'DC voltage range: {dc_range}')
    require(max(map(abs, omega_range)) < 0.05, f'Machine speed range: {omega_range}')
    require(v1_min < 0.7, f'Fault did not depress bus 1: {v1_min}')
    final_voltage = averaged[-1][1:1 + len(buses)]
    require(min(final_voltage) > 0.8 and max(final_voltage) < 1.2, f'Final voltage range: {min(final_voltage)}, {max(final_voltage)}')
    output.mkdir(parents=True, exist_ok=True)
    with (output / 'Hawaii.averaged.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time'] + channels)
        writer.writerows(averaged)
    steps, orders = [], collections.Counter()
    with step_path.open(newline='') as stream:
        for row in csv.DictReader(stream):
            steps.append(float(row['step']))
            orders[row['order']] += 1
    require(steps and all(math.isfinite(step) and step > 0 for step in steps), 'Invalid accepted solver steps')
    metrics.update({
        'samples': samples_seen, 'final_time_s': last_time, 'events': switch_changes,
        'initial_dispatch_max_error_pu_plant_base': pq_initial_error,
        'maximum_limited_current_ratio': imax_ratio,
        'limited_current_ratio_tolerance': current_limit_tolerance, 'dc_voltage_ratio_range': dc_range,
        'machine_speed_deviation_pu_range': omega_range, 'fault_bus_1_minimum_voltage_pu': v1_min,
        'final_cycle_voltage_pu_range': [min(final_voltage), max(final_voltage)],
        'accepted_steps': {'count': len(steps), 'minimum_s': min(steps), 'median_s': statistics.median(steps),
                           'maximum_s': max(steps), 'orders': dict(orders)},
        'averaging': {'cycles': 1, 'samples_per_cycle': window_count,
                      'description': '60 Hz demodulated voltage vector and P/Q/speed cycle means, labelled at window centre'},
    })
    metrics['averaged_sha256'] = hashlib.sha256((output / 'Hawaii.averaged.csv').read_bytes()).hexdigest()
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--reuse', type=Path, help='Completed run directory with immutable inputs and run.json')
    parser.add_argument('--tmax', type=float, default=1.5, help='CTest covers inception, clearing, and recovery; full study is 10 s')
    parser.add_argument('--rel-tol', type=float)
    parser.add_argument('--abs-tol', type=float)
    args = parser.parse_args()
    temporary = None
    if args.reuse:
        run = args.reuse.resolve()
        record = json.loads((run / 'run.json').read_text())
        for name, expected in record['input_sha256'].items():
            require(hashlib.sha256((run / name).read_bytes()).hexdigest() == expected,
                    f'Run input changed since simulation: {name}')
        require(hashlib.sha256((run / 'study.json').read_bytes()).hexdigest() == record['study_sha256'],
                'Study settings changed since simulation')
        study = json.loads((run / 'study.json').read_text())
    else:
        require(args.exe is not None, '--exe is required unless --reuse is supplied')
        require(math.isfinite(args.tmax) and args.tmax >= 1.2, 'Study must include fault recovery')
        if args.output:
            run = args.output.resolve()
            run.mkdir(parents=True, exist_ok=True)
        else:
            temporary = tempfile.TemporaryDirectory(prefix='gridkit-hawaii-')
            run = Path(temporary.name)
        require(run != ROOT, 'Use a separate run directory')
        record = {'input_sha256': {}}
        for name in ('Hawaii.case.json', 'Hawaii.state.json', 'conversion.json'):
            raw = (ROOT / name).read_bytes()
            (run / name).write_bytes(raw)
            record['input_sha256'][name] = hashlib.sha256(raw).hexdigest()
        study = json.loads((ROOT / 'Hawaii.solver.json').read_text())
        study.update(system_model_file='Hawaii.case.json', state_file='Hawaii.state.json', tmax=args.tmax)
        for key, value in [('rel_tol', args.rel_tol), ('abs_tol', args.abs_tol)]:
            if value is not None:
                require(math.isfinite(value) and value > 0, 'Tolerances must be finite and positive')
                study[key] = value
        (run / 'study.json').write_text(json.dumps(study, indent=2) + '\n')
        record['study_sha256'] = hashlib.sha256((run / 'study.json').read_bytes()).hexdigest()
        exe = args.exe.resolve()
        record['executable_sha256'] = hashlib.sha256(exe.read_bytes()).hexdigest()
        build = exe.parents[2]
        libraries = (build / 'GridKit').rglob('libgridkit*')
        record['library_sha256'] = {str(name.relative_to(build)): hashlib.sha256(name.read_bytes()).hexdigest()
                                   for name in libraries if name.is_file() and name.suffix in ('.so', '.dylib', '.dll')}
        (run / 'run.json').unlink(missing_ok=True)
        with (run / 'simulation.log').open('w') as log:
            subprocess.run([str(exe), str(run / 'study.json')], cwd=run, stdout=log, stderr=subprocess.STDOUT, check=True)
        record['cpu_time_s'] = float(re.search(r'Complete in (\S+) seconds', (run / 'simulation.log').read_text()).group(1))
        (run / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
    output = args.output.resolve() if args.output else run
    metrics = analyze(run / study['output_file'], run / study['step_output_file'], study,
                      json.loads((run / 'Hawaii.case.json').read_text()),
                      json.loads((run / 'Hawaii.state.json').read_text()),
                      json.loads((run / 'conversion.json').read_text()), record, output)
    (output / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
    print(json.dumps({key: value for key, value in metrics.items() if key not in ('study', 'library_sha256')}, indent=2))
    if temporary:
        temporary.cleanup()


if __name__ == '__main__':
    main()
