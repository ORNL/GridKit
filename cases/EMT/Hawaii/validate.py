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


def solve_pair(matrix, rhs):
    (a, b), (c, d) = matrix
    determinant = a * d - b * c
    return ((d * rhs[0] - b * rhs[1]) / determinant,
            (a * rhs[1] - c * rhs[0]) / determinant)


def decay_poles(matrix):
    (a, b), (c, d) = matrix
    trace, product = -(a + d), a * d - b * c
    require(trace > 0 and product > 0, 'Unstable rotor circuit')
    fast = (trace + math.sqrt(trace * trace - 4 * product)) / 2
    return [-product / fast, -fast]


def winding_checks(params, axis, data):
    keys = ('Lmd', 'Llfd', 'Ll1d', 'Rfd', 'R1d') if axis == 'd' else (
        'Lmq', 'Ll1q', 'Ll2q', 'R1q', 'R2q')
    lm, l1, l2, r1, r2 = (params[key] for key in keys)
    xl, omega = params['Ll'], 2 * math.pi * params['f']
    x, xp, xpp = (data[key] for key in ('source_X', 'source_Xp', 'effective_Xpp'))
    tp, tpp = data['source_time_parameters_s']
    require(all(value > 0 for value in (lm, l1, l2, r1, r2, tp, tpp)), 'Nonpositive winding parameter')
    require(abs(xl + lm - x) < 1e-12, 'Synchronous reactance conversion')
    require(abs(xl + 1 / (1 / lm + 1 / l1) - xp) < 1e-12, 'Transient reactance conversion')
    require(abs(xl + 1 / (1 / lm + 1 / l1 + 1 / l2) - xpp) < 1e-12,
            'Subtransient reactance conversion')

    # GENROU rotor equations at zero stator current, field voltage, and saturation.
    def genrou_matrix(subtransient):
        coupling = (x - xp) * (xp - subtransient) / (xp - xl)**2
        return [[-(1 + coupling) / tp, coupling / tp], [1 / tpp, -1 / tpp]]

    source_poles = decay_poles(genrou_matrix(data['source_Xpp']))
    effective = genrou_matrix(xpp)
    effective_poles = decay_poles(effective)
    rotor = [[lm + l1, lm], [lm, lm + l2]]
    winding_poles = decay_poles(list(zip(solve_pair(rotor, [-omega * r1, 0]),
                                        solve_pair(rotor, [0, -omega * r2]))))
    for label, poles in (('source', source_poles), ('effective', effective_poles), ('winding', winding_poles)):
        require(all(math.isclose(a, b, rel_tol=1e-12) for a, b in
                    zip(poles, data[label + '_open_poles_per_s'])), f'{label} pole report')
    pole_error = max(abs(a / b - 1) for a, b in zip(winding_poles, effective_poles))
    require(pole_error < 1e-12, 'Winding poles differ from coupled GENROU poles')

    # Compare operational reactance from the two independent state equations.
    coupling = (x - xp) * (xp - xpp) / (xp - xl)**2
    forcing = [(x - xp - coupling * (xp - xl)) / tp, (xp - xl) / tpp]
    weights = [(xpp - xl) / (xp - xl), (xp - xpp) / (xp - xl)]
    response_error = 0.0
    for s in (0.01 + 0.1j, 1 + 1j, 10 + 60j, 1 + 377j, 0.1 + 1e5j):
        source = [[s - effective[0][0], -effective[0][1]],
                  [-effective[1][0], s - effective[1][1]]]
        windings = [[s * rotor[0][0] + omega * r1, s * lm],
                    [s * lm, s * rotor[1][1] + omega * r2]]
        source_x = xpp + sum(a * b for a, b in zip(weights, solve_pair(source, forcing)))
        winding_x = xl + lm - s * lm**2 * sum(solve_pair(windings, [1, 1]))
        response_error = max(response_error, abs(winding_x / source_x - 1))
        if axis == 'd':
            source_field = sum(a * b for a, b in zip(weights, solve_pair(source, [1 / tp, 0])))
            winding_field = sum(solve_pair(windings, [omega * r1, 0]))
            response_error = max(response_error, abs(winding_field / source_field - 1))
    require(response_error < 1e-12, 'Winding transfer differs from unsaturated GENROU dynamics')
    return pole_error, response_error, max(abs(a / b - 1) for a, b in zip(effective_poles, source_poles))


def conversion_checks(case, state, report):
    counts = collections.Counter(d['class'] for d in case['devices'])
    expected = {'Machine': 30, 'Tgov1': 30, 'Ieeet1': 30, 'Ieeest': 14,
                'LineLumped': 77, 'Transformer': 12, 'LoadZ': 29, 'Switch': 2,
                'PLL': 9, 'OuterPowerControl': 9, 'InnerCurrentControl': 9,
                'PWM': 9, 'Converter': 9, 'Filter': 9, 'Park': 36}
    for kind, count in expected.items():
        require(counts[kind] == count, f'{kind} count: {counts[kind]} != {count}')
    require(not counts['Regfma'] and not counts['REGFMA'], 'Unexpected REGFMA replacement')
    require(not counts['DependentVoltageSource'], 'Inverter plants must use LCL Filters')
    devices = {d['id']: d for d in case['devices']}
    signals = {s['id']: s.get('value') for s in case['signals']}
    fault, discharge = devices['fault_load'], devices['fault_discharge_load']
    require(devices['fault_switch']['inputs'] == {'bus1': 'bus_1', 'bus2': 'fault_bus'}, 'Fault connection')
    require(devices['fault_discharge_switch']['inputs'] ==
            {'bus1': 'fault_bus', 'bus2': 'fault_discharge_bus'}, 'Discharge connection')
    require(fault['inputs']['bus'] == 'fault_bus' and discharge['inputs']['bus'] == 'fault_discharge_bus',
            'Fault and discharge loads')
    for i in range(3):
        for j in range(3):
            expected_x = 0.01 if i == j else 0.0
            require(fault['params']['R'][i][j] == 0, 'Fault must be purely inductive')
            require(abs(2 * math.pi * 60 * fault['params']['L'][i][j] / (138e3**2 / 1e8) - expected_x) < 1e-12,
                    'Fault reactance must match the phasor case')
            require(abs(discharge['params']['R'][i][j] / (138e3**2 / 1e8) - expected_x) < 1e-12,
                    'Discharge resistance')
    require(state['devices']['fault_switch']['open'] and
            not state['devices']['fault_discharge_switch']['open'], 'Initial fault switch states')
    require(all(value == 0 for value in state['devices']['fault_load'].values()), 'De-energized fault inductance')

    def phasor(values, prefix):
        a, b, c = (values[prefix + p] for p in 'abc')
        return complex(math.sqrt(2 / 3) * (a - (b + c) / 2), (b - c) / math.sqrt(2))

    for line in (d for d in case['devices'] if d['class'] == 'LineLumped'):
        p = line['params']
        require(math.isfinite(p['dx']) and p['dx'] > 0, 'Positive line length')
        for key in ('Rp', 'Lp', 'Gp', 'Cp'):
            matrix = p[key]
            require(len(matrix) == 3 and all(len(row) == 3 for row in matrix), 'Line matrix dimensions')
            require(all(math.isfinite(v) for row in matrix for v in row), 'Finite line matrix')
            diagonal, mutual = matrix[0][:2]
            require(diagonal - mutual >= 0 and diagonal + 2 * mutual >= 0, 'Passive transposed line')
            require(all(math.isclose(matrix[i][j], diagonal if i == j else mutual, rel_tol=1e-12)
                        for i in range(3) for j in range(3)), 'Transposed line symmetry')
        require(p['Lp'][0][1] > 0 and p['Cp'][0][1] < 0, 'Mutual line inductance and capacitance')
        v1, v2 = (phasor(state['buses'][line['inputs'][port]], 'v') for port in ('bus1', 'bus2'))
        current = phasor(state['devices'][line['id']], 'i12')
        z = p['dx'] * complex(p['Rp'][0][0] - p['Rp'][0][1],
                             2 * math.pi * 60 * (p['Lp'][0][0] - p['Lp'][0][1]))
        require(abs(v1 - v2 - z * current) / max(abs(v1), abs(v2)) < 1e-12, 'Initial coupled line KVL')

    for plant, data in report['inverters'].items():
        filt, pll, inner, outer, bridge = (devices[plant + suffix] for suffix in
                                             ('_filter', '_pll', '_inner', '_power', '_bridge'))
        terminal_voltage = devices[plant + '_terminal_voltage']
        voltage, current, grid_current, pwm = (devices[plant + suffix] for suffix in
                                               ('_voltage', '_current', '_grid_current', '_pwm'))
        require(filt['inputs']['e'] == bridge['outputs']['e'], 'Bridge voltage must drive its Filter')
        require(current['inputs']['input'] == filt['outputs']['i'], 'Inner-loop converter-current measurement')
        require(voltage['inputs']['input'] == filt['outputs']['vo'], 'Capacitor-voltage feedback')
        bus = devices[filt['inputs']['bus']]
        require([pll['inputs']['v' + p] for p in 'abc'] == [bus['outputs']['v' + p] for p in 'abc'],
                'PLL must read the terminal Bus voltage')
        require(terminal_voltage['inputs']['input'] == [bus['outputs']['v' + p] for p in 'abc'],
                'Outer-loop terminal Bus voltage measurement')
        require(outer['inputs']['v'] == terminal_voltage['outputs']['out'][:2],
                'Outer-loop terminal voltage feedback')
        require(outer['params']['Pref'] == data['dispatch_W']
                and outer['params']['Qref'] == data['dispatch_var'], 'Terminal power setpoints')
        require(grid_current['inputs']['input'] == filt['outputs']['ig'], 'Grid-current measurement')
        require(outer['inputs']['i'] == grid_current['outputs']['out'][:2], 'Outer-loop grid-current feedback')
        require(inner['inputs']['i'] == current['outputs']['out'][:2], 'Inner-loop converter-current feedback')
        require(inner['inputs']['v'] == voltage['outputs']['out'][:2], 'Inner-loop capacitor-voltage feedback')
        require(all(p['inputs']['theta'] == pll['outputs']['theta'] for p in (voltage, current, grid_current, terminal_voltage, pwm)),
                'Park transforms and PWM must share the PLL angle')
        require(pwm['inputs']['u'] == inner['outputs']['u'], 'PWM voltage command')
        require(bridge['inputs']['s'] == pwm['outputs']['s'], 'Bridge switching input')
        require(inner['inputs']['omega'] == pll['outputs']['omega'], 'Inner loop must use PLL frequency')
        require(outer['inputs']['ilim'] == inner['outputs']['ilim']
                and inner['inputs']['icmd'] == outer['outputs']['icmd'], 'Outer-loop anti-windup connection')
        require(pwm['inputs']['vdc'] == bridge['inputs']['vdc'], 'Shared DC voltage')
        require(signals[bridge['inputs']['vdc']] == data['vdc_V'], 'Constant DC voltage')
        require(inner['inputs']['ulim'] == pwm['outputs']['ulim'], 'Limited voltage feedback')
        initial = state['devices'][filt['id']]
        v = phasor(state['buses'][filt['inputs']['bus']], 'v')
        ig = phasor(initial, 'ig')
        power = v * ig.conjugate()
        require(abs(power - complex(data['dispatch_W'], data['dispatch_var'])) / data['rating_VA'] < 1e-12,
                'Initial terminal dispatch')
        require(set(initial) == {'iga', 'igb', 'igc'}, 'Filter state must prescribe only grid current')
        require(inner['id'] not in state['devices'] and outer['id'] not in state['devices'],
                'Controller commands must be derived by initialization')
    for name, change in report['exciter_adjustments'].items():
        require(change['Ke_original'] < 0 and devices[name]['params']['Ke'] == 0,
                'Hawaii exciter automatic Ke initialization')
    regularized = 0
    largest_pole_error = largest_response_error = largest_pole_change = 0.0
    for name, machine in report['machines'].items():
        for axis in ('d', 'q'):
            data = machine[axis]
            require(abs(data['reconstructed_Xp'] - data['source_Xp']) < 1e-12, 'Transient reactance conversion')
            require(abs(data['reconstructed_Xpp'] - data['effective_Xpp']) < 1e-12, 'Subtransient reactance conversion')
            pole_error, response_error, pole_change = winding_checks(devices[name]['params'], axis, data)
            largest_pole_error = max(largest_pole_error, pole_error)
            largest_response_error = max(largest_response_error, response_error)
            largest_pole_change = max(largest_pole_change, pole_change)
            regularized += data['effective_Xpp'] != data['source_Xpp']
    require(regularized == 18, 'Damper regularization count')
    require(report['power_flow']['max_KCL_pu'] < 1e-10, 'Initial network current balance')
    return {'LCL_plants': counts['Filter'], 'regularized_machines': regularized,
            'largest_open_pole_relative_error': largest_pole_error,
            'largest_rotor_transfer_relative_error': largest_response_error,
            'largest_regularized_open_pole_relative_change': largest_pole_change,
            'initial_max_KCL_pu': report['power_flow']['max_KCL_pu']}


def analyze(csv_path, step_path, study, case, state, report, record, output):
    metrics = {'conversion': conversion_checks(case, state, report), 'study': study,
               'source_revision': report['source_revision'], 'exciter_adjustments': report['exciter_adjustments'], 'cpu_time_s': record['cpu_time_s'],
               'study_sha256': record['study_sha256']}
    metrics['input_sha256'] = record['input_sha256']
    metrics['executable_sha256'] = record['executable_sha256']
    metrics['library_sha256'] = record['library_sha256']
    buses = [d for d in case['devices'] if d['class'] == 'Bus' and d['id'].startswith('bus_')]
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
    omega_range = [math.inf, -math.inf]
    pq_initial_error = 0.0
    v1_min = math.inf
    samples_seen = 0
    fault_current_base = 1e8 / 138e3
    fault_params = next(d['params'] for d in case['devices'] if d['id'] == 'fault_load')
    discharge_params = next(d['params'] for d in case['devices'] if d['id'] == 'fault_discharge_load')
    decay_rate = discharge_params['R'][0][0] / fault_params['L'][0][0]
    previous_fault = None
    clearing_current = None
    current_jump = open_current = decay_error = 0.0
    with csv_path.open(newline='') as stream:
        for raw in csv.DictReader(stream):
            row = {key: float(value) for key, value in raw.items()}
            require(all(math.isfinite(v) for v in row.values()), 'Nonfinite monitor value')
            time = row['t']
            require(time >= last_time - 1e-12, 'Monitor time runs backwards')
            opened = row['Switch_fault_switch_open'] > 0.5
            require(opened != (row['Switch_fault_discharge_switch_open'] > 0.5), 'Complementary fault switches')
            fault_current = [row[f'LoadZ_fault_load_i{p}'] for p in 'abc']
            if switch_state is not None and opened != switch_state:
                switch_changes.append([time, opened])
                require(abs(time - previous_fault[0]) < 1e-10, 'Fault event must record both sides')
                current_jump = max(current_jump, *(abs(a - b) / fault_current_base
                                                   for a, b in zip(fault_current, previous_fault[1])))
                if opened:
                    clearing_current = (time, fault_current)
            switch_state = opened
            previous_fault = (time, fault_current)
            if opened:
                open_current = max(open_current, *(abs(row[f'Switch_fault_switch_i12{p}']) / fault_current_base
                                                   for p in 'abc'))
            if clearing_current is not None:
                factor = math.exp(-decay_rate * (time - clearing_current[0]))
                scale = max(fault_current_base, *map(abs, clearing_current[1]))
                decay_error = max(decay_error, *(abs(a - factor * b) / scale
                                                for a, b in zip(fault_current, clearing_current[1])))
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
            if len(window) == window_count:
                centre = (window[0][0] + time) / 2
                fault_window = report['choices']['fault_on_s'] - 1 / 60 <= centre <= report['choices']['fault_off_s'] + 1 / 60
                if fault_window or samples_seen % decimation == 0:
                    average = [abs(value / window_count) if k < len(buses) else (value / window_count).real
                               for k, value in enumerate(sums)]
                    averaged.append([centre] + average)
    require(first is not None, 'Empty monitor output')
    require(abs(last_time - study['tmax']) < 1e-9, 'Study stopped before its final time')
    require(len(switch_changes) == 2, f'Expected two fault events, got {switch_changes}')
    expected_events = ((report['choices']['fault_on_s'], False), (report['choices']['fault_off_s'], True))
    for actual, expected in zip(switch_changes, expected_events):
        require(abs(actual[0] - expected[0]) < 1e-10 and actual[1] == expected[1], 'Fault event timing')
    require(current_jump < 1e-8, f'Fault inductor current changed at an event: {current_jump}')
    require(open_current < 1e-5, f'Open fault switch draws current: {open_current}')
    require(decay_error < 5e-4, f'Fault discharge differs from the analytical RL decay: {decay_error}')
    for m, data in report['machines'].items():
        for kind, unit in (('p', 'W'), ('q', 'var')):
            error = abs(first[f'Machine_{m}_{kind}'] - data[f'dispatch_{unit}']) / data['rating_VA']
            pq_initial_error = max(pq_initial_error, error)
    for plant, data in report['inverters'].items():
        filt = next(d for d in case['devices'] if d['id'] == plant + '_filter')
        v = [first[f'Bus_{filt["inputs"]["bus"]}_v{p}'] for p in 'abc']
        i = [first[f'Filter_{plant}_filter_ig{p}'] for p in 'abc']
        power = {'p': sum(a * b for a, b in zip(v, i)),
                 'q': ((v[1] - v[2]) * i[0] + (v[2] - v[0]) * i[1] + (v[0] - v[1]) * i[2]) / math.sqrt(3)}
        for kind, unit in (('p', 'W'), ('q', 'var')):
            error = abs(power[kind] - data[f'dispatch_{unit}']) / data['rating_VA']
            pq_initial_error = max(pq_initial_error, error)
    # Allow the study's relative accuracy in post-IC dispatch on each plant base.
    initial_dispatch_tolerance = study['rel_tol']
    require(pq_initial_error <= initial_dispatch_tolerance, f'Initial dispatch error: {pq_initial_error}')
    # Allow numerical interpolation error in the monitored algebraic limiter.
    current_limit_tolerance = 1e-4
    require(imax_ratio <= 1 + current_limit_tolerance, f'Current limiter exceeded: {imax_ratio}')
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
        'fault': {'inductance_H': fault_params['L'][0][0], 'discharge_resistance_ohm': discharge_params['R'][0][0],
                  'discharge_time_constant_s': 1 / decay_rate, 'event_current_jump_pu': current_jump,
                  'open_switch_current_pu': open_current, 'maximum_discharge_current_relative_error': decay_error},
        'initial_dispatch_max_error_pu_plant_base': pq_initial_error,
        'initial_dispatch_tolerance_pu_plant_base': initial_dispatch_tolerance,
        'maximum_limited_current_ratio': imax_ratio,
        'limited_current_ratio_tolerance': current_limit_tolerance,
        'machine_speed_deviation_pu_range': omega_range, 'fault_bus_1_minimum_voltage_pu': v1_min,
        'final_cycle_voltage_pu_range': [min(final_voltage), max(final_voltage)],
        'accepted_steps': {'count': len(steps), 'minimum_s': min(steps), 'median_s': statistics.median(steps),
                           'maximum_s': max(steps), 'orders': dict(orders)},
        'averaging': {'cycles': 1, 'samples_per_cycle': window_count,
                      'fault_output_interval_s': sample_period, 'outside_fault_output_interval_s': decimation * sample_period,
                      'description': '60 Hz demodulated voltage vector and P/Q/speed cycle means, labelled at window centre'},
    })
    metrics['averaged_sha256'] = hashlib.sha256((output / 'Hawaii.averaged.csv').read_bytes()).hexdigest()
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--reuse', type=Path, help='Completed run directory with immutable inputs and run.json')
    parser.add_argument('--tmax', type=float, default=1.5, help='CTest covers inception, clearing, and recovery; full study is 5 s')
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
