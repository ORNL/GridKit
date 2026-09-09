#!/usr/bin/env python3
"""Validate PLL-synchronized current and voltage control (stdlib only)."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent


def wiring_checks(case):
    """Keep the complete LCL feedback chain aligned with the component ports."""
    devices = {device['id']: device for device in case['devices']}
    assert not any(d['class'] in ('Angle', 'Modulation') for d in devices.values())
    filt, pll, pwm, bridge, inner = (devices[key] for key in
                                       ('filter', 'pll', 'pwm', 'bridge', 'current_control'))
    voltage, current, grid_current = (devices[key] for key in ('voltage', 'current', 'grid_current'))
    outer = devices.get('power_control', devices.get('voltage_control'))
    assert filt['inputs']['e'] == bridge['outputs']['e']
    assert current['inputs']['input'] == filt['outputs']['i']
    assert voltage['inputs']['input'] == filt['outputs']['vo']
    bus = devices[filt['inputs']['bus']]
    assert [pll['inputs']['v' + p] for p in 'abc'] == [bus['outputs']['v' + p] for p in 'abc']
    assert grid_current['inputs']['input'] == filt['outputs']['ig']
    assert inner['inputs']['i'] == current['outputs']['out'][:2]
    assert inner['inputs']['v'] == voltage['outputs']['out'][:2]
    grid_input = 'ig' if outer['class'] == 'OuterVoltageControl' else 'i'
    assert outer['inputs'][grid_input] == grid_current['outputs']['out'][:2]
    if outer['class'] == 'OuterVoltageControl':
        assert outer['inputs']['v'] == voltage['outputs']['out'][:2]
        assert outer['inputs']['omega'] == pll['outputs']['omega']
    else:
        terminal_voltage = devices['terminal_voltage']
        assert terminal_voltage['inputs']['input'] == [bus['outputs']['v' + p] for p in 'abc']
        assert outer['inputs']['v'] == terminal_voltage['outputs']['out'][:2]
        assert terminal_voltage['inputs']['theta'] == pll['outputs']['theta']
    assert all(d['inputs']['theta'] == pll['outputs']['theta'] for d in (voltage, current, grid_current, pwm))
    assert inner['inputs']['omega'] == pll['outputs']['omega']
    assert pwm['inputs']['u'] == inner['outputs']['u']
    assert inner['inputs']['ulim'] == pwm['outputs']['ulim']
    assert inner['inputs']['icmd'] == outer['outputs']['icmd']
    assert outer['inputs']['ilim'] == inner['outputs']['ilim']
    assert bridge['inputs']['s'] == pwm['outputs']['s']
    assert bridge['inputs']['vdc'] == pwm['inputs']['vdc']
    constants = {s['id']: s['value'] for s in case['signals'] if 'value' in s}
    assert constants[bridge['inputs']['vdc']] > 0


def steady_state(case, active, reactive):
    """Balanced LCL solution for the specified terminal powers."""
    devices = {d['id']: d for d in case['devices']}
    grid = devices['grid']['params']
    omega = grid['omega']
    source = math.sqrt(3) * grid['E'][0]
    filt = devices['filter']['params']
    z = complex(filt['Rg'][0][0], omega * filt['Lg'][0][0])
    current = complex(active, -reactive) / source
    drop = z * current
    voltage = source + drop
    c = filt['C'][0][0]
    inverter_current = current + 1j * omega * c * voltage
    return {'voltage_V': abs(voltage), 'vq_V': voltage.imag, 'id_A': inverter_current.real,
            'iq_A': inverter_current.imag, 'p_W': active,
            'q_var': reactive}


def measure(path, begin, end, final_time=None):
    integral = dict.fromkeys(['voltage_V', 'id_A', 'iq_A', 'p_W', 'q_var', 'vq_V', 'frequency_Hz'], 0.0)
    previous = None
    count = 0
    max_identity = 0.0
    max_limited = 0.0
    with path.open() as stream:
        for raw in csv.DictReader(stream):
            row = {key: float(value) for key, value in raw.items()}
            if not all(math.isfinite(value) for value in row.values()):
                raise AssertionError('Nonfinite monitor value')
            t = row['t']
            v = [row[f'Filter_filter_vo{p}'] for p in 'abc']
            i = [row[f'Filter_filter_ig{p}'] for p in 'abc']
            power = sum(a * b for a, b in zip(v, i))
            reactive = ((v[1]-v[2])*i[0] + (v[2]-v[0])*i[1] + (v[0]-v[1])*i[2]) / math.sqrt(3)
            vd, vq = row['Park_voltage_y1'], row['Park_voltage_y2']
            id_, iq = row['Park_grid_current_y1'], row['Park_grid_current_y2']
            max_identity = max(max_identity, abs(power - vd * id_ - vq * iq),
                               abs(reactive - vq * id_ + vd * iq))
            max_limited = max(max_limited, math.hypot(row['InnerCurrentControl_current_control_ilimd'],
                                                    row['InnerCurrentControl_current_control_ilimq']))
            bus_v = [row[f'Bus_terminal_v{p}'] for p in 'abc']
            power = sum(a * b for a, b in zip(bus_v, i))
            reactive = ((bus_v[1]-bus_v[2])*i[0] + (bus_v[2]-bus_v[0])*i[1]
                        + (bus_v[0]-bus_v[1])*i[2]) / math.sqrt(3)
            values = {'voltage_V': math.sqrt(sum(x*x for x in v)),
                      'id_A': row['Park_current_y1'], 'iq_A': row['Park_current_y2'],
                      'p_W': power, 'q_var': reactive, 'vq_V': vq,
                      'frequency_Hz': row['PLL_pll_omega'] / (2 * math.pi)}
            if previous:
                old_t, old = previous
                assert t >= old_t, 'Decreasing monitor time'
                lo, hi = max(begin, old_t), min(end, t)
                if hi > lo:
                    for key in integral:
                        slope = (values[key] - old[key]) / (t - old_t)
                        integral[key] += (hi-lo)*(old[key] + slope*((lo+hi)/2-old_t))
            previous = t, values
            count += 1
    assert abs(previous[0] - (end if final_time is None else final_time)) < 1e-10, 'Incomplete simulation'
    return {'mean': {key: value / (end-begin) for key, value in integral.items()},
            'samples': count, 'measurement_identity_error': max_identity,
            'maximum_limited_current_A': max_limited}


def validate(exe, output):
    solver = json.loads((HERE / 'GFL.solver.json').read_text())
    case_file = (HERE / solver['system_model_file']).resolve()
    state_file = (HERE / solver['state_file']).resolve()
    case = json.loads(case_file.read_text())
    wiring_checks(case)
    params = next(d['params'] for d in case['devices'] if d['id'] == 'power_control')
    report = {}
    for name, mu in [('steady_smooth', 240), ('steady_switching', 1e6)]:
        config = dict(solver, system_model_file=str(case_file), state_file=str(state_file),
                      events=[], tmax=0.3, dt_monitor=2e-6, mu=mu,
                      rel_tol=1e-8, abs_tol=1e-9, output_file=str(output / f'{name}.csv'))
        path = output / f'{name}.solver.json'
        path.write_text(json.dumps(config, indent=2) + '\n')
        with (output / f'{name}.log').open('w') as log:
            subprocess.run([str(exe), str(path)], cwd=output, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        measured = measure(output / f'{name}.csv', 0.25, 0.3)
        expected = steady_state(case, params['Pref'], params['Qref'])
        errors = {key: abs(measured['mean'][key] - value) for key, value in expected.items()}
        measured.update(expected=expected, absolute_errors=errors, mu=mu, mean_window_s=[0.25,0.3])
        report[name] = measured
        print(name, json.dumps(errors), flush=True)
        # The balanced oracle omits switching ripple.
        bounds = {'p_W': .001, 'q_var': .002, 'voltage_V': .0003,
                  'id_A': .00002, 'iq_A': .00002, 'vq_V': .0003}
        for key, bound in bounds.items():
            assert errors[key] < bound, (name, key, errors[key], bound)
        assert measured['measurement_identity_error'] < 1e-8
        assert measured['maximum_limited_current_A'] <= 30 + 1e-10
    (output / 'metrics.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def validate_voltage(exe, output):
    """Check phase-domain voltage tracking before, during, and after a reference step."""
    solver = json.loads((HERE / 'GFM.solver.json').read_text())
    final_time = max(solver['tmax'], 0.5)
    case_file = (HERE / solver['system_model_file']).resolve()
    state_file = (HERE / solver['state_file']).resolve()
    case = json.loads(case_file.read_text())
    wiring_checks(case)
    devices = {d['id']: d for d in case['devices']}
    constants = {s['id']: s['value'] for s in case['signals'] if 'value' in s}
    reference = [constants[name] for name in devices['voltage_control']['inputs']['vref']]
    frequency = devices['grid']['params']['omega'] / (2 * math.pi)
    # Each interval ends at an event or the final time. Use its last 20 ms
    # to measure tracking after the controller has responded.
    ends = sorted({event['time'] for event in solver['events']} | {solver['tmax'], final_time})
    references = []
    for end in ends:
        target = reference.copy()
        for event in solver['events']:
            if event['time'] < end:
                target['dq'.index(event['signal_id'][-1])] = event['value']
        references.append(target)
    report = {}
    for name, mu in [('voltage_smooth', 240), ('voltage_switching', 1e6)]:
        waveform = output / f'{name}.csv'
        config = dict(solver, system_model_file=str(case_file), state_file=str(state_file),
                      tmax=final_time, dt_monitor=2e-6, mu=mu, output_file=str(waveform))
        path = output / f'{name}.solver.json'
        path.write_text(json.dumps(config, indent=2) + '\n')
        with (output / f'{name}.log').open('w') as log:
            subprocess.run([str(exe), str(path)], cwd=output, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        intervals = []
        for end, target in zip(ends, references):
            measured = measure(waveform, end - .02, end, final_time)
            errors = {'voltage_V': abs(measured['mean']['voltage_V'] - math.hypot(*target)),
                      'vq_V': abs(measured['mean']['vq_V'] - target[1]),
                      'frequency_Hz': abs(measured['mean']['frequency_Hz'] - frequency)}
            measured.update(reference_V=target, absolute_errors=errors, mean_window_s=[end - .02, end])
            intervals.append(measured)
            print(name, end, json.dumps(errors), flush=True)
            # Allow 20% of the 1 V step during recovery; check settled tracking separately.
            voltage_bound = .05 if end == final_time else .2
            assert errors['voltage_V'] < voltage_bound, (name, end, errors)
            assert errors['vq_V'] < voltage_bound, (name, end, errors)
            assert errors['frequency_Hz'] < .05, (name, end, errors)
            assert measured['measurement_identity_error'] < 1e-8
            assert measured['maximum_limited_current_A'] <= devices['current_control']['params']['Imax'] + 1e-10
        report[name] = {'mu': mu, 'intervals': intervals}
    (output / 'metrics.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--scenario', choices=['GFL', 'GFM'], default='GFL')
    args = parser.parse_args()
    validator = validate if args.scenario == 'GFL' else validate_voltage
    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)
        validator(args.exe.resolve(), args.output.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='gridkit-current-control-') as directory:
            validator(args.exe.resolve(), Path(directory))


if __name__ == '__main__':
    main()
