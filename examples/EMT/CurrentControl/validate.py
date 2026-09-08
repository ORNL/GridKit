#!/usr/bin/env python3
"""Validate GFL current control against the balanced LCL steady state (stdlib only)."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent


def steady_state(case, irefd, irefq):
    """Balanced LCL solution with supplied current in the terminal-voltage frame."""
    devices = {d['id']: d for d in case['devices']}
    grid = devices['grid']['params']
    omega = grid['omega']
    source = math.sqrt(3) * grid['E'][0]
    line = devices['grid_filter']['params']
    z = complex(line['Rp'][0][0], omega * line['Lp'][0][0]) * line['dx']
    current = complex(irefd, irefq)
    drop = z * current
    voltage = math.sqrt(source * source - drop.imag * drop.imag) + drop.real
    c = devices['capacitor']['shunts']['C']['E'][0][0]
    inverter_current = current + 1j * omega * c * voltage
    return {'voltage_V': voltage, 'id_A': inverter_current.real,
            'iq_A': inverter_current.imag, 'p_W': voltage * irefd,
            'q_var': -voltage * irefq}


def measure(path, begin, end):
    integral = dict.fromkeys(['voltage_V', 'id_A', 'iq_A', 'p_W', 'q_var'], 0.0)
    previous = None
    count = 0
    max_identity = 0.0
    max_limited = 0.0
    minimum_dc = math.inf
    with path.open() as stream:
        for raw in csv.DictReader(stream):
            row = {key: float(value) for key, value in raw.items()}
            if not all(math.isfinite(value) for value in row.values()):
                raise AssertionError('Nonfinite monitor value')
            t = row['t']
            v = [row[f'Bus_capacitor_v{p}'] for p in 'abc']
            i = [row[f'LineLumped_grid_filter_i12{p}'] for p in 'abc']
            power = sum(a * b for a, b in zip(v, i))
            reactive = ((v[1]-v[2])*i[0] + (v[2]-v[0])*i[1] + (v[0]-v[1])*i[2]) / math.sqrt(3)
            vd, vq = row['Park_voltage_y1'], row['Park_voltage_y2']
            id_, iq = row['Park_grid_current_y1'], row['Park_grid_current_y2']
            max_identity = max(max_identity, abs(power - vd * id_ - vq * iq),
                               abs(reactive - vq * id_ + vd * iq))
            max_limited = max(max_limited, math.hypot(row['InnerCurrentControl_current_control_ilimd'],
                                                    row['InnerCurrentControl_current_control_ilimq']))
            minimum_dc = min(minimum_dc, row['DCLink_dc_vdc'])
            values = {'voltage_V': math.sqrt(sum(x*x for x in v)),
                      'id_A': row['Park_current_y1'], 'iq_A': row['Park_current_y2'],
                      'p_W': power, 'q_var': reactive}
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
    assert abs(previous[0] - end) < 1e-10, 'Incomplete simulation'
    return {'mean': {key: value / (end-begin) for key, value in integral.items()},
            'samples': count, 'measurement_identity_error': max_identity,
            'maximum_limited_current_A': max_limited, 'minimum_dc_voltage_V': minimum_dc}


def validate(exe, output):
    solver = json.loads((HERE / 'GFL.solver.json').read_text())
    case_file = (HERE / solver['system_model_file']).resolve()
    state_file = (HERE / solver['state_file']).resolve()
    case = json.loads(case_file.read_text())
    params = next(d['params'] for d in case['devices'] if d['id'] == 'power_control')
    irefd, irefq = params['Pref'] / params['V'], -params['Qref'] / params['V']
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
        expected = steady_state(case, irefd, irefq)
        errors = {key: abs(measured['mean'][key] - value) for key, value in expected.items()}
        measured.update(expected=expected, absolute_errors=errors, mu=mu, mean_window_s=[0.25,0.3])
        report[name] = measured
        print(name, json.dumps(errors), flush=True)
        # The balanced oracle omits switching ripple. Repeated current-reference
        # runs differ by 0.001194 var at resolved switching; allow 0.002 var.
        bounds = {'p_W': .001, 'q_var': .002, 'voltage_V': .0003,
                  'id_A': .00002, 'iq_A': .00002}
        for key, bound in bounds.items():
            assert errors[key] < bound, (name, key, errors[key], bound)
        assert measured['measurement_identity_error'] < 1e-8
        assert measured['maximum_limited_current_A'] <= 30 + 1e-10
        assert measured['minimum_dc_voltage_V'] > 0
    (output / 'metrics.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)
        validate(args.exe.resolve(), args.output.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='gridkit-gfl-') as directory:
            validate(args.exe.resolve(), Path(directory))


if __name__ == '__main__':
    main()
