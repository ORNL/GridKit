"""Convert the frozen PhasorDynamics Hawaii JSON using only the standard library."""

import argparse
import cmath
import collections
import hashlib
import json
import math
from pathlib import Path
import subprocess

SOURCE_REVISION = 'f16f3815e6f84c08bd0869773d3033b3dcb3dd0f'
SOURCE_PATH = 'cases/PhasorDynamics/Hawaii/Hawaii.case.json'
ROOT = Path(__file__).resolve().parent
SYSTEM_BASE = 100e6
FREQUENCY = 60.0
OMEGA = 2 * math.pi * FREQUENCY
CHOICES = {
    'damper_gap_fraction': 0.99,
    'transformer_I0': 0.001,
    'transformer_P0_W': 0.0,
    'transformer_knee': 1.2,
    'transformer_Lsat': 0.25,
    'filter_R_pu': 0.01,
    'filter_X_pu': 0.20,
    'inner_bandwidth_Hz': 300.0,
    'Mmax': 0.95,
    'PLL_Kp': 80.0,
    'PLL_Ki': 2500.0,
    'outer_Kp': 0.01,
    'outer_Ki': 40.0,
    'outer_Kaw': 200.0,
    'dc_voltage_ratio': 2.0,
    'dc_energy_seconds': 10.0,
    'carrier_Hz': 1800.0,
    'carrier_alignment': 0.5,
    'fault_R_pu': 0.01,
    'fault_on_s': 1.0,
    'fault_off_s': 1.1,
}


def diagonal(value):
    return [[value if i == j else 0.0 for j in range(3)] for i in range(3)]


def samples(value, prefix, scale=1.0):
    return {prefix + phase: scale * (value * cmath.exp(1j * angle)).real
            for phase, angle in zip('abc', (0.0, -2 * math.pi / 3, 2 * math.pi / 3))}


def solve_linear(matrix, rhs):
    """Pivoted elimination, used for the 72-variable rectangular power flow."""
    a = [list(row) + [value] for row, value in zip(matrix, rhs)]
    n = len(a)
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(a[i][k]))
        a[k], a[pivot] = a[pivot], a[k]
        if abs(a[k][k]) < 1e-14:
            raise ValueError('Singular conversion power-flow Jacobian')
        for i in range(k + 1, n):
            factor = a[i][k] / a[k][k]
            for j in range(k + 1, n + 1):
                a[i][j] -= factor * a[k][j]
    result = [0.0] * n
    for i in reversed(range(n)):
        result[i] = (a[i][n] - sum(a[i][j] * result[j] for j in range(i + 1, n))) / a[i][i]
    return result


def transformer_point(p, v1, v2):
    """Independent linear core/leakage circuit at the fundamental frequency."""
    ym = -0.5j * CHOICES['transformer_I0']
    ys = 1 / (1j * p['X'])
    r = p['R'] / 2
    a, b = 1 + r * (ym + ys), -r * ys
    e1 = (a * v1 - b * v2) / (a * a - b * b)
    e2 = (a * v2 - b * v1) / (a * a - b * b)
    i12 = ys * (e1 - e2)
    return e1, e2, i12, ym * e1 + i12, ym * e2 - i12


def axis_windings(p, axis):
    x, xp, xpp = p['X' + axis], p['X' + axis + 'p'], p['X' + axis + 'pp']
    xl = p['Xl']
    effective = xpp
    if xp == xpp:
        effective = xl + CHOICES['damper_gap_fraction'] * (xp - xl)
    if not x > xp > effective > xl > 0:
        raise ValueError(f'Unsupported {axis}-axis reactance ordering: {p}')
    lm = x - xl
    l1 = lm * (xp - xl) / (x - xp)
    l2 = (xp - xl) * (effective - xl) / (xp - effective)
    t1, t2 = p['T' + axis + 'op'], p['T' + axis + 'opp']
    r1 = (lm + l1) / (OMEGA * t1)
    r2 = (l2 + lm * l1 / (lm + l1)) / (OMEGA * t2)
    # Eigenvalues of -omega * inv(L_rotor) * R at zero stator current.
    determinant = (lm + l1) * (lm + l2) - lm * lm
    trace = OMEGA * ((lm + l2) * r1 + (lm + l1) * r2) / determinant
    product = OMEGA**2 * r1 * r2 / determinant
    fast = (trace + math.sqrt(trace * trace - 4 * product)) / 2
    slow = product / fast
    record = {
        'source_X': x, 'source_Xp': xp, 'source_Xpp': xpp,
        'effective_Xpp': effective, 'source_open_times_s': [t1, t2],
        'winding_open_times_s': [1 / slow, 1 / fast],
        'reconstructed_Xp': xl + 1 / (1 / lm + 1 / l1),
        'reconstructed_Xpp': xl + 1 / (1 / lm + 1 / l1 + 1 / l2),
    }
    return (lm, l1, l2, r1, r2), record


def operating_point(source):
    buses = {b['number']: b for b in source['buses']}
    numbers = list(buses)
    original = {n: complex(b['init']['Vr'], b['init']['Vi']) for n, b in buses.items()}
    y = {i: {j: 0j for j in numbers} for i in numbers}
    dispatch = {n: 0j for n in numbers}
    for d in source['devices']:
        p, ports = d['params'], d['ports']
        if d['class'] == 'Branch':
            i, j = ports['bus1'], ports['bus2']
            assert p['tap'] == 1 and p['phase'] == 0
            if d.get('extension', {}).get('xfmr', False):
                point1 = transformer_point(p, 1, 0)
                point2 = transformer_point(p, 0, 1)
                y[i][i] += point1[3]
                y[j][i] += point1[4]
                y[i][j] += point2[3]
                y[j][j] += point2[4]
            else:
                series = 1 / complex(p['R'], p['X'])
                y[i][i] += series + 0.5j * p['B']
                y[j][j] += series + 0.5j * p['B']
                y[i][j] -= series
                y[j][i] -= series
        elif d['class'] == 'LoadZIP':
            assert p['alphaI'] == p['alphaP'] == 0
            n = ports['bus']
            y[n][n] += complex(p['Pnom'], -p['Qnom']) / abs(original[n])**2
        elif d['class'] in ('Genrou', 'Regca'):
            dispatch[ports['bus']] += complex(p['p0'], p['q0'])
    unknown = [n for n in numbers if n != 23]
    voltage = dict(original)

    def residual(v):
        mismatch = [v[n] * sum(y[n][j] * v[j] for j in numbers).conjugate() - dispatch[n]
                    for n in unknown]
        return [part for value in mismatch for part in (value.real, value.imag)]

    for iteration in range(12):
        f = residual(voltage)
        if max(map(abs, f)) < 1e-11:
            break
        columns = []
        h = 1e-6
        for n in unknown:
            for direction in (1.0, 1j):
                perturbed = dict(voltage)
                perturbed[n] += h * direction
                columns.append([(a - b) / h for a, b in zip(residual(perturbed), f)])
        correction = solve_linear(list(zip(*columns)), [-value for value in f])
        for k, n in enumerate(unknown):
            voltage[n] += complex(correction[2 * k], correction[2 * k + 1])
    else:
        raise ValueError('Conversion power flow did not converge')
    powers = {n: voltage[n] * sum(y[n][j] * voltage[j] for j in numbers).conjugate()
              for n in numbers}
    delta = powers[23] - dispatch[23]
    slack = [d for d in source['devices'] if d['class'] == 'Genrou' and d['ports']['bus'] == 23]
    total = sum(d['params']['mva'] for d in slack)
    adjusted = {}
    for d in source['devices']:
        if d['class'] in ('Genrou', 'Regca'):
            p = d['params']
            adjusted[d['id']] = complex(p['p0'], p['q0'])
            if d in slack:
                adjusted[d['id']] += delta * p['mva'] / total
    net = {n: sum(adjusted[d['id']] for d in source['devices']
                  if d['id'] in adjusted and d['ports']['bus'] == n) for n in numbers}
    report = {
        'iterations': iteration,
        'max_voltage_change_pu': max(abs(voltage[n] - original[n]) for n in numbers),
        'slack_delta_W': delta.real * SYSTEM_BASE,
        'slack_delta_var': delta.imag * SYSTEM_BASE,
        'max_KCL_pu': max(abs((powers[n] - net[n]) / voltage[n]) for n in numbers),
        'bus_phasors_pu': {str(n): [voltage[n].real, voltage[n].imag] for n in numbers},
    }
    return buses, voltage, adjusted, report


def convert(source):
    assert source['params']['va_base'] == SYSTEM_BASE
    assert source['params']['freq_base'] == FREQUENCY
    buses, voltage, dispatch, flow_report = operating_point(source)
    devices, signals = [], {}
    state = {'header': {'version': 1, 'time': 0.0,
                        'description': 'Balanced fundamental operating point; see conversion.json'},
             'buses': {}, 'devices': {}}
    report = {'source_revision': SOURCE_REVISION, 'source_path': SOURCE_PATH,
              'choices': CHOICES, 'source_counts': dict(collections.Counter(d['class'] for d in source['devices'])),
              'power_flow': flow_report, 'machines': {}, 'inverters': {}}

    def signal(name, value=None):
        signals.setdefault(name, {'id': name})
        if value is not None:
            signals[name]['value'] = value
        return name

    def add(kind, name, params=None, inputs=None, outputs=None, mon=None, **extra):
        d = {'class': kind, 'id': name}
        for key, value in [('params', params), ('inputs', inputs), ('outputs', outputs), ('mon', mon)]:
            if value:
                d[key] = value
        d.update(extra)
        devices.append(d)
        return d

    bus_devices = {}
    for n, b in buses.items():
        vb = b['params']['kv'] * 1000
        name = f'bus_{n}'
        bus_devices[n] = add('Bus', name, outputs={f'v{p}': signal(f'b{n}_v{p}') for p in 'abc'},
                             mon=['va', 'vb', 'vc'])
        state['buses'][name] = samples(voltage[n], 'v', math.sqrt(2 / 3) * vb)
    for d in source['devices']:
        p, ports, name = d['params'], d['ports'], d['id']
        if d['class'] == 'Branch':
            i, j = ports['bus1'], ports['bus2']
            v1, v2 = buses[i]['params']['kv'] * 1000, buses[j]['params']['kv'] * 1000
            if d.get('extension', {}).get('xfmr', False):
                add('Transformer', name,
                    {'S': SYSTEM_BASE, 'V1': v1, 'V2': v2, 'f': FREQUENCY,
                     'R': p['R'], 'X': p['X'], 'I0': CHOICES['transformer_I0'],
                     'P0': CHOICES['transformer_P0_W'], 'knee': CHOICES['transformer_knee'],
                     'Lsat': CHOICES['transformer_Lsat'], 'tap': 1.0, 'split': 0.5},
                    {'bus1': f'bus_{i}', 'bus2': f'bus_{j}'})
                e1, e2, current, _, _ = transformer_point(p, voltage[i], voltage[j])
                state['devices'][name] = dict(samples(current, 'i12'),
                                              **samples(-1j * e1, 'psi1'), **samples(-1j * e2, 'psi2'))
            else:
                assert v1 == v2
                zb = v1**2 / SYSTEM_BASE
                add('LineLumped', name, {'dx': 1.0, 'Rp': diagonal(p['R'] * zb),
                                        'Lp': diagonal(p['X'] * zb / OMEGA),
                                        'Cp': diagonal(p['B'] / (OMEGA * zb))},
                    {'bus1': f'bus_{i}', 'bus2': f'bus_{j}'})
                current = (voltage[i] - voltage[j]) / complex(p['R'], p['X'])
                state['devices'][name] = samples(current, 'i12', math.sqrt(2 / 3) * SYSTEM_BASE / v1)
        elif d['class'] == 'LoadZIP':
            n = ports['bus']
            anchor = complex(buses[n]['init']['Vr'], buses[n]['init']['Vi'])
            vb = buses[n]['params']['kv'] * 1000
            add('LoadZ', name, {'R': diagonal(abs(anchor)**2 * vb**2 / (p['Pnom'] * SYSTEM_BASE))},
                {'bus': f'bus_{n}'})
            if p['Qnom']:
                assert p['Qnom'] < 0
                capacitance = -p['Qnom'] * SYSTEM_BASE / (OMEGA * abs(anchor)**2 * vb**2)
                bus_devices[n].setdefault('shunts', {})[name + '_capacitor'] = {
                    'rows': 3, 'cols': 3, 'E': diagonal(capacitance)}
        elif d['class'] == 'Genrou':
            n = ports['bus']
            vb, rating = buses[n]['params']['kv'] * 1000, p['mva'] * 1e6
            wd, rd = axis_windings(p, 'd')
            wq, rq = axis_windings(p, 'q')
            assert p['D'] == 0
            params = {'S': rating, 'V': vb, 'f': FREQUENCY, 'H': p['H'], 'F': 0.0,
                      'Rs': p['Ra'], 'Ll': p['Xl'], 'L0': p['Xl'],
                      'Lmd': wd[0], 'Llfd': wd[1], 'Ll1d': wd[2], 'Rfd': wd[3], 'R1d': wd[4],
                      'Lmq': wq[0], 'Ll1q': wq[1], 'Ll2q': wq[2], 'R1q': wq[3], 'R2q': wq[4],
                      'S10': p['S10'], 'S12': p['S12']}
            add('Machine', name, params,
                {'bus': f'bus_{n}', 'pm': signal(f's{ports["pmech"]}'), 'efd': signal(f's{ports["efd"]}')},
                {'speed': signal(f's{ports["speed"]}')}, ['omega', 'p', 'q'])
            current = (dispatch[name] * SYSTEM_BASE / (vb * voltage[n])).conjugate()
            state['devices'][name] = samples(current, 'i', math.sqrt(2 / 3))
            report['machines'][name] = {'d': rd, 'q': rq, 'rating_VA': rating,
                                        'dispatch_W': dispatch[name].real * SYSTEM_BASE,
                                        'dispatch_var': dispatch[name].imag * SYSTEM_BASE}
        elif d['class'] == 'Tgov1':
            gen = next(g for g in source['devices'] if g['class'] == 'Genrou'
                       and g['ports']['pmech'] == ports['pmech'])
            assert p['Trate'] == gen['params']['mva']
            add('Tgov1', name, p, {'speed': signal(f's{ports["speed"]}')},
                {'pmech': signal(f's{ports["pmech"]}')})
        elif d['class'] == 'Ieeet1':
            inputs = {'bus': f'bus_{ports["bus"]}', 'speed': signal(f's{ports["speed"]}')}
            for key in ('vs', 'vuel', 'voel'):
                if key in ports:
                    inputs[key] = signal(f's{ports[key]}')
            add('Ieeet1', name, dict(p, V=buses[ports['bus']]['params']['kv'] * 1000), inputs,
                {'efd': signal(f's{ports["efd"]}')})
        elif d['class'] == 'Ieeest':
            assert any(g['class'] == 'Genrou' and g['ports']['speed'] == ports['input']
                       for g in source['devices'])
            add('Ieeest', name, p, {'speed': signal(f's{ports["input"]}')},
                {'output': signal(f's{ports["output"]}')})

    signal('zero', 0.0)
    for d in source['devices']:
        if d['class'] != 'Regca':
            continue
        p, ports = d['params'], d['ports']
        prefix = d['id'].removesuffix('_regca')
        n = ports['bus']
        vb, rating = buses[n]['params']['kv'] * 1000, p['mva'] * 1e6
        r = CHOICES['filter_R_pu'] * vb**2 / rating
        inductance = CHOICES['filter_X_pu'] * vb**2 / (rating * OMEGA)
        wc = 2 * math.pi * CHOICES['inner_bandwidth_Hz']
        vdc = CHOICES['dc_voltage_ratio'] * vb
        capacitance = 2 * CHOICES['dc_energy_seconds'] * rating / vdc**2
        pq = dispatch[d['id']] * SYSTEM_BASE
        current = (pq / (vb * voltage[n])).conjugate()
        current_dq = current * cmath.exp(-1j * cmath.phase(voltage[n]))
        bridge_power = pq.real + r * abs(current)**2
        reecb = next(e for e in source['devices'] if e['id'] == prefix + '_reecb')
        imax = reecb['params']['Imax'] * rating / vb

        def s(key, value=None):
            return signal(prefix + '_' + key, value)

        def vector(key, phases):
            return [s(key + phase) for phase in phases]

        add('DCLink', prefix + '_dc', {'C': capacitance},
            {'isrc': s('isrc', bridge_power / vdc), 'idc': s('idc')}, {'vdc': s('vdc')}, ['vdc'])
        add('DependentVoltageSource', prefix + '_filter', {'Rs': diagonal(r), 'Ls': diagonal(inductance)},
            dict(zip(('ea', 'eb', 'ec'), vector('e', 'abc')), bus=f'bus_{n}'),
            dict(zip(('ia', 'ib', 'ic'), vector('i', 'abc'))))
        add('PLL', prefix + '_pll', {'V': vb, 'f': FREQUENCY, 'Kp': CHOICES['PLL_Kp'], 'Ki': CHOICES['PLL_Ki']},
            {'bus': f'bus_{n}'}, {'theta': s('theta'), 'omega': s('omega')}, ['omega'])
        add('Park', prefix + '_voltage', inputs={'input': [f'b{n}_v{phase}' for phase in 'abc'], 'theta': s('theta')},
            outputs={'out': vector('v', 'dq0')}, mon=['out'])
        add('Park', prefix + '_current', inputs={'input': vector('i', 'abc'), 'theta': s('theta')},
            outputs={'out': vector('i', 'dq0')}, mon=['out'])
        add('OuterPowerControl', prefix + '_power',
            {'V': vb, 'Pref': vb * current_dq.real, 'Qref': -vb * current_dq.imag,
             'Kp': CHOICES['outer_Kp'], 'Ki': CHOICES['outer_Ki'], 'Kaw': CHOICES['outer_Kaw']},
            {'i': vector('i', 'dq'), 'ilim': vector('ilim', 'dq')},
            {'icmd': vector('icmd', 'dq')}, ['icmd'])
        add('InnerCurrentControl', prefix + '_inner',
            {'L': inductance, 'Kp': inductance * wc, 'Ki': r * wc, 'Kaw': wc,
             'Imax': imax, 'Mmax': CHOICES['Mmax']},
            {'v': vector('v', 'dq'), 'i': vector('i', 'dq'), 'icmd': vector('icmd', 'dq'),
             'omega': s('omega'), 'vdc': s('vdc')},
            {'ilim': vector('ilim', 'dq'), 'u': vector('u', 'dq')}, ['ilim'])
        add('Park', prefix + '_inverse', {'inverse': True},
            {'input': vector('u', 'dq') + ['zero'], 'theta': s('theta')}, {'out': vector('u', 'abc')})
        add('Modulation', prefix + '_modulation', inputs={'u': vector('u', 'abc'), 'vdc': s('vdc')},
            outputs={'m': vector('m', 'abc')})
        add('PWM', prefix + '_pwm', {'fc': CHOICES['carrier_Hz'], 'alignment': CHOICES['carrier_alignment']},
            {'m': vector('m', 'abc')}, {'s': vector('s', 'abc')})
        add('Converter', prefix + '_bridge', inputs={'s': vector('s', 'abc'), 'vdc': s('vdc'), 'i': vector('i', 'abc')},
            outputs={'vo': vector('e', 'abc'), 'idc': s('idc')})
        state['devices'][prefix + '_dc'] = {'vdc': vdc}
        state['devices'][prefix + '_filter'] = samples(current, 'i', math.sqrt(2 / 3))
        state['devices'][prefix + '_power'] = {'icmdd': current_dq.real, 'icmdq': current_dq.imag}
        command = vb * abs(voltage[n]) + complex(r, OMEGA * inductance) * current_dq
        state['devices'][prefix + '_inner'] = {'ud': command.real, 'uq': command.imag}
        report['inverters'][prefix] = {'rating_VA': rating, 'voltage_V': vb, 'Imax_A': imax,
                                       'dispatch_W': pq.real, 'dispatch_var': pq.imag,
                                       'initial_vdc_V': vdc, 'dc_C_F': capacitance}
    add('Bus', 'fault_bus')
    add('Switch', 'fault_switch', {'open': True}, {'bus1': 'bus_1', 'bus2': 'fault_bus'}, mon=['open'])
    add('LoadZ', 'fault_load', {'R': diagonal(CHOICES['fault_R_pu'] * (buses[1]['params']['kv'] * 1000)**2 / SYSTEM_BASE)},
        {'bus': 'fault_bus'})
    state['devices']['fault_switch'] = {'open': True}
    case = {'header': {'case_name': 'Hawaii EMT',
                       'case_description': '37-bus network with winding machines and nine switching GFL plants',
                       'case_comments': 'Conversion choices and differences from PowerWorld are documented in README.md.'},
            'signals': list(signals.values()), 'devices': devices}
    report['emt_counts'] = dict(collections.Counter(d['class'] for d in devices))
    return case, state, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--output', type=Path, default=ROOT)
    args = parser.parse_args()
    if args.source:
        raw = args.source.read_bytes()
    else:
        raw = subprocess.check_output(['git', 'show', f'{SOURCE_REVISION}:{SOURCE_PATH}'], cwd=ROOT)
    case, state, report = convert(json.loads(raw))
    report['source_sha256'] = hashlib.sha256(raw).hexdigest()
    if args.source:
        report['source_revision'] = None
        report['source_path'] = args.source.name
    args.output.mkdir(parents=True, exist_ok=True)
    for name, data in [('Hawaii.case.json', case), ('Hawaii.state.json', state), ('conversion.json', report)]:
        (args.output / name).write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps(report['power_flow'], indent=2))


if __name__ == '__main__':
    main()
