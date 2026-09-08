#!/usr/bin/env python3
"""Build the ten-bus REGFMA study with physical, damped RLC filters."""
import cmath
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
OMEGA = 2 * math.pi * 60
POWER = 5e6
VOLTAGE = 13800.0
XL = 0.15
RL = 0.03
RESONANCE = 600.0
BASE_IMPEDANCE = VOLTAGE**2 / POWER
INDUCTANCE = XL * BASE_IMPEDANCE / OMEGA
CAPACITANCE = 1 / ((2 * math.pi * RESONANCE)**2 * INDUCTANCE)
DAMPING = 0.003 * BASE_IMPEDANCE


def diagonal(value):
    return [[value if row == column else 0.0 for column in range(3)] for row in range(3)]


def phasor(values, prefix):
    """Recover a phase RMS phasor from balanced instantaneous samples at t=0."""
    return math.sqrt(2) / 3 * sum(
        values[prefix + phase] * cmath.exp(2j * math.pi * n / 3)
        for n, phase in enumerate('abc'))


def phases(value, prefix):
    return {prefix + phase: math.sqrt(2) * (value * cmath.exp(-2j * math.pi * n / 3)).real
            for n, phase in enumerate('abc')}


def main():
    case = json.loads((HERE / 'TenBus.case.json').read_text())
    state = json.loads((HERE / 'TenBus.state.json').read_text())
    case['header'] = {
        'case_name': 'Ten-bus EMT grid with REGFMA',
        'case_description': 'Three governed synchronous machines and three REGFM_A1 sources with damped RLC filters',
        'case_comments': '13.8 kV, 60 Hz; the network and initial terminal injections follow TenBus.case.json.'}
    removed = {f'{kind}_{bus}' for bus in (4, 5, 6) for kind in ('pwm', 'converter', 'filter')}
    case['devices'] = [device for device in case['devices'] if device['id'] not in removed]
    for device in case['devices']:
        kind, name = device['class'], device['id']
        if kind == 'Bus':
            device['mon'] = ['va', 'vb', 'vc']
            if name in ('bus_4', 'bus_5', 'bus_6'):
                device['mon'] += ['i_sha', 'i_shb', 'i_shc']
        elif kind == 'Machine':
            device['mon'] = ['omega', 'p', 'q']
        elif kind == 'Switch':
            device['mon'] = ['open']
        elif name in ('line_4_7', 'line_5_8', 'line_6_8'):
            device['mon'] = ['i12a', 'i12b', 'i12c']
        elif name in ('load_4', 'load_5', 'load_6'):
            device['mon'] = ['ia', 'ib', 'ic']
        else:
            device.pop('mon', None)
    for bus in (4, 5, 6):
        # REGFMA owns the single physical RL coupling branch. The series RC
        # shunt is explicit so its capacitor voltage has a prescribed state.
        case['devices'] += [
            {'class': 'Bus', 'id': f'capacitor_{bus}',
             'shunts': {'C': {'E': diagonal(CAPACITANCE)}},
             'mon': ['va', 'vb', 'vc', 'i_sha', 'i_shb', 'i_shc']},
            {'class': 'LineLumped', 'id': f'damping_{bus}',
             'params': {'dx': 1.0, 'Rp': diagonal(DAMPING), 'Lp': diagonal(0.0)},
             'inputs': {'bus1': f'bus_{bus}', 'bus2': f'capacitor_{bus}'},
             'mon': ['i12a', 'i12b', 'i12c']}]
        case['devices'].append({
            'class': 'REGFMA', 'id': f'regfma_{bus}',
            'params': {'S': POWER, 'V': VOLTAGE, 'XL': XL, 'RL': RL,
                       'mp': 0.01, 'mq': 0.05, 'ImaxF': 2.0,
                       'VFlag': True, 'QVFlag': True},
            'inputs': {'bus': f'bus_{bus}'},
            'mon': ['i', 'e', 'pf', 'qf', 'vf', 'omega', 'edroop', 'delta', 'p', 'q', 'v']})
        voltage = phasor(state['buses'][f'bus_{bus}'], 'v')
        net_current = phasor(state['devices'].pop(f'filter_{bus}'), 'i')
        capacitor_voltage = voltage / (1 + 1j * OMEGA * DAMPING * CAPACITANCE)
        capacitor_current = 1j * OMEGA * CAPACITANCE * capacitor_voltage
        state['buses'][f'capacitor_{bus}'] = phases(capacitor_voltage, 'v')
        state['devices'][f'damping_{bus}'] = phases(capacitor_current, 'i12')
        state['devices'][f'regfma_{bus}'] = phases(net_current + capacitor_current, 'i')
    case['signals'] = [signal for signal in case['signals']
                       if signal['id'].startswith(('speed_', 'pmech_'))]
    state['header']['description'] = (
        'TenBus operating-point estimate with unchanged net network injections; '
        'REGFMA currents include the damped capacitor shunts at 60 Hz.')
    for name, data in (('REGFMA.case.json', case), ('REGFMA.state.json', state)):
        (HERE / name).write_text(json.dumps(data, indent=4) + '\n')
    print('Wrote REGFMA.case.json and REGFMA.state.json (10 network buses, 3 capacitor buses, 3 machines, 3 REGFMA sources).')


if __name__ == '__main__':
    main()
