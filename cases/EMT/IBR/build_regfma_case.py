#!/usr/bin/env python3
"""Replace the ten-bus PWM assemblies with REGFMA sources at the same operating point."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    case = json.loads((HERE / 'TenBus.case.json').read_text())
    state = json.loads((HERE / 'TenBus.state.json').read_text())
    case['header'] = {
        'case_name': 'Ten-bus EMT grid with REGFMA',
        'case_description': 'Three governed synchronous machines and three REGFM_A1 sources',
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
        case['devices'].append({
            'class': 'REGFMA', 'id': f'regfma_{bus}',
            'params': {'S': 5e6, 'V': 13800.0, 'XL': 0.15,
                       'mp': 0.01, 'mq': 0.05, 'ImaxF': 2.0,
                       'VFlag': True, 'QVFlag': True},
            'inputs': {'bus': f'bus_{bus}'},
            'mon': ['i', 'pf', 'qf', 'vf', 'omega', 'edroop', 'delta', 'p', 'q', 'v']})
        state['devices'][f'regfma_{bus}'] = state['devices'].pop(f'filter_{bus}')
    case['signals'] = [signal for signal in case['signals']
                       if signal['id'].startswith(('speed_', 'pmech_'))]
    state['header']['description'] = (
        'TenBus operating-point estimate with unchanged terminal currents; '
        'REGFMA initializes its internal voltage and references from these currents.')
    for name, data in (('REGFMA.case.json', case), ('REGFMA.state.json', state)):
        (HERE / name).write_text(json.dumps(data, indent=4) + '\n')
    print('Wrote REGFMA.case.json and REGFMA.state.json (10 buses, 3 machines, 3 REGFMA sources).')


if __name__ == '__main__':
    main()
