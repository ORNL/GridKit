"""Write Steady.state.json from the phasor solution of the transformer case."""
import cmath
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

ROTATION = (1.0, cmath.exp(-2j * math.pi / 3), cmath.exp(2j * math.pi / 3))


def solve(case):
    """Per-unit phasors of the bank with identity connection maps and the breaker closed."""
    devices = {device['id']: device for device in case['devices']}
    source, bank, load = (devices[name]['params'] for name in ('source', 'bank', 'load'))
    omega = source['omega']
    tap, split = bank.get('tap', 1.0), bank.get('split', 0.5)
    winding1, winding2 = bank['V1'] / math.sqrt(3), bank['V2'] / math.sqrt(3)
    base1, base2 = 3 * winding1 ** 2 / bank['S'], 3 * winding2 ** 2 / bank['S']
    conductance = bank.get('P0', 0.0) / bank['S']
    magnetizing = complex(conductance, -math.sqrt(bank['I0'] ** 2 - conductance ** 2))
    y1, y2 = split * magnetizing, (1 - split) * magnetizing
    r1 = r2 = bank['R'] / 2
    series = source['Rs'][0][0] / base1 + r1
    load_pu = load['R'][0][0] / base2
    k2 = y2 + 1 / (tap ** 2 * (load_pu + r2))
    e2 = (math.sqrt(2) * source['E'][0] / (math.sqrt(2) * winding1)) \
        / ((1 + 1j * bank['X'] * k2) * (1 + series * y1) + series * k2)
    e1 = e2 * (1 + 1j * bank['X'] * k2)
    i12 = k2 * e2
    return {
        'omega': omega,
        'i12': i12,
        'psi1': -1j * e1,
        'psi2': -1j * e2,
        'iw1': y1 * e1 + i12,
        'iw2': -e2 / (tap * (load_pu + r2)),
        'i_peak1': math.sqrt(2) * bank['S'] / (3 * winding1),
        'i_peak2': math.sqrt(2) * bank['S'] / (3 * winding2),
    }


def samples(phasor, time, omega):
    """Instantaneous phase a, b, and c values of a balanced positive-sequence phasor."""
    return [(phasor * rotation * cmath.exp(1j * omega * time)).real for rotation in ROTATION]


def main():
    case = json.loads((ROOT / 'Transformer.case.json').read_text())
    phasors = solve(case)
    bank = {}
    for name in ('i12', 'psi1', 'psi2'):
        for phase, value in zip('abc', samples(phasors[name], 0.0, phasors['omega'])):
            bank[name + phase] = value
    state = {
        'header': {'version': 1, 'time': 0.0,
                   'description': 'Phasor steady state with the breaker closed'},
        'devices': {'breaker': {'open': False}, 'bank': bank},
    }
    (ROOT / 'Steady.state.json').write_text(json.dumps(state, indent=2) + '\n')


if __name__ == '__main__':
    main()
