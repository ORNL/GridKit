"""Generate coupled 60 Hz line equivalents from GridWorkbench geometry."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import least_squares

from convert import FREQUENCY, ROOT, SYSTEM_BASE, read_source


def transposed(matrix):
    """Uniform equivalent of three equal cyclic phase transpositions."""
    diagonal = np.trace(matrix) / 3
    mutual = (matrix.sum() - np.trace(matrix)) / 6
    return (diagonal - mutual) * np.eye(3) + mutual * np.ones((3, 3))


def positive(matrix):
    return matrix[0, 0] - matrix[0, 1]


def generate(source, gridworkbench):
    sys.path.insert(0, str(gridworkbench))
    from gridworkbench.emt.parameters import EPS0, MU0, Earth, Line, sweep

    omega = 2 * np.pi * FREQUENCY
    buses = {bus['number']: bus['params']['kv'] for bus in source['buses']}

    def parameters(radius, spacing, kv):
        line = Line(
            x=np.array([-spacing, 0., spacing]),
            h=np.full(3, 14. if kv == 69 else 20.),
            r=np.full(3, radius), q=np.full(3, 0.3 * radius),
            sigma=np.full(3, 3.5e7), mu=np.full(3, MU0),
            phase=tuple('abc'), circuit=np.ones(3, dtype=int),
            length=1., earth=Earth(0.01, EPS0),
        )
        data = sweep(line, np.array([omega]))
        return {key: transposed(getattr(data, key)[0]) for key in ('R', 'L', 'G', 'C')}

    geometries, lines, fitted = {}, {}, {}
    for device in source['devices']:
        if device['class'] != 'Branch' or device.get('extension', {}).get('xfmr'):
            continue
        p, name = device['params'], device['id']
        kv = buses[device['ports']['bus1']]
        if kv not in (69, 138) or kv != buses[device['ports']['bus2']]:
            raise ValueError(f'Unsupported line voltage: {name}')
        key = (kv, p['R'], p['X'], p['B'])
        if key in fitted:
            lines[name] = fitted[key]
            continue
        zb = (kv * 1000)**2 / SYSTEM_BASE
        target = np.array([p['R'] * zb, p['X'] * zb / omega, p['B'] / (omega * zb)])
        if not np.all(np.isfinite(target) & (target > 0)):
            raise ValueError(f'Expected positive finite R, X, B: {name}')

        def residual(x):
            radius, spacing, length = np.exp(x)
            matrices = parameters(radius, spacing, kv)
            value = length * np.array([positive(matrices[k]) for k in ('R', 'L', 'C')])
            return np.log(value / target)

        fit = least_squares(
            residual, np.log([0.012, 3., 8000.]),
            bounds=(np.log([0.003, 0.75 if kv == 69 else 1.5, 100.]),
                    np.log([0.04, 15., 150000.])),
            ftol=1e-10, xtol=1e-10, gtol=1e-10,
        )
        if not fit.success:
            raise ValueError(f'Geometry fit failed: {name}: {fit.message}')
        radius, spacing, length = np.exp(fit.x)
        matrices = parameters(radius, spacing, kv)
        scale = {k: value / (length * positive(matrices[k]))
                 for k, value in zip(('R', 'L', 'C'), target)}
        for k, factor in scale.items():
            matrices[k] *= factor
        for k, matrix in matrices.items():
            if np.linalg.eigvalsh(matrix).min() < 0:
                raise ValueError(f'Nonpassive {k} matrix: {name}')
        geometries[name] = {
            'voltage_kV': kv, 'radius_m': float(radius), 'spacing_m': float(spacing),
            'calibration': scale,
            'params': {'dx': float(length), **{k + 'p': m.tolist() for k, m in matrices.items()}},
        }
        fitted[key] = lines[name] = name
    files = sorted((gridworkbench / 'gridworkbench/emt/parameters').glob('*.py'))
    return {
        'frequency_Hz': FREQUENCY,
        'description': 'Geometry-derived coupling; R, L, C scaled independently to preserve source positive sequence.',
        'geometry': {
            'arrangement': 'Ideally transposed horizontal three-wire overhead line; no shield wires.',
            'conductors': 'Synthetic tubular aluminum; inner/outer radius 0.3, conductivity 3.5e7 S/m, relative permeability 1.',
            'height_m': {'69': 14., '138': 20.}, 'earth_resistivity_ohm_m': 100.,
            'fit': 'Bounded least squares of logarithmic positive-sequence R, L, C errors.',
            'radius_bounds_m': [0.003, 0.04], 'spacing_bounds_m': {'69': [0.75, 15.], '138': [1.5, 15.]},
            'length_bounds_m': [100., 150000.],
        },
        'calculator_sha256': {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
        'geometries': geometries, 'lines': lines,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gridworkbench', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--output', type=Path, default=ROOT / 'line_parameters.json')
    args = parser.parse_args()
    raw = read_source(args.source)
    result = generate(json.loads(raw), args.gridworkbench.resolve())
    result['source_sha256'] = hashlib.sha256(raw).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    error = max(abs(1 / scale - 1) for g in result['geometries'].values() for scale in g['calibration'].values())
    print(f"Generated {len(result['lines'])} coupled lines from {len(result['geometries'])} geometries.")
    print(f'Maximum uncalibrated positive-sequence relative error: {error:.6g}')


if __name__ == '__main__':
    main()
