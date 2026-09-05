#!/usr/bin/env python3
"""Check rational network frequency response and the hybrid EMT trajectory."""

import argparse
import cmath
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile


HERE = Path(__file__).resolve().parent


def run(exe, directory, case, study):
    directory.mkdir(parents=True, exist_ok=True)
    case_path = directory / 'case.json'
    case_path.write_text(json.dumps(case) + '\n')
    study = dict(study, system_model_file=str(case_path), output_file='mon.csv')
    study_path = directory / 'solver.json'
    study_path.write_text(json.dumps(study) + '\n')
    with (directory / 'run.log').open('w') as log:
        result = subprocess.run([str(exe), str(study_path)], cwd=directory,
                                stdout=log, stderr=subprocess.STDOUT, timeout=150)
    if result.returncode:
        raise RuntimeError((directory / 'run.log').read_text())
    with (directory / 'mon.csv').open(newline='') as stream:
        rows = [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(stream)]
    if not rows or not all(math.isfinite(x) for row in rows for x in row.values()):
        raise AssertionError('Missing or nonfinite monitor output')
    if not math.isclose(rows[-1]['t'], study['tmax'], abs_tol=1e-10):
        raise AssertionError('Simulation did not reach its final time')
    if any(b['t'] <= a['t'] for a, b in zip(rows, rows[1:])):
        raise AssertionError('Monitor times are not increasing')
    expected = round(study['tmax'] / study['dt_monitor']) + 1
    if len(rows) != expected:
        raise AssertionError(f'Missing monitor samples: {len(rows)} != {expected}')
    return rows


def circuit(exe, results, frequency):
    case = json.loads((HERE / 'Circuit.case.json').read_text())
    omega = 2 * math.pi * frequency
    for device in case['devices']:
        if device['id'] == 'source':
            device['params']['omega'] = omega
    periods = max(12, math.ceil(.3 * frequency))
    study = dict(dt_monitor=1 / (96 * frequency), tmax=periods / frequency,
                 rel_tol=1e-8, abs_tol=1e-10, max_steps=1000000)
    rows = run(exe, results / f'{frequency}Hz', case, study)

    # Independent passive circuit relationships, with physical R/L/C values.
    # The oracle does not use the pole/residue evaluation in GridKit.
    s = 1j * omega
    ys = 1 / (10 + .01 * s) + 1 / (160 + .04 * s)
    zl = 1 + .004 * s + 1 / (1 / 3 + 1 / (.0015 * s))
    yh = (1e-5 + 2e-6 * s + 1 / (20000 + (40 / 3) * s)) / 2
    zz = 200 + .03 * s + 1 / (1 / 30 + 1 / (.06 * s))
    yl, yz = 1 / zl, 1 / zz
    a, b = ys + yh + yl, yz + yh + yl
    e = 13800 * math.sqrt(2 / 3)
    v1, v2 = ys * e * b / (a * b - yl * yl), ys * e * yl / (a * b - yl * yl)
    phasors = {'Bus_grid_v': v1, 'Bus_terminal_v': v2,
               'VoltageSource_source_i': ys * (e - v1),
               'LineLumped_line_i12': yl * (v1 - v2),
               'LineLumped_line_i_sh1': -yh * v1,
               'LineLumped_line_i_sh2': -yh * v2,
               'LoadZ_load_i': -yz * v2}
    error = 0.
    for row in rows[-4 * 96:]:
        for prefix, phasor in phasors.items():
            for phase, angle in zip('abc', (0, -2 * math.pi / 3, 2 * math.pi / 3)):
                expected = (phasor * cmath.exp(1j * (omega * row['t'] + angle))).real
                error = max(error, abs(row[prefix + phase] - expected) / max(1., abs(phasor)))
    # This includes integrator, CSV, and residual startup error, normalized by
    # each waveform's own peak. Four settled cycles check phase as well as gain.
    if error > 2e-4:
        raise AssertionError(f'{frequency} Hz circuit error {error:.3e} exceeds 2e-4')
    print(f'{frequency} Hz: maximum normalized waveform error {error:.3e}', flush=True)
    return error


def hybrid(exe, results):
    case = json.loads((HERE / 'Hybrid.case.json').read_text())
    study = json.loads((HERE / 'Hybrid.solver.json').read_text())
    rows = run(exe, results / 'hybrid', case, study)
    tight = dict(study, rel_tol=study['rel_tol'] / 10, abs_tol=study['abs_tol'] / 10)
    reference = run(exe, results / 'hybrid-tight', case, tight)
    if len(rows) != len(reference):
        raise AssertionError('Hybrid runs have different sample counts')
    scales = {key: max(1., max(abs(row[key]) for row in reference)) for key in reference[0]}
    error, kcl = 0., 0.
    for row, ref in zip(rows, reference):
        if not math.isclose(row['t'], ref['t'], abs_tol=1e-12):
            raise AssertionError('Hybrid sample times do not match')
        for key in ref:
            error = max(error, abs(row[key] - ref[key]) / scales[key])
        if abs(row['Machine_machine_efd'] - row['Ieeet1_exciter_efd']) > 1e-6:
            raise AssertionError('Machine field input differs from the IEEET1 output')
        for phase in 'abc':
            grid = (row['VoltageSource_source_i' + phase]
                    + row['LineLumped_line_i_sh1' + phase]
                    - row['LineLumped_line_i12' + phase])
            terminal = (row['Machine_machine_i' + phase]
                        + row['DependentVoltageSource_filter_i' + phase]
                        + row['LoadZ_load_i' + phase]
                        + row['LineLumped_line_i_sh2' + phase]
                        + row['LineLumped_line_i12' + phase])
            kcl = max(kcl, abs(grid), abs(terminal))
    if error > 2e-3 or kcl > .05:
        raise AssertionError(f'Hybrid convergence error={error:.3e}, KCL={kcl:.3e} A')
    if not all(.8 < row['Machine_machine_omega'] < 1.2 for row in rows):
        raise AssertionError('Hybrid machine speed left its operating range')
    print(f'Hybrid: tolerance refinement error {error:.3e}; maximum KCL error {kcl:.3e} A', flush=True)
    return dict(refinement=error, kcl_amps=kcl)


def validate(exe, results):
    summary = {'circuit': {f: circuit(exe, results, f) for f in (30, 60, 180)},
               'hybrid': hybrid(exe, results)}
    (results / 'validation.json').write_text(json.dumps(summary, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--results', type=Path, help='Keep logs and outputs in this directory')
    args = parser.parse_args()
    if args.results:
        validate(args.exe.resolve(), args.results.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='gridkit-frequency-') as directory:
            validate(args.exe.resolve(), Path(directory))


if __name__ == '__main__':
    main()
