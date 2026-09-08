"""Run the transformer scenarios, check the steady phasors and the inrush; optionally plot."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile

from steady_state import samples, solve

ROOT = Path(__file__).resolve().parent

# Measured floor at the 1e-8 IDA tolerances: 1.4e-7 relative to each
# monitor's phasor amplitude over the 0.1 s run.
STEADY_TOLERANCE = 2e-7
# The breaker closes at the phase-a voltage zero crossing; the measured flux
# peak is 1.98 per unit and the inrush peak is 4.9 times the steady peak.
INRUSH_RATIO = 4.0


def bank_parameters():
    case = json.loads((ROOT / 'Transformer.case.json').read_text())
    return {device['id']: device for device in case['devices']}['bank']['params']


def run(executable, scenario, output):
    output.mkdir(parents=True, exist_ok=True)
    result = subprocess.run([str(executable), str(ROOT / f'{scenario}.solver.json')],
                            cwd=output, text=True, capture_output=True, check=False)
    (output / 'simulation.log').write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    with (output / 'mon.csv').open(newline='') as stream:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    layout = json.loads((output / 'state.csv.json').read_text())['variables']
    for row in rows:
        assert all(math.isfinite(value) for value in row.values())
    return rows, layout


def expected_monitors(phasors, time):
    """Steady monitor values of the bank at one instant."""
    values = {'i12a': samples(phasors['i12'], time, phasors['omega'])[0]}
    for phase, value in zip('abc', samples(phasors['psi1'], time, phasors['omega'])):
        values['psi1' + phase] = value
    for phase, value in zip('abc', samples(phasors['iw1'], time, phasors['omega'])):
        values['i1' + phase] = -phasors['i_peak1'] * value
    values['i2a'] = -phasors['i_peak2'] * samples(phasors['iw2'], time, phasors['omega'])[0]
    return values


def scales(phasors):
    return {'i12a': abs(phasors['i12']),
            'psi1a': abs(phasors['psi1']), 'psi1b': abs(phasors['psi1']), 'psi1c': abs(phasors['psi1']),
            'i1a': phasors['i_peak1'] * abs(phasors['iw1']),
            'i1b': phasors['i_peak1'] * abs(phasors['iw1']),
            'i1c': phasors['i_peak1'] * abs(phasors['iw1']),
            'i2a': phasors['i_peak2'] * abs(phasors['iw2'])}


def validate_steady(executable, output, phasors):
    rows, layout = run(executable, 'Steady', output)
    differential = [item['component'] for item in layout if item['differential']]
    assert differential == ['bank'] * 9, differential
    solver = json.loads((ROOT / 'Steady.solver.json').read_text())
    assert len(rows) == round(solver['tmax'] / solver['dt_monitor']) + 1, len(rows)
    errors = {}
    for row in rows:
        for name, value in expected_monitors(phasors, row['t']).items():
            errors[name] = max(errors.get(name, 0.0), abs(row['Transformer_bank_' + name] - value))
    relative = {name: error / scale for (name, error), scale in zip(errors.items(), scales(phasors).values())}
    assert max(relative.values()) < STEADY_TOLERANCE, relative
    return rows, relative


def validate_energization(executable, output, phasors):
    rows, layout = run(executable, 'Energization', output)
    solver = json.loads((ROOT / 'Energization.solver.json').read_text())
    event = solver['events'][0]['time']
    before = [row for row in rows if row['t'] < event - 1e-12]
    after = [row for row in rows if row['t'] > event + 1e-12]
    at_event = [row for row in rows if abs(row['t'] - event) < 1e-12]
    assert before and after and len(at_event) == 2, (len(before), len(after), len(at_event))
    assert at_event[0]['Switch_breaker_open'] == 1.0 and at_event[1]['Switch_breaker_open'] == 0.0
    for row in before:
        assert row['Switch_breaker_open'] == 1.0
        assert abs(row['Transformer_bank_i1a']) < 1e-9
        assert abs(row['Transformer_bank_psi1a']) < 1e-12
    steady_peak = phasors['i_peak1'] * abs(phasors['iw1'])
    peak_flux = max(abs(row['Transformer_bank_psi1a']) for row in after)
    peak_current = max(abs(row['Transformer_bank_i1a']) for row in after)
    knee = bank_parameters()['knee']
    assert peak_flux > knee, peak_flux
    assert peak_current > INRUSH_RATIO * steady_peak, (peak_current, steady_peak)
    return rows, {'peak_flux': peak_flux, 'peak_current': peak_current, 'steady_peak': steady_peak}


def plot(output, steady, energization, phasors, knee):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    time = [row['t'] * 1e3 for row in steady]
    figure, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    for name, label in (('i1a', 'HV phase a'), ('i2a', 'LV phase a')):
        axes[0].plot(time, [row['Transformer_bank_' + name] for row in steady], label=label)
        axes[0].plot(time, [expected_monitors(phasors, row['t'])[name] for row in steady], 'k:', lw=0.8)
    axes[0].set_ylabel('A')
    axes[0].set_title('Steady start from the phasor state file (dotted: phasor solution)')
    axes[0].legend(loc='upper right')
    axes[1].plot(time, [row['Transformer_bank_psi1a'] for row in steady], label='psi1a')
    axes[1].plot(time, [row['Transformer_bank_i12a'] for row in steady], label='i12a')
    axes[1].set_ylabel('p.u.')
    axes[1].set_xlabel('ms')
    axes[1].legend(loc='upper right')
    for axis in axes:
        axis.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(output / 'steady.png', dpi=150)
    plt.close(figure)

    time = [row['t'] * 1e3 for row in energization]
    figure, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    for phase in 'abc':
        axes[0].plot(time, [row['Transformer_bank_i1' + phase] for row in energization], label='i1' + phase)
    axes[0].set_ylabel('A')
    axes[0].set_title('Energization at the phase-a voltage zero crossing')
    axes[0].legend(loc='upper right')
    for phase in 'abc':
        axes[1].plot(time, [row['Transformer_bank_psi1' + phase] for row in energization], label='psi1' + phase)
    for level in (knee, -knee):
        axes[1].axhline(level, color='k', lw=0.6, ls=':')
    axes[1].set_ylabel('p.u.')
    axes[1].legend(loc='upper right')
    axes[2].plot(time, [row['Transformer_bank_i2a'] for row in energization], label='i2a')
    axes[2].plot(time, [row['LoadZ_load_ia'] for row in energization], '--', label='load ia')
    axes[2].set_ylabel('A')
    axes[2].set_xlabel('ms')
    axes[2].legend(loc='upper right')
    for axis in axes:
        axis.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(output / 'energization.png', dpi=150)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', required=True, type=Path, help='EMTDynamicSimulation executable')
    parser.add_argument('--results', type=Path, help='directory for outputs; default is temporary')
    parser.add_argument('--plot', action='store_true', help='write steady.png and energization.png')
    arguments = parser.parse_args()
    case = json.loads((ROOT / 'Transformer.case.json').read_text())
    phasors = solve(case)
    knee = bank_parameters()['knee']

    executable = arguments.exe.resolve()

    def validate(output):
        steady, relative = validate_steady(executable, output / 'steady', phasors)
        energization, metrics = validate_energization(executable, output / 'energization', phasors)
        metrics['steady_relative_error'] = relative
        (output / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
        if arguments.plot:
            plot(output, steady, energization, phasors, knee)
        print(json.dumps(metrics, indent=2))

    if arguments.results:
        validate(arguments.results.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix='gridkit-emt-transformer-') as directory:
            validate(Path(directory))


if __name__ == '__main__':
    main()
