#!/usr/bin/env python3
"""Plot the ten-bus REGFMA response and check its recorded terminal quantities."""
import csv
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/gridkit-ibr-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
COLORS = ['#176b9c', '#e18727', '#2a9462', '#ae487b']
CLARKE = np.sqrt(2 / 3) * np.array([[1., -.5, -.5], [0., np.sqrt(3) / 2, -np.sqrt(3) / 2]])
plt.rcParams.update({'axes.prop_cycle': matplotlib.cycler(color=COLORS), 'axes.grid': True,
                     'grid.alpha': .22, 'font.size': 10, 'figure.dpi': 110, 'savefig.dpi': 145,
                     'axes.spines.top': False, 'axes.spines.right': False})


def read_csv(path):
    with path.open(newline='') as stream:
        header = next(csv.reader(stream))
    values = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    if values.shape[1] != len(header) or not np.isfinite(values).all():
        raise ValueError(f'Incomplete or nonfinite data: {path}')
    return {name: values[:, index] for index, name in enumerate(header)}


def phases(data, prefix, stem):
    return np.column_stack([data[f'{prefix}_{stem}{phase}'] for phase in 'abc'])


def main():
    directory = HERE / 'results'
    run = json.loads((directory / 'run.json').read_text())
    for name, digest in run['input_sha256'].items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Run input changed: {name}')
    study = json.loads((directory / 'FaultClearing.solver.json').read_text())
    case = json.loads((directory / study['system_model_file']).read_text())
    devices = {device['id']: device for device in case['devices']}
    data = read_csv(directory / study['output_file'])
    steps = read_csv(directory / study['step_output_file'])
    t = data['t']
    events = study['events']
    expected = round(study['tmax'] / study['dt_monitor']) + 1 + len(events)
    if (len(t) != expected or t[0] != 0 or abs(t[-1] - study['tmax']) > 1e-12
            or np.any(np.diff(t) < 0)):
        raise ValueError('Incomplete monitor timeline')
    for event in events:
        rows = np.flatnonzero(np.isclose(t, event['time'], rtol=0, atol=1e-12))
        command = float(event['open'])
        if len(rows) != 2 or not np.array_equal(data['Switch_fault_open'][rows], [1 - command, command]):
            raise ValueError(f'Incorrect event pair at {event["time"]}')
    if len(steps['time']) != run['ida']['steps'] or np.any(steps['step'] <= 0):
        raise ValueError('Accepted-step record does not match solver statistics')

    sources = {}
    checks = {'max_kcl_error_A': 0., 'max_zero_sequence_current_A': 0.,
              'max_active_power_error_W': 0., 'max_reactive_power_error_var': 0.}
    summary = {}
    for bus, remote in ((4, 7), (5, 8), (6, 8)):
        prefix = f'REGFMA_regfma_{bus}'
        params = devices[f'regfma_{bus}']['params']
        voltage = phases(data, f'Bus_bus_{bus}', 'v')
        current = phases(data, prefix, 'i')
        vab, iab = voltage @ CLARKE.T, current @ CLARKE.T
        magnitude = np.linalg.norm(iab, axis=1) / (params['S'] / params['V'])
        frequency = data[f'{prefix}_omega'] / (2 * np.pi)
        vpu = np.linalg.norm(vab, axis=1) / params['V']
        p = np.sum(voltage * current, axis=1)
        q = vab[:, 1] * iab[:, 0] - vab[:, 0] * iab[:, 1]
        kcl = (current + phases(data, f'LoadZ_load_{bus}', 'i')
               - phases(data, f'LineLumped_line_{bus}_{remote}', 'i12')
               - phases(data, f'Bus_bus_{bus}', 'i_sh'))
        for key, error in (
                ('max_kcl_error_A', np.max(np.abs(kcl))),
                ('max_zero_sequence_current_A', np.max(np.abs(current.sum(axis=1) / 3))),
                ('max_active_power_error_W', np.max(np.abs(p - data[f'{prefix}_p']))),
                ('max_reactive_power_error_var', np.max(np.abs(q - data[f'{prefix}_q'])))):
            checks[key] = max(checks[key], float(error))
        if magnitude.max() > params['ImaxF'] + 1e-5:
            raise ValueError(f'Bus {bus}: current exceeds ImaxF')
        sources[bus] = {'v': vpu, 'i': magnitude, 'f': frequency, 'p': p / 1e6, 'q': q / 1e6,
                        'pf': data[f'{prefix}_pf'] * params['S'] / 1e6,
                        'qf': data[f'{prefix}_qf'] * params['S'] / 1e6,
                        'e': data[f'{prefix}_edroop'], 'vabc': voltage, 'iabc': current}
        summary[f'regfma_{bus}'] = {
            'minimum_voltage_pu': float(vpu.min()), 'maximum_current_pu': float(magnitude.max()),
            'minimum_frequency_Hz': float(frequency.min()), 'maximum_frequency_Hz': float(frequency.max()),
            'final_frequency_Hz': float(frequency[-1]), 'final_voltage_pu': float(vpu[-1]),
            'final_active_power_MW': float(p[-1] / 1e6), 'final_reactive_power_Mvar': float(q[-1] / 1e6)}
    if (checks['max_kcl_error_A'] > .01 or checks['max_zero_sequence_current_A'] > .001
            or checks['max_active_power_error_W'] > .1 or checks['max_reactive_power_error_var'] > .1):
        raise ValueError(f'Terminal consistency check failed: {checks}')
    boundaries = [0.] + [event['time'] for event in events] + [study['tmax']]
    intervals = []
    for start, stop in zip(boundaries[:-1], boundaries[1:]):
        mask = (steps['time'] > start) & (steps['time'] <= stop)
        intervals.append({'start_s': start, 'stop_s': stop, 'accepted_steps': int(mask.sum()),
                          'median_step_s': float(np.median(steps['step'][mask]))})
    summary = {'run': run, 'monitor_samples': len(t), 'sources': summary, 'checks': checks,
               'step_intervals': intervals,
               'minimum_accepted_step_s': float(steps['step'].min()),
               'maximum_accepted_step_s': float(steps['step'].max())}
    (directory / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')

    plots = directory / 'plots'
    plots.mkdir(exist_ok=True)
    fault_start, fault_end = (event['time'] for event in events)

    def finish(fig, axes, name, limits):
        for ax in axes.flat:
            ax.axvspan(fault_start, fault_end, color='#555555', alpha=.10, linewidth=0)
            for event in events:
                ax.axvline(event['time'], color='#555555', lw=.7, ls=':')
            ax.set_xlim(*limits)
            ax.legend(loc='best', fontsize=8, ncol=2)
        for ax in axes[-1].flat:
            ax.set_xlabel('Time [s]')
        for extension in ('png', 'pdf'):
            fig.savefig(plots / f'{name}.{extension}')
        plt.close(fig)

    fig, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True, constrained_layout=True)
    fig.suptitle('Ten-bus REGFMA fault response\nThree-phase 2 Ω/phase shunt, applied at 1.00 s and cleared at 1.06 s')
    for color, (bus, source) in zip(COLORS, sources.items()):
        label = f'REGFMA at bus {bus}'
        for ax, key, ylabel in ((axes[0, 0], 'v', 'Terminal voltage [p.u.]'),
                                (axes[0, 1], 'f', 'Internal frequency [Hz]'),
                                (axes[1, 0], 'pf', 'Filtered active power [MW]'),
                                (axes[1, 1], 'qf', 'Filtered reactive power [Mvar]'),
                                (axes[2, 0], 'i', 'Terminal current [p.u.]'),
                                (axes[2, 1], 'e', 'Internal voltage command [p.u.]')):
            ax.plot(t, source[key], color=color, ls={4: '-', 5: '--', 6: '-.'}[bus], lw=1.2, label=label)
            ax.set_ylabel(ylabel)
    v7 = np.linalg.norm(phases(data, 'Bus_bus_7', 'v') @ CLARKE.T, axis=1) / 13800
    axes[0, 0].plot(t, v7, color=COLORS[3], lw=1, ls='--', label='Fault connection, bus 7')
    for machine, color in zip((1, 2, 3), COLORS):
        axes[0, 1].plot(t, data[f'Machine_machine_{machine}_omega'] * 60,
                        color=color, lw=.8, ls='--', label=f'Machine {machine}')
    axes[2, 0].axhline(devices['regfma_4']['params']['ImaxF'], color='#555555', ls='--', lw=1, label='Current limit')
    finish(fig, axes, 'response', (0, study['tmax']))

    fig, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True, constrained_layout=True)
    fig.suptitle('Ten-bus REGFMA fault detail\nRecorded terminal quantities; no waveform averaging')
    for color, (bus, source) in zip(COLORS, sources.items()):
        for ax, key, ylabel in ((axes[0, 0], 'v', 'Terminal voltage [p.u.]'),
                                (axes[0, 1], 'i', 'Terminal current [p.u.]'),
                                (axes[1, 0], 'p', 'Active power [MW]'),
                                (axes[1, 1], 'q', 'Reactive power [Mvar]')):
            ax.plot(t, source[key], color=color, ls={4: '-', 5: '--', 6: '-.'}[bus], lw=1.2,
                    label=f'REGFMA at bus {bus}')
            ax.set_ylabel(ylabel)
    for index, phase in enumerate('abc'):
        axes[2, 0].plot(t, sources[4]['vabc'][:, index] / 1000, lw=1, label=f'Phase {phase}')
        axes[2, 1].plot(t, sources[4]['iabc'][:, index], lw=1, label=f'Phase {phase}')
    axes[2, 0].set_ylabel('Bus 4 phase voltage [kV]')
    axes[2, 1].set_ylabel('REGFMA at bus 4 current [A]')
    axes[0, 1].axhline(devices['regfma_4']['params']['ImaxF'], color='#555555', ls='--', lw=1, label='Current limit')
    finish(fig, axes, 'fault', (fault_start - .05, fault_end + .14))

    (directory / 'index.html').write_text(f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><title>Ten-bus REGFMA response</title>
<style>body{{font:16px system-ui;max-width:1200px;margin:32px auto;padding:0 20px;color:#222}}
img{{width:100%;height:auto}}a{{color:#176b9c}}</style>
<h1>Ten-bus REGFMA response</h1>
<p>Three 5 MVA REGFMA sources and three synchronous machines. Three-phase 2 Ω/phase
fault from 1.00 to 1.06 s; 3 s simulation, μ = {study['mu']:g}.</p>
<p>{run['wall_seconds']:.3f} s wall time; {run['simulation_cpu_seconds']:.3f} s simulation CPU;
{run['ida']['steps']:,} accepted steps; {len(t):,} monitor samples.</p>
<p>Terminal consistency and event checks passed. These are simulation results,
not external model validation. See <a href="../README.md">study details</a>.</p>
<p><a href="summary.json">Metrics and checks</a> · <a href="FaultClearing.csv">Monitor CSV</a> ·
<a href="steps.csv">Accepted steps</a> · <a href="run.log">Solver log</a></p>
<h2>Response</h2><a href="plots/response.pdf">Vector PDF</a>
<img src="plots/response.png" alt="Voltage, frequency, power, current and voltage-command response">
<h2>Fault detail</h2><a href="plots/fault.pdf">Vector PDF</a>
<img src="plots/fault.png" alt="Fault response and three-phase waveforms at bus 4">
</html>\n''')
    print(f'Wrote {directory / "index.html"}; terminal consistency and event checks passed.')


if __name__ == '__main__':
    main()
