#!/usr/bin/env python3
"""Plot current-control runs with shared comparison axes and numerical summaries."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'gridkit-matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
NAMES = ['GFL', 'GFM']
BLUE, ORANGE, GREY = '#0072B2', '#D55E00', '#555555'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.alpha': .2, 'lines.linewidth': 1.2,
                     'savefig.dpi': 180, 'figure.constrained_layout.use': True})


def read_run(folder, name):
    record = json.loads((folder / f'{name}.run.json').read_text())
    for filename, expected in record['input_sha256'].items():
        if hashlib.sha256((folder / filename).read_bytes()).hexdigest() != expected:
            raise ValueError(f'{folder / filename}: input changed since simulation')
    solver = json.loads((folder / f'{name}.solver.json').read_text())
    case = json.loads((folder / solver['system_model_file']).read_text())
    csv = folder / solver['output_file']
    if csv.exists():
        raw = np.genfromtxt(csv, delimiter=',', names=True)
        data = {key: np.atleast_1d(raw[key]) for key in raw.dtype.names}
        np.savez_compressed(folder / f'{name}.npz', **data)
    else:
        with np.load(folder / f'{name}.npz') as raw:
            data = {key: raw[key] for key in raw.files}
    t = data['t']
    if len(t) < 2 or not all(np.isfinite(x).all() for x in data.values()) or (np.diff(t) < 0).any():
        raise ValueError(f'{folder}/{name}: invalid waveform data')
    if not np.isclose(t[-1], solver['tmax'], rtol=0, atol=1e-10):
        raise ValueError(f'{folder}/{name}: incomplete simulation')
    signals = {s['id']: s['value'] for s in case['signals'] if 'value' in s}
    signals.update(solver.get('signal_values', {}))
    devices = {d['id']: d for d in case['devices']}
    return {'name': name, 'folder': folder, 'data': data, 'solver': solver,
            'signals': signals, 'devices': devices, 'record': record,
            'fc': devices['pwm']['params']['fc'], 'mu': solver['mu'],
            'frequency': devices['grid']['params']['omega'] / (2 * np.pi)}


def integral(t, x, points):
    """Integrate linear segments, retaining both sides of duplicate event times."""
    dt = np.diff(t)
    cumulative = np.r_[0, np.cumsum(.5 * (x[:-1] + x[1:]) * dt)]
    slope = np.divide(np.diff(x), dt, out=np.zeros_like(x[:-1]), where=dt > 0)
    index = np.clip(np.searchsorted(t, points, side='right') - 1, 0, len(t) - 2)
    local = np.asarray(points) - t[index]
    return cumulative[index] + local * x[index] + .5 * local**2 * slope[index]


def mean(t, x, begin, end):
    return float(np.diff(integral(t, x, [begin, end]))[0] / (end - begin))


def carrier_mean(t, x, fc):
    boundaries = np.arange(int(np.floor(t[-1] * fc + 1e-9)) + 1) / fc
    return .5 * (boundaries[:-1] + boundaries[1:]), np.diff(integral(t, x, boundaries)) * fc


def trace(ax, run, key, label, color=BLUE):
    t, x = run['data']['t'], run['data'][key]
    ax.plot(t, x, color=color, alpha=.18, linewidth=.45)
    tm, xm = carrier_mean(t, x, run['fc'])
    ax.plot(tm, xm, color=color, label=label)


def duty(run):
    return .5 * (1 + np.column_stack([run['data'][f'PWM_pwm_m{p}'] for p in 'abc']))


def attenuation(x):
    return np.divide(2 * x * np.exp(-x), -np.expm1(-2 * x),
                     out=np.ones_like(x, dtype=float), where=x != 0)


def train(run, d, phase):
    """Periodic pulse function evaluated at the current duty and carrier phase."""
    fc, mu = run['fc'], run['mu']
    alignment = run['devices']['pwm']['params'].get('alignment', .5)
    on, off = alignment * (1 - d), alignment + (1 - alignment) * d
    tail = np.log(4 / np.finfo(float).eps)
    radius = int(np.ceil(tail * fc / mu))
    decay = 2 * np.pi**2 * fc / mu
    count = max(1, int(np.ceil(tail / decay)))
    # Use the shorter equivalent expansion for sharp edges or broad smoothing.
    if 2 * radius + 1 <= count:
        s = np.zeros_like(d)
        for offset in range(-radius, radius + 1):
            s += .5 * (np.tanh(.5 * mu / fc * (phase - on - offset))
                       - np.tanh(.5 * mu / fc * (phase - off - offset)))
        return s
    s = d.copy()
    center = .5 * (on + off)
    for n in range(1, count + 1):
        s += 2 * d * np.sinc(n * d) * attenuation(np.asarray(n * decay)) * np.cos(2 * np.pi * n * (phase - center))
    return s


def pulse_prediction(run):
    """Evaluate the continuous switching equation from monitored modulation."""
    phase = np.remainder(run['data']['t'] * run['fc'], 1)[:, None]
    return train(run, duty(run), phase)


def harmonics(run, end, prediction):
    """Compare Fourier amplitudes with the continuous PWM reference waveform."""
    fc, fm = run['fc'], run['frequency']
    last = int(np.floor(end * fc + 1e-9))
    count = min(last, int(round(fc / fm)))
    if count < 1:
        raise ValueError('The harmonic plot needs at least one complete carrier period')
    begin, end = (last - count) / fc, last / fc
    data = run['data']
    t, voltage = data['t'], data['Converter_bridge_ea']
    bridge = reference(run, run['devices']['bridge']['inputs']['vdc']) * (prediction[:, 0] - prediction.mean(axis=1))
    frequencies = np.array([fm, fc-2*fm, fc-fm, fc, fc+fm, fc+2*fm,
                            2*fc-fm, 2*fc+fm, 3*fc-2*fm, 3*fc+2*fm])
    measured, predicted = [], []
    for f in frequencies:
        phase = np.exp(-2j * np.pi * f * t)
        measured.append(abs(2 * np.diff(integral(t, voltage * phase, [begin, end]))[0] / (end - begin)))
        predicted.append(abs(2 * np.diff(integral(t, bridge * phase, [begin, end]))[0] / (end - begin)))
    return {'window_s': [begin, end], 'frequencies_hz': frequencies.tolist(),
            'measured_peak_V': measured, 'predicted_peak_V': predicted,
            'maximum_absolute_error_V': float(np.max(np.abs(np.array(measured) - predicted)))}


def reference(run, signal):
    values = np.full_like(run['data']['t'], run['signals'][signal])
    for event in run['solver']['events']:
        if event.get('signal_id') == signal:
            values[run['data']['t'] >= event['time']] = event['value']
    return values


def mu_label(mu):
    value = f'{mu:g}'
    if 'e' in value:
        coefficient, exponent = value.split('e')
        value = ('' if coefficient == '1' else coefficient + r'\times') + f'10^{{{int(exponent)}}}'
    return rf'$\mu={value}\,\mathrm{{s}}^{{-1}}$'


def time_axes(axes, run, end, begin=0):
    for ax in axes:
        ax.ticklabel_format(axis='both', style='plain', useOffset=False)
        ax.set_xlim(begin, end)
        ax.set_xticks(np.linspace(begin, end, 6))
        for event in run['solver']['events']:
            if begin <= event['time'] <= end:
                ax.axvline(event['time'], color=GREY, linewidth=.7, linestyle=':')
        ax.legend(loc='upper right', ncol=3, fontsize=9)
    axes[-1].set_xlabel('Time [s]')


def control_figure(run, end):
    name, d = run['name'], run['data']
    fig, ax = plt.subplots(4, 1, figsize=(10, 9), sharex=True)
    if name == 'GFL':
        trace(ax[0], run, 'Park_grid_current_y1', '$i_d$')
        trace(ax[1], run, 'Park_grid_current_y2', '$i_q$')
        params = run['devices']['power_control']['params']
        ax[0].axhline(params['Pref'] / params['V'], linestyle='--', color=ORANGE, label=r'$i_d^{\mathrm{ref}}$')
        ax[1].axhline(-params['Qref'] / params['V'], linestyle='--', color=ORANGE, label=r'$i_q^{\mathrm{ref}}$')
        ax[0].set_ylabel('d-axis current [A]')
        ax[1].set_ylabel('q-axis current [A]')
    else:
        trace(ax[0], run, 'Park_voltage_y1', '$v_d$')
        trace(ax[1], run, 'Park_voltage_y2', '$v_q$')
        vrefd, vrefq = run['devices']['voltage_control']['inputs']['vref']
        ax[0].plot(d['t'], reference(run, vrefd), '--', color=ORANGE, label=r'$v_d^{\mathrm{ref}}$')
        ax[1].plot(d['t'], reference(run, vrefq), '--', color=ORANGE, label=r'$v_q^{\mathrm{ref}}$')
        ax[0].set_ylabel('d-axis voltage [V]')
        ax[1].set_ylabel('q-axis voltage [V]')
    command = np.hypot(d['InnerCurrentControl_current_control_ud'], d['InnerCurrentControl_current_control_uq'])
    limit = np.sqrt(3 / 8) * run['devices']['pwm']['params']['Mmax'] * reference(run, run['devices']['bridge']['inputs']['vdc'])
    ax[2].plot(d['t'], command, color=BLUE, label=r'$\|\mathbf{u}\|_2$')
    ax[2].plot(d['t'], limit, color=ORANGE, linestyle='--', label='available voltage command')
    ax[2].set_ylabel('Voltage command [V]')
    ax[3].plot(d['t'], d['PLL_pll_omega'] / (2 * np.pi), color=BLUE, label='PLL frequency')
    ax[3].axhline(run['frequency'], color=ORANGE, linestyle='--', label='grid frequency')
    ax[3].set_ylabel('Frequency [Hz]')
    ax[0].set_title(f'{name}: {mu_label(run["mu"])}')
    time_axes(ax, run, end)
    ax[2].legend(loc='lower right', ncol=2, fontsize=9)
    fig.supxlabel('Feedback traces: faint instantaneous values; solid carrier-period means.', fontsize=9)
    return fig


def switching_figure(run, prediction, end):
    begin = max(0, end - 6 / run['fc'])
    mask = (run['data']['t'] >= begin) & (run['data']['t'] <= end)
    d = {key: value[mask] for key, value in run['data'].items()}
    fig, ax = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    ax[0].plot(d['t'], d['PWM_pwm_sa'], color=BLUE, label='$s_a$')
    ax[0].plot(d['t'], prediction[mask, 0], '--', color=ORANGE, label='continuous PWM reference')
    ax[0].set_ylabel('Switching function [−]')
    ax[0].set_title(f'{run["fc"]:g} Hz PWM, {run["name"]}, {mu_label(run["mu"])}')
    ax[1].plot(d['t'], d['Converter_bridge_ea'], color=BLUE, label=r'bridge $e_a$')
    ax[1].plot(d['t'], d['Filter_filter_voa'], color=ORANGE, label=r'capacitor $v_{\mathrm{o},a}$')
    ax[1].set_ylabel('Phase voltage [V]')
    ax[2].plot(d['t'], d['Filter_filter_ia'], color=BLUE, label='inverter-side $i_a$')
    ax[2].plot(d['t'], d['Filter_filter_iga'], color=ORANGE, label='grid-side $i_{g,a}$')
    ax[2].set_ylabel('Phase current [A]')
    time_axes(ax, run, end, begin)
    return fig


def summarize(run, prediction, spectrum):
    d, t = run['data'], run['data']['t']
    s = np.column_stack([d[f'PWM_pwm_s{p}'] for p in 'abc'])
    voltage = np.column_stack([d[f'Converter_bridge_e{p}'] for p in 'abc'])
    vdc = reference(run, run['devices']['bridge']['inputs']['vdc'])
    result = {'final_time_s': float(t[-1]), 'monitor_samples': len(t),
              'carrier_hz': run['fc'], 'mu': run['mu'],
              'mean_pll_frequency_Hz': mean(t, d['PLL_pll_omega'], *spectrum['window_s']) / (2 * np.pi),
              'maximum_pll_frequency_error_Hz': float(np.max(np.abs(d['PLL_pll_omega'] / (2 * np.pi) - run['frequency']))),
              'dc_voltage_V': float(vdc[0]),
              'edge_10_90_s': 2 * np.log(9) / run['mu'],
              'bridge_identity_max_error_V': float(np.max(np.abs(voltage - vdc[:, None] * (s - s.mean(axis=1)[:, None])))),
              'pulse_prediction_max_error': float(np.max(np.abs(prediction - s))),
              'maximum_reference_norm_A': float(np.max(np.hypot(d['InnerCurrentControl_current_control_ilimd'], d['InnerCurrentControl_current_control_ilimq']))),
              'maximum_voltage_command_V': float(np.max(np.hypot(d['InnerCurrentControl_current_control_ud'], d['InnerCurrentControl_current_control_uq']))),
              'tracking_window_s': spectrum['window_s'], 'harmonics': spectrum}
    for key, column in [('id_A', 'Park_current_y1'), ('iq_A', 'Park_current_y2'),
                        ('vd_V', 'Park_voltage_y1'), ('vq_V', 'Park_voltage_y2')]:
        result[f'mean_{key}'] = mean(t, d[column], *spectrum['window_s'])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('simulation'), help='Run directory, relative to this example')
    parser.add_argument('--compare', type=Path, nargs='+', default=[], help='Additional run directories with shared plot axes')
    parser.add_argument('--tmax', type=float, help='Plot end time; defaults to the common available interval')
    args = parser.parse_args()
    folders = [(HERE / path).resolve() for path in [args.output, *args.compare]]
    if len(set(folders)) != len(folders):
        parser.error('Run directories must be distinct')
    runs = {folder: {name: read_run(folder, name) for name in NAMES} for folder in folders}
    available = min(run['data']['t'][-1] for studies in runs.values() for run in studies.values())
    end = available if args.tmax is None else args.tmax
    if not np.isfinite(end) or end <= 0 or end > available + 1e-10:
        parser.error('tmax must be positive and within every simulation')
    figures = {name: [] for name in [*NAMES, 'switching', 'harmonics']}
    summaries = {}
    for folder, studies in runs.items():
        summary = {'time_window_s': [0, end], 'studies': {}}
        fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
        for ax, (name, run) in zip(axes, studies.items()):
            mask = run['data']['t'] <= end
            run['data'] = {key: value[mask] for key, value in run['data'].items()}
            prediction = pulse_prediction(run)
            spectrum = harmonics(run, end, prediction)
            summary['studies'][name] = summarize(run, prediction, spectrum)
            figures[name].append((folder, control_figure(run, end)))
            if name == 'GFM':
                figures['switching'].append((folder, switching_figure(run, prediction, end)))
            x = np.arange(len(spectrum['frequencies_hz']))
            ax.bar(x - .18, spectrum['measured_peak_V'], .36, color=BLUE, label='simulated bridge voltage')
            ax.bar(x + .18, spectrum['predicted_peak_V'], .36, color=ORANGE, label='continuous PWM reference')
            ax.set_yscale('log')
            ax.set_ylim(bottom=1e-3)
            ax.set_xlim(-.6, len(x) - .4)
            ax.set_ylabel('Peak amplitude [V]')
            ax.set_title(f'{name}: {mu_label(run["mu"])}')
            ax.legend(loc='upper right')
        axes[-1].set_xticks(x, [f'{f:g}' for f in spectrum['frequencies_hz']])
        axes[-1].set_xlabel('Frequency [Hz]')
        begin, finish = spectrum['window_s']
        fig.supxlabel(f'Fourier amplitudes over {begin:.5f}–{finish:.5f} s; the window includes transients.', fontsize=9)
        figures['harmonics'].append((folder, fig))
        summaries[folder] = summary
    limits = {}
    for name, group in figures.items():
        limits[name] = []
        for axes in zip(*(fig.axes for _, fig in group)):
            bounds = (min(ax.get_ylim()[0] for ax in axes), max(ax.get_ylim()[1] for ax in axes))
            for ax in axes:
                ax.set_ylim(*bounds)
            limits[name].append({'x': list(axes[0].get_xlim()), 'y': list(bounds)})
            if any(ax.get_xlim() != axes[0].get_xlim() for ax in axes):
                raise ValueError(f'{name}: incompatible comparison axes')
        for folder, fig in group:
            for extension in ['png', 'pdf']:
                fig.savefig(folder / f'{name}.{extension}')
            plt.close(fig)
    for folder, summary in summaries.items():
        summary['axis_limits'] = limits
        (folder / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(f'{folder}: four PNG/PDF figures and summary.json')


if __name__ == '__main__':
    main()
