#!/usr/bin/env python3
"""Plot actual EMT monitor/state CSVs and independently check network identities."""
import argparse
import csv
import html
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/gridkit-ibr-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle, Rectangle
import numpy as np

HERE = Path(__file__).resolve().parent

TITLES = {'01_Baseline': 'Baseline and startup', '02_LoadStep': '2 MW nominal load connection',
          '03_LoadPulse': 'Load connection and shedding', '04_FaultClearing': 'Three-phase resistive fault and clearing',
          '05_TieReclosing': 'Tie opening and reclosing',
          '06_FaultClearingSwitching': 'Fault clearing with resolved PWM switching'}
COLORS = ['#176b9c', '#e18727', '#2a9462', '#ae487b', '#6c60a8', '#697c32', '#b44137', '#218f9b', '#7e654d', '#454545']
plt.rcParams.update({'axes.prop_cycle': matplotlib.cycler(color=COLORS), 'axes.grid': True,
                     'grid.alpha': .22, 'font.size': 10, 'figure.dpi': 110, 'savefig.dpi': 145,
                     'axes.spines.top': False, 'axes.spines.right': False})


def read_csv(path):
    with path.open(newline='') as stream:
        header = next(csv.reader(stream))
    values = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    return values[:, 0], {key: values[:, i] for i, key in enumerate(header)}


def phase(data, prefix, stem):
    return np.column_stack([data[f'{prefix}_{stem}{p}'] for p in 'abc'])


def cycle_mean(t, value):
    """Trailing 1/60 s trapezoidal mean; startup windows are intentionally NaN."""
    cumulative = np.r_[0., np.cumsum(np.diff(t) * (value[:-1] + value[1:]) / 2)]
    result = 60 * (cumulative - np.interp(t - 1 / 60, t, cumulative))
    result[t < 1 / 60] = np.nan
    return result


def spectrum(signal, dt):
    frequency = np.fft.rfftfreq(len(signal), dt)
    amplitude = 2*np.abs(np.fft.rfft(signal))/len(signal)
    amplitude[0] *= .5
    if len(signal) % 2 == 0:
        amplitude[-1] *= .5
    return frequency, amplitude


def pwm_reference(t, pwm, mu, vdc):
    """Independent sigmoid sum with the current duty in every pulse replica."""
    phase = np.array([0, -2*np.pi/3, 2*np.pi/3])
    duty = (1 + pwm['M'] * np.sin(2*np.pi*pwm['fm']*t[:, None] + phase)) / 2
    on = pwm['alignment'] * (1 - duty)
    off = on + duty
    carrier = np.remainder(t * pwm['fc'], 1)[:, None]
    radius = int(np.ceil(np.log(4 / np.finfo(float).eps) * pwm['fc'] / mu))
    gate = np.zeros_like(duty)
    for k in range(-radius, radius + 1):
        gate += .5 * (np.tanh(.5 * mu / pwm['fc'] * (carrier - on - k))
                      - np.tanh(.5 * mu / pwm['fc'] * (carrier - off - k)))
    bridge = vdc * (gate - gate.mean(axis=1, keepdims=True))
    return gate[:, 0], bridge[:, 0]


def switching_detail(directory):
    """Inspect the switching transient and check continuous-PWM harmonics."""
    study = json.loads((directory / 'run.solver.json').read_text())
    case = json.loads(Path(study['system_model_file']).read_text())
    t, data = read_csv(directory / study['output_file'])
    dt = study['dt_monitor']
    mu = study.get('mu', 240.)
    mask = (t >= study['tmax'] - .5) & (t < study['tmax'])
    if (np.count_nonzero(mask) != round(.5 / dt)
            or not np.allclose(np.diff(t[mask]), dt, rtol=1e-7, atol=1e-12)):
        raise ValueError('Switching spectra require a complete uniform final 0.5 s window')
    constants = {s['id']: s['value'] for s in case['signals'] if 'value' in s}
    constants.update(study.get('signal_values', {}))
    devices = {d['id']: d for d in case['devices']}
    kcl_errors = []
    for b in (4, 5, 6):
        kcl = (phase(data, f'DependentVoltageSource_filter_{b}', 'i')
               - phase(data, f'Bus_bus_{b}', 'i_sh'))
        for device in case['devices']:
            if device['class'] == 'LoadZ' and device['inputs']['bus'] == f'bus_{b}':
                kcl += phase(data, 'LoadZ_' + device['id'], 'i')
            if device['class'] == 'LineLumped':
                for port, sign in (('bus1', -1), ('bus2', 1)):
                    if device['inputs'][port] == f'bus_{b}':
                        kcl += sign * phase(data, 'LineLumped_' + device['id'], 'i12')
        kcl_errors.append(float(np.max(np.abs(kcl))))
    checks = []
    fig, axs = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for b in (4, 5, 6):
        pwm = devices[f'pwm_{b}']['params']
        prediction = pwm_reference(t[mask], pwm, mu, constants[f'dc_{b}'])
        for ax, key, reference, atol in zip(
                axs, (f'PWM_pwm_{b}_sa', f'Converter_converter_{b}_voa'),
                prediction, (1e-8, 1e-3)):
            frequency, amplitude = spectrum(data[key][mask], dt)
            _, predicted = spectrum(reference, dt)
            indices = []
            for harmonic in range(1, 50):
                target = harmonic * pwm['fm']
                index = round(target * len(t[mask]) * dt)
                if index >= len(frequency) or not np.isclose(frequency[index], target):
                    raise ValueError('Predicted harmonic is not a resolved coherent FFT bin')
                indices.append(index)
                measured, expected = float(amplitude[index]), float(predicted[index])
                checks.append(dict(signal=key, frequency_hz=target,
                                   measured_peak=measured, predicted_peak=expected,
                                   absolute_error=abs(measured - expected),
                                   passed=abs(measured - expected) <= atol + 1e-5 * expected))
            if b == 4:
                ax.semilogy(frequency, np.maximum(amplitude, atol / 10), label='Simulation')
                visible = [index for index in indices if predicted[index] > atol]
                ax.scatter(frequency[visible], predicted[visible], facecolors='none',
                           edgecolors=COLORS[1], label='Continuous PWM reference', zorder=3)
                ax.set_xlim(0, 3000); ax.legend()
    report = dict(mu=mu, edge_10_90_us=2*np.log(9)/mu*1e6,
                  window_s=[study['tmax']-.5, study['tmax']],
                  relative_tolerance=1e-5, gate_absolute_tolerance=1e-8,
                  voltage_absolute_tolerance_v=1e-3,
                  max_ibr_bus_kcl_error_a=max(kcl_errors),
                  passed=max(kcl_errors) <= .02 and all(row['passed'] for row in checks),
                  harmonics=checks)
    (directory / 'switching_validation.json').write_text(json.dumps(report, indent=2)+'\n')
    plots = directory / 'plots'; plots.mkdir(exist_ok=True)
    axs[0].set_ylabel('Gate peak amplitude [−]')
    axs[1].set_ylabel('Bridge peak amplitude [V]'); axs[1].set_xlabel('Frequency [Hz]')
    save(fig, plots, 'switching_prediction', f'IBR 4 — continuous PWM check, μ={mu:g}', [], time_axis=False)
    transient_detail(t, data, study, plots)
    if not report['passed']:
        raise ValueError(f'Switching check failed: {directory / "switching_validation.json"}')
    print(f'{directory}: {len(checks)} harmonic checks and IBR bus KCL passed')


def transient_detail(t, data, study, plots):
    """Use matching fault windows and axes across smoothing settings."""
    plots.mkdir(exist_ok=True)
    events = study['events']
    if events:
        first, last = events[0]['time'], events[-1]['time']
        windows = [(first-.02, last+.04), (first-.002, first+.004),
                   (last-.002, last+.004)]
        titles = ('Disturbance and recovery', 'First event', 'Last event')
        if events[0]['element_id'] == 'fault':
            titles = ('Fault and recovery', 'Fault inception', 'Fault clearing')
        elif len(events) == 1:
            windows[2] = (first+.02, first+.026)
            titles = ('Disturbance and recovery', 'Event detail', 'Recovery detail')
    else:
        windows = [(0., study['tmax']), (0., .006), (study['tmax']-.006, study['tmax'])]
        titles = ('Startup and steady operation', 'Startup detail', 'Steady switching')
    fig, axs = plt.subplots(2, 3, figsize=(15, 7))
    for column, (start, end) in enumerate(windows):
        selected = (t >= start) & (t <= end)
        values = np.column_stack((t[selected]*1000,
                                  data['Converter_converter_4_voa'][selected]/1000,
                                  data['DependentVoltageSource_filter_4_ia'][selected]))
        np.savetxt(plots / f'switching_window_{column}.csv', values, delimiter=',',
                   header='time_ms,voa_kv,ia_a', comments='')
        for row, color in enumerate(('#0072B2', '#E69F00')):
            ax = axs[row, column]
            ax.plot(values[:, 0], values[:, row+1], color=color, lw=.9)
            for event in events:
                if start <= event['time'] <= end:
                    ax.axvline(event['time']*1000, color='#555555', ls=':', lw=.9)
            ax.set_xlim(start*1000, end*1000); ax.set_xlabel('Time [ms]')
        axs[0, column].set_title(titles[column])
        if events and events[0]['element_id'] == 'fault':
            axs[0, column].set_ylim(-21, 21)
            axs[1, column].set_ylim(((-1100, 1100), (-200, 850), (-550, 550))[column])
    axs[0, 0].set_ylabel('Bridge phase-a voltage [kV]')
    axs[1, 0].set_ylabel('Phase-a injection into bus 4 [A]')
    save(fig, plots, 'switching_transient',
         f'Ten-bus IBR 4 — μ={study.get("mu", 240.):g}, 900 Hz PWM, '
         f'{study["dt_monitor"]*1e6:g} μs samples', [], time_axis=False)


def powers(v, i):
    p = np.sum(v * i, axis=1)
    alpha_v = (2*v[:, 0] - v[:, 1] - v[:, 2]) / 3
    beta_v = (v[:, 1] - v[:, 2]) / np.sqrt(3)
    alpha_i = (2*i[:, 0] - i[:, 1] - i[:, 2]) / 3
    beta_i = (i[:, 1] - i[:, 2]) / np.sqrt(3)
    return p, 1.5 * (beta_v * alpha_i - alpha_v * beta_i)


def event_lines(axes, events):
    for ax in np.asarray(axes).ravel():
        for event in events:
            ax.axvline(event['time'], color='#555555', ls=':', lw=.8)
        ax.set_xlabel('Time [s]')


def save(fig, directory, name, title, gallery, events=(), time_axis=True):
    if time_axis:
        event_lines(fig.axes, events)
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, .965))
    fig.savefig(directory / (name + '.png'))
    fig.savefig(directory / (name + '.svg'))
    plt.close(fig)
    gallery.append((name, title))


def one_line(case, root):
    positions = {1: (0, 3), 2: (3, 3), 3: (6, 3), 4: (0, 0), 5: (3, 0), 6: (6, 0),
                 7: (1, 1.5), 8: (5, 1.5), 9: (8, 1.5), 10: (-2, 1.5)}
    fig, ax = plt.subplots(figsize=(15, 8))
    ax.set_aspect('equal')
    for device in case['devices']:
        if device['class'] not in ('LineLumped', 'Switch'):
            continue
        a = int(device['inputs']['bus1'].split('_')[-1])
        b = int(device['inputs']['bus2'].split('_')[-1])
        x1, y1 = positions[a]; x2, y2 = positions[b]
        if device['class'] == 'LineLumped' and (a, b) == (7, 8):
            ax.plot([x1, x1+1, x2-1, x2], [y1, y1+.42, y2+.42, y2], color='#65737b', lw=1.6)
            ax.text(3, 2.02, 'L7–8: parallel R–L path', ha='center', fontsize=9)
        else:
            switch = device['class'] == 'Switch'
            opened = switch and device['params']['open']
            ax.plot([x1, x2], [y1, y2], color='#b44137' if switch else '#65737b',
                    lw=1.6, ls='--' if opened else '-')
            if switch:
                mid = ((x1+x2)/2, (y1+y2)/2)
                ax.add_patch(Rectangle((mid[0]-.10, mid[1]-.07), .20, .14,
                                       facecolor='white' if opened else '#b44137', edgecolor='#b44137', zorder=4))
                ax.text(mid[0], mid[1]-.27, device['id'] + (' (open)' if opened else ' (closed)'), ha='center', fontsize=9)
    for bus, (x, y) in positions.items():
        ax.plot([x-.25, x+.25], [y, y], color='#172b3a', lw=5, zorder=5)
        ax.text(x+.32, y+.06, f'B{bus}', weight='bold')
        if bus <= 3:
            ax.plot([x, x], [y, y+.4], color=COLORS[0])
            ax.add_patch(Circle((x, y+.64), .23, fill=False, color=COLORS[0], lw=2))
            ax.text(x, y+.64, 'G', ha='center', va='center', color=COLORS[0], weight='bold')
            ax.text(x, y+1.02, f'M{bus}: 10 MVA + TGOV1', ha='center', fontsize=10)
        elif bus in (4, 5, 6):
            ax.plot([x, x], [y, y-.35], color=COLORS[2])
            ax.add_patch(Rectangle((x-.42, y-.8), .84, .45, fill=False, color=COLORS[2], lw=2))
            ax.text(x, y-.575, f'IBR {bus-3}', ha='center', va='center', color=COLORS[2], weight='bold')
            ax.text(x, y-1.02, 'DC → PWM → bridge → R–L', ha='center', fontsize=9)
    for device in case['devices']:
        if device['class'] == 'LoadZ':
            bus = int(device['inputs']['bus'].split('_')[-1]); x, y = positions[bus]
            nominal = 13800**2 / device['params']['R'][0][0] / 1e6
            label = f'{nominal:g} MW @ 1 pu' if bus != 10 else '2 Ω/phase fault resistor'
            ax.annotate(label, xy=(x, y), xytext=(x-.15, y-.46 if bus > 6 else y+.4),
                        ha='center', fontsize=9, arrowprops=dict(arrowstyle='->', color='#555555'))
    ax.set_xlim(-3.2, 9.5); ax.set_ylim(-1.5, 4.5)
    ax.axis('off')
    fig.suptitle('Synthetic 10-bus EMT grid · 13.8 kV · 60 Hz', fontsize=17)
    fig.text(.5, .03, 'Exactly 10 buses · 3 machines · 3 open-loop converters · 9 lines · 7 loads · initial switch states shown', ha='center')
    fig.tight_layout(rect=(0, .06, 1, .95))
    for extension in ('png', 'svg', 'pdf'):
        fig.savefig(root / ('one_line.' + extension))
    plt.close(fig)


def plot_scenario(name, case, results):
    directory = results / name
    study = json.loads((HERE / (name + '.solver.json')).read_text())
    t, data = read_csv(directory / study['output_file'])
    ts, states = read_csv(directory / study['state_output_file'])
    if not np.array_equal(t, ts):
        raise ValueError(f'{name}: state and monitor time grids differ')
    events = study['events']
    plots = directory / 'plots'; plots.mkdir(exist_ok=True)
    gallery = []
    title = TITLES[name]
    emit = lambda fig, key, label: save(fig, plots, key, title + ' — ' + label, gallery, events)
    voltage = {b: phase(data, f'Bus_bus_{b}', 'v') for b in range(1, 11)}
    rms = {b: np.sqrt(np.maximum(0, cycle_mean(t, np.mean((v - np.roll(v, -1, axis=1))**2, axis=1)))) / 13800 for b, v in voltage.items()}
    speed = {b: data[f'Machine_machine_{b}_omega'] * 60 for b in (1, 2, 3)}
    machine_p = {b: data[f'Machine_machine_{b}_p'] for b in (1, 2, 3)}
    machine_q = {b: data[f'Machine_machine_{b}_q'] for b in (1, 2, 3)}
    converter_i = {b: phase(data, f'DependentVoltageSource_filter_{b}', 'i') for b in (4, 5, 6)}
    converter_pq = {b: powers(voltage[b], converter_i[b]) for b in (4, 5, 6)}
    devices = case['devices']
    loads = [d for d in devices if d['class'] == 'LoadZ']
    lines = [d for d in devices if d['class'] == 'LineLumped']
    load_power = {d['id']: -np.sum(voltage[int(d['inputs']['bus'].split('_')[-1])] * phase(data, 'LoadZ_'+d['id'], 'i'), axis=1) for d in loads}

    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for b in range(1, 9): axs[0].plot(t, rms[b], label=f'B{b}')
    for b in (9, 10): axs[1].plot(t, rms[b], label=f'B{b}')
    for ax in axs: ax.set_ylabel('Line–line RMS [pu]'); ax.legend(ncol=4)
    emit(fig, '01_bus_rms', 'all bus voltages (trailing one-cycle RMS)')

    fig, axs = plt.subplots(5, 2, figsize=(14, 13), sharex=True)
    window = (t >= .98) & (t <= 1.10)
    for b, ax in zip(range(1, 11), axs.ravel()):
        for n, p in enumerate('abc'): ax.plot(t[window], voltage[b][window, n]/1000, label=p)
        ax.set_title(f'Bus {b}'); ax.set_ylabel('kV'); ax.set_xlim(.98, 1.10)
    axs[0, 0].legend(ncol=3)
    emit(fig, '02_bus_waveforms', 'three-phase event waveforms')

    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for b in (1, 2, 3):
        axs[0].plot(t, speed[b], label=f'M{b}')
        theta = data[f'Machine_machine_{b}_theta']
        axs[1].plot(t, np.rad2deg(theta - 2*np.pi*60*t - theta[0]), label=f'M{b}')
    axs[0].set_ylabel('Rotor electrical frequency [Hz]'); axs[1].set_ylabel('Angle deviation [deg]')
    axs[0].legend(ncol=3)
    emit(fig, '03_machine_motion', 'machine frequency and angle relative to 60 Hz')

    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for b in (1, 2, 3):
        axs[0].plot(t, cycle_mean(t, machine_p[b])/1e6, label=f'M{b}')
        axs[1].plot(t, cycle_mean(t, machine_q[b])/1e6, label=f'M{b}')
    axs[0].set_ylabel('Injected P [MW]'); axs[1].set_ylabel('Injected Q [Mvar]'); axs[0].legend(ncol=3)
    emit(fig, '04_machine_power', 'machine power (trailing one-cycle mean)')

    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for b in (4, 5, 6):
        for n, ax in enumerate(axs): ax.plot(t, cycle_mean(t, converter_pq[b][n])/1e6, label=f'IBR {b-3} at B{b}')
    axs[0].set_ylabel('Injected P [MW]'); axs[1].set_ylabel('Injected Q [Mvar]'); axs[0].legend(ncol=3)
    emit(fig, '05_converter_power', 'converter AC-terminal power (positive into grid)')

    fig, axs = plt.subplots(3, 2, figsize=(13, 9), sharex=True)
    window = (t >= .99) & (t <= 1.025)
    for b, row in zip((4, 5, 6), axs):
        for p in 'abc':
            row[0].plot(t[window], data[f'PWM_pwm_{b}_s{p}'][window], label=p)
            row[1].plot(t[window], data[f'Converter_converter_{b}_vo{p}'][window]/1000, label=p)
        row[0].set_ylabel(f'IBR {b-3}: s [–]'); row[1].set_ylabel('Bridge voltage [kV]')
        for ax in row: ax.set_xlim(.99, 1.025)
    axs[0, 0].legend(ncol=3)
    emit(fig, '06_pwm_bridge', 'actual smoothed PWM and bridge voltage')

    fig, axs = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    window = (t >= .97) & (t <= 1.12)
    for b, ax in zip((4, 5, 6), axs):
        for n, p in enumerate('abc'): ax.plot(t[window], converter_i[b][window, n], label=p)
        ax.set_ylabel(f'IBR {b-3} current [A]'); ax.set_xlim(.97, 1.12)
    axs[0].legend(ncol=3)
    emit(fig, '07_converter_currents', 'filter phase currents around the first event')

    fig, axs = plt.subplots(3, 3, figsize=(14, 10), sharex=True)
    for device, ax in zip(lines, axs.ravel()):
        current = phase(data, 'LineLumped_'+device['id'], 'i12')
        ax.plot(t, np.sqrt(np.maximum(0, cycle_mean(t, np.mean(current**2, axis=1)))))
        ax.set_title(device['id']); ax.set_ylabel('Phase RMS [A]')
    emit(fig, '08_line_currents', 'all nine line currents')

    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for key, power in load_power.items(): axs[1 if key == 'load_10' else 0].plot(t, cycle_mean(t, power)/1e6, label=key)
    for ax in axs: ax.set_ylabel('Consumed P [MW]'); ax.legend(ncol=3)
    emit(fig, '09_load_power', 'all loads and fault-resistor dissipation')

    fig, axs = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    for n, label in enumerate(('Turbine block output', 'Valve position', 'Mechanical power')):
        for b in (1, 2, 3): axs[n].plot(t, states[f'y:governor_{b}[{n}]'], label=f'GOV{b}')
        axs[n].set_ylabel(label+' [pu]')
    axs[0].legend(ncol=3)
    emit(fig, '10_governors', 'all TGOV1 states on the 10 MVA machine base')

    fig, axs = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
    for b, row in zip((1, 2, 3), axs):
        for n, label in ((2, 'psi_d'), (3, 'psi_q'), (5, 'psi_fd')):
            row[0].plot(t, states[f'y:machine_{b}[{n}]'], label=label)
        for key in ('ifd', 'ks', 'psi_at'): row[1].plot(t, data[f'Machine_machine_{b}_{key}'], label=key)
        row[0].set_ylabel(f'M{b} flux [pu]'); row[1].set_ylabel(f'M{b} [pu / –]')
    axs[0, 0].legend(ncol=3); axs[0, 1].legend(ncol=3)
    emit(fig, '11_machine_fields', 'rotor flux, field current and saturation')

    fig, axs = plt.subplots(3, 2, figsize=(13, 9), sharex=True)
    for switch, row in zip(('load_step', 'fault', 'tie'), axs):
        row[0].step(t, data[f'Switch_{switch}_open'], where='post'); row[0].set_ylabel(switch+' open [0/1]')
        for p in 'abc': row[1].plot(t, data[f'Switch_{switch}_i12{p}'], label=p)
        row[1].set_ylabel('Switch current [A]')
    axs[0, 1].legend(ncol=3)
    emit(fig, '12_switches', 'switch commands and phase currents')

    # Independent KCL reconstruction from recorded physical branch quantities.
    kcl = {b: -phase(data, f'Bus_bus_{b}', 'i_sh') for b in voltage}
    for b in (1, 2, 3): kcl[b] += phase(data, f'Machine_machine_{b}', 'i')
    for b in (4, 5, 6): kcl[b] += converter_i[b]
    for device in loads:
        b = int(device['inputs']['bus'].split('_')[-1])
        kcl[b] += phase(data, 'LoadZ_'+device['id'], 'i')
    losses = np.zeros_like(t); stored_rate = np.zeros_like(t)
    for device in lines:
        a, b = [int(device['inputs'][key].split('_')[-1]) for key in ('bus1', 'bus2')]
        prefix = 'LineLumped_'+device['id']; current = phase(data, prefix, 'i12')
        kcl[a] -= current
        kcl[b] += current
        r, l, c = [np.array(device['params'][key]) * device['params']['dx'] for key in ('Rp', 'Lp', 'Cp')]
        derivative = np.column_stack([states[f"yp:{device['id']}[{n}]"] for n in range(3)])
        losses += np.einsum('ni,ij,nj->n', current, r, current)
        stored_rate += np.einsum('ni,ij,nj->n', current, l, derivative)
        for bus in (a, b):
            derivative = np.column_stack([states[f'yp:bus_{bus}[{n}]'] for n in range(3)])
            stored_rate += .5 * np.einsum('ni,ij,nj->n', voltage[bus], c, derivative)
    for device in devices:
        if device['class'] == 'Switch':
            a, b = [int(device['inputs'][key].split('_')[-1]) for key in ('bus1', 'bus2')]
            current = phase(data, 'Switch_'+device['id'], 'i12')
            kcl[a] -= current; kcl[b] += current
    total_machine = sum(machine_p.values())
    total_converter = sum(pq[0] for pq in converter_pq.values())
    total_load = sum(load_power.values())
    imbalance = total_machine + total_converter - total_load - losses - stored_rate
    max_kcl = float(max(np.max(np.abs(value)) for value in kcl.values()))
    max_bridge = 0.; gate_min, gate_max = 1., 0.
    for b in (4, 5, 6):
        s = phase(data, f'PWM_pwm_{b}', 's'); e = phase(data, f'Converter_converter_{b}', 'vo')
        dc = study.get('signal_values', {}).get(f'dc_{b}', next(d['value'] for d in case['signals'] if d['id'] == f'dc_{b}'))
        expected = dc * (s - np.mean(s, axis=1, keepdims=True))
        max_bridge = max(max_bridge, float(np.max(np.abs(e - expected))))
        gate_min = min(gate_min, float(s.min())); gate_max = max(gate_max, float(s.max()))
    if max_kcl > .02 or max_bridge > 1e-7 or gate_min < 0 or gate_max > 1:
        raise ValueError(f'{name}: physical identity checks failed: KCL={max_kcl}, bridge={max_bridge}, gates={gate_min,gate_max}')
    fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    for value, label in ((total_machine, 'Machines'), (total_converter, 'Converters'), (total_load, 'Loads')):
        axs[0].plot(t, cycle_mean(t, value)/1e6, label=label)
    axs[0].set_ylabel('P [MW]'); axs[0].legend(ncol=3)
    for b in voltage: axs[1].plot(t, np.max(np.abs(kcl[b]), axis=1), label=f'B{b}')
    axs[1].set_ylabel('Max phase |KCL| [A]'); axs[1].set_yscale('symlog', linthresh=1e-9)
    axs[2].plot(t, imbalance); axs[2].set_ylabel('Power balance error [W]')
    emit(fig, '13_balance_checks', 'network power and independent conservation checks')

    # A coherent 0.5 s window at the end: no fabricated switching ripple.
    mask = (t >= 2.5) & (t < 3.)
    signal = data['Converter_converter_4_voa'][mask]
    frequency, amplitude = spectrum(signal, study['dt_monitor'])
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.semilogy(frequency, np.maximum(amplitude, 1e-12)); ax.set_xlim(0, 2000)
    ax.set_xlabel('Frequency [Hz]'); ax.set_ylabel('Peak voltage amplitude [V]')
    ax.axvline(60, color=COLORS[1], ls=':', label='60 Hz'); ax.axvline(900, color=COLORS[2], ls=':', label='900 Hz carrier'); ax.legend()
    save(fig, plots, '14_bridge_spectrum', title+' — bridge spectrum, 2.5–3.0 s', gallery, time_axis=False)

    # Every DAE y variable has a labeled trace, grouped in a multipage PDF too.
    schema = json.loads((directory / (study['state_output_file']+'.json')).read_text())
    variables = sorted(schema['variables'], key=lambda v: v['index'])
    with PdfPages(plots / 'all_dae_variables.pdf') as pdf:
        for first in range(0, len(variables), 18):
            fig, axs = plt.subplots(6, 3, figsize=(15, 17), sharex=True)
            for variable, ax in zip(variables[first:first+18], axs.ravel()):
                key = f"y:{variable['component']}[{variable['local_index']}]"
                ax.plot(t[::10], states[key][::10], lw=.8)
                ax.set_title(key + (' (differential)' if variable['differential'] else ' (algebraic)'), fontsize=9)
                ax.tick_params(labelsize=8)
            for ax in axs.ravel()[len(variables[first:first+18]):]: ax.axis('off')
            event_lines(axs, events)
            page = first//18 + 1
            fig.suptitle(f'{title} — all DAE variables, page {page}\nNative component units; display every tenth sample, raw CSV retains all samples', fontsize=13)
            fig.tight_layout(rect=(0, 0, 1, .965)); pdf.savefig(fig)
            key = f'states_{page:02d}'
            fig.savefig(plots / (key+'.png')); plt.close(fig)
            gallery.append((key, f'All DAE variables, page {page}'))
    after = t >= 1.; settled = t >= 2.5
    metrics = dict(scenario=name, title=title, mu=study.get('mu',240.0), dt_monitor=study['dt_monitor'], samples=len(t), monitor_columns=len(data), dae_variables=len(variables),
                   min_main_bus_rms_after_event_pu=float(min(np.nanmin(rms[b][after]) for b in range(1,9))),
                   min_machine_frequency_after_event_hz=float(min(np.min(speed[b][after]) for b in speed)),
                   max_machine_frequency_after_event_hz=float(max(np.max(speed[b][after]) for b in speed)),
                   final_mean_machine_mw=float(np.mean(total_machine[settled])/1e6),
                   final_mean_converter_mw=float(np.mean(total_converter[settled])/1e6),
                   max_kcl_error_a=max_kcl, max_bridge_identity_error_v=max_bridge,
                   max_power_balance_error_w=float(np.max(np.abs(imbalance))), pwm_range=[gate_min, gate_max])
    (directory / 'metrics.json').write_text(json.dumps(metrics, indent=2)+'\n')
    cards = ''.join(f'<figure><a href="plots/{key}.png"><img loading="lazy" src="plots/{key}.png" alt="{html.escape(label)}"></a><figcaption>{html.escape(label)}</figcaption></figure>' for key, label in gallery)
    events_text = ', '.join(f"{e['time']:g} s: {'open' if e['open'] else 'close'} {e['element_id']}" for e in events) or 'No switching events'
    (directory / 'index.html').write_text(page_html(title, f'<p><a href="../index.html">All scenarios</a> · {events_text}</p><p><a href="{study["output_file"]}">Model monitors CSV</a> · <a href="{study["state_output_file"]}">All states and derivatives CSV</a> · <a href="{study["state_output_file"]}.json">State index map</a> · <a href="plots/all_dae_variables.pdf">All DAE traces PDF</a> · <a href="metrics.json">Checks and metrics</a> · <a href="run.log">Solver log</a></p><p>Positive source P/Q means injection into the grid. RMS and displayed powers use a trailing 1/60 s window. Startup is retained. Dotted lines mark events.</p><div class="gallery">{cards}</div>'))
    print(name, json.dumps(metrics), flush=True)
    return metrics, comparison_trace(t, data)


def comparison_trace(t, data):
    v = phase(data, 'Bus_bus_7', 'v')
    rms = np.sqrt(np.maximum(0, cycle_mean(t, np.mean((v - np.roll(v, -1, axis=1))**2, axis=1)))) / 13800
    converter = sum(powers(phase(data, f'Bus_bus_{b}', 'v'),
                           phase(data, f'DependentVoltageSource_filter_{b}', 'i'))[0]
                    for b in (4, 5, 6))
    load = -sum(np.sum(phase(data, f'Bus_bus_{b}', 'v') *
                       phase(data, f'LoadZ_load_{b}', 'i'), axis=1)
                for b in range(4, 11))
    return dict(t=t[::10], voltage=rms[::10],
                frequency=data['Machine_machine_1_omega'][::10]*60,
                converter=cycle_mean(t, converter)[::10]/1e6,
                load=cycle_mean(t, load)[::10]/1e6)



def switching_comparison(results):
    series = {}
    for name in ('04_FaultClearing', '06_FaultClearingSwitching'):
        t, data = read_csv(results/name/(name+'.csv'))
        study = json.loads((HERE/(name+'.solver.json')).read_text())
        series[name] = (t, data, study)
    fig, axs = plt.subplots(3, 1, figsize=(13, 10), sharex=True)
    for name, (t, data, study) in series.items():
        mask = (t >= .970) & (t <= .976)
        label = f"mu={study.get('mu', 240):g}"
        for ax, key, scale, unit in zip(axs,
                ('PWM_pwm_4_sa', 'Converter_converter_4_voa', 'DependentVoltageSource_filter_4_ia'),
                (1, 1e-3, 1), ('Phase-a switching function', 'Bridge phase-a voltage [kV]', 'Filter phase-a current [A]')):
            ax.plot(t[mask], data[key][mask]*scale, label=label)
            ax.set_ylabel(unit)
    axs[0].legend()
    gallery = []
    save(fig, results, 'switching_waveform_comparison',
         'Fault-study pre-event waveforms — common DC voltage', gallery)
    fig, axs = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
    spectral_metrics = []
    for name, (t, data, study) in series.items():
        mask = (t >= 2.5) & (t < 3.0)
        label = f"mu={study.get('mu',240):g}"
        for ax, key, unit in zip(axs, ('PWM_pwm_4_sa','Converter_converter_4_voa'),
                                 ('Switching-function peak amplitude', 'Bridge-voltage peak amplitude [V]')):
            signal = data[key][mask]
            frequency, amplitude = spectrum(signal, study['dt_monitor'])
            ax.semilogy(frequency, np.maximum(amplitude,1e-14), label=label)
            ax.set_ylabel(unit); ax.set_xlim(0,3000); ax.set_ylim(bottom=1e-12)
            carrier = int(np.argmin(abs(frequency-900)))
            spectral_metrics.append(dict(scenario=name, signal=key, window_s=[2.5,3.0],
                                         carrier_frequency_hz=float(frequency[carrier]),
                                         carrier_peak_amplitude=float(amplitude[carrier]),
                                         sideband_780_peak=float(amplitude[np.argmin(abs(frequency-780))]),
                                         sideband_1020_peak=float(amplitude[np.argmin(abs(frequency-1020))])))
    for ax in axs:
        ax.set_xlabel('Frequency [Hz]'); ax.axvline(900,color='#777777',ls=':',label='900 Hz carrier'); ax.legend()
    save(fig, results, 'switching_spectrum_comparison',
         'Carrier and sidebands — coherent final 0.5 s, rectangular window', gallery, time_axis=False)
    (results/'switching_spectra.json').write_text(json.dumps(spectral_metrics,indent=2)+'\n')


def page_html(title, body):
    return f'<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title><style>body{{font:16px system-ui;margin:2rem auto;max-width:1500px;padding:0 1rem;color:#172b3a;background:#f5f7fa}}a{{color:#176b9c}}.gallery{{display:grid;grid-template-columns:repeat(auto-fit,minmax(450px,1fr));gap:1rem}}figure{{margin:0;background:white;padding:.7rem;border:1px solid #dce2e8}}img{{width:100%}}figcaption{{padding:.5rem}}table{{border-collapse:collapse;width:100%;background:white}}td,th{{padding:.6rem;border:1px solid #dce2e8;text-align:left}}h1{{font-size:1.8rem}}@media(max-width:600px){{.gallery{{display:block}}}}</style><h1>{html.escape(title)}</h1>{body}</html>'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=HERE/'results')
    parser.add_argument('--reuse-plots', action='store_true', help='Reuse existing per-scenario galleries for unchanged outputs')
    parser.add_argument('--switching-only', action='store_true', help='Plot and check only the resolved switching study')
    parser.add_argument('--transients-only', action='store_true', help='Plot matching voltage/current panels for available runs')
    args = parser.parse_args(); results = args.results.resolve(); results.mkdir(parents=True, exist_ok=True)
    if args.transients_only:
        cards = []
        for directory in sorted(results.iterdir()):
            if not (directory / 'run.solver.json').exists():
                continue
            study = json.loads((directory / 'run.solver.json').read_text())
            t, data = read_csv(directory / study['output_file'])
            transient_detail(t, data, study, directory / 'plots')
            cards.append(f'<figure><figcaption>{html.escape(directory.name)}</figcaption>'
                         f'<a href="{directory.name}/plots/switching_transient.svg">'
                         f'<img src="{directory.name}/plots/switching_transient.png"></a></figure>')
        (results / 'transients.html').write_text(page_html('Ten-bus IBR transients', ''.join(cards)))
        return
    if args.switching_only:
        switching_detail(results / '06_FaultClearingSwitching')
        return
    study = json.loads((HERE/'01_Baseline.solver.json').read_text())
    case = json.loads((HERE/study['system_model_file']).read_text())
    one_line(case, results)
    metrics, traces = [], {}
    for name in TITLES:
        if args.reuse_plots and (results/name/'index.html').exists():
            metric = json.loads((results/name/'metrics.json').read_text())
            t, data = read_csv(results/name/(name+'.csv'))
            trace = comparison_trace(t, data)
        else:
            metric, trace = plot_scenario(name, case, results)
        study = json.loads((HERE/(name+'.solver.json')).read_text())
        metric.setdefault('mu', study.get('mu',240.0))
        metric.setdefault('dt_monitor', study['dt_monitor'])
        metrics.append(metric); traces[name] = trace
    fig, axs = plt.subplots(4, 1, figsize=(13, 13), sharex=True)
    for name, trace in traces.items():
        for ax, key, label in zip(axs, ('voltage','frequency','converter','load'), ('B7 line–line RMS [pu]','M1 frequency [Hz]','Converter P [MW]','Load P [MW]')):
            ax.plot(trace['t'], trace[key], label=TITLES[name]); ax.set_ylabel(label)
    axs[0].legend(fontsize=9, ncol=2)
    gallery = []
    save(fig, results, 'scenario_comparison', 'Six sparse-solver scenarios — common network and event time base', gallery)
    fig, axs = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    baseline = traces['01_Baseline']
    for name, trace in traces.items():
        if name == '01_Baseline': continue
        for ax, key, label in zip(axs, ('frequency','voltage'), ('M1 frequency change vs baseline [Hz]', 'B7 RMS change vs baseline [pu]')):
            ax.plot(trace['t'], trace[key] - np.interp(trace['t'], baseline['t'], baseline[key]), label=TITLES[name]); ax.set_ylabel(label)
    axs[0].legend(fontsize=9,ncol=2)
    save(fig, results, 'baseline_differences', 'Disturbance response relative to the undisturbed run', gallery)
    switching_comparison(results)
    switching_detail(results / '06_FaultClearingSwitching')
    (results/'summary.json').write_text(json.dumps(metrics,indent=2)+'\n')
    with (results/'summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=sorted(metrics[0]));writer.writeheader();writer.writerows(metrics)
    rows=''.join(f'<tr><td><a href="{m["scenario"]}/index.html">{m["title"]}</a></td><td>{m["samples"]:,}</td><td>{m["min_main_bus_rms_after_event_pu"]:.4f}</td><td>{m["min_machine_frequency_after_event_hz"]:.4f}–{m["max_machine_frequency_after_event_hz"]:.4f}</td><td>{m["max_kcl_error_a"]:.2e}</td></tr>' for m in metrics)
    (results/'index.html').write_text(page_html('Synthetic EMT 10-bus simulation review', f'<p>Six actual GridKit/IDA sparse KLU runs, each 0–3 s; five at 100 µs monitoring and one at 10 µs. Three machines, three open-loop PWM converters and seven loads. Full startup and pre/post-event samples are retained.</p><p><b>Model scope:</b> ideal constant DC sources, fixed 60 Hz modulation, no PLL, current limit, DC-link dynamics or converter protection. The first five use μ=240 s⁻¹ and the sixth uses μ=50000 s⁻¹. All use the same 28.733 kV DC sources; larger μ reveals switching features. These are synthetic network-response demonstrations, not a hardware switching-loss or ride-through validation.</p><p><a href="../README.md">Case/scenario documentation</a> · <a href="summary.csv">Metrics CSV</a> · <a href="summary.json">Metrics JSON</a> · <a href="one_line.svg">One-line SVG</a> · <a href="one_line.pdf">One-line PDF</a></p><figure><img src="one_line.png" alt="Ten bus one-line diagram"></figure><table><tr><th>Scenario / plots and data</th><th>Samples</th><th>Minimum B1–B8 RMS after 1 s [pu]</th><th>Machine frequency range after 1 s [Hz]</th><th>Max KCL error [A]</th></tr>{rows}</table><p>RMS values use a trailing one-cycle window. Startup minima are excluded from this comparison; all raw startup data remain available.</p><div class="gallery"><figure><a href="switching_waveform_comparison.png"><img src="switching_waveform_comparison.png"></a></figure><figure><a href="switching_spectrum_comparison.png"><img src="switching_spectrum_comparison.png"></a></figure><figure><a href="scenario_comparison.png"><img src="scenario_comparison.png"></a></figure><figure><a href="baseline_differences.png"><img src="baseline_differences.png"></a></figure></div>'))


if __name__ == '__main__':
    main()
