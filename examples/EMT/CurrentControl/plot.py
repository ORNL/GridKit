#!/usr/bin/env python3
"""Plot resolved switching-control simulations and retain compressed waveforms."""
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
OUT = HERE / 'simulation'
FC = 6000.0
MU = 1e6
VDC = 400.0
UMAX = np.sqrt(3 / 8) * .95 * VDC
BLUE, ORANGE, GREY = '#0072B2', '#D55E00', '#555555'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.alpha': .2, 'lines.linewidth': 1.2,
                     'savefig.dpi': 180, 'figure.constrained_layout.use': True})


def read(name):
    csv = OUT / f'{name}.csv'
    if csv.exists():
        data = np.genfromtxt(csv, delimiter=',', names=True)
        data = {n: data[n] for n in data.dtype.names}
        np.savez_compressed(OUT / f'{name}.npz', **data)
        return data
    with np.load(OUT / f'{name}.npz') as data:
        return {n: data[n] for n in data.files}


def save(fig, name):
    for ext in ['png', 'pdf']:
        fig.savefig(OUT / f'{name}.{ext}')
    plt.close(fig)


def mean(t, x, begin, end):
    sel = (t >= begin) & (t <= end)
    return float(np.trapz(x[sel], t[sel]) / (t[sel][-1] - t[sel][0]))


def carrier_mean(t, x):
    """Time-weighted means for plotting, without changing the simulated model."""
    k = np.floor(t * FC).astype(int)
    dt = np.diff(t)
    weight = np.bincount(k[:-1], weights=dt)
    accum = np.bincount(k[:-1], weights=.5 * (x[:-1] + x[1:]) * dt)
    good = weight > 0
    return (np.arange(len(weight))[good] + .5) / FC, accum[good] / weight[good]


def trace(ax, t, x, label, color):
    ax.plot(t, x, color=color, alpha=.18, linewidth=.45)
    tm, xm = carrier_mean(t, x)
    ax.plot(tm, xm, color=color, label=label)


def predictor(d):
    """Reconstruct pulses from sampled commands and their documented edge times."""
    t = d['t']
    k = np.floor(t * FC + 1e-10).astype(int)
    sample_time = np.arange(k.max() + 1) / FC
    sample_indices = np.searchsorted(t, sample_time - 1e-13)
    if np.max(np.abs(t[sample_indices] - sample_time)) > 1e-10:
        raise ValueError('Carrier-boundary command samples are missing')
    modulation = np.column_stack([d[f'Modulation_modulation_m{p}'][sample_indices] for p in 'abc'])
    switches = np.zeros((len(t), 3))
    # Three nearby intervals cover all tails at this carrier and edge width.
    for offset in [-1, 0, 1]:
        interval = k + offset
        # All periodic replicas use the command active in interval k.
        idx = np.clip(k - 1, 0, len(modulation) - 1)
        duty = (1 + modulation[idx]) / 2
        on = (interval[:, None] + .5 * (1 - duty)) / FC
        off = (interval[:, None] + .5 * (1 + duty)) / FC
        sigmoid = lambda x: .5 * (1 + np.tanh(.5 * MU * x))
        switches += sigmoid(t[:, None] - on) - sigmoid(t[:, None] - off)
    return switches


def harmonics(d):
    """Integrate each held command's periodic PWM waveform analytically."""
    t = d['t']
    begin, end = .2 - 1 / 60, .2
    intervals = np.arange(round(begin * FC), round(end * FC))
    idx = np.searchsorted(t, (intervals - 1) / FC - 1e-13)
    m = np.column_stack([d[f'Modulation_modulation_m{p}'][idx] for p in 'abc'])
    duty = (1 + m) / 2
    # Fourier coefficients of the centered pulse, including sigmoid smoothing.
    decay = 2 * np.pi**2 * FC / MU
    count = max(1, int(np.ceil(np.log(4 / np.finfo(float).eps) / decay)))
    orders = np.arange(-count, count + 1)
    argument = decay * np.abs(orders)
    attenuation = np.ones(len(orders))
    nonzero = orders != 0
    x = argument[nonzero]
    attenuation[nonzero] = 2 * x * np.exp(-x) / (-np.expm1(-2 * x))
    coefficients = (duty[:, :, None] * np.sinc(duty[:, :, None] * orders)
                    * np.exp(-1j * np.pi * orders) * attenuation)
    bridge_coefficients = VDC * np.einsum('p,kpn->kn', [2/3, -1/3, -1/3], coefficients)
    inside = (t > begin) & (t < end)
    time = np.r_[begin, t[inside], end]
    voltage = np.interp(time, t, d['Converter_bridge_voa'])
    frequencies = np.array([60, 5880, 5940, 6000, 6060, 6120, 11940, 12060, 17880, 18120])
    measured, predicted = [], []
    for f in frequencies:
        omega = 2 * np.pi * f
        difference = orders - f / FC
        integral = np.exp(1j * np.pi * difference) * np.sinc(difference) / FC
        phase = np.exp(-1j * omega * intervals / FC)
        predicted.append(abs(2 * np.sum(phase[:, None] * bridge_coefficients * integral) / (end - begin)))
        measured.append(abs(2 * np.trapz(voltage * np.exp(-1j * omega * time), time) / (end - begin)))
    return frequencies, np.array(measured), np.array(predicted)


def main():
    data = {name: read(name) for name in ['GFL', 'GFM']}
    summary = {'carrier_hz': FC, 'edge_10_90_us': 2 * np.log(9) / MU * 1e6,
               'formal_tests_run': False, 'studies': {}}
    for name, d in data.items():
        t = d['t']
        if not all(np.isfinite(x).all() for x in d.values()) or (np.diff(t) < 0).any():
            raise ValueError(f'{name}: nonfinite or nonmonotonic output')
        s = np.column_stack([d[f'PWM_pwm_s{p}'] for p in 'abc'])
        bridge = np.column_stack([d[f'Converter_bridge_vo{p}'] for p in 'abc'])
        current = np.column_stack([d[f'DependentVoltageSource_filter_i{p}'] for p in 'abc'])
        ilim = np.hypot(d['InnerCurrentControl_current_control_ilimd'], d['InnerCurrentControl_current_control_ilimq'])
        u = np.hypot(d['InnerCurrentControl_current_control_ud'], d['InnerCurrentControl_current_control_uq'])
        pred = predictor(d)
        # The two observations at a sampling instant straddle the latch update.
        interior = np.abs(t * FC - np.round(t * FC)) > .05
        final = (.18, .2)
        summary['studies'][name] = {
            'final_time_s': float(t[-1]), 'samples': len(t),
            'maximum_reference_norm_A': float(ilim.max()), 'maximum_voltage_command_V': float(u.max()),
            'bridge_identity_max_error_V': float(np.max(np.abs(bridge - VDC * (s - s.mean(axis=1)[:, None])))),
            'bridge_power_max_error_W': float(np.max(np.abs((bridge * current).sum(axis=1) - VDC * d['Converter_bridge_idc']))),
            'pulse_reconstruction_max_error': float(np.max(np.abs(pred[interior] - s[interior]))),
            'final_mean_id_A': mean(t, d['Park_current_y1'], *final),
            'final_mean_iq_A': mean(t, d['Park_current_y2'], *final),
            'final_mean_vd_V': mean(t, d['Park_voltage_y1'], *final),
            'final_mean_vq_V': mean(t, d['Park_voltage_y2'], *final),
        }
        fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        if name == 'GFL':
            req = np.select([t < .04, t < .08, t < .12], [8, 16, 45], default=8)
            ax[0].plot(t, req, '--', color=GREY, label='requested $i_d$')
            ax[0].plot(t, d['InnerCurrentControl_current_control_ilimd'], color=ORANGE, label='limited $i_d$')
            trace(ax[0], t, d['Park_current_y1'], '$i_d$', BLUE)
            ax[0].set_ylabel('d-axis current [A]')
            trace(ax[1], t, d['Park_current_y2'], '$i_q$', BLUE)
            ax[1].axhline(0, color=ORANGE, linestyle='--', label='requested $i_q$')
            ax[1].set_ylabel('q-axis current [A]')
            ax[0].set_title('Current tracking, reference limiting, and recovery')
            summary['studies'][name]['limited_interval_mean_id_A'] = mean(t, d['Park_current_y1'], .11, .12)
        else:
            trace(ax[0], t, d['Park_voltage_y1'], '$v_d$', BLUE)
            ax[0].axhline(208, color=ORANGE, linestyle='--', label='$v_d^{ref}=208$ V')
            ax[0].set_ylabel('d-axis voltage [V]')
            trace(ax[1], t, d['Park_voltage_y2'], '$v_q$', BLUE)
            ax[1].axhline(0, color=ORANGE, linestyle='--', label='$v_q^{ref}=0$ V')
            ax[1].set_ylabel('q-axis voltage [V]')
            ax[0].set_title('Islanded voltage control: connect and disconnect a second load')
        ax[2].plot(t, u, color=BLUE, label=r'$\|\mathbf{u}\|_2$')
        ax[2].axhline(UMAX, color=ORANGE, linestyle='--', label='available voltage command')
        ax[2].set_ylabel('Voltage command [V]')
        ax[2].set_xlabel('Time [s]')
        for a in ax:
            for event in ([.04, .08, .12] if name == 'GFL' else [.04, .12]):
                a.axvline(event, color=GREY, linewidth=.7, linestyle=':')
            a.legend(loc='upper right', ncol=3, fontsize=9)
        ax[2].legend(loc='lower right', ncol=2, fontsize=9)
        fig.supxlabel('Feedback traces: faint instantaneous values; solid carrier-period means.', fontsize=9)
        save(fig, name)
    d = data['GFM']; t = d['t']; mask = (t >= .18) & (t <= .181)
    fig, ax = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    ax[0].plot(1e3 * t[mask], d['PWM_pwm_sa'][mask], color=BLUE, label='$s_a$')
    ax[0].plot(1e3 * t[mask], predictor(d)[mask, 0], '--', color=ORANGE, label='pulse-edge reconstruction')
    ax[0].set_ylabel('Gate state [−]'); ax[0].set_title('Resolved 6 kHz switching, islanded steady operation')
    ax[1].plot(1e3 * t[mask], d['Converter_bridge_voa'][mask], color=BLUE, label='bridge $e_a$')
    ax[1].plot(1e3 * t[mask], d['Bus_capacitor_va'][mask], color=ORANGE, label='capacitor $v_{o,a}$')
    ax[1].set_ylabel('Phase voltage [V]')
    ax[2].plot(1e3 * t[mask], d['DependentVoltageSource_filter_ia'][mask], color=BLUE, label='inverter-side $i_{i,a}$')
    ax[2].plot(1e3 * t[mask], d['LineLumped_grid_filter_i12a'][mask], color=ORANGE, label='grid-side $i_{g,a}$')
    ax[2].set_ylabel('Phase current [A]'); ax[2].set_xlabel('Time [ms]')
    for a in ax: a.legend(loc='best', ncol=2)
    save(fig, 'switching')
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for ax, (name, d) in zip(axes, data.items()):
        frequencies, measured, predicted = harmonics(d)
        x = np.arange(len(frequencies))
        ax.bar(x - .18, measured, .36, color=BLUE, label='simulated bridge voltage')
        ax.bar(x + .18, predicted, .36, color=ORANGE, label='analytic PWM prediction')
        ax.set_yscale('log')
        ax.set_ylabel('Peak amplitude [V]')
        ax.set_title(name)
        ax.legend(loc='upper right')
        summary['studies'][name]['harmonics'] = {
            'frequencies_hz': frequencies.tolist(), 'measured_peak_V': measured.tolist(),
            'predicted_peak_V': predicted.tolist(),
            'maximum_absolute_error_V': float(np.max(np.abs(measured - predicted))),
        }
    axes[-1].set_xticks(x, [str(f) for f in frequencies])
    axes[-1].set_xlabel('Frequency [Hz]')
    save(fig, 'harmonics')
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
