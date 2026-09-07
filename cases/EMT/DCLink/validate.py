"""Validate a DC capacitor and converter against an exact RC transient; optionally plot."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent


def validate(executable, output, plot=False):
    output.mkdir(parents=True, exist_ok=True)
    run = subprocess.run([str(executable), str(ROOT / 'DCLink.solver.json')],
                         cwd=output, text=True, capture_output=True, check=False)
    (output / 'simulation.log').write_text(run.stdout + run.stderr)
    if run.returncode:
        raise RuntimeError(run.stdout + run.stderr)
    with (output / 'mon.csv').open(newline='') as stream:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    with (output / 'state.csv').open(newline='') as stream:
        states = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    layout = json.loads((output / 'state.csv.json').read_text())['variables']
    assert [item['component'] for item in layout if item['differential']] == ['dc']
    assert len(rows) == len(states) == 3002  # Both sides of the event are sampled.
    capacitance, resistance, event = 0.02, 15.0, 0.5
    tau = resistance * capacitance
    at_event = 1200.0 - 600.0 * math.exp(-event / tau)
    exact = [1200.0 - 600.0 * math.exp(-r['t'] / tau) if r['t'] <= event
             else 300.0 + (at_event - 300.0) * math.exp(-(r['t'] - event) / tau)
             for r in rows]
    voltage, source, current, energy, power = [], [], [], [], []
    residual_error = power_error = 0.0
    for row, state in zip(rows, states):
        assert all(math.isfinite(value) for value in row.values())
        assert all(math.isfinite(value) for value in state.values())
        assert abs(row['t'] - state['time']) < 1e-12
        v, src, dc, e = (row['DCLink_dc_' + name] for name in ('vdc', 'isrc', 'idc', 'energy'))
        pac = sum(row['Converter_bridge_vo' + p] * row['DependentVoltageSource_filter_i' + p] for p in 'abc')
        power_error = max(power_error, abs(v * dc - pac))
        residual_error = max(residual_error, abs(capacitance * state['yp:dc[0]'] - (src - dc)))
        assert abs(dc - v / resistance) < 1e-8
        assert abs(e - 0.5 * capacitance * v * v) < 1e-8
        assert abs(v - state['y:dc[0]']) < 1e-10
        if abs(row['t'] - event) > 1e-12:
            assert src == (80.0 if row['t'] < event else 20.0)
        voltage.append(v)
        source.append(src)
        current.append(dc)
        energy.append(e)
        power.append(v * src - pac)
    error = max(abs(a - b) for a, b in zip(voltage, exact))
    assert error < 2e-6, error
    assert power_error < 1e-8, power_error
    assert residual_error < 2e-5, residual_error
    event_indices = [i for i, r in enumerate(rows) if abs(r['t'] - event) < 1e-12]
    assert len(event_indices) == 2
    left, right = event_indices
    assert source[left] == 80.0 and source[right] == 20.0
    assert abs(voltage[right] - voltage[left]) < 1e-10
    assert abs(states[right]['yp:dc[0]'] - states[left]['yp:dc[0]'] + 3000.0) < 2e-5
    # Trapezoidal quadrature never crosses the discontinuity: its two samples
    # have equal time, but distinct left/right power values.
    integrated = [energy[0]]
    for i in range(1, len(rows)):
        integrated.append(integrated[-1] + 0.5 * (power[i-1] + power[i]) * (rows[i]['t'] - rows[i-1]['t']))
    energy_error = max(abs(a - b) for a, b in zip(energy, integrated))
    assert energy_error < 0.02, energy_error  # O(dt_monitor^2) quadrature error, in joules.
    metrics = {'samples': len(rows), 'maximum_voltage_error_V': error,
               'maximum_converter_power_error_W': power_error,
               'maximum_capacitor_residual_A': residual_error,
               'maximum_integrated_energy_error_J': energy_error,
               'event_voltage_V': voltage[left], 'final_voltage_V': voltage[-1]}
    (output / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
    print(json.dumps(metrics, indent=2))
    if plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        time = [r['t'] for r in rows]
        fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True, layout='constrained')
        axes[0].plot(time, voltage, label='IDA simulation', linewidth=2)
        axes[0].plot(time[::60], exact[::60], 'o', fillstyle='none', markersize=4, label='Exact RC solution')
        axes[0].set_ylabel('DC voltage (V)')
        axes[0].set_title('DC-link charging and discharge — fixed bridge state (1, 0, 0)')
        axes[1].plot(time, source, label='Source current')
        axes[1].plot(time, current, label='Converter DC current')
        axes[1].set_ylabel('Current (A)')
        axes[2].plot(time, [e / 1000 for e in energy], label='Capacitor energy')
        axes[2].plot(time[::60], [e / 1000 for e in integrated[::60]], 'o', fillstyle='none', markersize=4, label='Integrated net power')
        axes[2].set_ylabel('Energy (kJ)')
        axes[2].set_xlabel('Time (s)')
        for ax in axes:
            ax.axvline(event, color='0.4', linestyle=':', linewidth=1)
            ax.grid(alpha=0.2)
            ax.legend(loc='best')
        for extension in ('png', 'svg'):
            fig.savefig(output / ('dc_link.' + extension), dpi=170)
        plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--results', type=Path)
    parser.add_argument('--plot', action='store_true', help='Requires matplotlib')
    args = parser.parse_args()
    if args.results:
        validate(args.exe.resolve(), args.results.resolve(), args.plot)
    else:
        with tempfile.TemporaryDirectory(prefix='gridkit-dc-link-') as directory:
            validate(args.exe.resolve(), Path(directory), args.plot)
