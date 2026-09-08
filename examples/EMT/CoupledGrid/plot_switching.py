"""Render one-cycle CoupledGrid switching traces at three smoothing settings."""
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess

EXAMPLE = Path(__file__).resolve().parent
HERE = EXAMPLE / 'results'
HERE.mkdir(parents=True, exist_ok=True)
ROOT = EXAMPLE.parents[2]
shutil.copyfile(EXAMPLE / 'switching.tex', HERE / 'switching.tex')
CASE = ROOT / 'cases/EMT/CoupledGrid'
PERIOD = 1 / 60
SETTINGS = [('low', 240), ('middle', 3600), ('high', 50000)]

for name in ('CoupledGrid.case.json', 'CoupledGrid.state.json'):
    shutil.copyfile(CASE / name, HERE / name)
base = json.loads((CASE / 'high.solver.json').read_text())
records = []
for name, mu in SETTINGS:
    waveform = HERE / f'mu-{mu}.csv'
    study = dict(base, mu=mu, tmax=PERIOD, dt_monitor=PERIOD / 10000,
                 system_model_file=str(HERE / 'CoupledGrid.case.json'),
                 state_file=str(HERE / 'CoupledGrid.state.json'),
                 events=[], output_file=str(waveform))
    config = HERE / f'mu-{mu}.solver.json'
    config.write_text(json.dumps(study, indent=2) + '\n')
    print(f'Plotting run: mu={mu}', flush=True)
    with (HERE / f'mu-{mu}.log').open('w') as log:
        subprocess.run([str(ROOT / 'build/application/EMT/EMTDynamicSimulation'), str(config)],
                       cwd=HERE, stdout=log, stderr=subprocess.STDOUT, check=True)
    values = []
    with waveform.open() as stream:
        for row in csv.DictReader(stream):
            values.append([1000 * float(row['t']), float(row['PWM_pwm_9_sa']),
                           float(row['Converter_converter_9_ea']) / 1000,
                           float(row['DependentVoltageSource_filter_9_ia'])])
    assert len(values) == 10001 and abs(values[-1][0] - 1000 * PERIOD) < 1e-10
    assert all(math.isfinite(x) for row in values for x in row)
    with (HERE / f'{name}.plot.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['time_ms', 's', 'e_kV', 'i_A'])
        writer.writerows(values)
    records.append({'mu': mu, 'dt_monitor_s': study['dt_monitor'], 'samples': len(values),
                    'gate_range': [min(row[1] for row in values), max(row[1] for row in values)],
                    'waveform': waveform.name})
(HERE / 'switching-runs.json').write_text(json.dumps(records, indent=2) + '\n')
for _ in range(2):
    with (HERE / 'latex.log').open('w') as log:
        subprocess.run(['pdflatex', '-cnf-line=extra_mem_top=20000000', '-cnf-line=extra_mem_bot=20000000', '-interaction=nonstopmode', '-halt-on-error', 'switching.tex'],
                       cwd=HERE, stdout=log, stderr=subprocess.STDOUT, check=True)
subprocess.run(['pdftoppm', '-png', '-r', '180', '-singlefile', 'switching.pdf', 'switching'], cwd=HERE, check=True)
print(HERE / 'switching.pdf', flush=True)
