"""Build editable PGFPlots comparisons against a fresh GridKit PhasorDynamics run."""

import argparse
import bisect
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent
LABELS = {'vmag': r'$|V|$ [p.u.]', 'omega': r'$\omega_r-1$ [p.u.]',
          'p': r'$P$ [p.u., 100 MVA]', 'q': r'$Q$ [p.u., 100 MVA]'}


def interpolate(time, values, point):
    k = min(max(bisect.bisect_right(time, point) - 1, 0), len(time) - 2)
    if time[k + 1] == time[k]:
        return values[k + 1]
    fraction = (point - time[k]) / (time[k + 1] - time[k])
    return values[k] + fraction * (values[k + 1] - values[k])


def document(kind, names, lower, upper, style_path, duration,
             emt_fault, phasor_fault, mu):
    """Use the EMT shared fonts/styles; data and size remain directly editable."""
    # All channels are retained. Highlight fault bus 1; other channels share a
    # neutral style so the comparison is readable without a 37-entry legend.
    count = len(names)
    machine_buses = sorted({name.split()[0] for name in names}, key=int)
    palette = ['emtBlue', 'referenceOrange', 'busGreen', 'busPurple', 'busSky', 'black']
    lines = [r'\documentclass[border=3pt]{standalone}',
             r'\usepackage{tikz}', r'\usepackage{pgfplots}',
             r'\input{' + style_path + '}',
             r'\usepgfplotslibrary{groupplots,fillbetween}', r'\pgfplotsset{compat=1.18}',
             r'\definecolor{emtBlue}{RGB}{0,114,178}', r'\definecolor{referenceOrange}{RGB}{213,94,0}',
             r'\definecolor{busGreen}{RGB}{0,158,115}', r'\definecolor{busPurple}{RGB}{204,121,167}',
             r'\definecolor{busSky}{RGB}{86,180,233}',
             '% TUNABLES: panel width, height, and shared axis limits',
             r'\def\figwidth{183mm}', r'\def\panelheight{49mm}',
             f'\\def\\ylower{{{lower:.8g}}}\\def\\yupper{{{upper:.8g}}}',
             r'\begin{document}', r'\begin{tikzpicture}[font=\small]',
             r'\begin{groupplot}[group style={group size=1 by 3,vertical sep=11mm,x descriptions at=edge bottom},',
             f'width=\\figwidth,height=\\panelheight,xmin=0,xmax={duration:g},grid=major,grid style={{gray!20}},',
             r'tick label style={font=\footnotesize},scaled y ticks=false,',
             r'xlabel={Time [s]},ylabel style={font=\small},legend style={draw=none,font=\footnotesize}]']
    for panel, color, title, fault in [
            ('e', 'emtBlue', f'(a) GridKit EMT ($\\mu={mu:g}$); fault cleared at {emt_fault[1]:.2f} s', emt_fault),
            ('r', 'referenceOrange', f'(b) GridKit PhasorDynamics; fault cleared at {phasor_fault[1]:.2f} s', phasor_fault)]:
        if kind != 'vmag':
            title += f' ({count} machines)'
        legend = (r',legend style={at={(0.98,0.98)},anchor=north east,legend columns=3,fill=white,draw=none,font=\scriptsize}'
                  if kind != 'vmag' else '')
        lines.append(f'\\nextgroupplot[title={{{title}}},ylabel={{{LABELS[kind]}}},ymin=\\ylower,ymax=\\yupper{legend}]')
        lines.append(f'\\path[fill=gray!15] (axis cs:{fault[0]},\\ylower) rectangle (axis cs:{fault[1]},\\yupper);')
        order = list(range(1, count)) + [0] if kind == 'vmag' else range(count)
        for k in order:
            style = f'{color}!65,line width=0.35pt'
            if kind == 'vmag' and k == 0:
                style = 'black,line width=0.8pt'
            elif kind != 'vmag':
                bus = names[k].split()[0]
                style = palette[machine_buses.index(bus)] + ',line width=0.45pt'
            lines.append(f'\\addplot[{style},no marks,forget plot] table[x=time,y={panel}{k},col sep=comma] {{Hawaii.{kind}.csv}};')
        if kind != 'vmag':
            for k, bus in enumerate(machine_buses):
                lines.extend([f'\\addlegendimage{{{palette[k]},line width=0.7pt}}', f'\\addlegendentry{{Bus {bus}}}'])
        if kind == 'vmag':
            lines.append(f'\\node[anchor=south east,font=\\footnotesize,fill=white,inner sep=2pt] at (rel axis cs:0.99,0.02) {{All {count} buses; bus 1 in black}};')
    lines.append(r'\nextgroupplot[title={(c) EMT minus PhasorDynamics: range across channels},ylabel={Difference [p.u.]},ymin=\ylower,ymax=\yupper]')
    lines.extend([
        f'\\addplot[name path=lo,draw=none] table[x=time,y=minimum,col sep=comma] {{Hawaii.{kind}.csv}};',
        f'\\addplot[name path=hi,draw=none] table[x=time,y=maximum,col sep=comma] {{Hawaii.{kind}.csv}};',
        r'\addplot[emtBlue!25] fill between[of=lo and hi];',
        f'\\addplot[emtBlue,line width=0.5pt] table[x=time,y=minimum,col sep=comma] {{Hawaii.{kind}.csv}};',
        f'\\addplot[emtBlue,line width=0.5pt] table[x=time,y=maximum,col sep=comma] {{Hawaii.{kind}.csv}};',
        f'\\addplot[gray,dashed] coordinates {{(0,0) ({duration:g},0)}};',
        r'\end{groupplot}', r'\end{tikzpicture}', r'\end{document}', ''])
    return '\n'.join(lines)


def render(output, name):
    command = ['pdflatex', '-cnf-line=extra_mem_top=20000000', '-cnf-line=extra_mem_bot=20000000',
               '-interaction=nonstopmode', '-halt-on-error', name + '.tex']
    for _ in range(2):
        result = subprocess.run(command, cwd=output, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(result.stdout.decode()[-5000:])
    subprocess.run(['pdftoppm', '-png', '-r', '600', '-singlefile', name + '.pdf', name], cwd=output, check=True)


def waveform_document(output, name, panels, mu, fault, columns=3, time_unit='s'):
    """Raw samples, with individual CSV tables and no RMS or cycle averaging."""
    style = os.path.relpath(ROOT.parents[2] / 'docs/Figures/EMT/diagram-style.tex', output)
    rows = math.ceil(len(panels) / columns)
    width = 58 if columns == 3 else 183
    lines = [r'\documentclass[border=3pt]{standalone}', r'\usepackage{pgfplots}',
             r'\usepgfplotslibrary{groupplots}', r'\pgfplotsset{compat=1.18}',
             r'\input{' + style + '}',
             r'\definecolor{phaseA}{HTML}{0072B2}', r'\definecolor{phaseB}{HTML}{D55E00}',
             r'\definecolor{phaseC}{HTML}{009E73}', r'\definecolor{commandQ}{HTML}{CC79A7}',
             '% TUNABLES: panel dimensions; CSV files contain instantaneous samples.',
             rf'\def\panelwidth{{{width}mm}}', r'\def\panelheight{47mm}',
             r'\begin{document}', r'\begin{tikzpicture}[font=\small]',
             rf'\begin{{groupplot}}[group style={{group size={columns} by {rows},horizontal sep=13mm,vertical sep=15mm}},',
             r'width=\panelwidth,height=\panelheight,grid=major,grid style={gray!15},',
             rf'xlabel={{Time [{time_unit}]}},tick label style={{font=\scriptsize}},title style={{font=\scriptsize}},',
             r'ylabel style={font=\scriptsize},scaled ticks=false,',
             r'legend style={font=\tiny,draw=none,fill=white,legend columns=2,cells={anchor=west}}]']
    for index, (title, ylabel, times, curves) in enumerate(panels):
        filename = f'{name}.{index}.csv'
        with (output / filename).open('w', newline='') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(['time'] + [f'y{k}' for k in range(len(curves))])
            writer.writerows([f'{t:.10g}'] + [f'{curve[2][n]:.9g}' for curve in curves]
                             for n, t in enumerate(times))
        title = title.replace('_', r'\_')
        lines.append(rf'\nextgroupplot[title={{{title}; $\mu={mu:g}$}},ylabel={{{ylabel}}},xmin={times[0]:.10g},xmax={times[-1]:.10g}]')
        for k, (label, style, _) in enumerate(curves):
            lines.append(rf'\addplot[{style},line width=0.5pt,no marks] table[x=time,y=y{k},col sep=comma] {{{filename}}};')
            if index == 0:
                lines.append(r'\addlegendentry{' + label + '}')
        for event in fault:
            if times[0] <= event <= times[-1]:
                x = (event - times[0]) / (times[-1] - times[0])
                lines.append(rf'\draw[gray,densely dashed] (rel axis cs:{x:.9g},0) -- (rel axis cs:{x:.9g},1);')
    lines.extend([r'\end{groupplot}', r'\end{tikzpicture}', r'\end{document}'])
    (output / (name + '.tex')).write_text('\n'.join(lines) + '\n')


def waveforms(run, output, metrics, switching):
    study = metrics['study']
    case = json.loads((run / 'Hawaii.case.json').read_text())
    report = json.loads((run / 'conversion.json').read_text())
    plants = list(report['inverters'])
    filters = {d['id'].removesuffix('_filter'): d for d in case['devices'] if d['class'] == 'Filter'}
    fault = sorted({e['time'] for e in study['events'] if e['element_id'] == 'fault_switch'})
    windows = [(fault[0] - .03, fault[0] + .03), (fault[1] - .03, fault[1] + .03),
               (study['tmax'] - .05, study['tmax'])]
    buses = ['1', '23', '33']
    selected_plants = ['23_10', '26_1', '33_1']
    keys = {f'Bus_bus_{bus}_v{p}' for bus in buses for p in 'abc'}
    for plant in plants:
        keys.update(f'Park_{plant}_grid_current_y{k}' for k in [1, 2])
        keys.update(f'Reecb_{plant}_electrical_icmd{q}' for q in 'dq')
        keys.add(f'PLL_{plant}_pll_omega')
        keys.update(f'Filter_{plant}_filter_ig{p}' for p in 'abc')
        keys.update(f'Bus_{filters[plant]["inputs"]["bus"]}_v{p}' for p in 'abc')
    for plant in selected_plants:
        keys.update(f'Filter_{plant}_filter_i{p}' for p in 'abc')
    raw, overview = [], []
    with (run / study['output_file']).open(newline='') as stream:
        for index, row in enumerate(csv.DictReader(stream)):
            t = float(row['t'])
            detailed = any(lo <= t <= hi for lo, hi in windows)
            sampled = index % (6 if .95 <= t <= 1.3 else 24) == 0 or t == study['tmax']
            if not (detailed or sampled):
                continue
            point = {'t': t, **{k: float(row[k]) for k in keys}}
            if detailed:
                raw.append(point)
            if sampled:
                overview.append(point)
    if abs(overview[-1]['t'] - study['tmax']) > 1e-9:
        raise ValueError('Waveform run is incomplete')
    phase_styles = ['phaseA', 'phaseB', 'phaseC']
    names = []
    panels = []
    for bus in buses:
        for lo, hi in windows:
            rows = [row for row in raw if lo <= row['t'] <= hi]
            curves = [(f'$v_{p}$', style, [row[f'Bus_bus_{bus}_v{p}'] / 1000 for row in rows])
                      for p, style in zip('abc', phase_styles)]
            panels.append((f'Bus {bus}', 'Voltage [kV]', [row['t'] for row in rows], curves))
    waveform_document(output, 'Hawaii.voltage-abc', panels, study['mu'], fault)
    names.append('Hawaii.voltage-abc')
    panels = []
    for plant in selected_plants:
        for lo, hi in windows:
            rows = [row for row in raw if lo <= row['t'] <= hi]
            curves = [(f'$i_{{{p}}}$', style, [row[f'Filter_{plant}_filter_i{p}'] / 1000 for row in rows])
                      for p, style in zip('abc', phase_styles)]
            curves += [(f'$i_{{g,{p}}}$', style + ',dashed', [row[f'Filter_{plant}_filter_ig{p}'] / 1000 for row in rows])
                       for p, style in zip('abc', phase_styles)]
            panels.append((f'Plant {plant}', 'Current [kA]', [row['t'] for row in rows], curves))
    waveform_document(output, 'Hawaii.current-abc', panels, study['mu'], fault)
    names.append('Hawaii.current-abc')
    for kind in ['current-dq', 'inverter-pq', 'pll']:
        panels = []
        for plant in plants:
            times = [row['t'] for row in overview]
            if kind == 'current-dq':
                curves = [(f'$i_{{g,{q}}}$', style, [row[f'Park_{plant}_grid_current_y{k}'] / 1000 for row in overview])
                          for q, k, style in [('d', 1, 'phaseA'), ('q', 2, 'phaseB')]]
                curves += [(rf'$i_{{g,{q}}}^{{\mathrm{{cmd}}}}$', style, [row[f'Reecb_{plant}_electrical_icmd{q}'] / 1000 for row in overview])
                           for q, style in [('d', 'black,dashed'), ('q', 'commandQ,dashed')]]
                ylabel = 'Current [kA]'
            elif kind == 'inverter-pq':
                p, q = [], []
                for row in overview:
                    v = [row[f'Bus_{filters[plant]["inputs"]["bus"]}_v{x}'] for x in 'abc']
                    i = [row[f'Filter_{plant}_filter_ig{x}'] for x in 'abc']
                    p.append(sum(a * b for a, b in zip(v, i)) / 1e6)
                    q.append(((v[1]-v[2])*i[0] + (v[2]-v[0])*i[1] + (v[0]-v[1])*i[2]) / math.sqrt(3) / 1e6)
                curves = [('$p$', 'phaseA', p), ('$q$', 'phaseB', q)]
                ylabel = 'Power [MW / Mvar]'
            else:
                curves = [(r'$f_{\mathrm{PLL}}$', 'phaseA', [row[f'PLL_{plant}_pll_omega'] / (2 * math.pi) for row in overview])]
                ylabel = 'Frequency [Hz]'
            panels.append((f'Plant {plant}', ylabel, times, curves))
        name = 'Hawaii.' + kind
        waveform_document(output, name, panels, study['mu'], fault)
        names.append(name)
    if switching:
        settings = json.loads((switching / 'study.json').read_text())
        if settings['mu'] != study['mu']:
            raise ValueError('Switching-detail mu differs from the fault run')
        with (switching / 'switching.csv').open(newline='') as stream:
            samples = list(csv.DictReader(stream))
        plant = selected_plants[0]
        end = settings['tmax']
        samples = [r for r in samples if float(r['t']) >= end - 6 / report['choices']['carrier_Hz']]
        times = [1000 * float(r['t']) for r in samples]
        panels = []
        for prefix, quantity, label, scale in [('PWM', 's', 'Switching function [--]', 1), ('Converter', 'e', 'Bridge voltage [kV]', 1000)]:
            device = 'pwm' if prefix == 'PWM' else 'bridge'
            curves = [(f'${quantity}_{p}$', style, [float(r[f'{prefix}_{plant}_{device}_{quantity}{p}']) / scale for r in samples])
                      for p, style in zip('abc', phase_styles)]
            panels.append((f'Plant {plant}: startup switching detail', label, times, curves))
        waveform_document(output, 'Hawaii.switching', panels, study['mu'], [], columns=1, time_unit='ms')
        names.append('Hawaii.switching')
    return names



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--averaged', type=Path, required=True)
    parser.add_argument('--reference', type=Path, default=ROOT / 'phasor-reference/monitored')
    parser.add_argument('--output', type=Path, default=ROOT / 'results')
    parser.add_argument('--no-render', action='store_true')
    parser.add_argument('--waveforms', action='store_true', help='Add instantaneous phase and inverter-control figures')
    parser.add_argument('--switching', type=Path, help='Completed high-resolution switching.py run directory')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    run_metrics = json.loads(args.averaged.with_name('metrics.json').read_text())
    averaged_hash = hashlib.sha256(args.averaged.read_bytes()).hexdigest()
    if averaged_hash != run_metrics['averaged_sha256']:
        raise ValueError('Averaged data do not match the validated run')
    reference_metrics = json.loads((args.reference / 'metrics.json').read_text())
    if run_metrics['exciter_adjustments'] != reference_metrics['exciter_adjustments']:
        raise ValueError('EMT and PhasorDynamics exciter adjustments differ')
    if run_metrics['source_revision'] != reference_metrics['source']['revision']:
        raise ValueError('EMT conversion and PhasorDynamics case must use the same source revision')
    if run_metrics['final_time_s'] < 1.5:
        raise ValueError('Comparison requires fault recovery through 1.5 s')
    with args.averaged.open(newline='') as stream:
        emt = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    emt_fault = [event['time'] for event in run_metrics['study']['events'] if event['element_id'] == 'fault_switch']
    phasor_fault = [event['time'] for event in reference_metrics['study']['events']]
    metrics = {'source_revision': run_metrics['source_revision'], 'averaged_sha256': averaged_hash,
               'input_sha256': run_metrics['input_sha256'],
               'reference': reference_metrics,
               'description': 'Descriptive differences between GridKit EMT and PhasorDynamics simulations.',
               'emt_fault_interval_s': emt_fault, 'reference_fault_interval_s': phasor_fault, 'channels': {}}
    figures = []
    for kind in LABELS:
        path = args.reference / f'Hawaii.{kind}.csv'
        if hashlib.sha256(path.read_bytes()).hexdigest() != reference_metrics['channels_sha256'][path.name]:
            raise ValueError(f'PhasorDynamics data changed since simulation: {path}')
        with path.open(newline='') as stream:
            reference = list(csv.DictReader(stream))
        names = [name for name in reference[0] if name != 'time']
        time = [float(row['time']) for row in reference]
        if time[0] > emt[0]['time'] or time[-1] < emt[-1]['time']:
            raise ValueError('PhasorDynamics simulation must cover the EMT comparison interval')
        values = {name: [float(row[name]) for row in reference] for name in names}
        aligned = [[interpolate(time, values[name], row['time']) for name in names] for row in emt]
        actual = [[row[f'{kind}:{name}'] for name in names] for row in emt]
        differences = [[a - b for a, b in zip(a_row, b_row)] for a_row, b_row in zip(actual, aligned)]
        channel_metrics = {}
        for k, name in enumerate(names):
            errors = [row[k] for row in differences]
            channel_metrics[name] = {'rmse_pu': math.sqrt(sum(e * e for e in errors) / len(errors)),
                                     'maximum_absolute_error_pu': max(map(abs, errors)),
                                     'final_error_pu': errors[-1]}
        metrics['channels'][kind] = channel_metrics
        metrics[kind] = {
            'maximum_absolute_error_pu': max(abs(e) for row in differences for e in row),
            'rmse_pu': math.sqrt(sum(e * e for row in differences for e in row) / (len(emt) * len(names))),
        }
        # Retain every available cycle average around the fault.
        selected = [k for k, row in enumerate(emt) if 0.9 <= row['time'] <= 1.4 or k % 4 == 0 or k == len(emt) - 1]
        with (output / f'Hawaii.{kind}.csv').open('w', newline='') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(['time'] + [f'e{k}' for k in range(len(names))]
                            + [f'r{k}' for k in range(len(names))] + ['minimum', 'maximum'])
            for k in selected:
                writer.writerow([f'{emt[k]["time"]:.9g}'] + [f'{x:.9g}' for x in actual[k] + aligned[k]]
                                + [f'{min(differences[k]):.9g}', f'{max(differences[k]):.9g}'])
        low = min(x for row in actual + aligned + differences for x in row)
        high = max(x for row in actual + aligned + differences for x in row)
        margin = max(1e-6, 0.07 * (high - low))
        name = f'Hawaii.{kind}'
        style = os.path.relpath(ROOT.parents[2] / 'docs/Figures/EMT/diagram-style.tex', output)
        (output / (name + '.tex')).write_text(document(kind, names, low - margin, high + margin,
                                                     style, run_metrics['final_time_s'],
                                                     emt_fault, phasor_fault, run_metrics['study']['mu']))
        figures.append(name)
    if args.waveforms:
        figures.extend(waveforms(args.averaged.resolve().parent, output, run_metrics, args.switching))
    if not args.no_render:
        with ThreadPoolExecutor(max_workers=4) as pool:
            list(pool.map(lambda name: render(output, name), figures))
    (output / 'comparison.json').write_text(json.dumps(metrics, indent=2) + '\n')


if __name__ == '__main__':
    main()
