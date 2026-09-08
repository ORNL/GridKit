"""Build editable PGFPlots comparisons against the frozen PowerWorld reference."""

import argparse
import bisect
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent
REVISION = 'f16f3815e6f84c08bd0869773d3033b3dcb3dd0f'
REFERENCE = 'examples/PhasorDynamics/Validation/Hawaii/reference'
LABELS = {'vmag': r'$|V|$ [p.u.]', 'omega': r'$\omega_r-1$ [p.u.]',
          'p': r'$P$ [p.u., 100 MVA]', 'q': r'$Q$ [p.u., 100 MVA]'}


def interpolate(time, values, point):
    k = min(max(bisect.bisect_right(time, point) - 1, 0), len(time) - 2)
    if time[k + 1] == time[k]:
        return values[k + 1]
    fraction = (point - time[k]) / (time[k + 1] - time[k])
    return values[k] + fraction * (values[k + 1] - values[k])


def document(kind, names, lower, upper, error_lower, error_upper, style_path, duration):
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
    for panel, color, title, clear in [('e', 'emtBlue', '(a) EMT: switching GFL plants; fault cleared at 1.10 s', 1.1),
                                      ('r', 'referenceOrange', '(b) PowerWorld reference; fault cleared at 1.15 s', 1.15)]:
        legend = (r',legend style={at={(0.02,0.98)},anchor=north west,legend columns=3,fill=white,draw=none,font=\scriptsize}'
                  if kind != 'vmag' else '')
        lines.append(f'\\nextgroupplot[title={{{title}}},ylabel={{{LABELS[kind]}}},ymin=\\ylower,ymax=\\yupper{legend}]')
        lines.append(f'\\path[fill=gray!15] (axis cs:1,\\ylower) rectangle (axis cs:{clear},\\yupper);')
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
        label = 'All 37 buses; bus 1 in black' if kind == 'vmag' else 'All 30 synchronous machines'
        lines.append(f'\\node[anchor=south west,font=\\footnotesize,fill=white,inner sep=2pt] at (rel axis cs:0.01,0.02) {{{label}}};')
    lines.append(f'\\nextgroupplot[title={{(c) EMT minus PowerWorld: range across channels}},ylabel={{Difference [p.u.]}},ymin={error_lower:.8g},ymax={error_upper:.8g}]')
    lines.extend([
        f'\\addplot[name path=lo,draw=none] table[x=time,y=minimum,col sep=comma] {{Hawaii.{kind}.csv}};',
        f'\\addplot[name path=hi,draw=none] table[x=time,y=maximum,col sep=comma] {{Hawaii.{kind}.csv}};',
        r'\addplot[emtBlue!25] fill between[of=lo and hi];',
        f'\\addplot[emtBlue,line width=0.5pt] table[x=time,y=minimum,col sep=comma] {{Hawaii.{kind}.csv}};',
        f'\\addplot[emtBlue,line width=0.5pt] table[x=time,y=maximum,col sep=comma] {{Hawaii.{kind}.csv}};',
        f'\\addplot[gray,dashed] coordinates {{(0,0) ({duration:g},0)}};',
        r'\end{groupplot}', r'\end{tikzpicture}', r'\end{document}', ''])
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--averaged', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT / 'results')
    parser.add_argument('--no-render', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    run_metrics = json.loads(args.averaged.with_name('metrics.json').read_text())
    averaged_hash = hashlib.sha256(args.averaged.read_bytes()).hexdigest()
    if averaged_hash != run_metrics['averaged_sha256']:
        raise ValueError('Averaged data do not match the validated run')
    if run_metrics['final_time_s'] < 1.5 or run_metrics['source_revision'] != REVISION:
        raise ValueError('Comparison requires fault recovery through 1.5 s at the frozen source revision')
    with args.averaged.open(newline='') as stream:
        emt = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    metrics = {'source_revision': REVISION, 'averaged_sha256': averaged_hash,
               'input_sha256': run_metrics['input_sha256'],
               'description': 'Descriptive errors between different models and fault clearing times; not validation tolerances.',
               'emt_fault_interval_s': [1.0, 1.1], 'reference_fault_interval_s': [1.0, 1.15], 'channels': {}}
    for kind in LABELS:
        raw = subprocess.check_output(['git', 'show', f'{REVISION}:{REFERENCE}/Hawaii.{kind}.ref.csv'], cwd=ROOT).decode()
        reference = list(csv.DictReader(io.StringIO(raw)))
        names = [name for name in reference[0] if name != 'time']
        time = [float(row['time']) for row in reference]
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
        # Retain 240 Hz around the fault and 60 Hz elsewhere in the editable plot data.
        selected = [k for k, row in enumerate(emt) if 0.9 <= row['time'] <= 1.4 or k % 4 == 0 or k == len(emt) - 1]
        with (output / f'Hawaii.{kind}.csv').open('w', newline='') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(['time'] + [f'e{k}' for k in range(len(names))]
                            + [f'r{k}' for k in range(len(names))] + ['minimum', 'maximum'])
            for k in selected:
                writer.writerow([f'{emt[k]["time"]:.9g}'] + [f'{x:.9g}' for x in actual[k] + aligned[k]]
                                + [f'{min(differences[k]):.9g}', f'{max(differences[k]):.9g}'])
        low = min(x for row in actual + aligned for x in row)
        high = max(x for row in actual + aligned for x in row)
        margin = max(1e-6, 0.07 * (high - low))
        elow, ehigh = min(min(row) for row in differences), max(max(row) for row in differences)
        emargin = max(1e-6, 0.1 * (ehigh - elow))
        name = f'Hawaii.{kind}'
        style = os.path.relpath(ROOT.parents[2] / 'docs/Figures/EMT/diagram-style.tex', output)
        (output / (name + '.tex')).write_text(document(kind, names, low - margin, high + margin,
                                                     elow - emargin, ehigh + emargin, style, run_metrics['final_time_s']))
        if not args.no_render:
            with tempfile.TemporaryDirectory(prefix='gridkit-hawaii-plot-') as temporary:
                command = ['pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-output-directory=' + temporary, name + '.tex']
                for _ in range(2):
                    result = subprocess.run(command, cwd=output, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    if result.returncode:
                        raise RuntimeError(result.stdout.decode()[-5000:])
                (output / (name + '.pdf')).write_bytes((Path(temporary) / (name + '.pdf')).read_bytes())
            subprocess.run(['pdftoppm', '-png', '-r', '600', '-singlefile', name + '.pdf', name], cwd=output, check=True)
    (output / 'comparison.json').write_text(json.dumps(metrics, indent=2) + '\n')


if __name__ == '__main__':
    main()
