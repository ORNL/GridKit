#!/usr/bin/env python3
"""Run the current-control example at a selected PWM resolution."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

HERE = Path(__file__).resolve().parent


def run(name, exe, output, overrides):
    solver = json.loads((HERE / f'{name}.solver.json').read_text())
    solver.update({key: value for key, value in overrides.items() if value is not None})
    solver['events'] = [event for event in solver['events'] if event['time'] <= solver['tmax']]
    solver['output_file'] = f'{name}.csv'
    inputs = []
    for key in ['system_model_file', 'state_file']:
        source = (HERE / solver[key]).resolve()
        destination = output / source.name
        if source == destination:
            raise ValueError('The output directory must differ from the input directory')
        shutil.copyfile(source, destination)
        solver[key] = destination.name
        inputs.append(destination)
    path = output / f'{name}.solver.json'
    path.write_text(json.dumps(solver, indent=2) + '\n')
    inputs.append(path)
    record_path = output / f'{name}.run.json'
    record_path.unlink(missing_ok=True)
    command = [str(exe), str(path)]
    start = time.perf_counter()
    log_path = output / f'{name}.log'
    with log_path.open('w') as log:
        subprocess.run(command, cwd=output, stdout=log, stderr=subprocess.STDOUT, check=True)
    wall_time = time.perf_counter() - start
    log = log_path.read_text()
    statistics = re.search(r'IDA statistics: (.+)', log).group(1)
    record = {
        'command': command,
        'executable_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
        'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'wall_time_s': wall_time,
        'cpu_time_s': float(re.search(r'Complete in (\S+) seconds', log).group(1)),
        'ida': {key: int(value) for key, value in re.findall(r'(\w+)=(\d+)', statistics)},
    }
    record_path.write_text(json.dumps(record, indent=2) + '\n')
    print(f'{name}: {output / solver["output_file"]} ({wall_time:.2f} s)', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, default=HERE.parents[2] / 'build/application/EMT/EMTDynamicSimulation')
    parser.add_argument('--output', type=Path, default=Path('simulation'), help='Output directory, relative to this example')
    parser.add_argument('--mu', type=float, help='PWM smoothing sharpness [1/s]')
    parser.add_argument('--tmax', type=float, help='Final time [s]')
    parser.add_argument('--dt-monitor', type=float, help='Monitor interval [s]')
    args = parser.parse_args()
    overrides = {'mu': args.mu, 'tmax': args.tmax, 'dt_monitor': args.dt_monitor}
    if any(value is not None and (not math.isfinite(value) or value <= 0) for value in overrides.values()):
        parser.error('mu, tmax, and dt-monitor must be finite and positive')
    output = (HERE / args.output).resolve()
    if output == HERE:
        parser.error('Use a separate output directory')
    output.mkdir(parents=True, exist_ok=True)
    run('GFL', args.exe.resolve(), output, overrides)


if __name__ == '__main__':
    main()
