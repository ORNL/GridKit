#!/usr/bin/env python3
"""Run the ten-bus REGFMA fault study and retain its inputs and solver statistics."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True, help='EMTDynamicSimulation executable')
    args = parser.parse_args()
    exe = args.exe.resolve()
    directory = HERE / 'results'
    directory.mkdir(exist_ok=True)
    study = json.loads((HERE / 'FaultClearing.solver.json').read_text())
    hashes = {}
    for key in ('system_model_file', 'state_file'):
        source = (HERE / study[key]).resolve()
        contents = source.read_bytes()
        (directory / source.name).write_bytes(contents)
        study[key] = source.name
        hashes[source.name] = hashlib.sha256(contents).hexdigest()
    solver = directory / 'FaultClearing.solver.json'
    solver.write_text(json.dumps(study, indent=4) + '\n')
    hashes[solver.name] = hashlib.sha256(solver.read_bytes()).hexdigest()
    (directory / 'run.json').unlink(missing_ok=True)
    start = time.perf_counter()
    with (directory / 'run.log').open('w') as log:
        subprocess.run([str(exe), solver.name], cwd=directory,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    elapsed = time.perf_counter() - start
    output = (directory / 'run.log').read_text()
    runtime = re.search(r'Complete in ([\d.eE+-]+) seconds', output)
    statistics = re.search(r'IDA statistics: (.+)', output)
    if runtime is None or statistics is None or 'SUNDIALS KLU (sparse)' not in output:
        raise ValueError('Missing successful sparse-solver run statistics; inspect results/run.log')
    stats = {name: int(value) for name, value in re.findall(r'(\w+)=(\d+)', statistics[1])}
    record = {'executable': str(exe), 'wall_seconds': elapsed,
              'simulation_cpu_seconds': float(runtime[1]), 'ida': stats,
              'input_sha256': hashes}
    (directory / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
    print(f'FaultClearing: {elapsed:.3f} s wall, {runtime[1]} s simulation CPU, {stats["steps"]} accepted steps')


if __name__ == '__main__':
    main()
