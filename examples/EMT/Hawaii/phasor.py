"""Run the copied Hawaii PhasorDynamics validation and export comparison traces."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[2]
BRANCH = 'lukel/cases-polish-dev'
VALIDATION = 'examples/PhasorDynamics/Validation/Hawaii'
CASE = 'cases/PhasorDynamics/Hawaii/Hawaii.case.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def export(directory, revision):
    paths = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', revision, VALIDATION],
                                    cwd=REPOSITORY, text=True).splitlines()
    paths = [path for path in paths if '/figures/' not in path] + [CASE]
    hashes = {}
    for path in paths:
        name = 'Hawaii.case.json' if path == CASE else str(Path(path).relative_to(VALIDATION))
        target = directory / name
        raw = subprocess.check_output(['git', 'show', f'{revision}:{path}'], cwd=REPOSITORY)
        if target.exists() and target.read_bytes() != raw:
            raise ValueError(f'Copied source differs: {target}; use a new output directory')
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
        hashes[name] = digest(target)
    return hashes


def run(exe, directory):
    command = [str(exe), 'Hawaii.solver.json']
    started = time.perf_counter()
    with (directory / 'simulation.log').open('w') as log:
        subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True)
    wall = time.perf_counter() - started
    log = (directory / 'simulation.log').read_text()
    result = {'command': command, 'wall_time_s': wall,
              'cpu_time_s': float(re.search(r'Complete in (\S+) seconds', log).group(1)),
              'input_sha256': {name: digest(directory / name)
                               for name in ('Hawaii.case.json', 'Hawaii.solver.json')}}
    write_json(directory / 'run.json', result)
    return result


def read_csv(path):
    with path.open(newline='') as stream:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    if not rows or not all(math.isfinite(value) for row in rows for value in row.values()):
        raise ValueError(f'Empty or nonfinite simulation: {path}')
    return rows


def export_channels(directory, case):
    rows = read_csv(directory / 'Hawaii.csv')
    channels = {'vmag': {str(bus['number']): f"Bus_{bus['name']}_Vm" for bus in case['buses']}}
    machines = [device for device in case['devices'] if device['class'] == 'Genrou']
    for kind in ('omega', 'p', 'q'):
        channels[kind] = {device['id'].removesuffix('_genrou').replace('_', ' '):
                          f"Genrou_{device['id']}_{kind}" for device in machines}
    hashes = {}
    for kind, columns in channels.items():
        path = directory / f'Hawaii.{kind}.csv'
        with path.open('w', newline='') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(['time', *columns])
            writer.writerows([row['t'], *(row[key] for key in columns.values())] for row in rows)
        hashes[path.name] = digest(path)
    return rows, hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT / 'phasor-reference')
    args = parser.parse_args()
    directory, exe = args.output.resolve(), args.exe.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    revision = subprocess.check_output(['git', 'rev-parse', BRANCH], cwd=REPOSITORY, text=True).strip()
    source_hashes = export(directory, revision)
    source = {'branch': BRANCH, 'revision': revision, 'input_sha256': source_hashes,
              'executable_sha256': digest(exe)}
    write_json(directory / 'source.json', source)

    # This is exactly the branch's Hawaii_validation command and input pair.
    validation = run(exe, directory)
    original = read_csv(directory / 'Hawaii.omega.csv')
    case = json.loads((directory / 'Hawaii.case.json').read_text())
    study = json.loads((directory / 'Hawaii.solver.json').read_text())
    for bus in case['buses']:
        bus['mon'] = ['Vm']
    for device in case['devices']:
        if device['class'] == 'Genrou':
            device['mon'] = ['omega', 'p', 'q']
    study['output_file'] = 'Hawaii.csv'
    for key in ('reference_file', 'error_type', 'error_tolerance', 'abs_err_threshold'):
        study.pop(key, None)
    # Added monitors preserve the source physical model and solver settings.
    target = directory / 'monitored'
    target.mkdir(exist_ok=True)
    write_json(target / 'Hawaii.case.json', case)
    write_json(target / 'Hawaii.solver.json', study)
    record = run(exe, target)
    rows, hashes = export_channels(target, case)
    if len(original) != len(rows):
        raise ValueError('Adding monitors changed the sample count')
    difference = max(abs(row[key] - ref[key]) for row, ref in zip(rows, original) for key in ref)
    if difference > 1e-10:
        raise ValueError(f'Adding monitors changed the speed trajectory: {difference}')
    record.update(source=source, validation=validation, channels_sha256=hashes,
                  maximum_validation_trace_difference_pu=difference,
                  final_time_s=rows[-1]['t'], study=study,
                  description='Original validation physics; only monitors and output selection changed.')
    write_json(target / 'metrics.json', record)
    print(f'monitored: {record["cpu_time_s"]:.6g} CPU seconds, final time {rows[-1]["t"]:g} s')
    print(f'Unmodified Hawaii_validation passed: {validation["cpu_time_s"]:.6g} CPU seconds')


if __name__ == '__main__':
    main()
