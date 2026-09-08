#!/usr/bin/env python3
"""Run the two switching-control examples without invoking the test suite."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, default=here.parents[2] / 'build/application/EMT/EMTDynamicSimulation')
    parser.add_argument('--scenario', choices=['GFL', 'GFM'])
    args = parser.parse_args()
    output = here / 'simulation'
    output.mkdir(exist_ok=True)
    for name in [args.scenario] if args.scenario else ['GFL', 'GFM']:
        with (output / f'{name}.log').open('w') as log:
            subprocess.run([str(args.exe.resolve()), str(here / f'{name}.solver.json')], cwd=output, stdout=log, stderr=subprocess.STDOUT, check=True)
        solver = json.loads((here / f'{name}.solver.json').read_text())
        paths = [args.exe.resolve(), here / f'{name}.solver.json']
        paths += [(here / solver[key]).resolve() for key in ['system_model_file', 'state_file']]
        record = {'command': [str(args.exe.resolve()), str(here / f'{name}.solver.json')],
                  'sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
        (output / f'{name}.run.json').write_text(json.dumps(record, indent=2) + '\n')
        print(f'{name}: {output / (name + ".csv")}', flush=True)


if __name__ == '__main__':
    main()
