"""Build and run the full GridKit ParaEMT 9-bus comparison (local outputs only)."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gridkit-build',type=Path,default=REPO/'build/emt-rational')
    parser.add_argument('--driver-build',type=Path,default=REPO/'build/paraemt-9bus')
    args = parser.parse_args()
    gridkit_build,driver_build=args.gridkit_build.resolve(),args.driver_build.resolve()
    # Use the same compiler/ABI as the configured GridKit build.
    cache=(gridkit_build/'CMakeCache.txt').read_text().splitlines()
    compiler=next(line.split('=',1)[1] for line in cache if line.startswith('CMAKE_CXX_COMPILER:'))
    subprocess.run([sys.executable,str(ROOT/'make_gridkit_case.py')],check=True)
    subprocess.run(['cmake','-S',str(ROOT),'-B',str(driver_build),f'-DGridKit_DIR={gridkit_build}',f'-DCMAKE_CXX_COMPILER={compiler}','-DCMAKE_BUILD_TYPE=Release'],check=True)
    subprocess.run(['cmake','--build',str(driver_build),'-j','10'],check=True)
    executable=driver_build/'paraemt_9bus'
    for event in ('governor_step','trip'):
      for label in ('tol1e-7','tol1e-8','tol1e-9'):
        output=ROOT/'results'/('gridkit' if event=='governor_step' else 'gridkit_trip')/label;output.mkdir(parents=True,exist_ok=True)
        command=[str(executable),str(ROOT/'gridkit'/f'{"trip-" if event=="trip" else ""}{label}.solver.json')]
        if event=='trip': command.append('trip')
        print('Running',event,label,flush=True)
        with (output/'run.log').open('w') as log: subprocess.run(command,cwd=output,stdout=log,stderr=subprocess.STDOUT,check=True)
        run=json.loads((output/'run.json').read_text())
        run['GridKit_revision']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
        run['command']=command
        run['executable_sha256']=hashlib.sha256(executable.read_bytes()).hexdigest()
        run['case_driver']='gridkit_9bus.cpp and GeneratorTerminalConstraint.hpp (local literature files)'
        (output/'run.json').write_text(json.dumps(run,indent=2)+'\n')
    subprocess.run([sys.executable,str(ROOT/'compare_gridkit.py')],check=True)


if __name__ == '__main__': main()
