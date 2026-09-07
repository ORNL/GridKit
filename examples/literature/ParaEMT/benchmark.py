"""Serial, CPU-pinned runtime trials; separately report capture-run costs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import shutil
import subprocess
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]
OUT=ROOT/'results/runtime'
CONFIGS={'ParaEMT':('dt50us','dt25us','dt12_5us'),'GridKit':('tol1e-7','tol1e-8','tol1e-9')}
ENV={'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','NUMBA_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}


def capture(event,sim,setting):
    folder=event if sim=='ParaEMT' else 'gridkit' if event=='governor_step' else 'gridkit_trip'
    return json.loads((ROOT/'results'/folder/setting/'run.json').read_text())


def run_trials(args):
    cpu=min(os.sched_getaffinity(0))
    environment=os.environ.copy();environment.update(ENV)
    environment.setdefault('MPLCONFIGDIR','/tmp/paraemt-matplotlib')
    machine={'platform':platform.platform(),'processor':next(line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),
             'logical_cpus':os.cpu_count(),'affinity_cpu':cpu,'environment':ENV,'repetitions':args.repetitions,
             'GridKit_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
             'execution':'Sequential fresh child processes pinned to one logical CPU; no concurrent benchmark jobs. Initialization includes ParaEMT preprocess/JIT; timed loops begin afterward.',
             'limits':'This compares implementation cost at the stated settings, not at matched waveform error. Trip G1 models differ. Python imports/process startup are included only in process wall time.'}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'environment.json').write_text(json.dumps(machine,indent=2)+'\n')
    for event in ('governor_step','trip'):
      for trial in range(1,args.repetitions+1):
       for sim,settings in CONFIGS.items():
        for setting in settings:
            destination=OUT/event/sim/setting/f'trial{trial}';destination.mkdir(parents=True,exist_ok=True)
            previous=ROOT/'results/runtime_cold'/event/sim/setting/f'trial{trial}'
            if sim=='GridKit' and (previous/'run.json').exists():
                # GridKit is compiled ahead of time; these already measured
                # trials are unaffected by correcting ParaEMT JIT warm-up.
                for name in ('run.json','run.log'):shutil.copyfile(previous/name,destination/name)
                run=json.loads((destination/'run.json').read_text())
                run['reused_from']=str(previous.relative_to(ROOT))
                (destination/'run.json').write_text(json.dumps(run,indent=2)+'\n')
                continue
            if sim=='GridKit':
                solver=ROOT/'gridkit'/f'{"trip-" if event=="trip" else ""}{setting}.solver.json'
                command=[str(args.executable.resolve()),str(solver)]+(['trip'] if event=='trip' else [])+['--benchmark']
                binary_hash=hashlib.sha256(args.executable.read_bytes()).hexdigest()
            else:
                step={'dt50us':50,'dt25us':25,'dt12_5us':12.5}[setting]
                command=[sys.executable,str(ROOT/'run_reference.py'),'--event',event.replace('_','-'),'--dt-us',str(step),
                         '--monitor-us','50','--benchmark','--output',str(destination)]
                binary_hash=hashlib.sha256((ROOT/'run_reference.py').read_bytes()).hexdigest()
            command=['taskset','-c',str(cpu)]+command
            print('Timing',event,sim,setting,trial,flush=True)
            begin=time.perf_counter()
            with (destination/'run.log').open('w') as log:
                subprocess.run(command,cwd=destination,env=environment,stdout=log,stderr=subprocess.STDOUT,check=True)
            elapsed=time.perf_counter()-begin
            run=json.loads((destination/'run.json').read_text())
            assert not run['result_capture']
            if sim=='ParaEMT':assert run['jit_specializations_stable']
            if sim=='GridKit': run['max_generator_kcl_mismatch_A']=None
            baseline=capture(event,sim,setting)
            np.testing.assert_allclose(run['final_state'],baseline['final_state'],rtol=1e-12,atol=1e-12,
                                       err_msg='Disabling result capture changed the trajectory')
            run.update(command=command,process_wall_s=elapsed,artifact_sha256=binary_hash,
                       final_state_matches_capture=True,GridKit_revision=machine['GridKit_revision'])
            (destination/'run.json').write_text(json.dumps(run,indent=2)+'\n')


def summarize():
    rows=[]
    for event in ('governor_step','trip'):
      for sim,settings in CONFIGS.items():
       for setting in settings:
        trials=[json.loads(p.read_text()) for p in sorted((OUT/event/sim/setting).glob('trial*/run.json'))]
        assert len(trials)>=3 and all(t['final_state_matches_capture'] for t in trials)
        steps=[t['steps'] if sim=='ParaEMT' else t['solver']['steps'] for t in trials]
        row={'event':event,'simulator':sim,'setting':setting,'trials':len(trials),
             'accepted_or_fixed_steps':statistics.median(steps),'duration_per_step_us':3e6/statistics.median(steps),
             'capture_loop_wall_s':capture(event,sim,setting)['loop_wall_s']}
        for key in ('loop_wall_s','loop_cpu_s','initialization_wall_s','process_wall_s'):
            values=[t[key] for t in trials]
            row[key]={'median':statistics.median(values),'min':min(values),'max':max(values)}
        row['jit_warmup_wall_s']=statistics.median(t.get('jit_warmup_wall_s',0) for t in trials)
        if sim=='GridKit':
            row['residual_evaluations']=statistics.median(t['solver']['residual_evaluations'] for t in trials)
            row['error_test_failures']=statistics.median(t['solver']['error_test_failures'] for t in trials)
        rows.append(row)
    (OUT/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    env=json.loads((OUT/'environment.json').read_text())
    lines=['# ParaEMT / GridKit runtimes','',
      'Three-second trajectories; medians of three processes per setting. Benchmark processes run sequentially on one pinned logical CPU with BLAS/OMP/Numba thread counts set to one. Result capture is disabled. Every benchmark final state agrees with its corresponding captured trajectory to `rtol=atol=1e-12`. ParaEMT stepping kernels are warmed on a disposable state before timing; their JIT signatures must remain unchanged throughout the timed loop. GridKit uses the previously measured trials of its unchanged compiled executable.','',
      f'Host: **{env["processor"]}**; pinned logical CPU {env["affinity_cpu"]}. GridKit revision `{env["GridKit_revision"]}`.','',
      '| Event | Simulator / setting | Loop wall, s | Loop CPU, s | Initialization, s | Process wall, s | Steps | Duration / steps, µs |',
      '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in rows:
        lines.append(f'| {r["event"]} | {r["simulator"]} {r["setting"]} | {r["loop_wall_s"]["median"]:.4f} | {r["loop_cpu_s"]["median"]:.4f} | {r["initialization_wall_s"]["median"]:.4f} | {r["process_wall_s"]["median"]:.4f} | {r["accepted_or_fixed_steps"]:,.0f} | {r["duration_per_step_us"]:.3f} |')
    lines += ['', 'ParaEMT initialization here includes explicit stepping-kernel warm-up in addition to network initialization. Its ordinary startup message does not mean all kernels are compiled: the first-use costs inside the original loops are retained in [RUNTIMES_COLD.md](RUNTIMES_COLD.md). Process wall time includes imports, loading, initialization, warm-up, the timed trajectory and metadata output. GridKit loop timing includes the event consistency solve. CPU time is measured inside each process.', '',
      '**Adaptive steps and output samples are different.** GridKit uses variable-step, variable-order IDA/BDF; the monitor interval is 50 µs. ParaEMT uses fixed 50, 25 or 12.5 µs integration steps. Both captured outputs use 50 µs spacing, giving 60,001 samples. The duration/steps column is an effective average, not a fixed GridKit step or a recorded step-size distribution. IDA may internally pass an output time and interpolate back.', '',
      'This is a cost comparison at the listed settings, not an equal-error benchmark. Use the waveform/refinement metrics alongside runtime. The governor event has closely mapped continuous models; the trip has different G1 post-event semantics, so trip speed ratios do not establish equivalent-model performance.', '',
      '## Result-capture runs', '',
      'These single runs generate the plotted data. GridKit streams monitor and full-state CSVs during its loop; ParaEMT copies samples into memory and exports files after its loop. Their capture-loop costs therefore include different output work and are not the primary timing comparison. They were not CPU-pinned. Compression and plotting are excluded from the benchmark loops.', '',
      '| Event | Simulator / setting | Capture loop wall, s | No-capture median loop wall, s |',
      '| --- | --- | ---: | ---: |']
    for r in rows:lines.append(f'| {r["event"]} | {r["simulator"]} {r["setting"]} | {r["capture_loop_wall_s"]:.4f} | {r["loop_wall_s"]["median"]:.4f} |')
    lines += ['', '[Raw trials, environment and ranges](results/runtime/) · [Runtime plot](plots/runtime_comparison.png) · [Governor adaptive work](plots/governor_step_adaptive_work.png) · [Trip adaptive work](plots/trip_adaptive_work.png)', '',
      'Reproduce with the configured Python environment: `python benchmark.py`. Use `python benchmark.py --summarize-only` to rebuild tables and plots from saved trials. Each raw trial records its command, executable/wrapper hash, settings, timing, final state and solver work.','']
    (ROOT/'RUNTIMES.md').write_text('\n'.join(lines))
    fig,axes=plt.subplots(2,1,figsize=(12,9),layout='constrained')
    for ax,event in zip(axes,('governor_step','trip')):
        selected=[r for r in rows if r['event']==event]
        values=[r['loop_wall_s']['median'] for r in selected]
        lower=[v-r['loop_wall_s']['min'] for v,r in zip(values,selected)]
        upper=[r['loop_wall_s']['max']-v for v,r in zip(values,selected)]
        ax.bar(range(6),values,yerr=[lower,upper],capsize=4,color=['#1768ac']*3+['#d45c21']*3)
        ax.set_xticks(range(6),[r['simulator']+'\n'+r['setting'] for r in selected]);ax.set_ylabel('Simulation-loop wall time (s)')
        ax.set_title(event.replace('_',' ')+' · median and min/max of three trials');ax.grid(axis='y',alpha=.2)
    fig.suptitle('Warmed simulation loops · pinned CPU · result capture disabled · 3 s simulated')
    fig.savefig(ROOT/'plots/runtime_comparison.png',dpi=160);plt.close(fig)
    for event in ('governor_step','trip'):
        fig,axes=plt.subplots(2,1,figsize=(12,9),sharex=True,layout='constrained')
        folder=ROOT/'results'/('gridkit' if event=='governor_step' else 'gridkit_trip')
        for setting in CONFIGS['GridKit']:
            data=pd.read_csv(folder/setting/'solver_work.csv')
            assert len(data)==301 and data.accepted_steps.is_monotonic_increasing
            delta=data.accepted_steps.diff();dt=data.time_s.diff()
            axes[0].plot(data.time_s,data.accepted_steps,label='GridKit '+setting)
            axes[1].plot(data.time_s,1e6*dt/delta,label='GridKit '+setting)
        for step,color in zip((50,25,12.5),('#d62728','#9467bd','#8c564b')):
            axes[0].plot([0,3],[0,3e6/step],ls='--',color=color,lw=1,label=f'ParaEMT {step:g} µs')
            axes[1].axhline(step,ls='--',color=color,lw=1,label=f'ParaEMT {step:g} µs')
        for ax in axes:ax.grid(alpha=.2);ax.axvline(1,color='.4',ls=':');ax.legend(fontsize=8,ncol=2)
        axes[0].set_ylabel('Cumulative accepted / fixed steps')
        axes[1].set_ylabel('10 ms / accepted steps in interval (µs)');axes[1].set_xlabel('Time (s)')
        fig.suptitle(event.replace('_',' ')+' · solver work through the event\nGridKit output interval: 50 µs; lower plot is an interval average, not each internal step')
        fig.savefig(ROOT/f'plots/{event}_adaptive_work.png',dpi=160);plt.close(fig)
    print(json.dumps(rows,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,default=REPO/'build/paraemt-9bus/paraemt_9bus')
    parser.add_argument('--repetitions',type=int,default=3)
    parser.add_argument('--summarize-only',action='store_true')
    args=parser.parse_args()
    if args.repetitions<3:parser.error('Use at least three trials')
    if not args.summarize_only:run_trials(args)
    summarize()


if __name__=='__main__':main()
