"""Normalize and compare actual GridKit and ParaEMT outputs for both disturbances."""
import gzip
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
PLOTS = ROOT / 'plots'
GROUPS = {'phase_voltage': lambda c: '_v' in c and c.endswith(('a_pu','b_pu','c_pu')),
          'voltage_magnitude': lambda c: '_vm_pu' in c,
          'rotor_speed': lambda c: '_speed_pu' in c,
          'mechanical_power': lambda c: '_pm_pu' in c,
          'field_voltage': lambda c: '_efd_pu' in c,
          'stabilizer_output': lambda c: '_pss_vs_pu' in c}


def normalized_gridkit(directory):
    raw = directory / 'mon.csv'
    if not raw.exists(): raw = directory / 'mon.csv.gz'
    df = pd.read_csv(raw)
    df['t'] = df.t.round(10)
    # The event restart records t=1 twice. Differential states are fixed;
    # both event limits are retained in raw monitoring. Keep the left limit
    # in this regular-grid CSV; the trip voltage impulse is not a finite sample.
    duplicates = df[df.t.duplicated(keep=False)]
    if len(duplicates) != 2 or not np.allclose(duplicates.t, 1.0):
        raise ValueError('Unexpected event records or discontinuous monitored output')
    (directory/'restart_monitor_adjustments.json').write_text((duplicates.iloc[1]-duplicates.iloc[0]).to_json(indent=2)+'\n')
    df = df.drop_duplicates('t',keep='first').reset_index(drop=True)
    out = {'time_s': df.t.to_numpy()}
    peak = np.sqrt(2/3)*230e3
    for bus in range(1,10):
        for phase in 'abc': out[f'bus{bus}_v{phase}_pu'] = df[f'Bus_bus_{bus}_v{phase}']/peak
        out[f'bus{bus}_vm_pu'] = np.sqrt(sum(out[f'bus{bus}_v{phase}_pu']**2 for phase in 'abc')*2/3)
    for gen in range(1,4):
        for key,column in [('speed_pu',f'Machine_gen_{gen}_omega'),('efd_pu',f'Machine_gen_{gen}_efd'),('pm_pu',f'GASTPTI_gov_{gen}_pmech'),('pss_vs_pu',f'IEEEST_stabilizer_{gen}_vss')]:
            out[f'gen{gen}_{key}'] = df[column]
    frame = pd.DataFrame(out)
    assert frame.shape == (6001,49) and np.isfinite(frame.to_numpy()).all()
    np.testing.assert_allclose(frame.time_s, np.arange(6001)*.0005,rtol=0,atol=1e-12)
    frame.to_csv(directory/'gridkit.csv.gz',index=False,float_format='%.12e',compression={'method':'gzip','mtime':0})
    schema = json.loads((RESULTS/'governor_step/dt25us/columns.json').read_text())
    for key in schema: schema[key] = schema[key].replace('ParaEMT','GridKit').replace('GAST ','GASTPTI ')
    (directory/'columns.json').write_text(json.dumps(schema,indent=2)+'\n')
    return frame


def metrics(a,b):
    np.testing.assert_allclose(a.time_s,b.time_s,rtol=0,atol=1e-12)
    result = {}
    for interval,mask in {'all': np.ones(len(a),dtype=bool),'before_event': a.time_s < 1,
                          'event_through_20ms': (a.time_s >= 1)&(a.time_s <= 1.02),
                          'after_event': a.time_s > 1,
                          'after_20ms': a.time_s > 1.02}.items():
        result[interval] = {}
        for group,select in GROUPS.items():
            columns = [c for c in a if select(c)]
            error = a.loc[mask,columns].to_numpy()-b.loc[mask,columns].to_numpy()
            row,col = np.unravel_index(np.argmax(abs(error)),error.shape)
            result[interval][group] = {'max_abs_pu': float(abs(error[row,col])), 'rms_pu': float(np.sqrt(np.mean(error**2))),
              'channel': columns[col], 'time_s': float(a.loc[mask,'time_s'].iloc[row])}
    return result


def compress(path):
    if not path.exists(): return
    with path.open('rb') as source, path.with_suffix(path.suffix+'.gz').open('wb') as target:
        with gzip.GzipFile(fileobj=target,mode='wb',mtime=0) as compressed: shutil.copyfileobj(source,compressed)
    path.unlink()


def compare_event(event):
    gridroot=RESULTS/('gridkit' if event == 'governor_step' else 'gridkit_trip')
    grids={label: normalized_gridkit(gridroot/label) for label in ('tol1e-7','tol1e-8','tol1e-9')}
    refs={label: pd.read_csv(RESULTS/event/label/'reference.csv.gz') for label in ('dt50us','dt25us','dt12_5us')}
    grid,ref=grids['tol1e-9'],refs['dt12_5us']
    summary={'experiment':event,
      'comparison':'Same 0.5 ms samples; no time/angle shifting, fitted gains, or baseline subtraction. GridKit t=1 is the left limit; ParaEMT t=1 includes its event step.',
      'gridkit_to_paraemt':{key:metrics(grid,value) for key,value in refs.items()},
      'gridkit_refinement':{'1e-7_to_1e-8':metrics(grids['tol1e-7'],grids['tol1e-8']),'1e-8_to_1e-9':metrics(grids['tol1e-8'],grid)},
      'paraemt_refinement':{'50_to_25us':metrics(refs['dt50us'],refs['dt25us']),'25_to_12.5us':metrics(refs['dt25us'],ref)},
      'initialization':json.loads((ROOT/'gridkit/initialization.json').read_text()),
      'per_channel':{}}
    for column in grid.columns.drop('time_s'):
        summary['per_channel'][column]={}
        for interval,mask in {'all':grid.time_s>=0,'after_event':grid.time_s>1,'after_20ms':grid.time_s>1.02}.items():
            error=grid.loc[mask,column]-ref.loc[mask,column]
            index=error.abs().idxmax()
            summary['per_channel'][column][interval]={'max_abs_pu':float(abs(error.loc[index])),
              'rms_pu':float(np.sqrt(np.mean(error**2))),'time_s':float(grid.time_s.loc[index])}
    if event == 'trip':
        summary['limitations']=[
          'GridKit uses an ideal terminal opening with explicit current-interruption projection. No finite breaker arc, snubber, or voltage impulse amplitude is modeled.',
          'ParaEMT removes G1 Norton injection and conductance, freezes its electrical history, but continues its state kernel. G1 post-trip machine/controller traces are not equivalent physical models.',
          'All traces and the nonconverging ParaEMT event sample are retained. Cross-simulator discrepancies are observations, not pass/fail errors.']
        connected=[c for c in grid if not c.startswith('gen1_')]
        summary['connected_generators_and_network']={key:metrics(grid[connected],value[connected]) for key,value in refs.items()}
        summary['paraemt_bus1_event_magnitude_pu']={key:float(value.loc[value.time_s.eq(1),'bus1_vm_pu'].iloc[0]) for key,value in refs.items()}
    target=RESULTS/('gridkit_comparison.json' if event=='governor_step' else 'gridkit_trip_comparison.json')
    target.write_text(json.dumps(summary,indent=2)+'\n')
    for label in grids:
        directory=gridroot/label
        compress(directory/'mon.csv')
        if label=='tol1e-9': compress(directory/'state.csv')
        else: (directory/'state.csv').unlink(missing_ok=True)
    print(event, json.dumps(summary['gridkit_to_paraemt']['dt12_5us']['after_20ms'],indent=2))
    return grids,refs


def main():
    for event in ('governor_step','trip'): compare_event(event)
    from plot_comparison import main as plot
    plot()


if __name__ == '__main__': main()
