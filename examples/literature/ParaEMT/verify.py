"""Verify artifact identity, complete channels, event invariants and current balance.

This checks integrity and internal consistency, not an externally specified
cross-simulator waveform acceptance tolerance.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parent


def check_event(directory,index,models):
    run=json.loads((directory/'run.json').read_text())
    limits=json.loads((directory/'event_state_limits.json').read_text())
    before,projected,after,tag=(limits[k] for k in ('before','after_projection','after_consistency','differential'))
    assert len(before)==len(projected)==len(after)==len(tag)==run['dae_variables']
    assert all(math.isfinite(x) for state in (before,projected,after) for x in state)
    projection=run['event']['explicit_state_projection']
    changed={p['index'] for p in projection}
    for i,differential in enumerate(tag):
        if i not in changed: assert before[i]==projected[i],(directory,i)
        if differential: assert projected[i]==after[i],(directory,i)
    assert run['event_differential_state_change']==0
    if run['event']['type']!='ideal-generator-terminal-opening': assert not changed;return
    expected={index[('gen_1',i)] for i in (2,3,4)}|{index[('xfmr_4_1',i)] for i in range(3)}
    assert changed==expected
    p=models['gen_1']['params']
    def y(i):return after[index[('gen_1',i)]]
    for axis,state,r1,r2,leak1,leak2 in [('d',2,5,6,'Llfd','Ll1d'),('q',3,7,8,'Ll1q','Ll2q')]:
        lm=p['Lm'+axis];l1,l2=p[leak1],p[leak2]
        # Rotor-flux equations at zero stator current, solved independently.
        determinant=(lm+l1)*(lm+l2)-lm*lm
        first=((lm+l2)*y(r1)-lm*y(r2))/determinant
        second=((lm+l1)*y(r2)-lm*y(r1))/determinant
        assert abs(y(state)-lm*(first+second))<1e-12
    assert max(abs(y(i)) for i in (9,10,11,21,22,23))<1e-10
    assert all(after[index[('xfmr_4_1',i)]]==0 for i in range(3))


def check_state(fine,models,trip):
    layout=json.loads((fine/'state.csv.json').read_text())['variables']
    index={(entry['component'],entry['local_index']):entry['index'] for entry in layout}
    pairs=[]
    for generator,transformer in zip(range(1,4),('xfmr_4_1','xfmr_7_2','xfmr_9_3')):
        p=models[f'gen_{generator}']['params'];scale=math.sqrt(2/3)*p['S']/p['V']
        for phase in range(3):pairs.append((generator,index[(f'gen_{generator}',21+phase)],index[(transformer,phase)],index[(transformer,6+phase)],scale))
    maximum=0
    with gzip.open(fine/'state.csv.gz','rt',newline='') as stream:
        reader=csv.reader(stream);header=next(reader)
        assert len(header)==1+2*len(layout)
        for count,row in enumerate(reader,1):
            assert len(row)==len(header)
            values=[float(v) for v in row];assert all(map(math.isfinite,values))
            assert abs(values[0]-(count-1)*.0005)<1e-11
            for generator,gen,transformer,shunt,scale in pairs:
                machine_current=scale*values[1+gen]
                opened=trip and generator==1 and values[0]>1
                maximum=max(maximum,abs((0 if opened else machine_current)+values[1+transformer]+values[1+shunt]))
                if opened:maximum=max(maximum,abs(machine_current))
        assert count==6001
    reported=json.loads((fine/'run.json').read_text())['max_generator_kcl_mismatch_A']
    assert abs(maximum-reported)<1e-10,(maximum,reported)
    for directory in fine.parent.glob('tol*'):check_event(directory,index,models)
    return maximum


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--regenerated',action='store_true',help='skip frozen local hashes after rerunning; retain source and numerical checks')
    args=parser.parse_args()
    manifest=json.loads((ROOT/'sources.json').read_text())
    for entry in manifest['artifacts']:
        path=ROOT/entry['path'];assert path.is_file(),path
        if not args.regenerated or entry.get('origin')!='local':
            raw=path.read_bytes()
            assert len(raw)==entry['bytes'],path
            assert hashlib.sha256(raw).hexdigest()==entry['sha256'],path
            if 'git_blob_sha1' in entry:assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==entry['git_blob_sha1'],path
    refs=sorted((ROOT/'results').glob('*/*/reference.csv.gz'))
    runs=refs+sorted((ROOT/'results').glob('gridkit*/*/gridkit.csv.gz'))
    assert len(refs)==6 and len(runs)==12,(len(refs),len(runs))
    expected=None
    for path in runs:
        schema=json.loads((path.parent/'columns.json').read_text())
        regular={}
        with gzip.open(path,'rt',newline='') as stream:
            reader=csv.DictReader(stream)
            assert len(reader.fieldnames)==49 and set(reader.fieldnames)==set(schema),path
            expected=set(schema)-{'time_s'}
            for count,row in enumerate(reader,1):
                assert all(math.isfinite(float(value)) for value in row.values()),(path,count)
                assert abs(float(row['time_s'])-(count-1)*.0005)<1e-11,(path,count)
                regular[round(float(row['time_s']),10)]=row
            assert count==6001,(path,count)
        metadata=json.loads((path.parent/'run.json').read_text())
        if metadata['simulator'].startswith('GridKit'):
            assert max(metadata['jacobian_check_max_scaled_difference'].values())<1e-4
        else:
            assert metadata['event']['applied_time_s']==1
            native=path.parent/'event_waveforms.csv.gz';dt=metadata['dt_s']
            with gzip.open(native,'rt',newline='') as stream:
                reader=csv.DictReader(stream);assert set(reader.fieldnames)==set(schema)
                for count,row in enumerate(reader,1):
                    assert all(math.isfinite(float(v)) for v in row.values())
                    t=round(float(row['time_s']),10)
                    assert abs(t-(.995+(count-1)*dt))<1e-11
                    if t in regular:
                        for col in schema:assert row[col]==regular[t][col],(native,t,col)
                assert count==round(.02/dt)+1
    case=json.loads((ROOT/'gridkit/9bus.case.json').read_text())
    models={entry['id']:entry for entry in case['devices']}
    drift={event:check_state(ROOT/'results'/folder/'tol1e-9',models,event=='trip') for event,folder in [('governor_step','gridkit'),('trip','gridkit_trip')]}
    reduction=json.loads((ROOT/'gridkit/index_reduction.json').read_text())
    for gen in range(1,4):
        p=models[f'gen_{gen}']['params'];d=reduction[str(gen)]
        for axis,leak1,leak2 in [('d','Llfd','Ll1d'),('q','Ll1q','Ll2q')]:
            lm=p['Lm'+axis]
            matrix=[[-lm-p['Ll'],lm,lm],[-lm,lm+p[leak1],lm],[-lm,lm,lm+p[leak2]]]
            for col in range(3):assert abs(sum(d[axis][row]*matrix[row][col] for row in range(3))-(1 if col==0 else 0))<1e-12
    coverage=json.loads((ROOT/'plots/coverage.json').read_text())
    for report in coverage:
        assert set(report['covered_channels'])==expected
        covered=set()
        for page in report['pages']:
            image=ROOT/'plots'/page['image'];assert image.read_bytes().startswith(b'\x89PNG')
            covered.update(page['channels'])
        assert covered==expected
        pdf=(ROOT/'plots'/report['pdf']).read_bytes()
        assert pdf.startswith(b'%PDF') and len(re.findall(rb'/Type /Page\b',pdf))==len(report['pages'])
    for link in re.findall(r'(?:href|src)="([^"]+)"',(ROOT/'plots/index.html').read_text()):
        if not link.startswith('#'):assert (ROOT/'plots'/link).is_file(),link
    print(f"Verified {len(manifest['artifacts'])} artifacts, 12 normalized runs, 6 native-step event exports, both event projections, full-state current balance, winding inverses, and all 48 channels in both PDF/gallery reports. KCL maxima (A): {drift}")


if __name__=='__main__': main()
