"""Re-run all gallery experiments against the configured GridKit build."""
from pathlib import Path
import gzip, hashlib, io, json, shlex, shutil, subprocess, tempfile, zipfile
import numpy as np
import csv

HERE=Path(__file__).resolve().parent
REPO=next(p for p in HERE.parents if (p/'GridKit/Solver').is_dir())
BUILD=REPO/'build'
def read_data(name):
    with zipfile.ZipFile(HERE/'data.zip') as archive:
        return io.BytesIO(archive.read(name))

def generate():
    link_dir=BUILD/'tests/IntegrationTests/PhasorDynamics'
    flags=(link_dir/'CMakeFiles/test_pd_signal_delay.dir/flags.make').read_text()
    includes=shlex.split(next(x.split('=',1)[1] for x in flags.splitlines() if x.startswith('CXX_INCLUDES')))
    link=shlex.split((link_dir/'CMakeFiles/test_pd_signal_delay.dir/link.txt').read_text())
    libraries=link[link.index('-o')+2:]
    with tempfile.TemporaryDirectory(prefix='gridkit-delay-gallery-') as tmp:
        tmp=Path(tmp)
        subprocess.run([link[0],'-std=c++20','-O2',*includes,str(HERE/'gallery.cpp'),'-o',str(tmp/'gallery'),*libraries],cwd=link_dir,check=True)
        with zipfile.ZipFile(HERE/'data.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
            for m in json.loads((HERE/'runs.json').read_text()):
                keys=['name','kind','tau','tf','rtol','compress','fixed_h','width','dt_monitor','history_rtol']
                with (tmp/(m['name']+'.log')).open('w') as log:
                    subprocess.run([str(tmp/'gallery'),*[str(m[k]) for k in keys]],cwd=tmp,stdout=log,stderr=log,check=True)
                for kind in ['samples','steps','attempts','midpoints']:
                    name=m['name']+'_'+kind+'.csv'
                    archive.write(tmp/name,name)
                archive.write(tmp/(m['name']+'.log'),m['name']+'.log')
                print(m['name'],flush=True)
        for name in ['app','grid']:
            target=tmp/name;target.mkdir()
            for file in ['case.json','solver.json']:shutil.copyfile(HERE/name/file,target/file)
            with (HERE/name/'run.log').open('w') as log:
                subprocess.run([str(BUILD/'application/PhasorDynamics/DynamicSimulation'),'solver.json'],cwd=target,stdout=log,stderr=log,check=True)
            with gzip.open(HERE/name/'output.csv.gz','wb') as f:f.write((target/'mon.csv').read_bytes())
    provenance=json.loads((HERE/'provenance.json').read_text())
    provenance['head']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
    provenance['diff_sha256']=hashlib.sha256(subprocess.check_output(['git','diff','HEAD'],cwd=REPO)).hexdigest()
    (HERE/'provenance.json').write_text(json.dumps(provenance,indent=2))

def summarize():
    rows=[]
    for m in json.loads((HERE/'runs.json').read_text()):
     if m['returncode']:continue
     a=np.genfromtxt(read_data(m['name']+'_samples.csv'),delimiter=',',names=True)
     s=np.genfromtxt(read_data(m['name']+'_steps.csv'),delimiter=',',names=True)
     attempts=np.genfromtxt(read_data(m['name']+'_attempts.csv'),delimiter=',',names=True)
     mid=np.genfromtxt(read_data(m['name']+'_midpoints.csv'),delimiter=',',names=True)
     e=np.r_[np.abs(a['z']-a['exact_z']),np.abs(mid['z']-mid['exact_z'])]
     # Exact at step endpoints, separately from IDA dense monitoring.
     q=s['t']-m['tau']
     if m['kind']=='feedback':ez=np.exp(-.4*q)
     elif m['kind']=='pulse':ez=1+.6*np.exp(-((q-2)/m['width'])**2)
     elif m['kind']=='chirp':ez=np.cos(q+.4*q*q)
     else:ez=np.cos(q)
     row={**m,'accepted':len(s),'rejected':int(sum(attempts['rejected'])),'max_abs_z_error':float(max(e)),
     'max_endpoint_z_error':float(max(abs(s['z']-ez))),
     'max_abs_u_error':float(max(abs(a['u']-a['exact_u']))),
     'min_h_tau':float(min(s['h_over_tau'])),'median_h_tau':float(np.median(s['h_over_tau'])),'max_h_tau':float(max(s['h_over_tau'])),
     'peak_segments':int(max(s['segments'])),'last_segments':int(s['segments'][-1]),
     'shift_signature_error':float(max(abs((a['u']-a['z'])-(a['exact_u']-a['exact_z'])))/m['tau'])}
     rows.append(row)
    with (HERE/'summary.csv').open('w') as f:
     w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (HERE/'summary.json').write_text(json.dumps(rows,indent=2))

if __name__=='__main__':
    generate()
    summarize()
    subprocess.run(['python3',str(HERE/'plot_gallery.py')],check=True)
