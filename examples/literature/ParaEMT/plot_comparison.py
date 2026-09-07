"""One response per figure: ParaEMT, GridKit, then error, with identical axes."""
import html
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
PLOTS=ROOT/'plots'
COLORS=('#1768ac','#d45c21','#298254')


def responses():
    for bus in range(1,10):
        yield f'bus{bus}_abc',f'Bus {bus} · three-phase voltage',[f'bus{bus}_v{p}_pu' for p in 'abc'],'Phase voltage (pu peak base)',0
        yield f'bus{bus}_vm',f'Bus {bus} · voltage magnitude',[f'bus{bus}_vm_pu'],'Three-phase magnitude (pu)',0
    for gen in range(1,4):
        for suffix,label,units,offset in [('speed_pu','rotor speed deviation','Speed − synchronous speed (pu)',1),
                                         ('pm_pu','mechanical power','Mechanical power (pu machine base)',0),
                                         ('efd_pu','field voltage','Field voltage (pu exciter base)',0),
                                         ('pss_vs_pu','stabilizer output','Stabilizer output (pu)',0)]:
            yield f'gen{gen}_{suffix}',f'Generator {gen} · {label}',[f'gen{gen}_{suffix}'],units,offset


def select(data,window):
    return data[(data.time_s>=window[0]-1e-10)&(data.time_s<=window[1]+1e-10)]


def generate(event):
    folder=PLOTS/event;folder.mkdir(parents=True,exist_ok=True)
    gridfolder=ROOT/'results'/('gridkit' if event=='governor_step' else 'gridkit_trip')/'tol1e-9'
    reffolder=ROOT/'results'/event/'dt12_5us'
    grid=pd.read_csv(gridfolder/'gridkit.csv.gz')
    ref=pd.read_csv(reffolder/'reference.csv.gz')
    native=pd.read_csv(reffolder/'event_waveforms.csv.gz')
    np.testing.assert_allclose(grid.time_s,ref.time_s,rtol=0,atol=1e-12)
    title='Governor reference −0.02 pu at 1 s' if event=='governor_step' else 'Generator 1 trip at 1 s'
    pages=[]
    with PdfPages(folder/'comparison.pdf',metadata={'Title':f'ParaEMT / GridKit 9-bus · {title}','Author':'GridKit case reproduction'}) as pdf:
        for name,label,columns,units,offset in responses():
            windows=[('full',(0.,3.)),('event',(.99,1.025) if name.startswith('bus') else (.95,1.35))]
            if name.endswith('_abc'):windows=[('event',(.99,1.025)),('late',(2.99,3.))]
            for view,window in windows:
                g,r=select(grid,window),select(ref,window)
                # Native-step data exposes the ParaEMT event spike. Error uses
                # the exact common 50 µs sample times, never interpolation.
                dense=view=='event' and name.startswith('bus')
                # Keep the short native-step interval and regular samples on
                # either side; no interpolation or missing ends of the window.
                top=select(pd.concat([ref,native]).drop_duplicates('time_s',keep='last').sort_values('time_s'),window) if dense else r
                error=g[columns].to_numpy()-r[columns].to_numpy()
                values=[top[columns].to_numpy()-offset,g[columns].to_numpy()-offset,error]
                assert all(np.isfinite(v).all() for v in values)
                low=min(float(v.min()) for v in values);high=max(float(v.max()) for v in values)
                margin=max((high-low)*.06,1e-10)
                limits=(low-margin,high+margin)
                fig,axes=plt.subplots(3,1,figsize=(12,10),sharex=True,sharey=True)
                for ax,data,time,panel in zip(axes,values,(top.time_s,g.time_s,g.time_s),
                    ('ParaEMT · 12.5 µs integration step'+(' · native samples near trip' if dense else ''),
                     'GridKit · IDA tolerance 1e−9','Error · GridKit − ParaEMT')):
                    for i,column in enumerate(columns):
                        ax.plot(time,data[:,i],color=COLORS[i],lw=1.25,label=f'Phase {"abc"[i]}' if len(columns)==3 else label)
                    ax.set_title(panel,loc='left',fontsize=11)
                    ax.set_ylabel(units,fontsize=10)
                    ax.set_xlim(*window);ax.set_ylim(*limits)
                    ax.grid(alpha=.22);ax.ticklabel_format(axis='y',style='plain',useOffset=False)
                    if window[0]<=1<=window[1]:ax.axvline(1,color='.45',ls=':',lw=.8)
                axes[2].axhline(0,color='.4',lw=.7)
                axes[2].set_xlabel('Time (s)')
                if len(columns)==3:
                    axes[0].legend(handles=[Line2D([0],[0],color=c,label=f'Phase {p}') for p,c in zip('abc',COLORS)],loc='upper right',ncol=3,fontsize=9)
                fig.suptitle(f'{title}\n{label} · {view.replace("event","event zoom").replace("late","late waveform")}',fontsize=14)
                note='Identical x/y scales in all three panels. Error at common 50 µs samples; no time shift or error magnification.'
                if dense:note+='\nParaEMT shows every integration step from 0.995 to 1.015 s, with 50 µs samples outside that interval.'
                if event=='trip':
                    note+='\nGridKit t=1 is the left limit; ParaEMT includes its trip step. GridKit does not assign a finite voltage-impulse amplitude.'
                    if name.startswith('gen1_'):note+='\nG1 post-trip models differ: GridKit isolates the machine; ParaEMT freezes electrical history but continues the state kernel.'
                fig.text(.5,.012,note,ha='center',va='bottom',fontsize=8)
                fig.tight_layout(rect=(0,.08,1,.975))
                target=folder/f'{name}_{view}.png'
                fig.savefig(target,dpi=140);pdf.savefig(fig);plt.close(fig)
                page={'image':str(target.relative_to(PLOTS)),'title':label,'response':name,'view':view,'channels':columns,
                      'x_limits_s':list(window),'y_limits':list(limits),'panels':['ParaEMT','GridKit','GridKit - ParaEMT'],
                      'shared_x_and_y':True,'response_offset_pu':offset,'paraemt_native_steps':dense}
                pages.append(page)
                print(f'Plotted {event}/{target.name}',flush=True)
    covered=set().union(*(set(p['channels']) for p in pages))
    assert covered==set(grid.columns)-{'time_s'}
    return {'event':event,'title':title,'pdf':f'{event}/comparison.pdf','pages':pages,'covered_channels':sorted(covered)}


def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'path.simplify':True,'path.simplify_threshold':.05})
    # Remove only the superseded figures named by this generator's old manifest.
    previous=PLOTS/'coverage.json'
    if previous.exists():
        for report in json.loads(previous.read_text()):
            for page in report['pages']:(PLOTS/page['image']).unlink(missing_ok=True)
    reports=[generate(event) for event in ('trip','governor_step')]
    previous.write_text(json.dumps(reports,indent=2)+'\n')
    document=['<!doctype html><html lang="en"><meta charset="utf-8"><title>GridKit / ParaEMT responses</title>',
      '<style>body{font:17px system-ui;max-width:1400px;margin:2rem auto;padding:0 1rem;color:#172a3a;background:#f6f8fa}a{color:#075b9d}nav{display:flex;gap:2rem}.gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:1rem}article{background:white;padding:1rem;border:1px solid #dce1e5;border-radius:6px}img{width:100%}h2{margin-top:3rem}table{border-collapse:collapse}td,th{padding:.5rem 1rem;border-bottom:1px solid #dce1e5}</style>',
      '<h1>ParaEMT / GridKit: each response separately</h1><nav><a href="#trip">Trip</a><a href="#governor_step">Governor step</a><a href="../RUNTIMES.md">Runtimes</a></nav>',
      '<p>Every figure: ParaEMT on top, GridKit in the middle, GridKit − ParaEMT on the bottom. All three panels use identical x/y scales. Three-phase signals stay together. Speed is shown as deviation from synchronous speed in pu. Errors use common output times without interpolation or phase alignment.</p>',
      '<p>Trip caveat: GridKit models an ideal terminal opening. ParaEMT freezes generator 1 electrical history but continues its state kernel. Its post-trip generator 1 response and trip spike cannot validate the disconnected GridKit machine.</p>']
    for report in reports:
        document += [f'<section id="{report["event"]}"><h2>{html.escape(report["title"])}</h2>',f'<p><a href="{report["pdf"]}">Complete PDF · {len(report["pages"])} figures</a></p>',
                     '<table><tr><th>Response</th><th>Full run</th><th>Event zoom</th><th>Late waveform</th></tr>']
        for name,label,_,_,_ in responses():
            row=[p for p in report['pages'] if p['response']==name]
            links={p['view']:f'<a href="{p["image"]}">{p["view"]}</a>' for p in row}
            document.append(f'<tr><td>{html.escape(label)}</td>'+''.join(f'<td>{links.get(v,"")}</td>' for v in ('full','event','late'))+'</tr>')
        document.append('</table><div class="gallery">')
        for page in report['pages']:
            if page['view']!='event':continue
            label=html.escape(page['title'])
            document.append(f'<article><h3>{label} · event zoom</h3><a href="{page["image"]}"><img loading="lazy" src="{page["image"]}" alt="{label}"></a></article>')
        document.append('</div></section>')
    document.append('</html>');(PLOTS/'index.html').write_text('\n'.join(document)+'\n')
    for target,source in {'gridkit_governor_step.png':'governor_step/gen1_speed_pu_full.png',
       'gridkit_voltages.png':'governor_step/bus1_abc_event.png','gridkit_stabilizers_and_refinement.png':'governor_step/gen1_pss_vs_pu_full.png',
       'governor_step_response.png':'governor_step/gen1_speed_pu_full.png','trip_response.png':'trip/gen2_speed_pu_full.png',
       'trip_waveforms.png':'trip/bus4_abc_event.png'}.items():shutil.copyfile(PLOTS/source,PLOTS/target)
    print('Saved reports:',[(r['event'],len(r['pages']),len(r['covered_channels'])) for r in reports],flush=True)


if __name__=='__main__':main()
