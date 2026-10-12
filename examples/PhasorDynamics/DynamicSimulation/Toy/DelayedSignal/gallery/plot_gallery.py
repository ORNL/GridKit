"""Render measured GridKit results. Run with: python3 plot_gallery.py."""
from pathlib import Path
import os, json, base64, html, io, zipfile
os.environ.setdefault('MPLCONFIGDIR', '/tmp/gridkit-delay-matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
ROWS=json.loads((ROOT/'summary.json').read_text());BY={r['name']:r for r in ROWS}
COLORS=['#166b9c','#d95b2b','#27836b','#8758a4','#b48113','#364152']
plt.rcParams.update({'font.size':10,'axes.titlesize':12,'axes.labelsize':10,'figure.titlesize':17,'axes.prop_cycle':plt.cycler(color=COLORS),'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2,'legend.frameon':False,'savefig.facecolor':'white','axes.formatter.useoffset':False})
PAGES=[]
def data(name,kind='samples'):
 with zipfile.ZipFile(ROOT/'data.zip') as archive:
  return np.genfromtxt(io.BytesIO(archive.read(f'{name}_{kind}.csv')),delimiter=',',names=True)
def app(name):return np.genfromtxt(ROOT/name/'output.csv.gz',delimiter=',',names=True)
def layout(title,nrows=2,ncols=1,height=7):
 fig,ax=plt.subplots(nrows,ncols,figsize=(11.5,height),squeeze=False)
 fig.suptitle(title,x=.07,ha='left',fontweight='bold');return fig,ax
pdf=PdfPages(ROOT/'gallery.pdf')
def save(fig,name,title,caption):
 fig.tight_layout(rect=[.015,.065,.99,.94]);fig.text(.07,.018,caption.split('\n')[0],fontsize=9,color='#444444')
 fig.savefig(OUT/(name+'.png'),dpi=155);fig.savefig(OUT/(name+'.svg'));pdf.savefig(fig)
 PAGES.append((name,title,caption));plt.close(fig)

# 1: Actual application output from a grid component, no manufactured generator model.
a=app('grid');t=a['t'];u=a['Genrou_DV1_omega'];taus=[.001,.1,.5]
fig,axes=layout('01  A generator signal delayed inside an actual grid case')
ax=axes[0,0];ax.plot(t,u*1000,color='#333333',lw=1.8,label='GENROU speed deviation')
for tau,col in zip(taus,['SignalDelay_speed_1ms_z','SignalDelay_speed_100ms_z','SignalDelay_speed_500ms_z']):
 ax.plot(t,a[col]*1000,lw=1.1,label=fr'$\tau={tau:g}$ s')
ax.set(ylabel='Speed deviation (10⁻³ pu)',title='TwoBusBasic + smooth mechanical power p(t) = cos(t), no switching events');ax.legend(ncol=4,fontsize=9)
ax=axes[1,0]
for tau,col in zip(taus,['SignalDelay_speed_1ms_z','SignalDelay_speed_100ms_z','SignalDelay_speed_500ms_z']):
 ref=np.interp(t-tau,t,u,left=u[0]);ax.plot(t,(a[col]-ref)*1e9,label=fr'$\tau={tau:g}$ s')
ax.set(xlabel='Time (s)',ylabel='Shift consistency difference (10⁻⁹ pu)',title='Compare delayed output with the same run’s recorded source shifted in time');ax.legend(ncol=3)
save(fig,'01_grid','Delayed GENROU speed','Actual DynamicSimulation CSV. Relative tolerance 10⁻⁸; monitor interval 0.5 ms.\nThe lower panel is a same-run consistency check, not an independent solution-error estimate. All three delays are integer multiples of the monitor interval. The source port is speed deviation (omega), not absolute speed (1 + omega).')

# 2: Application example with analytic reference and constant prehistory.
a=app('app');t=a['t'];fig,axes=layout('02  The JSON delay declarations, with an analytic reference')
axes[0,0].plot(t,np.cos(t),color='#333333',lw=1.5,label='Source cos(t), analytic')
for tau,col in zip([.001,.15,1.5],['SignalDelay_short_z','SignalDelay_long_z','SignalDelay_wide_z']):
 ref=np.cos(np.maximum(0,t-tau));axes[0,0].plot(t,a[col],label=fr'GridKit $\tau={tau:g}$ s')
 axes[1,0].plot(t,(a[col]-ref)*1e7,label=fr'$\tau={tau:g}$ s')
axes[0,0].set(ylabel='Signal value',title='Three taps share one source history; prehistory is constant at 1');axes[0,0].legend(ncol=4,fontsize=9)
axes[1,0].set(xlabel='Time (s)',ylabel='Output − exact delay (10⁻⁷)',title='Exact reference: 1 before the delay elapses, then cos(t − τ)');axes[1,0].legend(ncol=3)
save(fig,'02_cosine','Delayed cosine from JSON','Actual DynamicSimulation CSV. Relative tolerance 10⁻⁷; monitor interval 2 ms.\nThe large delay is visible as a waveform shift. The 1 ms delay is visually indistinguishable at this scale; the next figure expands its difference from the source. Constant prehistory meets cos(t) with a continuous first derivative at zero.')

# 3: Inspect the small delay directly instead of hiding it underneath the source.
fig,axes=layout('03  Make a 1 ms delay visible: inspect the difference from its source')
for h in [.0001,.001,.01,.1]:
 name=f'fixed_tau0.001_h{h:g}';a=data(name);t=a['t'];sig=(a['u']-a['z'])/.001
 axes[0,0].plot(t,sig,lw=1.15,label=fr'$h/\tau={h/.001:g}$ ($h={h:g}$ s)')
 axes[1,0].plot(t,sig-(a['exact_u']-a['exact_z'])/.001,lw=1.1)
a=data('fixed_tau0.001_h0.1');axes[0,0].plot(a['t'],(a['exact_u']-a['exact_z'])/.001,'k--',lw=1.2,label='Exact difference quotient')
axes[0,0].set(ylabel=r'$(u-z)/\tau$',title='The quotient approaches u′(t), and exposes whether the small shift is represented');axes[0,0].legend(ncol=3,fontsize=9)
axes[1,0].set(xlabel='Time (s)',ylabel='Difference-quotient error',yscale='symlog',title='Same runs, including low-order startup; symlog scale expands the small errors');axes[1,0].set_yscale('symlog',linthresh=1e-7)
save(fig,'03_small_delay','A small delay is not just overlapping traces','Fixed-step diagnostic sweep; loose history tolerance allows the selected steps.\nAll curves use the actual GridKit source and delayed algebraic variables. The analytic reference is [cos(t) − cos(t − τ)]/τ, with analytic prehistory. The difference quotient magnifies errors by 1/τ; visual overlap of raw waveforms alone does not validate a small delay.')

# 4: Separate h/tau from resolution of the signal itself.
fig,axes=layout('04  Step / delay ratio does not determine accuracy by itself',1,2,height=6.5)
for tau in [.001,.01,.1,1.]:
 rows=sorted([r for r in ROWS if r['fixed_h'] and r['tau']==tau],key=lambda r:r['fixed_h'])
 err=[r['max_abs_z_error'] for r in rows]
 axes[0,0].loglog([r['fixed_h']/tau for r in rows],err,'o-',ms=4,label=fr'$\tau={tau:g}$ s')
 axes[0,1].loglog([r['fixed_h']/(2*np.pi) for r in rows],err,'o-',ms=4,label=fr'$\tau={tau:g}$ s')
axes[0,0].axvline(1,color='#555555',ls=':',lw=1);axes[0,0].set(xlabel=r'$h/\tau$',ylabel='Maximum sampled |z − exact delay|',title='Crossing h/τ = 1 creates no accuracy cliff')
axes[0,1].set(xlabel=r'$h/T$, with signal period $T=2\pi$ s',ylabel='Maximum sampled |z − exact delay|',title='Signal resolution and startup matter');axes[0,1].legend()
save(fig,'04_step_sweep','28 fixed-step runs','Error includes uniform monitor samples AND the midpoint of every accepted step.\nFour delays × seven steps, h = 0.0001 … 0.5 s. Analytic cos(t) prehistory; compression disabled. Existing GridKit fixed-step mode limits IDA to order 2 and disables its integration error test. History rtol=1, atol=0.01 intentionally avoids accuracy-triggered rejection. These runs measure error; they do not certify tolerance compliance. The final shortened step is included.')

# 5: Adaptive time steps around smooth pulse.
fig,axes=layout('05  Adaptive steps follow a smooth pulse, without a shortest-delay cap',3,1,height=9)
for tau in [.001,.1,.5]:
 name=f'pulse_tau{tau:g}';a=data(name);s=data(name,'steps');at=data(name,'attempts')
 axes[0,0].plot(a['t'],a['z'],label=fr'$\tau={tau:g}$ s')
 axes[1,0].step(s['t'],s['h_over_tau'],where='pre',lw=1,label=fr'$\tau={tau:g}$ s')
 axes[2,0].step(s['t'],s['h'],where='pre',lw=1)
 rej=at[at['rejected']>0];axes[2,0].scatter(rej['t'],rej['h'],marker='x',s=35)
axes[0,0].plot(a['t'],a['exact_u'],'k--',lw=1,label='Analytic source');axes[0,0].set(ylabel='Signal value',title='u(t) = 1 + 0.6 exp[−((t − 2)/0.12)²]');axes[0,0].legend(ncol=4,fontsize=9)
axes[1,0].axhline(1,color='k',ls=':',lw=1);axes[1,0].set(yscale='log',ylabel=r'Accepted $h/\tau$',title='Dashed horizontal line: h = τ')
axes[2,0].set(xlabel='Time (s)',ylabel='Step size h (s)',yscale='log',title='Accepted steps; crosses show rejected candidate steps')
save(fig,'05_adaptive_pulse','Adaptive resolution of a smooth pulse','Relative tolerance 10⁻⁶; absolute tolerance 10⁻⁸; analytic prehistory; compression off.\nFor τ=1 ms, h reaches 0.512 s in a flat part but is at most 4.27 ms near the source pulse (t=1.8…2.2 s). This is h/τ up to 512, not a claim that the pulse itself is resolved by 0.512 s steps. Delayed outputs also affect integration accuracy and therefore the step sequence.')

# 6: Coupled differential equation, not only prescribed inputs.
fig,axes=layout('06  Delay inside feedback: check against an independent exact solution',2,2,height=8)
t=np.linspace(0,2,1000);axes[0,0].plot(t,np.exp(-.4*(t-.001)),'k--',label='Exact delayed state')
for tol in [1e-3,1e-5,1e-7]:
 name=f'feedback_{tol:g}';a=data(name);s=data(name,'steps')
 axes[0,0].plot(a['t'],a['z'],label=fr'rtol = {tol:g}')
 axes[0,1].plot(a['t'],a['z']-a['exact_z'],label=fr'rtol = {tol:g}')
 axes[1,0].step(s['t'],s['h_over_tau'],where='pre',lw=1,label=fr'rtol = {tol:g}')
r=sorted([r for r in ROWS if r['kind']=='feedback'],key=lambda r:r['rtol'])
axes[1,1].loglog([r['rtol'] for r in r],[r['max_abs_z_error'] for r in r],'o-',label='Measured delayed-output error')
axes[0,0].set(ylabel='Delayed state',title='x′ = −γx + 1.3x(t − τ), τ = 1 ms');axes[0,0].legend(fontsize=9)
axes[0,1].set(ylabel='z − exact delayed state',title='Exact x(t) = exp(−0.4t)',yscale='symlog');axes[0,1].set_yscale('symlog',linthresh=1e-9)
axes[1,0].axhline(1,color='k',ls=':',lw=1);axes[1,0].set(xlabel='Time (s)',ylabel=r'Accepted $h/\tau$',yscale='log')
axes[1,1].set(xlabel='Requested relative tolerance',ylabel='Maximum sampled delay error',title='Global error need not vary monotonically')
save(fig,'06_feedback','A true feedback connection','γ = 1.3 exp(0.4τ) + 0.4; exact exponential prehistory; compression off.\nAt rtol=10⁻⁷: maximum h/τ = 23.402, maximum sampled delayed-output error = 9.35×10⁻⁹. Tightening to 10⁻⁸ gives 1.59×10⁻⁸ in this run: local error control is not a strict global-error or monotonic-convergence guarantee. Six tolerances were run; all are retained in summary.csv.')

# 7: Optional compression evaluated honestly.
a=data('compression_0');b=data('compression_1');fig,axes=layout('07  Does history compression pay for itself here?',3,1,height=9)
axes[0,0].plot(a['t'],a['z'],label='Compression off');axes[0,0].plot(b['t'],b['z'],'--',label='Compression on');axes[0,0].set(ylabel='Delayed pulse',title='Same pulse, τ = 2 s; only the compression option changes');axes[0,0].legend()
for name,label in [('compression_0','off'),('compression_1','on')]:
 s=data(name,'steps');r=BY[name];axes[1,0].step(s['t'],s['segments'],where='post',label=f'{label}: peak {r["peak_segments"]} segments, {r["accepted"]} accepted steps')
axes[1,0].set(ylabel='Retained cubic segments',title='Window pruning remains enabled in both runs');axes[1,0].legend(fontsize=9)
axes[2,0].plot(a['t'],(b['z']-a['z'])*1e6);axes[2,0].set(xlabel='Time (s)',ylabel='On − off (10⁻⁶)',title='Trajectory difference; adaptive step sequences differ')
save(fig,'07_compression','Compression versus retention','Peak segments: 324 → 289 (10.8% reduction). Accepted steps: 538 → 564.\nSegment counts are not process-memory measurements. No speedup is established. Compression has a bounded change to stored polynomials; that is not a bound on the difference between two complete adaptive trajectories. Both runs keep a bounded delay window even with compression disabled.')

# 8: Increasing frequency makes restart costs visible.
fig,axes=layout('08  A smooth chirp exposes the cost of history-triggered retries',3,1,height=9)
a=data('chirp_1e-06');axes[0,0].plot(a['t'],a['z'],label='GridKit delayed signal');axes[0,0].plot(a['t'],a['exact_z'],'k--',lw=.9,label='Exact delay');axes[0,0].set(ylabel='Signal value',title='u(t) = cos(t + 0.4t²), τ = 1 ms');axes[0,0].legend()
for tol in [1e-4,1e-6,1e-8]:
 name=f'chirp_{tol:g}';s=data(name,'steps');at=data(name,'attempts');r=BY[name]
 axes[1,0].step(s['t'],s['h_over_tau'],where='pre',lw=.8,label=f'rtol={tol:g}: {r["accepted"]} accepted, {r["rejected"]} rejected')
 rej=at[at['rejected']>0];axes[2,0].step(rej['t'],np.arange(1,len(rej)+1),where='post',label=f'rtol={tol:g}')
axes[1,0].axhline(1,color='k',ls=':',lw=1);axes[1,0].set(ylabel=r'Accepted $h/\tau$',yscale='log');axes[1,0].legend(fontsize=9)
axes[2,0].set(xlabel='Time (s)',ylabel='Cumulative history rejections',title='Each rejection currently calls IDAReInit, discarding BDF history')
save(fig,'08_chirp','Retries on an increasingly oscillatory source','Smooth analytic prehistory; compression off; matched solver/history tolerances.\nAt rtol=10⁻⁶: 1,933 accepted steps, 42 history rejections, maximum sampled delay error 5.37×10⁻⁷. Rejections are the additional history-controller rejections, not IDA’s internal error-test or nonlinear-solve failures. This is evidence to improve the step predictor, not a runtime speed comparison.')
pdf.close()
# Self-contained HTML supports copying the gallery without broken image links.
intro='''<h1>GridKit delayed signals: measured results</h1><p>44 runs: two actual <code>DynamicSimulation</code> application cases and 42 experiments linked to the current GridKit libraries. All colored numerical curves are GridKit outputs. Dashed analytic curves and shifted-source comparisons are labeled.</p><p><strong>h is an accepted internal solver step.</strong> The monitor interval is a separate quantity. Large h/τ is possible; source resolution, interpolation and nonlinear/integration accuracy still limit h.</p><p><a href="gallery.pdf">Eight-page PDF</a> · <a href="summary.csv">42-run numerical summary</a> · <a href="README.md">Methods and reproduction</a> · <a href="provenance.json">Provenance</a></p>'''
parts=['<!doctype html><html lang="en"><meta charset="utf-8"><title>GridKit signal delay gallery</title><style>body{font:16px/1.55 system-ui,sans-serif;color:#202a35;max-width:1180px;margin:36px auto;padding:0 24px;background:#f6f8fb}h1{font-size:32px}h2{font-size:24px}a{color:#12649a}figure{margin:34px 0;padding:22px;background:white;border:1px solid #dce2e9;border-radius:9px}img{width:100%;height:auto}figcaption{color:#334155}nav{display:flex;gap:15px;flex-wrap:wrap}code{background:#e9edf3;padding:2px 5px}</style>',intro,'<nav>']
parts.extend(f'<a href="#{name}">{i+1}. {html.escape(title)}</a>' for i,(name,title,_) in enumerate(PAGES));parts.append('</nav>')
for name,title,caption in PAGES:
 encoded=base64.b64encode((OUT/(name+'.png')).read_bytes()).decode()
 parts.append(f'<figure id="{name}"><h2>{html.escape(title)}</h2><img alt="{html.escape(title)}" src="data:image/png;base64,{encoded}"><figcaption>{html.escape(caption).replace(chr(10),"<br><br>")}<p><a href="figures/{name}.png">PNG</a> · <a href="figures/{name}.svg">SVG</a></p></figcaption></figure>')
parts.append('</html>');(ROOT/'index.html').write_text('\n'.join(parts))
print(f'Wrote {len(PAGES)} figures, gallery.pdf and index.html in {ROOT}')
