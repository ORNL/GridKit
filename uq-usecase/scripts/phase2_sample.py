"""Phase 2: aleatoric baseline. Stratified sample over net load x commitment.

Design (see work-notes/datacenter_load_uq_plan.md section 7):

* net load is stratified into NET_BINS quantile bins of the empirical 8760
  joint (demand, wind), with HOURS_PER_BIN hours drawn from each
* commitment is NOT assumed to follow a distribution. At each hour it is
  stratified at fixed fractions of the measured feasible band
  [n_min, n_max], so results can be reported conditionally and reweighted
  later if a real UC distribution becomes available

Writes incrementally so a kill does not lose completed runs, and skips
conditions already present on resume.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp

DATA = "/kfs2/projects/scidac/scidac-data/ACTIVSg200/grid-v1"
WORK = Path("/tmp/phase2")
RESULTS = WORK / "phase2_results.csv"

NET_BINS = 5
HOURS_PER_BIN = 4
COMMIT_FRACTIONS = (0.0, 0.25, 0.50, 0.80)  # position within the feasible band
SEED = 0

cfg = cp.CampaignConfig(work=WORK)
order = cp.merit_order(cfg)

dem = pd.read_csv(f"{DATA}/demand.csv")
wnd = pd.read_csv(f"{DATA}/wind.csv")
DEMC = [c for c in dem.columns if c != "DateTime"]
WNDC = [c for c in wnd.columns if c != "DateTime"]
WG = {c: cfg.gid[int(c.replace("gen", "")) - 1] for c in WNDC}

net_all = dem[DEMC].sum(axis=1) - wnd[WNDC].sum(axis=1)
bins = pd.qcut(net_all, NET_BINS, labels=[f"q{i}" for i in range(NET_BINS)])

rng = np.random.default_rng(SEED)
conditions = []
for b in [f"q{i}" for i in range(NET_BINS)]:
    pool = net_all.index[bins == b].to_numpy()
    for h in rng.choice(pool, size=HOURS_PER_BIN, replace=False):
        h = int(h)
        lo, hi = cp.commitment_band(cfg, float(net_all[h]), order)
        if lo is None:
            print(f"  hour {h}: no feasible commitment, skipped")
            continue
        for f in COMMIT_FRACTIONS:
            n = int(round(lo + f * (hi - lo)))
            conditions.append({"bin": b, "hour": h, "frac": f, "n_commit": n})

C = pd.DataFrame(conditions).drop_duplicates(subset=["hour", "n_commit"])
print(f"{len(C)} conditions over {C.hour.nunique()} hours and {NET_BINS} net-load bins")
print(
    C.groupby("bin", observed=True)
    .agg(
        hours=("hour", "nunique"),
        runs=("n_commit", "size"),
        n_min=("n_commit", "min"),
        n_max=("n_commit", "max"),
    )
    .to_string()
)
print()

done = set()
if RESULTS.exists():
    prev = pd.read_csv(RESULTS)
    done = set(zip(prev.hour, prev.n_commit))
    print(f"resuming, {len(done)} conditions already done")

rows = []
for i, c in C.iterrows():
    key = (c.hour, c.n_commit)
    if key in done:
        continue
    tag = f"h{c.hour}_n{c.n_commit}"
    demand = {int(x.replace("bus", "")): float(dem.loc[c.hour, x]) for x in DEMC}
    wind = {WG[x]: float(wnd.loc[c.hour, x]) for x in WNDC}
    rec = cp.run_condition(cfg, tag, demand, wind, int(c.n_commit))
    rec.update({"bin": c.bin, "hour": int(c.hour), "frac": c.frac})
    rows.append(rec)
    pd.DataFrame(rows).to_csv(
        RESULTS,
        mode="a" if RESULTS.exists() else "w",
        header=not RESULTS.exists(),
        index=False,
    )
    rows = []
    print(
        f"  {tag:>14}  {rec['status']:<12} inertia={rec.get('inertia_MWs')}"
        f"  cct={rec.get('cct_s')}{'*' if rec.get('cct_censored') else ''}"
        f"  nadir={rec.get('coi_nadir_Hz')}"
        f"  {rec['wall_s']}s",
        flush=True,
    )

R = pd.read_csv(RESULTS)
print(
    f"\n{len(R)} total runs | ok={int((R.status=='ok').sum())} "
    f"| solve_failed={int((R.status=='solve_failed').sum())} "
    f"| sim_failed={int((R.status=='sim_failed').sum())}"
)
ok = R[R.status == "ok"]
if len(ok):
    print("\nQoI spread over the sample:")
    print(
        ok[
            [
                "inertia_MWs",
                "cct_s",
                "coi_nadir_Hz",
                "coi_max_rocof_Hz_s",
                "Vmin_post",
                "max_delta_spread_deg",
            ]
        ]
        .describe()
        .round(4)
        .to_string()
    )
    if "cct_censored" in ok:
        n_cens = int(ok.cct_censored.fillna(False).astype(bool).sum())
        print(
            f"\nCCT right-censored at {cfg.cct_max}s: {n_cens} of {len(ok)} "
            f"({100*n_cens/len(ok):.0f}%) -- lower bounds, not point estimates"
        )
