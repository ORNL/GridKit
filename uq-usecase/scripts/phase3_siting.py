"""Phase 3: datacenter siting. How much load can each bus take before the
transient stability margin (CCT) degrades?

Design deliberately small, to produce a result at modest cost:

    candidate buses  x  datacenter sizes  x  operating conditions

Candidate buses are stratified by hop distance from the fault bus, so the
sweep covers electrically near and far siting rather than clustering.
Operating conditions are drawn from the Phase 2 sample: one low, one mid, one
high net load, each at its own feasible commitment floor.

Baseline (dc_mw = 0) is run for each condition so the reported quantity is the
CHANGE in CCT caused by the datacenter, not its absolute value.
"""

import sys
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp

DATA = "/kfs2/projects/scidac/scidac-data/ACTIVSg200/grid-v1"
WORK = Path("/tmp/phase3")
RESULTS = WORK / "phase3_results.csv"

DC_SIZES_MW = (150.0, 400.0)
N_BUSES = 8
SEED = 0

cfg = cp.CampaignConfig(work=WORK)
order = cp.merit_order(cfg)

dem = pd.read_csv(f"{DATA}/demand.csv")
wnd = pd.read_csv(f"{DATA}/wind.csv")
DEMC = [c for c in dem.columns if c != "DateTime"]
WNDC = [c for c in wnd.columns if c != "DateTime"]
WG = {c: cfg.gid[int(c.replace("gen", "")) - 1] for c in WNDC}
net_all = dem[DEMC].sum(axis=1) - wnd[WNDC].sum(axis=1)

# --- operating conditions: low / mid / high net load ---
q = net_all.quantile([0.10, 0.50, 0.90]).values
HOURS = {
    lbl: int((net_all - v).abs().idxmin()) for lbl, v in zip(("lo", "mid", "hi"), q)
}

# --- candidate buses, stratified by hop distance from the fault ---
adj = {}
for d in cfg.base_case["devices"]:
    if d["class"] == "Branch":
        a, b = d["ports"]["bus1"], d["ports"]["bus2"]
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
dist, dq = {cfg.fault_bus: 0}, deque([cfg.fault_bus])
while dq:
    u = dq.popleft()
    for v in adj.get(u, ()):
        if v not in dist:
            dist[v] = dist[u] + 1
            dq.append(v)

load_buses = {
    d["ports"]["bus"] for d in cfg.base_case["devices"] if d["class"] == "LoadZIP"
}
cand = pd.DataFrame(
    [{"bus": b, "hops": dist.get(b, 99)} for b in sorted(load_buses) if b in dist]
)
rng = np.random.default_rng(SEED)
bands = pd.qcut(cand.hops, 4, labels=False, duplicates="drop")
picks = []
for band in sorted(pd.unique(bands)):
    pool = cand[bands == band].bus.to_numpy()
    k = min(N_BUSES // 4, len(pool))
    picks += list(rng.choice(pool, size=k, replace=False))
BUSES = sorted(int(b) for b in picks)

print(f"operating conditions: {HOURS}")
print(
    f"candidate buses (hops from fault bus {cfg.fault_bus}): "
    f"{[(b, dist[b]) for b in BUSES]}"
)

conditions = []
for lbl, h in HOURS.items():
    conditions.append({"cond": lbl, "hour": h, "bus": None, "mw": 0.0})
    for b in BUSES:
        for mw in DC_SIZES_MW:
            conditions.append({"cond": lbl, "hour": h, "bus": b, "mw": mw})
C = pd.DataFrame(conditions)
print(
    f"\n{len(C)} runs = {len(HOURS)} conditions x "
    f"({len(BUSES)} buses x {len(DC_SIZES_MW)} sizes + 1 baseline)"
)

done = set()
if RESULTS.exists():
    prev = pd.read_csv(RESULTS)
    done = set(zip(prev.cond, prev.dc_bus.fillna(-1), prev.dc_mw))
    print(f"resuming, {len(done)} already done")

for _, c in C.iterrows():
    # pandas stores a None bus as NaN, so test with isna, not `is None`
    is_base = pd.isna(c.bus)
    key = (c.cond, -1.0 if is_base else float(c.bus), float(c.mw))
    if key in done:
        continue
    h = int(c.hour)
    demand = {int(x.replace("bus", "")): float(dem.loc[h, x]) for x in DEMC}
    wind = {WG[x]: float(wnd.loc[h, x]) for x in WNDC}
    # commitment must resize with the datacenter, so recompute the band
    net = sum(demand.values()) + float(c.mw) - sum(wind.values())
    lo, hi = cp.commitment_band(cfg, net, order)
    if lo is None:
        print(f"  {c.cond} bus={c.bus} mw={c.mw}: no feasible commitment")
        continue
    tag = f"{c.cond}_b{0 if is_base else int(c.bus)}_mw{int(c.mw)}"
    rec = cp.run_condition(
        cfg,
        tag,
        demand,
        wind,
        lo,
        dc_bus=None if is_base else int(c.bus),
        dc_mw=float(c.mw),
    )
    rec.update(
        {
            "cond": c.cond,
            "hour": h,
            "hops": 0 if is_base else dist.get(int(c.bus)),
        }
    )
    pd.DataFrame([rec]).to_csv(
        RESULTS,
        mode="a" if RESULTS.exists() else "w",
        header=not RESULTS.exists(),
        index=False,
    )
    print(
        f"  {tag:>18} {rec['status']:<12} cct={rec.get('cct_s')}"
        f"{'*' if rec.get('cct_censored') else ''}  n={rec.get('n_commit')}"
        f"  {rec['wall_s']}s",
        flush=True,
    )

R = pd.read_csv(RESULTS)
ok = R[R.status == "ok"]
print(f"\n{len(R)} runs | ok={len(ok)}")
if len(ok):
    base = ok[ok.dc_mw == 0].set_index("cond")["cct_s"].to_dict()
    ok = ok.copy()
    ok["cct_base"] = ok["cond"].map(base)
    ok["d_cct"] = ok.cct_s - ok.cct_base
    print("\nchange in CCT caused by the datacenter (negative = worse):")
    print(
        ok[ok.dc_mw > 0]
        .pivot_table(
            index=["dc_bus", "hops"], columns=["cond", "dc_mw"], values="d_cct"
        )
        .round(4)
        .to_string()
    )
