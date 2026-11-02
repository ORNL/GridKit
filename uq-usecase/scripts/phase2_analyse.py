"""Phase 2 analysis: CCT distribution, drivers, and the peak-screening test."""

import numpy as np
import pandas as pd

R = pd.read_csv("/tmp/phase2/phase2_results.csv")
R = R[R.status == "ok"].copy()
R["cens"] = R.cct_censored.fillna(False).astype(bool)

print(
    f"{len(R)} conditions | censored at 0.5 s: {R.cens.sum()} ({100*R.cens.mean():.0f}%)\n"
)

print("CCT by net-load bin (censored = stable beyond any credible clearing time):")
g = (
    R.groupby("bin", observed=True)
    .agg(
        n=("cct_s", "size"),
        net_mean=("net_mw", "mean"),
        cct_min=("cct_s", "min"),
        cct_med=("cct_s", "median"),
        censored=("cens", "sum"),
        inertia_mean=("inertia_MWs", "mean"),
    )
    .round(3)
)
print(g.to_string())

unc = R[~R.cens]
print(f"\non UNCENSORED only (n={len(unc)}, biased toward the stressed tail):")
print(f"  corr(CCT, net load) = {unc.cct_s.corr(unc.net_mw):+.3f}")
print(f"  corr(CCT, inertia)  = {unc.cct_s.corr(unc.inertia_MWs):+.3f}")

spread = R.groupby("hour").cct_s.agg(["min", "max", "nunique"])
spread["range"] = spread["max"] - spread["min"]
print("\ndoes commitment change CCT at fixed hour?")
print(
    f"  hours where CCT is IDENTICAL across all commitments: "
    f"{int((spread['nunique'] == 1).sum())} of {len(spread)}"
)
print(f"  median CCT range within an hour: {spread['range'].median():.4f} s")
print(f"  max CCT range within an hour:    {spread['range'].max():.4f} s")

peak = R.loc[R.net_mw.idxmax()]
worse = R[R.cct_s < peak.cct_s]
print("\nPEAK-LOAD SCREENING TEST")
print(f"  highest-net-load condition: net={peak.net_mw:.0f} MW, CCT={peak.cct_s:.4f} s")
print(
    f"  conditions with CCT strictly BELOW that: {len(worse)} of {len(R)} "
    f"({100*len(worse)/len(R):.0f}%)"
)
print(
    f"  worst CCT overall: {R.cct_s.min():.4f} s "
    f"at net={R.loc[R.cct_s.idxmin(),'net_mw']:.0f} MW"
)

print("\nmargin against protection:")
for tc in (0.08, 0.10, 0.25):
    print(
        f"  t_clear={tc:.2f}s -> P(CCT < t_clear) = {100*(R.cct_s < tc).mean():5.1f}%"
        f" | worst-case margin {R.cct_s.min()/tc:.2f}x"
    )
