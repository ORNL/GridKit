"""Results for the two main campaigns used in report section 7.

Main campaigns (both with PV->PQ Q-limit switching, alphaP 0.8):
  aleatoric1_ap08_qlim   20 hours x 2 commitments x 4 buses x 5 sizes + baseline, fault 66
  faultscreen1_ap08_qlim 10 hours x 2 commitments x 17 faults x {baseline, bus 108, bus 148} at 200 MW
  aleatoric2_ap08_qlim   60 hours x 3 commitments x 4 buses x 5 sizes + baseline, fault 66 (figures 7.2b, 7.3)

Outcome: censored Kaplan-Meier over runs with status 'ok'. CCT is censored at the
0.5 s cap (cct_censored == True). Runs that failed (sim_failed, solve_failed) carry
no outcome and are reported as counts only.

Usage:
  python scripts/main_results.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/kfs2/projects/scidac/scidac-data/gridkit-campaigns")
FIG = Path("/home/isatkaus/gridkit/uq-usecase/work-notes/figures")
CAP = 0.5
THRESH = 0.25
rng = np.random.default_rng(0)


def load(name):
    p = ROOT / name / "results.parquet"
    if p.exists():
        R = pd.read_parquet(p)
    else:
        # Same table as sweep_collect.py, built in memory (nothing written to the store).
        R = pd.DataFrame([json.loads(f.read_text()) for f in sorted((ROOT / name / "records").glob("*.json"))])
    R["cens"] = R["cct_censored"].fillna(False).astype(bool)
    return R


def km(times, events):
    """Kaplan-Meier. times: CCT (s); events: True if CCT observed (not censored at CAP).
    Returns (t_grid, S(t)) with S = P(CCT > t)."""
    t = np.asarray(times, float)
    e = np.asarray(events, bool)
    order = np.argsort(t)
    t, e = t[order], e[order]
    uniq = np.unique(t[e]) if e.any() else np.array([])
    S = 1.0
    grid, surv = [0.0], [1.0]
    n_at_risk = len(t)
    for u in uniq:
        d = np.sum((t == u) & e)
        n_at_risk_u = np.sum(t >= u)
        S *= 1.0 - d / n_at_risk_u
        grid.append(u)
        surv.append(S)
    return np.array(grid), np.array(surv)


def plot_km(ax, times, cens, color, label):
    """KM step curve drawn out to the cap. Censored runs sit at CAP, so the curve
    holds its last value up to CAP; the tick marks the censored share."""
    cens = np.asarray(cens, bool)
    g, s = km(times, ~cens)
    g, s = np.append(g, CAP), np.append(s, s[-1])
    ax.step(g, s, where="post", color=color,
            label=f"{label}, {cens.mean():.0%} censored")
    ax.plot([CAP], [s[-1]], marker="|", ms=12, mew=2, color=color, ls="none")
    ax.set_xlim(0, CAP * 1.04)
    ax.set_ylim(0, 1)


def p_below(times, events, thresh=THRESH):
    """P(CCT < thresh) from the KM estimate."""
    g, s = km(times, events)
    return 1.0 - s[np.searchsorted(g, thresh, side="right") - 1]


def greenwood(times, events, thresh=THRESH):
    """P(CCT < thresh) from the KM estimate, with Greenwood's standard error."""
    t = np.asarray(times, float)
    e = np.asarray(events, bool)
    S, var = 1.0, 0.0
    for u in np.unique(t[e & (t <= thresh)]):
        d = np.sum((t == u) & e)
        n = np.sum(t >= u)
        S *= 1.0 - d / n
        var += d / (n * (n - d))
    return 1.0 - S, S * np.sqrt(var)


def km_median(times, events):
    g, s = km(times, events)
    idx = np.where(s <= 0.5)[0]
    return g[idx[0]] if len(idx) else np.nan  # nan = not reached before cap


def sample_space(name, R):
    print(f"\n== {name}")
    print("runs:", len(R), "| status:", R.status.value_counts().to_dict())
    base = R[R.dc_bus.isna()]
    sited = R[R.dc_bus.notna()]
    print("baseline runs:", len(base), base.status.value_counts().to_dict())
    print("sited runs:", len(sited), sited.status.value_counts().to_dict())
    print("hours:", R.hour.nunique(), "| bins:", R.groupby("bin").hour.nunique().to_dict())
    print("faults:", R.fault_bus.nunique(), "| placements:", sorted(R.dc_bus.dropna().unique().astype(int).tolist()))
    print("MW levels:", sorted(R.dc_mw.unique().tolist()))


def km_summary(label, R):
    ok = R[R.status == "ok"]
    ev = ~ok.cens
    print(f"{label}: n_ok={len(ok)}, censored={ok.cens.sum()} ({ok.cens.mean():.1%}), "
          f"P(CCT<{THRESH})={p_below(ok.cct_s, ev):.3f}, KM median={km_median(ok.cct_s, ev)}")


def aleatoric_block():
    A = load("aleatoric1_ap08_qlim")
    sample_space("aleatoric1_ap08_qlim", A)
    ok = A[A.status == "ok"]
    base_ok = ok[ok.dc_bus.isna()]
    sit_ok = ok[ok.dc_bus.notna()]
    km_summary("baseline", base_ok)
    km_summary("sited (all buses, all sizes)", sit_ok)

    # Simple KM figure: baseline vs sited at bus 108, 200 MW (fault 66 only)
    b108 = A[(A.dc_bus == 108) & (A.dc_mw == 200)]
    okb = b108[b108.status == "ok"]
    key = ["hour", "n_commit"]
    bm = A[A.dc_bus.isna()].merge(okb[key].drop_duplicates(), on=key)
    print("\nbus 108 200 MW: runs", len(b108), b108.status.value_counts().to_dict())
    km_summary("matched baseline (same conditions as bus 108 ok)", bm)
    km_summary("bus 108, 200 MW", okb)
    fig, ax = plt.subplots(figsize=(6, 4))
    plot_km(ax, bm.cct_s, bm.cens, "C0", f"baseline, same conditions, n={len(bm)}")
    plot_km(ax, okb.cct_s, okb.cens, "C1", f"bus 108, 200 MW, n={len(okb)}")
    ax.axvline(THRESH, ls=":", color="grey")
    ax.set_xlabel("CCT (s)")
    ax.set_ylabel("P(CCT > t)")
    ax.set_title("aleatoric1_ap08_qlim, fault 66, switching at alphaP 0.8")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "main_aleatoric_km_simple.png", dpi=150)
    plt.close(fig)

    # Net-load quintile KM (sited arm, pooled over buses and sizes)
    print("\nby net-load quintile (sited ok / baseline ok / status of sited):")
    fig, ax = plt.subplots(figsize=(6, 4))
    for i, b in enumerate(sorted(A.bin.unique())):
        s_ = sit_ok[sit_ok.bin == b]
        mean_net = A[A.bin == b].net_mw.mean()
        plot_km(ax, s_.cct_s, s_.cens, f"C{i}", f"{b}: mean net {mean_net:,.0f} MW, n={len(s_)}")
        print(f"  {b} mean net load {mean_net:,.1f} MW, net range "
              f"{A[A.bin == b].net_mw.min():,.1f}-{A[A.bin == b].net_mw.max():,.1f}")
    ax.axvline(THRESH, ls=":", color="grey")
    ax.set_xlabel("CCT (s)")
    ax.set_ylabel("P(CCT > t)")
    ax.set_title("aleatoric1_ap08_qlim, sited runs by net-load quintile")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "main_aleatoric_km_by_netload.png", dpi=150)
    plt.close(fig)
    print("no-outcome share by quintile:",
          (A[A.dc_bus.notna()].groupby("bin").status.apply(lambda x: (x != "ok").mean())).round(3).to_dict())
    print("sited dispatch failures by quintile:",
          A[A.dc_bus.notna() & (A.status == "solve_failed")].groupby("bin").size().to_dict())
    for b in sorted(A.bin.unique()):
        s = sit_ok[sit_ok.bin == b]
        bs = base_ok[base_ok.bin == b]
        sa = A[(A.bin == b) & A.dc_bus.notna()]
        print(f"  {b}: sited ok {len(s)}/{len(sa)} cens {s.cens.mean():.1%} "
              f"P(<{THRESH})={p_below(s.cct_s, ~s.cens):.3f} | baseline ok {len(bs)} "
              f"P(<{THRESH})={p_below(bs.cct_s, ~bs.cens) if len(bs) else np.nan:.3f}")

    # Damage at matched conditions: sited minus baseline, same hour, commitment, fault 66
    key = ["hour", "n_commit"]
    bb = A[A.dc_bus.isna()][key + ["cct_s", "cens", "status"]].rename(
        columns={"cct_s": "cct_base", "cens": "cens_base", "status": "st_base"})
    m = A[A.dc_bus.notna()].merge(bb, on=key)
    m = m[(m.status == "ok") & (m.st_base == "ok")]
    m = m[~m.cens & ~m.cens_base]
    m["dcct"] = m.cct_s - m.cct_base
    print("\nmatched damage (both measured), by MW:")
    print(m.groupby("dc_mw").dcct.agg(["count", "median", "min", "max"]).round(4))
    print("by bus:")
    print(m.groupby("dc_bus").dcct.agg(["count", "median", "min", "max"]).round(4))
    return A


def faultscreen_block():
    F = load("faultscreen1_ap08_qlim")
    sample_space("faultscreen1_ap08_qlim", F)
    ok = F[F.status == "ok"]
    km_summary("baseline", ok[ok.dc_bus.isna()])
    for b in [108, 148]:
        km_summary(f"bus {b} 200 MW", ok[ok.dc_bus == b])

    key = ["hour", "n_commit", "fault_bus"]
    base = F[F.dc_bus.isna()][key + ["cct_s", "cens", "status"]].rename(
        columns={"cct_s": "cct_base", "cens": "cens_base", "status": "st_base"})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, b in zip(axes, [108, 148]):
        sub = ok[ok.dc_bus == b]
        bm = ok[ok.dc_bus.isna()].merge(sub[key].drop_duplicates(), on=key)
        plot_km(ax, bm.cct_s, bm.cens, "C0", f"baseline, same conditions, n={len(bm)}")
        plot_km(ax, sub.cct_s, sub.cens, "C1", f"bus {b}, 200 MW, n={len(sub)}")
        ax.axvline(THRESH, ls=":", color="grey")
        ax.set_xlabel("CCT (s)")
        ax.set_title(f"faultscreen1_ap08_qlim, bus {b}")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("P(CCT > t)")
    fig.tight_layout()
    fig.savefig(FIG / "main_faultscreen_km.png", dpi=150)
    plt.close(fig)

    print("\nmatched pairs (same hour, commitment, fault) by placement:")
    for b in [108, 148]:
        m = F[F.dc_bus == b][key + ["cct_s", "cens", "status"]].merge(base, on=key)
        ok_pair = m[(m.status == "ok") & (m.st_base == "ok")]
        n_pairs = len(ok_pair)
        meas = ok_pair[~ok_pair.cens & ~ok_pair.cens_base]
        dcct = meas.cct_s - meas.cct_base
        print(f"  bus {b}: pairs with both ok {n_pairs}, both measured {len(meas)}, "
              f"sited-only failed {(m.status != 'ok').sum()}, median dCCT {dcct.median():.4f}, "
              f"min {dcct.min():.4f}, max {dcct.max():.4f}, "
              f"fraction dCCT<0 {(dcct < 0).mean():.3f}")

    print("\nper-fault baseline CCT (ok runs):")
    bo = ok[ok.dc_bus.isna()]
    tab = bo.groupby("fault_bus").agg(n=("cct_s", "size"), median=("cct_s", "median"),
                                      cens=("cens", "mean"))
    print(tab.round(3).to_string())
    return F


def pool_108(name):
    """Bus 108, 200 MW, ok runs, and the matched baseline on the same (hour, commitment) conditions."""
    R = load(name)
    sample_space(name, R)
    key = ["hour", "n_commit"]
    dc = R[(R.dc_bus == 108) & (R.dc_mw == 200) & (R.status == "ok")].reset_index(drop=True)
    base = R[R.dc_bus.isna()].merge(dc[key].drop_duplicates(), on=key)
    return dc, base


def km_37_vs_171_block():
    """Figure 7.2b: matched KM curves at 37 conditions (aleatoric1) and 171 (aleatoric2)."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, name in zip(axes, ["aleatoric1_ap08_qlim", "aleatoric2_ap08_qlim"]):
        dc, base = pool_108(name)
        km_summary(f"{name} matched baseline", base)
        km_summary(f"{name} bus 108, 200 MW", dc)
        for lab, R in [("baseline", base), ("bus 108, 200 MW", dc)]:
            p, se = greenwood(R.cct_s, ~R.cens)
            print(f"  {lab}: P(CCT<{THRESH})={p:.3f} Greenwood 95% [{p - 1.96 * se:.3f}, {p + 1.96 * se:.3f}]")
        plot_km(ax, base.cct_s, base.cens, "C0", f"baseline, n={len(base)}")
        plot_km(ax, dc.cct_s, dc.cens, "C1", f"bus 108, 200 MW, n={len(dc)}")
        ax.axvline(THRESH, ls=":", color="grey")
        ax.set_xlabel("CCT (s)")
        ax.set_title(f"{name}, fault 66")
        ax.legend(fontsize=8, loc="upper right")
    axes[0].set_ylabel("P(CCT > t)")
    fig.tight_layout()
    fig.savefig(FIG / "main_km_37_vs_171.png", dpi=150)
    plt.close(fig)


def random_orderings(pool, ns, B=2000, seed=0):
    """Random orderings of the pool (sampling without replacement). For each n, the KM estimate
    of P(CCT < THRESH) from the first n runs of each ordering. Returns the rows and the
    full-pool estimate."""
    cct = pool.cct_s.values
    ev = ~pool.cens.values
    rng = np.random.default_rng(seed)
    perms = [rng.permutation(len(pool)) for _ in range(B)]
    rows = []
    for n in ns:
        v = np.array([p_below(cct[p[:n]], ev[p[:n]]) for p in perms])
        rows.append((n, np.nanmedian(v), np.nanpercentile(v, 5), np.nanpercentile(v, 95), np.nanstd(v)))
    return np.array(rows)


def convergence_block():
    """Figure 7.3: spread of the KM estimate over random orderings, for the 37-condition pool
    (aleatoric1) and the 171-condition pool (aleatoric2), both at bus 108, 200 MW."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    specs = [("aleatoric1_ap08_qlim", "C2", [2, 3, 5, 10, 20, 30, 37]),
             ("aleatoric2_ap08_qlim", "C4", [2, 3, 5, 10, 20, 37, 60, 100, 171])]
    for name, color, ns in specs:
        dc, _ = pool_108(name)
        N = len(dc)
        p_full, se_full = greenwood(dc.cct_s, ~dc.cens)
        r = random_orderings(dc, ns)
        print(f"\nconvergence pool {name}: N={N}, full-pool P={p_full:.3f}, "
              f"Greenwood SE {se_full:.3f}, 95% [{p_full - 1.96 * se_full:.3f}, {p_full + 1.96 * se_full:.3f}]")
        for n, med, lo, hi, sd in r:
            theory = np.sqrt(p_full * (1 - p_full) / n * (N - n) / (N - 1))
            print(f"  n={int(n)}: median {med:.3f}, 5-95% [{lo:.3f}, {hi:.3f}] width {hi - lo:.3f}, "
                  f"sd {sd:.4f} (sampling-without-replacement formula {theory:.4f})")
        ax.fill_between(r[:, 0], r[:, 2], r[:, 3], color=color, alpha=0.25)
        ax.plot(r[:, 0], r[:, 1], "o-", color=color, ms=4,
                label=f"{name.split('_')[0]}, N={N}: median and 5-95% over orderings")
        ax.errorbar([N], [p_full], yerr=[[1.96 * se_full], [1.96 * se_full]], fmt="s", color="k", capsize=4,
                    label="Greenwood 95% CI, full pool" if N == 171 else None)
    ax.set_xlabel("number of operating conditions (hour x commitment)")
    ax.set_ylabel(f"KM estimate of P(CCT < {THRESH} s)")
    ax.set_ylim(0, 1)
    ax.set_title("bus 108, 200 MW, fault 66")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "main_convergence_37_vs_171.png", dpi=150)
    plt.close(fig)


def inflation(name):
    R = load(name)
    sited = R[(R.dc_bus.notna()) & (R.status == "ok")]
    n_tot = len(sited)
    n_unc = int((~sited.cens).sum())
    print(f"{name}: sited ok {n_tot}, uncensored {n_unc}, censored {sited.cens.mean():.1%}, "
          f"drop/keep inflation {n_tot / n_unc:.2f}")


if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    aleatoric_block()
    km_37_vs_171_block()
    convergence_block()
    faultscreen_block()
    print()
    inflation("aleatoric1_ap08_qlim")
    inflation("faultscreen1_ap08_qlim")
