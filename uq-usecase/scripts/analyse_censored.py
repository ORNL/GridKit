"""Censored-aware analysis of a campaign, plus the figures.

    python analyse_censored.py --name aleatoric1

CCT is right-censored at `cct_max` (0.5 s): a censored run means the system
survived the longest fault we are willing to consider, so the stored 0.5 is a
LOWER BOUND. Dropping those runs keeps only the stressed half of the operating
space and inflates every risk estimate; this script keeps them.
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402

CAP = 0.5


def km(times, events):
    """Kaplan-Meier S(t) = P(T > t) for right-censored data.

    Returns (t, S) step points. With all censoring at a single point beyond
    every observed event, as here, this reduces to the empirical estimator with
    the censored runs left in the denominator -- which is exactly the part that
    gets dropped by accident.
    """
    t = np.asarray(times, dtype=float)
    e = np.asarray(events, dtype=bool)
    out_t, out_s, S = [0.0], [1.0], 1.0
    for ti in np.unique(t[e]):
        at_risk = (t >= ti).sum()
        d = ((t == ti) & e).sum()
        if at_risk:
            S *= 1.0 - d / at_risk
        out_t.append(ti)
        out_s.append(S)
    return np.array(out_t), np.array(out_s)


def p_below(times, events, thr):
    t, s = km(times, events)
    return 1.0 - (s[t <= thr][-1] if (t <= thr).any() else 1.0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--name", default="aleatoric1")
    args = p.parse_args()

    cfg = cp.CampaignConfig(name=args.name)
    R = pd.read_parquet(cfg.store / "results.parquet")
    R["cens"] = R.cct_censored.fillna(False).astype(bool)
    ok = R[R.status == "ok"].copy()
    ok["evt"] = ~ok.cens  # event = CCT actually observed
    base = ok[ok.dc_bus.isna()]
    sited = ok[~ok.dc_bus.isna()]
    fig_dir = cfg.store / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"{args.name}: {len(R)} runs | ok {len(ok)} | censored {int(ok.cens.sum())} "
        f"({100*ok.cens.mean():.0f}%)"
    )

    # --- 1. the methods point: keeping vs dropping the censored runs ---------
    print("\n=== P(CCT < threshold), censored runs kept vs dropped ===")
    rows = []
    for thr in (0.08, 0.15, 0.20, 0.25, 0.30):
        keep = p_below(sited.cct_s, sited.evt, thr)
        drop = (sited[sited.evt].cct_s < thr).mean()
        rows.append(
            {
                "threshold_s": thr,
                "correct": round(keep, 3),
                "drop_censored": round(drop, 3),
                "inflation_x": round(drop / keep, 2) if keep else np.nan,
            }
        )
    T = pd.DataFrame(rows)
    print(T.to_string(index=False))
    T.to_csv(cfg.store / "p_below_threshold.csv", index=False)

    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(T.threshold_s, T.correct, "o-", label="censored runs kept (correct)")
    ax.plot(T.threshold_s, T.drop_censored, "s--", label="censored runs dropped")
    ax.set_xlabel("CCT threshold (s)")
    ax.set_ylabel("P(CCT < threshold)")
    ax.set_title(f"{args.name}: dropping censored runs inflates risk")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "p_below_threshold.png", dpi=140)
    plt.close(fig)

    # --- 2. survival curves, baseline vs sited ------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for lbl, g in (("no datacenter", base), ("with datacenter", sited)):
        t, s = km(g.cct_s, g.evt)
        ax.step(
            t,
            s,
            where="post",
            label=f"{lbl} (n={len(g)}, " f"{100*g.cens.mean():.0f}% censored)",
        )
    ax.axvline(CAP, color="k", ls=":", lw=1)
    ax.text(CAP, 0.95, " cap", fontsize=8, va="top")
    ax.set_xlabel("CCT (s)")
    ax.set_ylabel("P(CCT > t)")
    ax.set_title(f"{args.name}: Kaplan-Meier, right-censored at {CAP} s")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "km_baseline_vs_sited.png", dpi=140)
    plt.close(fig)

    # --- 3. by net-load quintile -------------------------------------------
    # qcut gives equal-count bins, so equal hours per bin is already
    # proportional allocation: no reweighting of the hour dimension is needed.
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for b in sorted(ok.bin.dropna().unique()):
        g = ok[ok.bin == b]
        t, s = km(g.cct_s, g.evt)
        ax.step(t, s, where="post", label=f"{b} (net {g.net_mw.mean():.0f} MW)")
    ax.axvline(CAP, color="k", ls=":", lw=1)
    ax.set_xlabel("CCT (s)")
    ax.set_ylabel("P(CCT > t)")
    ax.set_title(f"{args.name}: CCT by net-load quintile")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "km_by_netload.png", dpi=140)
    plt.close(fig)
    print("\n=== P(CCT < 0.25) by net-load quintile (censored-aware) ===")
    for b in sorted(ok.bin.dropna().unique()):
        g = ok[ok.bin == b]
        print(
            f"  {b}: net {g.net_mw.mean():7.0f} MW  n={len(g):>3}  "
            f"censored {100*g.cens.mean():3.0f}%  P(CCT<0.25)={p_below(g.cct_s, g.evt, 0.25):.3f}"
        )

    # --- 4. dCCT where both sides are measured ------------------------------
    # A baseline is per (hour, commitment, fault bus): with several
    # contingencies in play, joining on the hour alone fans out into a
    # cartesian product and silently inflates the pair count.
    key = ["hour", "n_commit", "fault_bus"]
    bb = base.set_index(key)[["cct_s", "cens"]].rename(
        columns={"cct_s": "cct_base", "cens": "base_cens"}
    )
    D = sited.join(bb, on=key)
    u = D[(~D.cens) & (~D.base_cens)].copy()
    u["d_cct"] = u.cct_s - u.cct_base
    if len(u):
        fig, ax = plt.subplots(figsize=(6.5, 4))
        for mw, g in u.groupby("dc_mw"):
            if len(g) >= 5:
                ax.scatter(g.net_mw, g.d_cct, s=18, alpha=0.7, label=f"{mw:.0f} MW")
        ax.axhline(0, color="k", lw=1)
        ax.set_xlabel("net load (MW)")
        ax.set_ylabel("dCCT (s), sited - baseline")
        ax.set_title(f"{args.name}: datacenter impact vs operating condition")
        ax.legend(fontsize=8, title="datacenter")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(fig_dir / "dcct_vs_netload.png", dpi=140)
        plt.close(fig)
        print(f"\n=== dCCT (uncensored pairs only, n={len(u)} of {len(D)}) ===")
        print(
            u.groupby("dc_mw").d_cct.agg(["count", "mean", "min"]).round(4).to_string()
        )
        print(
            "  NOTE: this subset is biased toward stressed conditions, because a pair"
        )
        print(
            "  survives only when BOTH sides are below the cap. Treat as conditional."
        )

    # --- 5. partial identification: use the one-sided pairs too -------------
    # A censored baseline is not uninformative. If baseline >= 0.5 and the sited
    # run measured 0.26, then damage >= 0.24 -- a bound in the safety-relevant
    # direction. Taking the baseline as exactly the cap gives a LOWER BOUND on
    # damage that is exact when both sides are measured.
    P = D[~D.cens].copy()  # sited measured; baseline may or may not be
    if len(P):
        P["dmg_lower"] = np.where(P.base_cens, CAP, P.cct_base) - P.cct_s
        exact = int((~P.base_cens).sum())
        bound = int(P.base_cens.sum())
        both_cens = int((D.cens & D.base_cens).sum())
        denom = len(P) + both_cens
        print(
            f"\n=== damage bounds: {len(P)} usable pairs "
            f"({exact} exact + {bound} one-sided) vs {exact} before ==="
        )
        print("  P(damage > x) is a LOWER bound: both-censored pairs carry no")
        print("  information, so they stay in the denominator only.")
        rows = [
            {
                "damage_gt_s": x,
                "n_pairs": int((P.dmg_lower > x).sum()),
                "P_lower": round(float((P.dmg_lower > x).sum()) / denom, 3),
            }
            for x in (0.0, 0.05, 0.10, 0.20, 0.30)
        ]
        B = pd.DataFrame(rows)
        print(B.to_string(index=False))
        B.to_csv(cfg.store / "damage_bounds.csv", index=False)
        print(
            f"  denominator {denom} = {len(P)} informative + {both_cens} both-censored"
        )
        print(f"  worst identified damage: {P.dmg_lower.max():.4f} s")

        fig, ax = plt.subplots(figsize=(6.5, 4))
        for lbl, g in (
            ("both measured (exact)", P[~P.base_cens]),
            ("baseline censored (lower bound)", P[P.base_cens]),
        ):
            if len(g):
                d = np.sort(g.dmg_lower.values)
                ax.step(
                    d,
                    np.arange(1, len(d) + 1) / denom,
                    where="post",
                    label=f"{lbl}, n={len(g)}",
                )
        d = np.sort(P.dmg_lower.values)
        ax.step(
            d,
            np.arange(1, len(d) + 1) / denom,
            where="post",
            color="k",
            lw=2,
            label=f"combined, n={len(P)}",
        )
        ax.set_xlabel("damage = baseline CCT - sited CCT (s)")
        ax.set_ylabel("lower bound on cumulative probability")
        ax.set_title(f"{args.name}: datacenter damage, partial identification")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(fig_dir / "damage_bounds.png", dpi=140)
        plt.close(fig)

    print(f"\nfigures -> {fig_dir}")


if __name__ == "__main__":
    main()
