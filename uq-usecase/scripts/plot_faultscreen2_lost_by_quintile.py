"""Figure 7.6: share of faultscreen-2 runs losing synchronism before the 0.5 s cap,
by net-load quintile, for each fault bus.

Left panel: no datacenter. Right panel: 100 MW datacenter at bus 108.
Denominator per point: runs with status 'ok' and a CCT (sim_failed and solve_failed
runs are excluded). Quintiles are the `bin` column (q0 lowest to q4 highest net load).
No title is drawn: the report caption carries it.

python plot_faultscreen2_lost_by_quintile.py [--out <png>]
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path("/kfs2/projects/scidac/scidac-data/gridkit-campaigns")
CAMPAIGN = "faultscreen2_ap08_qlim"
DEFAULT_OUT = Path("/home/isatkaus/gridkit/uq-usecase/work-notes/figures/main_faultscreen2_lost_by_quintile.png")
BINS = ["q0", "q1", "q2", "q3", "q4"]
NAMED = {66: "C0", 68: "C1", 94: "C2", 181: "C4"}  # colours as in the report figure
DC_BUS = 108
DC_MW = 100


def shares(ok, fault_buses):
    """Loss share per (fault bus, quintile): lost / runs with CCT."""
    S = ok[ok.fault_bus.isin(fault_buses)]
    g = S.groupby(["fault_bus", "bin"])["lost"].agg(["sum", "size"])
    return (g["sum"] / g["size"]).unstack("bin").reindex(columns=BINS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    R = pd.read_parquet(ROOT / CAMPAIGN / "results.parquet")
    ok = R[(R.status == "ok") & R.cct_s.notna()].copy()
    ok["lost"] = ~ok.cct_censored.astype(bool)

    base = ok[ok.dc_mw.fillna(0) == 0]
    dc = ok[(ok.dc_mw == DC_MW) & (ok.dc_bus == DC_BUS)]

    all_buses = sorted(R.fault_bus.dropna().astype(int).unique())
    others_base = [b for b in all_buses if b not in NAMED]
    panels = [("No datacenter", base, others_base), (f"{DC_MW} MW datacenter at bus {DC_BUS}", dc, None)]

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.7), dpi=100, sharey=True)
    for ax, (title, data, others) in zip(axes, panels):
        present = sorted(data.fault_bus.dropna().astype(int).unique())
        if others is None:
            others = [b for b in present if b not in NAMED]
        table = shares(data, present)
        x = range(len(BINS))
        n_others = 0
        for b in others:
            if b not in table.index:
                continue
            n_others += 1
            ax.plot(x, table.loc[b].values, color="#bbbbbb", lw=1.5)
        for b, col in NAMED.items():
            if b in table.index:
                ax.plot(x, table.loc[b].values, color=col, lw=2, marker="o", label=f"bus {b}")
        if n_others:
            ax.plot([], [], color="#bbbbbb", lw=1.5, label=f"{n_others} other fault buses")
        ax.set_xticks(list(x))
        ax.set_xticklabels([f"{q}\nlowest" if q == "q0" else f"{q}\nhighest" if q == "q4" else q for q in BINS])
        ax.set_xlabel("Net-load quintile")
        ax.set_title(title)
        ax.set_ylim(-0.05, 1.05)
        ax.legend(loc="center right", fontsize=9)
        print(f"{title}: {n_others} other fault buses with CCT data; "
              f"fault buses present: {present}")

    axes[0].set_ylabel("Fraction of runs losing synchronism\nbefore the 0.5 s cap")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=100)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
