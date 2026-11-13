"""Plot the fault-on nose trace from ida_nose_trace.py.

Reads /scratch/isatkaus/gridkit-campaigns/ida_zoom/nose_trace/trace.csv and writes
work-notes/figures/ida_nose_alpha08.png: |V| at the faulted bus vs fault strength s.
Left: full range 0 to 1. Right: zoom on the fold of bus 148.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

TRACE = Path("/scratch/isatkaus/gridkit-campaigns/ida_zoom/nose_trace/trace.csv")
FIG = Path(__file__).resolve().parents[1] / "work-notes" / "figures" / "ida_nose_alpha08.png"
STYLES = ((148, "o-", "bus 148 (fold)"), (3, "s-", "bus 3 (no fold)"))


def draw(ax, df, zoom=False):
    for bus, style, label in STYLES:
        d = df[df["bus"] == bus]
        ax.plot(d["s"], d["v_fault"], style, ms=3.5, lw=1.2, label=label)
    fold = df[df["bus"] == 148]
    s_stop = float(fold["s_stop"].iloc[0])
    ax.axvline(s_stop, color="C0", ls=":", lw=1)
    ax.axvline(1.0, color="k", ls="--", lw=1)
    if zoom:
        ax.set_xlim(0.30, 0.40)
        ax.set_ylim(0.15, 0.30)
        ax.set_title("zoom on the fold (bus 148)", fontsize=9)
    else:
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.05)
        ax.set_title("full range, upper branch", fontsize=9)
        ax.text(0.97, 0.9, "target s = 1\n(X = 0.01 pu)", fontsize=8, ha="right", va="top",
                transform=ax.transAxes)
        ax.annotate(f"fold: s_stop = {s_stop:.4f}", xy=(s_stop, 0.18), xytext=(0.45, 0.55),
                    textcoords="axes fraction", arrowprops=dict(arrowstyle="->", lw=0.8),
                    fontsize=8)
    ax.set_xlabel("fault strength s (fault shunt multiplier)")
    ax.set_ylabel("|V| at faulted bus (pu)")
    ax.grid(alpha=0.3)


def main():
    df = pd.read_csv(TRACE)
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4.2))
    draw(ax0, df)
    draw(ax1, df, zoom=True)
    ax1.annotate(f"X_crit = 0.01/s_stop = {0.01 / float(df[df['bus'] == 148]['s_stop'].iloc[0]):.4f} pu",
                 xy=(0.3672, 0.181), xytext=(0.315, 0.245), fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.8))
    ax0.legend(loc="lower left", fontsize=8)
    fig.suptitle("Upper branch of the fixed-state solution, alphaP 0.8, 200 MW, dispatch h1879",
                 fontsize=10)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=150)
    print(f"wrote {FIG}")


if __name__ == "__main__":
    main()
