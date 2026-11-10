"""Rebuild the alphaP frontier table from the per-condition records.

    python alphap_collect.py --name afront1

Reads store/records/af_*.json, which are the source of truth: one file per
(operating point, alphaP), written atomically. Answers the one question the
sweep exists for -- is the frontier a property of the bus, or of the operating
condition the dynamics initialise from?
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

STORE = Path("/kfs2/projects/scidac/scidac-data/gridkit-campaigns")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--name", default="afront1")
    args = p.parse_args()

    rec_dir = STORE / args.name / "records"
    rows = [json.loads(f.read_text()) for f in sorted(rec_dir.glob("af_*.json"))]
    if not rows:
        sys.exit(f"no records under {rec_dir}")
    df = pd.DataFrame(rows)
    print(f"{len(df)} records from {rec_dir}")

    sims = df[df.alpha_p_dc.notna()].copy()
    nodisp = df[df.alpha_p_dc.isna()]
    if len(nodisp):
        print(f"{len(nodisp)} dispatch-infeasible points (excluded from frontier)")

    sims["ok"] = sims.status.eq("ok")
    # fault_bus MUST be in the key: without it the max() below silently takes
    # the best fault location and every frontier reads 1.0.
    key = ["dc_bus", "dc_mw", "hour", "n_commit", "fault_bus"]

    # n_commit differs per (hour, MW), so pivoting on it is 78% structurally
    # empty. Rank it within each (hour, MW) instead: lo / hi commitment.
    lvl = (
        sims[["hour", "dc_mw", "n_commit"]]
        .drop_duplicates()
        .sort_values(["hour", "dc_mw", "n_commit"])
    )
    lvl["commit"] = (
        lvl.groupby(["hour", "dc_mw"])
        .n_commit.rank("dense")
        .map({1.0: "lo", 2.0: "hi"})
    )
    sims = sims.merge(lvl, on=["hour", "dc_mw", "n_commit"], how="left")

    # frontier = highest alphaP that ran; NaN if none did
    fr = sims[sims.ok].groupby(key)["alpha_p_dc"].max().rename("frontier").reset_index()
    allpts = sims[key + ["commit"]].drop_duplicates()
    fr = allpts.merge(fr, on=key, how="left")

    print("\n=== frontier (highest alphaP that simulated) ===")
    print("NaN = no alphaP ran at all.  '-' columns do not exist by design.")
    pv = fr.pivot_table(
        index=["dc_bus", "dc_mw"],
        columns=["hour", "commit"],
        values="frontier",
        dropna=False,
    )
    print(pv.to_string())

    spread = (pv.max(axis=1) - pv.min(axis=1)).rename("spread")
    print("\n=== does the operating condition move the frontier? ===")
    print(spread.to_string())
    moved = spread[spread > 0]
    print(
        f"\n{len(moved)} of {len(spread)} (bus, MW) cells have a frontier that "
        f"moves with the operating condition."
    )
    if len(moved) == 0:
        print(
            "VERDICT: frontier is a property of BUS and SIZE only. The AC "
            "solution and commitment do not move it."
        )
    else:
        print(
            "VERDICT: the operating condition DOES move the frontier, so it is "
            "not a pure bus property. Cells affected:"
        )
        print(moved.to_string())

    print("\n=== failure modes ===")
    bad = sims[~sims.ok]
    if len(bad):
        print(bad.groupby(["fail_mode"], dropna=False).size().to_string())
        print("\nfailure mode by alphaP:")
        print(pd.crosstab(bad.alpha_p_dc, bad.fail_mode.fillna("?")).to_string())
    else:
        print("none")

    if sims.fault_bus.nunique() > 1:
        fault_location(sims)

    out = STORE / args.name / "frontier.csv"
    fr.to_csv(out, index=False)
    sims.to_csv(STORE / args.name / "frontier_runs.csv", index=False)
    print(f"\nwrote {out}")


def fault_location(sims):
    """Does moving the fault change the frontier?

    IDACalcIC failures happen before the fault exists, so a dependence here
    would mean the BusFault device participates in the initial-condition
    solve -- which is the point of the experiment.
    """
    print("\n=== does FAULT LOCATION move the frontier? ===")
    key = ["dc_bus", "dc_mw", "hour", "n_commit"]
    f = (
        sims[sims.ok]
        .groupby(key + ["fault_bus"])["alpha_p_dc"]
        .max()
        .rename("frontier")
        .reset_index()
    )
    allp = sims[key + ["fault_bus"]].drop_duplicates()
    f = allp.merge(f, on=key + ["fault_bus"], how="left")
    pv = f.pivot_table(index=key, columns="fault_bus", values="frontier", dropna=False)
    spread = pv.max(axis=1) - pv.min(axis=1)
    moved = spread[spread > 0]
    print(f"{len(pv)} (point) cells x {pv.shape[1]} fault locations")
    print(f"cells whose frontier changes with fault location: {len(moved)}")
    if len(moved) == 0:
        print(
            "VERDICT: fault location does NOT move the frontier. Consistent with\n"
            "         the failure being an initial-condition problem that the\n"
            "         BusFault device plays no part in."
        )
    else:
        print("VERDICT: fault location DOES move the frontier. Affected cells:")
        print(pv.loc[moved.index].to_string())

    print("\nfrontier by DC bus and size, averaged over conditions:")
    print(
        f.pivot_table(
            index=["dc_bus", "dc_mw"], columns="fault_bus", values="frontier"
        ).to_string()
    )

    bad = sims[~sims.ok]
    if len(bad):
        print("\nfailure count by fault location (should be flat if pre-fault):")
        print(pd.crosstab(bad.fault_bus, bad.fail_mode.fillna("?")).to_string())


if __name__ == "__main__":
    main()
