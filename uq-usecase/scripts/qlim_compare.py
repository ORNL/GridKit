"""Compare CCT between unlimited and Q-limited campaigns, matched by tag.

usage: python qlim_compare.py <unlimited-campaign> <qlim-campaign>
Reads records/*.json from the kfs2 campaign store. Prints counts, censoring,
and paired CCT differences for conditions ok in both campaigns, split into
datacenter (dc_mw > 0) and baseline (dc_mw == 0) rows.
"""

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/kfs2/projects/scidac/scidac-data/gridkit-campaigns")


def load(name):
    rows = [json.load(open(p)) for p in sorted((ROOT / name / "records").glob("*.json"))]
    df = pd.DataFrame(rows)
    df["campaign"] = name
    return df


def main(unl, qlim):
    a, b = load(unl), load(qlim)
    print(f"{unl}: {len(a)} records; {qlim}: {len(b)} records")
    for name, df in ((unl, a), (qlim, b)):
        print(f"  {name}: status counts {df['status'].value_counts().to_dict()}")

    ok_a = a[a["status"] == "ok"].set_index("tag")
    ok_b = b[b["status"] == "ok"].set_index("tag")
    both = ok_a.index.intersection(ok_b.index)
    only_a = ok_a.index.difference(ok_b.index)
    only_b = ok_b.index.difference(ok_a.index)
    print(f"ok in both: {len(both)}; ok only unlimited: {len(only_a)}; ok only qlim: {len(only_b)}")

    pair = pd.DataFrame({
        "cct_unl": ok_a.loc[both, "cct_s"],
        "cct_qlim": ok_b.loc[both, "cct_s"],
        "cens_unl": ok_a.loc[both, "cct_censored"],
        "cens_qlim": ok_b.loc[both, "cct_censored"],
        "dc_mw": ok_b.loc[both, "dc_mw"],
        "bin": ok_b.loc[both, "bin"],
    })
    pair[["cens_unl", "cens_qlim"]] = pair[["cens_unl", "cens_qlim"]].astype(bool)
    pair["d_cct"] = pair["cct_qlim"] - pair["cct_unl"]

    for label, sub in (
        ("datacenter rows (dc_mw > 0)", pair[pair["dc_mw"] > 0]),
        ("baseline rows (dc_mw == 0)", pair[pair["dc_mw"] == 0]),
    ):
        print(f"\n== {label}: n = {len(sub)}")
        if sub.empty:
            continue
        print(f"  CCT unlimited: median {sub['cct_unl'].median():.4f} s, "
              f"censored {int(sub['cens_unl'].sum())}")
        print(f"  CCT qlim:      median {sub['cct_qlim'].median():.4f} s, "
              f"censored {int(sub['cens_qlim'].sum())}")
        print(f"  d_cct (qlim - unl): mean {sub['d_cct'].mean():+.4f} s, "
              f"min {sub['d_cct'].min():+.4f}, max {sub['d_cct'].max():+.4f}")
        print(f"  changed (|d| > 1e-6 s): {int((sub['d_cct'].abs() > 1e-6).sum())} of {len(sub)}")
        print(f"  qlim lower than unlimited: {int((sub['d_cct'] < -1e-6).sum())}; "
              f"higher: {int((sub['d_cct'] > 1e-6).sum())}")

    unc = pair[(~pair["cens_unl"]) & (~pair["cens_qlim"])]
    print(f"\n== uncensored in both (bisected, not at cap): n = {len(unc)}")
    if not unc.empty:
        print(f"  d_cct (qlim - unl): mean {unc['d_cct'].mean():+.4f} s, "
              f"min {unc['d_cct'].min():+.4f}, max {unc['d_cct'].max():+.4f}")
        print(f"  changed (|d| > 1e-6 s): {int((unc['d_cct'].abs() > 1e-6).sum())} of {len(unc)}")

    dc = pair[pair["dc_mw"] > 0]
    if not dc.empty:
        print("\n== datacenter rows, by bin (median CCT in s)")
        print(dc.groupby("bin")[["cct_unl", "cct_qlim", "d_cct"]].median().round(4).to_string())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
