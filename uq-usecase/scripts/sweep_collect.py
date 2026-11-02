"""Gather per-condition records into one table and report what is missing.

    python sweep_collect.py --name smoke

Reads `store/records/*.json`, which is the only shared state the workers
produce, and writes `store/results.parquet` plus a csv for eyeballing.
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--name", required=True)
    args = p.parse_args()

    cfg = cp.CampaignConfig(name=args.name)
    files = sorted((cfg.store / "records").glob("*.json"))
    if not files:
        print(f"no records under {cfg.store/'records'}")
        return

    R = pd.DataFrame([json.loads(f.read_text()) for f in files])
    R.to_parquet(cfg.store / "results.parquet", compression="zstd", index=False)
    R.to_csv(cfg.store / "results.csv", index=False)
    print(f"{len(R)} records -> {cfg.store/'results.parquet'}")

    cond = cfg.store / "conditions.csv"
    if cond.exists():
        C = pd.read_csv(cond)
        missing = set(C.tag) - set(R.tag)
        print(f"planned {len(C)}, collected {len(R)}, MISSING {len(missing)}")
        for t in sorted(missing)[:10]:
            print("   missing:", t)
    print("\nstatus:", R.status.value_counts().to_dict())

    if "fail_reason" in R:
        fr = R.fail_reason.dropna()
        if len(fr):
            print("fail reasons:", fr.value_counts().to_dict())

    ok = R[R.status == "ok"]
    if len(ok) and "cct_s" in ok:
        c = ok.cct_s.dropna()
        print(
            f"\ncct: n={len(c)} min={c.min():.4f} median={c.median():.4f} "
            f"max={c.max():.4f}  censored={int(ok.cct_censored.sum())}"
        )
    if len(ok) and "wall_s" in ok:
        print(
            f"wall per condition: median {ok.wall_s.median():.1f}s "
            f"total {R.wall_s.sum():.0f}s"
        )
    # which host ran what is not recorded; node coverage is checked via the
    # slurm chunk logs instead.


if __name__ == "__main__":
    main()
