"""Enumerate a sweep design and write the sbatch scripts for it. No simulation.

    python sweep_build.py --name smoke --nodes 1 --workers 4 --hours 2 --submit

Produces, under the campaign store:

    conditions.csv   one row per condition, the unit of work
    meta.json        provenance: git sha, design, seed, config
    slurm/chunk_NN.sh

A condition is one `OperatingPoint` crossed with one `Disturbance`. The two are
kept apart because stage A costs a Julia solve and stage B/C cost ~0.4 s, so
the worker sorts by operating point and solves each one once. See
work-notes/datacenter_load_uq_plan.md, "Parallel sampling runs".
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402
import slurm_utils as su  # noqa: E402

REPO = Path("/home/isatkaus/gridkit/uq-usecase")
PY = "/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python"
DATA = "/kfs2/projects/scidac/scidac-data/ACTIVSg200/grid-v1"


def build_conditions(
    cfg,
    hours_per_bin,
    net_bins,
    commit_fracs,
    fault_buses,
    durations,
    h_scales,
    seed,
    dc_buses=(),
    dc_mws=(),
):
    dem = pd.read_csv(f"{DATA}/demand.csv")
    wnd = pd.read_csv(f"{DATA}/wind.csv")
    demc = [c for c in dem.columns if c != "DateTime"]
    wndc = [c for c in wnd.columns if c != "DateTime"]
    net_all = dem[demc].sum(axis=1) - wnd[wndc].sum(axis=1)
    bins = pd.qcut(net_all, net_bins, labels=[f"q{i}" for i in range(net_bins)])

    # (None, 0.0) is the no-datacenter baseline every sited case is compared to,
    # so it is always present rather than being a separate run.
    siting = [(None, 0.0)] + [(int(b), float(m)) for b in dc_buses for m in dc_mws]

    rng = np.random.default_rng(seed)
    order = cp.merit_order(cfg)
    rows = []
    for b in [f"q{i}" for i in range(net_bins)]:
        pool = net_all.index[bins == b].to_numpy()
        take = min(hours_per_bin, len(pool))
        for h in rng.choice(pool, size=take, replace=False):
            h = int(h)
            lo, hi = cp.commitment_band(cfg, float(net_all[h]), order)
            if lo is None:
                continue
            for f in commit_fracs:
                n = int(round(lo + f * (hi - lo)))
                for dcb, dcm in siting:
                    op = cp.OperatingPoint(hour=h, n_commit=n, dc_bus=dcb, dc_mw=dcm)
                    for fb in fault_buses:
                        for dur in durations:
                            for hs in h_scales:
                                d = cp.Disturbance(
                                    fault_bus=fb, duration=dur, h_scale=hs
                                )
                                rows.append(
                                    {
                                        "tag": f"{op.tag}__{d.tag}",
                                        "bin": b,
                                        "hour": h,
                                        "n_commit": n,
                                        "dc_bus": op.dc_bus,
                                        "dc_mw": op.dc_mw,
                                        "fault_bus": fb,
                                        "duration": dur,
                                        "h_scale": hs,
                                    }
                                )
    C = pd.DataFrame(rows).drop_duplicates(subset="tag")
    # Sorted by operating point so a slice can be cut on point boundaries:
    # splitting one point across workers would re-solve stage A for each.
    C = C.sort_values(
        ["hour", "n_commit", "dc_bus", "dc_mw"], na_position="first"
    ).reset_index(drop=True)
    C.insert(0, "idx", C.index)
    return C


def slice_on_points(C, n_chunks):
    """Contiguous (start, count) ranges that never split an operating point."""
    sizes = (
        C.groupby(["hour", "n_commit", "dc_bus", "dc_mw"], dropna=False, sort=False)
        .size()
        .tolist()
    )
    target = len(C) / max(n_chunks, 1)
    chunks, start, cur = [], 0, 0
    for s in sizes:
        cur += s
        # close the chunk once it reaches its share, unless chunks remain scarce
        if cur >= target and len(chunks) < n_chunks - 1:
            chunks.append((start, cur))
            start += cur
            cur = 0
    if cur:
        chunks.append((start, cur))
    return chunks


def sbatch_script(node_i, sub, args, store, slurm_dir):
    """One node. Each worker gets a contiguous slice so it can reuse a solve."""
    launches = "\n".join(
        f"{PY} {REPO}/scripts/sweep_worker.py --conditions {store}/conditions.csv "
        f"--name {args.name} --start {s} --count {c} "
        f"> {slurm_dir}/w{node_i:02d}_{w:02d}.log 2>&1 &"
        for w, (s, c) in enumerate(sub)
    )
    qos = f"#SBATCH --qos={args.qos}" if args.qos else ""
    return f"""#!/bin/bash
#SBATCH --job-name=sweep_{node_i:02d}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={len(sub)}
#SBATCH --time={args.walltime}
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --output={slurm_dir}/chunk_{node_i:02d}.log
{qos}

set -uo pipefail
module load julia/1.12.1 2>/dev/null || true
export PYTHONPATH={REPO}/py-utils:${{PYTHONPATH:-}}

echo "node $(hostname)  workers {len(sub)}  $(date)"
{launches}
wait
echo "done $(date)"
"""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--name", required=True)
    p.add_argument("--nodes", type=int, default=1)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--hours", type=int, default=2, help="hours per net-load bin")
    p.add_argument("--net-bins", type=int, default=2)
    p.add_argument("--commit-fracs", default="0.5")
    p.add_argument("--fault-buses", default="66")
    p.add_argument("--durations", default="0.1")
    p.add_argument("--h-scales", default="1.0")
    p.add_argument("--dc-buses", default="", help="datacenter candidate buses")
    p.add_argument("--dc-mw", default="", help="datacenter sizes in MW")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--walltime", default="01:00:00")
    p.add_argument("--account", default="scidac")
    p.add_argument("--partition", default="short")
    p.add_argument("--qos", default=None)
    p.add_argument("--submit", action="store_true")
    args = p.parse_args()

    flist = lambda s: [float(x) for x in s.split(",")]  # noqa: E731
    ilist = lambda s: [int(x) for x in s.split(",")]  # noqa: E731

    cfg = cp.CampaignConfig(name=args.name)
    store, slurm_dir = cfg.store, cfg.store / "slurm"
    slurm_dir.mkdir(parents=True, exist_ok=True)

    C = build_conditions(
        cfg,
        args.hours,
        args.net_bins,
        flist(args.commit_fracs),
        ilist(args.fault_buses),
        flist(args.durations),
        flist(args.h_scales),
        args.seed,
        dc_buses=ilist(args.dc_buses) if args.dc_buses else (),
        dc_mws=flist(args.dc_mw) if args.dc_mw else (),
    )
    C.to_csv(store / "conditions.csv", index=False)
    print(
        f"{len(C)} conditions over {C.hour.nunique()} hours, "
        f"{C.groupby(['hour','n_commit','dc_bus','dc_mw'], dropna=False).ngroups} "
        f"distinct operating points"
    )

    sha = subprocess.run(
        ["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True
    ).stdout.strip()
    (store / "meta.json").write_text(
        json.dumps(
            {
                "name": args.name,
                "git_sha": sha,
                "seed": args.seed,
                "n_conditions": len(C),
                "nodes": args.nodes,
                "workers": args.workers,
                "design": {
                    k: getattr(args, k)
                    for k in (
                        "hours",
                        "net_bins",
                        "commit_fracs",
                        "fault_buses",
                        "durations",
                        "h_scales",
                        "dc_buses",
                        "dc_mw",
                    )
                },
                "scratch": str(cfg.scratch),
                "store": str(store),
            },
            indent=1,
        )
    )

    node_chunks = slice_on_points(C, args.nodes)
    paths = []
    for i, (lo, n) in enumerate(node_chunks):
        sub = slice_on_points(C.iloc[lo : lo + n], args.workers)
        sub = [(lo + s, c) for s, c in sub]
        sp = slurm_dir / f"chunk_{i:02d}.sh"
        sp.write_text(sbatch_script(i, sub, args, store, slurm_dir))
        sp.chmod(0o755)
        paths.append(str(sp))
        print(
            f"  chunk_{i:02d}: conditions {lo}-{lo+n-1} over {len(sub)} workers "
            f"{[(s, c) for s, c in sub]}"
        )

    if args.submit:
        ids = su.submit_jobs_batch(paths)
        (
            su.save_job_metadata(
                str(store / "job_metadata.json"), ids, len(C), args.nodes
            )
            if hasattr(su, "save_job_metadata")
            else None
        )
        print("submitted:", ids)
    else:
        print(
            f"\nnot submitted. to run:\n  for f in {slurm_dir}/chunk_*.sh; do sbatch $f; done"
        )


if __name__ == "__main__":
    main()
