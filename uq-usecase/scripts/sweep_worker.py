"""Run a contiguous slice of conditions. One record file per condition.

    python sweep_worker.py --conditions C.csv --name smoke --start 0 --count 8

Sorts its slice by operating point so one dispatch solve serves every
disturbance on that point: stage A is ~29 s, stage B/C ~0.4 s. Never appends to
a shared file; a condition is done if its record exists, which is also how
resume works.
"""

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402

DATA = "/kfs2/projects/scidac/scidac-data/ACTIVSg200/grid-v1"


def profiles(cfg):
    dem = pd.read_csv(f"{DATA}/demand.csv")
    wnd = pd.read_csv(f"{DATA}/wind.csv")
    demc = [c for c in dem.columns if c != "DateTime"]
    wndc = [c for c in wnd.columns if c != "DateTime"]
    wg = {c: cfg.gid[int(c.replace("gen", "")) - 1] for c in wndc}

    def at(h):
        return (
            {int(x.replace("bus", "")): float(dem.loc[h, x]) for x in demc},
            {wg[x]: float(wnd.loc[h, x]) for x in wndc},
        )

    return at


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--conditions", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--count", type=int, required=True)
    p.add_argument("--no-cct", action="store_true")
    # datacenter constant-power share; the campaign default is 0.8
    p.add_argument("--dc-alpha-p", type=float, default=0.8)
    # PV->PQ switching of generator Q limits in the AC PF. Default: on, unless
    # --ac-opf. --no-enforce-q-lims reproduces the unlimited-Q campaigns.
    qg = p.add_mutually_exclusive_group()
    qg.add_argument("--enforce-q-lims", dest="enforce_q_lims", action="store_true")
    qg.add_argument("--no-enforce-q-lims", dest="enforce_q_lims", action="store_false")
    p.set_defaults(enforce_q_lims=None)
    p.add_argument("--ac-opf", action="store_true")
    # control for --ac-opf: generator Q limits widened so they never bind
    p.add_argument("--ac-opf-relax-q", action="store_true")
    args = p.parse_args()

    cfg = cp.CampaignConfig(
        name=args.name,
        dc_alpha_p=args.dc_alpha_p,
        enforce_q_lims=args.enforce_q_lims,
        ac_opf=args.ac_opf,
        ac_opf_relax_q=args.ac_opf_relax_q,
    )
    missing = cp.check_runner(cfg)
    if missing:
        print(
            f"FATAL: DynamicSimulation cannot load its libraries on {__import__('socket').gethostname()}:",
            flush=True,
        )
        for m in missing:
            print("   ", m, flush=True)
        sys.exit(2)
    at = profiles(cfg)
    C = pd.read_csv(args.conditions)
    mine = C.iloc[args.start : args.start + args.count]
    # group by operating point: the expensive stage is solved once per group
    mine = mine.sort_values(["hour", "n_commit", "dc_bus", "dc_mw"])
    print(
        f"worker start={args.start} count={len(mine)} "
        f"points={mine.groupby(['hour','n_commit','dc_bus','dc_mw'], dropna=False).ngroups}",
        flush=True,
    )

    t_all = time.time()
    done = skipped = failed = 0

    # Collect every operating point this worker owns, then solve them all in
    # one batched Julia pass: launch cost (~10 s) dwarfs a warm solve (~0.1 s).
    groups, specs = [], []
    for key, grp in mine.groupby(
        ["hour", "n_commit", "dc_bus", "dc_mw"], dropna=False, sort=False
    ):
        hour, n_commit, dc_bus, dc_mw = key
        todo = [
            r
            for r in grp.itertuples()
            if not (cfg.store / "records" / f"{r.tag}.json").exists()
        ]
        skipped += len(grp) - len(todo)
        if not todo:
            continue
        demand, wind = at(int(hour))
        groups.append((key, todo))
        # Name by the FULL operating point: siting variants share an hour and a
        # commitment, so omitting dc_bus/dc_mw makes them overwrite each other's
        # dcin.m/dc.m/ac.m and read the wrong dispatch.
        op = cp.OperatingPoint(
            hour=int(hour),
            n_commit=int(n_commit),
            dc_bus=None if pd.isna(dc_bus) else int(dc_bus),
            dc_mw=0.0 if pd.isna(dc_mw) else float(dc_mw),
        )
        specs.append(
            {
                "demand": demand,
                "wind": wind,
                "n_commit": int(n_commit),
                "workdir": cfg.scratch / op.tag,
                "dc_bus": op.dc_bus,
                "dc_mw": op.dc_mw,
            }
        )

    if not specs:
        print(f"worker done: 0 run, {skipped} already present, 0.0s", flush=True)
        return

    t0 = time.time()
    solved = cp.solve_points_batch(cfg, specs)
    t_solve = time.time() - t0
    print(
        f"  batched stage A: {len(specs)} points in {t_solve:.1f}s "
        f"({t_solve/len(specs):.1f}s each)",
        flush=True,
    )

    for (key, todo), (case, base), spec in zip(groups, solved, specs):
        hour = key[0]
        point_dir = spec["workdir"]
        if case is None:
            # Dispatch infeasible: that is a result for every disturbance here.
            for r in todo:
                rec = dict(
                    base, tag=r.tag, bin=r.bin, hour=int(hour), wall_s=round(t_solve, 1)
                )
                cp.persist(cfg, r.tag, rec)
                failed += 1
            print(
                f"  {key} SOLVE FAILED {base.get('fail_reason')} "
                f"({len(todo)} conditions)",
                flush=True,
            )
            continue

        for r in todo:
            dist = cp.Disturbance(
                fault_bus=int(r.fault_bus),
                duration=float(r.duration),
                h_scale=float(r.h_scale),
            )
            t1 = time.time()
            rec = dict(base)
            rec.update(
                cp.run_disturbance(
                    cfg, case, dist, point_dir / dist.tag, with_cct=not args.no_cct
                )
            )
            rec.update(
                tag=r.tag, bin=r.bin, hour=int(hour), wall_s=round(time.time() - t1, 1)
            )
            cp.persist(cfg, r.tag, rec, rec.pop("_csv", None))
            done += 1
            print(
                f"  {r.tag:<28} {rec['status']:<12} "
                f"cct={rec.get('cct_s')} {rec['wall_s']}s",
                flush=True,
            )

    print(
        f"worker done: {done} run, {skipped} already present, "
        f"{failed} solve-failed, {time.time() - t_all:.1f}s total",
        flush=True,
    )


if __name__ == "__main__":
    main()
