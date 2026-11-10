"""CCT at two alphaP values over the alphap_frontier design.

Same dispatches as alphap_frontier.py (same buses, hours, commitment fractions
and MW levels, same defaults), but every run carries CCT (`with_cct=True`) and
only the alphaP values given by --alphas are run. Use --start/--count to split
the design across workers, as afront3.sh does.

    python cct_check.py --name cct181 --fault-buses 181 --start 0 --count 17
"""

import argparse
import copy
import sys
import time

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")

import alphap_frontier as af  # noqa: E402
import campaign as cp  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--name", required=True)
    p.add_argument("--fault-buses", default="94,68")
    p.add_argument("--alphas", default="0.2,1.0")
    p.add_argument("--dc-mw", default="100,200")
    p.add_argument("--commit-fracs", default="0.25,0.75")
    p.add_argument("--n-buses", type=int, default=16)
    p.add_argument("--duration", type=float, default=0.1)
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--count", type=int, default=10**9)
    p.add_argument("--list", action="store_true")
    args = p.parse_args()

    alphas = [float(x) for x in args.alphas.split(",")]
    fbs = [int(x) for x in args.fault_buses.split(",")]
    cfg = cp.CampaignConfig(name=args.name)
    missing = cp.check_runner(cfg)
    if missing:
        print("FATAL: runner cannot load its libraries:", *missing, sep="\n  ")
        sys.exit(2)

    at = af.profiles(cfg)
    cand = cp.candidate_buses(cfg, n=args.n_buses, require_load=True, keep=af.KEEP)
    buses = sorted(set(cand.bus.astype(int)) | {34})
    fracs = [float(x) for x in args.commit_fracs.split(",")]
    mws = [float(x) for x in args.dc_mw.split(",")]
    D = af.build_design(cfg, at, buses, af.HOURS, fracs, mws)
    mine = D.iloc[args.start : args.start + args.count]
    print(
        f"design: {len(D)} dispatches, slice {args.start}..{args.start+len(mine)-1}, "
        f"x faults {fbs} x alphaP {alphas}",
        flush=True,
    )
    if args.list:
        print(mine.to_string(index=False))
        return

    specs = []
    for r in mine.itertuples():
        demand, wind = at(int(r.hour))
        op = cp.OperatingPoint(
            hour=int(r.hour),
            n_commit=int(r.n_commit),
            dc_bus=int(r.dc_bus),
            dc_mw=float(r.dc_mw),
        )
        specs.append(
            {
                "demand": demand,
                "wind": wind,
                "n_commit": int(r.n_commit),
                "workdir": cfg.scratch / op.tag,
                "dc_bus": op.dc_bus,
                "dc_mw": op.dc_mw,
            }
        )

    t0 = time.time()
    solved = cp.solve_points_batch(cfg, specs)
    print(f"stage A: {len(specs)} points in {time.time()-t0:.1f}s", flush=True)

    for r, (case, base), spec in zip(mine.itertuples(), solved, specs):
        stem = f"h{r.hour}_n{r.n_commit}_b{int(r.dc_bus)}_m{int(r.dc_mw)}"
        if case is None:
            print(f"  {stem:<28} DISPATCH {base.get('fail_reason')}", flush=True)
            continue
        for fb in fbs:
            dist = cp.Disturbance(fault_bus=fb, duration=args.duration)
            for a in alphas:
                tag = f"cct_{stem}_f{fb}_a{a:g}"
                c = copy.deepcopy(case)
                for d in c["devices"]:
                    if d["id"] == f"datacenter_{int(r.dc_bus)}":
                        d["params"]["alphaP"] = a
                        d["params"]["alphaI"] = 1.0 - a
                wd = spec["workdir"] / f"f{fb}_a{a:g}"
                rec = dict(base)
                rec.update(cp.run_disturbance(cfg, c, dist, wd, with_cct=True))
                rec.update(
                    tag=tag,
                    hour=int(r.hour),
                    alpha_p_dc=a,
                    dc_bus=int(r.dc_bus),
                    n_commit=int(r.n_commit),
                )
                if rec["status"] != "ok":
                    rec["fail_mode"] = af.fail_mode((wd / "sim.log").read_text())
                cp.persist(cfg, tag, rec, rec.pop("_csv", None))
                print(
                    f"  {stem:<28} f{fb:<4} a{a:g} {rec['status']} "
                    f"cct={rec.get('cct_s')} cens={rec.get('cct_censored')}",
                    flush=True,
                )
    print(f"done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
