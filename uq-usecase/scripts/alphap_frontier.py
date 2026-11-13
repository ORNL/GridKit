"""Map the feasible-alphaP frontier for an added datacenter load.

    python alphap_frontier.py --name afront1            # build design + run all
    python alphap_frontier.py --name afront1 --list     # design only, no run

Separates two explanations for why a high constant-power share is sometimes
unsimulable:

  * a BUS property  -> the frontier is the same at every operating condition
  * a CONDITION property -> the frontier moves with hour or commitment, which
    would point at the AC solution / UC decision the dynamics initialise from

alphaP does not enter the power flow at all: at nominal voltage every ZIP mix
draws the same P, so one dispatch serves the whole ladder and the AC solution
is identical across it. That makes this a pure stage-B fan-out, and it is also
what makes the comparison clean.
"""

import argparse
import copy
import json
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import campaign as cp  # noqa: E402

DATA = "/kfs2/projects/scidac/scidac-data/ACTIVSg200/grid-v1"
HOURS = (7549, 5572, 6164, 4737)
ALPHAS = (0.2, 0.5, 0.7, 0.85, 1.0)
KEEP = (93, 108, 129, 158, 172)


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


def fail_mode(log: str) -> str:
    """Which of the two observed failures this is, from the runner log."""
    if "IDACalcIC" in log:
        return "IDACalcIC"
    if "IDASolve" in log:
        return "IDASolve"
    return "other"


def build_design(cfg, at, buses, hours, fracs, mws):
    order = cp.merit_order(cfg)
    rows = []
    for h in hours:
        demand, wind = at(h)
        frames = cp.apply_hour(cfg, demand, wind)
        net = frames["bus"]["PD"].sum() - frames["gen"].loc[cfg.nd, "PG"].sum()
        for mw in mws:
            lo, hi = cp.commitment_band(cfg, float(net) + mw, order)
            for f in fracs:
                n = int(round(lo + f * (hi - lo)))
                for b in buses:
                    rows.append(
                        {
                            "hour": h,
                            "n_commit": n,
                            "dc_bus": int(b),
                            "dc_mw": float(mw),
                        }
                    )
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--name", required=True)
    p.add_argument("--dc-mw", default="100,200")
    p.add_argument("--n-buses", type=int, default=16)
    p.add_argument("--fault-buses", default="66")
    p.add_argument("--hours", default=",".join(str(h) for h in HOURS))
    p.add_argument("--duration", type=float, default=0.1)
    p.add_argument("--commit-fracs", default="0.25,0.75")
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--count", type=int, default=10**9)
    p.add_argument("--list", action="store_true")
    args = p.parse_args()

    # Unlimited Q, as these campaigns were run (the config default is switched).
    cfg = cp.CampaignConfig(name=args.name, enforce_q_lims=False)
    missing = cp.check_runner(cfg)
    if missing:
        print("FATAL: runner cannot load its libraries:", *missing, sep="\n  ")
        sys.exit(2)

    at = profiles(cfg)
    cand = cp.candidate_buses(cfg, n=args.n_buses, require_load=True, keep=KEEP)
    buses = sorted(set(cand.bus.astype(int)) | {34})  # 34 is greenfield, no host load
    fracs = [float(x) for x in args.commit_fracs.split(",")]
    mws = [float(x) for x in args.dc_mw.split(",")]
    fbs = [int(x) for x in args.fault_buses.split(",")]
    hours = [int(x) for x in args.hours.split(",")]
    D = build_design(cfg, at, buses, hours, fracs, mws)

    print(
        f"design: {len(buses)} buses x {len(hours)} hours x {len(fracs)} commit "
        f"x {len(mws)} MW = {len(D)} dispatches, "
        f"x {len(fbs)} fault x {len(ALPHAS)} alphaP "
        f"= {len(D)*len(fbs)*len(ALPHAS)} sims"
    )
    print("buses:", buses)
    if args.list:
        print(D.to_string(index=False))
        return

    mine = D.iloc[args.start : args.start + args.count]
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

    out = []
    for r, (case, base), spec in zip(mine.itertuples(), solved, specs):
        stem = f"h{r.hour}_n{r.n_commit}_b{int(r.dc_bus)}_m{int(r.dc_mw)}"
        if case is None:
            rec = dict(base, tag=stem, hour=int(r.hour), alpha_p_dc=None)
            cp.persist(cfg, f"af_{stem}", rec)
            out.append(dict(rec, frontier=None, note="dispatch infeasible"))
            print(f"  {stem:<28} DISPATCH {base.get('fail_reason')}", flush=True)
            continue

        for fb in fbs:
            dist = cp.Disturbance(fault_bus=fb, duration=args.duration)
            results = {}
            for a in ALPHAS:
                tag = f"af_{stem}_f{fb}_a{a:g}"
                fj = cfg.store / "records" / f"{tag}.json"
                if fj.exists():
                    results[a] = json.loads(fj.read_text()).get("status")
                    continue
                c = copy.deepcopy(case)
                for d in c["devices"]:
                    if d["id"] == f"datacenter_{int(r.dc_bus)}":
                        d["params"]["alphaP"] = a
                        d["params"]["alphaI"] = 1.0 - a
                wd = spec["workdir"] / f"f{fb}_a{a:g}"
                rec = dict(base)
                rec.update(cp.run_disturbance(cfg, c, dist, wd, with_cct=False))
                rec.update(tag=tag, hour=int(r.hour), alpha_p_dc=a)
                if rec["status"] != "ok":
                    rec["fail_mode"] = fail_mode((wd / "sim.log").read_text())
                cp.persist(cfg, tag, rec, rec.pop("_csv", None))
                results[a] = rec["status"]

            ok = [a for a in ALPHAS if results.get(a) == "ok"]
            frontier = max(ok) if ok else None
            out.append(
                dict(
                    hour=int(r.hour),
                    n_commit=int(r.n_commit),
                    dc_bus=int(r.dc_bus),
                    dc_mw=float(r.dc_mw),
                    fault_bus=fb,
                    net_mw=base.get("net_mw"),
                    frontier=frontier,
                    **{f"a{a:g}": results.get(a) for a in ALPHAS},
                )
            )
            print(
                f"  {stem:<28} f{fb:<4} frontier={frontier} "
                + " ".join(f"{a:g}:{results.get(a)}" for a in ALPHAS),
                flush=True,
            )

    df = pd.DataFrame(out)
    # per-slice name: workers share the store, a fixed filename is a race
    dest = cfg.store / f"frontier_{args.start:04d}.csv"
    df.to_csv(dest, index=False)
    print(f"\nwrote {dest}  ({len(df)} rows, {time.time()-t0:.1f}s)")

    if "frontier" in df and df.frontier.notna().any():
        print("\nfrontier by bus and size (cols = operating condition):")
        pv = df.pivot_table(
            index=["dc_bus", "dc_mw"], columns=["hour", "n_commit"], values="frontier"
        )
        print(pv.to_string())
        spread = (pv.max(axis=1) - pv.min(axis=1)).rename("spread")
        print("\nspread of frontier across operating conditions:")
        print(spread.to_string())
        print(
            "\nall spreads 0  -> the frontier is a BUS (and size) property;\n"
            "any spread > 0 -> the operating condition moves it, which points\n"
            "                 at the AC solution / commitment the IC starts from."
        )


if __name__ == "__main__":
    main()
