"""Committed-capacity margin for each dispatch, to test DC OPF infeasibility.

Usage:
    python dispatch_margin.py <campaign>

For each distinct dispatch (hour, n_commit, dc_bus, dc_mw) in the campaign
records, prints the thermal PMAX committed by the merit-order commitment minus
the thermal energy the dispatch must supply (net_mw, which already includes
the datacenter MW because add_datacenter runs before the net is taken). A negative
margin means no thermal dispatch can serve the load, so DC OPF is infeasible
on capacity alone. Line limits are not tested here.
"""

import collections
import glob
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py-utils"))
import campaign as cp  # noqa: E402

STORE = Path("/kfs2/projects/scidac/scidac-data/gridkit-campaigns")


def main(name):
    cfg = cp.CampaignConfig(name=name, dc_alpha_p=0.8)
    order = cp.merit_order(cfg)
    gen = cfg.seed["gen"]
    pmax = dict(zip(cfg.gid, gen["PMAX"].astype(float)))
    nd = set(cfg.nd)

    disp = {}
    for f in glob.glob(str(STORE / name / "records" / "*.json")):
        r = json.load(open(f))
        key = (r["hour"], r["n_commit"], r.get("dc_bus"), r.get("dc_mw", 0.0))
        disp.setdefault(key, {"status": set(), "net": r["net_mw"]})
        disp[key]["status"].add(r["status"])
        if r["status"] == "solve_failed":
            disp[key]["stage"] = f"{r['fail_stage']} {r['fail_reason']}"

    rows = []
    for (hour, n, dc_bus, dc_mw), d in disp.items():
        on = set(order["gid"].iloc[:n])
        committed = sum(pmax[g] for g in cfg.gid if g in on and g not in nd)
        need = d["net"]
        status = "solve_failed" if "solve_failed" in d["status"] else "solved"
        rows.append((hour, n, dc_bus, dc_mw, round(committed, 1), round(need, 1),
                     round(committed - need, 1), status, d.get("stage", "")))

    print("hour n_commit dc_bus dc_mw committed_pmax need_mw margin_mw status stage")
    for r in sorted(rows, key=lambda x: x[6]):
        if r[6] < 150 or r[7] == "solve_failed":
            print(*r)
    print("\nsolve_failed margins:", sorted(r[6] for r in rows if r[7] == "solve_failed"))
    solved = [r[6] for r in rows if r[7] == "solved" and r[3] > 0]
    print("solved datacenter dispatches: n=%d min margin=%s" % (len(solved), min(solved)))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
