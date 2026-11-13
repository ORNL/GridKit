"""Tally dispatch-stage (solve_failed) and sim_failed rows in a campaign.

    python failure_audit.py <campaign> [--ref <campaign>] [--out <file>]

Reads records/*.json under the campaign store and prints:
  - solve_failed by (fail_stage, fail_reason), then each failing dispatch
    (hour, n_commit, dc_bus, dc_mw) with its stage and fault-row count. The
    fault rows of one dispatch share the same operating point.
  - sim_failed by (sim_rc, sim_errors). Each run's scratch sim.log is read and
    classified from the solver message: IDACalcIC (initial conditions, at the
    event) or IDASolve (integration; the failure time is printed).
  - with --ref, sim_failed rows are split into those that also fail in the
    reference campaign and those that do not.

sim.log lives in scratch and may be purged; the record counts do not need it.
"""

import argparse
import collections
import json
import re
from pathlib import Path

STORE = Path("/kfs2/projects/scidac/scidac-data/gridkit-campaigns")
SCRATCH = Path("/scratch/isatkaus/gridkit-campaigns")


def load(name):
    return {
        r["tag"]: r
        for r in (json.load(open(f)) for f in sorted((STORE / name / "records").glob("*.json")))
    }


def sim_kind(campaign, tag):
    """Classify one sim_failed run from its scratch sim.log."""
    disp, fault = tag.split("__")
    p = SCRATCH / campaign / disp / fault / "sim.log"
    if not p.exists():
        return "no log", None
    text = p.read_text()
    if "IDASolve" in text:
        m = re.search(r"At t = ([0-9.eE+-]+)", text)
        return "IDASolve", float(m.group(1)) if m else None
    if "IDACalcIC" in text:
        return "IDACalcIC", None
    return "other", None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("campaign")
    ap.add_argument("--ref", help="reference campaign for the sim_failed overlap")
    ap.add_argument("--out", help="also write the report to this file")
    a = ap.parse_args()

    rows = load(a.campaign)
    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out(f"== {a.campaign}: {len(rows)} records")
    out("status: " + str(dict(collections.Counter(r["status"] for r in rows.values()))))

    solve = {t: r for t, r in rows.items() if r["status"] == "solve_failed"}
    out("solve_failed by (stage, reason): " + str(dict(collections.Counter(
        (r["fail_stage"], r["fail_reason"]) for r in solve.values()))))
    disp = collections.defaultdict(list)
    for r in solve.values():
        disp[(r["hour"], r["n_commit"], r["dc_bus"], r["dc_mw"])].append(r["fail_stage"])
    out(f"solve_failed dispatches: {len(disp)}")
    for k in sorted(disp):
        out(f"  hour={k[0]} n_commit={k[1]} dc_bus={k[2]} dc_mw={k[3]} "
            f"stage={disp[k][0]} rows={len(disp[k])}")

    sim = {t: r for t, r in rows.items() if r["status"] == "sim_failed"}
    out("sim_failed by (sim_rc, sim_errors): " + str(dict(collections.Counter(
        (r["sim_rc"], r["sim_errors"]) for r in sim.values()))))

    if a.ref:
        ref = {t for t, r in load(a.ref).items() if r["status"] == "sim_failed"}
        common = {t for t in sim if t in ref}
        groups = [("also in " + a.ref, common), ("only here", set(sim) - common)]
        out(f"sim_failed also in {a.ref}: {len(common)}; only in {a.campaign}: {len(set(sim) - common)}")
    else:
        groups = [("all", set(sim))]

    for label, S in groups:
        kinds = collections.Counter()
        solve_times = []
        for t in sorted(S):
            k, tt = sim_kind(a.campaign, t)
            kinds[k] += 1
            if k == "IDASolve":
                solve_times.append((t, tt))
        out(f"{label}: {len(S)} rows; log kind {dict(kinds)}")
        out(f"  dc_bus: {dict(collections.Counter(sim[t]['dc_bus'] for t in S))}")
        out(f"  fault_bus (top 6): {dict(collections.Counter(sim[t]['fault_bus'] for t in S).most_common(6))}")
        out(f"  dispatches: {len({t.split('__')[0] for t in S})}")
        for t, tt in solve_times:
            out(f"  IDASolve {t} at t={tt}")

    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        Path(a.out).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
