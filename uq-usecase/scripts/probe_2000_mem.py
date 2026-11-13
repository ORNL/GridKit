"""Peak memory of one ACTIVSg2000 dynamic run and one stage-A batch.

    /usr/bin/time -v python probe_2000_mem.py dyn
    /usr/bin/time -v python probe_2000_mem.py stagea

Run under /usr/bin/time -v: its "Maximum resident set size" covers this process
and its waited-for children (the DynamicSimulation runner, or the Julia batch).
One case, one run each. This is a quick peek for the per-worker memory number,
not a campaign measurement.
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/scripts")
import pm_utils  # noqa: E402
import scale_probe as sp  # noqa: E402

M2000 = "/kfs2/projects/scidac/scidac-data/ACTIVSg2000/grid-v1/case_ACTIVSg2000.m"
PROJ = Path("/home/isatkaus/gridkit/uq-usecase/pm-solver")


def dyn():
    # Full monitors, tmax 10 s, the same settings as the scale probe in plan §scaling.
    row = sp.run(sp.CASES["ACTIVSg2000"], "ACTIVSg2000", True)
    print(row, flush=True)


def stagea():
    jb = pm_utils.resolve_julia()
    w = Path(f"/scratch/{os.environ['USER']}/scale-probe/stageA2000")
    w.mkdir(parents=True, exist_ok=True)
    cases = [
        {"input": M2000, "output": str(w / f"{pf}.m"), "pf_type": pf}
        for pf in ("dc-opf", "ac")
    ]
    recs, err, rc = pm_utils.run_pm_solve_batch(jb, PROJ, PROJ / "pm_solve.jl", cases, w)
    print(f"rc={rc}", [(r.get("status"), r.get("solve_time")) for r in recs], flush=True)
    if not recs:
        print(err[-2000:])


if __name__ == "__main__":
    {"dyn": dyn, "stagea": stagea}[sys.argv[1]]()
