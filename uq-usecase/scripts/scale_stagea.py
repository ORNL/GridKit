"""Warm per-case stage-A cost: ACTIVSg2000 vs ACTIVSg200.

    python scale_stagea.py

Runs both cases through one batched Julia process so the ~6 s launch and JIT
cost is paid once and the per-case `solve_time` the solver reports is the
quantity compared. Cold single-shot timings are launch-dominated and useless
for projection.
"""

import json
import os
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, "/home/isatkaus/gridkit/uq-usecase/py-utils")
import pm_utils  # noqa: E402

M = {
    "ACTIVSg200": "/kfs2/projects/scidac/scidac-data/ACTIVSg200/grid-v1/case_ACTIVSg200.m",
    "ACTIVSg2000": "/kfs2/projects/scidac/scidac-data/ACTIVSg2000/grid-v1/case_ACTIVSg2000.m",
}
REPEATS = 3
PROJ = Path("/home/isatkaus/gridkit/uq-usecase/pm-solver")
W = Path(f"/scratch/{os.environ['USER']}/scale-probe/stageA")


def main():
    W.mkdir(parents=True, exist_ok=True)
    jb = pm_utils.resolve_julia()
    cases, meta = [], []
    for name, m in M.items():
        if not Path(m).exists():
            print(f"MISSING {m}")
            continue
        for pf in ("dc-opf", "ac"):
            for k in range(REPEATS):
                cases.append(
                    {
                        "input": m,
                        "output": str(W / f"{name}.{pf}.{k}.m"),
                        "pf_type": pf,
                    }
                )
                meta.append((name, pf))

    recs, err, rc = pm_utils.run_pm_solve_batch(
        jb, PROJ, PROJ / "pm_solve.jl", cases, W
    )
    print(f"batch rc={rc}, {len(recs)} results")
    if not recs:
        print(err[-2000:])
        return

    agg = {}
    for (name, pf), r in zip(meta, recs):
        agg.setdefault((name, pf), []).append(
            (r.get("solve_time"), r.get("ok"), r.get("status"))
        )

    print(f"\n{'case':<12} {'stage':<8} {'ok':<4} {'median solve_s':>15}")
    out = {}
    for (name, pf), v in agg.items():
        ts = [t for t, ok, _ in v if t is not None and ok]
        med = st.median(ts) if ts else float("nan")
        out[(name, pf)] = med
        print(f"{name:<12} {pf:<8} {str(all(ok for _, ok, _ in v)):<4} {med:>15.3f}")

    print("\nratio 2000/200:")
    for pf in ("dc-opf", "ac"):
        a, b = out.get(("ACTIVSg200", pf)), out.get(("ACTIVSg2000", pf))
        if a and b:
            print(f"  {pf:<8} {b/a:5.1f}x   ({a:.3f} s -> {b:.3f} s)")
    tot200 = sum(out.get(("ACTIVSg200", p), 0) for p in ("dc-opf", "ac"))
    tot2000 = sum(out.get(("ACTIVSg2000", p), 0) for p in ("dc-opf", "ac"))
    if tot200:
        print(
            f"\n  full stage A chain: {tot200:.3f} s -> {tot2000:.3f} s "
            f"({tot2000/tot200:.1f}x)"
        )

    json.dump(
        {f"{k[0]}|{k[1]}": v for k, v in out.items()},
        open(W / "stagea_times.json", "w"),
        indent=2,
    )


if __name__ == "__main__":
    main()
