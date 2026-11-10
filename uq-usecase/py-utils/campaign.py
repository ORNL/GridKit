"""Campaign machinery: operating conditions -> solved case -> dynamics -> QoI.

Factored out of the Phase 1 scripts so Phase 2 (aleatoric sampling) and Phase 3
(datacenter siting) share one implementation.

The pipeline for one condition:

    hour profile  ->  frames (load + wind)
                  ->  merit-order commitment
                  ->  loss-aware DC OPF -> AC PF
                  ->  case.json  (wind converted to a non-inertial stand-in)
                  ->  DynamicSimulation
                  ->  scalar QoI

See work-notes/datacenter_load_uq_plan.md for why each step is the way it is.
"""

from __future__ import annotations

import copy
import json
import math
import os
import subprocess
import time
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd

import aux_to_m
import m_to_case_utils as mu
import perturb
import pm_utils
import solver_utils as su

__all__ = [
    "CampaignConfig",
    "OperatingPoint",
    "Disturbance",
    "prepare_storage",
    "check_runner",
    "merit_order",
    "commitment_band",
    "commit_n",
    "apply_hour",
    "solve_chain",
    "solve_chain_batch",
    "wind_to_negative_load",
    "reduce_qoi",
    "simulate",
    "is_unstable",
    "find_cct",
    "write_timeseries",
    "persist",
    "add_datacenter",
    "split_out_datacenter",
    "apply_disturbance",
    "solve_point",
    "solve_points_batch",
    "candidate_buses",
    "run_disturbance",
    "run_condition",
]

BASEMVA = 100.0

# DynamicSimulation needs the libraries it was built against (libomp lives
# here, and a compute node has no /lib64/libomp.so even though a login node
# may). This MUST stay scoped to that child process: putting it on the job's
# LD_LIBRARY_PATH makes Julia load this SuiteSparse instead of its own and
# PowerModels dies in SuiteSparse_jll init.
_GRIDKIT_DEPS = "/nopt/nrel/apps/cpu_stack/software/gridkit"
RUNNER_LD = ":".join(
    [
        f"{_GRIDKIT_DEPS}/llvm/16.0.6/lib",
        f"{_GRIDKIT_DEPS}/suitesparse/7.8.2/lib64",
        f"{_GRIDKIT_DEPS}/ipopt/3.14.14/lib",
        f"{_GRIDKIT_DEPS}/sundials/develop/lib64",
    ]
)


def runner_env() -> dict:
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = RUNNER_LD + ":" + env.get("LD_LIBRARY_PATH", "")
    return env


def check_runner(cfg: CampaignConfig) -> list[str]:
    """Unresolved shared libraries for the DynamicSimulation binary.

    Without this a worker discovers a missing libomp only after paying a ~30 s
    dispatch solve, and the failure then looks like a modelling problem rather
    than a missing library.
    """
    if not Path(cfg.runner).exists():
        return [f"runner not found: {cfg.runner}"]
    r = subprocess.run(
        ["ldd", str(cfg.runner)], capture_output=True, text=True, env=runner_env()
    )
    return [ln.strip() for ln in r.stdout.splitlines() if "not found" in ln]


# --- storage -------------------------------------------------------------
# Scratch holds per-run churn and is purged; store holds the two artifacts
# worth keeping. Splitting them keeps a small-file storm off ProjectFS.
SCRATCH_ROOT = Path(f"/scratch/{os.environ.get('USER', 'nobody')}/gridkit-campaigns")
STORE_ROOT = Path("/kfs2/projects/scidac/scidac-data/gridkit-campaigns")

# Stripe narrow: the workload is file-per-process and no file exceeds ~27 MB,
# so throughput comes from many files spread over OSTs, not from splitting one
# file. The scratch default (c=8 on the disk pool past 2 MB) puts a 25 MB file
# on 9 OST objects; these put it on 2.
_LAYOUT_SCRATCH = "-E 2M -S 1M -c 1 -p flash -E -1 -S 4M -c 1 -p flash"
_LAYOUT_SERIES = "-E 2M -S 1M -c 1 -p disk -E -1 -S 4M -c 1 -p disk"
# records/ holds only ~2 KB json, so Data-on-MDT keeps it off the OSTs entirely.
# Safe ONLY because nothing large lives here: on a mixed directory every file
# over 256K would also burn 256K of MDT for no benefit.
_LAYOUT_RECORDS = "-E 256K -L mdt -E -1 -S 1M -c 1 -p disk"


def _setstripe(path: Path, layout: str) -> bool:
    """Apply a Lustre layout to a directory. Returns False off Lustre."""
    path.mkdir(parents=True, exist_ok=True)
    try:
        r = subprocess.run(
            ["lfs", "setstripe", *layout.split(), str(path)],
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, OSError):
        return False
    return r.returncode == 0


def prepare_storage(scratch: Path, store: Path) -> dict:
    """Create the campaign trees with their Lustre layouts. Idempotent.

    A layout only affects files created afterwards, so this must run before the
    sweep, not after.
    """
    return {
        "scratch": _setstripe(Path(scratch), _LAYOUT_SCRATCH),
        "series": _setstripe(Path(store) / "series", _LAYOUT_SERIES),
        "records": _setstripe(Path(store) / "records", _LAYOUT_RECORDS),
    }


# --- the sample, split by where each parameter binds ----------------------
# Stage A costs a Julia solve; stages B and C are a JSON edit and a re-run.
# Keeping them apart is what lets one solved point fan out over many
# disturbances. Both are frozen so a point can be deduplicated or cached.


@dataclass(frozen=True)
class OperatingPoint:
    """Stage A: anything that changes the dispatch, so needs OPF + PF."""

    hour: int
    n_commit: int
    dc_bus: int | None = None
    dc_mw: float = 0.0
    load_scale: float = 1.0
    wind_scale: float = 1.0

    @property
    def tag(self) -> str:
        p = [f"h{self.hour}", f"n{self.n_commit}"]
        if self.dc_bus is not None and self.dc_mw:
            p.append(f"dc{self.dc_bus}-{self.dc_mw:g}")
        if self.load_scale != 1.0:
            p.append(f"l{self.load_scale:g}")
        if self.wind_scale != 1.0:
            p.append(f"w{self.wind_scale:g}")
        return "_".join(p)


@dataclass(frozen=True)
class Disturbance:
    """Stages B and C: everything variable without re-solving the dispatch."""

    fault_bus: int = 66
    duration: float | None = 0.1  # None -> bisect for CCT
    t_on: float = 1.0
    tmax: float = 10.0
    h_scale: float = 1.0
    alpha_p: float | None = None
    bus_mon: tuple[str, ...] = ("Vr", "Vi")
    genrou_mon: tuple[str, ...] = ("delta", "omega")

    @property
    def t_off(self) -> float:
        if self.duration is None:
            raise ValueError("t_off is undefined while duration is None")
        return self.t_on + self.duration

    @property
    def tag(self) -> str:
        p = [f"b{self.fault_bus}"]
        p.append("cct" if self.duration is None else f"d{self.duration:g}")
        if self.h_scale != 1.0:
            p.append(f"h{self.h_scale:g}")
        if self.alpha_p is not None:
            p.append(f"ap{self.alpha_p:g}")
        return "_".join(p)


class CampaignConfig:
    """Everything a campaign run needs that does not change between conditions."""

    def __init__(
        self,
        spec=None,
        name: str = "campaign",
        scratch: str | Path | None = None,
        store: str | Path | None = None,
        work: str | Path | None = None,
        fault_bus: int = 66,
        t_on: float = 1.0,
        t_off: float = 1.1,
        tmax: float = 10.0,
        reserve: float = 1.10,
        min_thermal: int = 8,
        loss_iters: int = 1,
        cct_max: float = 0.5,
        keep_timeseries: bool = True,
        dc_alpha_p: float = 0.2,
    ):
        self.spec = spec or mu.ACTIVSG200
        self.name = name
        # `work` is the old single-directory argument; it now means scratch.
        self.scratch = Path(scratch or work or SCRATCH_ROOT / name)
        self.store = Path(store or STORE_ROOT / name)
        self.work = self.scratch
        self.keep_timeseries = keep_timeseries
        # alphaP=1.0 is not simulable under a fault once the datacenter is large
        # relative to its host load; 0.2 is, and measured CCT is indistinguishable
        # wherever both run (see plan doc, section 4).
        self.dc_alpha_p = dc_alpha_p
        prepare_storage(self.scratch, self.store)
        self.fault_bus = fault_bus
        self.t_on, self.t_off, self.tmax = t_on, t_off, tmax
        # Campaign-wide default; a sweep overrides it per sample.
        self.disturbance = Disturbance(
            fault_bus=fault_bus, duration=t_off - t_on, t_on=t_on, tmax=tmax
        )
        self.reserve, self.min_thermal = reserve, min_thermal
        self.loss_iters = loss_iters
        # Capped at a credible protection time, not extended until instability.
        self.cct_max = cct_max

        self.runner = (
            Path.home() / "gridkit/build/application/PhasorDynamics/DynamicSimulation"
        )
        shipped = mu.load_case(self.spec.case_json_path)
        aux, dyr = mu.parse_aux(self.spec.aux), mu.parse_dyr(self.spec.dyr)
        rebuilt, _ = mu.build_case_from_aux_dyr(self.spec, shipped, aux, dyr)
        overlay = mu.make_creator_overlay(shipped, rebuilt, tol=1e-4)
        self.base_case, _ = mu.build_case_from_aux_dyr(
            self.spec, shipped, aux, dyr, creator_overlay=overlay
        )
        self.dyr = dyr
        self.base_solver = json.loads(Path(self.spec.solver_json_path).read_text())
        self.f_base = float(self.base_case["params"]["freq_base"])

        # Every worker derives the same seed.m from the same aux, but they all
        # share this path on Lustre, so a plain write interleaves. Write
        # per-process and rename: os.replace is atomic, so a reader always sees
        # a complete file.
        seed_m = self.work / "seed.m"
        tmp = self.work / f"seed.{os.getpid()}.tmp"
        aux_to_m.write_matpower_m(aux, tmp, case_name="seed")
        os.replace(tmp, seed_m)
        self.seed = aux_to_m.m_to_frames(seed_m)
        self.gid = list(perturb.gen_index(self.seed))
        self.fuel = [str(x).lower() for x in self.seed["genfuel"]]
        self.nd = np.array([f in perturb.NONDISPATCHABLE_FUELS for f in self.fuel])
        self.hmva = {
            mu.gen_id_of(d["id"], self.spec): d["params"]["H"] * d["params"]["mva"]
            for d in self.base_case["devices"]
            if d["class"] == "Genrou"
        }
        self.bus_of = {
            g: int(self.seed["gen"]["GEN_BUS"].iloc[i]) for i, g in enumerate(self.gid)
        }

    def inertia_of(self, case) -> float:
        """Fleet H*MVA read from the case, so Disturbance.h_scale is reflected."""
        return sum(
            d["params"]["H"] * d["params"]["mva"]
            for d in case["devices"]
            if d["class"] == "Genrou"
        )


def merit_order(cfg: CampaignConfig) -> pd.DataFrame:
    """Thermal units cheapest-first, nuclear must-run at the head."""
    g = cfg.seed["gen"].copy()
    g["gid"], g["fuel"] = cfg.gid, cfg.fuel
    g["c1"] = cfg.seed["gencost"]["C1"].values
    g["hmva"] = [cfg.hmva.get(x, 0.0) for x in cfg.gid]
    th = g[~cfg.nd].sort_values(["c1", "PMAX"])
    return pd.concat([th[th.fuel == "nuclear"], th[th.fuel != "nuclear"]])


def commitment_band(cfg: CampaignConfig, net_mw: float, order=None):
    """Feasible count of committed thermal units at this net load.

    Lower bound: adequacy (capacity >= net * reserve) and the `min_thermal`
    floor, below which the network loses voltage support and IDACalcIC fails.
    Upper bound: sum(PMIN) <= net, above which the OPF is infeasible.

    This is a necessary condition only. The OPF solve is the real test.
    """
    order = merit_order(cfg) if order is None else order
    cap = order["PMAX"].cumsum().values
    pmin = order["PMIN"].cumsum().values
    ok = [
        n
        for n in range(1, len(order) + 1)
        if cap[n - 1] >= net_mw * cfg.reserve
        and pmin[n - 1] <= net_mw
        and n >= cfg.min_thermal
    ]
    return (ok[0], ok[-1]) if ok else (None, None)


def commit_n(cfg: CampaignConfig, frames: dict, n: int, order=None) -> dict:
    """Commit the n cheapest thermal units. Non-dispatchable units stay on."""
    order = merit_order(cfg) if order is None else order
    on = set(order["gid"].iloc[:n])
    frames["gen"]["GEN_STATUS"] = [
        1 if (nd or g in on) else 0 for g, nd in zip(cfg.gid, cfg.nd)
    ]
    return frames


def apply_hour(cfg: CampaignConfig, demand: dict, wind: dict) -> dict:
    """Frames with this hour's bus loads and wind output. Wind is pinned."""
    f = perturb.copy_frames(cfg.seed)
    cur = f["bus"].set_index(f["bus"].BUS_I.astype(int))["PD"]
    f = perturb.scale_load(
        f, {b: (mw / cur[b] if cur[b] else 0.0) for b, mw in demand.items()}
    )
    f = perturb.set_gen_output(f, pg_mw=wind)
    g = f["gen"]
    g.loc[cfg.nd, "PMIN"] = g.loc[cfg.nd, "PG"]
    g.loc[cfg.nd, "PMAX"] = g.loc[cfg.nd, "PG"]
    return f


def _fail(workdir: Path, stage: str, it: int, err: str) -> dict:
    """Record a solver failure without putting its log in the returned dict.

    Solver logs contain newlines, and writing one into a CSV field corrupts
    every subsequent row. The log goes to a file; the record keeps a one-line
    reason that is safe to serialise.
    """
    (Path(workdir) / f"fail_{stage}_{it}.log").write_text(err)
    reason = "unknown"
    for key in (
        "LOCALLY_INFEASIBLE",
        "INFEASIBLE",
        "ITERATION_LIMIT",
        "TIME_LIMIT",
        "NUMERICAL_ERROR",
    ):
        if key in err:
            reason = key
            break
    return {"fail_stage": f"{stage}@{it}", "fail_reason": reason}


def solve_chain(cfg: CampaignConfig, frames: dict, workdir: Path):
    """Loss-aware DC OPF then AC PF. Returns (solved frames | None, info)."""
    workdir.mkdir(parents=True, exist_ok=True)
    pd_total = frames["bus"]["PD"].sum()
    loss, solved, info = 0.0, None, {}
    for i in range(cfg.loss_iters + 1):
        # DC OPF is lossless, so inflate load by the measured loss and let the
        # next AC PF land with a near-zero slack correction.
        infl = perturb.scale_load(perturb.copy_frames(frames), 1.0 + loss / pd_total)
        aux_to_m.write_matpower_m(infl, workdir / "dcin.m", case_name="c")
        _, err, rc = pm_utils.solve_pf(
            workdir / "dcin.m", workdir / "dc.m", pf_type="dc-opf"
        )
        if rc != 0:
            return None, _fail(workdir, "dc-opf", i, err)
        f2 = perturb.copy_frames(frames)  # original load, OPF dispatch
        f2["gen"]["PG"] = aux_to_m.m_to_frames(workdir / "dc.m")["gen"]["PG"].values
        aux_to_m.write_matpower_m(f2, workdir / "acin.m", case_name="c")
        _, err, rc = pm_utils.solve_pf(
            workdir / "acin.m", workdir / "ac.m", pf_type="ac"
        )
        if rc != 0:
            return None, _fail(workdir, "ac-pf", i, err)
        solved = aux_to_m.m_to_frames(workdir / "ac.m")
        on = solved["gen"][solved["gen"].GEN_STATUS == 1]
        loss = on["PG"].sum() - pd_total
        info = {
            "losses_mw": round(float(loss), 2),
            "max_over_pmax": round(
                float((on["PG"] - on["PMAX"]).clip(lower=0).max()), 3
            ),
            "vmin": round(float(solved["bus"]["VM"].min()), 4),
            "vmax": round(float(solved["bus"]["VM"].max()), 4),
        }
    return solved, info


def solve_chain_batch(cfg: CampaignConfig, items: list) -> list:
    """`solve_chain` for many conditions, in one Julia process per stage.

    Within a condition the DC OPF and AC PF are strictly sequential, so they
    cannot be batched together. The *same stage across conditions* is
    independent, so the conditions advance in lockstep: one batched DC OPF for
    all of them, then one batched AC PF, repeated per loss iteration. That is
    `2*(loss_iters+1)` Julia launches for the whole slice instead of per
    condition, and launch cost dominates (~10 s versus ~0.1 s to solve).

    `items` is a list of `(frames, workdir)`. Returns a list of
    `(solved_frames | None, info)` in the same order.
    """
    items = [(f, Path(w)) for f, w in items]
    for _, w in items:
        w.mkdir(parents=True, exist_ok=True)
    n = len(items)
    pd_tot = [float(f["bus"]["PD"].sum()) for f, _ in items]
    loss = [0.0] * n
    out = [None] * n
    info = [{} for _ in range(n)]
    dead = [None] * n  # fail dict once a condition drops out

    jb = pm_utils.resolve_julia()
    proj = pm_utils.PF_BACKENDS["pm"]["project"]
    script = pm_utils.PF_BACKENDS["pm"]["script"]
    # Per-process: workers share the campaign scratch root, so a fixed name
    # here would have every worker overwrite the others' manifest and results.
    batch_dir = items[0][1].parent / f"_batch_{os.getpid()}"

    def run_stage(stage, pf_type, prepare):
        """prepare(k) writes the input .m and returns (input, output) paths."""
        cases, idx = [], []
        for k in range(n):
            if dead[k] is not None:
                continue
            src, dst = prepare(k)
            cases.append({"input": str(src), "output": str(dst), "pf_type": pf_type})
            idx.append(k)
        if not cases:
            return
        recs, err, rc = pm_utils.run_pm_solve_batch(jb, proj, script, cases, batch_dir)
        if rc != 0:
            for k in idx:
                dead[k] = _fail(items[k][1], stage, it, err)
            return
        if len(recs) != len(idx):
            # Never silently zip-truncate: a condition with no record would
            # advance to the next stage with no input file.
            for k in idx:
                dead[k] = _fail(
                    items[k][1],
                    stage,
                    it,
                    f"batch returned {len(recs)} records for "
                    f"{len(idx)} cases\n{err}",
                )
            return
        for rec, k in zip(recs, idx):
            if not rec.get("ok"):
                dead[k] = _fail(
                    items[k][1],
                    stage,
                    it,
                    str(rec.get("status", "")) + "\n" + str(rec.get("error", "")),
                )

    for it in range(cfg.loss_iters + 1):

        def prep_dc(k):
            frames, wd = items[k]
            # DC OPF is lossless, so inflate load by the measured loss and let
            # the next AC PF land with a near-zero slack correction.
            infl = perturb.scale_load(
                perturb.copy_frames(frames), 1.0 + loss[k] / pd_tot[k]
            )
            aux_to_m.write_matpower_m(infl, wd / "dcin.m", case_name="c")
            return wd / "dcin.m", wd / "dc.m"

        run_stage("dc-opf", "dc-opf", prep_dc)

        def prep_ac(k):
            frames, wd = items[k]
            f2 = perturb.copy_frames(frames)  # original load, OPF dispatch
            f2["gen"]["PG"] = aux_to_m.m_to_frames(wd / "dc.m")["gen"]["PG"].values
            aux_to_m.write_matpower_m(f2, wd / "acin.m", case_name="c")
            return wd / "acin.m", wd / "ac.m"

        run_stage("ac-pf", "ac", prep_ac)

        for k in range(n):
            if dead[k] is not None:
                continue
            solved = aux_to_m.m_to_frames(items[k][1] / "ac.m")
            on = solved["gen"][solved["gen"].GEN_STATUS == 1]
            loss[k] = float(on["PG"].sum() - pd_tot[k])
            out[k] = solved
            info[k] = {
                "losses_mw": round(loss[k], 2),
                "max_over_pmax": round(
                    float((on["PG"] - on["PMAX"]).clip(lower=0).max()), 3
                ),
                "vmin": round(float(solved["bus"]["VM"].min()), 4),
                "vmax": round(float(solved["bus"]["VM"].max()), 4),
            }

    return [
        (None, dead[k]) if dead[k] is not None else (out[k], info[k]) for k in range(n)
    ]


def wind_to_negative_load(cfg: CampaignConfig, case: dict, frames: dict):
    """Replace wind machines with a non-inertial constant-current injection.

    alphaI=1, not alphaP=1: a constant-power negative injection makes
    IDACalcIC fail at high wind penetration regardless of how many synchronous
    machines are committed (measured). Constant current is also the better
    analogue of a current-controlled converter.
    """
    case = copy.deepcopy(case)
    g = frames["gen"]
    wind = {
        gid: (float(g["PG"].iloc[i]), float(g["QG"].iloc[i]))
        for i, (gid, nd) in enumerate(zip(cfg.gid, cfg.nd))
        if nd and g["GEN_STATUS"].iloc[i] == 1
    }
    drop = {gid + s for gid in wind for s in ("_genrou", "_tgov1", "_sexs_pti")}
    case["devices"] = [d for d in case["devices"] if d["id"] not in drop]
    for gid, (pg, qg) in wind.items():
        case["devices"].append(
            {
                "class": "LoadZIP",
                "ports": {"bus": cfg.bus_of[gid]},
                "id": f"{gid}_ibr_stand_in",
                "params": {
                    "Pnom": -pg / BASEMVA,
                    "Qnom": -qg / BASEMVA,
                    "alphaI": 1.0,
                    "alphaP": 0.0,
                },
            }
        )
    return case, len(wind)


def reduce_qoi(
    cfg: CampaignConfig, csv: Path, dist: Disturbance, blank: float = 0.05
) -> dict:
    """One monitor csv -> scalars. Independent of machine count.

    Two frequency views, not interchangeable:

    * `coi_*` uses the mean bus frequency, a proxy for the centre of inertia.
      This is the system-wide quantity and the one to report.
    * `local_max_rocof_Hz_s` is the max over buses and post-blank time. Kept
      for continuity but it is **not** a system ROCOF: measured on this case it
      lands at the faulted bus at `t_off + blank` every time, so it tracks the
      tail of the clearing transient rather than an inertial response.

    A self-clearing bus fault leaves generation and load balanced, so none of
    these are frequency-event metrics in the UFLS sense. GridKit exposes only
    FAULT_ON/FAULT_OFF, so a true power-imbalance event is not available.
    """
    df = pd.read_csv(csv)
    df.columns = [c.strip() for c in df.columns]
    # Events emit two rows at the same t; a zero dt would give inf gradients.
    df = df.drop_duplicates(subset="t", keep="last").reset_index(drop=True)
    t = df["t"].to_numpy()
    vr = df[[c for c in df.columns if c.endswith("_Vr")]].to_numpy()
    vi = df[[c for c in df.columns if c.endswith("_Vi")]].to_numpy()
    w = df[[c for c in df.columns if c.endswith("_omega")]].to_numpy()
    dl = np.rad2deg(df[[c for c in df.columns if c.endswith("_delta")]].to_numpy())
    if vr.size == 0:
        raise RuntimeError(f"{csv}: no bus voltage channels, monitors not enabled")
    vmag = np.hypot(vr, vi)
    va = np.unwrap(np.arctan2(vi, vr), axis=0)
    during = (t >= dist.t_on) & (t <= dist.t_off)
    post = t > dist.t_off
    meas = t > (dist.t_off + blank)
    busf = cfg.f_base + np.gradient(va, t, axis=0) / (2.0 * np.pi)
    coi = busf.mean(axis=1)
    coi_rocof = np.abs(np.gradient(coi, t))
    spread = np.rad2deg(va.max(axis=1) - va.min(axis=1))
    at_nadir = busf[meas][np.argmin(coi[meas])]
    return {
        "Vmin_fault": round(float(vmag[during].min()), 4),
        "Vmin_post": round(float(vmag[post].min()), 4),
        "coi_nadir_Hz": round(float(coi[meas].min()), 4),
        "coi_max_rocof_Hz_s": round(float(coi_rocof[meas].max()), 4),
        "freq_spread_at_nadir_Hz": round(float(at_nadir.max() - at_nadir.min()), 4),
        "local_max_rocof_Hz_s": round(
            float(np.abs(np.gradient(busf, t, axis=0))[meas].max()), 3
        ),
        "nadir_mach_Hz": round(float((cfg.f_base * (1 + w)).min()), 4),
        "max_va_spread_deg": round(float(spread.max()), 2),
        "final_va_spread_deg": round(float(spread[-1]), 2),
        "max_delta_spread_deg": round(
            float((dl.max(axis=1) - dl.min(axis=1)).max()), 2
        ),
        "n_machines": int(w.shape[1]),
        "any_nan": bool(np.isnan(vmag).any() or np.isnan(w).any()),
        "t_end": float(t[-1]),
    }


def simulate(cfg: CampaignConfig, workdir: Path, dist: Disturbance):
    """One dynamic run of a given disturbance. Returns (rc, csv, log).

    Expects `case.json` to already be in `workdir`, with the stage-B parts of
    `dist` applied. Clears stale csv first so the caller cannot pick up a
    previous run's output.
    """
    for f in workdir.glob("*.csv"):
        f.unlink()
    case = mu.load_case(workdir / "case.json")
    eid = su.fault_element_id(case, dist.fault_bus)
    solver = su.apply_solver_overrides(
        cfg.base_solver,
        {
            "tmax": dist.tmax,
            "events": [
                {"time": dist.t_on, "type": "fault_on", "element_id": eid},
                {
                    "time": dist.t_off,
                    "type": "fault_off",
                    "element_id": eid,
                },
            ],
        },
    )
    solver["system_model_file"] = "case.json"
    su.write_solver(solver, workdir / "solver.json")
    r = subprocess.run(
        [str(cfg.runner), "solver.json"],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=runner_env(),
    )
    log = r.stdout + r.stderr
    (workdir / "sim.log").write_text(log)
    return r.returncode, next(iter(workdir.glob("*.csv")), None), log


def _coi_weights(cfg: CampaignConfig, cols) -> np.ndarray:
    """Inertia weights for the centre-of-inertia angle, from device H*MVA."""
    w = []
    for c in cols:
        gid = mu.gen_id_of(c[len("Genrou_") : -len("_delta")], cfg.spec)
        w.append(cfg.hmva.get(gid, 0.0))
    w = np.asarray(w, dtype=float)
    return w / w.sum() if w.sum() > 0 else np.full(len(w), 1.0 / max(len(w), 1))


def is_unstable(cfg: CampaignConfig, csv, pole_slip_deg: float = 180.0) -> bool:
    """Transient instability: a machine pole-slips relative to the COI.

    Also treats a missing output, NaN, or a dead run as unstable, because the
    integrator giving up is itself evidence the point was not survivable.

    Deliberately NOT used: "final spread is close to max spread". That was
    tried and gives false positives, because the post-fault angle spread
    settles to a nonzero value, so a mild fault looks like sustained swinging.
    """
    if csv is None or not Path(csv).exists():
        return True
    df = pd.read_csv(csv)
    df.columns = [c.strip() for c in df.columns]
    df = df.drop_duplicates(subset="t", keep="last")
    cols = [c for c in df.columns if c.endswith("_delta")]
    if not cols:
        return True
    dl = np.rad2deg(df[cols].to_numpy())
    if np.isnan(dl).any():
        return True
    coi = dl @ _coi_weights(cfg, cols)
    return bool(np.abs(dl - coi[:, None]).max() > pole_slip_deg)


def find_cct(
    cfg: CampaignConfig,
    workdir: Path,
    dist: Disturbance,
    lo: float = 0.05,
    hi: float | None = None,
    tol: float = 0.01,
):
    """Critical clearing time by bisection. Returns (cct_s, censored, n_sims).

    `censored` means the system was still stable at `hi`, so `cct_s` is a
    LOWER BOUND, not a point estimate, and must be treated as right-censored
    when fitting a distribution.

    `hi` defaults to `cfg.cct_max`, deliberately capped at a credible
    protection time rather than extended until instability appears. EHV primary
    clearing is ~0.05-0.10 s and breaker-failure backup ~0.25-0.50 s; a
    multi-second fault is not a contingency any protection scheme permits, so a
    CCT beyond the cap has no engineering meaning and is reported as censored.
    """
    hi = cfg.cct_max if hi is None else hi
    n = 0
    # Own subdirectory: the bisection would otherwise overwrite the base-case
    # monitor csv, leaving a probe run on disk where the recorded QoI came from.
    d = Path(workdir) / "cct"
    d.mkdir(parents=True, exist_ok=True)
    # A probe is asked one question by is_unstable: did any machine pole-slip?
    # That needs t + delta, 21 of 441 columns. Writing the 400 bus voltage
    # columns here was half the campaign's csv volume and was never read.
    probe = su.enable_monitors(
        mu.load_case(Path(workdir) / "case.json"), bus=(), genrou=("delta",)
    )
    mu.write_case(probe, d / "case.json")

    def unstable_at(dur):
        nonlocal n
        n += 1
        rc, csv, _ = simulate(cfg, d, replace(dist, duration=dur))
        return rc != 0 or is_unstable(cfg, csv)

    if unstable_at(lo):
        return float("nan"), False, n  # unstable even at the shortest fault
    if not unstable_at(hi):
        return float(hi), True, n  # still stable at the cap: censored
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if unstable_at(mid):
            hi = mid
        else:
            lo = mid
    return round(float(lo), 4), False, n


def write_timeseries(csv, out: Path, float32: bool = True) -> int:
    """Monitor csv -> parquet. Returns rows written.

    Drops the duplicate timestamps GridKit emits at each event here, once, so
    no downstream consumer has to rediscover that `np.gradient` returns inf
    across a zero dt. `t` stays float64; the channels do not need it.
    """
    df = pd.read_csv(csv)
    df.columns = [c.strip() for c in df.columns]
    df = df.drop_duplicates(subset="t", keep="last")
    if float32:
        df = df.astype({c: "float32" for c in df.columns if c != "t"})
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, compression="zstd", index=False)
    return len(df)


def _jsonable(v):
    if isinstance(v, (np.floating, np.integer, np.bool_)):
        v = v.item()
    if isinstance(v, float) and not math.isfinite(v):
        return None  # JSON has no NaN literal; null is the honest encoding
    return v


def persist(cfg: CampaignConfig, tag: str, rec: dict, csv=None) -> dict:
    """Write the durable artifacts for one condition to `cfg.store`.

    Everything else a run produces lives on scratch and dies with the job.
    Written straight to /projects rather than staged and copied, so there is no
    migration step to forget before the scratch purge.

    A failed condition is persisted too: an infeasible point is a result.
    """
    if cfg.keep_timeseries and csv is not None and Path(csv).exists():
        series = cfg.store / "series" / f"{tag}.parquet"
        rec["series_rows"] = write_timeseries(csv, series)
        rec["series"] = str(series)

    record = cfg.store / "records" / f"{tag}.json"
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_text(json.dumps({k: _jsonable(v) for k, v in rec.items()}, indent=1))
    return rec


def add_datacenter(
    cfg: CampaignConfig, frames: dict, bus: int, mw: float, power_factor: float = 0.98
) -> dict:
    """Add datacenter load at a bus, before dispatch. Returns modified frames.

    Added to the `.m` so the OPF has to serve it. The commitment is not resized
    here: `commit_n` runs afterwards with the same `n`, so a dispatch the
    committed fleet cannot serve is a failure, not a re-commit.
    `power_factor` 0.98 lagging is an assumption, not a sourced value.
    """
    b = frames["bus"]
    row = b["BUS_I"].astype(int) == int(bus)
    if not row.any():
        raise KeyError(f"bus {bus} not in case")
    mvar = mw * np.tan(np.arccos(power_factor))
    b.loc[row, "PD"] = b.loc[row, "PD"] + mw
    b.loc[row, "QD"] = b.loc[row, "QD"] + mvar
    return frames


def split_out_datacenter(
    cfg: CampaignConfig,
    case: dict,
    bus: int,
    mw: float,
    power_factor: float = 0.98,
    alpha_p: float = 1.0,
    alpha_i: float = 0.0,
):
    """Separate the datacenter from the bus load so it can have its own ZIP mix.

    `build_case_from_m` folds all demand at a bus into one LoadZIP, which would
    force the datacenter to share the surrounding load's constant-impedance
    behaviour. A real datacenter is rectifier-fed and closer to constant power,
    which does NOT self-shed during a dip, so it must be a separate device.
    A bus with no load in the seed has no LoadZIP in the case template, even
    though `add_datacenter` put the MW into the dispatch. There is then nothing
    to split: the new device carries the whole datacenter, which is what keeps
    the dynamic model's load equal to the dispatched load. Greenfield siting at
    an unloaded bus is realistic, so this is a supported case, not an error.
    """
    case = copy.deepcopy(case)
    mvar = mw * np.tan(np.arccos(power_factor))
    host = [
        d
        for d in case["devices"]
        if d["class"] == "LoadZIP"
        and d["ports"]["bus"] == int(bus)
        and not d["id"].endswith("_ibr_stand_in")
    ]
    if host:
        h = host[0]
        h["params"]["Pnom"] -= mw / BASEMVA
        h["params"]["Qnom"] -= mvar / BASEMVA
    case["devices"].append(
        {
            "class": "LoadZIP",
            "ports": {"bus": int(bus)},
            "id": f"datacenter_{bus}",
            "params": {
                "Pnom": mw / BASEMVA,
                "Qnom": mvar / BASEMVA,
                "alphaI": alpha_i,
                "alphaP": alpha_p,
            },
        }
    )
    return case


def run_condition(
    cfg: CampaignConfig,
    tag: str,
    demand: dict,
    wind: dict,
    n_commit: int,
    with_cct: bool = True,
    dc_bus: int | None = None,
    dc_mw: float = 0.0,
    dist: Disturbance | None = None,
):
    """One operating condition end to end. Returns a flat record.

    Thin wrapper over `solve_point` + `run_disturbance`, kept because the phase
    scripts use it. A sweep that fans several disturbances over one solved
    point should call the two stages itself rather than re-solving each time.
    """
    t0 = time.time()
    d = cfg.work / tag
    dist = cfg.disturbance if dist is None else dist
    case, rec = solve_point(cfg, demand, wind, n_commit, d, dc_bus, dc_mw)
    rec["tag"] = tag
    if case is None:
        rec["wall_s"] = round(time.time() - t0, 1)
        return persist(cfg, tag, rec)

    rec.update(run_disturbance(cfg, case, dist, d, with_cct=with_cct))
    rec["wall_s"] = round(time.time() - t0, 1)
    return persist(cfg, tag, rec, rec.pop("_csv", None))


def apply_disturbance(cfg: CampaignConfig, case: dict, dist: Disturbance) -> dict:
    """Stage B: reshape case.json without touching the dispatch.

    Deep-copies, because the point of the split is to reuse one solved case
    across many disturbances; mutating it in place would leak H scaling and
    monitor choices from one sample into the next.
    """
    case = copy.deepcopy(case)
    case = su.set_fault_bus(case, dist.fault_bus)
    case = su.enable_monitors(case, bus=dist.bus_mon, genrou=dist.genrou_mon)
    if dist.h_scale != 1.0:
        for dev in case["devices"]:
            if dev["class"] == "Genrou":
                dev["params"]["H"] *= dist.h_scale
    if dist.alpha_p is not None:
        # Original bus loads only: the wind stand-in needs alphaI=1 and the
        # datacenter alphaP=1, and overwriting either changes a different model.
        for dev in case["devices"]:
            if (
                dev["class"] == "LoadZIP"
                and not dev["id"].endswith("_ibr_stand_in")
                and not dev["id"].startswith("datacenter_")
            ):
                dev["params"]["alphaP"] = dist.alpha_p
                dev["params"]["alphaI"] = 0.0
    return case


def solve_point(
    cfg: CampaignConfig,
    demand: dict,
    wind: dict,
    n_commit: int,
    workdir: Path,
    dc_bus: int | None = None,
    dc_mw: float = 0.0,
):
    """Stage A: hour and commitment -> solved case, before any disturbance.

    Returns `(case, rec)`. `case` is None when the dispatch did not solve, and
    the caller must not run dynamics on it.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    frames = apply_hour(cfg, demand, wind)
    if dc_bus is not None and dc_mw > 0:
        frames = add_datacenter(cfg, frames, dc_bus, dc_mw)
    net = frames["bus"]["PD"].sum() - frames["gen"].loc[cfg.nd, "PG"].sum()
    frames = commit_n(cfg, frames, n_commit, merit_order(cfg))
    rec = {
        "demand_mw": round(float(frames["bus"]["PD"].sum()), 1),
        "wind_mw": round(float(frames["gen"].loc[cfg.nd, "PG"].sum()), 1),
        "net_mw": round(float(net), 1),
        "n_commit": n_commit,
        "dc_bus": dc_bus,
        "dc_mw": dc_mw,
    }

    solved, info = solve_chain(cfg, frames, workdir)
    rec.update(info)
    if solved is None:
        rec["status"] = "solve_failed"
        return None, rec

    case, _ = mu.build_case_from_m(
        cfg.spec, cfg.base_case, workdir / "ac.m", dyr=cfg.dyr
    )
    case, n_wind = wind_to_negative_load(cfg, case, solved)
    if dc_bus is not None and dc_mw > 0:
        case = split_out_datacenter(
            cfg,
            case,
            dc_bus,
            dc_mw,
            alpha_p=cfg.dc_alpha_p,
            alpha_i=1.0 - cfg.dc_alpha_p,
        )
    rec["wind_as_load"] = n_wind
    rec["dc_alpha_p"] = cfg.dc_alpha_p
    return case, rec


def solve_points_batch(cfg: CampaignConfig, specs: list) -> list:
    """Stage A for many operating points, sharing one Julia process per stage.

    `specs` is a list of dicts with keys `demand`, `wind`, `n_commit`,
    `workdir`, and optionally `dc_bus`, `dc_mw`. Returns `(case | None, rec)`
    per spec, in order, matching `solve_point`.
    """
    prepared, items = [], []
    order = merit_order(cfg)
    for s in specs:
        frames = apply_hour(cfg, s["demand"], s["wind"])
        dc_bus, dc_mw = s.get("dc_bus"), float(s.get("dc_mw") or 0.0)
        if dc_bus is not None and dc_mw > 0:
            frames = add_datacenter(cfg, frames, int(dc_bus), dc_mw)
        net = frames["bus"]["PD"].sum() - frames["gen"].loc[cfg.nd, "PG"].sum()
        frames = commit_n(cfg, frames, int(s["n_commit"]), order)
        prepared.append((frames, net, dc_bus, dc_mw))
        items.append((frames, Path(s["workdir"])))

    results = solve_chain_batch(cfg, items)

    out = []
    for (frames, net, dc_bus, dc_mw), (solved, info), s in zip(
        prepared, results, specs
    ):
        rec = {
            "demand_mw": round(float(frames["bus"]["PD"].sum()), 1),
            "wind_mw": round(float(frames["gen"].loc[cfg.nd, "PG"].sum()), 1),
            "net_mw": round(float(net), 1),
            "n_commit": int(s["n_commit"]),
            "dc_bus": dc_bus,
            "dc_mw": dc_mw,
        }
        rec.update(info)
        if solved is None:
            rec["status"] = "solve_failed"
            out.append((None, rec))
            continue
        wd = Path(s["workdir"])
        case, _ = mu.build_case_from_m(
            cfg.spec, cfg.base_case, wd / "ac.m", dyr=cfg.dyr
        )
        case, n_wind = wind_to_negative_load(cfg, case, solved)
        if dc_bus is not None and dc_mw > 0:
            case = split_out_datacenter(
                cfg,
                case,
                int(dc_bus),
                dc_mw,
                alpha_p=cfg.dc_alpha_p,
                alpha_i=1.0 - cfg.dc_alpha_p,
            )
        rec["wind_as_load"] = n_wind
        rec["dc_alpha_p"] = cfg.dc_alpha_p
        out.append((case, rec))
    return out


def candidate_buses(
    cfg: CampaignConfig,
    n: int = 16,
    from_bus: int = 66,
    require_load: bool = False,
    keep: tuple = (),
    seed: int = 0,
) -> pd.DataFrame:
    """Sample `n` buses stratified by hop distance from `from_bus`.

    Returns bus, hops, branches, sum_rate (MVA) and PD, which is what a siting
    or contingency design needs to be reported as stratified rather than
    arbitrary.

    `require_load=True` keeps only buses that already carry load: a datacenter
    is a constant-power device, and alone at an unloaded bus it prevents the
    integrator from initialising (see plan doc, section 7). `keep` forces buses
    in regardless, which is how deliberately weak buses stay in the sample --
    their infeasibility is a result, not a failure.
    """
    br, bus = cfg.seed["branch"], cfg.seed["bus"]
    adj: dict = {}
    nbr: dict = {}
    rate: dict = {}
    for f, t, r in zip(
        br.F_BUS.astype(int), br.T_BUS.astype(int), br.RATE_A.astype(float)
    ):
        adj.setdefault(f, []).append(t)
        adj.setdefault(t, []).append(f)
        nbr[f] = nbr.get(f, 0) + 1
        nbr[t] = nbr.get(t, 0) + 1
        rate[f] = rate.get(f, 0.0) + r
        rate[t] = rate.get(t, 0.0) + r

    hops = {int(from_bus): 0}
    queue = [int(from_bus)]
    while queue:
        u = queue.pop(0)
        for v in adj.get(u, []):
            if v not in hops:
                hops[v] = hops[u] + 1
                queue.append(v)

    pdm = dict(zip(bus.BUS_I.astype(int), bus.PD.astype(float)))
    df = pd.DataFrame({"bus": sorted(hops)})
    df["hops"] = df.bus.map(hops)
    df["branches"] = df.bus.map(nbr).fillna(0).astype(int)
    df["sum_rate"] = df.bus.map(rate).round(1)
    df["PD"] = df.bus.map(pdm).round(2)

    pool = df[df.bus != int(from_bus)]
    if require_load:
        pool = pool[pool.PD > 0.0]
    keep = [int(b) for b in keep]
    rng = np.random.default_rng(seed)
    sel = list(keep)
    rings = {h: g for h, g in pool[~pool.bus.isin(keep)].groupby("hops")}
    ladder = sorted(rings)
    i = 0
    while len(sel) < n and ladder and i < 1000:
        cand = [b for b in rings[ladder[i % len(ladder)]].bus.tolist() if b not in sel]
        if cand:
            sel.append(int(rng.choice(cand)))
        i += 1
    return df[df.bus.isin(sel)].sort_values("hops").reset_index(drop=True)


def run_disturbance(
    cfg: CampaignConfig,
    case: dict,
    dist: Disturbance,
    workdir: Path,
    with_cct: bool = True,
) -> dict:
    """Stages B and C: apply a disturbance to a solved case and measure it.

    Cheap relative to `solve_point`, which is the whole reason for the split.
    The monitor csv is returned under `_csv` for the caller to persist.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    case = apply_disturbance(cfg, case, dist)
    rec = {
        "fault_bus": dist.fault_bus,
        "fault_duration_s": dist.duration,
        "h_scale": dist.h_scale,
        "alpha_p": dist.alpha_p,
        "n_genrou": sum(1 for x in case["devices"] if x["class"] == "Genrou"),
        "inertia_MWs": round(cfg.inertia_of(case)),
    }

    for f in workdir.glob("*.csv"):
        f.unlink()
    mu.write_case(case, workdir / "case.json")
    rc, csv, log = simulate(cfg, workdir, dist)
    rec["sim_rc"] = rc
    rec["sim_errors"] = sum("ERROR" in ln for ln in log.splitlines())
    if rc == 0 and csv is not None:
        rec.update(reduce_qoi(cfg, csv, dist))
        rec["status"] = "ok"
        rec["_csv"] = csv
    else:
        rec["status"] = "sim_failed"

    if with_cct and rec["status"] == "ok":
        cct, censored, nsim = find_cct(cfg, workdir, dist)
        rec["cct_s"] = cct
        rec["cct_censored"] = censored
        rec["cct_sims"] = nsim
    return rec
