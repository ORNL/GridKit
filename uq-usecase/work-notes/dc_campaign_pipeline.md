# Datacenter UQ campaign: code pipeline

**Audience:** someone (human or agent) who must modify or extend this code
without having run the study. It maps every code path from raw case files to a
QoI record.

| document | answers | read it when |
|---|---|---|
| [datacenter_pilot_report.md](datacenter_pilot_report.md) | *what did we find, and how, in general terms* | you want the study. Stands alone |
| [datacenter_load_uq_plan.md](datacenter_load_uq_plan.md) | *why each choice was made; every number and its provenance* | you want to question a decision |
| **this file** | *which function does what* | you are changing the code |

This file is the *what* and *where*; the plan is the *why*. If the two
disagree about a settled parameter, the
[settled design table](datacenter_load_uq_plan.md#settled-design-authoritative)
in the plan wins.

**Repo root:** `~/gridkit/uq-usecase`
**Python:** `/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python`
**Store:** `/kfs2/projects/scidac/scidac-data/gridkit-campaigns/<name>/`

---

## 0. Places

Where everything lives. Paths are absolute because workers run from varying
directories.

### Code

| what | where |
|---|---|
| this project | `~/gridkit/uq-usecase/` |
| pipeline library | `~/gridkit/uq-usecase/py-utils/` |
| Julia PF/OPF solver | `~/gridkit/uq-usecase/pm-solver/pm_solve.jl` (+ `Project.toml`) |
| campaign + analysis scripts | `~/gridkit/uq-usecase/scripts/` |
| notes, plan, report | `~/gridkit/uq-usecase/work-notes/` |
| figures | `~/gridkit/uq-usecase/figs/` |
| GridKit source | `~/gridkit/` |
| **the runner binary** | `~/gridkit/build/application/PhasorDynamics/DynamicSimulation` |
| shipped cases | `~/gridkit/cases/PhasorDynamics/<CASE>/<CASE>.case.json` |
| validation examples | `~/gridkit/examples/PhasorDynamics/Validation/<CASE>/` |

### Input data (read-only)

| what | where |
|---|---|
| TAMU raw data root | `/kfs2/projects/scidac/scidac-data/` |
| ACTIVSg200 raw (`.aux`, `.dyr`, `.m`) | `.../ACTIVSg200/raw-tamu-data/` |
| ACTIVSg200 profiles | `.../ACTIVSg200/grid-v1/{demand,wind}.csv` |
| ACTIVSg2000 raw | `.../ACTIVSg2000/raw-tamu-data/` |
| ACTIVSg2000 profiles | `.../ACTIVSg2000/grid-v1/{demand,wind,solar}.csv` |
| ACTIVSg2000 template `.m` | `.../ACTIVSg2000/grid-v1/case_ACTIVSg2000.m` |

### Experiment output

| what | where | lifetime |
|---|---|---|
| **durable results** | `/kfs2/projects/scidac/scidac-data/gridkit-campaigns/<name>/` | kept, **not backed up** |
| per-condition records | `<campaign>/records/<tag>.json` | source of truth |
| monitor timeseries | `<campaign>/series/<tag>.parquet` | |
| collected tables | `<campaign>/{results.parquet,results.csv,frontier.csv}` | |
| sbatch logs | `<campaign>/slurm/` | |
| **per-run churn** | `/scratch/$USER/gridkit-campaigns/<op tag>/` | **purged after 28 days idle** |
| scaling probes | `/scratch/$USER/scale-probe/` | scratch |

Campaigns run so far: `phase3b`, `aleatoric1`, `faultscreen1`, `afront1`,
`afront2`, `afront3`, `cct_check`, `cct181`, `faultscreen1_ap1`,
`aleatoric1_ap1`, `faultscreen1_ap08`, `aleatoric1_ap08`.

**Never put campaign data in `/tmp`** on these nodes: it is tmpfs, i.e. RAM.

### Environment

```bash
PY=/home/isatkaus/flash-w-dom/conda-envs/h-py312-basic/bin/python
module load julia/1.12.6          # `module avail julia` lists only 1.12.6 as a module; pm_utils.resolve_julia() uses the 1.12.1 binary under /nopt/nrel/apps/cpu_stack/software/julia; PATH may pick ~/bin/julia (1.12.5); ~/bin/julia112 is broken, do not use
sbatch --account=scidac --partition=short
```

`campaign.runner_env()` scopes GridKit's `libomp` to the simulation subprocess.
Do **not** export it globally: that makes Julia load the wrong SuiteSparse.

---

## 1. Module map

| file | lines | role |
|---|---|---|
| `py-utils/campaign.py` | 1277 | **the pipeline.** Config, stages, QoI, CCT, storage, PV→PQ Q-limit switching (§4.2a) |
| `py-utils/m_to_case_utils.py` | 868 | `.aux`/`.dyr`/`.m` → `case.json`; device construction |
| `py-utils/aux_to_m.py` | 381 | `.aux` ↔ MATPOWER `.m` frames |
| `py-utils/perturb.py` | 327 | frame-level edits: load scaling, commitment, gen output |
| `py-utils/pm_utils.py` | 400 | PowerModels.jl driver (incl. batch mode) |
| `py-utils/solver_utils.py` | 351 | `solver.json`: faults, monitors, events |
| `py-utils/slurm_utils.py` | 1058 | SLURM submit/monitor. **Only the submit half is used** |
| `pm-solver/pm_solve.jl` | ~480 | Julia: DC OPF / AC PF, single and `--batch` |
| `scripts/sweep_build.py` | — | design → `conditions.csv` + sbatch |
| `scripts/sweep_worker.py` | 181 | run a slice of conditions; `--dc-alpha-p` sets the datacenter `alphaP` (default 0.8); `--no-enforce-q-lims` turns PV→PQ switching off (on by default, §4.2a); `--enforce-q-lims` is still accepted; `--ac-opf` replaces DC OPF and AC PF with a PM.jl AC OPF (§4.2b); `--ac-opf-relax-q` is that route's no-Q-limit control |
| `scripts/highp_rerun.sh` | 36 | SLURM: re-run the sited arm of a campaign at another `alphaP` (`sbatch highp_rerun.sh <src> <dst> <alphaP> <workers>`); baseline rows dropped. **Not valid for Q-limited campaigns**, where the baseline changes too |
| `scripts/qlim_campaign.sh` | — | SLURM: full design (baseline and sited rows) at `alphaP` with `--enforce-q-lims`. `sbatch --account=scidac --partition=short qlim_campaign.sh <src> <dst> <alphaP> <workers>`. Used for `faultscreen1_ap08_qlim` (job 19797428) and `faultscreen1_qlim` (job 19797515). The fifth argument is the method: `switch` (default), `acopf` (§4.2b), or `acopf_noq` (control). Used for `faultscreen1_{ap08_,}acopf` (jobs 19797889/890, rerun 19797987/988) and `faultscreen1_{ap08_,}acopfnoq` (jobs 19798072/073, rerun 19798109/110) |
| `scripts/qlim_pf_check.py` | — | Five-point check of §4.2a: unlimited vs switched AC PF (Q violations, units fixed, vmin, losses). Writes under `/scratch/.../qlim_check/` |
| `scripts/qlim_dyn_check.py` | — | Pre-fault drift check on one sited point with the default disturbance, unlimited vs switched |
| `scripts/acopf_pilot.py` | — | Five-point check of §4.2b: AC OPF vs switched PF (losses, dispatch change, slack at PMAX, Vmax, warm solve time, pre-fault drift). `--design --points --dc-alpha-p --dyn`. Writes under `/scratch/.../qlim_check/acopf_pilot/` (copy at `/kfs2/.../gridkit-campaigns/qlim_check/acopf_pilot/`) |
| `scripts/qlim_compare.py` | — | Paired CCT comparison of two campaigns by tag: censoring, mean/max ΔCCT, uncensored-in-both subset. `python qlim_compare.py <first> <second>`; the output labels the first arm "unlimited" whatever it is, so read the directories |
| `scripts/sweep_collect.py` | — | records → `results.parquet` |
| `scripts/analyse_censored.py` | — | censored-aware stats + figures |
| `scripts/alphap_frontier.py` | — | feasible-`alphaP` surface over buses / sizes / faults |
| `scripts/alphap_collect.py` | — | frontier tables + failure-mode breakdown |
| `scripts/check_zip_composition.py` | — | verifies co-located `LoadZIP`s compose correctly |
| `scripts/cct_check.py` | 121 | CCT at `alphaP` 0.2 vs 1.0 over the `alphap_frontier` design (`--name --fault-buses --alphas --start --count`; `with_cct=True`) |
| `scripts/afront3.sh` | 27 | SLURM wrapper: `alphap_frontier.py` over faults 94, 68, 181 (campaign `afront3`) |
| `scripts/cctcheck.sh` | 31 | SLURM wrapper: `sbatch cctcheck.sh <name> <faults>`; 16 workers over the 272-dispatch design (campaigns `cct_check`, `cct181`) |
| `scripts/ida_zoom.py` | 252 | IDA initial-condition zoom-in on one failing case (`--stage bus\|solver\|mw`, `--src --name --alphas --buses --fault-bus --workers`); writes `/scratch/.../gridkit-campaigns/ida_zoom/<name>/results.csv` |
| `scripts/ida_discriminate.py` | 82 | Event-schedule discriminator: imports `ida_zoom`, adds in-process modes `open_ev_t02` (fault_on at 1.0, tmax 0.2) and `open_noev_t12` (no event, tmax 1.2), runs `none/open/open_ev_t02/open_noev_t12/event/closed` over `--buses` at `--alpha`; writes `/scratch/.../gridkit-campaigns/ida_zoom/<name>/results.csv`. Does not modify `ida_zoom.py`. Submitted with `sbatch --wrap` |
| `scripts/pm_fault_check.py` | 187 | PowerModels.jl AC PF over fault bus x datacenter MW (`pm_utils.run_pm_solve_batch`); builds one `.m` per case from the base `ac.m`; prints PF convergence next to IDA IC status |
| `scripts/ida_xsweep.py` | 73 | IDA fault-strength sweep (test 3, `ida_ic_issue.md` §7): `ida_zoom.run_one` with `fault_x` over X in {1.0 … 0.01} at buses 1, 148, 94, `alphaP` 1.0, closed mode. Output `ida_zoom/<name>/results.csv` (adds `fault_x`). Run `python ida_xsweep.py --name xsweep181 --workers 16` |
| `scripts/ida_event_repro.py` | 110 | Self-contained event-time reproducer (test 4): writes case dirs (`case.json`, `solver.json`, `sim.log`) for event times {0.5, 1.0, 2.0} at buses 1 and 94, plus no-event controls. Output `ida_zoom/event_repro/results.csv`. Uses the plain campaign runner, not `ida_zoom`'s global `FAULT_ON_T` |
| `scripts/ida_pvpq.jl` | 142 | PowerModels.jl PF with PV→PQ Q-limit switching (test 2). Separate from `pm-solver/pm_solve.jl` (not modified). Switched gens are deactivated and replaced by a load of (−Pg, −Qlim); one-way switch, slack never switched. Run as `julia --project=pm-solver ida_pvpq.jl manifest.json results.json` |
| `scripts/ida_pvpq.py` | 114 | Driver for `ida_pvpq.jl` (test 2): builds cases with `pm_fault_check.build_case` (warm and flat starts, buses 1, 148, 94, 181, none × `MW_GRID`), calls `pm_utils.resolve_julia()`, prints P/S/X tables next to IDA IC status. Output `ida_zoom/pvpq/results.csv` |
| `scripts/ida_mechanism.py` | 87 | Mechanism tests on `h1879` (`ida_ic_issue.md` §7 item 6), event mode, 913 runs via `ida_zoom.run_one`: `nodc` (`dc_zero`, all 200 buses), `mw` (MW ∈ {25…200} × alphaP ∈ {0.2, 0.8, 1.0} for the 23 buses failing at 1.0), `xsweep` (X ∈ {0.5…0.01} × alphaP ∈ {0.2, 0.8}). Output `ida_zoom/mechanism/results.csv`. Run `python ida_mechanism.py --workers 16` (SLURM ~1 min) |
| `scripts/ida_seed_test.py` | 122 | IDA seed test (`ida_ic_issue.md` §7 item 7): solves fault-on AC PF with `pm_fault_check.build_case` and `pm_utils.run_pm_solve_batch` (PM.jl) at buses 64, 146, 1, 94 and datacenter 0 / 200 MW; writes the solved bus voltages into a copy of `case.json` as the bus `init` (the IDA t = 0 guess, `Bus::initialize`); runs `ida_zoom.run_one` in `closed` mode with the original and seeded guess. Output `ida_zoom/seed_test/results.csv`. Run `python ida_seed_test.py` on the current exclusive compute node (check with `hostname`; it changes from day to day), never the login node (~1 min) |
| `scripts/failure_audit.py` | 111 | Tally of `solve_failed` by stage and reason, and of `sim_failed` by `sim_rc` and class (IDACalcIC vs IDASolve, read from scratch `sim.log`). `python failure_audit.py <campaign> [--ref <campaign>] [--out <file>]`; `--ref` splits the sim-failed overlap. Issue doc §10.11 |
| `scripts/dispatch_margin.py` | 63 | Committed thermal PMAX minus `net_mw` per dispatch (net_mw includes the datacenter MW), to test DC OPF infeasibility. `python dispatch_margin.py <campaign>` |
| `scripts/dispatch_pf_probe.py` | 82 | Solves four `solve_failed` operating points from `faultscreen1` with Q limits off and on, logging each Julia batch pass. Writes only to scratch `dispatch_probe/`. `python dispatch_pf_probe.py` |
| `scripts/fold_check_h1879.py` | 72 | Checks h1879 sim-failed rows against the fixed-state fold sets (`ida_zoom/fixed_state/alpha_sweep.csv`, IDACalcIC at 200 MW). Assumes the sweep used the same datacenter bus (not checked). `python fold_check_h1879.py` |
| `scripts/check_docs.py` | — | validates work-notes tables, fences, anchors, links |

`slurm_utils.generate_sim_script` is **not** used: it xargs `DynamicSimulation`
over pre-built directories with no Python in the loop, which does not fit this
pipeline. Only `compute_slices`, `submit_jobs_batch` and the `sacct` helpers are
reused.

---

## 2. End-to-end flow

```
hour (load, wind)  ──┐
commitment n       ──┼─► apply_hour ─► add_datacenter ─► commit_n ─► frames (dict of DataFrames)
datacenter bus/MW  ──┘                                                   │
                                                                         ▼
                                            solve_chain_batch  (Julia, batched)
                                              DC OPF ─► AC PF ─► loss loop
                                                                         │
                                                      solved frames + ac.m
                                                                         ▼
                                        build_case_from_m ─► case.json skeleton
                                        wind_to_negative_load
                                        split_out_datacenter
                                                                         │
                          ═════════ end of STAGE A (expensive, ~2.5 s batched) ═════════
                                                                         ▼
                                        apply_disturbance  (deepcopy!)
                                          set_fault_bus, enable_monitors,
                                          h_scale, alpha_p
                                                                         ▼
                                        simulate ─► DynamicSimulation ─► *.csv
                                                                         ▼
                                        reduce_qoi ─► scalars
                                        find_cct   ─► bisection on duration
                                                                         ▼
                                        persist ─► record.json + mon.parquet
                          ═════════ STAGE B/C (cheap, ~0.4–5 s) ═════════
```

**The stage boundary is the single most important fact in this codebase.**
Stage A costs a Julia solve; stage B/C is a JSON edit and a re-run. Measured:
Julia launch 5.9 s, warm solve 0.04 s, one dynamic run 0.4–0.8 s (full monitors,
`tmax` 10 s; 0.28 s with `tmax` 5 s).

---

## 3. Data structures

### `frames` — MATPOWER tables
`dict` of DataFrames: `bus`, `branch`, `gen`, `gencost`, `genfuel`. Column names
are MATPOWER (`BUS_I`, `PD`, `QD`, `GEN_BUS`, `PG`, `PMAX`, `PMIN`,
`GEN_STATUS`, `F_BUS`, `T_BUS`, `BR_R`, `BR_X`, `RATE_A`). Produced by
`aux_to_m.m_to_frames(path)`, written by `aux_to_m.write_matpower_m`.

**Always copy before editing:** `perturb.copy_frames(frames)`.

### `case` — GridKit `case.json`
```python
{"params": {"freq_base": 60.0, ...},
 "buses":   [{"id": .., "mon": ["Vr","Vi"], ...}],
 "devices": [{"class": "Genrou"|"LoadZIP"|"BusFault"|"SexsPti"|"Tgov1",
              "id": "...", "ports": {"bus": 129}, "params": {...},
              "mon": ["delta","omega"]}]}
```
Load devices: `LoadZIP` with `Pnom`, `Qnom` **in pu on 100 MVA**, plus `alphaI`,
`alphaP`; `alphaZ = 1 - alphaI - alphaP` is implicit and **not settable**.

Device id conventions that code depends on:

| suffix / prefix | meaning |
|---|---|
| `<gid>_genrou` | synchronous machine |
| `<gid>_ibr_stand_in` | wind, converted to negative `LoadZIP` |
| `datacenter_<bus>` | the added datacenter load |

`apply_disturbance` **skips both of the latter** when applying a global
`alpha_p`, because the wind stand-in needs `alphaI=1` and the datacenter carries
its own mix.

### `rec` — the record
Flat `dict`, one per condition, serialised to `store/records/<tag>.json`.
Non-finite floats become `null` (JSON has no `NaN`). Keys include `tag`,
`status`, `net_mw`, `n_commit`, `dc_bus`, `dc_mw`, `fault_bus`, `inertia_MWs`,
`cct_s`, `cct_censored`, the `reduce_qoi` scalars, and on failure `fail_stage`
and `fail_reason`.

`status` is one of `ok`, `solve_failed`, `sim_failed`.

---

## 4. Stage reference

### 4.1 Case construction (once, in `CampaignConfig.__init__`)

| step | function |
|---|---|
| parse PowerWorld `.aux` and `.dyr` | `mu.parse_aux`, `mu.parse_dyr` |
| build reference case from source | `mu.build_case_from_aux_dyr` |
| diff against the shipped `case.json` | `mu.make_creator_overlay` |
| rebuild with the overlay applied | `mu.build_case_from_aux_dyr(..., creator_overlay=)` |
| write seed `.m` and load frames | `aux_to_m.write_matpower_m`, `m_to_frames` |

The overlay exists because the shipped case differs from a clean rebuild in one
machine. Verified: 0 diffs with overlay `{197_1_genrou: 14}`.

`seed.m` is written **per process** (`seed.<pid>.tmp` then `os.replace`) because
workers share the campaign scratch on Lustre and a plain write interleaves.

### 4.2 Operating point (stage A)

| function | does |
|---|---|
| `apply_hour(cfg, demand, wind)` | per-bus loads for the hour; pins wind `PMIN=PMAX=PG` |
| `add_datacenter(cfg, frames, bus, mw, pf=0.98)` | adds MW/MVAr to `bus` **before** dispatch, so the OPF must serve it |
| `merit_order(cfg)` | nuclear must-run first, then by `c1` |
| `commitment_band(cfg, net_mw, order)` | `(n_min, n_max)`; necessary condition only |
| `commit_n(cfg, frames, n, order)` | commits the `n` cheapest; non-dispatchables stay on |
| `solve_chain(cfg, frames, workdir)` | DC OPF → AC PF with loss loop, **one condition**; a one-item call to `solve_chain_batch` |
| `solve_chain_batch(cfg, items)` | same, **many conditions in lockstep**, 2·(loss_iters+1) Julia launches total; with `enforce_q_lims` each AC step is the switching loop of §4.2a |
| `solve_point` / `solve_points_batch` | the above plus case construction; return `(case, rec)` |

`solve_*` returns `case is None` on dispatch failure — **callers must check**,
and a failed case writes no `ac.m`, so nothing stale can be picked up.

### 4.2a Reactive limits: PV-to-PQ switching

On by default: `CampaignConfig(enforce_q_lims=None)` resolves to on unless `ac_opf`. Pass `enforce_q_lims=False` or
`sweep_worker.py --no-enforce-q-lims` to reproduce the unlimited campaigns exactly.
The physics, the OPF/PF split and the measured effect are in
[plan §10](datacenter_load_uq_plan.md#reactive-limits-in-the-ac-pf-pv-to-pq-switching)
and [pilot §4.2a](datacenter_pilot_report.md#42a-reactive-limits-pv-to-pq-switching).
This section is the code.

Module-level helpers (before `solve_chain`), constants `Q_TOL_MVAR` 0.1,
`V_TOL_PU` 1e-4, `Q_MAX_ITER` 30, `Q_MAX_RETURNS` 2:

| function | does |
|---|---|
| `_q_units(base)` | in-service, non-slack gens → `{gen row: (bus, QMAX, QMIN, VG)}` |
| `_q_switched(base, state, vm)` | copy of the frames with each switched unit moved into its bus: `PD -= PG`, `QD -= Qfix`, `GEN_STATUS = 0`. Bus becomes PQ (type 1) unless an in-service unit still holds it (type 2). `vm` warm-starts the voltages |
| `_q_step(base, solved, state, returns)` | updates `state` in place from one solved case. A violating unit is fixed at its limit (`("max", QMAX)` or `("min", QMIN)`). A fixed unit whose voltage passes its setpoint returns to PV, counted in `returns`, at most `Q_MAX_RETURNS` times. Returns the number of changes |
| `_q_restore(base, solved, state)` | the solved frames, with bus `PD`/`QD`/`BUS_TYPE` from `base`, and each switched unit back in service with `PG` from `base` and `QG = Qfix` |

The batch driver is the closure `ac_switched(it)` inside `solve_chain_batch`. It
replaces the single `run_stage("ac-pf", ...)` when the flag is on:

1. `ac_frames(k)`: original load with the DC OPF `PG`. Same start as the unswitched path.
2. Up to `Q_MAX_ITER` passes. Each pass runs `_q_switched`, writes `qin.m`, runs one
   batched `run_stage("ac-pf")` over the conditions still switching, reads `qsol.m`,
   then `_q_step`. A condition with no changes is done: `_q_restore` writes `ac.m`.
3. A condition still changing at the cap is marked dead, `fail_reason` `ac-pf-q`.

Per-condition state: `q_state[k]` (fixed units), `q_returns[k]`, `q_vm[k]` (warm
start), `q_iters[k]` (AC solves in the last loss pass). These feed the record fields
`q_lims`, `q_fixed` and `q_solves`. `run_stage` gained `skip=` so that settled and
dead conditions are not re-run.

`solve_chain` is now a one-item call to `solve_chain_batch`, so there is one
implementation.

**Invariant.** The dynamic build reads only the final `ac.m`
(`build_case_from_m`, §4.3). A switched unit is therefore never absent from `ac.m`:
it is restored in service at its fixed `QG`. `qin.m` and `qsol.m` are intermediate
(scratch) files only.

### 4.2b AC OPF path (PM.jl), second Q-limit route

Opt-in: `CampaignConfig(ac_opf=True)`, `sweep_worker.py --ac-opf`, or launcher
method `acopf`. It replaces the whole DC OPF plus AC PF chain with one PM.jl
`ac-opf` solve. The reasoning and the pilot are in
[plan §10](datacenter_load_uq_plan.md#ac-opf-as-a-second-q-limit-route-pmjl) and
[pilot §4.2b](datacenter_pilot_report.md#42b-ac-opf-as-a-second-q-limit-route).

Flow in `solve_chain_batch`, when `cfg.ac_opf` is set:

1. `prep_opf(k)` writes the original frames (load, commitment, no switching) to
   `opfin.m`. `run_stage("ac-opf", "ac-opf", prep_opf)` makes one batched
   `pm_solve.jl --pf-type ac-opf` launch. `solve_ac_opf` carries `QMIN`/`QMAX` as
   constraints. `write_solved_m` returns `ac.m` with only `Vm`/`Va` and `Pg`/`Qg`
   patched in.
2. `collect()` reads `ac.m` and fills `losses_mw` (in-service PG minus total PD),
   `max_over_pmax`, `vmin`, `vmax`, and `q_viol`.
3. The loss loop runs zero times (`range(0)`), because AC OPF includes losses.

`q_viol` counts in-service non-slack units outside their limits by more than
`Q_TOL_MVAR`, using the limits the solve used. In the relaxed control the solve used
±9999, so `collect()` scores against the original limits from the frames
(`q_src`). Records carry `ac_opf`, `ac_opf_relax_q`, and `q_lims`, which is true for
switched campaigns and for AC OPF without relaxation.

**Dispatch differs.** AC OPF re-dispatches P, unlike the DC OPF path. On the five
pilot points the summed change was about 280 MW, and one unit moved 53 MW. Rows with
no datacenter also shift. So an AC OPF campaign differs from a switched campaign in
dispatch as well as Q. The `acopf_noq` control isolates the Q part (§4.2b in the
pilot report).

**Pitfall.** Do not read `sacct` COMPLETED as success for this path (§6).

### 4.3 Case shaping

| function | does |
|---|---|
| `build_case_from_m(spec, base_case, ac.m, dyr=)` | solved `.m` → `case.json`. **Template-driven**: a bus with no load in `base_case` gets no `LoadZIP` |
| `wind_to_negative_load(cfg, case, frames)` | drops wind `Genrou` triplets, adds negative `LoadZIP` with `alphaI=1` at the same bus |
| `split_out_datacenter(cfg, case, bus, mw, alpha_p=, alpha_i=)` | reduces the host load, appends `datacenter_<bus>`. If no host exists, appends only — keeping dynamic load equal to dispatched load |

### 4.4 Disturbance (stage B/C)

| function | does |
|---|---|
| `apply_disturbance(cfg, case, dist)` | **deepcopies**, then sets fault bus, monitors, `h_scale`, `alpha_p` |
| `simulate(cfg, workdir, dist)` | writes `solver.json`, runs `DynamicSimulation`, returns `(rc, csv, log)` |
| `reduce_qoi(cfg, csv, dist)` | csv → 13 scalars |
| `is_unstable(cfg, csv)` | pole slip vs centre of inertia, 180° |
| `find_cct(cfg, workdir, dist)` | bisection on duration; probes use `delta`-only monitors |
| `run_disturbance(cfg, case, dist, workdir)` | the above, returns a record fragment |

`apply_disturbance` **must** deepcopy: the whole point of the split is reusing
one solved case across many disturbances.

### 4.5 Persist

| function | does |
|---|---|
| `write_timeseries(csv, out)` | → parquet (zstd, float32); drops duplicate `t` **once** |
| `persist(cfg, tag, rec, csv)` | `store/records/<tag>.json` + `store/series/<tag>.parquet` |
| `prepare_storage(scratch, store)` | creates trees, applies Lustre layouts |

Failed conditions are persisted too — infeasibility is a result, and a record on
disk is how resume works.

---

## 5. Configuration

`CampaignConfig(name=, scratch=, store=, ...)`. Constructor is cheap (0.2 s).

| knob | default | note |
|---|---|---|
| `dc_alpha_p` | `0.8` | datacenter constant-power share; `alphaI` takes the remainder, so the shipped mix is `P/I/Z = 80/20/0` (the first campaigns ran `20/80/0`, `alphaP` 0.2). `1.0` is the intended physics but has more IDA initial-condition folds (issue doc §10.2) |
| `cct_max` | `0.5` | censoring cap |
| `fault_bus`, `t_on`, `t_off`, `tmax` | 66, 1.0, 1.1, 10.0 | build the default `cfg.disturbance`; a sweep overrides per sample |
| `min_thermal` | `8` | commitment floor |
| `reserve` | `1.10` | adequacy margin |
| `loss_iters` | `1` | DC OPF loss iterations |
| `enforce_q_lims` | `None` | PV→PQ switching of generator Q limits in the AC PF (§4.2a). `None` means on unless `ac_opf`. `False` reproduces the unlimited campaigns |
| `ac_opf` | `False` | PM.jl AC OPF replaces DC OPF and AC PF; loss loop skipped (§4.2b). Mutually exclusive with `enforce_q_lims` (ValueError) |
| `ac_opf_relax_q` | `False` | with `ac_opf`: QMIN/QMAX widened to ±9999 on the OPF input only, as a control. Requires `ac_opf` (ValueError otherwise) |
| `keep_timeseries` | `True` | write `mon.parquet` |

`cfg.work` is a legacy alias for `cfg.scratch`.

---

## 6. Traps

Each of these has already caused a real failure here.

| trap | consequence | guard |
|---|---|---|
| **Fixed filename under shared scratch** | three separate races: `seed.m`, the Julia batch manifest, and per-condition workdirs that omitted `dc_bus`/`dc_mw` | PID-scope or include the full `OperatingPoint.tag` |
| **Solver log into a CSV field** | newlines corrupt every later row; cost 29 of 51 runs once | logs to file, records keep one-line `fail_reason` |
| **`dict.get(key, 0.0)` in a physics sum** | missing machines silently contributed zero inertia, 1.4 % low | read from the object |
| **Global `LD_LIBRARY_PATH` for the runner** | GridKit needs `libomp` from its build prefix, but putting it on the job env makes Julia load the wrong SuiteSparse | `campaign.runner_env()`, child process only |
| **`zip()` over solver results** | a short result list silently truncates and a case advances with no input | explicit length check |
| **Pandas `None` → `NaN`** | `is None` tests fail on `dc_bus` | `pd.isna` |
| **Aggregating over a swept dimension** | omitting `fault_bus` from a `groupby` key made `max()` pick the best fault location, so every frontier read `1.0` and the real effect vanished | put every swept dimension in the key |
| **Monitors default to almost nothing** | `"mon": []` on every bus | `su.enable_monitors` |
| **Event duplicate timestamps** | two rows at the same `t`; `np.gradient` → `inf` | `drop_duplicates(subset="t")`, done once at parquet write |
| **`dc_alpha_p` on `cct_check` records** | `cct_check.py` sets `alphaP` by editing the case dict, so the record field `dc_alpha_p` stays at the config default for both arms (0.2 when `cct_check` ran; 0.8 now); the true value is in `alpha_p_dc` | read `alpha_p_dc` for `cct_check` / `cct181`; `dc_alpha_p` for the `highp_rerun` campaigns |
| **`;;` at the end of gen rows in `ac.m`** | the base `ac.m` of `faultscreen1_ap1` has `;;` on gen rows; PowerModels fails with `matrix parsing error, inconsistent number of items in each row` on every case | `build_case` in `pm_fault_check.py` replaces `;;` with `;` before writing |
| **Docstring before `using` in a `.jl` script** | Julia errors `LoadError: cannot document the following expression: using PowerModels`; the Python caller raises `RuntimeError` with that text | write the header as `#` comments (as in `ida_pvpq.jl`) |
| **ZIP mix is `alphaI = 1 − alphaP`, so `alphaZ` is always 0** | `split_out_datacenter` callers in `campaign.py` hard-code `alpha_i=1.0 - cfg.dc_alpha_p`; `dc_alpha_p` 0.8 gives 80/20/0, 0.7 gives 70/30/0, and 70/20/10 is not reachable | needs a new `dc_alpha_i` parameter (not yet added); `alphaZ = 1 − alphaI − alphaP` is derived in GridKit `LoadZIP` |
| **16 workers started together race on `seed.m`** | `faultscreen1_ap08` worker 9 crashed at `CampaignConfig` init with `FileNotFoundError: seed.m` while other workers were starting; 43 conditions were missing | rerun the crashed worker alone with its `--start`/`--count` (existing records are skipped); check `records` count against `conditions.csv` after every `highp_rerun.sh` job |
| **`\|` inside a markdown table cell** | breaks the table even in backticks | avoid in docs |
| **Dropping a switched unit from the final `ac.m`** | `build_case_from_m` drops every gen with `GEN_STATUS = 0`, so the machine leaves the dynamic case and its inertia and dynamics change with no error | `_q_restore` returns every unit to service at its fixed `QG`. The restore check: 49 gens, 33 in service, same `PD`/`QD` as the unlimited case |
| **Loop variable read in a skipped loop** | `run_stage` reads the enclosing `it` on its failure path. The AC OPF branch skips the loop that assigns `it`, so 3 of 16 workers per campaign crashed with `NameError` and 188 of 1020 rows were missing. `sacct` still said COMPLETED | assign `it = 0` before the branch. Check the record count against `conditions.csv` and grep `slurm/*.log` for `Traceback`, not `sacct` |
| **Control scored against its own widened limits** | the relaxed AC OPF solve writes ±9999 into `ac.m`, so `q_viol` read zero at every point | `collect()` scores against the original frames (`q_src`). Records from before this fix are invalid for the control, and the control was rerun |
| **`highp_rerun.sh` on a Q-limited campaign** | it drops the baseline rows, but with Q limits the baseline solution changes too (2 → 0 violations), so the arms are not comparable | use `qlim_campaign.sh`, which keeps every row |

---

## 7. Verifying a change

Cheapest first. All are things that have actually caught bugs here.

1. **Import and construct**: `CampaignConfig(name="x", scratch="/tmp/x")` —
   catches syntax and ordering errors in ~1 s.
2. **One condition end to end**: `run_condition` on hour 6164, mid commitment.
   Known-good: `cct_s = 0.2891`, `inertia_MWs = 11864`.
3. **Equivalence, not just success**: when optimising, compare against the slow
   path on `net_mw`, `losses_mw`, `vmin`, `vmax`, `inertia_MWs`. The Julia
   batching was validated this way (0 mismatches).
4. **Tolerances**: the integrator noise floor is ~1e-8 pu. A flat run measured
   `max|omega| = 7.2e-09`. Do not assert equality tighter than that — a
   composition check failed spuriously at a 1e-9 threshold.
5. **Parallel-only bugs**: two of the three races above were invisible to
   single-worker tests. Run ≥2 workers before trusting a sweep.
6. **Markdown**: tables and fences validate with the snippet in
   `work-notes/` history; broken tables render silently wrong.
7. **Q-limit switching** (only when `enforce_q_lims` is touched): run
   `scripts/qlim_pf_check.py` (five points, unlimited vs switched: Q violations,
   units fixed, vmin, losses) and `scripts/qlim_dyn_check.py` (pre-fault drift with
   the default disturbance). Known-good 2026-11-09: switched violations 0 at all five
   points; pre-fault drift 1.3e-6 unlimited, 2.2e-6 switched.
8. **AC OPF path** (only when `ac_opf` or `ac_opf_relax_q` is touched): run
   `scripts/acopf_pilot.py` (five points, AC OPF vs switched). Then check that the
   `acopf_noq` control shows `q_viol` above zero against the real limits. Known-good
   2026-11-09: 0 violations at all points with limits, and 952 of 952 control rows
   violate (max 16 Mvar). Finally confirm the record count equals `conditions.csv`
   and that `slurm/*.log` has no `Traceback`.

---

## 8. Known-good reference values

Hour 6164, `n_commit = 33`, fault bus 66, ACTIVSg200. All verified 2026-11-02:

| configuration | `net_mw` | `cct_s` | `inertia_MWs` |
|---|---|---|---|
| baseline, no datacenter | 1243.2 | **0.2891** | 11,864 |
| 200 MW at bus 108 | 1443.2 | 0.1836 | |
| 150 MW at bus 56 | 1393.2 | 0.1836 | |

The two sited configurations give the *same* CCT: the bisection tolerance is
0.01 s, so they are indistinguishable. Use `net_mw` to tell them apart —
quoting 0.1836 without the configuration is ambiguous.

Case facts: 200 buses, 246 branches, 49 gens (25 coal, 17 gas, 6 wind,
1 nuclear), seed demand 2,178 MW (annual peak), 40 `Genrou`, 109 `LoadZIP`.

---

## 9. Working with ACTIVSg2000

`CaseSpec` is `mu.ACTIVSG2000`. The 2000-bus case is structurally richer than
ACTIVSg200 and four assumptions baked into the 200-bus path had to be lifted.

| assumption | why it fails at 2000 buses |
|---|---|
| one machine class (`Genrou`) | 314 `Genrou` **+ 20 `Gensal`** |
| one exciter, one governor class | `Ieeet1` (324) + `Esdc1a` (10); `Tgov1` (314) + `Hygov` (20); plus `Ieeest` stabiliser (333) |
| every online generator has a dynamic model | **101 online units have none** (renewables, small plant) |
| unit ids are positional within a bus | 13 buses carry 2–11 generators and the numbering shifts |

What was added to `m_to_case_utils`:

* `CaseSpec.machine_classes`, `.machine_suffixes`, `.unit_classes`,
  `.all_suffixes` — multi-class machines and companions.
* `machines(case, spec)` — all machine devices regardless of class.
* `match_gens_by_mva(case, spec, gen)` and `CaseSpec.match_by_mva` — see below.
* `build_case_from_m(..., strict=False)` now reports `counts["gen_static"]`,
  the online units with no dynamic model, instead of only raising.

**The positional-id trap, and how it presents.** At bus 4192 the `.m` lists 10
generators including one that is *offline*; the case carries 9 machines and
skips it, so every later unit at that bus is off by one. The machine then gets
another unit's dispatch, and the simulation dies at initialisation with

```
[ERROR] Tgov1: initial valve position is outside [Pvmin, Pvmax]
[ERROR] System model initialization failed with status 2
```

Only 2 of the 334 machines were affected and both were at multi-generator buses.
Matching on `(bus, MBASE)` fixes it: ratings are stable under the shift, and
units sharing a bus *and* a rating are interchangeable. **If you see a `Tgov1`
valve-position error, suspect generator matching before suspecting the
dispatch.**

**The 101 static units** must be given to `static_to_negative_load` (see
`scripts/test_2k_rebuild.py`) or the dynamic model will be short ~9.8 GW of
injection and fail `IDACalcIC`. This is the same negative-load stand-in used
for wind on ACTIVSg200, `alphaI = 1`.

**`omega` in the monitor output is a deviation from synchronous speed**, not
absolute pu. A healthy machine sits near `0.0`. Testing `abs(omega - 1)` looks
like a 1.0 pu excursion on every healthy run.

**`enable_monitors` used to target `Genrou` only** *(fixed 2026-11-02)*. On
ACTIVSg2000 that left the 20 `Gensal` machines without a `delta` channel, and
since the stability test is COI pole slip on `delta`, those 20 machines were
**invisible to it** — a silently incomplete stability criterion rather than an
error. The signature was in the column arithmetic: 4,649 monitor columns
decompose as `1 + 2x2000 buses + 314 Genrou x2 + 20 Gensal x1`, the lone `x1`
being the missing `delta`. `enable_monitors` now takes `machine_classes`,
default `("Genrou", "Gensal")`, and all 331 machines carry both channels.

> Monitor column counts are a useful integrity check precisely because they are
> exactly predictable: `1 + 2 x n_bus + 2 x n_monitored_machines`. A count that
> does not factor cleanly means some device group was skipped.

**Not yet ported.** `CampaignConfig` is still ACTIVSg200-specific: merit order,
commitment band, `apply_hour` and the wind-column mapping all assume that
case's generator table and profile format. A `faultscreen`-shaped campaign
needs none of it and runs today; an `aleatoric`-shaped campaign over hours and
commitments does.

### Measured costs

| quantity | ACTIVSg200 | ACTIVSg2000 |
|---|---|---|
| `build_case_from_m` | — | 0.12 s |
| DC OPF + AC PF, warm batched | 0.029 s | 0.543 s |
| one dynamic run (`delta`+`omega`, `tmax` 5 s) | ~0.3 s | **5.95 s** |
| monitor output per run | 1.9 MB | 18.2 MB |

**The ratio holds; the absolute cost does not.** A dynamic run is ~10x the warm
dispatch solve at *both* scales (0.3 s vs 0.029 s at 200 buses; 5.95 s vs
0.543 s at 2000), so the stage-split fan-out still pays. What changes is the
absolute number: an 8-probe CCT bisection goes from ~2.5 s to ~48 s, so the
thing worth minimising is the **number of dynamic runs**.

Note the two "stage A" figures mean different things: `0.029 s` is the solver's
warm `solve_time`, `~2.8 s per point` is what a batched campaign actually pays
including `.m` I/O and case building.

### Scripts

| script | does |
|---|---|
| `scripts/scale_probe.py` | dynamic-run cost and monitor size, both cases |
| `scripts/scale_stagea.py` | warm batched DC OPF / AC PF cost, both cases |
| `scripts/test_2k_rebuild.py` | end-to-end rebuild + dynamics on ACTIVSg2000 |
| `scripts/test_2k_screen.py` | `faultscreen`-shaped fan-out on ACTIVSg2000 |
| `scripts/bisect_2k.py` | isolates which rebuild step breaks initialisation |
