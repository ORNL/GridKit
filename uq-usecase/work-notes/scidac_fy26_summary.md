# SciDAC FY26 Summary: GridKit UQ Use-Case

Consolidated summary of the work documented in `~/gridkit/uq-usecase/`. Every number
here is copied from a source document in this repo, and each claim carries a pointer to
where it came from. Nothing is estimated or rounded for effect.

**Documented period:** earliest dated note 2026-06-17 ([`notes_1.md`](notes_1.md)),
latest 2026-09-30 ([`m_to_case_helper.md`](m_to_case_helper.md) and its companion notebook,
both revised that day).
**Summary written:** 2026-09-30.

**Source documents:** 20 markdown files. 10 in `work-notes/` (4,806 lines total, of which
2 are superseded drafts), 6 in `cases/`, plus `contents.md`, `kestrel_install.md`, and
`notebooks/setup.md` / `setup_env.md`. The per-file index is
[`contents.md`](../contents.md); this document is the narrative on top of it.

---

## Contents

- [1. What the project is](#1-what-the-project-is)
- [2. Workstreams](#2-workstreams)
  - [2.1 Kestrel install and build engineering](#21-kestrel-install-and-build-engineering)
  - [2.2 Epistemic UQ campaigns](#22-epistemic-uq-campaigns)
  - [2.3 Hypothesis-driven sensitivity design](#23-hypothesis-driven-sensitivity-design-hawaii-v6v10)
  - [2.4 Power-flow solver investigation](#24-power-flow-solver-investigation)
  - [2.5 MATPOWER to case.json pipeline](#25-matpower-to-casejson-pipeline)
  - [2.6 Grid visualization](#26-grid-visualization)
  - [2.7 Case reference documentation](#27-case-reference-documentation)
  - [2.8 Forward-looking work](#28-forward-looking-work)
- [3. Headline results](#3-headline-results)
- [4. Deliverables](#4-deliverables)
- [5. Status and open items](#5-status-and-open-items)
- [6. Known documentation inconsistencies](#6-known-documentation-inconsistencies)

---

## 1. What the project is

Uncertainty quantification for GridKit phasor-dynamic simulation, on two TAMU synthetic
grids: **Hawaii40** (37 bus) and **Illinois / ACTIVSg200** (200 bus).

[`uq_plan.md`](uq_plan.md) defines two independent tracks:

| Track | Type | Source of uncertainty | Sample basis |
|---|---|---|---|
| Epistemic | Dynamic parameter uncertainty | LHS / Morris over Genrou `H`, `D` | per-campaign, 20 to 16,000 |
| Aleatoric | Load / dispatch variability | Hourly PCM `.m` solutions | 8,760 hours |

The two can be crossed into a 2D scenario-by-parameter study. In practice the **epistemic
track is mature and has produced the bulk of the data**, while the aleatoric track is still
in the tooling phase, because it needs a working `.m` to `case.json` pipeline and a
trustworthy power-flow solver first. Those two prerequisites (Task 0 in the plan) became
substantial workstreams in their own right, covered in sections 2.4 and 2.5.

---

## 2. Workstreams

### 2.1 Kestrel install and build engineering

Source: [`kestrel_install.md`](../kestrel_install.md)

GridKit and its full dependency stack were built from source on Kestrel (RHEL 8.8, Intel
Sapphire Rapids), installed to a shared, reusable prefix at
`/nopt/nrel/apps/cpu_stack/software/gridkit/`.

Six-stage build, scripted as `scripts/1_build_llvm.sh` through `6_build_gridkit.sh`:

| Stage | Component | Version / note |
|---|---|---|
| 1 | LLVM | 16.0.6, `LLVM_ENABLE_RTTI=ON` |
| 2 | Enzyme | v0.0.206 |
| 3 | SuiteSparse / KLU | 7.8.2 / 2.3.4 |
| 4 | SUNDIALS | `develop` branch, required for an adjoint memory fix |
| 5 | IPOPT | 3.14.14 with MUMPS |
| 6 | GridKit | branch `isatkaus/uq-usecase` |

The build gotchas documented are the substantive part of this workstream:

- **`RelWithDebInfo`, not `Release`.** `-O3` over-optimizes the LLVM IR and makes
  `__enzyme_todense` untraceable. Documented as "the root cause of all Enzyme failures."
- **`--gcc-toolchain` on clang only**, never on gcc. Clang otherwise picks up GCC 12 on
  Kestrel, while GridKit needs GCC 13 for the `<format>` header and `GLIBCXX_3.4.29`.
- **No `cray-libsci` module**, since it pins the GNU compiler generation back to 12.1.
- IPOPT moved from HSL MA27/MA57 (broken, `libhsl.so` does not exist on the system) to
  open-source MUMPS via ThirdParty-Mumps.

Also built out here: the MkDocs documentation site (`build_mkdocs.sh`, with a generated
`hooks.py` that converts ` ```math ` fences to `$$` blocks for arithmatex), the
`rebuild_if_updated.sh` daily-maintenance script that merges `origin/develop`, rebuilds,
runs `ctest` and pushes back, and the `wip.sh` staging script.

**Upstream contribution.** The `m_viz` work was contributed to ORNL/GridKit. The original
PR #490 did not follow GridKit directory conventions and lacked CMake wiring and unit
tests; the PI reorganized it into **PR #516**, which `kestrel_install.md` records as the
reference example of a compliant PR. Conventions captured from it: per-feature directories
under `application/notebooks/`, `pyproject.toml` install rather than `sys.path.insert`,
env vars instead of hardcoded paths, pytest suite under `tests/UnitTests/Python/`, and
CMake gating behind `GridKit_ENABLE_NOTEBOOKS`.

### 2.2 Epistemic UQ campaigns

Source: [`all_epistemic_runs.md`](all_epistemic_runs.md), [`notes_4.md`](notes_4.md)

**19 completed sweeps** across the two grids, sampling Genrou inertia `H`.

| Grid | Runs | Sample counts |
|---|---|---|
| Hawaii40 | v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, ian-csv | 20 to 16,000 |
| Illinois (ACTIVSg200) | v1, v2, v6-il, v7-il, v8-il, v9-il, v10-il, ian-csv | 1,000 to 4,000 |

Production runs are **hawaii-v5** (16,000 samples, Gaussian sigma=12%, 4 generators at hop
distances 1 to 3+, 182G) and **illinois-v2** (4,000 samples, Gaussian sigma=12%, 4
generators at hop distances 4 to 9, 144G).

Methodology evolved across versions, visible in the run table: uniform +/-10% (v1, v2),
then a switch to Gaussian sigma=12% as the "utility database" epistemic marginal (v3 to
v5), then uniform +/-20% for the hypothesis-testing series (v6 to v10). Two Morris designs
(`hawaii-ian-csv`, `illinois-ian-csv`) came from an external collaborator and sample all
39 fleet `H` parameters simultaneously rather than 4.

**Totals as of 2026-09-28**, verified by re-summing the run table in
[`all_epistemic_runs.md`](all_epistemic_runs.md):

- **69,045 dynamic simulations** (43,020 Hawaii + 26,025 Illinois)
- **~1.4 TB of Parquet results**, roughly one file per sample

**These numbers are a moving target.** New epistemic sweeps land regularly, so treat any
total quoted here as a snapshot with a date attached, and take
[`all_epistemic_runs.md`](all_epistemic_runs.md) as the live source of truth. Its run
table is the authoritative record; re-sum the N-samples and disk columns for a current
figure rather than quoting this document.

The growth so far, as captured by successive notes:

| Snapshot | Date | Simulations | Data |
|---|---|---:|---:|
| [`notes_2.md`](notes_2.md) exec summary | 2026-08-26 | ~43,000 | not stated |
| [`notes_4.md`](notes_4.md) | 2026-09-11 | 48,020 | ~678 GB |
| [`all_epistemic_runs.md`](all_epistemic_runs.md) | 2026-09-28 | 69,045 | ~1.4 TB |

**Throughput** (from SLURM accounting, `all_epistemic_runs.md`):

- Hawaii40 v5: 16,000 simulations in **~5.2 min compute** on 4 nodes x 4 workers, which is
  3,050 samples/min aggregate and 190/min per worker. Collection added ~6.3 min, so
  **~11.6 min end to end**.
- Illinois v2: each 1,000-sample chunk in ~1 min 40 to 47 sec, 590 samples/min aggregate,
  145/min per worker. One chunk (job 15338688) failed after 1 second and was resubmitted
  successfully; `notes_4.md` states the root cause was not re-diagnosed.

Excluded from totals: `threebusbasic-uq/` and `threebusbasic-v1/`, 77M each, early
smoke tests on the 3-bus toy case with no Parquet results.

### 2.3 Hypothesis-driven sensitivity design (Hawaii v6-v10)

Source: [`notes_3.md`](notes_3.md), [`hawaii_case_tables.md`](../cases/hawaii_case_tables.md)

Rather than sampling more of the same, v6 to v10 were designed to disentangle *why* `H`
sensitivity is or is not detectable. The organizing quantity is **H x mva** (MWs), the
actual stored rotational energy, as distinct from `H` alone (seconds).

Hawaii40 fleet, 39 Genrou: total **H x mva = 6,685.72 MWs**, ranging from 10.3 MWs (`2_1`,
ALOHA69) to 648.9 MWs (`35_8`, KALAELOA138). Illinois, 40 Genrou: total
**15,693.99 MWs**, with `H` from 1.4636 s (`69_1`) to 9.3080 s (`167_1`) and machine
ratings from 2.04 to 682.98 MVA ([`illinois_case_tables.md`](../cases/illinois_case_tables.md)).

Each experiment is N=4,000, LHS, uniform +/-20%, 4 generators:

| Run | Selection criterion | Hypothesis being tested |
|---|---|---|
| v6 | smallest `H` | weakest intrinsic response |
| v7 | farthest from fault | electrical distance masks the signal |
| v8 | smallest `mva` | machine size, confounded with dispatch |
| v9 | largest `H x mva` | **positive control**: if undetectable here, the problem is the metric, not the physics |
| v10 | small `p0`, moderate `H x mva` | lightly dispatched machines contribute little to the swing |

The v9 positive control is the methodological keystone: it makes the whole series
falsifiable. The v10 note is explicit that `p0` and `H x mva` remain correlated in a single
operating-point snapshot, so it is a **partial** isolation test, and a true `p0` isolation
needs multiple operating points, which is the aleatoric track.

The sampling-range rationale is also documented rather than assumed: +/-2 to 5% is
manufacturer nameplate, +/-10 to 15% is the NERC/IEEE database assumption, +/-20% is a
defensible "fairly uncertain database value" prior, and +/-50% is a stress test rather than
a calibration prior.

### 2.4 Power-flow solver investigation

Source: [`pf_helper.md`](pf_helper.md), 924 lines, the single largest investigation

**Central finding: GridKit `solve_pf` has no PV to PQ switching.** Verified across six
independent code paths (`BusPV.cpp`, `GeneratorPV.cpp`, `BusFactory.hpp`, `Kinsol.cpp`,
`solve_pf.cpp`) and then confirmed behaviorally under stress. **PowerModels.jl is therefore
the production PF solver** for all `case.json` generation.

Three-way base-case comparison on ACTIVSg200:

| Comparison | max\|dV\| |
|---|---|
| PM.jl vs TAMU/PowerWorld reference | 0.000011 pu |
| GridKit vs TAMU/PowerWorld reference | 0.029730 pu |
| PM.jl vs GridKit | 0.029728 pu |

PM.jl matches the reference to numerical noise (also max\|dTheta\| = 9.3e-4 deg); the
~0.030 pu deviation is entirely on GridKit's side.

A **safe envelope** was mapped where GridKit's offset stays constant and it remains usable
as an angle-and-flow sanity checker: 0.5x to 1.5x uniform load, N<=2 largest generators
off. Outside it the error grows to 0.089 pu at 2.0x load and 0.279 pu at N=25 generators
off, and the voltage-violation counts go qualitatively wrong: 0 versus 44 to 85 over-voltage
buses at 0.1x to 0.5x load, 0 versus 85 to 177 under-voltage buses at N=10 to 25 generators off.

Two explanatory results worth carrying forward:

- **Why voltage is stiff.** The regulated-bus fraction is 24.5% for ACTIVSg200 and 24.3%
  for Hawaii40. With roughly a quarter of buses holding voltage, the network responds to
  perturbation almost entirely through angle. [`notes_2.md`](notes_2.md) states it as:
  across load +/-80%, wind curtailment, and up to 10 generators offline, "voltage magnitude
  barely moves (< 0.01 pu)."
- **Why it matters for UQ, not just for PF.** Using GridKit PF outside the envelope does
  not crash and produces no error signal. It silently poisons the ensemble initial
  conditions, biasing Sobol/PCE sensitivity indices toward zero for any generator whose Q
  limit was active, which would make `H` variation look *less* detectable than it is. This
  is a direct threat to the section 2.3 conclusions, which is why it was chased down.

Cost of the chosen solver: PM.jl is **~8.5 s per call**, of which ~5.5 s is Julia startup
and ~3 s is the Ipopt solve. For 8,760 hourly scenarios the documented mitigation is a
batched-manifest mode looping inside one Julia process, plus sbatch parallelism.

A pluggable dispatcher `run_pf_solve(solver="pm"|"gridkit")` exists in `pf_utils.py`.
Closing the GridKit gap is scoped as adding post-solve PG/QG writeback to `solve_pf.cpp`,
which is a missing *output* feature rather than a missing solve capability, plus the PV to
PQ switching itself.

### 2.5 MATPOWER to case.json pipeline

Source: [`m_to_case_helper.md`](m_to_case_helper.md), 937 lines, last revised 2026-09-30.
Companion notebook `m_to_case_helper.ipynb`, modules `m_to_case_utils.py`, `aux_to_m.py`,
`perturb.py`, `solver_utils.py`, `pm_utils.py`. This is the most recently active
workstream.

The bridge from a solved MATPOWER `.m` to a GridKit `case.json`, and the prerequisite for
the aleatoric track.

#### Architecture: bootstrap once, produce forever

The organizing idea, and the reason the design is defensible:

```
BOOTSTRAP (once per case, to establish trust)
  .aux + .dyr ──build_case_from_aux_dyr──> rebuild of illinois.case.json
                     └─ diff vs shipped file ──> 0 diffs (with creator overlay)
  .aux        ──write_matpower_m────────> seed .m matching the shipped operating point
  ── PowerWorld exits here, permanently ──

PRODUCTION (every scenario)
  seed .m ──(fix commitment)──> base .m                          [once]
       └─ hourly load / wind ──(PM.jl)──> solved .m              [x8760]
              └─ build_case_from_m ──> case.json ──> DynamicSimulation
```

Both production inputs are open standards (MATPOWER `.m`, which Sienna, PowerModels.jl and
every PCM tool speak, plus PSS/E `.dyr`). Because commitment is fixed during bootstrap,
the 8,760-scenario loop never adds or removes devices; it only patches bus `init.Vr/Vi`,
Genrou `p0/q0` and `LoadZIP.{Pnom,Qnom}`, which is the narrowest and best-tested code path.

#### Why not TAMU's shipped `.m`

`write_matpower_m` reads **no template `.m`**; the table layout is a documented constant in
`aux_to_m.py` and every value comes from the `.aux`. Three measured reasons for rejecting
`case_ACTIVSg200.m` as the seed:

1. It is a **different operating point** from `illinois.case.json`.
2. It is **missing branch (64, 82)**, which exists in the `.aux` and in the shipped case as
   `BR_82_64_1`. The `.m` has 245 branches; the `.aux` and `case.json` have 246.
3. It is a **rounded export**: over the 245 shared branches, `max|dR| = 5.0e-07`,
   `max|dX| = 5.0e-07`, 245/245 within 1e-6. The `.aux` carries ~15 significant figures
   where the `.m` stores ~7.

Two `.aux` parsing traps documented: the `Branch` object appears in **two `DATA` blocks
with different schemas** (180 Lines using `LineR/LineX/LineC`, 66 Transformers using
`LineR:1/LineX:1/LineC:1`, so reading `LineR` for a transformer silently yields `NaN`), and
all 66 transformers are `TAP = 1.0`, `SHIFT = 0.0`, meaning ACTIVSg200 has no off-nominal
taps or phase shifters. `gencost` is fully derivable from PowerWorld's cubic IO curve
scaled by fuel price (`c0 = GenFixedCost`, 49/49 exact). Fields genuinely absent from a
`.aux` are reported through `AUX_ABSENT` rather than silently defaulted; note `VMAX/VMIN`
fall back to 1.1/0.9, which matters if OPF with voltage constraints is ever added.

#### Validated against three independent oracles

Nothing rests on "the mapping looks right":

| Oracle | Check | Result (pure `.aux` build) |
|---|---|---|
| 1 | AC PF mismatch `S = V conj(Y_bus V)` vs scheduled injections | max\|dP\| 0.0063 MW, max\|dQ\| 0.0018 MVAr |
| 1 | network losses vs `PG - PD` | 24.749 / 24.749 MW |
| 2 | PM.jl re-solve of the generated `.m` | 4 iterations, constraint violation 1.6e-13, `LOCALLY_SOLVED` |
| 3 | `.aux` + `.dyr` to `case.json`, with creator overlay | **0 diffs** vs shipped |
| 3 | same, without overlay | only `197_1` differs above 1e-4 (14 params) |
| 3 | seed `.m` to `case.json` | **0 diffs @ 1e-6**, 1 @ 1e-7, 13 @ 1e-8 |

Oracle 1's residual is write precision, not modelling error, and it is also what confirms
branch (64, 82) belongs: omit it and buses 64 and 82 show large mismatch. Converging in 4
iterations at 1e-13 is what an already-solved input looks like. The sub-1e-4 rows in
oracle 3 are `.dyr` print rounding: the `.dyr` prints 4 decimals while the JSON keeps full
precision (`49_1` `H` is 5.2115 in the `.dyr`, 5.21145296 in the JSON).

**Provenance conclusion:** `illinois.case.json` comes from the PowerWorld `.aux` plus the
PSS/E `.dyr`, not from `case_ACTIVSg200.m`. Field-level evidence:

| Check | Result |
|---|---|
| `Genrou.p0/q0/mva` vs `.aux` Gen | 40/40 exact |
| Bus `init.Vr/Vi` vs `.aux` Bus | 300/400 bit-exact, 100/400 within 1e-7, 0 mismatches |
| `LoadZIP.Pnom/Qnom` vs `.aux` Load | 320/320 agree, worst residual 5.0e-09 |
| Shunt `Qnom = -(SSNMVR/100)*\|V\|^2` | 4/4, residuals 6.8e-9 to 2.7e-8 |
| `Genrou` dynamics vs `.dyr` | 39/40 exact at the `.dyr` 4-dp print precision |

The single exception, `197_1`, is handled by a **creator overlay generated by diffing
rather than hand-written**, so it stays correct if GridKit reships the case. The notebook
health check re-derives it on every run and asserts the rebuild still gives 0 diffs.

#### Operating point established by the build

| Quantity | Value |
|---|---|
| total PD / QD | 2178.00 MW / 620.73 MVAr |
| total PG (online) | 2202.75 MW |
| losses | 24.749 MW |
| gens online / total | 40 / 49 |
| bus types | 151 PQ, 48 PV, 1 slack (bus 189) |
| branches (transformers) | 246 (66) |
| LoadZIP devices | 164 = 160 real loads + 4 shunts |
| load buses / multi-ZIP | 108 / 36, distribution `{1:72, 2:22, 3:12, 4:2}` |

Zero `PG > PMAX` violations against the `.aux` limits. Five `PG < PMIN` rows remain, all
wind (buses 65, 104, 105, 114, 147) curtailed below economic minimum, which is fine for PF
and relevant only if OPF is added.

#### Test suite

Beyond the identity round trip, eight perturbation-driven tests (4a to 4h) exercise each
code path in isolation: bus init only, Genrou dispatch only, offline-gen removal (including
the `strict=True` raise on the rebuild-required direction), single-bus LoadZIP patch,
multi-ZIP proportional scaling (sum matches `PD/baseMVA` exactly, inter-device ratio
preserved to 1e-12), shunt exclusion across all four shunts, zero-load bus, and a full
end-to-end 0.8x load perturbation through PM.jl.

The tolerance philosophy is stated explicitly and is worth keeping: the identity round trip
does **not** compare against mathematical truth, it compares against a *different-solver
equilibrium*, since the base JSON came from PowerWorld and the round trip uses PM.jl. So
LoadZIP values are *expected* to change; the invariant is that the multi-ZIP proportion
within a bus is preserved.

#### End-to-end run, verified 2026-09-30

The notebook now executes the complete loop (perturb, write `.m`, PM.jl solve, build
`case.json`, write `solver.json`, run `DynamicSimulation`, reduce to QoI) across five smoke
scenarios. All cells executed successfully; results from the notebook outputs:

| Scenario | PF status | PF iters | total PD (MW) | gens online | Genrou in case | sim rc | wall (s) |
|---|---|---:|---:|---:|---:|---:|---:|
| base | LOCALLY_SOLVED | 4 | 2178.00 | 40 | 40 | 0 | 35.7 |
| load_90 | LOCALLY_SOLVED | 4 | 1960.20 | 40 | 40 | 0 | 9.3 |
| load_110 | LOCALLY_SOLVED | 4 | 2395.80 | 40 | 40 | 0 | 9.0 |
| wind_50 | LOCALLY_SOLVED | 4 | 2178.00 | 40 | 40 | 0 | 9.1 |
| gen_off | LOCALLY_SOLVED | 4 | 2178.00 | 39 | 39 | 0 | 8.9 |

Every scenario converged, every dynamic run exited 0 with zero `[ERROR]` lines, and the
device count matched the requested commitment in all five. Per-scenario cost after the
first is **~9 s**, dominated by Julia startup.

Quantities of interest, fault at bus 2 for 1.0 to 1.1 s, `freq_base = 60.0` Hz read from
the case rather than assumed:

| Scenario | machines | nadir f (Hz) | peak f (Hz) | max angle spread (deg) | Vmin during fault |
|---|---:|---:|---:|---:|---:|
| base | 40 | 59.645 | 60.478 | 87.24 | 0.0212 |
| load_90 | 40 | 59.622 | 60.471 | 88.94 | 0.0212 |
| load_110 | 40 | 59.669 | 60.485 | 85.46 | 0.0212 |
| wind_50 | 40 | 59.645 | 60.478 | 87.72 | 0.0212 |
| gen_off | 39 | 59.661 | 60.481 | 77.66 | 0.0212 |

No NaNs, all runs reached `t_end = 10.0` s. These are **smoke tests proving the pipeline
reacts and stays numerically sound**, not a physics study: a uniform "all load x 0.90" is
not a meaningful operating condition, and the notebook says so explicitly. Production
driving will be per-bus from hourly profiles.

#### Two correctness traps documented in the notebook

- **`element_id` is not a bus number.** It is the 0-based ordinal in which GridKit
  registers `BusFault` devices while reading `case.json`, and for Illinois that ordering is
  not the bus ordering: bus 1 is `element_id` 1 while bus 2 is `element_id` 0. Resolved by
  `su.fault_element_id(case, bus)`, never by hand.
- **`omega` is a deviation, not an absolute speed.** Verified in `GenrouImpl.hpp`: the
  state initialises to `omega(0.0)` (line 465) and enters the algebraic equations as
  `(ONE + omega)` (line 481). So synchronous speed is `omega == 0` and frequency is
  `freq_base * (1 + omega)`. Getting this backwards would corrupt every frequency QoI.

#### Bugs and limitations found

- **The `Vnom` non-parameter.** The superseded `m_to_case.py` writes `Vnom` at five call
  sites, but GridKit's LoadZIP has no such parameter. Verified consequence: the simulator
  ignores it, `mon.csv` is byte-identical, and it costs one spurious `[ERROR]` line per
  device.
- **`m_to_case.py` was never runnable** in `h-py312-basic`: it does
  `from matpowercaseframes import Case`, but the installed v2.1.0 exports `CaseFrames`
  only, and the `try/except` rewrites the failure into a misleading "pip install" message
  for a package that is already installed. `m_to_case_utils.py` replaces it. Only
  `m_to_case_helper.ipynb` imported the old module, so **the completed epistemic runs are
  unaffected**.
- **Multi-unit buses are guarded, not supported.** `parse_dyr` keys on `(bus, model)` and
  drops the PSS/E unit id, and `gen_ids` uses positional rank. Verified harmless on
  ACTIVSg200 (49 generators on 49 distinct buses, every unit id `'1'`, zero collisions) but
  **will bite on ACTIVSg2000 and ACTIVSg10k**. Guards: `parse_dyr` raises on collision,
  `gen_ids` warns.
- **Adding a generator is untested.** `_add_gens_from_dyr()` exists but every round trip so
  far reports `gen_added: 0`, so that path has never executed against real data.
- **`LoadZ` (Hawaii) not implemented.** Hawaii uses `LoadZ` with params `R`, `X`, requiring
  a different patch algorithm than Illinois's `LoadZIP`.
- **Slack has no headroom.** Bus 189 runs at PG 569.141817 against PMAX 569.150019, i.e.
  **0.008202 MW**. Switching off anything larger than a few MW needs an explicit
  redistribution rule, which `perturb.set_gen_status(..., redistribute="headroom")` now
  provides. The documented probe unit is `65_1` (PAXTON 1 1): H x mva = 1135 MWs, 7.2% of
  the 15,694 MWs fleet total, but only 4 MW of dispatch, so it changes inertia
  substantially while barely moving the power balance.
- **GridKit cannot enable or disable a generator at runtime**, so `GEN_STATUS=0` machines
  must be absent from `case.json`. This is why commitment is frozen before the hourly loop.
- **ZIP mix is preserved verbatim**; `alphaI`/`alphaP` are never modified.

### 2.6 Grid visualization

Source: [`notes_1.md`](notes_1.md), [`notes_2.md`](notes_2.md)

`m_viz.ipynb` plus `m_viz_utils.py`: interactive geographic Plotly maps of any MATPOWER
`.m` case showing buses, branches, generators and branch loading, with GridKit `case.json`
device IDs cross-referenced in hover tooltips. Tested on **4 TAMU grids from 37 to ~10,000
buses** (Hawaii40, ACTIVSg200, ACTIVSg2000, ACTIVSg10k). Geo source falls back from `.gic`
to `.AUX` depending on case. This is the work that became upstream PR #516.

### 2.7 Case reference documentation

Source: `cases/` (6 files)

Per-grid reference docs, several machine-generated from source data by
`gen_hawaii_tables.py` so they cannot drift:

| File | Contents |
|---|---|
| [`hawaii.md`](../cases/hawaii.md) | 37 buses, 45 gens, 89 branches; `.m` to `hawaii.json` field mapping; UQ parameter selection v1 to v5 |
| [`hawaii_m_tables.md`](../cases/hawaii_m_tables.md) | Parsed MATPOWER tables. Total PD 1136.29 MW, PG online 1154.66 MW, imbalance 18.37 MW |
| [`hawaii_case_tables.md`](../cases/hawaii_case_tables.md) | Genrou params ranked by H x mva |
| [`illinois.md`](../cases/illinois.md) | 49 gens (25 coal, 17 ng, 6 wind, 1 nuclear); patching spec; UQ selection v1/v2 |
| [`illinois_m_tables.md`](../cases/illinois_m_tables.md) | 200 buses (151 PQ, 48 PV, 1 slack); PD 1475.69 MW, QD 420.55 MVAr; V in [1.010, 1.056] pu; wind PG 536.20 MW |
| [`illinois_case_tables.md`](../cases/illinois_case_tables.md) | 40 Genrou ranked by H x mva |

Two modelling facts recorded here matter downstream: `hawaii.json` carries 39 of 45
generators (6 offline large-MVA units excluded, each modeled as a Genrou + Tgov1 + Ieeet1
triplet), and `illinois.json` carries 40 of 49, of which **two (buses 161 and 197) are kept
despite `GEN_STATUS=0`**, flagged as an open question for the GridKit dev team. Also
recorded: JSON `Pnom` runs ~1.49x higher than `PD/100` from the TAMU base case, ratio range
1.36 to 1.60.

### 2.8 Forward-looking work

**Aleatoric track** ([`uq_plan.md`](uq_plan.md)): 8,760 hourly PCM operating points. Task 0
(scenario generation) is essentially complete now that sections 2.4 and 2.5 landed. Task 0b,
a fault-response sensitivity prototype to check the working hypothesis that flat pre-fault
`|V|` does not imply flat fault response, is recorded as **planning note only, not yet
implemented**.

**Graph transformers** ([`gt_ideas.md`](gt_ideas.md), 399 lines): an explicitly labelled
idea survey rather than a committed plan. Five candidate problems (inverse parameter
inference, fast forward surrogate, fault localization, identifiability-driven sensor
placement, operating-point-conditioned dynamics), backbone survey (Graphormer, GraphGPS,
Exphormer, TokenGT), and three sketched data runs (R1 fault-location sweep, R2 expanded
parameter set, R3 aleatoric coupling). Carries its own caveat that attention is quadratic
in node count and that any result must be benchmarked against MPNN/GAT baselines.

---

## 3. Headline results

1. **69,045 dynamic simulations, ~1.4 TB**, across 19 sweeps on 2 grids, as of 2026-09-28
   and still growing. See [`all_epistemic_runs.md`](all_epistemic_runs.md) for the current
   count.
2. **16,000-sample sweep in ~11.6 minutes end to end** on 4 nodes (5.2 min compute,
   6.3 min collect), which is 190 samples/min per worker.
3. **GridKit `solve_pf` lacks PV to PQ switching**, proven at source level across six code
   paths and confirmed behaviorally; PM.jl reproduces the TAMU/PowerWorld reference to
   1.1e-5 pu where GridKit is off by 0.030 pu. PM.jl adopted as production.
4. **A safe envelope for GridKit PF was quantified** (0.5x to 1.5x load, N<=2 gens off)
   together with the reason it matters: outside it, ensemble initial conditions are
   silently poisoned and sensitivity indices are biased toward zero.
5. **Provenance of `illinois.case.json` established** as `.aux` + `.dyr`, validated against
   three independent oracles (AC power-flow mismatch, PM.jl re-solve, round trip to the
   shipped file at 0 diffs), with a diff-generated overlay for the single divergent machine
   `197_1`. The whole `.m` to `case.json` chain is reproducible and self-checking.
6. **The scenario pipeline runs end to end** (perturb, PM.jl, `case.json`,
   `DynamicSimulation`, QoI): 5 scenarios, all `LOCALLY_SOLVED` in 4 iterations, all
   dynamic runs exit 0 with zero error lines, ~9 s per scenario after the first.
7. **Full GridKit dependency stack built on Kestrel**, with the Enzyme/`-O3` root cause
   and the GCC 13 toolchain requirements documented.
8. **Upstream contribution merged path**: `m_viz` via ORNL/GridKit PR #516.

---

## 4. Deliverables

**Notebooks** (`notebooks/`): `uq_setup.ipynb` (main production sweep: sampling, SLURM
submission, collection, QC), `m_viz.ipynb`, `gridkit_viz.ipynb`, `pf_helper.ipynb`,
`pm_helper.ipynb`, `m_to_case_helper.ipynb`, `gridkit_helper.ipynb`.

**Python helpers** (`py-utils/`): `gridkit_utils.py` (sampling, run dirs, collection),
`slurm_utils.py` (submission, status, failure recovery), `m_to_case_utils.py`,
`aux_to_m.py`, `perturb.py`, `solver_utils.py`, `pf_utils.py`, `pm_utils.py`,
`m_viz_utils.py`, `gen_hawaii_tables.py`.

**Scripts** (`scripts/`): the six-stage dependency build, `rebuild_if_updated.sh`,
`build_mkdocs.sh`, `wip.sh`, `create_h-py312-basic.sh`.

**Data** (`/kfs2/projects/scidac/scidac-data/`): `gridkit-runs/` (all epistemic output),
`pcm-runs/` (hourly `.m` inputs for the aleatoric track).

**Documentation**: this repo's `work-notes/` and `cases/`, indexed by
[`contents.md`](../contents.md) and published through MkDocs.

---

## 5. Status and open items

| Item | Status |
|---|---|
| Epistemic pipeline (sample, run, collect, QC) | Production, 19 sweeps completed |
| v6-v10 sensitivity **analysis** | Runs complete; ranking analysis is the stated next step in [`notes_2.md`](notes_2.md) |
| PF solver choice | Decided: PM.jl. GridKit parity plan scoped but not implemented |
| `.m` to `case.json` pipeline | Provenance and mapping verified against 3 oracles; full scenario loop runs end to end (2026-09-30). 6 notebook-dependent checks still marked `[TODO]` in [`m_to_case_helper.md`](m_to_case_helper.md) |
| Aleatoric 8,760-hour track | Prerequisites done; `run_scenario` is already array-job shaped, but the production notebook (`aleatoric_setup`) and the hourly-profile loader are not yet written |
| Task 0b fault-response prototype | Planning note only |
| Graph transformer work | Idea survey only |
| Fault-location sweep | Identified as next data run in [`notes_2.md`](notes_2.md) and `gt_ideas.md` R1 |

Known limitations carried in the docs, worth restating because they bite at scale:

1. **Multi-unit buses.** `parse_dyr` raises if a bus carries more than one unit per model.
   ACTIVSg200 has 49 generators on 49 distinct buses so it never fires, but ACTIVSg2000 and
   larger will need `(bus, model, unit)` keying.
2. **Constant-impedance loads.** All LoadZIP devices have `alphaI = alphaP = 0`, hence
   `alphaZ = 1`. Voltage-dependent load response is 100% constant impedance everywhere.
3. **Slack headroom.** 0.008202 MW on bus 189, as above.
4. **The two anomalous Illinois generators** (buses 161, 197) kept despite `GEN_STATUS=0`.
5. **Generator addition is untested** (`gen_added: 0` on every round trip so far), and
   Hawaii's `LoadZ` patch path is not implemented.

---

## 6. Known documentation inconsistencies

Flagged rather than silently corrected, since resolving them needs a judgement call:

- **Simulation and data totals are living numbers.** Every count in this document is a
  dated snapshot, and new epistemic sweeps land often. Where two documents disagree on a
  total, the one with the later date is usually just newer rather than more correct.
  [`all_epistemic_runs.md`](all_epistemic_runs.md) is the live table; its Headline metrics
  section is maintained, while its Appendix ("Bragging points") lags behind and should not
  be quoted.
- [`notes_4.md`](notes_4.md) (48,020 / ~678 GB) and [`notes_2.md`](notes_2.md) (~43,000) are
  correct **for their dates** and should be read as snapshots, not as contradictions.
- `cases/illinois.md` is stale in several places, superseded by the verified counts in
  [`m_to_case_helper.md`](m_to_case_helper.md): load-bus counts are 108 buses / 36 multi-ZIP
  (not 109 / 38, and multi-ZIP does not mean exactly 2), there are 4 shunts and not 1, and
  it lists `Vnom` as a `LoadZIP` parameter when GridKit has no such parameter.
- `cases/hawaii.md` names Hawaii's load class `Load`; it is actually `LoadZ`, with params
  `R` and `X` rather than `Pnom`/`Qnom`.
- `pf_helper_old.md` and `pf_helper.md.bak` are superseded drafts. Per
  [`contents.md`](../contents.md), do not cite them.
