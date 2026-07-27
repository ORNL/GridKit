# UQ Plan: GridKit Aleatoric + Epistemic Uncertainty Quantification

## Overview

Two independent UQ tracks applied to the **Illinois (ACTIVSg200) case**:

| Track | Type | Source of uncertainty | Samples |
|-------|------|----------------------|---------|
| Aleatoric | Load/dispatch variability | Hourly PCM (.m solutions) | 8760 |
| Epistemic | Dynamic parameter uncertainty | LHS over Genrou H, D, ... | ~100-500 |

Tracks can be combined for a 2D study (scenario x parameter sample) or run independently.

---

## Task 0: Scenario generation pipeline (PF solver, prerequisite for Track 1)

### Goal

The ACTIVSg200 system is likely overbuilt: too much inertia, so fault disturbances are well
damped and not very interesting for dynamic UQ. The plan is to manually take the base `.m`
file, turn a few spinning generators off (reducing system inertia and making it dynamically
more fragile), then vary demand/wind values to produce 8760 realistic operating-point scenarios.
Each scenario must have a valid AC power-flow solution before it can become a `case.json`.

This is **not** full generator commitment optimization (that would require DCOPF via a PCM tool
like Sienna; deferred to later). The goals are:

1. **Vary demand and wind generation** in the base `.m` file to create a range of operating conditions (as opposed to using the single static base-case dispatch).
2. **Thin the generator fleet** by setting selected spinning generators to `GEN_STATUS=0`, reducing system inertia to make fault responses more interesting for dynamic UQ.

After making these changes, a valid AC power-flow solution is needed — not an optimal dispatch, just a feasible solution that satisfies network equations and can serve as initial conditions for `case.json`.

### Workflow

```
base case_ACTIVSg200.m
  |
  +-- manually set selected gens to GEN_STATUS=0 (reduce inertia)
  |
  +-- substitute hourly demand/wind values (8760 scenarios)
  |
  +-- solve PF for each scenario  <-- Task 0 focus
  |
  +-- export solved .m files (same format as pcm-runs/ dataset)
  |
  +-- m_to_case.py: convert each solved .m → case.json
        (offline gens = GEN_STATUS=0 must be ABSENT from JSON,
         since GridKit has no mechanism to disable a Genrou device at runtime)
```

### Off-generator rule (important)

GridKit has no way to turn a generator on or off during simulation. Therefore:
- If `GEN_STATUS=0` in the `.m` file, that generator must be **absent** from `case.json`
- This means `m_to_case.py` must rebuild the Genrou/Tgov1/SexsPti device triplet list
  from scratch for each scenario, not just patch `p0`/`q0` on a fixed base JSON
- The current `illinois.json` (with its 2 anomalous offline Genrous at buses 161, 197)
  is the base for the epistemic track only; the aleatoric track needs scenario-specific JSONs

### PF solver options

See **[`pf_helper.md`](pf_helper.md)** for detailed notes, design decisions, and implementation plan.

Brief summary of options considered:

#### Option A: GridKit PowerFlow
Check if GridKit's `SteadyState/Kinsol` infrastructure supports standalone AC power flow
on a MATPOWER-format case. Pros: no external dependency. Cons: unclear scope, needs investigation.

#### Option B: PowerModels.jl
Use PowerModels.jl to solve ACOPF/PF and export the solution as a `.m` file, exactly as was
done to produce the existing `pcm-runs/ACTIVSg200_wind_demand/` dataset.
Pros: already proven for this case. Cons: Julia dependency.

### Decision

TBD. Option B (PowerModels.jl) is the safe choice since it is already proven. Option A worth
investigating first to avoid the Julia dependency. Full generator commitment optimization
(e.g. DCOPF via Sienna PCM) is out of scope for now.

---

## Track 1: Aleatoric UQ (PCM operating conditions)

### Goal
Run GridKit dynamic simulation for each of the 8760 hourly power-flow solutions
from the PCM (ACTIVSg200), using the Illinois case (`illinois.json`, 200-bus),
capturing how the system dynamic response varies with real operating conditions.

### Data
- `.m` solution files: `/kfs2/projects/scidac/scidac-data/pcm-runs/ACTIVSg200_wind_demand/matpower/ACTIVSg200_sol_hour_*.m`
- 8760 files (one per hour), sorted by hour index
- Each file is a full MATPOWER case with a converged AC PF solution: `mpc.bus` (Vm, Va per bus), `mpc.gen` (PG, QG per gen), branches, etc.
- The PCM dataset is `ACTIVSg200_wind_demand`: both load and wind generation vary across hours.

### What the PCM `.m` solution files contain

Each `.m` file stores the result of a PCM dispatch plus an AC PF solve for one operating
hour. The fields relevant to GridKit initialization are:

| MATPOWER field | Column | Units | Meaning |
|---|---|---|---|
| `mpc.bus` | `VM` (col 8) | pu | Voltage magnitude at each bus |
| `mpc.bus` | `VA` (col 9) | degrees | Voltage angle at each bus |
| `mpc.bus` | `PD` (col 3) | MW | Active load at each bus |
| `mpc.bus` | `QD` (col 4) | MVAr | Reactive load at each bus |
| `mpc.gen` | `PG` (col 2) | MW | Active power output per generator |
| `mpc.gen` | `QG` (col 3) | MVAr | Reactive power output per generator |
| `mpc.gen` | `GEN_STATUS` (col 8) | 0/1 | Generator online/offline |

### How a MATPOWER `.m` solution maps to a GridKit `case.json`

GridKit's `case.json` has two per-scenario fields that must be updated from the `.m` file:

**1. Bus initial conditions (`init.Vr`, `init.Vi`)**

Each bus entry in `case.json` has:
```json
{ "number": 42, "class": "bus", "init": {"Vr": 1.02, "Vi": -0.12}, ... }
```
`Vr` and `Vi` are the rectangular components of the bus voltage phasor. They come from
the polar form in the `.m` file:
```
Vr = Vm * cos(Va_rad)
Vi = Vm * sin(Va_rad)     where Va_rad = Va_deg * pi / 180
```
The `number` field in `case.json` equals `BUS_I` in `mpc.bus` (confirmed by `attach_json_ids`).

**2. Generator initial dispatch (`params.p0`, `params.q0`)**

Each Genrou device entry in `case.json` has:
```json
{ "class": "Genrou", "id": "126_1", "params": {"p0": 0.39, "q0": 0.04, ...} }
```
`p0` and `q0` are in per-unit on the system MVA base (100 MVA for ACTIVSg200):
```
p0 = PG / baseMVA
q0 = QG / baseMVA
```
The Genrou `id` is `"{bus}_{rank}"` where `rank` is the 1-based position of that generator
among all generators at that bus in the `mpc.gen` matrix (same convention as `attach_json_ids`).

**3. Load (NOT yet patched)**

The `case.json` contains `LoadZIP` devices with constant-power parameters (`P0`, `Q0`). The
PCM `.m` file has per-hour `PD`/`QD` values per bus in `mpc.bus`. Currently `m_to_case.py` does
**not** patch load parameters. For Track 1, whether load variation is needed depends on
whether the dynamic simulation is sensitive to the load model (likely yes for longer sims).
This is an open question — see [Open Questions](#open-questions), item 4.

### What is and is not changing across PCM scenarios

The 8760 scenarios differ in:
- **Wind generation**: wind PG varies hour by hour (total wind 0–536 MW range in ACTIVSg200)
- **Load**: total PD and its distribution across buses varies with time of day and season
- **Generator commitment**: which ng generators are online varies (GEN_STATUS per hour)
- **Bus voltages**: the AC PF solution (Vm, Va) is unique to each dispatch point

The 8760 scenarios share (held constant in the base `illinois.json`):
- Network topology (branches, their R/X/B parameters)
- Generator dynamic parameters (H, D, Xd, Td', ...)
- Fault specification in `solver.json`

### Implementation

**`m_to_case.py`** (DONE) in `py-utils/`:
- `patch_case_from_m(base_case_json, m_path, output_path)`:
  - Patches bus `init.Vr` / `init.Vi` from `mpc.bus` VM/VA columns
  - Patches Genrou `params.p0` / `params.q0` from `mpc.gen` PG/QG columns
  - Does NOT patch load, dynamic params, or solver settings
- `patch_cases_from_m_list(base_case_json, m_paths, output_dir)`: batch version
- `scenario_summary(paths)`: returns DataFrame of (scenario, gen_id, p0, q0) for sanity checks

**`uq_setup.ipynb`** (TODO):
- Select subset or all 8760 scenarios
- Generate scenario case JSONs via `patch_cases_from_m_list`
- For each scenario: `run_sample(scenario_dir, runner, solver_fn)`
- `collect_and_save` with extra `scenario_id` column

### Key ID mapping
- `MATPOWER BUS_I` == GridKit bus `number` field (confirmed by `attach_json_ids`)
- GridKit Genrou `id` = `"{bus}_{rank}"` (rank = 1-based position among generators at that bus in `mpc.gen`)
- baseMVA = 100 MVA for ACTIVSg200 (stored in each `.m` file header)

---

## Track 2: Epistemic UQ (dynamic parameter uncertainty)

### Goal
Quantify sensitivity of dynamic response to uncertain Genrou parameters using
Latin Hypercube Sampling. 

### Sampling method: LHS vs plain random

For a given marginal distribution (uniform or Gaussian), **Latin Hypercube Sampling (LHS)**
stratifies the CDF into N equal-probability intervals and draws exactly one sample per interval
per dimension, then randomly pairs the draws across dimensions. This guarantees uniform coverage
of each marginal with no gaps or clusters, unlike plain random sampling which can leave large
empty regions by chance. For N=1000 in 4 dimensions, LHS produces a sample set that is much
more space-filling than 1000 i.i.d. draws, which matters when estimating statistics or training
a surrogate from a limited simulation budget.

### Completed runs

All runs sample **H (inertia constant)** only, one Genrou per bus. See case docs for full parameter tables and rationale.

| Run | Case | N | Dist | Generators | Solver overrides | Data path |
|-----|------|---|------|------------|-----------------|-----------|
| threebusbasic-v1 | ThreeBus/Basic (3 bus) | 100 | uniform ±10% | `genrou_2_1`, `genrou_3_1` | none | `.../threebusbasic-v1/` |
| hawaii-v1 | Hawaii (37 bus) | 20 | uniform ±10% | `2_1`, `23_1`, `34_1`, `35_1` | none | `.../hawaii-v1/` |
| hawaii-v2 | Hawaii (37 bus) | 1K | uniform ±10% | `2_1`, `23_1`, `34_1`, `35_1` | none | `.../hawaii-v2/` |
| hawaii-v3 | Hawaii (37 bus) | 1K | normal std=12% | `2_1`, `23_1`, `34_1`, `35_1` | none | `.../hawaii-v3/` |
| hawaii-v4 | Hawaii (37 bus) | 4K | normal std=12% | `2_1`, `23_1` | none | `.../hawaii-v4/` |
| **hawaii-v5** | Hawaii (37 bus) | **16K** | normal std=12% | `2_1`, `23_1`, `34_1`, `35_1` | none | `.../hawaii-v5/` |
| illinois-v1 | Illinois (200 bus) | 1K | normal std=12% | `126_1`, `135_1`, `115_1`, `189_1` | tmax=10s, fault t=1.0/1.1s | `.../illinois-v1/` |
| **illinois-v2** | Illinois (200 bus) | **4K** | normal std=12% | `126_1`, `135_1`, `115_1`, `189_1` | tmax=10s, fault t=1.0/1.1s | `.../illinois-v2/` |

All paths are under `/kfs2/projects/scidac/scidac-data/gridkit-runs/`. Bold = current best version per case.

Case documentation:
- [`cases/hawaii.md`](../cases/hawaii.md) — bus/gen tables, parameter selection rationale, version history (v1–v5)
- [`cases/illinois.md`](../cases/illinois.md) — bus/gen tables, hop distances from fault, version history (v1–v2)

Utility functions in `gridkit_utils.py`:
- `generate_samples(param_specs, N, seed, method)` — LHS or random sampling over param specs
- `make_run_dir(...)` — creates `run_NNN/`, copies and patches `case.json` + `solver.json`; overwrites `mon[]` arrays in each bus/device entry with `MONITORS_BY_CLASS` values
- `run_sample(...)` — subprocess call to `DynamicSimulation`
- `collect_parallel(run_root, out_path, n_workers)` — `ThreadPoolExecutor` collect of `mon.csv` → `run_NNN.parquet`; I/O-bound, GIL released during reads/writes; safe range on Kestrel Lustre: 32 (conservative) to 64 (aggressive)
- `collect_and_save(...)` — stacked mode (single `results.parquet`); use only for N < 500

#### [`gridkit_helper.ipynb`](../notebooks/gridkit_helper.ipynb) — single test run

Use to validate a case config before launching a full UQ sweep.

- **runner setup**: locates `DynamicSimulation` binary, adds to `PATH`
- **simple run config**: `BASE_CASE_DIR`, `CASE_JSON`, `SOLVER_JSON`, `RUN_DIR`, `MONITORS_BY_CLASS`, `SOLVER_OVERRIDES`
- **simple run execute**: copies case files, patches solver if `SOLVER_OVERRIDES` set, runs binary, loads and shows `mon.csv`
- **plot cell**: time-series of all monitored signals

#### [`uq_setup.ipynb`](../notebooks/uq_setup.ipynb) — full epistemic UQ sweep

Workflow: imports → runner setup → ONE case config cell → shared experiment cells.

**Case config cells** (run exactly one):
- `### hawaii setup` — Hawaii medium case; current: v5, 16K samples, 4 generators (`2_1`, `23_1`, `34_1`, `35_1`), Gaussian std=12%, no solver overrides, monitors: `Vm`/`Va` (buses), `delta`/`omega` (genrou)
- `### illinois setup` — Illinois large case; current: v2, 4K samples, 4 generators at increasing hop distance from fault (`126_1` hop 4, `135_1` hop 5, `115_1` hop 7, `189_1` hop 9), Gaussian std=12%, solver override: tmax=10 s, fault at t=1.0/1.1 s

**Shared experiment cells** (run in order after config):
1. Write `meta.yml` to `UQ_RUN_ROOT` (records case, sampling params, solver overrides, monitor config)
2. Generate `samples.csv` via LHS (or random; `SAMPLE_METHOD`)
3. Create `run_NNN/` dirs — patched `case.json` (params + mon arrays) and `solver.json`; 3-digit zero-padded format
4. **Run simulations** — choose ONE:
   - **Option A: serial** — loop on interactive node; suitable for N ≤ ~200 or quick tests
   - **Option B: SLURM parallel** — one sbatch script per node, each runs a sample-index slice via `xargs -P WORKERS_PER_NODE`
5. Collect results (`per_run` mode default: one `run_NNN.parquet` per run) via `collect_parallel()`; or SLURM collect job if `COLLECT_IN_SLURM=True`
6. Results size + data quality check (file count, NaN scan, column consistency, soft physical bounds)
7. Sanity check plots (random subset of runs, one signal per variable)

**SLURM config** (`SLURM_ACCOUNT`, `SLURM_PARTITION`, `WALLTIME`, `N_NODES`, `WORKERS_PER_NODE`, `SLURM_QOS`, `SLURM_MODULES`):
- `COLLECT_IN_SLURM=True` — writes and auto-submits `slurm/collect.sh` with `--dependency=afterok` on all sim jobs; collect job calls `collect_parallel()` via inline Python heredoc
- `COLLECT_IN_SLURM=False` — collect manually via notebook collect cell after sim jobs complete
- Status cell: polls `sacct` for all sim + collect job IDs; scans `chunk_NN_failed.txt`; only prints "safe to proceed" when no RUNNING/PENDING jobs remain

**Data quality check** (results size cell):
- File count vs `N_SAMPLES`
- Column consistency across sampled files
- NaN scan over `QC_SAMPLE_N=100` randomly chosen runs (fixed or random seed toggle)
- Soft physical bounds: `Vm` [0, 1.6], `Va` [-4, 4], `omega` [-1, 1] (deviation, not absolute), `delta` [-20, 20]

---

## Directory Structure

```
/kfs2/projects/scidac/scidac-data/gridkit-runs/

  # Track 2 (epistemic UQ) — completed runs
  hawaii-v5/                    # 16K samples, 4 gens, Gaussian std=12%
    meta.yml
    samples.csv
    runs/
      run_000/
        hawaii.json             # patched case (params + mon arrays)
        hawaii.solver.json
        mon.csv                 # simulator output
        stdout.txt
        stderr.txt
      run_001/
      ...
      run_000.parquet           # collected output (one per run)
      run_001.parquet
      ...
    slurm/
      chunk_00.sh  chunk_01.sh  ...  # sim sbatch scripts
      collect.sh                    # collect job (COLLECT_IN_SLURM=True)
      chunk_00_NNN.out  ...         # SLURM stdout/stderr

  illinois-v2/                  # 4K samples, 4 gens, Gaussian std=12%
    meta.yml
    samples.csv
    runs/
      run_000/  ...  run_000.parquet  ...
    slurm/

  # Track 1 (aleatoric UQ) — planned
  illinois-aleatoric-v1/
    scenario_cases/             # patched case JSONs from m_to_case
      scenario_0000.json  ...
    scenario_0000/
      run_000/  mon.csv  ...
      results.parquet
    ...
    summary.parquet
```

---

## Notebooks

| Notebook | Purpose | Status |
|----------|---------|--------|
| [`gridkit_helper.ipynb`](../notebooks/gridkit_helper.ipynb) | Single-case runner + solver override support; validate case before UQ sweep | Done |
| [`uq_setup.ipynb`](../notebooks/uq_setup.ipynb) | Epistemic UQ sweep (Track 2); SLURM parallel + collect; data quality check; sanity plots | Done (hawaii-v5, illinois-v2) |
| `uq_viz.ipynb` | Visualize ensemble results (time series, statistics, scatter) | TODO |
| `m_viz.ipynb` | Visualize PCM (.m) data | Exists |

---

## Open Questions

### Case correspondence (must resolve before any util work)

The `.m <-> .json` mapping has several subtleties that need to be studied on the **base raw case**
(`/kfs2/projects/scidac/scidac-data/ACTIVSg200/raw-tamu-data/case_ACTIVSg200.m` vs `illinois.json`)
before writing or modifying any utility functions:

1. **Off-generator question (critical)**: The illinois.json has 40 Genrou devices, while the .m
   has 49 generators (11 offline). 9 offline ng units are simply absent from the JSON. 2 offline
   units (buses 161, 197) ARE in the JSON with p0=q0=0. The question is: does GridKit require that
   a Genrou device be absent when offline, or can it be present with p0=q0=0? The existing JSON
   encodes two different answers for the same case. Need to understand the "art" behind which
   offline gens are kept vs dropped before deciding how `m_to_case` should handle them.

2. **Scenario hours where excluded gens may become online**: In the 8760 PCM run, some of the
   9 excluded generators (e.g. the Champaign ng units) may be dispatched online in some hours.
   Those hours cannot be correctly represented by just patching p0/q0. Need to understand the
   full extent of this: how many hours, which units, and whether those scenarios must be skipped,
   handled differently, or whether the JSON must be rebuilt for them.

3. **Near-zero dispatch**: Many generators in hour_1 have PG ~ -8e-7 (numerically zero from the
   solver). Clamping strategy TBD.

4. **Load variation**: `m_to_case` currently does not update LoadZIP parameters from `mpc.bus`
   PD/QD. Whether load variation is needed for Track 1 is TBD.

5. **Scenario filter**: The PCM dataset is `ACTIVSg200_wind_demand` (combined). No separate
   wind-only or demand-only `.m` files exist in `pcm-runs/`. If a simpler study is desired,
   subsampling by season or wind level is needed.

### Simulation setup (deferred)

6. **Solver JSON**: Which solver fn (e.g. fault scenario) should Track 1 use? Same as epistemic runs?

7. **Parallelism**: 8760 scenarios requires SLURM array on Kestrel, not in-notebook.

---

## Existing Assets to Reuse

From `gridkit_utils.py`:
- `generate_samples` (LHS, for Track 2)
- `make_run_dir` (patches case.json params, creates run_XXX/ dirs)
- `run_sample` (subprocess runner for DynamicSimulation)
- `collect_and_save` (aggregates mon.csv -> Parquet)
- `attach_json_ids` (MATPOWER <-> GridKit ID bridge, used internally by m_to_case)

From `m_viz_utils.py`:
- `.m` file parsing helpers (matpowercaseframes wrapper)
- `make_bus_trace`, `make_gen_trace` for input visualization

---

## Next Steps

1. **Study base raw case**: Compare `raw-tamu-data/case_ACTIVSg200.m` against `illinois.json`
   side by side — understand exactly how the JSON was built from the .m, especially the logic
   for which offline generators are included vs excluded. Document findings in `cases/illinois.md`.

2. **Audit 8760 hours for excluded-gen conflicts**: For each of the 9 excluded generators,
   check how many of the 8760 `.m` solution files have that generator online (PG > 0,
   GEN_STATUS=1). This determines whether scenario filtering is needed and at what scale.

3. **Prototype single-hour patch**: Test `m_to_case.py` on one representative `.m` solution
   (e.g. hour_1 or a peak-load hour); inspect patched JSON and run GridKit to confirm
   convergence before scaling.

4. **Decide on initial scale**: Start with a small subset (e.g. 24 or 168 hours) before
   committing to 8760.
4. Extend `collect_and_save` (or write a new helper) to add `scenario_id` column for
   combined aleatoric/epistemic results
5. Build `uq_setup.ipynb` (combined runner)
6. Build `uq_viz.ipynb` (ensemble plots, statistics, sensitivity)

---

## References

- MATPOWER column definitions: https://matpower.org/docs/
- GridKit case JSON format: `~/gridkit/build/examples/PhasorDynamics/Large/Illinois/illinois.json`
- PCM data: `/kfs2/projects/scidac/scidac-data/pcm-runs/ACTIVSg200_wind_demand/matpower/`
