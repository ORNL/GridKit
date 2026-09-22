# Work Notes — Progress Summary (2026-07-27)

## `m_viz.ipynb` — MATPOWER case geovisualization

**Summary**: Interactive notebook that loads MATPOWER `.m` case files and produces a
geographic map of the grid — buses, branches, and generators — with hover details and
optional overlays. Tested on Hawaii40, ACTIVSg200 (Illinois), ACTIVSg2000 (Texas), and
ACTIVSg10k (WECC); built generically against the MATPOWER `.m` format so it should work
for other TAMU synthetic grid cases with a `.gic`/`.AUX` geo file, without code changes.

**Basic functionality:**
- Parses `.m` files into normalized `bus` / `gen` / `branch` tables.
- Loads geographic coordinates from `.gic` or `.AUX` files.
- Plots an interactive Plotly map: bus markers sized by load, generator markers sized by
  dispatch and colored by fuel type, branches colored by loading percentage.
- Exports to standalone HTML for sharing or embedding (e.g. in MkDocs).

**Key features:**
- **Branch loading colormap** — viridis-colored by % of RATE_A, computed from PF/QF/PT/QT.
- **Colocated bus splitting** — buses sharing identical coordinates are spread on a small
  circle so each remains individually hoverable.
- **Generator fan-out** — multiple generators at one bus are spread radially with connector
  lines back to the bus (handles Hawaii's 10-generator buses, etc.).
- **Fault bus overlay** — highlight a bus with a marker, by MATPOWER bus number, bus name,
  or GridKit JSON device ID.
- **GridKit ID hover augmentation** — cross-references a GridKit `case.json` so hover
  tooltips show GridKit device IDs alongside MATPOWER identifiers.

**Utilities (`m_viz_utils.py`, `gridkit_utils.py`):** `read_matpower_case`, `summarize_case`,
`validation_report_df`, `attach_geo_to_case`, `plot_grid`, `lookup_fault_bus`,
`attach_json_ids`, `compute_branch_loading`.

**Next steps:**
- Extend the map to visualize GridKit `case.json` dynamic device models directly —
  Genrou, SexsPti, Ieeet1, Tgov1, and other device parameters overlaid on the same
  geographic layout, not just the MATPOWER `.m` static fields.
- Add dynamic simulation result visualization: plot `mon.csv` / collected Parquet
  time-series (voltage, angle, frequency, rotor angle) directly on the map or as linked
  time-series views, to connect UQ run outputs back to grid geography.

---

## `uq_setup.ipynb` — sampling-based UQ workflow (epistemic track)

**Summary**: End-to-end pipeline to run large ensembles of GridKit dynamic simulations
with sampled generator parameter uncertainty (currently: inertia constant H), parallelized
on Kestrel via SLURM, with automated result collection and QC.

**Workflow stages:**
1. Sample parameter sets via Latin Hypercube Sampling (uniform or Gaussian marginals).
2. Generate one patched `case.json` + `solver.json` per sample (`run_NNN/` directories).
3. Execute simulations in parallel across SLURM nodes (array of sbatch scripts, one per
   node, workers fanned out with `xargs`).
4. Collect per-run outputs into Parquet files (parallelized reader, safe on Lustre).
5. Automated data quality checks (file counts, NaN scan, physical bounds) and sanity plots.

**Companion notebook**: `gridkit_helper.ipynb` — single-run validator for testing a case
configuration (including solver overrides like fault timing) before launching a full sweep.

**Parallelization**: two levels — across SLURM nodes (one sbatch script per node, each
covering a slice of the sample indices) and within a node (multiple worker processes per
node via `xargs -P`), so the sample budget scales with both node count and per-node core
count.

**Two result serialization modes** — chosen per run:
- **Stacked** (`results.parquet`) — all runs combined into a single file; simplest to load,
  used for small runs (e.g. the 20-sample Hawaii v1 run).
- **Per-run** (`runs/run_NNN.parquet`) — one file per sample; scales to large N without a
  single huge in-memory concat step; used for all runs from ~1K samples upward.

**Sample of completed runs** (`/kfs2/projects/scidac/scidac-data/gridkit-runs/`; full run
history recorded in `cases/hawaii.md` and `cases/illinois.md`):

| Case | Run | Samples | Distribution | Generators sampled | Serialization |
|---|---|---:|---|---|---|
| Hawaii40 | v2 | 1,000 | uniform ±10% | 4 | per-run |
| Hawaii40 | **v5** | **16,000** | Gaussian σ=12% | 4 (hop distances 1–3+ from fault) | per-run |
| Illinois (ACTIVSg200) | **v2** | **4,000** | Gaussian σ=12% | 4 (hop distances 4–9 from fault) | per-run |

Both uniform and Gaussian sampling have been exercised end-to-end; Gaussian σ=12% is the
current standard for production epistemic runs on both grids. These represent initial
epistemic UQ runs (~22,000 dynamic simulations total) sampling inertia (H) uncertainty on
a fixed subset of generators per grid. Next planned expansions: broaden parameter
coverage to additional generators/device parameters beyond H, and sweep fault location
across multiple buses rather than a single fixed fault point per case.

**Utilities (`gridkit_utils.py`):** `generate_samples`, `make_run_dir`, `run_sample`,
`collect_parallel`, `collect_and_save`, `attach_json_ids`.


---

## PF solver work (`solve_pf`) and MATPOWER↔GridKit case tooling

**Summary**: A thin standalone wrapper around GridKit's existing Newton/KINSOL power-flow
solver, exposed as a CLI tool (`solve_pf`) that reads a MATPOWER `.m` file directly and
produces a converged bus voltage/angle solution. No new solver was written — the work was
(1) wrapping GridKit's existing PF machinery for standalone `.m`-file use, (2) adjusting
the local MATPOWER parser to handle unit conversions and ACTIVSg-extended fields correctly,
and (3) running a series of perturbation experiments against it. Purpose: generate
realistic operating-point scenarios for GridKit dynamic simulation.

**Basic functionality:**
- Reads `.m` files, warm-starts from stored Vm/Va, solves standard AC power flow using
  GridKit's existing solver internals.
- Can optionally write back a solved `.m` file with updated Vm/Va (topology unchanged).
- Verified convergence on 3-bus, Hawaii40 (37-bus), and ACTIVSg200 (200-bus) cases.

**Parser fixes made along the way**: two unit-conversion bugs in the local MATPOWER
parser were found and fixed — bus angles stored in degrees needed conversion to radians,
and generator/load powers stored in MW/MVAr needed conversion to per-unit.

**Supporting Python utilities (`pf_utils.py`):** run/parse the solver, diff two solutions,
and generate perturbed `.m` scenarios (load scaling, wind curtailment, generators offline)
for stress-testing the solver.

**Key finding**: across all tested perturbations (load ±80%, wind curtailment, up to
10 generators offline), power-flow solutions respond almost entirely through bus **angle**
— voltage magnitude barely moves (< 0.01 pu even under the largest perturbation tested).
This is explained by the fraction of voltage-regulated (PV/slack) buses — about 24% in
both Illinois and Hawaii — which is consistent with typical real transmission networks.
Practical implication: for UQ studies built on this solver, bus angle and branch flow are
the more informative quantities of interest; voltage magnitude is comparatively insensitive
to steady-state operating point.

**`m_to_case.py`** — converts a solved `.m` file into a patched GridKit `case.json`:
updates bus voltage initial conditions and generator dispatch (P/Q) to match a given
power-flow solution. This is the bridge that will let many operating-point scenarios
(e.g. 8760 hourly PCM solutions) each become a valid GridKit simulation starting point.

**In progress / planned**: extending this pipeline to generate a full year of hourly
operating-point scenarios (aleatoric UQ track), including load-model patching and a
small prototype study to confirm that fault-response outputs (rotor angle swings,
frequency nadir) actually vary meaningfully across operating points before committing to
the full-scale run.

---

## Executive Summary (updated 2026-08-26)

**1. Grid visualization (`m_viz`)**
- Interactive geographic map of any MATPOWER `.m` case — buses, branches, generators,
  branch loading — tested on 4 TAMU grids from 37 to ~10,000 buses.
- Cross-references GridKit `case.json` device IDs directly in the map hover tooltips.
- Next: add dynamic device model (Genrou, SexsPti, Tgov1) and simulation result viz.

**2. Epistemic UQ pipeline (`uq_setup`)**
- HPC-parallelized (across nodes + within node): sampling → simulation → collection → QC.
- ~43,000 epistemic dynamic simulations completed across 2 grids (Hawaii40, Illinois).
- Hawaii: systematic hypothesis-driven experiments (v6–v10) varying 4 generators each
  with N=4000 LHS samples, targeting H-sensitivity detectability vs machine size,
  distance, and dispatch. Plus a 39-parameter Morris sensitivity design (N=1000)
  sampling all fleet H parameters simultaneously.
- Next: analyze v6–v10 results for sensitivity ranking, then sweep fault location.

**Sample of completed runs** (`/kfs2/projects/scidac/scidac-data/gridkit-runs/`):

| Case | Run | Samples | Distribution | Gens sampled | Hypothesis / purpose |
|---|---|---:|---|---:|---|
| Hawaii40 | v2 | 1,000 | uniform ±10% | 4 | initial bus-diverse set |
| Hawaii40 | v5 | 16,000 | Gaussian σ=12% | 4 | production epistemic (hop distances 1–3+) |
| Hawaii40 | v6 | 4,000 | uniform ±20% | 4 | smallest H |
| Hawaii40 | v7 | 4,000 | uniform ±20% | 4 | farthest from fault |
| Hawaii40 | v8 | 4,000 | uniform ±20% | 4 | smallest mva + dispatch |
| Hawaii40 | v9 | 4,000 | uniform ±20% | 4 | largest H×mva (positive control) |
| Hawaii40 | v10 | 4,000 | uniform ±20% | 4 | small p0, moderate H×mva |
| Hawaii40 | ian-csv | 1,000 | Morris | 39 | all fleet H, sensitivity |
| Illinois | v2 | 4,000 | Gaussian σ=12% | 4 | hop distances 4–9 from fault |

See [`notes_3.md`](notes_3.md) for v6–v10 experiment designs and fleet analysis.

**3. Power-flow tooling (GridKit `solve_pf` + PM.jl)**
- Two PF solvers tested on ACTIVSg200: GridKit `solve_pf` (Newton/KINSOL) and PM.jl
  (Ipopt, industry-standard).
- **PM.jl is the production PF solver** for all case.json generation: correct PV→PQ
  switching, writes post-solve PG/QG, reproduces TAMU/PowerWorld reference to 1e-5 pu.
- GridKit `solve_pf` has no PV→PQ switching; ~0.030 pu voltage bias vs PM.jl. Usable
  only inside a "safe envelope" (±80% per-bus load, N≤2 gens off) for relative angle
  comparisons. See [`pf_helper.md`](pf_helper.md) Sections 1 and 11.
- Pluggable solver interface designed: `run_pf_solve(solver="pm"|"gridkit")` dispatcher
  in `pf_utils.py`. GridKit swap requires adding post-solve PG/QG writeback to
  `solve_pf.cpp` (missing output feature, not missing solve capability).
- `.m` → `case.json` bridge (`m_to_case.py`) implemented and documented: handles bus
  init, gen dispatch, LoadZIP patching, offline-gen removal. **TODO**: run full
  notebook test suite (`m_to_case_helper.ipynb` has 4 known bugs to fix first).
  See [`m_to_case_helper.md`](m_to_case_helper.md).
