# SCIDAC / GridKit UQ-Usecase: Contents Index

One place to find "what file explains X" and "what path is Y." This is a **committed,
git-tracked doc** (unlike `/memories/repo/`, which is personal and local — see
[`~/cmdlines/copilot/project_dev.md`](~/cmdlines/copilot/project_dev.md) for that
distinction). Update this file whenever a new work-note, case doc, or notebook is added.

---

## Main paths

| Path | What it is |
|---|---|
| `~/gridkit/uq-usecase/` | Project root for all UQ work |
| `~/gridkit/uq-usecase/work-notes/` | Design docs, investigation logs, the main plan |
| `~/gridkit/uq-usecase/cases/` | Per-grid reference docs (Illinois, Hawaii40) |
| `~/gridkit/uq-usecase/notebooks/` | All notebooks (viz, PF solving, UQ sweeps) |
| `~/gridkit/uq-usecase/py-utils/` | Shared Python helpers imported by the notebooks |
| `~/gridkit/uq-usecase/scripts/` | Build/rebuild/wip shell scripts |
| `~/gridkit/uq-usecase/pf-solver/`, `pm-solver/` | GridKit `solve_pf` and PowerModels.jl driver code |
| `~/gridkit/uq-usecase/kestrel_install.md` | Kestrel HPC environment setup / build instructions |
| `/kfs2/projects/scidac/scidac-data/gridkit-runs/` | All epistemic UQ run output (Parquet + SLURM logs) |
| `/kfs2/projects/scidac/scidac-data/pcm-runs/` | PCM hourly `.m` solution files (aleatoric track input) |
| `/nopt/nrel/apps/cpu_stack/software/gridkit/` | Installed GridKit dependencies on Kestrel |

---

## `work-notes/` — design docs and investigation logs

| File | What it covers |
|---|---|
| [`uq_plan.md`](work-notes/uq_plan.md) | **The main plan.** Aleatoric + epistemic track definitions, Task 0/0b, directory structure, notebook workflow, open questions. Has its own Contents + "Related investigative documents" section linking everything below. Start here. |
| [`pf_helper.md`](work-notes/pf_helper.md) | GridKit `solve_pf` vs PowerModels.jl power-flow solver comparison; why PM.jl is the production choice; safe-envelope findings |
| [`m_to_case_helper.md`](work-notes/m_to_case_helper.md) | MATPOWER `.m` → GridKit `case.json` mapping spec: `patch_case_from_m` / `build_case_from_solved_m`, field-by-field rules, test plan |
| [`notes_1.md`](work-notes/notes_1.md) | Early `m_viz` workflow summary (2026-06-17) |
| [`notes_2.md`](work-notes/notes_2.md) | Progress summary (2026-07-27): `m_viz`, `uq_setup` epistemic pipeline, PF solver work; executive summary |
| [`notes_3.md`](work-notes/notes_3.md) | Hawaii v6–v10 epistemic sensitivity run designs (H parameter, hypotheses, fleet tables) |
| [`notes_4.md`](work-notes/notes_4.md) | Full epistemic run inventory: sample counts, disk usage (GB), SLURM timing, bragging-point metrics |
| [`gt_ideas.md`](work-notes/gt_ideas.md) | Exploratory: graph transformers / ML on the UQ dataset (idea survey, not a fixed plan) |
| `pf_helper_old.md`, `pf_helper.md.bak` | Superseded drafts — **not current**, kept for history only, don't cite these |

## `cases/` — per-grid reference docs

| File | What it covers |
|---|---|
| [`illinois.md`](cases/illinois.md) | Illinois (ACTIVSg200, 200 bus): `.m` gen/bus tables, `illinois.json` field mapping, PF-solution patching guide, UQ parameter selection (v1/v2) |
| [`illinois_m_tables.md`](cases/illinois_m_tables.md) | Full parsed MATPOWER tables for Illinois (`mpc.bus`, `mpc.branch`, `mpc.gen`, `mpc.gencost`) |
| [`hawaii.md`](cases/hawaii.md) | Hawaii40 (37 bus): `.m` gen/bus/branch tables, `hawaii.json` field mapping, UQ parameter selection (v1–v5) |
| [`hawaii_m_tables.md`](cases/hawaii_m_tables.md) | Full parsed MATPOWER tables for Hawaii40 |
| [`hawaii_case_tables.md`](cases/hawaii_case_tables.md) | Hawaii `hawaii.json` device tables parsed directly from the JSON (Genrou params, H×mva ranking) |

## Root-level and setup docs

| File | What it covers |
|---|---|
| [`kestrel_install.md`](kestrel_install.md) | Kestrel HPC environment setup, module loads, build instructions |
| [`notebooks/setup.md`](notebooks/setup.md), [`setup_env.md`](notebooks/setup_env.md) | Notebook/conda environment setup |

## Notebooks (`notebooks/`)

| Notebook | Purpose |
|---|---|
| `m_viz.ipynb` | Interactive geographic grid visualization (any MATPOWER `.m` case) |
| `gridkit_viz.ipynb` | GridKit `case.json` device-level visualization (companion to `m_viz`) |
| `pf_helper.ipynb` | GridKit `solve_pf` experiments (companion to `work-notes/pf_helper.md`) |
| `pm_helper.ipynb` | PowerModels.jl experiments (companion to `work-notes/pf_helper.md`) |
| `m_to_case_helper.ipynb` | Validates the `.m` → `case.json` mapping (companion to `work-notes/m_to_case_helper.md`) |
| `gridkit_helper.ipynb` | Single-run validator for a case config before launching a full UQ sweep |
| `uq_setup.ipynb` | Full epistemic UQ sweep: sampling → SLURM run → collect → QC (the main production notebook) |
| `work.ipynb` | Scratch / in-progress notebook, not a stable reference |

## `py-utils/` — shared Python helpers

| File | What it provides |
|---|---|
| `gridkit_utils.py` | LHS/random sampling, run-dir generation, SLURM run/collect helpers, dispatch readers |
| `m_to_case.py` | `patch_case_from_m`, `build_case_from_solved_m`, batch helpers (documented in `m_to_case_helper.md`) |
| `m_viz_utils.py` | MATPOWER `.m` parsing, geo attachment, Plotly map building |
| `pf_utils.py` | GridKit `solve_pf` wrapper, perturbed-`.m` generators |
| `pm_utils.py` | PowerModels.jl (`pm_solve.jl`) wrapper |
| `gen_hawaii_tables.py` | Generates `cases/hawaii_*_tables.md` from source data |

## `scripts/` — build and workflow scripts

| Script | Purpose |
|---|---|
| `1_build_llvm.sh` … `6_build_gridkit.sh` | Full dependency-stack build, in order (rarely needed) |
| `rebuild_if_updated.sh` | Checks for upstream changes, rebuilds only if needed (the usual daily command) |
| `wip.sh` | Commits work-in-progress to the `isatkaus/uq-usecase` branch |
| `build_mkdocs.sh` | Builds the mkdocs site from `work-notes/` + `cases/` + `notebooks/` |
| `create_h-py312-basic.sh` | Creates the `h-py312-basic` conda environment |

---

## Relationship to `/memories/repo/`

This file is the **git-tracked, always-visible** index — anyone who clones the repo sees it.
`/memories/repo/scidac_shared_context.md` (see `project_dev.md` for its real path) is a
**personal, local-only** summary the agent maintains for cross-session continuity (current
run status, recent decisions). They're complementary, not duplicates: this file answers
"where is the doc that explains X," the memory file answers "what did we decide / where did
we leave off."
