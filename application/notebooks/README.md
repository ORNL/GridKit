# application/notebooks

Python notebooks and supporting packages for exploratory analysis and
visualization of GridKit cases. Not part of the C++ build; enabled via the
`GridKit_ENABLE_NOTEBOOKS` CMake option (see below) or usable standalone
with `pip`.

## Packages

| Package | Purpose |
|---|---|
| `m_viz/` | MATPOWER `.m` case parsing, geo-coordinate attachment, and Plotly geo plotting. Used by `m_viz.ipynb`. |
| `gridkit_common/` | Small helpers shared across notebooks (linking MATPOWER tables back to GridKit JSON case IDs). |
| `uq_sampling/` | Parameter sampling (Latin Hypercube / random), run-directory orchestration, and Genrou dispatch patching for uncertainty-quantification sweeps. Not yet exercised by a notebook in this repo — added ahead of the UQ sweep notebook that will use it. |

## Environment setup

Developed and tested with Python 3.12; 3.11+ should also work.

```sh
conda create --prefix ./my-env python=3.12 --yes
conda activate ./my-env
pip install -e .[notebook,test]
```

Or with a plain venv:

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[notebook,test]
```

This installs `m_viz`, `gridkit_common`, and `uq_sampling` as regular
importable packages — no `sys.path` manipulation needed in notebooks.

## Notebooks

### m_viz

`m_viz.ipynb` parses MATPOWER `.m` case files (ACOPF or PF solutions) into normalized pandas
DataFrames and produces interactive Plotly geo maps. Geographic coordinates come from a `.gic`
or `.AUX` file. The map shows bus loads (PD), dispatched generation (fuel-colored, sized by MW),
and branch loading percentage (viridis color scale).

Optionally, set `GRIDKIT_REPO` to the root of a local GridKit repository clone to augment hover
labels with GridKit-assigned bus and branch IDs (`gridkit_common.ids`, requires `scipy`). Set
`GRIDKIT_REPO = None` to skip this.

Supported cases: **Hawaii40**, **Illinois (ACTIVSg200)**, **Texas (ACTIVSg2000)**, **WECC (ACTIVSg10k)**.

#### Input data

Download TAMU synthetic grid cases from
[Texas A&M Electric Grid Test Cases](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/).
After unzipping, set `CASE_DATA_DIR` (or the `GRIDKIT_CASE_DATA_DIR` environment variable) to the
folder containing the extracted files and `CASE_NAME` to match. Everything else is auto-detected.

##### Hawaii40

| File | Used for |
|---|---|
| `Hawaii40_20231026.m` | MATPOWER case (buses, generators, branches) |
| `Hawaii40_GIC_data.gic` | Geographic coordinates (preferred) |
| `Hawaii40_20231026.AUX` | Geographic coordinates (fallback if no `.gic`) |

##### Illinois (ACTIVSg200)

| File | Used for |
|---|---|
| `case_ACTIVSg200.m` | MATPOWER case (buses, generators, branches) |
| `ACTIVSg200_GIC_data.gic` | Geographic coordinates (preferred) |
| `ACTIVSg200.AUX` | Geographic coordinates (fallback if no `.gic`) |

##### Texas (ACTIVSg2000)

| File | Used for |
|---|---|
| `case_ACTIVSg2000.m` | MATPOWER case (buses, generators, branches) |
| `ACTIVSg2000_GIC_data.gic` | Geographic coordinates (preferred) |
| `ACTIVSg2000_dynamics.AUX` | Geographic coordinates (fallback if no `.gic`) |

##### WECC (ACTIVSg10k)

| File | Used for |
|---|---|
| `case_ACTIVSg10k.m` | MATPOWER case (buses, generators, branches) |
| `ACTIVSg10k_GIC_data.gic` | Geographic coordinates |

More notebooks will be added in future PRs.
