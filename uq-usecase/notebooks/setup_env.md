# Conda environment setup

The notebooks were developed and tested with Python 3.12. Other versions (3.11+) should also work.

More notebooks will be added in future PRs. The environment below covers all currently included
notebooks.

## Packages

| Channel | Packages |
|---|---|
| conda | `python=3.12 numpy pandas plotly scipy jupyter` |
| pip | `matpowercaseframes==2.1.0` |

```bash
conda create --prefix ./my-env python=3.12 numpy pandas plotly scipy jupyter --yes
conda activate ./my-env
pip install matpowercaseframes==2.1.0
```

## Notebooks

### m_viz

`m_viz.ipynb` parses MATPOWER `.m` case files (ACOPF or PF solutions) into normalized pandas
DataFrames and produces interactive Plotly geo maps. Geographic coordinates come from a `.gic`
or `.AUX` file. The map shows bus loads (PD), dispatched generation (fuel-colored, sized by MW),
and branch loading percentage (viridis color scale).

Optionally, set `GRIDKIT_REPO` to the root of a local GridKit repository clone to augment hover
labels with GridKit-assigned bus and branch IDs (`gridkit_utils.py`, requires `scipy`). Set
`GRIDKIT_REPO = None` to skip this.

Supported cases: **Hawaii40**, **Illinois (ACTIVSg200)**, **Texas (ACTIVSg2000)**, **WECC (ACTIVSg10k)**.

#### Input data

Download TAMU synthetic grid cases from
[Texas A&M Electric Grid Test Cases](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/).
After unzipping, set `CASE_DATA_DIR` to the folder containing the extracted files and `CASE_NAME`
to match. Everything else is auto-detected.

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

#### Supporting files

| File | Purpose |
|---|---|
| `py-utils/m_viz_utils.py` | MATPOWER `.m` parsing, geo merge, colocated bus splitting, generator fan-out, plotting |
| `py-utils/gridkit_utils.py` | GridKit JSON ID augmentation for hover labels (optional) |
