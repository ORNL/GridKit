# Notebooks: Setup and Usage

## Conda Environment

The notebooks were developed and tested with a Python 3.12 conda environment with the following packages.

### Required packages

| Channel | Packages |
|---|---|
| conda | `python=3.12 numpy pandas pyyaml jupyter plotly h5py poppler pillow scipy pyarrow` |
| pip | `tables matpowercaseframes==2.1.0` |

Create it manually (any system, adjust `--prefix` or `--name` to your setup):
```bash
conda create --prefix ./my-env python=3.12 numpy pandas pyyaml jupyter plotly h5py poppler pillow scipy pyarrow --yes
conda activate ./my-env
pip install tables matpowercaseframes==2.1.0
```

### Build with the provided script (Kestrel HPC only)

A convenience script is provided for Kestrel at `../scripts/create_h-py312-basic.sh`.
It is **Kestrel-specific**: the `envs_dir` path and module commands (`ml mamba`, `ml conda`)
must be adjusted for other systems. Edit `envs_dir` near the top of the script before running.

Build the environment (run on a compute node)

```bash
salloc --time=1:00:00 --account=scidac --nodes=1 --partition=debug
bash ~/gridkit/uq-usecase/scripts/create_h-py312-basic.sh
```

The environment installs to `~/flash-w-dom/conda-envs/h-py312-basic/`
(flash filesystem with 256K DoM for fast conda env access).

## Notebooks

| Notebook | Description |
|---|---|
| `gridkit_helper.ipynb` | GridKit run setup and result parsing helpers |
| `gridkit_viz.ipynb` | GridKit result visualization |
| `m_viz.ipynb` | Geovisualization of MATPOWER `.m` cases with optional GridKit JSON ID augmentation |

### m_viz notebook

`m_viz.ipynb` parses MATPOWER `.m` case files (ACOPF or PF solutions) into normalized pandas
DataFrames (bus, gen, branch) and produces interactive Plotly geo maps. Geographic coordinates
come from a `.gic` or `.AUX` file. Generator fuel labels and multi-circuit branch data are read
directly from the `.m` file.

The map visualizes the solution stored in the `.m` file: bus loads (PD), dispatched generation
(marker size proportional to MW, fuel-colored), and branch loading percentage (viridis color
scale, requires PF/QF/PT/QT/RATE_A columns).

Optionally, set `GRIDKIT_REPO` to the root of a local GridKit repository clone to augment hover
labels with GridKit-assigned bus and branch IDs. JSON case files are read from
`GRIDKIT_REPO/examples/PhasorDynamics/`. Set `GRIDKIT_REPO = None` to skip this augmentation.

Key processing steps:

- **Colocated bus splitting**: buses sharing identical lat/lon are spread on a small circle so
  each bus is individually hoverable on the map.
- **Generator fan-out**: when a bus has multiple generators, their markers are spread radially
  so each is individually hoverable. Thin connector lines link markers back to the bus (visual
  only). Optionally, single-generator buses can also be offset (`FAN_SINGLE_GEN = True`) to
  keep the bus marker visible underneath.
- **Fault bus overlay**: optionally highlights a faulted bus with a red marker.
- **Branch loading coloring**: branches colored by loading percentage (viridis) when flow data
  is available.

Currently supported cases: **Hawaii40**, **Illinois (ACTIVSg200)**, **Texas (ACTIVSg2000)**,
**WECC (ACTIVSg10k)**.

## Supporting Files

- `../py-utils/gridkit_utils.py` — Python utilities for GridKit runs
- `../py-utils/m_viz_utils.py` — Python utilities for MATPOWER `.m` case geovisualization (`m_viz` notebook)
- `../context/` — input context data files used by notebooks
- `../figs/` — figures (tracked on branch, excluded from PRs)
- `../mds/` — markdown/doc assets
