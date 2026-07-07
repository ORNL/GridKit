# Notebooks: Setup and Usage

## Conda Environment

The notebooks require the `h-py312-basic` conda environment.
The creation script is in `../scripts/create_h-py312-basic.sh`.

### Build the environment (run on a compute node)

```bash
salloc --time=1:00:00 --account=scidac --nodes=1 --partition=debug
bash ~/gridkit/uq-usecase/scripts/create_h-py312-basic.sh
```

The environment installs to `~/flash-w-dom/conda-envs/h-py312-basic/`
(flash filesystem with 256K DoM for fast conda env access).

### Register the Jupyter kernel (one-time, run once after building the env)

```bash
# Register h-py312-basic as a persistent Jupyter kernel so that it survives VSCode window reloads
conda run -n h-py312-basic python -m ipykernel install --user --name h-py312-basic --display-name "h-py312-basic"
```

After registration, the kernel appears permanently in VSCode's kernel picker
without needing to re-add it after each window reload.

## Notebooks

| Notebook | Description |
|---|---|
| `gridkit_helper.ipynb` | GridKit run setup and result parsing helpers |
| `gridkit_viz.ipynb` | GridKit result visualization |
| `m_viz.ipynb` | Monitoring/signal visualization |

## Supporting Files

- `../py-utils/gridkit_utils.py` — Python utilities for GridKit runs
- `../py-utils/m_viz_utils.py` — Python utilities for monitoring visualization
- `../context/` — input context data files used by notebooks
- `../figs/` — figures (tracked on branch, excluded from PRs)
- `../mds/` — markdown/doc assets
