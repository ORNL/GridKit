# Work Notes

## m_viz workflow summary (2026-06-17)

### What `m_viz.ipynb` does

- Loads MATPOWER `.m` cases into normalized tables (`bus`, `gen`, `branch`).
- Runs quick integrity checks and previews.
- Resolves geo source files with fallback logic:
  - prefer `.gic` when available,
  - otherwise use `.AUX`.
- Builds first map visualization with buses, branches, and generators.

### Utilities added in `m_viz_utils.py`

- MATPOWER ingestion and normalization:
  - `read_matpower_case`
  - `normalize_bus_df`
  - `normalize_gen_df`
  - `normalize_branch_df`
  - `summarize_case`
  - `validate_case_tables`
  - `validation_report_df`
- Geo parsing:
  - `create_bus_geo_df_from_gic`
  - `create_bus_geo_df_from_aux`
  - `create_bus_geo_df` (dispatches by file type)
- Geo attachment and preprocessing:
  - `attach_geo_to_case`
  - `split_points` (bus split for colocated buses)
  - `split_generator_points` (generator fan-out around bus)
- Plotting helpers:
  - `plot_first_geo`
  - `_build_branch_line_arrays`
  - `_build_gen_connector_line_arrays`
- Branch metrics:
  - `compute_branch_loading`

### Geo behavior differences: Hawaii vs ACTIVSg

#### Hawaii (`Hawaii40_20231026`)

- Typical input has no matching `.gic`; notebook prints a warning and continues with `.AUX`.
- Geo extraction comes from AUX `Substation` + `Bus` table mapping.
- Generator overlap is handled with Hawaii-scoped fan-out in the notebook:
  - generators on the same bus are displayed radially around bus center,
  - optional connector overlay is shown as `fake gen connections`.
- Goal: preserve generator visibility and hoverability when many units share one bus.

#### ACTIVSg cases

- Usually have `.gic` available and use direct GIC-based bus coordinates.
- Existing bus `split_points` behavior remains the primary de-overlap mechanism.
- No Hawaii-only generator fan-out toggles are applied by default.
- Result: closer to original ACTIVSg plotting style while retaining current first-geo map flow.

### Notes

- Generator fan-out and connector lines are visualization-only and do not modify electrical topology.
- Branch lines remain the real network branches; `fake gen connections` are just aid lines for display.
