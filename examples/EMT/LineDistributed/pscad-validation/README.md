# LineDistributed against PSCAD

Validation of the `LineDistributed` example (adaptive implicit integration that may
step over the line delay) against PSCAD 5.1 from the same source data: tower
geometries and conductor catalogs. GridWorkbench computes the line constants and
this folder fits them for GridKit; PSCAD receives the same tower and fits its own
frequency-dependent phase model. Both tools integrate identical forcing: the
example's finite C2 transitions are reproduced in PSCAD by a scripted component.

## Layout

| Path | Contents |
| --- | --- |
| `lines/*.line.json`, `lines/catalog.json` | Source data: four single-circuit lines (69 kV 16 km, 138 kV 40 km, 345 kV 100 km twin bundle, 765 kV 322 km quad bundle) and their conductor/tower catalog |
| `lines/<line>/` | Generated: `constants.csv` (modal z, y), `sequence.json` (60 Hz data, surge impedances, delays), `tower.json` (PSCAD tower), `fd.fit.json` (GridKit input), `bergeron.fit.json`, `fit-report.json`, `fits.png`, `constants-comparison.{png,json}` (PSCAD LCP against GridWorkbench), `pscad-fit.json` (PSCAD's own fit summary) |
| `scenarios/scenarios.json` | Every run object: energize (open end), pulse (5 us rise), fault (receiving end, 2 ms clearing ramp), sag |
| `scenarios/gridkit/<line>__<model>.json` | `line_distributed` suites |
| `scenarios/pscad/*.pscx` | One PSCAD project per run, plus half-step `_h` refinements |
| `results/` (ignored, not kept) | Written by the runs: `gridkit/<line>__<model>/<scenario>/` with `samples.csv`, `steps.csv` (accepted steps, h/tau, history error, BDF order) and `statistics.json`; `pscad/<project>/` with the CSV on the study clock, `.psout`, PSCAD messages and line-constants outputs |
| `plots/<line>/<model>_<scenario>.png`, `plots/summary.png` | Time-domain comparisons and the error heatmap |
| `plots/steps/` | Accepted step over the propagation delay against time, per run |
| `validation.md`, `validation.json` | Error table (peak and RMS, normalized by the PSCAD signal) and step statistics |
| `diagnostics/` | Step-size study: line variants (2 to 200 km, 4x losses, 1000 ohm m soil) with R-L and resistive terminations (`step_study.png/json`), suites for the tolerance, order, pole-cap and algebraic-suppression experiments (`runs_*/`), PSCAD runtime cases (`pscad/S_*.pscx`) and the runtime summary `timing.json`; their raw output (`results*/`, `pscad-results/`) is ignored |
| `lib/` | `linedata.py` (constants, modes, delays), `vectorfit.py`, `scenarios.py`, `pscad_circuit.py` (PSCAD project through GridWorkbench's authoring internals), `pscad_run.py` |

## Reproduce

GridWorkbench environment (Windows, PSCAD 5.1, mhi-pscad, mhi-psout) for the
Python steps; the gridkit-sdk container with the IDAS step-check patch for GridKit.

```sh
python make_lines.py                 # lines/<line>/...
python make_scenarios.py             # suites and PSCAD projects (--gridkit-only skips PSCAD)
python run_pscad.py                  # results/pscad
docker exec <container> bash -lc 'cd /work/examples/EMT/LineDistributed/pscad-validation && \
  bash run_gridkit.sh /home/gridkit/build/examples/EMT/LineDistributed/line_distributed'
python compare.py                    # plots, validation.md
python compare_constants.py          # lines/<line>/constants-comparison.*
python plot_steps.py                 # plots/steps
```

The container must mount this worktree at `/work`; `run_gridkit.sh` writes beside the
suites. The GridKit build used here was the example's targets against SUNDIALS 7.8.0
with `buildsystem/spack_repo/gridkit/packages/sundials/idas-step-check.patch`.

## State of the results

- Raw simulation output is not kept in the repository (see `.gitignore`); `plots/`,
  `validation.*` and the diagnostics summaries are what remains of each run, and a
  rerun of the commands above reproduces them.
- `plots/` and `validation.*` are the first full pass: frequency-dependent model, all
  four lines, four scenarios, with the line data produced by averaging the
  capacitance matrix under ideal transposition. PSCAD averages the Maxwell
  potential-coefficient matrix instead; `lib/linedata.py` has since been switched to
  that convention (it reproduces PSCAD's capacitances to 0.000%), but `lines/`,
  `scenarios/` and the plots have not been regenerated under it. On the flat bundled
  towers the old convention puts the aerial velocity 0.4% below PSCAD's, which is
  the arrival-time offset visible in the 345 kV and 765 kV pulse and energization
  plots. The short lines agree to 1 to 2% peak on voltages and better than 0.2% on
  the receiving-end current. Rerunning the four commands above refreshes everything
  under the PSCAD convention.
- `diagnostics/` was produced after that switch.
- Lossless Bergeron runs exist for the record but are not reported.

## Findings about the implementation

- The adaptive step exceeds the propagation delay only with a resistive load and
  lines short enough that the delay is below the solver's steady-state step: at
  rel_tol 1e-6 that step is 30 to 50 us for every length from 2 to 200 km, so
  h/tau > 1 for most of the time at 2 to 8 km and 22% of the time at 16 km. A series
  R-L load holds every length at 1 to 6 us steps because both terminations reflect
  high frequencies fully and only the line losses damp them; longer lines or higher
  losses help gradually but never reach the resistive case. The example's own 80 km
  R-L fault case behaves the same way.
- The steady-state step ceiling is not set by BDF order, history tolerance,
  compression or fast fit poles; the worst-weighted local errors sit on the wave
  variables fed by the delayed signals, consistent with the C1 cubic history putting
  second-derivative kinks into the delayed input. Excluding algebraic variables from
  the error test raises the ceiling modestly.
- Without history compression, fronts recirculating between two near-reflecting
  terminations (open end, lossless line) spawn a history segment per accepted step
  each round trip and the step collapses; all suites therefore run with
  `compress: true`. Clearing a fault with a 50 us conductance ramp forces the source
  inductance current to zero and produces megavolt spikes in both tools, so the
  fault scenarios clear over 2 ms.
- PSCAD runtime for the step-over-delay cases: GridKit 0.11 to 0.14 s per 100 ms
  (0.04 to 0.06 s at rel_tol 1e-4) against 0.9 to 1.0 s for PSCAD at the 2 us step
  that matches GridKit's accuracy; PSCAD's coarsest admissible steps are faster but
  3 to 13% off its own fine solution.

## Uncommitted source changes that go with this folder

- `examples/EMT/LineDistributed/LineCircuit.cpp`: current tolerance scale includes
  the characteristic admittance (open-ended lines no longer get a 1e-10 A absolute
  tolerance on wave currents); voltage scale includes pulse amplitudes.
- `Input.hpp`, `main.cpp`: optional solver fields `max_order` and
  `suppress_algebraic`; `steps.csv` gains a BDF `order` column; `errors.csv` lists the
  variables dominating IDAS's weighted local error every 50 steps.
- `LineDistributed.hpp`, `LineCircuit.hpp`: variable names for those diagnostics.
- `GridKit/Solver/Dynamic/Ida.{hpp,cpp}`: `getLastOrder`, `getEstimatedLocalErrors`,
  `getErrorWeights` accessors.
