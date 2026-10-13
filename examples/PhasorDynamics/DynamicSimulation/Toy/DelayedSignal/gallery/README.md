# Signal delay gallery

Open [index.html](index.html). All seven figures are PNGs in `figures/`.
Every waveform comparison includes the actual input signal.

| Figure | What changes |
| --- | --- |
| 01 — smooth ramp | Delays from 10 ms to 1.5 s |
| 02 — pulse widths | Gaussian width from 25 to 400 ms, with a fixed 250 ms delay |
| 03 — multiple timescales | A broad pulse followed by a narrow pulse |
| 04 — oscillatory packet | A smooth envelope containing oscillations |
| 05 — adaptive pulse | Delays of 1, 100, and 500 ms; accepted step sizes and step/delay ratios |
| 06 — long history | A 120 ms feature delayed by up to 8 s |
| 07 — compression | Input, both delayed outputs, history counts, and numerical errors |

These are 20 runs of the actual GridKit `SystemModel`, `SignalDelay`,
`SignalHistory`, and IDA wrapper linked from `build/`. The small driver
[gallery.cpp](gallery.cpp) supplies prescribed smooth inputs, evaluated at the
queried time like `FunctionSignalSource`, and records diagnostics. It contains
no separate delay implementation.

`h` is an accepted internal solver step. The monitor interval is separate.
All runs use relative tolerance `1e-6`, absolute tolerance `1e-8`, `h_min = 1e-6` s, and analytic
prehistory. References evaluate the input function at `t - tau`. Gaussian
widths are scale parameters in `exp(-((t-center)/width)^2)`, not full widths
at half maximum. Compression is disabled except in its explicit comparison.
Window pruning remains enabled in every run.

## Data and interpretation

- `runs.json` contains all experiment settings.
- `summary.csv` and `summary.json` contain errors, accepted steps, history
  rejections, native error/nonlinear failures, and retained segment counts.
- `data.zip` contains the complete CSVs and logs: `*_samples.csv` at uniform
  monitor times, `*_midpoints.csv` at every accepted step midpoint,
  `*_steps.csv` at accepted endpoints, and `*_attempts.csv` for history checks.
- `provenance.json` identifies the source revision and local source diff.
- `performance/` contains separate `DynamicSimulation` application benchmarks,
  their inputs, and before/after timings.

Errors are sampled maxima, not continuous error bounds. Generation checks for
nonfinite data, inconsistent step sizes, failed acceptance checks, and missed
features. Sparse samples cannot guarantee detection of every arbitrarily narrow
smooth input. Source resolution still matters even though there is no shortest-delay
step cap. Every crossed history polynomial is checked against the candidate IDAS dense
polynomial, including features between its sample times. `h_min` is a failure
floor, not a source sampling guarantee. Rejected trials use native IDAS rollback,
which preserves the BDF history and order; curves show actual accepted steps,
without smoothing.

Compression can change the adaptive step sequence. Segment counts are not total
process memory. History rejection counts exclude IDA's internal error-test and
nonlinear-solve failures.

## Reproduce

From the repository root, with the configured CMake Makefiles build, the
[IDAS step-check extension](../../../../../../buildsystem/spack_repo/gridkit/packages/sundials/README.md),
and NumPy/Matplotlib:

```bash
cmake --build build --target DynamicSimulation test_pd_signal_delay -j 10
python3 -B examples/PhasorDynamics/DynamicSimulation/Toy/DelayedSignal/gallery/generate.py
```

Compilation and raw intermediates use a temporary directory. Results and PNGs
are written here. To redraw the saved results without another simulation:

```bash
python3 -B examples/PhasorDynamics/DynamicSimulation/Toy/DelayedSignal/gallery/plot_gallery.py
```

To refresh the application timings with the current build:

```bash
python3 examples/PhasorDynamics/DynamicSimulation/Toy/DelayedSignal/gallery/performance/measure.py
```

The benchmark cases use one 1 ms tap (`short`), 17 taps from 1 ms to 4 s
(`many`), and a 20 s tap with compression enabled (`long`). They retain the
original and preceding-checkpoint timings for comparison. These are local CPU
measurements; speedups are not uniform across cases.
