# Signal delay gallery

Open [index.html](index.html). Fourteen sections contain 36 PNG views in
`figures/`. Choose a variant beneath each heading to see one experiment at a
time. Every variant also opens as a standalone PNG, including without JavaScript.
The editable source is `plot_gallery.py`; recorded CSV data are in `data.zip`.
Only PNGs are generated.

| Figure | Coverage |
| --- | --- |
| 01 | Smooth rising transition with several delays |
| 02 | Gaussian pulse widths |
| 03 | Broad and narrow features in one source |
| 04 | Oscillatory packet |
| 05 | Adaptive pulse propagation |
| 06 | Long retained history |
| 07 | Compression on and off |
| 08 | Sharp rising/falling edges and an asymmetric pulse |
| 09 | Bipolar sine waves at 1, 5, and 20 Hz |
| 10 | Increasing-frequency chirp |
| 11 | Continuous ramp and triangle signals with derivative corners |
| 12 | Narrow Gaussian pulses at different offsets |
| 13 | An intentionally missed pulse beside the source-resolved result |
| 14 | Accepted source knots and reconstructed history near an edge/corner |

The driver uses the actual GridKit history, delay component, and IDAS wrapper.
Analytic references evaluate the source at `t - tau`; prehistory provides both
value and derivative. Signal values are continuous. No value jumps are tested.

Ordinary views contain two aligned panels: the input and delayed output above,
and the accepted step-to-delay ratio `h/tau` below. The input is always dashed
grey and the delayed output solid blue. A dashed line at `h/tau = 1` marks a step
equal to the delay. Each view uses the same type, line weights, axes, and palette.
There are no error panels, retry markers, or duplicate absolute-step panels.

Sine views show at most three cycles. Short time windows use milliseconds with
an explicit origin; other views use seconds. The plots retain the recorded
samples and step sequence without smoothing or averaging.

Compression has a dedicated comparison of retained segments and `h/tau`.
The source-resolution figure compares missed and resolved output against sparse
analytic reference markers. The history detail shows actual accepted values
and their interpolation polynomials near a smooth edge and a derivative corner.
Error measurements and retry counts remain in the data and the gallery's
collapsed numerical-results section.

`h` is the accepted internal step, independent of the monitor interval. `h_min`
is a failure floor. `h_max` limits source sampling gaps and must be chosen for
the feature being resolved. The uncapped 1 ms Gaussian intentionally demonstrates
that an unobserved source pulse can be missed despite passing history checks.
It is explicitly marked `expected_missed`; every other run must pass its sampled
analytic-error check. Sparse probing cannot guarantee resolution of arbitrary
unknown features.

Smooth edge `width` is the 10–90% transition time for the new rise/fall cases.
The original `ramp` uses a tanh scale parameter. Gaussian `width` is the scale in
`exp(-((t-center)/width)^2)`, not FWHM. Frequencies are in Hz. The corner cases
compare adaptive refinement with explicit source-stencil restarts at known
corner times. Compression preserves marked corners.

## Data

- `runs.json`: complete experiment settings, including per-run step limits.
- `summary.csv` and `summary.json`: sampled errors, steps, retries, and storage.
- `data.zip`: monitor samples, accepted-step midpoints, steps, history-check
  attempts, accepted source-polynomial coefficients, and logs.
- `provenance.json`: the original simulation source/build record, plus a separate
  rendering record with plot-source and data hashes, Matplotlib version, and PNG list.
- `performance/`: separate `DynamicSimulation` timings using identical inputs
  and the median of three app-reported CPU times.

Errors are sampled maxima. Segment counts exclude container overhead and other
process memory. Compression bounds concern stored history; they do not bound
the difference between complete numerical trajectories. Compression may change
the adaptive step sequence.

## Reproduce

From the repository root, with the configured build, the
[IDAS step-check extension](../../../../../../buildsystem/spack_repo/gridkit/packages/sundials/README.md),
and NumPy/Matplotlib:

```bash
cmake -S . -B build
python3 -B examples/PhasorDynamics/DynamicSimulation/Toy/DelayedSignal/gallery/generate.py
```

Generation builds the CMake `signal_delay_gallery` target, runs all experiments,
validates the data, and redraws the figures. Use `--build-dir` for another
configured build. Raw intermediates use a temporary directory.

To redraw all PNGs from saved data:

```bash
python3 -B examples/PhasorDynamics/DynamicSimulation/Toy/DelayedSignal/gallery/plot_gallery.py
```

Application timings can be refreshed separately:

```bash
cmake --build build --target DynamicSimulation -j 10
python3 -B examples/PhasorDynamics/DynamicSimulation/Toy/DelayedSignal/gallery/performance/measure.py
```
