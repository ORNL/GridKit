# Signal delay results

Open [index.html](index.html) for the self-contained gallery, or [gallery.pdf](gallery.pdf)
for eight pages of plots. PNG and SVG versions are in `figures/`.

The gallery contains 44 runs against the local GridKit implementation:

1. Actual `DynamicSimulation` output from a two-bus GENROU case. Smooth mechanical
   power forcing, `pmech = cos(t)`, excites a machine-speed signal with three delay taps.
2. Actual `DynamicSimulation` output from the JSON cosine-source example, with
   three delays and constant prehistory.
3. Twenty-eight fixed-step runs: four delays and seven step sizes.
4. Six feedback runs with an independent exact exponential solution.
5. Three smooth-pulse runs with different delays.
6. Two runs comparing history compression on and off.
7. Three smooth-chirp runs with different error tolerances.

The 42 instrumented runs use the actual GridKit `SystemModel`, `SignalDelay`,
`SignalHistory`, and IDA wrapper, linked from `build/`. The standalone
[gallery.cpp](gallery.cpp) adds only manufactured source/feedback components and
records solver diagnostics. It does not implement another delay algorithm.

The complete instrumented CSV files and logs are grouped in `data.zip` to keep
the directory compact. The two application CSVs are compressed separately.

`h` is an actual accepted internal solver step. `dt_monitor` only determines the
uniform output sampling. Data includes every accepted step midpoint, so fixed
steps that align with monitor times do not hide interpolation errors.

## Reading the results

- `summary.csv` and `summary.json`: measured errors, steps, rejections, and segment counts.
- `runs.json`: all 42 harness experiment settings. Absolute tolerance is relative
  tolerance times 0.01, for both the solver and the history controller.
- `data.zip: *_samples.csv`: uniform callback output plus the initial condition.
- `data.zip: *_midpoints.csv`: IDA dense output at every accepted step midpoint.
- `data.zip: *_steps.csv`: accepted endpoint values, step sizes, history counts,
  and the history error estimate.
- `data.zip: *_attempts.csv`: candidate steps and additional history-control rejections.
- `app/` and `grid/`: JSON inputs, application logs, and compressed original monitor CSVs.
- `provenance.json`: source revision, local diff fingerprint, build description, and caveats.

Maximum errors are sampled maxima, not continuous error bounds. Analytic
references are used for all manufactured problems. The generator plot instead
compares with its own recorded source shifted in time; it is a consistency check.

The fixed-step sweep deliberately uses loose history tolerances (`rtol=1`,
`atol=0.01`) so that selected coarse steps are measured rather than rejected.
GridKit's existing fixed-step mode disables the IDA integration error test and
caps BDF order at two. Adaptive experiments use matched solver/history tolerances.
All harness cases use analytic prehistory; the JSON application cases use default
constant prehistory.

Compression-off still prunes history outside the largest required delay window.
Segment counts describe retained cubic segments, not total process memory.
No runtime speedup is claimed. Rejection counts exclude IDA's own internal
error-test and nonlinear-solve failures.

## Reproduce

From the repository root, with the existing CMake Makefiles build and NumPy/Matplotlib:

```bash
cmake --build build --target DynamicSimulation test_pd_signal_delay -j 10
python3 examples/PhasorDynamics/DynamicSimulation/Toy/DelayedSignal/gallery/generate.py
```

Compilation and raw intermediate output use a temporary directory. Compressed
results and figures are written here. To redraw the saved results without running
GridKit again:

```bash
python3 examples/PhasorDynamics/DynamicSimulation/Toy/DelayedSignal/gallery/plot_gallery.py
```
