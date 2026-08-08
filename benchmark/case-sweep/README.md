# Case sweep benchmarks

Use the canonical PhasorDynamics cases and their existing example wrappers:

- [`ACTIVSg200`](../../cases/PhasorDynamics/ACTIVSg200) with its
  [`Validation`](../../examples/PhasorDynamics/Validation/ACTIVSg200) wrapper
- [`ACTIVSg500`](../../cases/PhasorDynamics/ACTIVSg500) with its
  [`Validation`](../../examples/PhasorDynamics/Validation/ACTIVSg500) wrapper
- [`WECC240`](../../cases/PhasorDynamics/WECC240) with its
  [`Validation`](../../examples/PhasorDynamics/Validation/WECC240) wrapper
- [`ACTIVSg10k`](../../cases/PhasorDynamics/ACTIVSg10k) with its
  [`Large`](../../examples/PhasorDynamics/Large/ACTIVSg10k) wrapper

Use identical, explicit solver settings for each comparison. Record the
`DynamicSimulation` application's `Complete in` time and report the median of
three runs. ACTIVSg200, ACTIVSg500, and WECC240 provide correctness coverage;
ACTIVSg10k provides the large scaling point.

The [`2026-08-26-current-cases`](2026-08-26-current-cases) study contains small
solver definitions that reference the canonical cases directly. It uses a
10-second horizon, a 100 ms fault, relative tolerance `1e-5`, absolute
tolerance `1e-7`, AMD ordering, and no periodic monitor output. Run it with:

```bash
ruby benchmark/case-sweep/2026-08-26-current-cases/run_bench.rb
```

Do not copy case JSON files into this directory. Generated monitor CSV files,
standard output, and plots remain untracked.
