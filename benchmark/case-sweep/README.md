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

Generate temporary solver variants, standard output, timing records, monitor
CSV files, and summary tables under `build/` or `/tmp`. Do not copy case or
solver JSON files into this directory.
