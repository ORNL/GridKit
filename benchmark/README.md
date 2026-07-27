# GridKit benchmarks

This directory contains portable benchmark methodology and drivers. Benchmark
inputs come from the canonical cases in [`../cases/PhasorDynamics`](../cases/PhasorDynamics);
do not copy case or solver JSON files into this directory.

Generated standard output, `/usr/bin/time` records, monitor CSV files, profiles,
and plots belong under `build/` or `/tmp`, not in source control.

Current benchmark families:

- [`adaptive-step/`](adaptive-step/): PhasorDynamics adaptive-step simulation.
