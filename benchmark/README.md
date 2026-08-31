# GridKit benchmarks

This directory contains reproducible benchmark drivers and methodology. The
drivers use canonical case inputs without modifying them; generated inputs,
raw output, timing records, profiles, and plots belong under `build/` or
`/tmp`, not in source control.

Current benchmark families:

- [`network-current/`](network-current/): compares direct PhasorDynamics
  network-current injection against solver-owned current variables and provides
  paired timing and FlameGraph drivers.
