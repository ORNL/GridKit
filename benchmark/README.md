# GridKit benchmarks

This directory contains reproducible benchmark drivers and methodology. The
drivers use canonical case inputs without modifying them; generated inputs,
raw output, timing records, profiles, and plots belong under `build/` or
`/tmp`, not in source control.

Current benchmark families:

- [`branch-model/`](branch-model/): compares the PhasorDynamics Branch model
  with direct bus-current coupling against four solver-owned terminal-current
  algebraic variables.
