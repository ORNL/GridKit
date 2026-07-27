# Adaptive-step simulation benchmarks

Adaptive-step benchmarks run `DynamicSimulation` with the repository's canonical
PhasorDynamics cases and their existing solver wrappers:

| Case data | Solver wrapper |
|---|---|
| `cases/PhasorDynamics/ACTIVSg200/ACTIVSg200.case.json` | `examples/PhasorDynamics/Validation/ACTIVSg200/ACTIVSg200.solver.json` |
| `cases/PhasorDynamics/ACTIVSg500/ACTIVSg500.case.json` | `examples/PhasorDynamics/Validation/ACTIVSg500/ACTIVSg500.solver.json` |
| `cases/PhasorDynamics/WECC240/WECC240.case.json` | `examples/PhasorDynamics/Validation/WECC240/WECC240.solver.json` |
| `cases/PhasorDynamics/ACTIVSg10k/ACTIVSg10k.case.json` | `examples/PhasorDynamics/Large/ACTIVSg10k/ACTIVSg10k.solver.json` |

CMake stages these canonical inputs in the corresponding build-tree example
directories. This generated staging is not a second tracked copy of a case.

Build the application before collecting timings:

```bash
cmake --build build --target DynamicSimulation -j 10
```

Use the application-reported `Complete in` time and report the median of three
runs made with identical solver settings. Put generated solver variations and
all resulting output under `build/` or `/tmp`.
