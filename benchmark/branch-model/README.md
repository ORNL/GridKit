# PhasorDynamics Branch formulation benchmark

This study compares two formulations using the same develop-based code,
instrumentation, canonical model data, solver controls, compiler settings, and
host:

1. Branch terminal currents are evaluated directly into bus current-balance
   residuals.
2. Each Branch owns four algebraic variables for its two complex terminal
   currents, and the bus equations couple to those variables.

The representative systems are the canonical ACTIVSg10k case and the Texas
ACTIVSg2000 example. The runner creates monitor-free case copies under the
requested build output directory. It never rewrites the canonical JSON files.

## Build

Use a separate build tree for each branch so both binaries remain unchanged
while trials are interleaved:

```bash
git switch lukel/branch-bench-dev
cmake --preset enzyme -B build/branch-bench-dev
cmake --build build/branch-bench-dev --target DynamicSimulation -j 10

git switch lukel/branch-bench-dev2
cmake --preset enzyme -B build/branch-bench-dev2
cmake --build build/branch-bench-dev2 --target DynamicSimulation -j 10
```

## Run

The default protocol performs one discarded warm-up per branch/case cell,
then five interleaved trials pinned to logical CPU 4:

```bash
ruby benchmark/branch-model/run_bench.rb \
  --baseline-binary build/branch-bench-dev/application/PhasorDynamics/DynamicSimulation \
  --variant-binary build/branch-bench-dev2/application/PhasorDynamics/DynamicSimulation \
  --output build/benchmark/branch-model/2026-08-27
```

The benchmark solver definitions use adaptive IDA, relative tolerance `1e-5`,
absolute tolerance `1e-7`, `max_steps=1000000`, and the develop branch's
default KLU ordering. Periodic and model-defined monitoring are disabled. The
same event schedule is generated for both formulations.

Raw stdout/stderr, GNU `time` records, generated monitor-free cases, trial JSON
and TSV, and a long-form summary TSV are written below `--output`.

## Interpreting timers

The timers are hierarchical. Do not sum nested buckets as independent shares.

- Application phases use `application_wall_seconds` as their parent.
- `consistent_ic_seconds` and `ida_solve_seconds` describe the solver envelope.
- Residual, Jacobian, KLU setup, and KLU solve timers are nested inside that
  envelope.
- Residual input/model/output timers subdivide `residual_seconds`.
- Model-family residual timers subdivide `residual_model_seconds`.
- Jacobian input/model/zero/structure/value-copy timers subdivide
  `jacobian_seconds`.
- Model-family Jacobian timers subdivide `jacobian_model_seconds`.
- System Jacobian zero/scatter timers are nested in the system/model Jacobian
  timers.

Compare both total work and cost per call. In particular, report state count,
Jacobian nonzeros, peak RSS, integration counters, KLU setup and solve time,
residual/Jacobian time, the Branch family buckets, and trial spread.
