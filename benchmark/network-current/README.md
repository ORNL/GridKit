# PhasorDynamics network-current formulation benchmark

This study compares two aggregate formulations using the same model data,
solver controls, instrumentation, compiler settings, and host:

1. `direct_injection`: Branch, LoadZ, LoadZIP, GENROU, and GENSAL evaluate
   their network-current contributions directly into bus current-balance
   residuals.
2. `internal_current_variables`: those five model families own algebraic
   network-current variables coupled to the bus equations.

REGCA retains its public `IR` and `II` signals in both formulations. Its 37
instances, or 74 current states, in WECC240 are therefore not part of the state
delta. REPCA consumes those currents as public signals, and the current
single-column `SignalNode` contract cannot derive them from the two bus-voltage
columns without an interface and Jacobian redesign. GenClassical is also
unchanged, but has its own timing and count bucket.

The runner copies the canonical cases without monitoring and writes generated
solver files below the requested build output. Solver files contain absolute
paths to those generated cases, and each simulation runs from its generated
input directory so relative study paths retain their normal interpretation.
Canonical JSON files are never rewritten.

| Case | Branch | LoadZ | LoadZIP | GENROU | GENSAL | Expected state increase |
|---|---:|---:|---:|---:|---:|---:|
| ACTIVSg10k | 12,706 | 0 | 4,722 | 926 | 530 | 63,180 |
| ACTIVSg2000 | 3,206 | 1,507 | 0 | 432 | 0 | 16,702 |
| WECC240 | 447 | 0 | 141 | 103 | 0 | 2,276 |

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

## Run the timing benchmark

The default protocol performs one discarded warm-up per formulation/case cell,
then five interleaved trials pinned to logical CPU 4:

```bash
ruby benchmark/network-current/run_bench.rb \
  --baseline-binary build/branch-bench-dev/application/PhasorDynamics/DynamicSimulation \
  --variant-binary build/branch-bench-dev2/application/PhasorDynamics/DynamicSimulation \
  --baseline-ref lukel/branch-bench-dev \
  --variant-ref lukel/branch-bench-dev2 \
  --output build/benchmark/network-current/2026-08-31
```

The generated studies use adaptive IDA, relative tolerance `1e-5`, absolute
tolerance `1e-7`, `max_steps=1000000`, no periodic or model-defined monitoring,
and the same 10-second fault schedule. Every measured pair must match the state
increase in the table; model counts emitted by each binary must match the
canonical inputs.

Raw stdout/stderr, GNU `time` records, generated inputs, and the following
machine-readable reports are written below `--output`:

- `trials.json` records binary/input SHA-256 hashes and nested metrics.
- `trials.tsv` contains one wide record per measured run.
- `summary.tsv` reports median, MAD, minimum, and maximum for every metric.

The reports include state and Jacobian sizes, peak RSS, IDA counters, residual
and Jacobian time, KLU setup/solve totals and per-call costs, and separate
REGCA, Branch, LoadZ, LoadZIP, GENROU, GENSAL, and GenClassical timers. Timers
are hierarchical, so nested buckets must not be summed as independent shares.
Sparse-factor fill is not measured and should be described only as an
inference from the KLU results.

## Produce paired FlameGraphs

FlameGraph is an external prerequisite and is not vendored. Use a checkout at
the pinned upstream commit:

```bash
git clone https://github.com/brendangregg/FlameGraph.git /path/to/FlameGraph
git -C /path/to/FlameGraph checkout 41fee1f99f9276008b7cd112fca19dc3ea84ac32
```

The host must also provide a `perf` executable compatible with its running
Linux kernel. On WSL2 this commonly requires matching WSL kernel tools. The
driver does not install packages or change `perf_event_paranoid`; it fails
before warm-ups or captures if `perf` is absent, and runs a sampling smoke test
before starting the profile matrix.

```bash
ruby benchmark/network-current/run_flamegraphs.rb \
  --trials build/benchmark/network-current/2026-08-31/trials.json \
  --flamegraph-dir /path/to/FlameGraph
```

The driver verifies the recorded binary, generated-case, and solver hashes,
then uses those exact binaries and inputs. It performs one unprofiled warm-up
and three interleaved captures per formulation/case with:

```text
perf record -F 99 -e cpu-clock:u --call-graph dwarf,16384
```

Each capture is folded with `stackcollapse-perf.pl --all`. The three captures
for a cell are summed before rendering a stable `--hash` FlameGraph. For each
case, `difffolded.pl -n` also produces a normalized differential whose widths
show the internal-current-variable formulation. The driver rejects empty
samples, unresolved GridKit/IDA/KLU profiles, noninteractive SVGs, and an
above-5-percent share of samples containing any unknown frame by default.

All `perf.data`, unfolded and folded stacks, logs, metadata, and SVGs remain
under `build/benchmark/network-current/`, which is ignored by Git. Do not
commit generated profiles or plots.
