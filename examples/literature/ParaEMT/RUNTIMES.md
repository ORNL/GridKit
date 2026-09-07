# ParaEMT / GridKit runtimes

Three-second trajectories; medians of three processes per setting. Benchmark processes run sequentially on one pinned logical CPU with BLAS/OMP/Numba thread counts set to one. Result capture is disabled. Every benchmark final state agrees with its corresponding captured trajectory to `rtol=atol=1e-12`. ParaEMT stepping kernels are warmed on a disposable state before timing; their JIT signatures must remain unchanged throughout the timed loop. GridKit uses the previously measured trials of its unchanged compiled executable.

Host: **12th Gen Intel(R) Core(TM) i9-12900K**; pinned logical CPU 0. GridKit revision `257ccb317a7b5d4a530ba6b96688c6ff3c6f5c49`.

| Event | Simulator / setting | Loop wall, s | Loop CPU, s | Initialization, s | Process wall, s | Steps | Duration / steps, µs |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| governor_step | ParaEMT dt50us | 7.0475 | 6.3753 | 17.6761 | 25.6580 | 60,000 | 50.000 |
| governor_step | ParaEMT dt25us | 12.3489 | 11.6747 | 18.1733 | 31.1864 | 120,000 | 25.000 |
| governor_step | ParaEMT dt12_5us | 29.3164 | 27.0875 | 17.9584 | 48.7472 | 240,000 | 12.500 |
| governor_step | GridKit tol1e-7 | 1.0623 | 1.0192 | 0.0197 | 1.1034 | 72,480 | 41.391 |
| governor_step | GridKit tol1e-8 | 1.4668 | 1.3794 | 0.0222 | 1.5097 | 98,255 | 30.533 |
| governor_step | GridKit tol1e-9 | 2.1716 | 2.0417 | 0.0145 | 2.2076 | 133,682 | 22.441 |
| trip | ParaEMT dt50us | 5.8573 | 5.6747 | 15.3386 | 21.9829 | 60,000 | 50.000 |
| trip | ParaEMT dt25us | 11.2496 | 10.9187 | 15.5911 | 27.7022 | 120,000 | 25.000 |
| trip | ParaEMT dt12_5us | 25.7674 | 24.3915 | 15.3789 | 41.9238 | 240,000 | 12.500 |
| trip | GridKit tol1e-7 | 1.2474 | 1.2065 | 0.0158 | 1.2871 | 95,002 | 31.578 |
| trip | GridKit tol1e-8 | 1.8089 | 1.7352 | 0.0168 | 1.8449 | 137,888 | 21.757 |
| trip | GridKit tol1e-9 | 2.6526 | 2.5607 | 0.0145 | 2.6928 | 194,152 | 15.452 |

ParaEMT initialization here includes explicit stepping-kernel warm-up in addition to network initialization. Its ordinary startup message does not mean all kernels are compiled: the first-use costs inside the original loops are retained in [RUNTIMES_COLD.md](RUNTIMES_COLD.md). Process wall time includes imports, loading, initialization, warm-up, the timed trajectory and metadata output. GridKit loop timing includes the event consistency solve. CPU time is measured inside each process.

**Adaptive steps and output samples are different.** GridKit uses variable-step, variable-order IDA/BDF; the monitor interval is 50 µs. ParaEMT uses fixed 50, 25 or 12.5 µs integration steps. Both captured outputs use 50 µs spacing, giving 60,001 samples. The duration/steps column is an effective average, not a fixed GridKit step or a recorded step-size distribution. IDA may internally pass an output time and interpolate back.

This is a cost comparison at the listed settings, not an equal-error benchmark. Use the waveform/refinement metrics alongside runtime. The governor event has closely mapped continuous models; the trip has different G1 post-event semantics, so trip speed ratios do not establish equivalent-model performance.

## Result-capture runs

These single runs generate the plotted data. GridKit streams monitor and full-state CSVs during its loop; ParaEMT copies samples into memory and exports files after its loop. Their capture-loop costs therefore include different output work and are not the primary timing comparison. They were not CPU-pinned. Compression and plotting are excluded from the benchmark loops.

| Event | Simulator / setting | Capture loop wall, s | No-capture median loop wall, s |
| --- | --- | ---: | ---: |
| governor_step | ParaEMT dt50us | 16.4186 | 7.0475 |
| governor_step | ParaEMT dt25us | 22.1297 | 12.3489 |
| governor_step | ParaEMT dt12_5us | 35.4347 | 29.3164 |
| governor_step | GridKit tol1e-7 | 10.7793 | 1.0623 |
| governor_step | GridKit tol1e-8 | 11.0520 | 1.4668 |
| governor_step | GridKit tol1e-9 | 11.7881 | 2.1716 |
| trip | ParaEMT dt50us | 15.7052 | 5.8573 |
| trip | ParaEMT dt25us | 21.3726 | 11.2496 |
| trip | ParaEMT dt12_5us | 33.5419 | 25.7674 |
| trip | GridKit tol1e-7 | 10.9549 | 1.2474 |
| trip | GridKit tol1e-8 | 12.1528 | 1.8089 |
| trip | GridKit tol1e-9 | 12.6822 | 2.6526 |

[Raw trials, environment and ranges](results/runtime/) · [Runtime plot](plots/runtime_comparison.png) · [Governor adaptive work](plots/governor_step_adaptive_work.png) · [Trip adaptive work](plots/trip_adaptive_work.png)

Reproduce with the configured Python environment: `python benchmark.py`. Use `python benchmark.py --summarize-only` to rebuild tables and plots from saved trials. Each raw trial records its command, executable/wrapper hash, settings, timing, final state and solver work.
