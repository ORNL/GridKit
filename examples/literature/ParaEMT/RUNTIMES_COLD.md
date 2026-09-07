# Cold-process runtimes

Three-second trajectories; medians of three fresh processes per setting. All benchmark processes run sequentially on one pinned logical CPU with BLAS/OMP/Numba thread counts set to one. Result capture is disabled. Every benchmark final state agrees with its corresponding captured trajectory to `rtol=atol=1e-12`.

Host: **12th Gen Intel(R) Core(TM) i9-12900K**; pinned logical CPU 0. GridKit revision `257ccb317a7b5d4a530ba6b96688c6ff3c6f5c49`.

| Event | Simulator / setting | Loop wall, s | Loop CPU, s | Initialization, s | Process wall, s | Steps | Duration / steps, µs |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| governor_step | ParaEMT dt50us | 20.3402 | 18.4711 | 6.1526 | 27.3904 | 60,000 | 50.000 |
| governor_step | ParaEMT dt25us | 25.2189 | 23.4087 | 6.7371 | 33.0259 | 120,000 | 25.000 |
| governor_step | ParaEMT dt12_5us | 34.8503 | 33.3000 | 5.9369 | 42.0295 | 240,000 | 12.500 |
| governor_step | GridKit tol1e-7 | 1.0623 | 1.0192 | 0.0197 | 1.1034 | 72,480 | 41.391 |
| governor_step | GridKit tol1e-8 | 1.4668 | 1.3794 | 0.0222 | 1.5097 | 98,255 | 30.533 |
| governor_step | GridKit tol1e-9 | 2.1716 | 2.0417 | 0.0145 | 2.2076 | 133,682 | 22.441 |
| trip | ParaEMT dt50us | 17.8341 | 16.8877 | 5.5528 | 24.0667 | 60,000 | 50.000 |
| trip | ParaEMT dt25us | 24.5051 | 23.3151 | 5.4582 | 32.1707 | 120,000 | 25.000 |
| trip | ParaEMT dt12_5us | 37.2346 | 35.1962 | 6.1495 | 44.2438 | 240,000 | 12.500 |
| trip | GridKit tol1e-7 | 1.2474 | 1.2065 | 0.0158 | 1.2871 | 95,002 | 31.578 |
| trip | GridKit tol1e-8 | 1.8089 | 1.7352 | 0.0168 | 1.8449 | 137,888 | 21.757 |
| trip | GridKit tol1e-9 | 2.6526 | 2.5607 | 0.0145 | 2.6928 | 194,152 | 15.452 |

**The ParaEMT loop includes lazy compilation of stepping kernels on first use.** Initialization compiles the network kernel but does not warm every stepping kernel. These are cold-process costs, not warmed integration throughput. Process wall time additionally includes Python imports, library loading and result-metadata output. GridKit loop timing includes the event consistency solve. CPU time is measured inside each process.

**Adaptive steps and output samples are different.** GridKit uses variable-step, variable-order IDA/BDF; the monitor interval is 50 µs. ParaEMT uses fixed 50, 25 or 12.5 µs integration steps. Both captured outputs use 50 µs spacing, giving 60,001 samples. The duration/steps column is an effective average, not a fixed GridKit step or a recorded step-size distribution. IDA may internally pass an output time and interpolate back.

This is a cost comparison at the listed settings, not an equal-error benchmark. Use the waveform/refinement metrics alongside runtime. The governor event has closely mapped continuous models; the trip has different G1 post-event semantics, so trip speed ratios do not establish equivalent-model performance.

## Result-capture runs

These single runs generate the plotted data. GridKit streams monitor and full-state CSVs during its loop; ParaEMT copies samples into memory and exports files after its loop. Their capture-loop costs therefore include different output work and are not the primary timing comparison. They were not CPU-pinned. Compression and plotting are excluded from the benchmark loops.

| Event | Simulator / setting | Capture loop wall, s | No-capture median loop wall, s |
| --- | --- | ---: | ---: |
| governor_step | ParaEMT dt50us | 16.4186 | 20.3402 |
| governor_step | ParaEMT dt25us | 22.1297 | 25.2189 |
| governor_step | ParaEMT dt12_5us | 35.4347 | 34.8503 |
| governor_step | GridKit tol1e-7 | 10.7793 | 1.0623 |
| governor_step | GridKit tol1e-8 | 11.0520 | 1.4668 |
| governor_step | GridKit tol1e-9 | 11.7881 | 2.1716 |
| trip | ParaEMT dt50us | 15.7052 | 17.8341 |
| trip | ParaEMT dt25us | 21.3726 | 24.5051 |
| trip | ParaEMT dt12_5us | 33.5419 | 37.2346 |
| trip | GridKit tol1e-7 | 10.9549 | 1.2474 |
| trip | GridKit tol1e-8 | 12.1528 | 1.8089 |
| trip | GridKit tol1e-9 | 12.6822 | 2.6526 |

[Raw trials, environment and ranges](results/runtime_cold/) · [Runtime plot](plots/runtime_cold/runtime_comparison.png) · [Governor adaptive work](plots/runtime_cold/governor_step_adaptive_work.png) · [Trip adaptive work](plots/runtime_cold/trip_adaptive_work.png)

To repeat a cold ParaEMT trial using the current wrapper, append `--cold-start` to its saved `--benchmark` command. The original commands refer to the earlier wrapper hash recorded in each trial. GridKit's `--benchmark` command is unchanged. For example:

```bash
python run_reference.py --event trip --dt-us 12.5 --benchmark --cold-start --output /tmp/paraemt-cold-trip
```
