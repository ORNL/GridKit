# Algebraic Collapse Benchmark

`develop` (`f559083d8`) vs `lukel/algebraic-collapse-dev` (`379909220`).

## Setup

- 10 s simulation, bus fault 1.0 to 1.1 s, adaptive IDA, `rel_tol=1e-5`,
  `abs_tol=1e-7`, `max_steps=1e6`, no monitors (`dt_monitor=0`, one `IDASolve`
  per event segment).
- Clang 16 `RelWithDebInfo`, SUNDIALS 7.8.0 IDA with KLU, i9-12900K under WSL2,
  pinned to one logical CPU.
- One discarded warm-up, then 5 interleaved trials per cell. Time is the
  app-reported `Complete in`; peak RSS is from GNU `time`.
- Both builds run from binary and library snapshots, so neither is rebuilt
  between trials.

## Results

| Case | Unknowns | develop (s) | collapse (s) | Speedup | Paired ratio | Peak RSS (MiB) |
|---|---|---:|---:|---:|---|---|
| Hawaii[^hawaii] | 1,654 → 854 | 0.103 ± 0.005 | 0.055 ± 0.002 | 1.86× | 0.56 [0.49, 0.62] | 12.1 → 11.2 |
| IEEE39 | 596 → 366 | 0.016 ± 0.001 | 0.015 ± 0.001 | 1.09× | 0.89 [0.66, 1.10] | 10.2 → 10.0 |
| ACTIVSg200 | 1,540 → 882 | 0.056 ± 0.005 | 0.047 ± 0.006 | 1.21× | 0.83 [0.81, 0.95] | 12.3 → 11.7 |
| ACTIVSg500 | 2,806 → 1,722 | 0.116 ± 0.006 | 0.086 ± 0.006 | 1.35× | 0.74 [0.65, 0.81] | 15.0 → 13.9 |
| WECC240 | 5,670 → 2,839 | 0.183 ± 0.025 | 0.135 ± 0.001 | 1.35× | 0.70 [0.63, 0.93] | 19.1 → 15.9 |
| ACTIVSg2000 | 20,314 → 11,746 | 3.79 ± 0.52 | 2.31 ± 0.07 | 1.64× | 0.61 [0.40, 0.74] | 50.4 → 40.5 |
| ACTIVSg10k | 92,836 → 55,203 | 31.9 ± 3.0 | 26.6 ± 1.0 | 1.20× | 0.84 [0.72, 0.88] | 190.9 → 144.9 |

- Times are median ± median absolute deviation; speedup is the ratio of medians.
- Paired ratio is collapse/develop within each interleaved round, given as the
  median [min, max]. It cancels round-to-round host noise.
- Host noise on this run is higher than in earlier studies; trust the paired
  ratios over single medians for IEEE39, ACTIVSg2000, and ACTIVSg10k.

[^hawaii]: Both builds use `abs_tol=1e-5`, the Hawaii validation tolerance.
    At `abs_tol=1e-7`, develop fails `IDACalcIC` (`IDA_LINESEARCH_FAIL`) and the
    collapsed branch completes.
