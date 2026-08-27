# Branch Current Benchmark

Two designs considered:
- **A**: Branch has no internal variables. Terminal current is computed from bus voltages and added to the bus residuals.
- **B**: Branch has four internal algebraic terminal-current variables and four equations.

B is 1.89× slower on ACTIVSg10k and 3.91× slower on ACTIVSg2000. Three quarters of the added time is KLU setup and solve.

## Setup

| Case | Buses | Branches | Other models |
|---|---:|---:|---|
| ACTIVSg10k | 10,000 | 12,706 | 4,722 LoadZIP, 1,456 IEEEST, 1,456 machines, 1 BusFault |
| ACTIVSg2000 | 2,000 | 3,206 | 1,507 LoadZ, 432 GENROU, 432 IEEET1, 414 TGOV1, 2,000 BusFault |

10second simulation, adpative stepping, `rel_tol=1e-5`, `abs_tol=1e-7`, fault 1.0 to 1.1 s, no monitors. Clang 16 `RelWithDebInfo`, SUNDIALS 7.7.0 IDA with KLU, i9-12900K pinned to one core, five trials per cell, medians are reported here.

## Wall time

| Case | A (s) | B (s) | B/A |
|---|---:|---:|---:|
| 10k | 27.30 ± 0.13 | 51.63 ± 0.08 | 1.89× |
| 2k | 1.71 ± 0.08 | 6.70 ± 0.18 | 3.91× |

## Size and memory

| Case | Scenario | States | Differential | Algebraic | Jacobian nnz | nnz/state | Peak RSS (KiB) |
|---|---|---:|---:|---:|---:|---:|---:|
| 10k | A | 95,480 | 27,391 | 68,089 | 384,472 | 4.0267 | 204,272 |
| 10k | B | 146,304 | 27,391 | 118,913 | 591,680 | 4.0442 | 252,488 |
| 2k | A | 24,352 | 5,148 | 19,204 | 100,448 | 4.1248 | 59,560 |
| 2k | B | 37,176 | 5,148 | 32,028 | 156,056 | 4.1978 | 77,840 |

## Solver time by category (s)

| Case | Scenario | Total | Residual | Jacobian | KLU setup | KLU solve | Other |
|---|---|---:|---:|---:|---:|---:|---:|
| 10k | A | 26.98685 | 13.94361 | 2.53241 | 1.38485 | 4.88375 | 4.14112 |
| 10k | B | 51.29529 | 16.40101 | 3.21097 | 11.49818 | 13.93948 | 6.17889 |
| 2k | A | 1.64811 | 0.65635 | 0.26427 | 0.19110 | 0.28221 | 0.25418 |
| 2k | B | 6.62895 | 1.49944 | 0.27648 | 1.97699 | 2.23319 | 0.62316 |

## Added solver time, B − A

Five-run mean differences.

| Bucket | 10k (s) | 10k share | 2k (s) | 2k share |
|---|---:|---:|---:|---:|
| Residual | +2.772767 | 11.13% | +0.881901 | 17.76% |
| Jacobian | +0.724118 | 2.91% | +0.013749 | 0.28% |
| KLU setup | +10.147749 | 40.74% | +1.782893 | 35.91% |
| KLU solve | +9.124068 | 36.63% | +1.922413 | 38.72% |
| Other | +2.137470 | 8.58% | +0.364518 | 7.34% |
| **Total** | **+24.906171** | **100.00%** | **+4.965474** | **100.00%** |

## Per-call cost

| Case | Operation | A calls | B calls | B/A calls | A µs/call | B µs/call | B/A µs/call | B/A total time |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 10k | Residual | 2,320 | 2,152 | 0.928× | 6,010.18 | 7,621.29 | 1.268× | 1.176× |
| 10k | Jacobian | 174 | 189 | 1.086× | 14,554.09 | 16,989.26 | 1.167× | 1.268× |
| 10k | KLU setup | 174 | 189 | 1.086× | 7,958.92 | 60,836.91 | 7.644× | 8.303× |
| 10k | KLU solve | 2,320 | 2,152 | 0.928× | 2,105.07 | 6,477.45 | 3.077× | 2.854× |
| 2k | Residual | 1,078 | 1,164 | 1.080× | 608.86 | 1,288.18 | 2.116× | 2.285× |
| 2k | Jacobian | 96 | 79 | 0.823× | 2,752.86 | 3,499.73 | 1.271× | 1.046× |
| 2k | KLU setup | 96 | 79 | 0.823× | 1,990.60 | 25,025.20 | 12.572× | 10.345× |
| 2k | KLU solve | 1,078 | 1,164 | 1.080× | 261.79 | 1,918.54 | 7.329× | 7.913× |

## Branch residual

Share is of total model residual time.

| Case | A (s) | B (s) | B/A | A share | B share |
|---|---:|---:|---:|---:|---:|
| 10k | 5.654278 | 7.283599 | 1.288× | 47.09% | 51.48% |
| 2k | 0.195917 | 0.587551 | 2.999× | 35.22% | 45.72% |
