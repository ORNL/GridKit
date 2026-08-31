# Network Current Benchmark

Two designs considered:
- **A**: Branch, LoadZ, LoadZIP, GENROU, and GENSAL have no internal network-current variables. Current is computed from bus voltage and added to the bus residuals.
- **B**: Branch has four internal algebraic terminal-current variables, while LoadZ, LoadZIP, GENROU, and GENSAL each have two internal algebraic network-current variables.

REGCA is unchanged: its two public current signals remain in both designs.

B is 1.90× slower on ACTIVSg10k, 3.88× slower on ACTIVSg2000, and 1.71× slower on WECC240. Three quarters of the added time is KLU setup and solve.

## Setup

| Case | Buses | Branches | Other models |
|---|---:|---:|---|
| ACTIVSg10k | 10,000 | 12,706 | 4,722 LoadZIP, 1,456 IEEEST, 1,456 machines, 1 BusFault |
| ACTIVSg2000 | 2,000 | 3,206 | 1,507 LoadZ, 432 GENROU, 432 IEEET1, 414 TGOV1, 2,000 BusFault |
| WECC240 | 243 | 447 | 141 LoadZIP, 103 GENROU, 37 REGCA |

10second simulation, adpative stepping, `rel_tol=1e-5`, `abs_tol=1e-7`, fault 1.0 to 1.1 s, no monitors. Clang 16 `RelWithDebInfo`, SUNDIALS 7.7.0 IDA with KLU, i9-12900K pinned to one core, five trials per cell, medians are reported here.

## Wall time

| Case | A (s) | B (s) | B/A |
|---|---:|---:|---:|
| 10k | 29.61 ± 0.94 | 56.36 ± 1.86 | 1.90× |
| 2k | 2.06 ± 0.31 | 7.98 ± 0.37 | 3.88× |
| WECC | 0.166 ± 0.006 | 0.284 ± 0.013 | 1.71× |

## Size and memory

| Case | Scenario | States | Differential | Algebraic | Jacobian nnz | nnz/state | Peak RSS (KiB) |
|---|---|---:|---:|---:|---:|---:|---:|
| 10k | A | 83,124 | 27,391 | 55,733 | 340,872 | 4.1008 | 194,212 |
| 10k | B | 146,304 | 27,391 | 118,913 | 591,680 | 4.0442 | 252,584 |
| 2k | A | 20,474 | 5,148 | 15,326 | 86,664 | 4.2329 | 56,716 |
| 2k | B | 37,176 | 5,148 | 32,028 | 156,056 | 4.1978 | 77,820 |
| WECC | A | 5,425 | 1,804 | 3,621 | 19,686 | 3.6288 | 19,580 |
| WECC | B | 7,701 | 1,804 | 5,897 | 29,146 | 3.7847 | 21,808 |

## Solver time by category (s)

| Case | Scenario | Total | Residual | Jacobian | KLU setup | KLU solve | Other |
|---|---|---:|---:|---:|---:|---:|---:|
| 10k | A | 29.30364 | 15.45922 | 3.17890 | 1.60592 | 5.08457 | 3.93142 |
| 10k | B | 56.03708 | 18.32613 | 3.45491 | 11.97431 | 15.35627 | 6.92547 |
| 2k | A | 1.99391 | 0.88292 | 0.28773 | 0.18517 | 0.36195 | 0.27048 |
| 2k | B | 7.89756 | 1.93148 | 0.31752 | 2.07130 | 2.83228 | 0.74502 |
| WECC | A | 0.15353 | 0.06115 | 0.02501 | 0.01305 | 0.02390 | 0.03090 |
| WECC | B | 0.27061 | 0.07250 | 0.02545 | 0.05673 | 0.07298 | 0.04659 |

## Added solver time, B − A

Five-run mean differences.

| Bucket | 10k (s) | 10k share | 2k (s) | 2k share | WECC (s) | WECC share |
|---|---:|---:|---:|---:|---:|---:|
| Residual | +2.914612 | 10.48% | +1.071108 | 17.83% | +0.009630 | 8.19% |
| Jacobian | +0.362189 | 1.30% | +0.019318 | 0.32% | −0.000039 | −0.03% |
| KLU setup | +10.667856 | 38.34% | +1.892193 | 31.49% | +0.043211 | 36.74% |
| KLU solve | +10.685275 | 38.41% | +2.545384 | 42.36% | +0.049212 | 41.84% |
| Other | +3.192198 | 11.47% | +0.481014 | 8.00% | +0.015606 | 13.27% |
| **Total** | **+27.822130** | **100.00%** | **+6.009017** | **100.00%** | **+0.117620** | **100.00%** |

## Per-call cost

| Case | Operation | A calls | B calls | B/A calls | A µs/call | B µs/call | B/A µs/call | B/A total time |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 10k | Residual | 2,390 | 2,186 | 0.915× | 6,468.29 | 8,383.41 | 1.296× | 1.185× |
| 10k | Jacobian | 191 | 175 | 0.916× | 16,643.45 | 19,742.34 | 1.186× | 1.087× |
| 10k | KLU setup | 191 | 175 | 0.916× | 8,407.96 | 68,424.63 | 8.138× | 7.456× |
| 10k | KLU solve | 2,390 | 2,186 | 0.915× | 2,127.43 | 7,024.83 | 3.302× | 3.020× |
| 2k | Residual | 1,110 | 1,164 | 1.049× | 795.42 | 1,659.35 | 2.086× | 2.188× |
| 2k | Jacobian | 91 | 79 | 0.868× | 3,161.84 | 4,019.23 | 1.271× | 1.104× |
| 2k | KLU setup | 91 | 79 | 0.868× | 2,034.86 | 26,218.94 | 12.885× | 11.186× |
| 2k | KLU solve | 1,110 | 1,164 | 1.049× | 326.08 | 2,433.23 | 7.462× | 7.825× |
| WECC | Residual | 799 | 839 | 1.050× | 76.54 | 86.41 | 1.129× | 1.186× |
| WECC | Jacobian | 65 | 58 | 0.892× | 384.78 | 438.83 | 1.140× | 1.018× |
| WECC | KLU setup | 65 | 58 | 0.892× | 200.79 | 978.01 | 4.871× | 4.346× |
| WECC | KLU solve | 799 | 839 | 1.050× | 29.91 | 86.98 | 2.908× | 3.054× |

## Branch residual

Share is of total model residual time.

| Case | A (s) | B (s) | B/A | A share | B share |
|---|---:|---:|---:|---:|---:|
| 10k | 6.584044 | 8.082607 | 1.228× | 49.29% | 51.00% |
| 2k | 0.278229 | 0.739956 | 2.660× | 38.05% | 45.19% |
| WECC | 0.015481 | 0.017768 | 1.148× | 29.10% | 30.07% |
