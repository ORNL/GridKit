# Network-current formulation benchmark

The aggregate `internal_current_variables` formulation is slower than
`direct_injection` in all three cases: 1.910x on ACTIVSg10k, 3.951x on
ACTIVSg2000, and 1.782x on WECC240. All 30 measured runs succeeded, and the
required state increase is exact in all 15 measured pairs.

The largest measured change is in KLU. Using differences between median timing
buckets, added KLU setup and solve time represents 77.1%, 73.4%, and 77.4% of
the added application-reported CPU time, respectively. KLU setup cost per call
grows by 8.14x, 12.89x, and 4.87x; solve cost per call grows by 3.30x, 7.46x,
and 2.91x.

## Formulations

- `direct_injection` at `f503b3057bc8caa16fde0a12434c83ba0f1dcb3a`:
  Branch, LoadZ, LoadZIP, GENROU, and GENSAL contribute network current
  directly to the bus current balance.
- `internal_current_variables` at
  `3a5953959857a49817fec8b850a62ef63d6bea2c`: Branch owns four terminal-current
  variables, and LoadZ, LoadZIP, GENROU, and GENSAL each own two network-current
  variables.

REGCA is deliberately identical in both endpoints. Its public `IR` and `II`
signals are consumed by REPCA and cannot currently be derived through the
single-column `SignalNode` contract from both bus-voltage columns. WECC240 has
37 REGCA instances, so its 74 REGCA current states are present in both
formulations and are excluded from the state delta. GenClassical is also
unchanged and has no instances in these three cases.

## Protocol

- Separate, complete `RelWithDebInfo` build trees; Clang 16.0.6 with
  `-march=haswell -mtune=native -O2 -g -DNDEBUG`.
- SUNDIALS 7.7.0 IDA with KLU.
- Intel Core i9-12900K under WSL2 kernel
  `6.18.33.2-microsoft-standard-WSL2`, pinned to logical CPU 4.
- One discarded warm-up for each formulation/case cell, then five
  rotated/interleaved measured trials.
- Adaptive IDA, `rel_tol=1e-5`, `abs_tol=1e-7`,
  `max_steps=1000000`, no monitors, and a 10-second simulation with the fault
  applied from 1.0 to 1.1 seconds.
- Headline time is the application's `Complete in` CPU time. MAD is the
  unscaled median absolute deviation. Peak RSS is from GNU `time`.

The frozen binaries and SHA-256 hashes were:

| Formulation | Binary SHA-256 |
|---|---|
| `direct_injection` | `a3595023354791903217dae7037021925218680550f7fcc6df960c20a2e05f46` |
| `internal_current_variables` | `a2ef4db194251b5a2fc30b15b69c74140375f9f2c52077c77b67ecbd79ef3715` |

The ignored raw records, generated inputs, `trials.json`, and `summary.tsv`
are under `build/benchmark/network-current/2026-08-31-aggregate/`. An
independent recomputation of all 672 summary rows from `trials.json` matched
the recorded median, MAD, minimum, maximum, and sample count.

## Device counts and state accounting

| Case | Buses | Branch | REGCA | LoadZ | LoadZIP | GENROU | GENSAL | GenClassical |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ACTIVSg10k | 10,000 | 12,706 | 0 | 0 | 4,722 | 926 | 530 | 0 |
| ACTIVSg2000 | 2,000 | 3,206 | 0 | 1,507 | 0 | 432 | 0 | 0 |
| WECC240 | 243 | 447 | 37 | 0 | 141 | 103 | 0 | 0 |

Every differential-variable count is unchanged. All added states are
algebraic current variables.

| Case | Formulation | States | Differential | Algebraic | Jacobian nnz | nnz/state | Peak RSS KiB, median ± MAD [min, max] |
|---|---|---:|---:|---:|---:|---:|---:|
| ACTIVSg10k | direct | 83,124 | 27,391 | 55,733 | 340,872 | 4.1008 | 194,212 ± 28 [194,152, 194,276] |
| ACTIVSg10k | internal | 146,304 | 27,391 | 118,913 | 591,680 | 4.0442 | 252,584 ± 56 [252,504, 252,640] |
| ACTIVSg2000 | direct | 20,474 | 5,148 | 15,326 | 86,664 | 4.2329 | 56,716 ± 48 [56,624, 56,888] |
| ACTIVSg2000 | internal | 37,176 | 5,148 | 32,028 | 156,056 | 4.1978 | 77,820 ± 0 [77,816, 77,832] |
| WECC240 | direct | 5,425 | 1,804 | 3,621 | 19,686 | 3.6288 | 19,580 ± 12 [19,504, 19,620] |
| WECC240 | internal | 7,701 | 1,804 | 5,897 | 29,146 | 3.7847 | 21,808 ± 92 [21,716, 21,936] |

| Case | Measured state change | Required state change | Result |
|---|---:|---:|---|
| ACTIVSg10k | 83,124 -> 146,304 (+63,180) | +63,180 | 5/5 pairs pass |
| ACTIVSg2000 | 20,474 -> 37,176 (+16,702) | +16,702 | 5/5 pairs pass |
| WECC240 | 5,425 -> 7,701 (+2,276) | +2,276 | 5/5 pairs pass |

The earlier absolute estimates of 86,036 -> 149,216 and 21,338 -> 38,040
were each high by one two-current pair per affected machine: 2,912 states for
the 926 GENROU plus 530 GENSAL instances in ACTIVSg10k, and 864 states for the
432 GENROU instances in ACTIVSg2000. Those estimates retained the pre-reduction
machine-current pairs in the direct endpoint while also counting their
restoration in the delta. The measured deltas agree exactly with the intended
aggregate formulations.

## Application CPU time

| Case | Direct seconds, median ± MAD [min, max] | Internal seconds, median ± MAD [min, max] | Internal/direct |
|---|---:|---:|---:|
| ACTIVSg10k | 29.396600 ± 0.993600 [26.935400, 37.121800] | 56.154100 ± 1.892700 [54.261400, 65.524600] | 1.910x |
| ACTIVSg2000 | 2.010770 ± 0.340920 [1.558150, 2.953170] | 7.945210 ± 0.406890 [7.538320, 8.680020] | 3.951x |
| WECC240 | 0.153168 ± 0.004574 [0.148594, 0.171137] | 0.273016 ± 0.011158 [0.259367, 0.302899] | 1.782x |

The short WECC240 runs are especially sensitive to fixed overhead and host
jitter. The spread is reported rather than treating any single run as
representative.

## IDA work counters

The counters are identical across all five trials within each cell, but the
two formulations do not always follow the same nonlinear path.

| Case | Formulation | Steps | Residual evals | Jacobian evals | Linear setups | Error-test failures | Nonlinear iterations | Nonlinear convergence failures |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| ACTIVSg10k | direct | 1,508 | 2,390 | 191 | 191 | 67 | 2,384 | 50 |
| ACTIVSg10k | internal | 1,391 | 2,186 | 175 | 175 | 53 | 2,180 | 45 |
| ACTIVSg2000 | direct | 768 | 1,110 | 91 | 91 | 29 | 1,104 | 8 |
| ACTIVSg2000 | internal | 822 | 1,164 | 79 | 79 | 22 | 1,158 | 9 |
| WECC240 | direct | 602 | 799 | 65 | 65 | 23 | 793 | 1 |
| WECC240 | internal | 610 | 839 | 58 | 58 | 17 | 833 | 0 |

## Callback and KLU costs

| Case | Formulation | Operation | Calls | Total seconds, median ± MAD [min, max] | Microseconds/call, median ± MAD [min, max] |
|---|---|---|---:|---:|---:|
| ACTIVSg10k | direct | Residual | 2,390 | 15.459221 ± 0.663896 [14.025067, 20.212904] | 6,468.29 ± 277.78 [5,868.23, 8,457.28] |
| ACTIVSg10k | direct | Jacobian | 191 | 3.178898 ± 0.192608 [2.765985, 3.762069] | 16,643.45 ± 1,008.42 [14,481.60, 19,696.69] |
| ACTIVSg10k | direct | KLU setup | 191 | 1.605921 ± 0.068904 [1.537017, 1.884052] | 8,407.96 ± 360.75 [8,047.21, 9,864.15] |
| ACTIVSg10k | direct | KLU solve | 2,390 | 5.084565 ± 0.134349 [4.817115, 6.448539] | 2,127.43 ± 56.21 [2,015.53, 2,698.13] |
| ACTIVSg10k | internal | Residual | 2,186 | 18.326125 ± 0.744064 [17.582061, 21.473219] | 8,383.41 ± 340.38 [8,043.03, 9,823.06] |
| ACTIVSg10k | internal | Jacobian | 175 | 3.454909 ± 0.172307 [3.282603, 4.012000] | 19,742.34 ± 984.61 [18,757.73, 22,925.71] |
| ACTIVSg10k | internal | KLU setup | 175 | 11.974310 ± 0.525988 [11.448322, 14.023384] | 68,424.63 ± 3,005.65 [65,418.98, 80,133.62] |
| ACTIVSg10k | internal | KLU solve | 2,186 | 15.356269 ± 0.387695 [14.968574, 17.776812] | 7,024.83 ± 177.35 [6,847.47, 8,132.12] |
| ACTIVSg2000 | direct | Residual | 1,110 | 0.882922 ± 0.166823 [0.643211, 1.407096] | 795.42 ± 150.29 [579.47, 1,267.65] |
| ACTIVSg2000 | direct | Jacobian | 91 | 0.287728 ± 0.035707 [0.238589, 0.404712] | 3,161.84 ± 392.39 [2,621.85, 4,447.38] |
| ACTIVSg2000 | direct | KLU setup | 91 | 0.185172 ± 0.004246 [0.180926, 0.208227] | 2,034.86 ± 46.66 [1,988.19, 2,288.21] |
| ACTIVSg2000 | direct | KLU solve | 1,110 | 0.361953 ± 0.060535 [0.271638, 0.550776] | 326.08 ± 54.54 [244.72, 496.19] |
| ACTIVSg2000 | internal | Residual | 1,164 | 1.931482 ± 0.134183 [1.797300, 2.203542] | 1,659.35 ± 115.28 [1,544.07, 1,893.08] |
| ACTIVSg2000 | internal | Jacobian | 79 | 0.317519 ± 0.009150 [0.298685, 0.333003] | 4,019.23 ± 115.82 [3,780.83, 4,215.23] |
| ACTIVSg2000 | internal | KLU setup | 79 | 2.071297 ± 0.002634 [2.068663, 2.116257] | 26,218.94 ± 33.34 [26,185.61, 26,788.06] |
| ACTIVSg2000 | internal | KLU solve | 1,164 | 2.832281 ± 0.175797 [2.656484, 3.253970] | 2,433.23 ± 151.03 [2,282.20, 2,795.51] |
| WECC240 | direct | Residual | 799 | 0.061153 ± 0.002501 [0.058653, 0.066028] | 76.54 ± 3.13 [73.41, 82.64] |
| WECC240 | direct | Jacobian | 65 | 0.025010 ± 0.000958 [0.024052, 0.032454] | 384.78 ± 14.74 [370.03, 499.30] |
| WECC240 | direct | KLU setup | 65 | 0.013052 ± 0.000490 [0.012123, 0.013617] | 200.79 ± 7.54 [186.51, 209.49] |
| WECC240 | direct | KLU solve | 799 | 0.023897 ± 0.001233 [0.022551, 0.025548] | 29.91 ± 1.54 [28.22, 31.98] |
| WECC240 | internal | Residual | 839 | 0.072499 ± 0.003257 [0.066316, 0.078960] | 86.41 ± 3.88 [79.04, 94.11] |
| WECC240 | internal | Jacobian | 58 | 0.025452 ± 0.002308 [0.023144, 0.033302] | 438.83 ± 39.80 [399.04, 574.17] |
| WECC240 | internal | KLU setup | 58 | 0.056725 ± 0.001688 [0.054137, 0.058413] | 978.01 ± 29.11 [933.39, 1,007.12] |
| WECC240 | internal | KLU solve | 839 | 0.072977 ± 0.001699 [0.070039, 0.078917] | 86.98 ± 2.03 [83.48, 94.06] |

Difference-of-medians decomposition:

| Case | Added CPU seconds | Residual | Jacobian | KLU setup | KLU solve | KLU share of added CPU |
|---|---:|---:|---:|---:|---:|---:|
| ACTIVSg10k | 26.757500 | +2.866904 | +0.276011 | +10.368389 | +10.271704 | 77.14% |
| ACTIVSg2000 | 5.934440 | +1.048561 | +0.029791 | +1.886125 | +2.470328 | 73.41% |
| WECC240 | 0.119848 | +0.011345 | +0.000442 | +0.043673 | +0.049080 | 77.39% |

The internal formulation is slower even where it performs fewer KLU setups.
This separates the cost-per-call effect from changes in IDA work counts.

## Model-family timers

Times are median ± MAD [min, max] seconds. A zero means median, MAD, minimum,
and maximum are all zero because that family has no instances in the case. The
same REGCA implementation has different total WECC240 time because the
formulations have different residual and Jacobian call counts.

### ACTIVSg10k

| Family | Direct residual | Internal residual | Direct Jacobian | Internal Jacobian |
|---|---:|---:|---:|---:|
| REGCA | 0 | 0 | 0 | 0 |
| Branch | 6.584044 ± 0.287110 [6.020928, 8.502310] | 8.082607 ± 0.345771 [7.736836, 9.407399] | 0.797707 ± 0.057283 [0.678671, 0.945491] | 0.777629 ± 0.047298 [0.727599, 0.894155] |
| LoadZ | 0 | 0 | 0 | 0 |
| LoadZIP | 1.438129 ± 0.078507 [1.242806, 1.966515] | 2.242297 ± 0.103155 [2.139142, 2.736322] | 0.227297 ± 0.025283 [0.175414, 0.292481] | 0.268062 ± 0.020836 [0.245412, 0.327194] |
| GENROU | 0.816060 ± 0.031100 [0.730815, 1.100289] | 0.867960 ± 0.034227 [0.833733, 1.037006] | 0.149818 ± 0.007717 [0.136189, 0.175033] | 0.232186 ± 0.004934 [0.226686, 0.253722] |
| GENSAL | 0.485305 ± 0.024765 [0.432864, 0.641248] | 0.498455 ± 0.009211 [0.489244, 0.597655] | 0.092332 ± 0.003899 [0.083916, 0.105604] | 0.073176 ± 0.003841 [0.069334, 0.086569] |
| GenClassical | 0 | 0 | 0 | 0 |

### ACTIVSg2000

| Family | Direct residual | Internal residual | Direct Jacobian | Internal Jacobian |
|---|---:|---:|---:|---:|
| REGCA | 0 | 0 | 0 | 0 |
| Branch | 0.278229 ± 0.048057 [0.205525, 0.454328] | 0.739956 ± 0.048155 [0.691800, 0.832735] | 0.056527 ± 0.005349 [0.048098, 0.080095] | 0.055801 ± 0.001921 [0.051854, 0.058375] |
| LoadZ | 0.081499 ± 0.012419 [0.063659, 0.116506] | 0.232486 ± 0.016873 [0.215613, 0.263926] | 0.014484 ± 0.001845 [0.011149, 0.020134] | 0.019230 ± 0.000857 [0.018064, 0.021445] |
| LoadZIP | 0 | 0 | 0 | 0 |
| GENROU | 0.080506 ± 0.014760 [0.057980, 0.131740] | 0.149540 ± 0.011149 [0.138391, 0.173586] | 0.029161 ± 0.002013 [0.026679, 0.035212] | 0.028771 ± 0.000477 [0.027744, 0.029248] |
| GENSAL | 0 | 0 | 0 | 0 |
| GenClassical | 0 | 0 | 0 | 0 |

### WECC240

| Family | Direct residual | Internal residual | Direct Jacobian | Internal Jacobian |
|---|---:|---:|---:|---:|
| REGCA | 0.004722 ± 0.000143 [0.004579, 0.005634] | 0.005132 ± 0.000202 [0.004848, 0.005763] | 0.001928 ± 0.000044 [0.001884, 0.002317] | 0.001763 ± 0.000093 [0.001669, 0.002137] |
| Branch | 0.015481 ± 0.000387 [0.014482, 0.015868] | 0.017768 ± 0.000622 [0.016617, 0.020002] | 0.003590 ± 0.000399 [0.003151, 0.004006] | 0.003134 ± 0.000129 [0.003005, 0.003712] |
| LoadZ | 0 | 0 | 0 | 0 |
| LoadZIP | 0.002011 ± 0.000138 [0.001851, 0.002394] | 0.003611 ± 0.000079 [0.003532, 0.003959] | 0.000494 ± 0.000023 [0.000472, 0.000657] | 0.000798 ± 0.000016 [0.000782, 0.000913] |
| GENROU | 0.005927 ± 0.000304 [0.005624, 0.006333] | 0.006583 ± 0.000459 [0.006118, 0.007597] | 0.003865 ± 0.000260 [0.003606, 0.004678] | 0.003314 ± 0.000223 [0.003091, 0.004283] |
| GENSAL | 0 | 0 | 0 | 0 |
| GenClassical | 0 | 0 | 0 | 0 |

## Interpretation

The internal-current formulation increases state count by 76.0% on
ACTIVSg10k, 81.6% on ACTIVSg2000, and 42.0% on WECC240. Jacobian nonzeros rise
by 73.6%, 80.1%, and 48.1%; median peak RSS rises by 30.1%, 37.2%, and 11.4%.
Residual and Jacobian callback cost also rises, but the largest aggregate
increase is KLU setup and solve time.

The larger sparse systems and higher KLU per-call costs are consistent with
more expensive factorization and triangular solves. Sparse-factor fill was not
measured, however, so fill and cache pressure remain hypotheses rather than
established causes. Confirming them requires KLU symbolic/numeric factor
statistics or equivalent sparse-factor instrumentation. Timers are
hierarchical and must not be summed as independent shares.

WSL2 introduced visible spread, especially in the high ACTIVSg10k maxima and
the direct ACTIVSg2000 cell. Ratios above use medians from the five measured
trials. The model-family totals also reflect differing callback counts and are
not pure per-model microbenchmarks.

## FlameGraph status

`run_flamegraphs.rb` implements the upstream
[capture, fold, and render workflow](https://github.com/brendangregg/FlameGraph/blob/master/README.md)
and the normalized
[differential FlameGraph method](https://www.brendangregg.com/blog/2014-11-09/differential-flame-graphs.html).
It verifies the recorded binary and input hashes, performs one unprofiled
warm-up and three interleaved captures per formulation/case, captures
user-space samples with
`perf record -F 99 -e cpu-clock:u --call-graph dwarf,16384`, folds with
`stackcollapse-perf.pl --all`, renders stable `--hash` colors, and uses
`difffolded.pl -n` so differential widths are normalized to the internal-current
variant.

FlameGraph itself is not vendored. The driver requires an external checkout at
`41fee1f99f9276008b7cd112fca19dc3ea84ac32`.

No profiles or SVGs were produced on this host because `perf` is absent. The
prerequisite check failed before any capture with:

```text
perf executable not found; install a perf build matching this Linux/WSL kernel or pass --perf PATH. No profile captures were started.
```

Use matching WSL perf tooling or a native Linux host without changing the
recorded binaries or generated inputs. The driver will then reject empty
captures, unresolved GridKit/IDA/KLU frames, noninteractive SVGs, or more than
the configured share of samples containing unknown frames.
