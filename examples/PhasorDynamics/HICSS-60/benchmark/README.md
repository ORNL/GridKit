| Setting | Value |
|---|---|
| Duration | 10 s |
| Fault R | 0 p.u. |
| Fault X | 0.01 p.u. |
| Fault on | 1 s |
| Fault off | 1.15 s |
| Timing | Median of 7 CPU times |
| Warm-up | 1 per point |
| Monitoring | Off |
| CPU affinity | 14 |

| Case | `mu` | `rel_tol` | `abs_tol` |
|---|---|---|---|
| NE | 240 | 3e-5 | 2e-5 |
| Illinois | 240 | 3e-5 | 2e-6 |
| Hawaii | 240 | 1e-6 | 7e-5 |
| WECC | 240 | 2e-6 | 1e-6 |
| Texas | 240 | 1e-4 | 1e-6 |

## Results

| | Variables | | Runtime (s) | | | $\epsilon_{\mathrm{RMSE}}^{\text{abs}}$ | | $\epsilon_{\infty}^{\text{rel}}$ | |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| Case | Diff. | Alg. | PowerWorld | GridKit | Speedup | $\omega$ | $\lvert V\rvert$ | $\omega$ | $\lvert V\rvert$ |
| IEEE N.E. | 180 | 266 | 6.969 | 0.011 | 615.0x | 8.52e-6 | 1.56e-4 | 5.78e-3 | 1.14e-3 |
| Synth. Illinois | 400 | 482 | 4.032 | 0.029 | 137.6x | 5.82e-6 | 3.91e-5 | 4.08e-3 | 7.90e-4 |
| Synth. Hawaii | 602 | 668 | 2.172 | 0.066 | 32.8x | 2.95e-5 | 2.16e-4 | 6.19e-3 | 4.44e-3 |
| WECC | 1804 | 2451 | 4.375 | 0.115 | 38.1x | 9.26e-7 | 1.06e-5 | 5.70e-3 | 1.25e-4 |
| Synth. Texas | 6056 | 8314 | 14.703 | 0.993 | 14.8x | 6.04e-6 | 1.91e-4 | 1.43e-2 | 3.12e-3 |
