# DOE

| Case | GridKit | Fixed-Step | Speedup | $\epsilon_{\mathrm{RMSE}}^{\text{abs}}$ | $\epsilon_{\infty}^{\text{rel}}$ |
|:--|--:|--:|--:|--:|--:|
| Hawaii | 0.055 [s] | 2.172 [s] | 39.8× | 2.73e-05 | 5.71e-03 |
| IEEE39 | 0.010 [s] | 6.969 [s] | 730.6× | 8.55e-06 | 5.77e-03 |
| ACTIVSg200 | 0.021 [s] | 4.032 [s] | 191.8× | 5.60e-06 | 4.04e-03 |
| WECC240 | 0.101 [s] | 4.375 [s] | 43.4× | 8.20e-07 | 2.79e-03 |
| ACTIVSg2000 | 0.654 [s] | 14.703 [s] | 22.5× | 6.04e-06 | 1.40e-02 |
| ACTIVSg10k | 6.989 [s] | 62.735 [s] | 9.0× | 9.39e-06 | 3.23e-02 |
| ACTIVSg25k | 16.013 [s] | 109.500 [s] | 6.8× | - | - |
| ACTIVSg70k | 59.442 [s] | 385.562 [s] | 6.5× | - | - |

GridKit: CPU time of one run with monitoring off. Study settings are in [solvers](solvers/).

Errors: generator speed against PowerWorld from [Validation](../Validation/).
