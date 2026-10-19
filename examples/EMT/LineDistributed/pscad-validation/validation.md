# GridKit LineDistributed against PSCAD: validation summary

Peak errors are max |GridKit - PSCAD| over the three phases and the study window, normalized by the peak |PSCAD| of the same quantity. RMS errors are normalized by the PSCAD RMS. GridKit steps are the accepted adaptive steps; h/tau is the accepted step over the shortest modal delay (median) and 'time above all delays' is the fraction of elapsed time with steps longer than every modal delay. PSCAD steps are tf over the fixed step.

| Line | Model | Scenario | v1 peak % | v2 peak % | source i peak % | i2 peak % | v2 RMS % | i2 RMS % | PSCAD h/2 effect on v2 % | GridKit steps | median h/tau | time above all delays | GridKit s | PSCAD steps |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 69kv-wood-pole | fd | energize | 1.342 | 1.462 | 1.787 | 0.876 | 0.599 | 0.831 | 0.744 | 66618 | 0.01 | 3% | 2.5 | 31000 |
| 69kv-wood-pole | fd | pulse | 18.356 | 13.742 | 18.469 | 0.668 | 4.800 | 0.240 | 7.095 | 6435 | 0.00 | 27% | 0.3 | 7000 |
| 69kv-wood-pole | fd | fault | 5.423 | 2.016 | 0.464 | 0.187 | 0.848 | 0.148 | 1.046 | 109974 | 0.01 | 0% | 2.7 | 60000 |
| 69kv-wood-pole | fd | sag | 1.335 | 1.967 | 1.910 | 0.040 | 0.196 | 0.013 | - | 108341 | 0.01 | 0% | 2.7 | 50000 |
| 138kv-delta | fd | energize | 1.163 | 1.270 | 1.273 | 0.614 | 0.475 | 0.576 | 0.648 | 41476 | 0.00 | 2% | 1.9 | 31000 |
| 138kv-delta | fd | pulse | 18.347 | 11.383 | 18.457 | 0.575 | 4.483 | 0.216 | 6.070 | 8417 | 0.00 | 21% | 0.4 | 7000 |
| 138kv-delta | fd | fault | 15.298 | 1.813 | 0.615 | 0.234 | 1.086 | 0.078 | 1.012 | 94599 | 0.01 | 0% | 2.7 | 60000 |
| 138kv-delta | fd | sag | 1.138 | 1.742 | 1.685 | 0.032 | 0.135 | 0.010 | - | 68162 | 0.01 | 0% | 2.0 | 50000 |
| 345kv-h-frame | fd | energize | 12.346 | 14.012 | 19.331 | 14.195 | 11.200 | 11.365 | 0.315 | 29211 | 0.00 | 2% | 2.0 | 31000 |
| 345kv-h-frame | fd | pulse | 18.375 | 22.025 | 18.695 | 31.390 | 8.165 | 12.073 | 5.557 | 7401 | 0.00 | 14% | 0.3 | 7000 |
| 345kv-h-frame | fd | fault | 9.477 | 6.149 | 3.425 | 0.644 | 5.863 | 0.405 | 0.587 | 44944 | 0.00 | 0% | 1.7 | 60000 |
| 345kv-h-frame | fd | sag | 3.751 | 4.622 | 6.042 | 0.814 | 1.265 | 0.195 | - | 47618 | 0.00 | 0% | 1.7 | 50000 |
| 765kv-horizontal | fd | energize | 6.925 | 7.507 | 8.629 | 7.567 | 5.324 | 5.400 | 0.144 | 12211 | 0.00 | 0% | 1.6 | 31000 |
| 765kv-horizontal | fd | pulse | 18.384 | 53.273 | 18.679 | 57.955 | 23.802 | 26.341 | 4.223 | 8118 | 0.00 | 0% | 0.9 | 17563 |
| 765kv-horizontal | fd | fault | 27.970 | 8.258 | 3.330 | 17.784 | 6.203 | 3.594 | 3.141 | 49034 | 0.00 | 0% | 2.3 | 60000 |
| 765kv-horizontal | fd | sag | 1.891 | 2.127 | 2.678 | 0.598 | 0.696 | 0.195 | - | 25011 | 0.00 | 0% | 1.6 | 50000 |

## Pre-event steady state (sinusoidal cases, 0 to 20 ms)

| Line | Model | Scenario | v2 peak % | i2 peak % | source i peak % |
|---|---|---|---:|---:|---:|
| 69kv-wood-pole | fd | fault | 0.074 | 0.001 | 0.076 |
| 69kv-wood-pole | fd | sag | 0.074 | 0.001 | 0.076 |
| 138kv-delta | fd | fault | 0.073 | 0.003 | 0.075 |
| 138kv-delta | fd | sag | 0.073 | 0.003 | 0.075 |
| 345kv-h-frame | fd | fault | 0.080 | 0.015 | 0.339 |
| 345kv-h-frame | fd | sag | 0.080 | 0.015 | 0.339 |
| 765kv-horizontal | fd | fault | 0.088 | 0.117 | 0.928 |
| 765kv-horizontal | fd | sag | 0.088 | 0.117 | 0.928 |
