# Current canonical case sweep, 2026-08-26

This study runs the canonical ACTIVSg200, ACTIVSg500, WECC240, and ACTIVSg10k
cases without copying their case JSON files.

All four solver definitions use a 10-second horizon, a fault from 1.0 to 1.1
seconds, relative tolerance `1e-5`, absolute tolerance `1e-7`, AMD ordering,
and no periodic monitor output. The runner reports the median of three
interleaved trials using the application's `Complete in` time.

```bash
ruby benchmark/case-sweep/2026-08-26-current-cases/run_bench.rb
```

## Results

These are median application-reported CPU times from three interleaved trials.

| Case | Buses | States | Jacobian nonzeros | `Complete in` | Steps |
|---|---:|---:|---:|---:|---:|
| ACTIVSg200 | 200 | 882 | 5,616 | 0.029631 s | 860 |
| ACTIVSg500 | 500 | 1,722 | 10,956 | 0.045682 s | 656 |
| WECC240 | 243 | 4,292 | 18,656 | 0.111918 s | 609 |
| ACTIVSg10k | 10,000 | 73,998 | 333,732 | 11.7716 s | 1,532 |

The raw local trials are written to `results.tsv` and should remain untracked.

## ACTIVSg10k load representation

All 4,722 ACTIVSg10k `LoadZIP` devices had `alphaI=alphaP=0`, making them pure
constant-impedance loads. The conversion script replaces them with equivalent
`LoadZ` devices using the initialized bus-voltage magnitude:

```text
R = Pnom * V0^2 / (Pnom^2 + Qnom^2)
X = Qnom * V0^2 / (Pnom^2 + Qnom^2)
```

The median ACTIVSg10k runtime fell from 14.4478 seconds to 11.7716 seconds, a
1.23x speedup. In a representative profile, the load residual bucket fell from
0.977731 seconds to 0.000039 seconds.
