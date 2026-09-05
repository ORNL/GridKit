# CoupledGrid EMT case

A second, larger frequency dependent EMT example: 14 buses at 13.8 kV and
60 Hz, thirteen 300–600 m overhead-line sections, two synchronous machines
with TGOV1 governors and IEEET1 exciters, three PWM–Converter–DVS sources,
seven unbalanced RL loads, and one switched resistive load. The eighth ring
bus connects back to the grid source; five spurs connect remote loads and
converters. Bus 14 connects to bus 10 at 0.5 s and disconnects at 0.8 s.
The simulation ends at 1.2 s.

The topology and dispatch are synthetic. The line parameters come from a
physical, untransposed ACSR conductor geometry calculated with the local
OpenLine tool. See [line provenance and validation](lines/README.md).
Every line uses full 3-by-3 rational impedance and admittance matrices,
including mutual impedance and capacitance. The eight shared poles describe
skin effect and lossy earth return; the Maxwell capacitance is represented
exactly. The fitted impedance is passive by construction. Its maximum
held-out entrywise error is 0.2292% over 1 Hz–30 kHz.

The grid source and converter reactors use exact RL admittances expressed as
VectorFit submodels. The source has R = 0.12 ohm and L = 1 mH per phase.
Each converter reactor has R = 5.0784 ohm and L = 0.101031558 H, corresponding
to 0.02 pu resistance and 0.15 pu reactance on a 0.75 MVA base.
Machine ratings are 5 and 4 MVA, with initial active dispatch of 1.2 and 1 MW.
Connected loads total 4.1 MW at nominal voltage and 0.97 power factor; their
phase powers vary by up to 6%. The switched load adds nominally 0.8 MW.
Actual load power varies with voltage because these are impedance loads.

## Mu comparison

The common carrier frequency is 900 Hz, modulation index is 0.8, and the
DC voltage is fixed at 30,306.481793 V in every run. The three settings are

| Study | Mu | Definition |
| --- | ---: | --- |
| `low.solver.json` | 240 | 4 times 60 Hz |
| `middle.solver.json` | 929.516003 | Geometric midpoint |
| `high.solver.json` | 3600 | 4 times the 900 Hz carrier |

The shared CommonMath parameter also changes smoothing in machine saturation
and controller limit equations. In PWM, it acts on time differences in seconds.
Changing mu changes the differential equations, including the applied inverter
fundamental voltage; it is not a solver tolerance. DC voltage is deliberately
held fixed to expose this effect. The high-mu case is a comparison baseline,
not the hard-switching limit or an accuracy reference.

The exact ideal-pulse Fourier coefficient is attenuated by

```math
A(f,\mu)=\frac{2\pi^2 f/\mu}{\sinh(2\pi^2 f/\mu)}.
```

At 60 Hz the three gains are 0.070985, 0.773148, and 0.982186.
The independent calculation in `pwm_analysis.py` includes sampling and pulse
alignment. Converter common-mode removal cancels the 900 Hz triplen carrier;
780 Hz and 1020 Hz sidebands remain useful comparisons.

These are open-loop ideal-DC converter sources with RL reactors. There is no
PLL, current control, DC-link dynamics, or protection. At low mu their reduced
internal fundamental causes large reactive absorption. This demonstrates the
smoothing parameter's effect, not credible protected IBR operation at that
setting. The untransposed line and unequal loads excite phase coupling.

## Run and render

Use an EMT build with Enzyme, SUNDIALS, and sparse KLU. NumPy, SciPy, and
Matplotlib are needed for case generation, analysis, and plotting.

```bash
cmake --build build/emt-rational --target EMTDynamicSimulation -j 10
python3 cases/EMT/CoupledGrid/build_case.py
python3 cases/EMT/CoupledGrid/run.py
python3 cases/EMT/CoupledGrid/analyze.py --case cases/EMT/CoupledGrid/CoupledGrid.case.json --results build/emt-coupled-results
python3 cases/EMT/CoupledGrid/report.py --results build/emt-coupled-results
python3 cases/EMT/CoupledGrid/draw_oneline.py --case cases/EMT/CoupledGrid/CoupledGrid.case.json --output build/emt-coupled-results --line-data cases/EMT/CoupledGrid/lines/response.csv
python3 cases/EMT/CoupledGrid/plot_results.py --case cases/EMT/CoupledGrid/CoupledGrid.case.json --results build/emt-coupled-results
```

`run.py --binary /path/to/EMTDynamicSimulation --output /path/to/results`
selects another build or output location. It performs three sequential trials
per mu in rotated order, then repeats each mu with tenfold tighter tolerances.
Repeated primary CSVs must be byte-identical. Logs, solver files, raw monitor
CSVs, executable/library/case/state hashes, application CPU time, wall time, and accumulated
IDA counters are retained. The application reports counters summed across
event reinitializations; `linear_setups` counts IDA linear-solver setups.
Tight runs are excluded from runtime medians.

For a single run, invoke the EMT application with one of the supplied solver
files. Physical case, initial state, and numerical study configuration are
separate JSON files. `build_case.py` regenerates them from `lines/model.json`
and a full three-phase 60 Hz network solve. No propagation or delay model is
required by these lumped pi sections.

## Interpreting the results

Primary tolerances are relative `1e-6`, absolute `1e-7`; tighter runs use
`1e-7` and `1e-8`. Monitoring is at 60 kHz, independently of IDA's adaptive
internal steps. Coherent 0.1 s measurement windows are before switching
`[0.4, 0.5)`, with the load connected `[0.7, 0.8)`, and after opening
`[1.1, 1.2)`. These are finite-window measurements, not claims of equilibrium.

All studies share the same fundamental initial estimate. Machine dispatch
and bus voltages are initialized, while internal line/filter states use their
model defaults. Consequently startup contains artificial energization
transients. The small physical shunt capacitances give network modes above
the line-fit band and monitor Nyquist frequency. IDA resolves those modes
internally, but the recorded startup cannot validate their waveforms.
Full-response plots mark the first 0.2 s; comparisons of operating quantities
use the later windows. Ideal switch edges also excite frequencies beyond the
validated band, so sub-millisecond switching peaks are model results without
a physical high-frequency accuracy claim.

The result directory contains a rendered one-line, physical geometry and
matrix plots, overlays using one color per mu, waveform differences relative
to mu = 3600, harmonic spectra, timing/solver-work comparisons, and the raw
simulation evidence. `index.html` links the gallery; `metrics.json` records
KCL, independent PWM checks, operating quantities, and tolerance comparisons.
Generated CSVs and rendered results belong in the build directory.
