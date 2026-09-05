# Ten-bus EMT converter scenarios

These six studies share [one case](../../../cases/EMT/IBR/TenBus.case.json) and
[one initial-state estimate](../../../cases/EMT/IBR/TenBus.state.json).
See the [case description](../../../cases/EMT/IBR/README.md) for topology,
parameters, initialization, and the explicitly compensated PWM smoothing.

| Solver file | Events |
|---|---|
| `01_Baseline.solver.json` | No events; startup and undisturbed operation |
| `02_LoadStep.solver.json` | Connect the 2 MW nominal load at 1.0 s |
| `03_LoadPulse.solver.json` | Connect that load at 1.0 s; shed it at 1.5 s |
| `04_FaultClearing.solver.json` | Apply a three-phase 2 Ω/phase shunt at 1.0 s; clear at 1.06 s |
| `05_TieReclosing.solver.json` | Open the 7–8 tie at 1.0 s; reclose at 1.8 s |
| `06_FaultClearingSwitching.solver.json` | Same fault/clearing times as 04, with `mu=50000`, 10 µs monitoring and explicitly adjusted DC constants |

Every run spans 0–3 s with adaptive IDA stepping and sparse KLU. The first
five use 100 µs monitor spacing and default `mu=240`; the sixth uses 10 µs.
All use relative tolerance 1e-7 and absolute tolerance 1e-9. Monitoring spacing is not
the internal integration step. Ideal-switch events retain both pre-event and
post-event samples at the same timestamp.

## Run and review

Build `EMTDynamicSimulation` with Enzyme and SUNDIALS KLU enabled, then run
from the repository root:

```bash
python3 examples/EMT/IBR/run.py --exe build/application/EMT/EMTDynamicSimulation --jobs 4
python3 examples/EMT/IBR/plot.py
```

Use an absolute executable path when outside the repository. The runner resolves
case/state paths relative to the example files and isolates each scenario's
outputs. `--scenario 04_FaultClearing` selects a run, `--check-only` validates
existing outputs, and `--results PATH` selects another output directory (pass
it to both scripts). The runner requires only Python's standard library;
plotting and case regeneration require NumPy and Matplotlib.

Open [results/index.html](results/index.html) for the one-line diagram, scenario
comparison, and links to each scenario's plots and raw data. The HTML is local
and requires no server or internet connection. PNG figures are accompanied by
SVG for the summary plots and PDF for the one-line and complete DAE traces.
Generated data and plots live under ignored `results/`; no generated output is
part of the source case. Expect roughly 6 GB of CSV data for all six runs.

A solver JSON also runs directly with `EMTDynamicSimulation`; its output paths
are relative to the current working directory. Build and install trees preserve
the same relative case paths. CTest `EMTTenBusIBR` runs all six scenarios with
a shortened time horizon and checks their events and output completeness.

## Recorded data and checks

Each full run records 222 model quantities plus time, including all available
Bus, Machine, PWM, Converter, dependent-source, line, load and switch monitors.
The complete state CSV separately records all 231 DAE variables (78 differential, 153 algebraic) and all 231
derivatives, including the three internal states of each TGOV1 governor.
The companion `.states.csv.json` maps global columns to component classes,
paths, local indices, and differential/algebraic classification. Consult each
component's internal-variable enum for local-index meaning and native units.
Algebraic derivatives are solver interpolants, not independent physical states.

The runner rejects nonfinite values, incomplete rows, missing final times,
incomplete state index maps, and incorrect pre/post-event switch commands.
It records input hashes, timings and sample counts in `validation.json`.
The plotting script independently reconstructs every bus's phase KCL and
checks the converter identity `vo = vdc (s - mean(s))` and PWM bounds. It also
plots network power balance including line resistance losses and line magnetic
and capacitive energy rates from recorded states/derivatives. Converter power
is measured at the AC bus, so filter losses lie upstream of that boundary.
These checks assess numerical/model consistency, not external validation.

Displayed RMS and power traces use a trailing 1/60 s trapezoidal window; incomplete
startup windows are left blank. Positive source P/Q denotes injection into the
grid. LoadZ records current injection (negative for consumption); load-power
plots reverse that sign to show positive consumption. Reactive power uses the amplitude-invariant Clarke transform. Frequency
plots show machine rotor speed times 60 Hz, not PLL-estimated bus frequency.
Rotor-angle deviation is relative to the initialized angle advancing at 60 Hz.
Spectra use a coherent final 0.5 s window at each study's actual sample rate.
The overview includes direct PWM, bridge-voltage, filter-current, and spectral
comparisons between the two fault studies. Overview metrics exclude startup
before 1 s; the gallery and raw data retain the complete startup transient.

Every DAE variable has a trace in the per-scenario state pages and multipage
PDF. Those dense pages display every tenth recorded sample for readability;
the original CSVs retain every sample and every derivative.

## Resolved switching study

`06_FaultClearingSwitching.solver.json` sets `mu=50000` before model
construction. The logistic 10–90% transition width is `2 ln(9)/mu`, about
87.9 µs, resolved by the 10 µs output interval. The carrier remains 900 Hz.

Its `signal_values` override sets each DC source to 28918.846170570516 V.
This replaces the default case's 407357 V smoothing compensation and preserves
the 1.02 pu open-circuit AC fundamental; it is not an automatic rescaling by
the solver. Changing only mu while retaining 407357 V would apply excessive
AC voltage. The shared mu also changes machine saturation and governor
limiter smoothing, so this is a study-wide smoothing comparison.

All full simulations use sparse KLU. Earlier short dense-versus-Enzyme
cross-checks are retained as historical verification data; the EMT driver now
rejects dense fallback, and subsequent simulations use sparse exclusively.
