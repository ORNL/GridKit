# Hawaii EMT fault study

The [case README](../../../cases/EMT/Hawaii/README.md) documents the conversion,
fabricated data, and deviations from the validated PhasorDynamics case.
The EMT fault runs from 1.0 to 1.10 s. The PowerWorld reference fault runs from
1.0 to 1.15 s. Both intervals are labelled in the plots.

The four figures contain all 37 buses or all 30 synchronous machines, followed
by the range of differences across matching channels. Bus voltage uses a
one-cycle positive-sequence estimate. Machine speed, active power, and reactive
power use cycle means; their colours identify the six machine buses. The time
label is the centre of the averaging window;
the approximately half-cycle transition spreading is a measurement effect.
Powers use the 100 MVA system base. These are comparisons between different
dynamic models, not a new PowerWorld validation.

The inverter plants use the same LCL Filter wiring as the GFL control example.
The generator preserves initial terminal dispatch, includes both reactor losses
in DC power, and initializes capacitor voltage and both inductor currents.
The filter and controller parameters are synthetic, as listed in the case README.

Generate current results using the commands below. Plots, editable TeX/data,
and measured validation/comparison statistics are written under ignored
`results/`; no numerical snapshot is assumed to describe a regenerated case.

Differences from the reference include the unequal fault intervals. The
resistive EMT fault produces an active-power surge and initial machine deceleration; the inductive
reference fault produces a different initial response. The comparison therefore
does not establish equivalent disturbance dynamics, despite similar recovery
trends. The [deviation list](../../../cases/EMT/Hawaii/README.md#fault-and-comparison)
also covers the winding approximation and replacement inverter controls.

`results/Hawaii.{vmag,omega,p,q}.pdf` contains the four vector figures;
matching PNGs are rendered at 600 DPI. Each figure retains its editable
`.tex` source and `.csv` data. The figures remain on the TeX editing path.
`results/comparison.json` records per-channel errors against the frozen
PowerWorld reference. The validated run's `metrics.json` records physical
checks, input/executable/library hashes, and accepted-step statistics.

`convergence.py` compares the full run with a tenfold tighter run through
1.5 s, covering the disturbance and early recovery. It requires matching
physical inputs and records maximum differences in `results/convergence.json`.

`switching.py` checks all nine bridges in a separate one-cycle run sampled at
720 kHz, using the same physical case and solver tolerances. It compares the
switching functions and bridge harmonics with an independent periodic sigmoid
sum at the instantaneous duty, and checks bridge AC/DC power balance using
Filter's converter-side current. `results/switching.json` records fundamental,
carrier, and sideband amplitudes. The window includes startup; these checks
establish the continuous switching equations and DC conversion, not a
steady-state harmonic-performance specification.

From the repository root, with Enzyme and SUNDIALS KLU enabled:

```bash
python3 cases/EMT/Hawaii/convert.py
cmake --build build --target EMTDynamicSimulation -j 10
python3 cases/EMT/Hawaii/validate.py \
  --exe build/application/EMT/EMTDynamicSimulation \
  --tmax 10 --output /tmp/gridkit-hawaii-run
python3 examples/EMT/Hawaii/plot.py \
  --averaged /tmp/gridkit-hawaii-run/Hawaii.averaged.csv
python3 cases/EMT/Hawaii/validate.py \
  --exe build/application/EMT/EMTDynamicSimulation \
  --tmax 1.5 --rel-tol 1e-6 --abs-tol 1e-7 --output /tmp/gridkit-hawaii-refined
python3 examples/EMT/Hawaii/convergence.py \
  --baseline /tmp/gridkit-hawaii-run --refined /tmp/gridkit-hawaii-refined
python3 examples/EMT/Hawaii/switching.py \
  --exe build/application/EMT/EMTDynamicSimulation
ctest --test-dir build -R EMTHawaiiCase --output-on-failure
```

The CTest study ends at 1.5 s and covers fault inception, clearing, and early
recovery. The ten-second run checks the longer response separately.
For a quicker plot, use `--tmax 1.5` in the first validation command; the
plotter follows the validated run duration.
Plot generation reads the PowerWorld CSVs with `git show` at the same frozen
source revision used by the converter; no checkout is needed. It requires
PGFPlots, pdfLaTeX, and Poppler. Simulation and validation require only Python's
standard library in addition to the EMT executable.
