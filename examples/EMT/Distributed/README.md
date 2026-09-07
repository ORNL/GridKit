# Distributed EMT line studies

These scripts run `EMTDynamicSimulation` with adaptive IDA stepping and sparse
KLU. They compare a lossless Bergeron line with one and ten π sections, then
compare frequency-dependent distributed lines with 60 Hz π equivalents in
synthetic 8- and 20-bus networks. Generated cases, raw waveforms, accepted-step
logs, resource measurements, and plots remain in `local/emt-distributed`.

## Square-pulse example

The focused pulse study uses the same 60 km overhead geometry with a 1 kV,
100 µs square source-voltage pulse, a 100 Ω source resistance, and a 600 Ω receiving
resistance per phase. Equal phase excitation isolates the zero-sequence
response. The line starts unenergized, and the simulation retains 2 ms after
pulse launch.

The dashed trace is the prescribed source voltage $e(t)$. The current into
the line is $i(t)=[e(t)-v_s(t)]/(100\,\Omega)$ per phase, where $v_s$ is the
computed sending-terminal voltage.

After building the application and generating the `fits-final` coefficients
as described below, run:

```sh
python3 examples/EMT/Distributed/pulse.py
```

Open `local/emt-distributed/pulse/index.html` for three figures: dimensioned
geometry and sag, the three line discretizations, and sending/receiving
pulses. A table reports the response peaks and spatial convergence. Each
figure is saved as SVG, PDF, and PNG. The geometry and circuits also include
editable TeX sources and data, using the shared EMT README diagram style.
Rendering those diagrams requires `pdflatex`, `circuitikz`, `dvisvgm`, and
`pdftoppm`. Raw cases, solver logs, coefficient files, waveforms,
metrics, validation details, and input hashes remain in the same directory.
Use `--output` to retain a separate study, `--fits` to select its source fits,
and `--stage plot` to regenerate figures and reference checks from retained
simulations. Add `--gridworkbench /path/to/GridWorkbench` to recompute the
geometry on an independent frequency grid and check the new series fits
between the original fitting samples.

Use `--stage diagrams` to redraw only the geometry, circuits, and gallery
from retained results, preserving the waveform figure and simulation data.
The circuit figure shows `n=1,5,20` on a common 60 km horizontal scale,
with section widths proportional to `1/n`. Its series R(ω) and L(ω) symbols
represent the fitted frequency-dependent impedance. End shunts are C/2 and
interior shunts are C; each lower node denotes the common reference potential.
The waveform figure compares the same three resolutions with the distributed
model over 0–1000 µs after launch. The complete five-model data and full 2 ms
validation window remain in the results files.

The distributed model uses the existing `Yc` and `H` fits. Its 1-, 5-, 20-, and
100-section lumped comparisons use frequency-dependent `Zp` and `Yp`, rather
than the older study's 60 Hz π equivalents. Series impedance is fitted from
the same `parameters.npz` samples as a sum of passive modal RL branches,
with positive modal resistance and inductance. The shunt capacitance is exact
for this geometry. The script checks both modal series fits against a 0.02%
maximum relative-error limit over the supplied frequency samples.

All five EMT circuits run at nominal and ten-times tighter tolerances using
adaptive IDA and sparse KLU. Independent Laplace-domain circuit solutions
check each simulation; halving the inverse-transform spacing checks the
references. The RMS spatial comparison uses the complete 2 ms response
window. Numerical checks report errors at the source discontinuities
separately. An ideal square pulse has unlimited bandwidth, whereas the source
geometry fits cover 0.01 Hz–10 MHz. See the generated report for measured
errors and the inherited distributed fit's passivity limitations.

## Run

Build the EMT application with Enzyme and sparse SUNDIALS enabled:

```sh
cmake --build build --target EMTDynamicSimulation -j 10
python3 examples/EMT/Distributed/run.py --trials 3
python3 examples/EMT/Distributed/run.py --sine --short --models distributed distributed_capped pi_1 --trials 3
```

The lossless studies need no external fitting package. Install the packages
in `requirements.txt` to generate frequency-dependent fits and plots. Supply
an accessible GridWorkbench checkout:

```sh
python3 examples/EMT/Distributed/fit.py --gridworkbench /path/to/GridWorkbench --output local/emt-distributed/fits-final
python3 examples/EMT/Distributed/run.py --family frequency --trials 1
python3 examples/EMT/Distributed/run.py --family network --models distributed pi_1 --buses 8 20 --horizon .2 --event load --trials 3
python3 examples/EMT/Distributed/run.py --family network --models distributed pi_1 --buses 8 20 --horizon .2 --event fault --trials 3
python3 examples/EMT/Distributed/run.py --family network --models distributed pi_1 --buses 8 20 --horizon .05 --event load --hybrid --trials 1
python3 examples/EMT/Distributed/run.py --family network --models distributed pi_1 --buses 8 20 --horizon .05 --event fault --hybrid --trials 1
python3 examples/EMT/Distributed/plot.py --results local/emt-distributed/studies
```

Generate the lossless tolerance-refinement curves and a more accurate
frequency-dependent line response with:

```sh
python3 examples/EMT/Distributed/run.py --sine --short --models distributed distributed_capped --rtol 1e-5 --atol 1e-6 --trials 1 --results local/emt-distributed/refinement/coarse
python3 examples/EMT/Distributed/run.py --sine --short --models distributed distributed_capped --rtol 1e-7 --atol 1e-8 --trials 1 --results local/emt-distributed/refinement/fine
python3 examples/EMT/Distributed/run.py --family frequency --models distributed --rtol 1e-7 --atol 1e-9 --trials 1 --results local/emt-distributed/refinement/frequency-fine
python3 examples/EMT/Distributed/plot.py --results local/emt-distributed/studies
```

The plotter reads the nominal lossless runs and the `coarse` and `fine`
subdirectories of `--refinement`, then computes the analytic errors directly
from their waveforms. It writes the measurements to `plots/refinement.json`.
To plot the refined frequency-dependent response separately, pass its results
directory with `--results` and choose a separate `--output` directory.

`--event fault_a` applies a phase-a-to-ground fault near the phase-a voltage
crest; `fault` is a balanced three-phase shunt fault. `baseline` runs only
initial energization. Choose another `--results` directory to retain separate
tolerance refinements. A `STOP` file in the results directory stops a batch
at the next completed-trial boundary; remove it before resuming. `--exe`, `--fits`, `--rtol`, and `--atol` override the
executable, coefficient directory, and solver tolerances. No study sets
`dt_fixed`. `distributed_capped` bounds the adaptive step by the shortest
transport delay; it does not request fixed stepping.

## Cases

The lossless line has characteristic impedance 300 Ω, source resistance
100 Ω, and receiving resistance 600 Ω. Its delay is 300 µs for a 1,000 V
source step at 1 ms, or 20 µs for a 60 Hz sinusoid with 1,000 V peak.
The receiving reflection ratio is −1/6 per round trip. The scripts compare
both the reflection staircase and the analytic steady sinusoidal solution.
The `frequency` family energizes a 60 km fitted line with an equal three-phase
step to excite its zero-sequence mode. The plotter independently inverts the
fitted terminal transfer on a Laplace contour; it uses neither the EMT state
equations nor the delay-history interpolation to compute that reference.

Each network contains a ring with additional chords and two switched shunt
buses, included in the stated bus count. Main buses carry 250 kW resistive
loads at 13.8 kV. The 8-bus case has two governed 10 MVA machines and eight
lines; the 20-bus case has three machines and 23 lines. Line lengths repeat
5, 20, 40, and 60 km. A 500 kW load is connected halfway through the study
and removed at three quarters, or a 10 Ω fault is applied halfway through
and cleared 20 ms later.

The hybrid cases add one converter to the 8-bus case and three to the 20-bus
case. They use 900 Hz PWM, modulation index 0.8, an imposed 28.919 kV DC
voltage, and 1 Ω / 16 mH interface filters. `mu=50000` resolves the existing
smoothed PWM transitions. These are open-loop switching demonstrations with
an ideal DC supply; they do not model a controlled inverter plant or DC-link
energy balance.

Machine and controller parameters come from the existing
[IBR example](../IBR/README.md). The overhead geometry is illustrative 345 kV
geometry used here at 13.8 kV to exercise EMT behavior. These are synthetic
numerical studies, not calibrated utility networks.

Line prehistories and rational bus states start at zero. Machine currents
are initialized from an equal-share estimate of the load; IDA computes
consistent algebraic variables and derivatives. The initial transient is
therefore part of each study, not a solved network steady operating point.

## Fitting

[Geometry and provenance](../../../cases/EMT/Distributed/README.md) identify
the GridWorkbench data. The script computes skin-effect and earth-return
frequency dependence using GridWorkbench, then cyclically transposes the
physical series-impedance and shunt-admittance matrices. Constant real
positive-sequence and zero-sequence projectors give reciprocal modal fits.

`scikit-rf` fits real, stable rational coefficients. Characteristic admittance
uses a constant term. Propagation uses strictly proper rational parts after
removing only the vacuum flight time, `length / 299792458`. Modal dispersion
remains in the rational part; the explicit delay is not the entire 60 Hz phase
shift. Both modal groups can have the same explicit delay.

The fitting band is 0.01 Hz–10 MHz. The script checks the Hermitian part of
the full line's common and differential terminal admittances at 10,003
frequencies from DC to 10 GHz. It applies the smallest sampled attenuation
correction needed to remove negative conductance, rejecting corrections
larger than 0.1%. `summary.json` records errors before correction and the
corrected full propagation error. This sampled check is not an all-frequency
passivity proof.

## Measurements

Each trial records process CPU time reported by the application, elapsed
wall time, peak resident memory, DAE variables, Jacobian nonzeros, accepted
steps, residual evaluations, linear setups, error-test failures, nonlinear
iterations, convergence failures, and the executable SHA-256. The plotter
reports median timings, min/median/95th-percentile/max accepted step sizes,
and the fraction of steps exceeding the shortest and longest delay.

`response.csv` is sampled at the monitor interval (2 or 10 µs for the
lossless cases, 20 µs for networks). `steps.csv` contains actual accepted
internal steps. Monitor spacing is not a solver step-size restriction.
The network phase norm is $\sqrt{(v_a^2+v_b^2+v_c^2)/3}$; it equals
phase RMS voltage for a balanced sinusoid. Plots of fast PWM waveforms are
sampled observations, so they need not
capture every switching-edge extremum.

A π/distributed difference measures a model difference. Use the analytic
lossless reference, tighter tolerances, and a delay-capped distributed run to
assess numerical error separately. Implicit delay overlap allows steps larger
than a delay, but its interpolation error is not independently controlled by
IDA's ordinary error estimator.
