# Inverter current and voltage control

`GFL` tracks parameter-derived dq grid-current references against a stiff grid using
a terminal-voltage PLL and cascaded grid-side and inverter-side current control.
`GFM` regulates capacitor voltage against a stiff grid through
cascaded voltage and current control. Both cases obtain frame angle and
frequency from the capacitor-voltage PLL through signal ports.
Both use the same continuous PWM, converter, and physical LCL filter.

The bridge has a 20 mF DC link initialized at 400 V and 6 kHz centered PWM.
Constant source current matches the initial fundamental power balance. The
filter uses 2 mH / 0.2 Ω, 100 µF, and 1 mH / 0.1 Ω. Nominal voltage is 208 V
line-to-line RMS at 60 Hz; dq quantities use the power-invariant Park transform.
The current loop has a nominal 400 Hz bandwidth. The grid-connected voltage
loop uses `Kp=0.1507964474` S, `Ki=2.131834551` S/s, and `Kaw=376.9911184` s⁻¹.

The current targets are fixed derived parameters `Pref/V` and `-Qref/V`,
chosen to match the initial measured grid current. The q-axis target is zero
up to roundoff, and the inverter-current limit is 30 A. The PLL gains are
80 rad/s and 2500 rad/s²; the outer loop uses `Kp=0.01`, `Ki=40`, and `Kaw=200`.
The voltage study steps the d-axis capacitor-voltage reference from 208 V
to 209 V at 0.04 s and restores 208 V at 0.12 s. Its grid voltage matches
the initial terminal phasor (206.93 V line-to-line RMS).
Initial states are fundamental operating-point estimates; switching ripple
develops during startup. See the [cases](../../../cases/EMT/CurrentControl/README.md)
for connections and initialization.

## Run and plot

Defaults are 0.2 s, 2 µs monitoring, and `mu=1000000`.

```bash
cmake --build build --target EMTDynamicSimulation -j 10
python3 examples/EMT/CurrentControl/run.py
python3 examples/EMT/CurrentControl/plot.py
```

Compare resolved switching and broad smoothing over the same 0.2 s window,
including the voltage-reference step and return:

```bash
python3 examples/EMT/CurrentControl/run.py --tmax 0.2
python3 examples/EMT/CurrentControl/run.py --mu 240 --tmax 0.2 --dt-monitor 1e-5 --output simulation-smooth
python3 examples/EMT/CurrentControl/plot.py --compare simulation-smooth
```

`mu` changes switching resolution while preserving the fixed-duty
carrier-period mean. The logistic 10–90% width is 4.39 µs at `1000000`
and 18.3 ms at `240`. The circuit and controller parameters are unchanged.

`--scenario GFL` or `--scenario GFM` selects one study. `--exe` selects the
executable; output and comparison paths are relative to this example.
The plotter uses the common available interval, or its `--tmax` option,
and gives corresponding panels identical x- and y-axis limits.

Each run retains its effective inputs, hashes, waveforms, log, and IDA statistics.
The plotter writes PNG/PDF figures, compressed waveforms, and `summary.json`
with tracking, bridge power balance, DC-link energy balance, and continuous
PWM comparisons. Generated inputs and data are ignored; PNG/PDF figures remain available for review.
CSV waveforms may be removed after plotting; the plotter also reads NPZ files.

Plot | Resolved switching | Broad smoothing
---- | ------------------ | ---------------
Current tracking | [GFL](simulation/GFL.png) | [GFL](simulation-smooth/GFL.png)
Voltage control | [GFM](simulation/GFM.png) | [GFM](simulation-smooth/GFM.png)
PWM and bridge | [Switching](simulation/switching.png) | [Switching](simulation-smooth/switching.png)
Fourier amplitudes | [Harmonics](simulation/harmonics.png) | [Harmonics](simulation-smooth/harmonics.png)

## Regression coverage

The continuous-PWM unit fixture checks duty mean, the smoothing-to-switching
transition, and input values and gradients at the same evaluation time. These
examples report closed-loop behavior without prescribing exact transients or solver steps.

`validate.py` runs the GFL model against the analytic high-voltage LCL solution
for prescribed dq grid current. It checks broad smoothing and resolved
switching, and is registered as `EMTGflPowerControl`.
`--scenario GFM` checks phase-domain voltage tracking, PLL
alignment and frequency, and current limiting before, during, and after the
reference step, in both PWM resolutions. It is registered as
`EMTPLLVoltageControl`. Both validators use only the Python standard library.

```bash
python3 examples/EMT/CurrentControl/validate.py --exe build/application/EMT/EMTDynamicSimulation
python3 examples/EMT/CurrentControl/validate.py --exe build/application/EMT/EMTDynamicSimulation --scenario GFM
```

The validator uses the balanced LCL solution for the supplied
dq current targets and reconstructs P and Q independently from phase samples
and Park measurements. The reactive-power bound is 0.002 var: repeated resolved-switching runs
differed from the balanced analytic solution by 0.001194 var. The current
reference fixes dq current; it does not impose exact average reactive power
in the presence of switching ripple.
