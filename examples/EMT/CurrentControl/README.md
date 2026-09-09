# Inverter power and current control

`GFL` regulates terminal active and reactive power against a stiff grid using
a terminal-voltage PLL and an inner inverter-side current controller.
Frame angle and frequency come from the terminal-voltage PLL through signal ports.
The case uses continuous PWM, a converter, and a physical LCL filter.
The [Filter](../../../GridKit/Model/EMT/Component/Filter/README.md) component
connects the converter to the terminal Bus and supplies converter-side current,
capacitor voltage, and grid-side current through its `i`, `vo`, and `ig` outputs.

A constant 400 V signal feeds PWM and Converter, with 6 kHz centered PWM.
The filter uses 2 mH / 0.2 Ω, 100 µF, and 1 mH / 0.1 Ω. Nominal voltage is 208 V
line-to-line RMS at 60 Hz; dq quantities use the power-invariant Park transform.
The current loop has a nominal 400 Hz bandwidth.

The power loop measures terminal Bus voltage and Filter `ig` in the same
PLL frame. Its setpoints match the initial terminal P/Q, and its power errors
are normalized by the rated voltage. The inverter-current limit is 30 A.
The PLL gains are
80 rad/s and 2500 rad/s²; the outer loop uses `Kp=0.01`, `Ki=40`, and `Kaw=200`.
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

Compare resolved switching and broad smoothing over the same 0.2 s window:

```bash
python3 examples/EMT/CurrentControl/run.py --tmax 0.2
python3 examples/EMT/CurrentControl/run.py --mu 240 --tmax 0.2 --dt-monitor 1e-5 --output simulation-smooth
python3 examples/EMT/CurrentControl/plot.py --compare simulation-smooth
```

`mu` changes switching resolution while preserving the fixed-duty
carrier-period mean. The logistic 10–90% width is 4.39 µs at `1000000`
and 18.3 ms at `240`. The circuit and controller parameters are unchanged.

`--exe` selects the executable; output and comparison paths are relative to
this example.
The plotter uses the common available interval, or its `--tmax` option,
and gives corresponding panels identical x- and y-axis limits.
The PWM and bridge figure shows the final six carrier periods so individual
switching pulses are visible.

Each run retains its effective inputs, hashes, waveforms, log, and IDA statistics.
The plotter writes PNG/PDF figures, compressed waveforms, and `summary.json`
with tracking, bridge voltage projection, and continuous PWM comparisons. Generated inputs and data are ignored; PNG/PDF figures remain available for review.
CSV waveforms may be removed after plotting; the plotter also reads NPZ files.

Plot | Resolved switching | Broad smoothing
---- | ------------------ | ---------------
Current tracking | [GFL](simulation/GFL.png) | [GFL](simulation-smooth/GFL.png)
PWM and bridge | [Switching](simulation/switching.png) | [Switching](simulation-smooth/switching.png)
Fourier amplitudes | [Harmonics](simulation/harmonics.png) | [Harmonics](simulation-smooth/harmonics.png)

## Regression coverage

The continuous-PWM unit fixture checks duty mean, the smoothing-to-switching
transition, and input values and gradients at the same evaluation time. This
example reports closed-loop behavior without prescribing exact transients or solver steps.

`validate.py` runs the GFL model against the balanced LCL solution for prescribed
terminal P/Q. It checks broad smoothing and resolved switching, uses only the
Python standard library, and is registered as `EMTGflPowerControl`.

```bash
python3 examples/EMT/CurrentControl/validate.py --exe build/application/EMT/EMTDynamicSimulation
```

The validator uses the balanced LCL solution for the specified terminal P/Q
and reconstructs power independently from phase samples and Park measurements.
The reactive-power bound is 0.002 var.
