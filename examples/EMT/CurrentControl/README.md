# Inverter current and voltage control

`GFL` tracks a current reference against a stiff grid with a known angle.
`GFM` supplies an islanded load through cascaded voltage and current control.
Both use the same continuous PWM, converter, and physical LCL filter.

The bridge has a 20 mF DC link initialized at 400 V and 6 kHz centered PWM.
Constant source current matches the initial fundamental power balance. The
filter uses 2 mH / 0.2 Ω, 100 µF, and 1 mH / 0.1 Ω. Nominal voltage is 208 V
line-to-line RMS at 60 Hz; dq quantities use the power-invariant Park transform.
Current and voltage loop bandwidths are nominally 400 Hz and 60 Hz.

The current reference steps through 8, 16, 45, and 8 A, with a 30 A limit.
The voltage study connects and disconnects a second 20 Ω/phase load.
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

Compare resolved switching and broad smoothing over the same 0.1 s window:

```bash
python3 examples/EMT/CurrentControl/run.py --tmax 0.1
python3 examples/EMT/CurrentControl/run.py --mu 240 --tmax 0.1 --dt-monitor 1e-5 --output simulation-smooth
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
