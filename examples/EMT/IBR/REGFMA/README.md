# Ten-bus EMT REGFMA fault response

This study replaces the three PWM/Converter/filter assemblies in the
[ten-bus case](../../../../cases/EMT/IBR/README.md) with three 5 MVA
[REGFMA (REGFM_A1)](../../../../GridKit/Model/EMT/Component/Source/REGFMA/README.md)
sources at buses 4, 5, and 6. The ten network buses, three governed synchronous
machines, lines, loads, switches, and initial net network injections are retained.
Each source uses 13.8 kV, `XL=0.15`, `RL=0.03`, `mp=0.01`, `mq=0.05`, `ImaxF=2`,
`VFlag=true`, and `QVFlag=true`; other parameters use the model defaults.
Unconnected power and voltage references are initialized from the operating point.

The physical filter uses the source's single RL coupling branch and a series RC
shunt at its bus. On the 5 MVA, 13.8 kV base, its per-phase values are

Quantity | Value
-------- | -----
Series resistance | 1.14264 Ω
Series inductance | 15.1547 mH
Shunt capacitance | 4.64290 µF
Capacitor damping resistance | 0.114264 Ω
Undamped LC resonance | 600 Hz

These filter parameters follow the main-circuit defaults in the
[PNNL PSCAD reference](https://github.com/pnnl/PSCAD-and-PSSE-Version-of-WECC-Grid-Forming-Inverter-Models/releases/tag/V1).
Its separate virtual-admittance and inner-current controls are not reproduced
by this averaged EMT realization. No additional series inductance is added to
`XL`. Three capacitor buses and three resistive `LineLumped` branches realize
the shunts with existing components. The capacitor voltages and currents are
initialized at 60 Hz; source currents include the shunt currents so the net
network injections match the original operating point.

`FaultClearing.solver.json` applies the existing three-phase 2 Ω/phase shunt
at 1.00 s and clears it at 1.06 s. The simulation spans 0–3 s with adaptive
IDA stepping, sparse KLU, `mu=240`, relative tolerance 1e-7, and absolute
tolerance 1e-6. Monitoring uses 100 µs spacing, with both pre-event and
post-event samples retained. Monitor spacing is not the accepted solver step.

## Run and review

Build `EMTDynamicSimulation` with Enzyme and SUNDIALS KLU enabled, then run
from the repository root:

```bash
python3 examples/EMT/IBR/REGFMA/run.py --exe build/application/EMT/EMTDynamicSimulation
python3 examples/EMT/IBR/REGFMA/plot.py
```

The runner requires Python's standard library; plotting also requires NumPy
and Matplotlib. To regenerate the case and initial state from the existing
ten-bus files, run `python3 cases/EMT/IBR/build_regfma_case.py`.

Open [the local gallery](results/index.html) to inspect the
[full response](results/plots/response.png) and
[fault detail](results/plots/fault.png). Each PNG has a companion vector PDF.
Generated files remain in ignored `results/`, alongside input snapshots,
`run.log`, runtime and input hashes in `run.json`, monitor data in
`FaultClearing.csv`, accepted solver steps in `steps.csv`, and measured extrema
and numerical checks in `summary.json`. The study records selected monitors
without a complete DAE state dump.

## Recorded data and checks

The response figure shows terminal voltage, REGFMA internal frequency and
machine rotor frequency, filtered active and reactive power, terminal current,
and internal voltage command. The detail figure shows instantaneous power and
the recorded three-phase voltage and current at bus 4 around the fault.
The response figure's voltage inset retains the full event peaks.
Shading marks the fault interval. Positive source power denotes injection.
Buses 4 and 5 have identical source and feeder parameters; their responses
coincide with the 7–8 tie closed.

Voltage and current magnitudes use the power-invariant Clarke transform.
For balanced sinusoidal waveforms, the voltage magnitude equals line-to-line
RMS voltage; its base is 13.8 kV. Current magnitude uses the source base
`S/V`; 1 pu corresponds to `S/(sqrt(3) V)` amperes of phase RMS current.
The plotted magnitudes are instantaneous αβ norms, without cycle averaging.
Filtered power uses the model's `pf` and `qf` states; terminal power in the
detail figure is reconstructed directly from phase quantities. REGFMA
frequency is `omega/(2 pi)`; machine rotor frequency is `60 omega`.

The plotting script verifies input hashes, finite data, the monitor timeline,
pre/post-event switch commands, and accepted-step counts. It independently
reconstructs phase current balance at buses 4–6 and the capacitor buses, and
each REGFMA's terminal active and reactive power. It checks the damping-resistor
voltage drops, zero-sequence current, inductor-current continuity across events,
and the radial current-reference limit reconstructed from the applied source
voltage. Actual current peaks are recorded separately; physical inductor current
can overshoot the reference limit. These are numerical consistency checks,
not external validation.

REGFMA uses WECC droop and limiting controls with a physical RL branch in
balanced EMT coordinates. The surrounding network retains its line dynamics.
This example demonstrates that averaged model's fault response; it does not
represent PWM switching or hardware protection.
