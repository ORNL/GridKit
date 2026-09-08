# Synthetic 10-bus EMT grid

`TenBus.case.json` contains exactly ten buses, three synchronous machines,
three PWM/Converter assemblies, nine three-phase lines, seven resistive loads,
and three ideal switches. `TenBus.state.json` supplies a fundamental-frequency
operating-point estimate. This is a synthetic demonstration, not a utility case.

All buses use 13.8 kV line-to-line RMS and 60 Hz nominal values. Electrical
network parameters are SI; machine winding parameters and internal variables
follow the Machine model's per-unit convention. There is no global power base.

| Buses | Equipment |
|---|---|
| 1, 2, 3 | 10 MVA machines with TGOV1 governors; inertia constants 3.7, 4.5, 3.0 s |
| 4, 5, 6 | PWM → Converter → dependent voltage source with 1 Ω + 16 mH series filter; 0.5 MW local load each |
| 7 | 6 MW nominal load and fault-switch connection |
| 8 | 5 MW nominal load and switched-load connection |
| 9 | Initially disconnected 2 MW nominal load |
| 10 | Initially disconnected 2 Ω/phase fault resistor |

Load powers are three-phase values at nominal voltage, not constant-power
loads. The fault is a finite-resistance three-phase shunt, not a bolted fault.

Lines connect 1–7, 2–7, 3–8, 4–7, 5–8, 6–8, 1–2, 2–3, and 7–8.
The `tie` switch initially joins buses 7 and 8 directly, in parallel with the
7–8 line. Opening it transfers current to the finite-impedance paths; it does
not island the network. `load_step` joins buses 8 and 9. `fault` joins buses 7
and 10. Open load/fault buses are anchored to ground by their resistors.

The 1–2 and 2–3 lines each have 1 µF total shunt capacitance, split equally
between their ends. These capacitors provide the machine-terminal voltage
states needed during initialization. Other line shunts are zero. Series
matrices are diagonal, with no mutual coupling. `dx=1 m` makes the per-length
matrices numerically equal to the total line parameters.

## Converter scope and smoothing

Each PWM uses M=0.8, fm=60 Hz, fc=900 Hz, and centered pulses. The default
shared `GridKit::Math::MU<double>` is 240 s⁻¹. Continuous PWM preserves the
fixed-duty carrier-period mean. `build_case.py` therefore sets
`Vdc = 2 sqrt(2/3) (1.02 × 13800) / M = 28732.515 V` for a 1.02 pu
mean open-circuit line-to-line voltage. All studies use this same DC voltage;
changing mu requires no DC compensation.

At `mu=240`, the 18.3 ms logistic edge width exceeds the 1.11 ms carrier
period; this setting suppresses switching. At `mu=50000`, the width is
87.9 µs and the carrier switching is resolved.

The converters have fixed-frequency modulation and ideal unlimited DC sources.
There is no PLL, current controller, current limit, DC energy storage, exciter,
or converter protection. Fault current and ride-through traces describe these
model equations; they are not hardware capability predictions. Machine field
voltage is held at its initialized value; governors regulate mechanical power.

## Initialization and reproduction

Run `python3 cases/EMT/IBR/build_case.py` from the repository root (requires
NumPy) to regenerate the case and state files. The script solves a linear
fundamental-frequency network with specified machine terminal voltages and
converter source voltages, then derives machine P/Q, instantaneous bus
voltages, line currents and filter currents. It does not synthesize any
simulation results.

IDA computes consistent algebraic variables and derivatives while retaining
these differential initial states. This is a fundamental-frequency estimate;
switching ripple and model initialization can still produce startup transients.
Full studies apply disturbances at 1 s. Compare each event run against the
undisturbed baseline.

A sixth solver study uses runtime `mu=50000` to resolve carrier switching. The shared
scale affects machine saturation and governor limiters as well as PWM.

See [the six scenarios and plotting workflow](../../../examples/EMT/IBR/README.md).

## REGFMA variant

`REGFMA.case.json` replaces the three PWM/Converter/filter assemblies with
5 MVA REGFMA (REGFM_A1) sources at buses 4–6. The network, machines, loads,
switches, and nominal voltage are retained. Each source has physical RL coupling
and a damped capacitor shunt. `REGFMA.state.json` preserves the initial network
bus voltages and net injections; source currents also supply the initialized
capacitor shunts. Each source derives its internal voltage and unattached
references from that operating-point estimate.

Run `python3 cases/EMT/IBR/build_regfma_case.py` to regenerate these two files
from the existing ten-bus case and state. See the
[REGFMA fault study](../../../examples/EMT/IBR/REGFMA/README.md) for source
parameters, simulation commands, plots, and numerical checks.
