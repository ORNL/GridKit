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
shared `GridKit::Math::MU<double>` is 240 s⁻¹. For the documented logistic
smoothing, the fundamental attenuation is

`(π ω / μ) / sinh(π ω / μ) = 0.07098472`, with `ω=2π·60`.

`build_case.py` computes the exact fundamental Fourier coefficient of the
sampled ideal pulse train, applies this attenuation, and selects an effective
constant DC value of about 407.357 kV to give a 1.02 pu open-circuit AC
fundamental. This is explicit mathematical compensation for broad smoothing,
not a realistic DC-link rating or an implicit transformer. The actual gate
signals are correspondingly smooth and the carrier ripple is strongly suppressed.
Changing shared MU requires regenerating and rechecking this case.

The converters have fixed-frequency modulation and ideal unlimited DC sources.
There is no PLL, current controller, current limit, DC energy storage, exciter,
or converter protection. Fault current and ride-through traces describe these
model equations; they are not hardware capability predictions. Machine field
voltage is held at its initialized value; governors regulate mechanical power.

## Initialization and reproduction

Run `python3 cases/EMT/IBR/build_case.py` from the repository root (requires
NumPy) to regenerate the case and state files. The script solves a linear
fundamental-frequency network with specified machine terminal voltages and
converter source voltages, then derives machine P/Q and instantaneous bus
voltages. It does not synthesize any simulation results.

Line and filter currents still start at their model defaults. IDA computes
consistent algebraic variables and derivatives while retaining differential
states, so an energization transient remains. All examples retain it and apply
disturbances at 1 s. Compare each event run against the undisturbed baseline.

A sixth solver study uses runtime `mu=50000` and explicit DC overrides to
resolve carrier switching while preserving the AC fundamental. The shared
scale affects machine saturation and governor limiters as well as PWM.

See [the six scenarios and plotting workflow](../../../examples/EMT/IBR/README.md).
