# Frequency dependent EMT network

`Hybrid.case.json` connects a synchronous Machine with IEEET1 and a
PWM → Converter → DependentVoltageSource assembly to a sinusoidal source
through a lumped line. The source and converter filter admittances, line
series impedance and shunt admittances, and load impedance contain real
poles. All electrical coefficients use SI units; `dx=1 m` makes the line's
per-length coefficients equal to its total coefficients. The two line
terminals use independent states for the shared `Yp` coefficient set.

The network uses the existing LineLumped model. Propagation, Delay, and
LineDistributed are not required. The Machine retains its nonlinear winding
equations and receives the exciter's per-unit field voltage.

`Circuit.case.json` contains the source, line, and load alone. Each phase has
the following passive realizations, where `||` denotes parallel connection:

| Element | Transfer function |
| --- | --- |
| Source admittance | `1/(10 + 0.01s) + 1/(160 + 0.04s)` |
| Line series impedance | `1 + 0.004s + (3 || 0.0015s)` |
| Total line shunt admittance | `1e-5 + 2e-6s + 1/(20000 + (40/3)s)` |
| Load impedance | `200 + 0.03s + (30 || 0.06s)` |
| Converter filter admittance | `1/(8 + 0.04s) + 1/(400 + 0.2s)` |

The shunt admittance is split equally between the line terminals. These are
synthetic coefficients derived from passive circuit elements, not measured
fits. The source and load impedance values therefore vary with frequency
beyond a constant R/L approximation.

Run the validation from the repository root, using an Enzyme and SUNDIALS/KLU
build:

```bash
python3 cases/EMT/FrequencyDependent/validate.py \
  --exe build/emt-rational/application/EMT/EMTDynamicSimulation
```

The script checks settled voltage and current waveforms at 30, 60, and
180 Hz against an independent complex circuit solve. It also runs the hybrid
case twice with a factor of ten tighter solver tolerances, checks every
monitored waveform, and verifies Kirchhoff's current law at both buses.
It requires only Python's standard library and uses temporary output files.
Pass `--results /tmp/emt-frequency-results` to retain logs and measured errors.
The same checks are registered as `EMTFrequencyDependent` in CTest.

To run the hybrid study directly, pass `Hybrid.solver.json` to
`EMTDynamicSimulation`. It records 0.2 s at 100 µs intervals. Initial bus
voltages and machine dispatch are estimates; the simulation includes the
network energization transient. PWM runs at 900 Hz with `mu=50000` and an
ideal 28.919 kV DC supply, selected for a 1.02 pu fundamental open-circuit
voltage. The sinusoidal source and initial voltages share its phase. This
case exercises the switching and field-control
equations; the converter has fixed modulation and no current or DC-link
controller.
