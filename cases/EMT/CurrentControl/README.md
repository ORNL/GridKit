# Switching inverter controls

`GFL.case.json` and `GFM.case.json` share an ideal two-level switching
bridge, constant DC voltage, and physical LCL filter connected to a stiff grid.
The converter voltage feeds a `Filter` component, which owns both inductor
currents and the capacitor voltage and injects grid-side current into the
terminal Bus. Its `i`, `vo`, and `ig` signal outputs supply the controllers
and reference-frame operators.

Both cases connect `PLL` to the terminal Bus voltage outputs. Its `theta` output supplies every
`Park` transform, and its `omega` output supplies the controller frequency
inputs. `PWM` limits the dq voltage command to the available DC voltage,
returns the limited command to `InnerCurrentControl`, and applies the inverse
transform before carrier comparison.

`GFL` uses `OuterPowerControl` with grid-current targets derived as `Pref/V`
and `-Qref/V`. `GFM` uses `OuterVoltageControl` to regulate the
capacitor voltage, with a magnitude reference step from 208 V to 209 V at 0.04 s
and back at 0.12 s. It is a PLL-synchronized, grid-connected voltage-control
study. Both outer loops supply `icmd` to `InnerCurrentControl` and receive
its `ilim` output for anti-windup.

A constant 400 V signal supplies `vdc` to PWM and Converter. The GFM
reference step scales both components of the initial capacitor-voltage
reference in the terminal-voltage frame.

Initial states are balanced fundamental operating-point estimates. The PLL
infers its initial angle from terminal Bus voltage; controller outputs determine
the PI integral contributions. IDA preserves differential states and obtains
consistent algebraic values and derivatives. Switching ripple develops during
startup. In `GFM`, the stiff-grid magnitude and phase match the
initialized terminal phasor, preserving the initial filter operating point.

Run and review the plots in the [example directory](../../../examples/EMT/CurrentControl/README.md).
