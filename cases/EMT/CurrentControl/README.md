# Switching inverter controls

`GFL.case.json` uses an ideal two-level switching
bridge, constant DC voltage, and physical LCL filter connected to a stiff grid.
The converter voltage feeds a `Filter` component, which owns both inductor
currents and the capacitor voltage and injects grid-side current into the
terminal Bus. Its `i`, `vo`, and `ig` signal outputs supply the controllers
and reference-frame operators.

The case connects `PLL` to the terminal Bus voltage outputs. Its `theta` output supplies every
`Park` transform, and its `omega` output supplies the controller frequency
inputs. `PWM` limits the dq voltage command to the available DC voltage,
returns the limited command to `InnerCurrentControl`, and applies the inverse
transform before carrier comparison.

`GFL` uses `OuterPowerControl` to regulate terminal active and reactive power
against external constant `Pref` and `Qref` signals.
A separate Park transform supplies terminal Bus voltage; grid-side current
comes from Filter `ig`. The outer loop supplies `icmd` to `InnerCurrentControl`
and receives its `ilim` output for anti-windup.

A constant 400 V signal supplies `vdc` to PWM and Converter.

The state file prescribes frequency, terminal Bus voltage, and Filter grid
current. Initialization derives the remaining Filter outputs and controller
commands through their signal connections, including the PI integral contributions. IDA preserves differential states and obtains
consistent algebraic values and derivatives. Switching ripple develops during
startup.

Run and review the plots in the [example directory](../../../examples/EMT/CurrentControl/README.md).
