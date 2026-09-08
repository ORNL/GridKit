# Switching inverter controls

`GFL.case.json` and `GFM.case.json` share an ideal two-level switching
bridge, dynamic DC link, and physical LCL filter connected to a stiff grid.
The converter voltage feeds a `DependentVoltageSource` with series resistance
and inductance. The capacitor is a Bus shunt, and `LineLumped` carries grid-side
current from the capacitor bus to the terminal bus.

Both cases use `PLL` at the capacitor bus. Its `theta` output supplies every
`Park` transform, and its `omega` output supplies the controller frequency
inputs. The inverse transform and `Modulation` convert the limited dq voltage
command to the phase modulation inputs of `PWM`.

`GFL` uses `OuterPowerControl` with grid-current targets derived as `Pref/V`
and `-Qref/V`. `GFM` uses `OuterVoltageControl` to regulate the
capacitor voltage, with a d-axis reference step from 208 V to 209 V at 0.04 s
and back at 0.12 s. It is a PLL-synchronized, grid-connected voltage-control
study. Both outer loops supply `icmd` to `InnerCurrentControl` and receive
its `ilim` output for anti-windup.

`DCLink` uses 20 mF with an initial voltage of 400 V. Constant source current
matches the initial fundamental power balance; subsequent power imbalance
changes the DC voltage. Its voltage feeds current control, Modulation, and
Converter; the converter DC current returns to the capacitor.

Initial states are balanced fundamental operating-point estimates. The PLL
infers its initial angle from capacitor voltage; controller outputs determine
the PI integral contributions. IDA preserves differential states and obtains
consistent algebraic values and derivatives. Switching ripple develops during
startup. In `GFM`, the stiff-grid magnitude and phase match the
initialized terminal phasor, preserving the initial filter operating point.

Run and review the plots in the [example directory](../../../examples/EMT/CurrentControl/README.md).
