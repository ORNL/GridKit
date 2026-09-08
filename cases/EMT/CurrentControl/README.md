# Switching inverter controls

`GFL.case.json` and `GFM.case.json` share an ideal two-level switching bridge,
dynamic DC link, and LCL filter. The converter voltage feeds a
`DependentVoltageSource` with series resistance and inductance. The capacitor is a Bus shunt, and
`LineLumped` carries grid-side current from the capacitor bus to the terminal
bus. This preserves all physical inductor currents and capacitor voltages.

`Park` transforms the inverter current, capacitor voltage, and grid-side
current. The inverse transform and `Modulation` convert the limited dq
voltage command to the phase modulation inputs of `PWM`. GFM uses `Angle`
and `OuterVoltageControl`. GFL uses `PLL` at the capacitor terminal and
`OuterPowerControl` with current targets derived as `Pref/V` and `-Qref/V`.
The power-reference parameters use rated voltage and the initial measured dq
current to preserve the operating point. The outer loop reads
Park-transformed grid-side current, supplies `icmd` to
InnerCurrentControl, and receives its `ilim` output for anti-windup.

`DCLink` uses 20 mF with an initial voltage of 400 V. Constant source current
matches the initial fundamental power balance; subsequent power imbalance
changes the DC voltage. Its voltage feeds current control, Modulation, and
Converter; the converter DC current returns to the capacitor.

Initial state files contain balanced fundamental operating-point estimates,
with controller output values that determine the PI integral contributions.
IDA preserves differential states and obtains consistent algebraic values and
derivatives. Switching ripple develops
from this estimate during startup; the initial point is not a periodic switching
orbit. GFL starts from 1664 W and 0 var exported at the capacitor terminal;
its PLL infers the angle from the initialized bus voltage. GFM uses a fixed
voltage/frequency primary reference.

Run and review the plots in the [example directory](../../../examples/EMT/CurrentControl/README.md).
