# Switching inverter controls

`GFL.case.json` and `GFM.case.json` share an ideal two-level switching bridge
and LCL filter. The converter voltage feeds a `DependentVoltageSource` with
series resistance and inductance. The capacitor is a Bus shunt, and
`LineLumped` carries grid-side current from the capacitor bus to the terminal
bus. This preserves all physical inductor currents and capacitor voltages.

`Park` transforms the inverter current, capacitor voltage, and grid-side
current. The inverse transform and `Modulation` convert the limited dq
voltage command to the phase modulation inputs of `PWM`. `Angle` integrates
the supplied angular frequency. The GFM case adds `OuterVoltageControl`;
the GFL case supplies current references directly.

Initial state files contain balanced fundamental operating-point estimates,
including PI integral contributions. IDA preserves differential states and
obtains consistent algebraic values and derivatives. Switching ripple develops
from this estimate during startup; the initial point is not a periodic switching
orbit. GFL uses a known grid angle; GFM uses a fixed voltage/frequency primary
reference.

Run and review the plots in the [example directory](../../../examples/EMT/CurrentControl/README.md).
