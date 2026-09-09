# Controller Models

EMT controller models exchange signals with other components. In the switching
inverter examples, PLL supplies angle to the Park operators and frequency to
the controllers through signal ports.

Arrows indicate signal flow; dashed blue paths denote measurements and switching signals.
The converter current $\mathbf{i}$ is positive
out of the bridge.
The power-invariant transforms are implicit at the $abc$/$dq0$ boundary;
the implemented controllers use the $d$ and $q$ components. PWM applies the
inverse transform internally.

Symbol | Producer | Consumer | Coordinates
------ | -------- | -------- | -----------
$\mathbf{e}$ | `Converter.e` | `Filter.e` | $abc$
$\mathbf{v}_g$ | Terminal Bus | Grid-voltage Park input | $abc$
$\mathbf{v}_{\mathrm{o}}$ | `Filter.vo` | Voltage Park input | $abc$
$\mathbf{i}$ | `Filter.i` | Current Park input | $abc$
$\mathbf{i}_g$ | `Filter.ig` | Terminal Bus and grid-current Park input | $abc$
$v_{\mathrm{dc}}$ | External constant | Converter and PWM | Scalar
$\mathbf{u}^{\mathrm{cmd}}$ | `InnerCurrentControl.u` | `PWM.u` | $dq$
$\mathbf{u}^{\mathrm{lim}}$ | `PWM.ulim` | `InnerCurrentControl.ulim` | $dq$

The current and voltage controllers retain their local voltage-input name
`v`: it receives the $dq$ components of transformed $\mathbf{v}_{\mathrm{o}}$.
`Filter.v` is the separate $abc$ voltage input from the terminal Bus.
The filter is drawn as a Z–C–Z circuit, with series impedances
$\mathbf Z_{\mathrm{s}}=\mathbf R_{\mathrm{s}}+s\mathbf L_{\mathrm{s}}$ and
$\mathbf Z_g=\mathbf R_g+s\mathbf L_g$, and shunt capacitance $\mathbf C$.

## Grid Following

![Grid-following inverter wiring](../../../../../docs/Figures/EMT/Controller/diagram_gfl.png)

### REGCA correspondence

[REGCA](../../../PhasorDynamics/Converter/REGCA/README.md) is a controlled
current source. Its two low-voltage functions have no block here.

REGCA | Grid following
----- | --------------
LVPL, the $I_L(V_M)$ ceiling on active current | Not modeled. The only current bound is the fixed [InnerCurrentControl](InnerCurrentControl/README.md) limit $I^{\max}$ on $\|\mathbf{i}^{\mathrm{cmd}}\|_2$.
LVACM, the $\text{linseg}(V_T;V_{A0},V_{A1},1)$ factor on injected active current | Not modeled. REGCA scales the injection because it has no converter or current loop. Here the injected current follows the current loop, PWM, Converter, and Filter.

## GFM Voltage Control

![GFM control schematic](../../../../../docs/Figures/EMT/Controller/diagram_gfm.png)

The $\mathbf v_g$ input to Voltage Control is proposed in this general schematic.
The implemented OuterVoltageControl uses $\mathbf v_{\mathrm{o}}$ and
$\mathbf i_g$; its measurement interface is unchanged.

## Models

- [InnerCurrentControl](InnerCurrentControl/README.md): converter current control in $dq$ coordinates.
- [REECB](REECB/README.md): terminal power, voltage, and current-priority control.
- [REPCA](REPCA/README.md): renewable plant voltage, reactive-power, and frequency control.
- [OuterPowerControl](OuterPowerControl/README.md): active and reactive power control in $dq$ coordinates.
- [OuterVoltageControl](OuterVoltageControl/README.md): filter-capacitor voltage control in $dq$ coordinates.
- [IEEET1](IEEET1/README.md)
- [PWM](PWM/README.md)
- [TGOV1](TGOV1/README.md)
- [SEXS-PTI](SEXS-PTI/README.md): simplified excitation system.
- [GASTPTI](GASTPTI/README.md): gas turbine governor with exhaust-temperature limiting.
- [IEEEST](IEEEST/README.md): power system stabilizer.
