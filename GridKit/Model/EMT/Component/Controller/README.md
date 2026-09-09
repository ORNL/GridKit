# Controller Models

EMT controller models exchange signals with other components. Reference-frame
angle and frequency are supplied through signal ports.

Arrows indicate signal flow; dashed blue paths denote measurements and switching signals.
Measured currents are positive from the converter toward the grid.
Electrical-side voltages are phase-to-neutral vectors.
The power-invariant transforms are implicit at the $abc$/$dq0$ boundary;
the implemented controllers use the $d$ and $q$ components. PWM applies the
inverse transform internally.

Symbol | Producer | Consumer | Coordinates
------ | -------- | -------- | -----------
$\mathbf{e}$ | `Converter.e` | `Filter.e` | $abc$
$\mathbf{v}_{\mathrm{t}}$ | Inverter terminal Bus | `Filter.v` and terminal-voltage Park input | $abc$
$\mathbf{v}_{\mathrm{o}}$ | `Filter.vo` | Voltage Park input | $abc$
$\mathbf{i}$ | `Filter.i` | Current Park input | $abc$
$\mathbf{i}_g$ | `Filter.ig` | Terminal Bus and grid-current Park input | $abc$
$\mathbf{v}_{\mathrm{PCC}}$ | Grid-side PCC Bus (GFM schematic) | Plant Control voltage measurement | $abc$
$\mathbf{i}_{\mathrm{PCC}}$ | PCC disconnect current (GFM schematic) | Plant Control current measurement | $abc$
$v_{\mathrm{dc}}$ | External constant | Converter and PWM | Scalar
$\mathbf{u}^{\mathrm{cmd}}$ | `InnerCurrentControl.u` | `PWM.u` | $dq$
$\mathbf{u}^{\mathrm{lim}}$ | `PWM.ulim` | `InnerCurrentControl.ulim` | $dq$

The current and voltage controllers retain their local voltage-input name
`v`: it receives the $dq$ components of transformed $\mathbf{v}_{\mathrm{o}}$.
`Filter.v` is the separate $abc$ terminal voltage $\mathbf v_{\mathrm{t}}$.
The LCL filter has series impedances
$\mathbf Z_{\mathrm{s}}=\mathbf R_{\mathrm{s}}+s\mathbf L_{\mathrm{s}}$ and
$\mathbf Z_g=\mathbf R_g+s\mathbf L_g$, and shunt capacitance $\mathbf C$.
The GFM schematic shows the series resistors and inductors explicitly.

## Grid Following

![Grid-following inverter wiring](../../../../../docs/Figures/EMT/Controller/diagram_gfl.png)

### REGCA correspondence

[REGCA](../../../PhasorDynamics/Converter/REGCA/README.md) is a controlled
current source. Its two low-voltage functions have no block here.

REGCA | Grid following
----- | --------------
LVPL, the $I_L(V_M)$ ceiling on active current | Not modeled. The only current bound is the fixed [InnerCurrentControl](InnerCurrentControl/README.md) limit $I^{\max}$ on $\|\mathbf{i}^{\mathrm{cmd}}\|_2$.
LVACM, the $\text{linseg}(V_T;V_{A0},V_{A1},1)$ factor on injected active current | Not modeled. REGCA scales the injection because it has no converter or current loop. Here the injected current follows the current loop, PWM, Converter, and Filter.

## Grid Forming

![GFM control schematic](../../../../../docs/Figures/EMT/Controller/diagram_gfm.png)

PCC labels the grid-side connection point, not the disconnect. This assumes
the connection is the Local EPS–Area EPS boundary in IEEE 1547 terminology;
an individual inverter terminal is not necessarily the PCC.[^pcc]
The subscripts identify physical locations; they are not mandated by the standard.
For transmission-connected plants, IEEE 2800 distinguishes the unit point of
connection (POC), plant point of measurement (POM), and point of interconnection
(POI); these need not coincide.[^ieee2800]

$\mathbf i$ is the converter-side filter current, $\mathbf v_{\mathrm o}$ the
capacitor voltage, $\mathbf i_g$ the grid-side filter current, and
$\mathbf v_{\mathrm{t}}$ the inverter terminal voltage before the disconnect.
For the unbranched circuit shown, $\mathbf i_g=\mathbf i_{\mathrm{PCC}}$ and
$\mathbf i-\mathbf i_g=\mathbf C\,\mathrm d\mathbf v_{\mathrm o}/\mathrm dt$.
An ideal closed disconnect gives $\mathbf v_{\mathrm{t}}=\mathbf v_{\mathrm{PCC}}$;
when open, $\mathbf i_{\mathrm{PCC}}=0$ and the two voltages can differ.

Primary Control follows the UNIFI cascade: local $\mathbf v_{\mathrm{t}}$ and
$\mathbf i_g$ determine power feedback; $P^{\mathrm{ref}}$ and $Q^{\mathrm{ref}}$
set the operating target. Droop, VSM, or another primary law generates
$\mathbf v^{\mathrm{ref}}$, $\omega$, and $\theta$, with
$\mathrm d\theta/\mathrm dt=\omega$.[^unifi] Nominal voltage and frequency are
configured references. The shared $\theta$ supplies all implicit Park transforms
and PWM; $\omega$ supplies current- and voltage-loop decoupling.
The inner outline groups PWM, Current Control, and Voltage Control in this
shared reference frame; angle and frequency connections are implicit.
Voltage Control receives $\mathbf v^{\mathrm{ref}}$ and measures
$\mathbf v_{\mathrm o}$ and $\mathbf i_g$, matching OuterVoltageControl.

Plant Control regulates PCC exchange and supplies $P^{\mathrm{ref}}$ and
$Q^{\mathrm{ref}}$ to Primary Control. REPCA is one plant-controller model;
its `pext` and `qext` are power commands, not voltage-loop inputs. The diagram
shows active/reactive-power dispatch; plant voltage and frequency-response
modes require their corresponding references and measurements.

Primary Control and Plant Control are architectural blocks here, not a wired
switching example. Islanded operation needs a local network and load, an
autonomous primary controller, and coordination of plant references when the
PCC opens. Reconnection also requires synchronization across the disconnect.

[^pcc]: NREL, [*Clause-by-Clause Summary of Requirements in IEEE Standard 1547-2018*](https://www.nlr.gov/docs/fy20osti/75184.pdf), 2020, Sections 1–2; EPRI, [*Definitions, Acronyms, and References*](https://der-interconnection.epri.com/hub/tiir/en/ch2-definitions.html).

[^unifi]: UNIFI Consortium, [*UNIFI's Grid-Forming Inverter Reference Design*](https://docs.nlr.gov/docs/fy25osti/92994.pdf), 2025, Sections 2.1–2.2, Figures 2–4.

[^ieee2800]: NERC–EPRI–NATF, [*2023 Planning and Modeling Virtual Seminar*](https://www.nerc.com/globalassets/who-we-are/standing-committees/rstc/irps/2023_natf-epri-nerc_pm_virtual_seminar-day_1.pdf#page=135), IEEE 2800 reference points of applicability.

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
