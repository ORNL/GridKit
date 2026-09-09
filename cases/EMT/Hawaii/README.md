# Hawaii EMT case

This case ports the 37-bus Hawaii system at `f16f3815e6f84c08bd0869773d3033b3dcb3dd0f`
on `lukel/cases-polish-dev`. It retains 30 synchronous machines, 30 TGOV1
governors, 30 IEEET1 exciters, 14 IEEEST stabilizers, nine inverter plants,
nine REPCA plant controllers, 27 constant-impedance loads, and 89 branches, including 12 transformers.
The source system base is 100 MVA and 60 Hz. Bus voltage bases are 138 or 69 kV.

`convert.py` reads the source JSON with `git show`; it does not check out the
source branch. It writes `Hawaii.case.json`, `Hawaii.state.json`, and
`conversion.json`. The latter records every machine conversion and operating
point adjustment. To regenerate from an exported source JSON, pass `--source`.

## Network and operating point

The 77 lines are coupled, ideally transposed three-phase pi sections.
`line_parameters.py` uses GridWorkbench's overhead-line calculations at
60 Hz, including conductor internal impedance, earth return, and mutual
capacitance. It fits synthetic conductor radius, phase spacing, and length
to each distinct source R/X/B combination. The horizontal three-wire
geometry has no shield wires; heights are 14 m at 69 kV and 20 m at 138 kV.
`line_parameters.json` records the assumptions, bounds, calculator hashes,
geometries, calibration factors, and full per-metre matrices.

The source's positive-sequence values do not determine a physical geometry.
Many branches cannot be matched within the chosen overhead-geometry bounds.
The generator therefore scales each complete R, L, and C matrix by a positive
factor to preserve the source positive sequence exactly. These are calibrated
lumped equivalents, not raw geometry results or surveyed Hawaii lines.
The largest uncalibrated relative error is 243%; the calibration must not be
interpreted as an independently validated zero-sequence or broadband model.

For each per-metre matrix $\mathbf M$, ideal transposition gives equal self
terms $M_s$ and mutual terms $M_m$. Its positive-sequence value is $M_s-M_m$;
its zero-sequence value is $M_s+2M_m$. Calibration preserves the geometry's
ratio between these values and its positive semidefiniteness. With
$Z_\mathrm{b}=V_\mathrm{b}^2/S_\mathrm{b}$, the total positive-sequence values
remain $R=R_\mathrm{pu}Z_\mathrm{b}$,
$L=X_\mathrm{pu}Z_\mathrm{b}/\omega_\mathrm{b}$, and
$C=B_\mathrm{pu}/(\omega_\mathrm{b}Z_\mathrm{b})$.
`LineLumped` applies the physical fitted length and half-length terminal shunts.
This case is for a balanced fault.

All source ZIP fractions have `alphaI = alphaP = 0`. Each load resistance
preserves its initial active power at the source voltage magnitude. The
negative reactive load at bus 16 is a separate bus capacitance with the same
initial reactive power. Line shunts belong to the terminal buses.

Transformer resistance and leakage reactance retain the source values on a
100 MVA bank base. All source taps are one and phase shifts are zero. The
EMT banks use grounded-wye identity connection maps and equal core splits.
The missing no-load and saturation data are explicitly fabricated:
`I0 = 0.001` p.u., `P0 = 0` W, `knee = 1.2` p.u., and `Lsat = 0.25` p.u.
The rating is a conversion base, not a thermal rating. The no-load branch adds
reactive demand absent from the source.

A rectangular-coordinate power flow includes these core branches and the
original load admittances. Bus 23 retains its source complex voltage; all
other buses retain their specified generator powers. The small slack-power
adjustment is shared among the eight synchronous machines at bus 23 in
proportion to rating. Inverter dispatch is unchanged. `conversion.json`
records the voltage changes, slack adjustment, and current-balance residual.
Balanced phasors initialize bus voltages, line currents, transformer leakage
currents and core fluxes, and machine current injections. Governors and
exciters use the normal machine operating-point reconciliation.

## Machine winding conversion

The source GENROU parameters describe transient and subtransient responses;
the EMT Machine requires physical winding circuits. The conversion preserves
the unsaturated coupled rotor dynamics on each axis, subject to the
reactance regularization below. Let $X$, $X'$, $X''$, and $X_l$ denote the
source reactances and $T'_0$, $T''_0$ its time parameters:

```math
\begin{aligned}
L_m &= X-X_l, &
L_1 &= \frac{L_m(X'-X_l)}{X-X'}, &
L_2 &= \frac{(X'-X_l)(X''-X_l)}{X'-X''}, \\
R_1 &= \frac{L_m+L_1}{\omega_\mathrm{b}T'_0}, &
R_2 &= \frac{L_2+L_m L_1/(L_m+L_1)}{\omega_\mathrm{b}T''_0}.
\end{aligned}
```

On the d axis, winding 1 is the field and winding 2 is the damper. On the q
axis they are the two dampers. `Ll = L0 = Xl`; rating, inertia, stator
resistance, and saturation data are copied. All source damping coefficients
are zero, so EMT friction is zero. Differential leakage is zero as required
by the existing Machine model. The conversion records the two unsaturated
open-circuit poles separately from the time parameters. GENROU retains rotor
coupling, so its poles are the eigenvalues of

```math
A=\begin{bmatrix}
-(1+a)/T'_0 & a/T'_0 \\
1/T''_0 & -1/T''_0
\end{bmatrix},\qquad
a=\frac{(X-X')(X'-X'')}{(X'-X_l)^2}.
```

`conversion.json` distinguishes the source poles, the effective poles after
regularization, and the winding poles, in inverse seconds from slow to fast.
The validator checks the actual winding parameters against these poles and
independently compares operational reactance and d-axis field response with
the unsaturated GENROU equations.

For the 18 units with $X'_d=X''_d$, the effective subtransient reactance is
$X''_{d,\mathrm{eff}}=X_l+0.99(X'_d-X_l)$. This changes only one percent of
the gap above leakage and gives finite positive damper leakage. Each changed
value appears in `conversion.json`. EMT saturation uses air-gap flux and
scales both magnetizing inductances; GENROU uses subtransient flux. EMT also
retains stator transients and fault-induced DC current offsets. All source
stator resistances are zero and are retained; the phasor stator equations
omit these transient currents. Exciter sensing uses the instantaneous balanced voltage
norm. Existing controller time-constant floors and shared smooth limiters
remain those documented by their model READMEs. In particular, IEEET1 raises
the source's zero voltage-sensing time constant to 0.001 s and reports that
change at startup.

The local phasor study and EMT conversion set the 22 negative IEEET1 `Ke`
parameters to zero, selecting the model's existing automatic initialization.
At this operating point those exciters are below the saturation knee; retaining
negative `Ke` leaves an unstable field-voltage difference mode between colocated
machines. The eight positive `Ke` values are preserved. This is a documented
study adjustment, not recovered equipment data; `conversion.json` records each
changed value. The original exported validation inputs remain intact.

The largest regularization changes $X''_d$ by 0.25 percent and an open-circuit
pole by 0.371 percent. The effective GENROU and winding poles agree to
roundoff. This rotor correspondence does not establish equivalence of the
full fault trajectories, which include the other model differences below.

## Inverter plants

Each source REGCA unit is replaced by the
[GFL switching-inverter arrangement](../../../GridKit/Model/EMT/Component/Controller/README.md#grid-following): PLL,
voltage and current Park transforms, REECB, InnerCurrentControl,
PWM, Converter, and a physical LCL Filter.
`Converter.e` drives `Filter.e`; `Filter.ig` injects into the original bus.
`Filter.i` supplies the inner current loop, while `Filter.ig`
supplies the outer loop through a separate Park transform. PLL reads terminal
Bus voltage. Separate voltage Park transforms supply Bus voltage to
REECB and `Filter.vo` to InnerCurrentControl. All Park transforms
and PWM share PLL's `theta`, and the inner controller receives PLL's `omega`.
PWM limits the dq voltage command and returns it to the inner controller.
REPCA reads the same terminal voltage and grid current as REECB.
Its source parameters, selectors, and plant rating are retained. All nine
plants regulate terminal voltage (`RefFlag = true`); their initialized voltage
references are latched locally. REPCA `qext` supplies REECB `Qref`,
while `Pref` is an external constant at the source active dispatch. REECB retains the source measurement lags, voltage-response logic, and current
priority. All nine plants use reactive-power references and Q priority.
Power-invariant current base is $S/V$.

The following plant data are fabricated because REGCA supplies no bridge,
filter, or switching-control parameters:

Parameter | Value
--------- | -----
Converter-side resistance and reactance | 0.01 and 0.20 p.u. on the plant rating
Grid-side resistance and reactance | 0.005 and 0.10 p.u. on the plant rating
Capacitor susceptance | 0.10 p.u.; $C=0.10/(\omega_\mathrm{b}Z_\mathrm{b})$
Undamped LCL resonance | 734.85 Hz; $\sqrt{(L_{\mathrm{s}}+L_g)/(L_{\mathrm{s}}L_g C)}/(2\pi)$
Inner-loop bandwidth | $2\pi\,300$ rad/s; $K_P=L_{\mathrm{s}}\omega_c$, $K_I=R_{\mathrm{s}}\omega_c$, $K_{\mathrm{aw}}=\omega_c$
Capacitor-voltage measurement | 5 ms lag for current-reference compensation
Bridge current limit | $1.5S/V$; separate from source REECB terminal-current limit $1.3S/V$
Modulation limit | 0.95
PLL | $K_P=80$ rad/s, $K_I=2500$ rad/s$^2$
DC voltage | Constant $2V_\mathrm{b}$, shared by PWM and Converter
Carrier | 1800 Hz, centre aligned

The default study uses $\mu=240$ s$^{-1}$ and BDF order at most 2.
It retains `rel_tol = 1e-5` and uses model-scaled `abs_tol = 1e-6` with
`scaled_abs_tol = true`. Conversion supplies nominal Bus voltages, network
phase-current bases on 100 MVA, and inverter phase-current ratings on each
plant base. The monitor interval is $1/7200$ s; accepted steps are recorded
separately.

The separate switching check uses $\mu=50000$ s$^{-1}$, order 5, and scalar
`abs_tol = 1e-6`. Its logistic 10--90 percent edge width is
$2\ln(9)/\mu=87.9$ microseconds.

With power-invariant balanced phasors and the original terminal injection $I_g$, initialization uses

```math
\begin{aligned}
V_{\mathrm{o}} &= V+(R_g+\mathrm{j}\omega_{\mathrm{b}}L_g)I_g, \\
I &= I_g+\mathrm{j}\omega_{\mathrm{b}}CV_{\mathrm{o}}, \\
E &= V_{\mathrm{o}}+(R_{\mathrm{s}}+\mathrm{j}\omega_{\mathrm{b}}L_{\mathrm{s}})I, \\
P_{\mathrm{bridge}} &= P_{\mathrm{grid}}+R_g|I_g|^2+R_{\mathrm{s}}|I|^2.
\end{aligned}
```

The nine Filter states follow these phasors. PLL aligns to terminal Bus voltage;
inner-loop voltage output initializes to $E$ in that frame. Its capacitor
compensation maps REECB terminal-current commands to converter-current commands
using $\mathbf i^\star=\mathbf i_g^{\mathrm{cmd}}+\omega C\mathbf J\mathbf v_f$.
Here $\mathbf v_f$ is the filtered capacitor voltage, equal to $\mathbf v_o$
in steady state. Initialization reverses that map and requests terminal current from REECB.
REECB derives its states and requests terminal Q from REPCA, which derives its
PI/lead-lag states and voltage reference through the normal initialization procedure. Initial conditions describe
the fundamental operating point, not the periodic switching orbit.

The replacement omits REGCA current-source lag, LVPL, high-voltage reactive
current logic, and source reactive-current rate limits. REPCA and REECB control
equations and parameters are retained; REPCA frequency control is disabled in
every source plant. The explicit PLL, inner current regulator, LCL filter, and
switching bridge supply additional dynamics. Capacitor-reference compensation
is a balanced-fundamental relation; capacitor transients remain physical.
Matching initial dispatch does not establish complete REGCA equivalence. See
[REGCA correspondence](../../../GridKit/Model/EMT/Component/Controller/README.md#regca-correspondence)
for LVPL and LVACM.

## Fault and comparison

At 1.0 s, `fault_switch` connects a three-phase grounded inductive `LoadZ` at
bus 1; it opens at 1.15 s. The fault uses zero resistance and diagonal
inductance $L_f=0.01 V_{\mathrm{b},1}^2/(\omega_\mathrm{b}S_\mathrm{b})
=5.05158$ mH per phase, matching the source $\mathrm{j}0.01$ p.u. impedance
at 60 Hz. The inductors start de-energized. No transmission branch trips.

`fault_discharge_switch` connects a separate grounded resistive `LoadZ` to
`fault_bus` whenever the grid fault switch is open. It opens at fault inception
and closes at clearing, in the same event groups as `fault_switch`.
The discharge resistance is $R_d=1.9044$ ohm per phase; it is disconnected
throughout the applied fault. After clearing, the isolated inductor currents
remain continuous and decay according to $L_f\dot{\mathbf{i}}+R_d\mathbf{i}=0$,
with $L_f/R_d=2.65258$ ms. This discharge circuit is an explicit idealized
clearing assumption; its stored energy is dissipated locally after the grid
connection opens.

The source PhasorDynamics validation solver and EMT solver both clear at
**1.15 s**. The comparison plots use freshly simulated GridKit PhasorDynamics
trajectories and label the actual intervals recorded by each run. Differences include
machine conversion, core branches, and replacement
inverter controls; the error statistics describe these different models.
The phasor comparison uses the same exciter adjustment as EMT. Both faults have the same fundamental
impedance and clearing time; EMT additionally retains the inductor transient.
Machine speed is plotted as $\omega_r-1$, powers on the 100 MVA system base,
and bus voltage as the positive-sequence magnitude on its local voltage base.

Deviation | EMT choice and consequence
--------- | --------------------------
GENROU reduction | Corresponding unsaturated rotor circuits, retained stator transients, shared-axis saturation
Equal d-axis reactances | 18 positive-damper regularizations, at most 0.25 percent change in $X''_d$
Transformer no-load and connection data | Fabricated core branches and grounded-wye banks; source voltage bases retain nominal ratio
Initial power flow | At bus 23, total synchronous dispatch changes by -0.152426 MW and +1.179855 Mvar; maximum bus phasor change is 0.000156170 p.u.
Positive-sequence network | Coupled transposed equivalents calibrated to source R/X/B; zero sequence inferred from synthetic geometry
Exciter sensing | Existing 1 ms sensing floor replaces source zero lag
REGCA | Explicit switching plants with fabricated filter/inner-control/DC data and the omitted functions listed above
REECB | Source electrical-control equations and parameters retained; SI terminal inputs and current commands
REPCA | Source plant-control equations and parameters retained; SI terminal measurements and power outputs
Fault impedance | 5.05158 mH per phase reproduces the source inductive impedance at 60 Hz
Fault clearing circuit | Isolated switched discharge resistors preserve inductor current and dissipate stored energy
Fault clearing | 1.15 s, matching the PhasorDynamics validation solver
Measurements | Centred cycle averages remove switching ripple from the phasor comparison

## Run and validation

From this directory:

```bash
python3 convert.py
../../../build/application/EMT/EMTDynamicSimulation Hawaii.solver.json
python3 validate.py --exe ../../../build/application/EMT/EMTDynamicSimulation
```

Ordinary conversion uses the stored `line_parameters.json` and requires only
Python's standard library. To regenerate those parameters, use NumPy, SciPy,
and a GridWorkbench checkout containing `gridworkbench.emt.parameters`:

```bash
python3 line_parameters.py --gridworkbench /path/to/GridWorkbench
python3 convert.py
```

Pass the same `--source` to both scripts when using an exported source case.
The converter rejects coefficients generated from a different source file.

The validator uses only the Python standard library. It checks conversion
invariants, LCL wiring and initial KVL/current/DC-power balances, finite
trajectories, paired event timing, fault-current continuity and analytical
discharge decay, terminal dispatch, current-reference limits,
DC voltage, and the fault/recovery response. The monitored current-limit
ratio allows $10^{-4}$ numerical interpolation error. Solver monitor cadence
does not measure accepted solver steps. The full five-second study and
GridKit PhasorDynamics comparison scripts live in `examples/EMT/Hawaii`; the four comparison figures are tracked, while metrics, raw simulation CSVs,
and logs remain in ignored result directories.
