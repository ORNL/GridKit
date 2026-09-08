# Hawaii EMT case

This case ports the 37-bus Hawaii system at `f16f3815e6f84c08bd0869773d3033b3dcb3dd0f`
on `lukel/cases-polish-dev`. It retains 30 synchronous machines, 30 TGOV1
governors, 30 IEEET1 exciters, 14 IEEEST stabilizers, nine inverter plants,
27 constant-impedance loads, and 89 branches, including 12 transformers.
The source system base is 100 MVA and 60 Hz. Bus voltage bases are 138 or 69 kV.

`convert.py` reads the source JSON with `git show`; it does not check out the
source branch. It writes `Hawaii.case.json`, `Hawaii.state.json`, and
`conversion.json`. The latter records every machine conversion and operating
point adjustment. To regenerate from an exported source JSON, pass `--source`.

## Network and operating point

The 77 lines are balanced, uncoupled three-phase pi sections. With
$Z_\mathrm{b}=V_\mathrm{b}^2/S_\mathrm{b}$, their phase parameters are
$R=R_\mathrm{pu}Z_\mathrm{b}$, $L=X_\mathrm{pu}Z_\mathrm{b}/\omega_\mathrm{b}$,
and total $C=B_\mathrm{pu}/(\omega_\mathrm{b}Z_\mathrm{b})$. Length is one metre;
these are lumped equivalent coefficients, not inferred conductor geometry.
Positive- and zero-sequence impedances are identical because the source has
only positive-sequence data. This case is for a balanced fault.

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
the EMT Machine requires physical winding circuits. Fundamental and standard
parameter sets are distinct, as described in the
[MathWorks parameterization documentation](https://www.mathworks.com/help/sps/ug/machine-parameterization.html).
The following classical separated-time-scale conversion is an approximation,
not an exact GENROU reduction. For either axis, let $X$, $X'$, $X''$, and $X_l$
denote the source reactances and $T'_0$, $T''_0$ its open-circuit time constants:

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
open-circuit poles of each coupled winding circuit; these differ from the
separated-time-scale inputs. The reactance limits are checked independently.

For the 18 units with $X'_d=X''_d$, the effective subtransient reactance is
$X''_{d,\mathrm{eff}}=X_l+0.99(X'_d-X_l)$. This changes only one percent of
the gap above leakage and gives finite positive damper leakage. Each changed
value appears in `conversion.json`. Saturation acts on both EMT magnetizing
axes, and stator transients are retained; those dynamics are absent from the
phasor reduction. Exciter sensing uses the instantaneous balanced voltage
norm. Existing controller time-constant floors and shared smooth limiters
remain those documented by their model READMEs. In particular, IEEET1 raises
the source's zero voltage-sensing time constant to 0.001 s and reports that
change at startup.

The largest regularization changes $X''_d$ by 0.25 percent. The largest
time-constant discrepancy is much larger: units `37_3` and `37_5` have a
coupled q-axis slow time constant of 1.04168 s instead of 0.58 s (79.6 percent),
and a fast time constant of 0.03898 s instead of 0.07 s. These are material
model deviations. The approximate winding conversion should not be used as
evidence of matching GENROU fault dynamics.

## Inverter plants

Each source REGCA unit is replaced by a PLL, voltage and current Park
transforms, OuterPowerControl, InnerCurrentControl, inverse Park, Modulation,
PWM, Converter, DCLink, and a DependentVoltageSource with an RL filter.
The filter connects directly to the original bus. Its measured injection is
the outer loop's measured dq-current feedback. Current targets `Pref/V` and
`-Qref/V` are derived from power-reference parameters and rated voltage.
Those parameters are chosen to reproduce the initialized dq injection.
The targets remain fixed as voltage changes; this replacement does not
regulate constant P/Q.
Power-invariant current base is $S/V$.

The following plant data are fabricated because REGCA supplies no bridge,
filter, DC energy, or switching-control parameters:

Parameter | Value
--------- | -----
Filter resistance and reactance | 0.01 and 0.20 p.u. on the plant rating
Inner-loop bandwidth | $2\pi\,300$ rad/s; $K_P=L\omega_c$, $K_I=R\omega_c$, $K_{\mathrm{aw}}=\omega_c$
Current limit | Source REECB `Imax` times $S/V$
Modulation limit | 0.95
PLL | $K_P=80$ rad/s, $K_I=2500$ rad/s$^2$
Outer current loop | $K_P=0.01$, $K_I=40$ s$^{-1}$, $K_{\mathrm{aw}}=200$ s$^{-1}$
DC voltage | $2V_\mathrm{b}$
DC energy | $H_\mathrm{dc}=10$ s; $C=2H_\mathrm{dc}S/V_\mathrm{dc}^2$
DC source | Constant current equal to initial bridge power divided by initial DC voltage
Carrier | 1800 Hz, centre aligned

The study uses shared $\mu=50000$ s$^{-1}$, giving a logistic 10--90 percent
edge width of $2\ln(9)/\mu=87.9$ microseconds. Relative and absolute solver
tolerances are $10^{-5}$ and $10^{-6}$; the monitor interval is $1/7200$ s.
Adaptive accepted steps resolve the edges more finely than the monitor
interval. Accepted-step statistics are reported separately. This study does
not use the low-$\mu$ fundamental-only setting.

The large DC energy represents an aggregate energy buffer; it is not a
manufacturer capacitor value or a DC-voltage regulator. Initial current-loop
integrals supply filter resistive drop, and outer-loop integrals supply the
dispatched current. Initial conditions describe the fundamental operating
point, not the periodic switching orbit.

The replacement omits REGCA current-source lag, LVPL, high-voltage reactive
current logic, and source ramp-rate behavior; REECB voltage-dip reactive
injection, P/Q priority, voltage/Q control flags, measurement lags, and
current-order recovery; and REPCA plant voltage/reactive-power regulation,
line-drop compensation, voltage freeze, deadbands, and frequency dispatch.
Its PLL, outer current PI loop, circular current limiter, and physical bridge are
different dynamics. Matching initial dispatch does not validate those omitted
controls.

## Fault and comparison

At 1.0 s, `fault_switch` closes a three-phase grounded resistive `LoadZ` at
bus 1; it opens at 1.10 s. Each phase resistance is
$0.01 V_{\mathrm{b},1}^2/S_\mathrm{b}=1.9044$ ohm, matching the magnitude of
the source inductive fault impedance. This resistive choice changes the fault current angle and dissipation. It permits
commanded clearing without stranding inductive energy in an ideal open switch.
No transmission branch trips.

The supplied PowerWorld reference and source validation solver clear at
**1.15 s**, whereas this request specifies **1.10 s**. The four comparison
plots label both intervals. Differences include timing, fault impedance,
machine conversion, core branches, and replacement inverter controls; their
error statistics are descriptive and are not a PowerWorld validation claim.
Machine speed is plotted as $\omega_r-1$, powers on the 100 MVA system base,
and bus voltage as the positive-sequence magnitude on its local voltage base.

Deviation | EMT choice and consequence
--------- | --------------------------
GENROU reduction | Fundamental winding circuits, stator transients, shared-axis saturation, approximate open-circuit times
Equal d-axis reactances | 18 positive-damper regularizations, at most 0.25 percent change in $X''_d$
Transformer no-load and connection data | Fabricated core branches and grounded-wye banks; source voltage bases retain nominal ratio
Initial power flow | At bus 23, total synchronous dispatch changes by -0.152426 MW and +1.179855 Mvar; maximum bus phasor change is 0.000156170 p.u.
Positive-sequence network | Uncoupled balanced three-phase equivalents; no source-derived zero-sequence or frequency-dependent data
Exciter sensing | Existing 1 ms sensing floor replaces source zero lag
REGCA/REECB/REPCA | Explicit switching plants with fabricated filter/control/DC data and the omitted functions listed above
Fault impedance | Resistive magnitude equivalent replaces inductance
Fault clearing | Requested 1.10 s replaces reference 1.15 s
Measurements | Centred cycle averages remove switching ripple from the phasor comparison

## Run and validation

From this directory:

```bash
python3 convert.py
../../../build/application/EMT/EMTDynamicSimulation Hawaii.solver.json
python3 validate.py --exe ../../../build/application/EMT/EMTDynamicSimulation
```

The validator uses only the Python standard library. It checks conversion
invariants, finite trajectories, event timing, initial dispatch, current-reference
limits, DC voltage, and the fault/recovery response. The monitored current-limit
ratio allows $10^{-4}$ numerical error: the full baseline run exceeds unity
by $5.543\times10^{-5}$ at 2.000556 s. Through 1.5 s, the maximum excess is
$6.813\times10^{-6}$, decreasing to $3.561\times10^{-11}$ with tenfold tighter
solver tolerances. This refinement covers early recovery only. Solver monitor cadence
does not measure accepted solver steps. The full ten-second study and
PowerWorld comparison artifacts live in `examples/EMT/Hawaii`; raw simulation
CSVs and logs belong in a local run directory.
