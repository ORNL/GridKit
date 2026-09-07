# GridKit equivalent 9-bus case

The complete 3-second governor-step and ideal-opening cases have been simulated in GridKit
with three machines, three SEXS-PTI exciters, three GASTPTI governors,
three IEEEST stabilizers, six pi lines, three transformer series branches,
and three series RL loads. The scripts and results are in this directory.
Both output files now use a 50 microsecond monitor interval. ParaEMT
uses the number of saved snapshots to select one-, two-, then
three-point predictor startup in `predictX`; this wrapper retains that
behavior. The no-capture timing path preserves the same snapshot-count
transitions without retaining result arrays, and its final state must
match the capture run. This output-cadence dependence is an upstream
numerical detail, not a GridKit adaptive-step restriction.

The reference is the authors' pinned ParaEMT revision
[`d79d735a`](https://github.com/NatLabRockies/ParaEMT_public/tree/d79d735a4a587d56c5b88187d1a499195b6b2b84),
not an arbitrary IEEE 9-bus parameter set.

The normal `EMTDynamicSimulation` application cannot run this case by
itself. [gridkit_9bus.cpp](gridkit_9bus.cpp) supplies steady-state network
currents/derivatives, writable governor references, the scheduled step,
and the exact generator-terminal constraint reformulation below. All
machine and controller residuals come from the actual GridKit libraries.
No controllers are frozen or replaced by constant outputs.

| Part | GridKit implementation |
| --- | --- |
| GENROU equivalent circuits | Existing EMT `Machine`, with converted fundamental winding parameters |
| SEXS | EMT `SEXS-PTI`, direct PhasorDynamics port, plus the actual 0.02 s terminal-magnitude measurement lag |
| GAST | EMT `GASTPTI`, full droop, valve antiwindup, fuel-flow and temperature lags, temperature selection, damping, rating conversion, and initialization |
| IEEEST | EMT `IEEEST`, complete notch/lead–lag/washout/output-limit cascade; absolute rotor speed is converted to deviation inside the model |
| Lines and transformers | Existing `LineLumped`, diagonal phase matrices, transformer shunts zero |
| Loads | Existing `LoadZ`, series RL, using ParaEMT load option 1 |
| Disturbances | Bus-1 governor reference increment −0.02 pu, or an ideal generator-terminal opening, at 1 s |

The three model commits are `bef5be65` (`SEXS-PTI Implementaiton`),
`17e4f901` (`GASTPTI Implementaiton`), and `34a1f088`
(`IEEEST Implementaiton`). They include JSON/builder registration,
monitoring, analytic/dependency-tracking and sparse Enzyme Jacobians,
documentation, and tests against the PhasorDynamics implementations.
IEEEST also implements compensated-voltage cutout using CommonMath
smoothing. The separate PSLF transport-delay extension is explicitly
rejected for nonzero `Tdelay`; the reference case uses no delay and
has both voltage cutouts disabled.

## Machine parameters

[make_gridkit_case.py](make_gridkit_case.py) reads the original worksheet
exports, solved JSON, and saved `machine_parameters.json`. The latter
contains the actual matrices assembled by ParaEMT's `ToEquiCirData` and
`MergeMacG` routines in [Lib_BW.py](upstream/Lib_BW.py).

| GridKit parameter | ParaEMT expression |
| --- | --- |
| `Ll`, `Lmd`, `Lmq`, `L0`, `Rs` | `ec_Ll`, `ec_Lad`, `ec_Laq`, `ec_L0`, `ec_Ra` |
| `Llfd` | `ec_Lffd - ec_Lad` |
| `Ll1d` | `ec_L11d - ec_Lad` |
| `Ll1q` | `ec_L11q - ec_Laq` |
| `Ll2q` | `ec_L22q - ec_Laq` |
| `Rfd`, `R1d`, `R1q`, `R2q` | Corresponding `ec_R* / ws` |
| `S`, `V`, `f` | Original 150/250/100 MVA ratings in VA; referred 230 kV LL RMS; 60 Hz |
| `H`, `F`, `S10`, `S12` | Workbook `H`; zero damping and saturation for this case |

GridKit's rotor equations use `(1/ws)*d(psi)/dt`; ParaEMT uses an unscaled
rotor flux derivative. The resistance and exciter conversions therefore
change together by `ws`. Compare the published exciter-base `efd` output,
not the differently scaled internal rotor voltage.

`Ll1d` is **0.1083333333**, not the intermediate ParaEMT `ec_L1d=0.2`.
`Ll2q` is **0.1086956522**, not `ec_L2q=0.125`. These values reproduce
the assembled winding self-inductances. The comparison uses external
voltages, speed, power, field voltage, and stabilizer output; internal
rotor angles are not compared across different Park conventions.

## Network and operating point

ParaEMT uses dimensionless bus voltages and series RL transformer
branches with unit taps on their respective voltage bases. Referring all
bus and machine voltages together to 230 kV preserves this network. The
network base is 100 MVA, so `Z_base=529 ohm`. Physical terminal voltages
can be recovered on the original 16.5/18/13.8/230 kV bases from the saved
per-unit waveforms.

For each line, `R=Re(Zpu)*Z_base`, `L=Im(Zpu)*Z_base/ws`, and total
`C=Bpu/(ws*Z_base)`, split equally between the terminals. Transformer
branches have their original series impedance and zero shunt. Loads use
`Zpu=Vm^2/conjugate((P+jQ)/100 MVA)`, retaining the **original** load
impedances when refining the operating point.

The supplied JSON rounds voltages and dispatch independently. Its initial
KCL mismatch prevents a consistent inductive start. The converter solves
the passive network by Schur reduction, fixing all three generator
voltage magnitudes, generator-2/3 real power, and generator-1 angle.
[initialization.json](gridkit/initialization.json) records every change:
maximum voltage-phasor adjustment **4.8466e-6 pu**, slack generation
adjustment **0.00156184 MW**, and maximum reactive-generation adjustment
**0.00069636 Mvar**. No dynamic parameter is fitted to the reference.
The resulting passive-node current mismatch is about `1.1e-14 pu`.

## Generator-terminal current constraints

At buses 1–3, the machine and transformer impose a constraint on
inductive states. This connection has a higher-index DAE; direct IDA
consistency calculation and direct integration failed, even after the
balanced initial residual was reduced to `5.82e-11` in mixed physical
units. IDA's documented consistency algorithm targets index-one systems:
[SUNDIALS mathematical considerations](https://sundials.readthedocs.io/en/latest/ida/Mathematics_link.html#initial-condition).

[GeneratorTerminalConstraint.hpp](GeneratorTerminalConstraint.hpp)
differentiates only these three-phase KCL constraints, using the exact
unsaturated winding equations. Write

```math
\begin{bmatrix}\psi_d\\\psi_{fd}\\\psi_{1d}\end{bmatrix}
= M_d\begin{bmatrix}i_d\\i_{fd}\\i_{1d}\end{bmatrix},\qquad
\begin{bmatrix}\psi_q\\\psi_{1q}\\\psi_{2q}\end{bmatrix}
= M_q\begin{bmatrix}i_q\\i_{1q}\\i_{2q}\end{bmatrix}.
```

The first rows of `inverse(M_d)` and `inverse(M_q)` give `di_d/dt`
and `di_q/dt` from the corresponding flux derivatives. For phase angle
`theta_k`, the exact current derivative is

```math
\dot i_k = \cos\theta_k\,\dot i_d-\sin\theta_k\,\dot i_q
-\omega_b\omega(\sin\theta_k\,i_d+\cos\theta_k\,i_q)
-\dot\psi_0/L_0.
```

Original terminal KCL is `I_base*i_k+i_transformer,k=0`. Its residual is
replaced by `(I_base*di_k/dt+di_transformer,k/dt)/omega_b=0`; the original
zero initial KCL is preserved in exact arithmetic. The transformer shunt
currents are identically zero. The helper adds **no states** and no
physical admittance. It substitutes the flux expressions above rather
than treating algebraic current derivatives as extra independent states.
The driver rejects saturated machines or nonzero transformer shunts.
This is a case-specific reformulation, not a general index-reduction
facility in the EMT application.

The governor-step system has **258 variables** and **1086 Jacobian entries**.
The driver checks the whole Jacobian against central differences at the
initial and final states, including `dF/dy + alpha*dF/dyp`. Maximum scaled
difference is `3.81e-6` (the denominator is `1+abs(FD)`, with an absolute
perturbation of `1e-6`); this checks differentiation, not solution accuracy.
It also verifies that event consistency calculation changes **no
differential state**. `run.json` records these checks.

Original generator-terminal KCL is monitored independently throughout
the run. Maximum sampled mismatch falls from **0.009896 A** to
**0.000609 A** to **0.000162 A** as GridKit tolerance goes from `1e-7`
to `1e-8` to `1e-9`. These are currents on the referred 230 kV base.
They quantify numerical drift in the differentiated constraint.

## Ideal generator-terminal opening

The trip case adds a three-variable `machine_bus_1`. Before the event,
its three equations enforce equality with network `bus_1`, so the circuit
is the same. At 1 s the helper separates these terminals and uses

```math
\dot i_{\mathrm{transformer},k}/\omega_b=0,\qquad
I_{\mathrm{base}}\dot i_{\mathrm{machine},k}/\omega_b=0
```

at the network and machine terminal, respectively. The currents are
explicitly projected to zero at the opening, so these differentiated
constraints preserve the open-circuit KCL. The original full machine and
controller equations remain active, with the exciter sensing the isolated
machine terminal. The helper retains a fixed sparse structure in both
configurations. The trip system has **261 variables** and **1128 Jacobian
entries**, with assembled-Jacobian finite-difference checks before and
after opening.

Rotor winding fluxes cannot jump without impulsive rotor voltages.
For the d axis, the open-circuit projection therefore solves

```math
\begin{bmatrix}L_{md}+L_{lfd}&L_{md}\\L_{md}&L_{md}+L_{l1d}\end{bmatrix}
\begin{bmatrix}i_{fd}^{+}\\i_{1d}^{+}\end{bmatrix}
=\begin{bmatrix}\psi_{fd}^{-}\\\psi_{1d}^{-}\end{bmatrix},\qquad
\psi_d^{+}=L_{md}(i_{fd}^{+}+i_{1d}^{+}).
```

The q axis uses the corresponding two damper windings; `psi_0` becomes
zero. Only these three stator fluxes and the three transformer currents
are projected. Rotor fluxes, angle, speed, all controller states, other
inductor currents, and network capacitor voltages retain their left limits.
The subsequent IDA consistency solve changes no differential state from
that explicit projection. `event_state_limits.json` saves all before,
projected, and consistent values; `run.json` identifies the six affected
indices. `verify.py` checks this independently against the winding
parameters and verifies zero stator current after consistency.

This ideal interruption is a finite pre/post-event construction. A voltage
impulse is required to interrupt stored inductive current instantaneously;
its finite peak and waveform require a specified arc/snubber/breaker model.
No such model, pulse width, or artificial admittance has been added here.
Continuous controller equations are integrated on each side of the event;
this construction does not claim a resolved measurement-chain response
to an unspecified voltage impulse.

ParaEMT's event is different. `Lib_BW.GenTrip` zeros G1 injection and
rebuilds the network without its Norton conductance. `lib_numba.numba_updateIg`
skips G1 afterward, leaving its electrical history unchanged, while
`numba_updateX` still updates its states using that history and network
bus voltage (the skip is commented out). `Re_Init` also changes network
history. G1's resulting internal trajectory does not represent the same
isolated full machine as GridKit. It is plotted explicitly, with that
limitation; it is not used as an equivalent-machine answer key.

The ParaEMT bus-1 event-step magnitude grows from **3.595605** to
**8.062099** to **17.045898 pu** at 50, 25 and 12.5 microseconds. Every
integration step from 0.995 to 1.015 s is retained in the native event
CSV, including the spike. On the common 50 µs grid, GridKit retains its
left-limit sample at 1 s, whereas ParaEMT includes its event step. Raw
GridKit monitoring retains both restart limits; neither those records
nor the ParaEMT spike are silently removed from the archived data.

## Agreement and remaining limits

The three GridKit tolerance runs and three ParaEMT time-step runs are
compared in [gridkit_comparison.json](results/gridkit_comparison.json).
No waveform shifting or fitted gain is applied. Most discrepancies
approximately halve as ParaEMT's step is halved; GridKit's own refinement
differences are much smaller. This supports the parameter/equation
mapping for this balanced governor-step experiment.

The maximum voltage-magnitude discrepancy at 12.5 µs is `3.03192e-4 pu`,
at bus 3 at **50 microseconds**, during ParaEMT's startup transient. It does not
follow the same full-run halving trend. Its individual causes have not
been separated. ParaEMT has explicit controller updates, time-step-dependent
network damping, and fixed machine companion damping (`99/101`);
GridKit solves the continuous DAE with IDA. The original rounded initial
state also differs slightly. These are measured comparisons, not an
externally established acceptance tolerance or a reproduction of a
specific published figure.

These experiments do not validate unbalanced transformer winding behavior,
saturation, PLL dynamics, or active voltage cutouts. The trip is an
ideal-opening comparison with different post-trip G1 semantics; its
nonconverged ParaEMT voltage spike remains unsuitable as an answer key. General application support for consistent
inductive initialization, scheduled writable setpoints, and these
terminal constraints remains to be integrated outside this case driver.
