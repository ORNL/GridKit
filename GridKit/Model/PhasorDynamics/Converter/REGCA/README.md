# **Renewable Energy Generator/Converter Model (REGCA)**

REGCA is a first-generation WECC renewable generator/converter model for
inverter-coupled resources.

## Notes

None.

## Block Diagram

![REGCA generator/converter block diagram](../../../../../docs/Figures/PhasorDynamics_REGCA_Diagram.png)

Figure 1: REGCA generator/converter model. Figure courtesy of the
[PowerWorld REGC_A model reference](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Machine%20Model%20REGC_A.htm).

## Model Parameters

Symbol                           | Units    | JSON     | Description                                           | Typical Value | Note
---------------------------------|----------|----------|-------------------------------------------------------|---------------|------
$P_0$                            | [p.u.]   | `p0`     | Initial active power injection                        | 1.0           | System base; required initialization source
$Q_0$                            | [p.u.]   | `q0`     | Initial reactive power injection                      | 0.0           | System base; required initialization source
$S^\mathrm{base}$                | [MVA]    | `mva`    | REGCA component power base                            | 100.0         |
$T_\mathrm{g}$                   | [sec]    | `Tg`     | Converter current-control lag time constant           | 0.02          | Block name: `Tg`
$T_M$                            | [sec]    | `TM`     | Terminal voltage sensor time constant                 | 0.02          | Block name: `Tfltr`
$R_q^{\max}$                     | [p.u./s] | `Rqmax`  | Reactive-current recovery positive rate limit         | 999.0         | Block name: `Iqrmax`; disabled when $R_q^{\max}\le 0$
$R_q^{\min}$                     | [p.u./s] | `Rqmin`  | Reactive-current recovery negative rate limit         | -999.0        | Block name: `Iqrmin`; disabled when $R_q^{\min}\ge 0$
$R_p^{\max}$                     | [p.u./s] | `Rpmax`  | Active-current magnitude recovery rate limit          | 999.0         | Block name: `rrpwr`; must be nonnegative
$s_L$                            | [binary] | `sL`     | LVPL switch                                           | 1             | Block name: `LPVLSW`
$I_{L1}$                         | [p.u.]   | `IL1`    | LVPL upper-current ceiling                            | 1.1           | Block name: `LVPL1`
$V_{L0}$                         | [p.u.]   | `VL0`    | LVPL zero-crossing voltage                            | 0.4           | Block name: `zerox`
$V_{L1}$                         | [p.u.]   | `VL1`    | LVPL upper breakpoint voltage                         | 0.9           | Block name: `brkpt`
$V_{A0}$                         | [p.u.]   | `VA0`    | LVACM lower breakpoint voltage                        | 0.4           | Block name: `LVPnt0`
$V_{A1}$                         | [p.u.]   | `VA1`    | LVACM upper breakpoint voltage                        | 0.9           | Block name: `LVPnt1`
$V_\mathrm{hv}^{\max}$           | [p.u.]   | `Vhvmax` | HV reactive management activation threshold           | 1.2           | Block name: `VLim`
$Q^{\min}$                       | [p.u.]   | `Qmin`   | PowerWorld compatibility field                         |               | Optional; accepted but unused
$K_\mathrm{hv}$                  | [p.u.]   | `Khv`    | HV reactive management gain                           | 0.7           | Optional; defaults to 0.7; block name: `Khv`
$X_\mathrm{e}$                   | [p.u.]   | `Xe`     | PowerWorld compatibility field                         |               | Optional; accepted but unused

All listed JSON parameters are required unless marked optional.

### Parameter Validation

Invalid REGCA parameter sets are rejected by the following checks. Let $\epsilon_T=10^{-3}$.
Time constants below $\epsilon_T$ are raised to $\epsilon_T$ and logged as a warning,
every other condition is a configuration error.

```math
\begin{aligned}
  T &\leftarrow \max(T, \epsilon_T)
    \quad T\in\{T_\mathrm{g},T_M\} \\
  S^\mathrm{base}
    &> 0 \\
  R_p^{\max}
    &\ge 0 \\
  I_{L1}
    &\ge 0 \\
  s_L
    &\in \{0,1\} \\
  0
    &\le V_{L0} < V_{L1} \\
  0
    &\le V_{A0} < V_{A1} < V_\mathrm{hv}^{\max} \\
  0 \le K_\mathrm{hv}
    &< \infty
\end{aligned}
```

### Model Derived Parameters

```math
\begin{aligned}
  s_L^\mathrm{off}
    &= 1 - s_L \\
  k_\mathrm{base}
    &= \dfrac{S^\mathrm{sys}}{S^\mathrm{base}} \\
  K_L
    &= 100
\end{aligned}
```

The fixed slope $K_L$ [p.u./p.u.] approximates the unbounded LVPL
characteristic above $V_{L1}$.

## Model Ports

Name       | Port   | Init    | Description
-----------|--------|---------|------
`bus`      | Bus    | Known   | Terminal bus voltage
`ipcmd`    | Input  | Unknown | Active-current command input
`iqcmd`    | Input  | Unknown | Reactive-current command input
`ibranchr` | Output | Known   | Branch-current real-component output
`ibranchi` | Output | Known   | Branch-current imaginary-component output
`pbranch`  | Output | Known   | Branch active-power output
`qbranch`  | Output | Known   | Branch reactive-power output

## Model Variables

### Internal Variables

#### Differential

Symbol                | Units  | Description               | Note
----------------------|--------|---------------------------|------
$V_M$                 | [p.u.] | Filtered terminal voltage | State 3 in Fig. 1
$I_q$                 | [p.u.] | Reactive-current state    | State 1 in Fig. 1 before the `-1` block; component base
$I_p$                 | [p.u.] | Active-current state      | State 2 in Fig. 1; component base

#### Algebraic

Symbol                     | Units    | Description                                                           | Note
---------------------------|----------|-----------------------------------------------------------------------|------
$V_T$                      | [p.u.]   | Terminal voltage magnitude                                            |
$I_\mathrm{r}$             | [p.u.]   | Branch-current real component                                         | System base
$I_\mathrm{i}$             | [p.u.]   | Branch-current imaginary component                                    | System base
$I_q^\mathrm{extra}$       | [p.u.]   | Extra inductive current from high-voltage reactive current management | Component base
$I_L$                      | [p.u.]   | LVPL upper-limit current curve                                        | Component base; function of $V_M$
$P^\mathrm{br}$            | [p.u.]   | Branch active power                                                   | System base
$Q^\mathrm{br}$            | [p.u.]   | Branch reactive power                                                 | System base

### External Variables

#### Differential
None.

#### Algebraic

Symbol                          | Units  | Init    | Description                                                      | Note
--------------------------------|--------|---------|------------------------------------------------------------------|------
$V_\mathrm{r}$                  | [p.u.] | Known   | Terminal voltage, real component                                 | Bus input
$V_\mathrm{i}$                  | [p.u.] | Known   | Terminal voltage, imaginary component                            | Bus input
$I_p^\mathrm{cmd}$              | [p.u.] | Unknown | Active-current command in the terminal-voltage reference frame   | Optional signal port `ipcmd`; system base
$I_q^\mathrm{cmd}$              | [p.u.] | Unknown | Reactive-current command in the terminal-voltage reference frame | Optional signal port `iqcmd`; system base

## Model Equations

Define the pre-limit current derivatives:

```math
\begin{aligned}
  f_\mathrm{q} &= \dfrac{1}{T_\mathrm{g}} (k_\mathrm{base} I_q^\mathrm{cmd} - I_q) \\
  f_\mathrm{p} &= \dfrac{1}{T_\mathrm{g}} (k_\mathrm{base} I_p^\mathrm{cmd} - I_p)
\end{aligned}
```

Figure 1 places LVPL on the active-current integrator ceiling, realized by
the `awmax` gate in the differential equations. The ceiling moves with the
sensed voltage, and a pinned $I_p$ tracks it as a non-windup limit. `rrpwr`
applies the active-current recovery rule according to the sign of $I_p$.

The limited active-current integrator drive applies the recovery rate rule
of [Appendix A](#appendix-a-rrpwr):

```math
f_\mathrm{p}^{\lim}
  = \text{rrpwr}(I_p, f_\mathrm{p}; R_p^{\max}).
```

### Differential Equations

The $I_q$ limiter branch is selected by the initial reactive power $Q_0$ and
the sign that enables the corresponding limit.

```math
\begin{aligned}
  0 &= -\dot V_M + \dfrac{1}{T_M} (V_T - V_M) \\
  0 &= -\dot I_q +
    \begin{cases}
      \text{min}(f_\mathrm{q}, R_q^{\max})
        & Q_0 > 0 \land R_q^{\max} > 0 \\
      \text{max}(f_\mathrm{q}, R_q^{\min})
        & Q_0 < 0 \land R_q^{\min} < 0 \\
      f_\mathrm{q} & \text{otherwise}
    \end{cases} \\
  0 &= -\dot I_p +
    \begin{cases}
      f_\mathrm{p}^{\lim} & s_L = 0 \\
      \text{awmax}(I_p, f_\mathrm{p}^{\lim}; I_L, \dot I_L) & s_L = 1
    \end{cases}
\end{aligned}
```


### Algebraic Equations

```math
\begin{aligned}
  0 &= -V_T^2 + V_\mathrm{r}^2 + V_\mathrm{i}^2 \\
  0 &= -k_\mathrm{base} V_T I_\mathrm{r}
       + V_\mathrm{i}(I_q - I_q^\mathrm{extra})
       + V_\mathrm{r} I_p\,\text{linseg}(V_T; V_{A0}, V_{A1}, 1) \\
  0 &= -k_\mathrm{base} V_T I_\mathrm{i}
       - V_\mathrm{r}(I_q - I_q^\mathrm{extra})
       + V_\mathrm{i} I_p\,\text{linseg}(V_T; V_{A0}, V_{A1}, 1) \\
  0 &= -I_q^\mathrm{extra}
       + K_\mathrm{hv}\,\text{ramp}(V_T - V_\mathrm{hv}^{\max}) \\
  0 &= -I_L
       + \text{linseg}(V_M; V_{L0}, V_{L1}, I_{L1})
       + K_L\,\text{ramp}(V_M - V_{L1}) \\
  0 &= -P^\mathrm{br}
       + V_\mathrm{r} I_\mathrm{r} + V_\mathrm{i} I_\mathrm{i} \\
  0 &= -Q^\mathrm{br}
       + V_\mathrm{i} I_\mathrm{r} - V_\mathrm{r} I_\mathrm{i}
\end{aligned}
```

CommonMath defines the [primitives](../../../../CommonMath.md#primitives) and
[derived functions](../../../../CommonMath.md#derived-functions) used above.

## Network Interface

```math
\begin{aligned}
  I_\mathrm{r}^\mathrm{inj} &:= I_\mathrm{r} \\
  I_\mathrm{i}^\mathrm{inj} &:= I_\mathrm{i}
\end{aligned}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
  V_\mathrm{r}, V_\mathrm{i}
    &\leftarrow \text{terminal-bus voltage} \\
  P_0, Q_0
    &\leftarrow \text{power-flow injection on system base}
\end{aligned}
```

### Internal Initialization

REGCA requires $V_{A1} \le V_{T,0}$, which excludes initialization below the
nominal upper LVACM breakpoint.

With LVPL enabled, REGCA additionally requires $I_{p,0} \le I_{L,0}$.
Initialization rejects an operating point above the active-current integrator
ceiling.

Subscript $0$ denotes initial values; all internal derivatives are initialized
to zero:

```math
\begin{aligned}
  V_{T,0}
    &= \sqrt{V_{\mathrm{r},0}^2 + V_{\mathrm{i},0}^2} \\
  V_{M,0}
    &= V_{T,0} \\
  A_0^\mathrm{LVACM}
    &= \text{linseg}(V_{T,0}; V_{A0}, V_{A1}, 1) \\
  I_{L,0}
    &= \text{linseg}(V_{T,0}; V_{L0}, V_{L1}, I_{L1})
       + K_L\,\text{ramp}(V_{T,0} - V_{L1}) \\
  I_{p,0}
    &= \dfrac{k_\mathrm{base}P_0}{V_{T,0}A_0^\mathrm{LVACM}} \\
  k_\mathrm{base} I_{p,0}^\mathrm{cmd}
    &= I_{p,0} \\
  I_{q,0}^\mathrm{extra}
    &= K_\mathrm{hv}\,\text{ramp}(V_{T,0} - V_\mathrm{hv}^{\max}) \\
  I_{q,0}^\mathrm{cmd}
    &= \dfrac{Q_0}{V_{T,0}}
       + \dfrac{I_{q,0}^\mathrm{extra}}{k_\mathrm{base}} \\
  I_{q,0}
    &= k_\mathrm{base} I_{q,0}^\mathrm{cmd}
\end{aligned}
```

The remaining algebraic quantities are then initialized as follows:

```math
\begin{aligned}
  I_{\mathrm{r},0}
    &= \dfrac{V_{\mathrm{r},0}P_0 + V_{\mathrm{i},0}Q_0}{V_{T,0}^2} \\
  I_{\mathrm{i},0}
    &= \dfrac{V_{\mathrm{i},0}P_0 - V_{\mathrm{r},0}Q_0}{V_{T,0}^2} \\
  P_0^\mathrm{br}
    &= P_0 \\
  Q_0^\mathrm{br}
    &= Q_0
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
  I_p^\mathrm{cmd}
    &\leftarrow I_{p,0}^\mathrm{cmd} \\
  I_q^\mathrm{cmd}
    &\leftarrow I_{q,0}^\mathrm{cmd}
\end{aligned}
```

## Monitorable Outputs

Output | Units  | Description                 | Note
-------|--------|-----------------------------|------
`ir`   | [p.u.] | Real current injection      | System base; exported through `ibranchr` when assigned
`ii`   | [p.u.] | Imaginary current injection | System base; exported through `ibranchi` when assigned
`p`    | [p.u.] | Active-power output         | System base; exported through `pbranch` when assigned
`q`    | [p.u.] | Reactive-power output       | System base; exported through `qbranch` when assigned

## Testing

- `validation()` checks construction, monitor creation, parameter validation, bus and signal configuration, and minimum time-constant handling.
- `initializationAndSignals()` checks power-flow initialization, base conversion, signal publication, bus injection, and unattached-command latching.
- `initializationDomain()` checks rejected and accepted voltage and LVPL initialization boundaries.
- `residualEquations()` checks every model residual at a hand-computable midpoint state.
- `activeCurrentControl()` checks `rrpwr`, enabled and bypassed LVPL behavior, and tracking of a moving LVPL ceiling.
- `reactiveCurrentControl()` checks the positive, negative, and unrestricted reactive-current recovery-rate branches.
- `highVoltageManagement()` checks HVRCM initialization, $K_\mathrm{hv}$ loading, threshold behavior, and its local derivative.
- `jacobian()` compares the dependency-tracking and Enzyme Jacobians for enabled and bypassed LVPL configurations when Enzyme support is enabled.

Tests use model identities rather than frozen smoothing decimals and reserve
$100 \epsilon$ for roundoff checks.

## Appendix A: `rrpwr`

The exact active-current rate-limit rule is

```math
\text{rrpwr}(x, f; r) =
  \begin{cases}
    \text{max}(f, -r) & x < 0 \\
    \text{clamp}(f; -r, r) & x = 0 \\
    \text{min}(f, r) & x > 0
  \end{cases}
```

The model evaluates this rule with the following continuously differentiable
($C^1$) approximation:

```math
\begin{aligned}
  t(x) &= \tanh\!\left(\dfrac{\mu x}{2}\right)=2\sigma(x)-1 \\[0pt]
  w_+(x) &= \dfrac{t(x)^2+t(x)\lvert t(x)\rvert}{2} \\[0pt]
  w_-(x) &= \dfrac{t(x)^2-t(x)\lvert t(x)\rvert}{2} \\[0pt]
  \text{rrpwr}(x,f;r)
    &\approx f
      +\left[1-w_+(x)\right]\text{ramp}(-f-r)
      -\left[1-w_-(x)\right]\text{ramp}(f-r).
\end{aligned}
```

The one-sided weights and their first derivatives vanish at $x=0$. The
approximation therefore equals `slew` exactly at zero and preserves the
outward rate limit for finite $\mu$ while gradually releasing restoring motion.

## Appendix B: `awmax`

The exact anti-windup rule under a moving upper bound $u$ with rate $\dot u$
is

```math
\text{awmax}(x, f; u, \dot u) =
  \begin{cases}
    f & x < u \\
    \text{min}(f, \dot u) & x \ge u
  \end{cases}
```

Below the bound the unconstrained derivative passes. Pinned at the bound, the
state tracks $\text{min}(f, \dot u)$, so a falling bound drags the state down
with it. The rule is fixed-bound anti-windup on the gap $g = x - u$ held below
zero: pinned, $\dot g = \text{min}(f - \dot u, 0) \le 0$, so the gap cannot
grow and closes whenever $f < \dot u$.

The model evaluates this rule with the following smooth approximation:

```math
\text{awmax}(x, f; u, \dot u)
  \approx \dot u
    + \left[\sigma(u-x)+\left(1-\sigma(u-x)\right)\sigma(\dot u-f)\right]
      (f - \dot u)
```

With a stationary bound ($\dot u = 0$) this reduces to the `antiwindup` of
CommonMath restricted to its upper limit, which admits an algebraic bound $u$.
