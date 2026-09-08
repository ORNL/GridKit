# GASTPTI Model

GASTPTI is a gas turbine-governor model with speed-droop fuel control and an
exhaust-temperature low-value selector.

## Notes

- PowerWorld caps its load reference at $A_T$ during transient simulation;
  GridKit does not.
- The GASTD-only `dbL`/`dbH` speed deadband is not part of GASTPTI.
- Unlike PowerWorld, GridKit rejects
  rather than swaps reversed $V^{\min}$ and $V^{\max}$ values.

> [!WARNING]
> GridKit does not yet apply the associated generator's Governor Response Limits
> modes `Down Only` and `Fixed` to GASTPTI. Normal response is always used.

## Block Diagram

![GASTPTI governor block diagram](../../../../../../docs/Figures/PhasorDynamics/GASTPTI/diagram.png)

Figure 1: GASTPTI governor model. Figure courtesy of the
[PowerWorld GAST_PTI model reference](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Governor%20GAST_PTI%20and%20GASTD.htm).

The PhasorDynamics GASTPTI equations, smooth low-value selector, valve
anti-windup, thermal feedback, and initialization are retained. EMT `speed`
is absolute rotor speed (one at synchronous), and `pmech`/`pref` use the
connected machine power base. Required `S` supplies that base in VA.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$S^{\mathrm{machine}}$ | [VA] | `S` | Connected-machine power base | Required
$R$ | [p.u.] | `R` | Permanent speed droop | Default 0.05; Speed deviation per component-base power deviation
$T_1$ | [s] | `T1` | Fuel-valve time constant | Default 0.4
$T_2$ | [s] | `T2` | Fuel-flow time constant | Default 0.1
$T_3$ | [s] | `T3` | Exhaust-temperature time constant | Default 3.0
$A_T$ | [p.u.] | `At` | Ambient-temperature load limit | Default 1.0; Component base
$K_T$ | [p.u.] | `Kt` | Exhaust-temperature feedback gain | Default 2.0
$V^{\max}$ | [p.u.] | `Vmax` | Upper valve response limit | Default 1.0; Component base
$V^{\min}$ | [p.u.] | `Vmin` | Lower valve response limit | Default 0.0; Component base
$D^\mathrm{turb}$ | [p.u.] | `Dturb` | Turbine damping coefficient | Default 0.0; Component-base power per speed deviation
$T^\mathrm{rate}$ | [MW] | `Trate` | Turbine rating | Default Machine base; Same-valued MVA component base when provided; GridKit addition

### Parameter Validation

```math
\begin{aligned}
  R &> 0 \\
  T_1,T_2,T_3 &\ge 0 \\
  A_T,K_T,D^\mathrm{turb} &\ge 0 \\
  T^\mathrm{rate} &> 0 \quad \text{when provided} \\
  V^{\min} &\le V^{\max}
\end{aligned}
```

### Derived Parameters

Let $\epsilon_T=10^{-3}\ \mathrm{s}$. Accepted time constants below
$\epsilon_T$ are raised to that floor in place:

```math
\begin{aligned}
  T_x &\leftarrow \max\!\left(T_x,\epsilon_T\right),
    && x\in\{1,2,3\} \\
  S^{\mathrm{base}}
    &\leftarrow
      \begin{cases}
        10^6 T^\mathrm{rate} & T^\mathrm{rate}\text{ provided} \\
        S^{\mathrm{machine}} & T^\mathrm{rate}\text{ omitted}
      \end{cases} \\
  k_{\mathrm{base}}
    &= \dfrac{S^{\mathrm{machine}}}{S^{\mathrm{base}}}
\end{aligned}
```

Multiplication by $k_{\mathrm{base}}$ converts machine-base power to component
base. $S^{\mathrm{machine}}$ and $S^{\mathrm{base}}$ are stored in VA.

`S` must be finite and positive. The optional turbine rating `Trate` remains
in MW; when omitted it uses the machine base `S`. In equations below, the
speed-deviation variable is the EMT input `speed` minus one.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\omega_r$ | `speed` | Input | [p.u.] | Machine rotor speed | Optional, defaults to one
$P^\mathrm{ref}$ | `pref` | Input | [p.u.] | Active-power reference | Machine base; inferred when unattached
$P_{\mathrm{m}}$ | `pmech` | Output | [p.u.] | Mechanical power | Machine base; seeded by the machine

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units  | Description                        | Note
-------|--------|------------------------------------|-----
$x_V$  | [p.u.] | Fuel-valve state                   | State 1 in Fig. 1; component base
$x_F$  | [p.u.] | Fuel-flow state                    | State 2 in Fig. 1; component base
$x_T$  | [p.u.] | Exhaust-temperature feedback state | State 3 in Fig. 1; component base

#### Algebraic

Symbol           | Units  | Description                   | Note
-----------------|--------|-------------------------------|-----
$V_D$            | [p.u.] | Speed/load fuel demand        | Component base
$V_T$            | [p.u.] | Temperature-limit fuel demand | Component base
$V$              | [p.u.] | Low-value selector output     | Component base
$P_{\mathrm{m}}$   | [p.u.] | Mechanical power output       | Machine base

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\omega_r$ | [p.u.] | Machine rotor speed | Optional `speed`; defaults to one
$P^\mathrm{ref}$ | [p.u.] | Active-power reference | Optional `pref`; machine base

## Model Equations

The speed deviation is $\omega = \omega_r - 1$.

### Internal Equations

#### Differential

The initialized constants $V_{\mathrm{resp}}^{\min}$,
$V_{\mathrm{resp}}^{\max}$, and $s^{\mathrm{valve}}$ are defined under
[Internal Initialization](#internal-initialization).

```math
\begin{aligned}
  0 &=
    -\dfrac{\mathrm{d}x_V}{\mathrm{d}t}
    + \dfrac{s^{\mathrm{valve}}}{T_1}
      \mathrm{antiwindup}\!\left(
        x_V,V-x_V;
        V_{\mathrm{resp}}^{\min},V_{\mathrm{resp}}^{\max}
      \right) \\
  0 &=
    -\dfrac{\mathrm{d}x_F}{\mathrm{d}t}
    + \dfrac{1}{T_2}\left(-x_F+x_V\right) \\
  0 &=
    -\dfrac{\mathrm{d}x_T}{\mathrm{d}t}
    + \dfrac{1}{T_3}\left(-x_T+x_F\right).
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
  0 &= - \omega + R(k_{\mathrm{base}}P^\mathrm{ref}-V_D) \\
  0 &= -V_T + A_T+K_T(A_T-x_T) \\
  0 &= -V + \min(V_D,V_T) \\
  0 &= -k_{\mathrm{base}}P_{\mathrm{m}} + x_F-D^\mathrm{turb}\omega.
\end{aligned}
```

CommonMath defines the [`antiwindup`](../../../../../CommonMath.md#antiwindup)
and [`min`](../../../../../CommonMath.md#minimum) targets and smooth approximations.

### External Equations

None.

## Initialization

### Input Initialization

```math
\begin{aligned}
  \omega &\leftarrow \text{machine speed deviation} \\
  P_{\mathrm{m}} &\leftarrow \text{machine mechanical power}
\end{aligned}
```

### Internal Initialization

```math
\begin{aligned}
  x_F
    &\leftarrow k_{\mathrm{base}}P_{\mathrm{m}}
       +D^\mathrm{turb}\omega \\
  x_V,x_T
    &\leftarrow x_F \\
  V_T
    &\leftarrow A_T+K_T\left(A_T-x_F\right) \\
  m_T
    &\leftarrow V_T-x_F \\
  \left(V_{\mathrm{resp}}^{\min},V_{\mathrm{resp}}^{\max}\right)
    &\leftarrow
      \left(\min(V^{\min},x_F),\max(V^{\max},x_F)\right)
\end{aligned}
```

```math
s^{\mathrm{valve}}
  \leftarrow
  \begin{cases}
    1 & V_{\mathrm{resp}}^{\min}\lt V_{\mathrm{resp}}^{\max} \\
    0 & V_{\mathrm{resp}}^{\min}=V_{\mathrm{resp}}^{\max}
  \end{cases}
```

```math
\left(V_{D},V\right)
  \leftarrow
  \begin{cases}
    \left(
      V_T-\mathrm{iramp}\!\left(m_T\right),
      x_F
    \right)
      & s^{\mathrm{valve}}=1 \\
    \left(
      x_F,
      \min\!\left(x_F,V_T\right)
    \right)
      & s^{\mathrm{valve}}=0
  \end{cases}
```

When $s^{\mathrm{valve}}=1$, initialization requires a finite positive
$m_T$ so `iramp` is defined. All candidates and response bounds are validated
before state, derivatives, or signals are changed; failed initialization is
atomic.

### Reference Initialization

For an unattached `pref` input, the controller latches the inferred reference:

```math
P^\mathrm{ref}
  \leftarrow
  \dfrac{1}{k_{\mathrm{base}}}
  \left(V_D+\dfrac{\omega}{R}\right).
```

Initialization uses the resolved machine-base $P_{\mathrm{m}}$ and preserves
any supplied `pref` input. If that input differs from the inferred reference,
consistent initialization resolves algebraic values and derivatives while
retaining the initialized differential states.

## Monitors

Monitor  | Units  | Description                        | Note
---------|--------|------------------------------------|-----
`pmech`  | [p.u.] | Mechanical-power output            | $P_{\mathrm{m}}$; machine base
`xvalve` | [p.u.] | Fuel-valve state                   | $x_V$; component base
`xflow`  | [p.u.] | Fuel-flow state                    | $x_F$; component base
`xtemp`  | [p.u.] | Exhaust-temperature feedback state | $x_T$; component base
`vload`  | [p.u.] | Speed/load fuel demand             | $V_D$; component base
`vtemp`  | [p.u.] | Temperature-limit fuel demand      | $V_T$; component base

## Appendix A: `iramp`

For a positive smooth-ramp output $v>0$ and CommonMath smoothing parameter
$\mu$,

```math
\mathrm{iramp}(v) = v+\dfrac{1}{\mu}\log\left(1-e^{-\mu v}\right).
```

This is the positive-range inverse of GridKit's smooth
[`ramp`](../../../../../CommonMath.md#ramp).
