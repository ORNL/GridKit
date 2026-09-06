# GASTPTI

GASTPTI is a gas turbine-governor model with speed-droop fuel control and an
exhaust-temperature low-value selector.

> [!WARNING]
> GridKit does not yet apply the associated generator's Governor Response Limits
> modes `Down Only` and `Fixed` to GASTPTI. Normal response is always used.

## Notes

- PowerWorld caps its load reference at $A_T$ during transient simulation;
  GridKit does not.
- The GASTD-only `dbL`/`dbH` speed deadband is not part of GASTPTI.
- Unlike PowerWorld, GridKit rejects
  rather than swaps reversed $V^{\min}$ and $V^{\max}$ values.

## Block Diagram

![GASTPTI governor block diagram](../../../../../docs/Figures/PhasorDynamics/GASTPTI/diagram.png)

Figure 1: GASTPTI governor model. Figure courtesy of the
[PowerWorld GAST_PTI model reference](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Governor%20GAST_PTI%20and%20GASTD.htm).

## Model Parameters

Symbol            | Units  | JSON    | Description                       | Typical Value | Note
------------------|--------|---------|-----------------------------------|---------------|---------------------------------------------------------------
$R$               | [p.u.] | `R`     | Permanent speed droop             | 0.05          | Speed deviation per component-base power deviation
$T_1$             | [s]    | `T1`    | Fuel-valve time constant          | 0.4           |
$T_2$             | [s]    | `T2`    | Fuel-flow time constant           | 0.1           |
$T_3$             | [s]    | `T3`    | Exhaust-temperature time constant | 3.0           |
$A_T$             | [p.u.] | `At`    | Ambient-temperature load limit    | 1.0           | Component base
$K_T$             | [p.u.] | `Kt`    | Exhaust-temperature feedback gain | 2.0           |
$V^{\max}$        | [p.u.] | `Vmax`  | Upper valve response limit        | 1.0           | Component base
$V^{\min}$        | [p.u.] | `Vmin`  | Lower valve response limit        | 0.0           | Component base
$D^\mathrm{turb}$ | [p.u.] | `Dturb` | Turbine damping coefficient       | 0.0           | Component-base power per speed deviation
$T^\mathrm{rate}$ | [MW]   | `Trate` | Turbine rating                    | Required      | Same-valued MVA component base; GridKit addition

### Parameter Validation

A valid GASTPTI parameter set must satisfy the following conditions:

```math
\begin{aligned}
  R &> 0 \\
  T_1,T_2,T_3 &\ge 0 \\
  A_T,K_T,D^\mathrm{turb} &\ge 0 \\
  T^\mathrm{rate} &> 0 \\
  V^{\min} &\le V^{\max}
\end{aligned}
```

### Model Derived Parameters

Let $\epsilon_T=10^{-3}\ \mathrm{s}$. Accepted time constants below
$\epsilon_T$ are raised to that floor in place:

```math
\begin{aligned}
  T_x &\leftarrow \max\!(T_x,\epsilon_T),
    && x\in\{1,2,3\} \\
  S^{\mathrm{base}}
    &\leftarrow 10^6 T^\mathrm{rate} \\
  k_{\mathrm{base}}
    &= \dfrac{S^{\mathrm{sys}}}{S^{\mathrm{base}}}
\end{aligned}
```

Multiplication by $k_{\mathrm{base}}$ converts system-base power to component
base. $S^{\mathrm{sys}}$ and $S^{\mathrm{base}}$ are stored in VA.

## Model Ports

Name    | Port   | Init    | Description
--------|--------|---------|------------
`speed` | Input  | Known   | Optional machine speed deviation; defaults to zero
`pref`  | Input  | Unknown | Optional active-power/load reference; latches its initialized value when unattached
`pmech` | Output | Known   | Required mechanical power output

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
$P_\mathrm{m}$   | [p.u.] | Mechanical power output       | System base

### External Variables

#### Differential

None.

#### Algebraic

Symbol           | Units  | Description                 | Note
-----------------|--------|-----------------------------|-----------------------------------
$\omega$         | [p.u.] | Machine speed deviation     | Optional `speed`; defaults to zero
$P^\mathrm{ref}$ | [p.u.] | Active-power/load reference | Optional `pref`; system base

## Model Equations

Smooth functions: [`antiwindup`](../../../../CommonMath.md#antiwindup), [`min`](../../../../CommonMath.md#minimum).

### Internal Equations

#### Differential

The initialized constants $V_{\mathrm{resp}}^{\min}$,
$V_{\mathrm{resp}}^{\max}$, and $s^{\mathrm{valve}}$ are defined under
[Internal Initialization](#internal-initialization).

```math
\begin{aligned}
  0 &=
    -\dot{x}_V
    + \dfrac{s^{\mathrm{valve}}}{T_1}
      \text{antiwindup}\!(
        x_V,V-x_V;
        V_{\mathrm{resp}}^{\min},V_{\mathrm{resp}}^{\max}
      ) \\
  0 &=
    -\dot{x}_F
    + \dfrac{1}{T_2}(-x_F+x_V) \\
  0 &=
    -\dot{x}_T
    + \dfrac{1}{T_3}(-x_T+x_F)
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
  0 &= - \omega + R(k_{\mathrm{base}}P^\mathrm{ref}-V_D) \\
  0 &= -V_T + A_T+K_T(A_T-x_T) \\
  0 &= -V + \text{min}(V_D,V_T) \\
  0 &= -k_{\mathrm{base}}P_\mathrm{m} + x_F-D^\mathrm{turb}\omega
\end{aligned}
```

### External Equations

None.

## Initialization

### Input Initialization

```math
\begin{aligned}
  \omega &\leftarrow \text{machine speed deviation} \\
  P_\mathrm{m} &\leftarrow \text{machine mechanical power}
\end{aligned}
```

### Internal Initialization

```math
\begin{aligned}
  x_F
    &\leftarrow k_{\mathrm{base}}P_\mathrm{m}
       +D^\mathrm{turb}\omega \\
  x_V,x_T
    &\leftarrow x_F \\
  V_T
    &\leftarrow A_T+K_T(A_T-x_F) \\
  m_T
    &\leftarrow V_T-x_F \\
  (V_{\mathrm{resp}}^{\min},V_{\mathrm{resp}}^{\max})
    &\leftarrow
      (\min(V^{\min},x_F),\max(V^{\max},x_F))
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
(V_{D},V)
  \leftarrow
  \begin{cases}
    (
      V_T-\text{iramp}\!(m_T),
      x_F,
    )
      & s^{\mathrm{valve}}=1 \\
    (
      x_F,
      \text{min}\!(x_F,V_T)
    )
      & s^{\mathrm{valve}}=0
  \end{cases}
```

When $s^{\mathrm{valve}}=1$, initialization requires a finite positive
$m_T$ so `iramp` is defined. All candidates and response bounds are validated
before state, derivatives, or signals are changed; failed initialization is
atomic.

A response limit that excludes the initialized fuel flow is widened to include
it, and a warning is logged. This matches PowerWorld's default
`Modify Limits and Run` treatment of initial limit violations.

### Output Initialization

```math
P^\mathrm{ref}
  \leftarrow
  \dfrac{1}{k_{\mathrm{base}}}
  \left(V_D+\dfrac{\omega}{R}\right)
```

Initialization preserves the machine-seeded system-base $P_\mathrm{m}$. An
attached `pref` signal receives the initialized reference; an unattached port
latches that value for subsequent residual evaluations.

## Monitors

Monitor  | Units  | Description                        | Note
---------|--------|------------------------------------|-----
`pmech`  | [p.u.] | Mechanical-power output            | $P_\mathrm{m}$; system base
`xvalve` | [p.u.] | Fuel-valve state                   | $x_V$; component base
`xflow`  | [p.u.] | Fuel-flow state                    | $x_F$; component base
`xtemp`  | [p.u.] | Exhaust-temperature feedback state | $x_T$; component base
`vload`  | [p.u.] | Speed/load fuel demand             | $V_D$; component base
`vtemp`  | [p.u.] | Temperature-limit fuel demand      | $V_T$; component base

## Appendix A: `iramp`

For a positive smooth-ramp output $v>0$ and CommonMath smoothing parameter
$\mu$,

```math
\text{iramp}(v) = v+\dfrac{1}{\mu}\log(1-\exp(-\mu v))
```

This is the positive-range inverse of GridKit's smooth
[`ramp`](../../../../CommonMath.md#ramp).
