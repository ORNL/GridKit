# **Gas Turbine-Governor Model (GASTPTI)**

GASTPTI is a gas turbine-governor model with speed-droop fuel control and an
exhaust-temperature low-value selector.

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

Symbol            | Units     | JSON    | Description                           | Typical Value | Note
------------------|-----------|---------|---------------------------------------|---------------|-----
$R$               | [p.u.]    | `R`     | Permanent speed droop                 | 0.05          | Speed deviation per component-base power deviation
$T_1$             | [sec]     | `T1`    | Fuel-valve time constant              | 0.4           |
$T_2$             | [sec]     | `T2`    | Fuel-flow time constant               | 0.1           |
$T_3$             | [sec]     | `T3`    | Exhaust-temperature time constant     | 3.0           |
$A_T$             | [p.u.]    | `At`    | Ambient-temperature load limit        | 1.0           | Component base
$K_T$             | [p.u.]    | `Kt`    | Exhaust-temperature feedback gain     | 2.0           |
$V^{\max}$        | [p.u.]    | `Vmax`  | Upper valve response limit            | 1.0           | Component base
$V^{\min}$        | [p.u.]    | `Vmin`  | Lower valve response limit            | 0.0           | Component base
$D^\mathrm{turb}$ | [p.u.]  | `Dturb` | Turbine damping coefficient           | 0.0           | Component-base power per speed deviation
$T^\mathrm{rate}$ | [MW]   | `Trate` | Turbine rating                        | 100.0         | Same-valued MVA component base; GridKit addition
$\mathrm{mode}$   | [integer] | `mode`  | Governor response-limit mode          | 0             | 0 = Normal, 1 = Down Only, 2 = Fixed

### Parameter Validation

```math
\begin{aligned}
  R &> 0 \\
  T_1,T_2,T_3 &\ge 0 \\
  A_T,K_T,D^\mathrm{turb} &\ge 0 \\
  T^\mathrm{rate} &> 0 \\
  V^{\min} &\le V^{\max} \\
  \mathrm{mode} &\in\{0,1,2\}
\end{aligned}
```

The system and component power bases and both reciprocal conversion ratios must
also be finite and positive; port requirements are defined under Model Ports.

### Model Derived Parameters

Let $\epsilon_T=10^{-3}\ \mathrm{s}$. Accepted time constants below
$\epsilon_T$ are raised to that floor in place:

```math
\begin{aligned}
  T_x &\leftarrow \max\!\left(T_x,\epsilon_T\right),
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
`speed` | Input  | Known   | Machine speed deviation
`pref`  | Input  | Unknown | Active-power/load reference
`pmech` | Output | Known   | Mechanical power output

The `pmech` output must be assigned. Inputs are optional, but attached `speed`
and `pref` ports must be linked. Unattached `speed` is zero; initialization
publishes the resolved reference through attached `pref` or latches it otherwise.
Indexed `speed`, `pref`, and `pmech` ports must be distinct.

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
$P_{\text{m}}$   | [p.u.] | Mechanical power output       | System base

### External Variables

#### Differential

None.

#### Algebraic

Symbol           | Units  | Init    | Description                 | Note
-----------------|--------|---------|-----------------------------|-----
$\omega$         | [p.u.] | Known   | Machine speed deviation     | Optional `speed`; defaults to zero
$P^\mathrm{ref}$ | [p.u.] | Unknown | Active-power/load reference | Optional `pref`; system base

## Model Equations

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
      \text{antiwindup}\!\left(
        x_V,V-x_V;
        V_{\mathrm{resp}}^{\min},V_{\mathrm{resp}}^{\max}
      \right) \\
  0 &=
    -\dot{x}_F
    + \dfrac{1}{T_2}\left(-x_F+x_V\right) \\
  0 &=
    -\dot{x}_T
    + \dfrac{1}{T_3}\left(-x_T+x_F\right).
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
  0 &= - \omega + R(k_\mathrm{base}P^\mathrm{ref}-V_D) \\
  0 &= -V_T + A_T+K_T(A_T-x_T) \\
  0 &= -V + \text{min}(V_D,V_T) \\
  0 &= -k_\mathrm{base}P_\text{m} + x_F-D^\mathrm{turb}\omega.
\end{aligned}
```

CommonMath defines the [`antiwindup`](../../../../CommonMath.md#antiwindup)
and [`min`](../../../../CommonMath.md#min) targets and smooth approximations.

### External Equations

#### Differential

None.

#### Algebraic

None.

## Initialization

### Input Initialization

```math
\begin{aligned}
  \omega &\leftarrow \text{machine speed deviation} \\
  P_\text{m} &\leftarrow \text{machine mechanical power}
\end{aligned}
```

### Internal Initialization

```math
\begin{aligned}
  x_{F}
    &\leftarrow k_{\mathrm{base}}P_{\text{m}}
       +D^\mathrm{turb}\omega \\
  x_{V},x_{T}
    &\leftarrow x_{F} \\
  V_{T}
    &\leftarrow A_T+K_T\left(A_T-x_{F}\right) \\
  m_{T}
    &\leftarrow V_{T}-x_{F} \\
  \left(V_{\mathrm{resp}}^{\min},V_{\mathrm{resp}}^{\max}\right)
    &\leftarrow
      \begin{cases}
        \left(\min(V^{\min},x_{F}),\max(V^{\max},x_{F})\right)
          & \mathrm{mode}=0 \\
        \left(\min(V^{\min},x_{F}),x_{F}\right)
          & \mathrm{mode}=1 \\
        \left(x_{F},x_{F}\right)
          & \mathrm{mode}=2
      \end{cases} \\
  s^{\mathrm{valve}}
    &\leftarrow
      \begin{cases}
        1 & V_{\mathrm{resp}}^{\min}<V_{\mathrm{resp}}^{\max} \\
        0 & V_{\mathrm{resp}}^{\min}=V_{\mathrm{resp}}^{\max}
      \end{cases} \\
  \left(V_{D},V\right)
    &\leftarrow
      \begin{cases}
        \left(
          V_{T}-\text{iramp}\!\left(m_{T}\right),
          x_{F}
        \right)
          & s^{\mathrm{valve}}=1 \\
        \left(
          x_{F},
          \text{min}\!\left(x_{F},V_{T}\right)
        \right)
          & s^{\mathrm{valve}}=0
      \end{cases}
\end{aligned}
```

### Output Initialization

```math
P^\mathrm{ref}
  \leftarrow
  \dfrac{1}{k_{\mathrm{base}}}
  \left(V_{D,0}+\dfrac{\omega}{R}\right).
```

## Monitorable Outputs

Output   | Units  | Description                        | Note
---------|--------|------------------------------------|-----
`pmech`  | [p.u.] | Mechanical-power output            | $P_{\text{m}}$; system base
`xvalve` | [p.u.] | Fuel-valve state                   | $x_V$; component base
`xflow`  | [p.u.] | Fuel-flow state                    | $x_F$; component base
`xtemp`  | [p.u.] | Exhaust-temperature feedback state | $x_T$; component base
`vload`  | [p.u.] | Speed/load fuel demand             | $V_D$; component base
`vtemp`  | [p.u.] | Temperature-limit fuel demand      | $V_T$; component base

## Testing

- `validation()` checks defaults, parameter domains, response modes, signal
  configuration, and time-constant floors.
- `initializationAndSignals()` checks base conversion, signal initialization,
  monitor values, and unattached-reference latching.
- `initializationDomain()` checks accepted and rejected operating points for
  each response mode.
- `initializationExactness()` checks the smooth-selector inverse.
- `residualEquations()` checks every residual against a fixed answer key.
- `governorControl()` checks droop, damping, response limits, and anti-windup.
- `temperatureLimiting()` checks the low-value selector.
- `jacobian()` checks fixed DependencyTracking oracles and Enzyme agreement to
  $10^{-9}$ when enabled.

## Appendix A: `iramp`

For a positive smooth-ramp output $v>0$,

```math
\text{iramp}(v) = v+\dfrac{1}{\mu}\log\left(1-e^{-\mu v}\right).
```

This is the positive-range inverse of GridKit's smooth
[`ramp`](../../../../CommonMath.md#rho-ramp).
