# **Gas Turbine-Governor Model (GASTPTI)**

GASTPTI is a gas turbine-governor model for thermal generating units with a
speed-droop fuel valve and an exhaust-temperature load limit.

## Notes

- Power signals and the `pmech` monitor output are on system base.
- Internal fuel-demand, fuel-flow, and temperature quantities are on GASTPTI
  component base.
- GASTPTI uses $T^\mathrm{rate}$, loaded from `Trate`, as its component power base.
- `mode = 2` (Down Only) is accepted with a warning and currently simulated as
  `mode = 0` (Normal).
- The GASTD-only speed deadband (`dbL`/`dbH`) shown in Figure 1 is not part of
  GASTPTI.
- Narrow valve limits attenuate the smooth anti-windup gate in Normal mode;
  verification warns when the maximum gate falls below 0.99.

## Block Diagram

![GASTPTI governor block diagram](../../../../../docs/Figures/PhasorDynamics/GASTPTI/diagram.png)

Figure 1: GASTPTI governor model. Figure courtesy of the
[PowerWorld GAST_PTI model reference](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Governor%20GAST_PTI%20and%20GASTD.htm).

## Model Parameters

Symbol             | Units     | JSON    | Description                            | Typical Value | Note
-------------------|-----------|---------|----------------------------------------|---------------|------
$R$                | [p.u.]    | `R`     | Permanent droop                        | 0.05          |
$T_1$              | [sec]     | `T1`    | Fuel-valve time constant               | 0.4           | State 1; raised to the minimum-time floor
$T_2$              | [sec]     | `T2`    | Fuel-flow time constant                | 0.1           | State 2; raised to the minimum-time floor
$T_3$              | [sec]     | `T3`    | Exhaust-temperature time constant      | 3.0           | State 3; raised to the minimum-time floor
$A_T$              | [p.u.]    | `At`    | Ambient-temperature load limit         | 1.0           |
$K_T$              | [p.u.]    | `Kt`    | Exhaust-temperature feedback gain      | 2.0           |
$V^{\max}$         | [p.u.]    | `Vmax`  | Maximum fuel-valve/turbine-power limit | 1.0           |
$V^{\min}$         | [p.u.]    | `Vmin`  | Minimum fuel-valve/turbine-power limit | 0.0           |
$D^\mathrm{turb}$  | [p.u.]    | `Dturb` | Turbine damping coefficient            | 0.0           | Multiplied by speed deviation
$T^\mathrm{rate}$  | [MW]      | `Trate` | Turbine-rating power base              | 100.0         |
$\mathrm{mode}$    | [integer] | `mode`  | Governor response mode                 | 0             | Source label: `TSGovRespLimit`; 0 = Normal, 1 = Fixed, 2 = Down Only

### Parameter Validation

Invalid GASTPTI parameter sets are rejected by the following checks:

```math
\begin{aligned}
  R &> 0 \\
  T_1, T_2, T_3
    &\ge 0 \\
  A_T
    &\ge 0 \\
  K_T
    &\ge 0 \\
  V^{\min}
    &\le V^{\max}
    \quad \mathrm{mode}=1 \\
  V^{\min}
    &< V^{\max}
    \quad \mathrm{mode}\in\{0,2\} \\
  D^\mathrm{turb}
    &\ge 0 \\
  T^\mathrm{rate}
    &> 0 \\
  \mathrm{mode}
    &\in \{0,1,2\}
\end{aligned}
```

### Model Derived Parameters

Let $\epsilon_T=10^{-3}\ \mathrm{s}$. A time constant below $\epsilon_T$ is
raised to that floor in place, so every equation below uses the raised value:

```math
\begin{aligned}
  T_x
    &\leftarrow \max\!\left(T_x,\epsilon_T\right),
       \quad x\in\{1,2,3\} \\
  k_{\mathrm{base}}
    &= \dfrac{S^\mathrm{sys}}{T^\mathrm{rate}} \\
  s^\mathrm{fix}
    &=
    \begin{cases}
      0 & \mathrm{mode}=1 \\
      1 & \mathrm{mode}\in\{0,2\}
    \end{cases}
\end{aligned}
```

Multiplying by $k_\mathrm{base}$ converts system base to component base.

## Model Ports

Name    | Port   | Init    | Description
--------|--------|---------|------
`speed` | Input  | Known   | Machine speed deviation
`pref`  | Input  | Unknown | Active-power/load reference
`pmech` | Output | Known   | Mechanical power output

`Known` ports are seeded before `initialize()` and preserved by it. `Unknown`
inputs are resolved during initialization and written to attached signal
storage, or retained as constant inputs when the port is unattached.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units  | Description                        | Note
-------|--------|------------------------------------|------
$x_V$  | [p.u.] | Fuel-valve state                   | State 1 in Fig. 1
$x_F$  | [p.u.] | Fuel-flow state                    | State 2 in Fig. 1
$x_T$  | [p.u.] | Exhaust-temperature feedback state | State 3 in Fig. 1

#### Algebraic

Symbol           | Units  | Description                   | Note
-----------------|--------|-------------------------------|------
$V_D$            | [p.u.] | Speed/load fuel demand        | LV gate input
$V_T$            | [p.u.] | Temperature-limit fuel demand | LV gate input
$V$              | [p.u.] | LV gate output                | Smooth minimum of $V_D$ and $V_T$
$P_{\mathrm{m}}$ | [p.u.] | Mechanical power to generator | System base; assigned to `pmech`

### External Variables

#### Differential

None.

#### Algebraic

Symbol           | Units  | Init    | Description                 | Note
-----------------|--------|---------|-----------------------------|------
$\omega$         | [p.u.] | Known   | Machine speed deviation     | Optional signal port `speed`; defaults to zero
$P^\mathrm{ref}$ | [p.u.] | Unknown | Active-power/load reference | Optional signal port `pref`; system base

## Model Equations

### Differential Equations

```math
\begin{aligned}
  0 &=
    -\dot{x}_V
    + \dfrac{s^\mathrm{fix}}{T_1}\,
      \text{antiwindup}
      \left(x_V, V - x_V;\, V^{\min}, V^{\max}\right) \\
  0 &=
    -\dot{x}_F
    + \dfrac{s^\mathrm{fix}}{T_2}
      \left(-x_F + x_V\right) \\
  0 &=
    -\dot{x}_T
    + \dfrac{s^\mathrm{fix}}{T_3}
      \left(-x_T + x_F\right)
\end{aligned}
```

CommonMath defines the [`antiwindup`](../../../../CommonMath.md#antiwindup)
target and smooth approximation.

### Algebraic Equations

```math
\begin{aligned}
  0 &=
    -\omega
    + R\left(k_{\mathrm{base}}P^\mathrm{ref} - V_D\right) \\
  0 &=
    -V_T
    + A_T
    + K_T\left(A_T - x_T\right) \\
  0 &=
    -V
    + \text{min}\left(V_D, V_T\right) \\
  0 &=
    -k_{\mathrm{base}}P_{\mathrm{m}}
    + x_F
    - D^\mathrm{turb}\omega
\end{aligned}
```

CommonMath defines helper targets and smooth approximations for
[min](../../../../CommonMath.md#derived-functions).

## Initialization

### Input Initialization

```math
\begin{aligned}
  \omega
    &\leftarrow \text{machine speed deviation} \\
  P_{\mathrm{m}}
    &\leftarrow \text{machine mechanical-power seed on system base}
\end{aligned}
```

Initialization never replaces the system-base value held in $P_{\mathrm{m}}$.

### Internal Initialization

Subscript $0$ denotes initial values; all internal derivatives start at zero.
$\rho^{-1}$ is the exact inverse of the smooth CommonMath
[`ramp`](../../../../CommonMath.md#primitives) (`iramp`), so the load demand
seats behind the smooth LV gate exactly:

```math
\begin{aligned}
  x_{F,0}
    &= k_{\mathrm{base}}P_{\mathrm{m},0} + D^\mathrm{turb}\omega_0 \\
  x_{V,0}
    &= x_{F,0} \\
  x_{T,0}
    &= x_{F,0} \\
  V_0
    &= x_{F,0} \\
  V_{T,0}
    &= A_T + K_T\left(A_T - x_{T,0}\right) \\
  V_{D,0}
    &= V_{T,0} - \rho^{-1}\!\left(V_{T,0} - x_{F,0}\right)
\end{aligned}
```

Initialization rejects an operating point when the temperature-gate margin
$V_{T,0} - x_{F,0}$ is not positive. A non-Fixed-mode fuel flow outside
$[V^{\min}, V^{\max}]$ produces a warning and initializes at the dispatched
value.

Every check resolves before any storage is written, so a rejected
initialization leaves state, the `pmech` seed, and external signals unchanged.

### Output Initialization

```math
\begin{aligned}
  P^\mathrm{ref}
    &\leftarrow
      \dfrac{1}{k_{\mathrm{base}}}
      \left(V_{D,0} + \dfrac{\omega_0}{R}\right)
\end{aligned}
```

GASTPTI writes the resolved active-power/load reference to an attached `pref`
signal input. If no controller is connected, that value is used as a constant
reference input.

## Monitorable Outputs

Output   | Units  | Description                        | Note
---------|--------|------------------------------------|------
`pmech`  | [p.u.] | Mechanical-power output            | $P_{\mathrm{m}}$ (system base)
`xvalve` | [p.u.] | Fuel-valve state                   | $x_V$ (component base)
`xflow`  | [p.u.] | Fuel-flow state                    | $x_F$ (component base)
`xtemp`  | [p.u.] | Exhaust-temperature feedback state | $x_T$ (component base)
`vload`  | [p.u.] | Speed/load fuel demand             | $V_D$ (component base)
`vtemp`  | [p.u.] | Temperature-limit fuel demand      | $V_T$ (component base)
