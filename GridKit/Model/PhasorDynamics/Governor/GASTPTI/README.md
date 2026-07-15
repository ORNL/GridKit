# **Gas Turbine-Governor Model (GASTPTI)**

GASTPTI is a gas turbine-governor model for thermal generating units.

## Notes

> [!NOTE]
> `mode = 2` (Down Only) is accepted with a warning and currently simulated as
> `mode = 0` (Normal).

- The GASTD-only speed deadband (`dbL`/`dbH`) shown in Figure 1 is not part of
  GASTPTI.

## Block Diagram

![](../../../../../docs/Figures/PhasorDynamics/GASTPTI/diagram.png)

Figure 1: GASTPTI block diagram. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol                          | Units    | JSON     | Description                                  | Typical Value | Note
--------------------------------|----------|----------|----------------------------------------------|---------------|------
$R$                             | [p.u.]   | `R`      | Permanent droop                              | 0.05          |
$T_1$                           | [sec]    | `T1`     | Fuel-valve time constant                     | 0.4           |
$T_2$                           | [sec]    | `T2`     | Fuel-flow time constant                      | 0.1           |
$T_3$                           | [sec]    | `T3`     | Exhaust-temperature time constant            | 3.0           |
$A_T$                           | [p.u.]   | `At`     | Ambient-temperature load limit               | 1.0           |
$K_T$                           | [p.u.]   | `Kt`     | Exhaust-temperature feedback gain            | 2.0           |
$V^{\max}$                      | [p.u.]   | `Vmax`   | Maximum fuel-valve/turbine-power limit       | 1.0           |
$V^{\min}$                      | [p.u.]   | `Vmin`   | Minimum fuel-valve/turbine-power limit       | 0.0           |
$D^\mathrm{turb}$               | [p.u.]   | `Dturb`  | Turbine damping coefficient                  | 0.0           |
$T^\mathrm{rate}$               | [MW]     | `Trate`  | Turbine-rating power base                    | 100.0         |
$\mathrm{mode}$                  | [-]      | `mode`   | Governor response mode                       | 0             | PowerWorld: `TSGovRespLimit`; 0: Normal, 1: Fixed, 2: Down Only

### Parameter Validation

> [!NOTE]
> The current $T \leftarrow \max(T,10^{-3})$ treatment for
> $T\in\{T_1,T_2,T_3\}$ and the governor `mode` behaviour should be replaced by a structural template-based
> approach in the near future.

> [!WARNING]
> Until $\mu$ is configurable, narrow limits can interfere with the smooth anti-windup gate below in Normal mode (windowed limits), [TGOV1](../Tgov1/README.md),
> [IEEET1](../../Exciter/IEEET1/README.md), and
> [SEXS-PTI](../../Exciter/SEXS-PTI/README.md).
GASTPTI applies the following parameter rules with $\epsilon_T=10^{-3}$:

```math
\begin{aligned}
  T &\leftarrow \max\!\left(T, \epsilon_T\right)
    \quad T\in\{T_1,T_2,T_3\} \\
  R, T^\mathrm{rate} &> 0 \\
  A_T, K_T, D^\mathrm{turb} &\ge 0 \\
  \mathrm{mode} &\in \{0,1,2\} \\
  V^{\min} &\le V^{\max}
    \quad \mathrm{mode}=1 \\
  V^{\min} &< V^{\max}
    \quad \mathrm{mode}\in\{0,2\}
\end{aligned}
```

### Model Derived Parameters

```math
\begin{aligned}
  k^\mathrm{base}
    &= \dfrac{S^\mathrm{sys}}{T^\mathrm{rate}} \\
  s^\text{fix}
    &=
    \begin{cases}
      0 & \mathrm{mode}=1 \\
      1 & \mathrm{mode}\in\{0,2\}
    \end{cases}
\end{aligned}
```

## Model Ports

Name    | Port   | Init    | Description
--------|--------|---------|------
`speed` | Input  | Known   | Machine speed deviation
`pref`  | Input  | Unknown | Active-power/load reference
`pmech` | Output | Unknown | Mechanical power output

## Model Variables

### Internal Variables

#### Differential

Symbol                  | Units  | Description                        | Note
------------------------|--------|------------------------------------|------
$x_V$                   | [p.u.] | Fuel-valve state                   | State 1 in Fig. 1
$x_F$                   | [p.u.] | Fuel-flow state                    | State 2 in Fig. 1
$x_T$                   | [p.u.] | Exhaust-temperature feedback state | State 3 in Fig. 1

#### Algebraic

Symbol                          | Units  | Description                          | Note
--------------------------------|--------|--------------------------------------|------
$V_D$                           | [p.u.] | Speed/load fuel demand               | LV gate input
$V_T$                           | [p.u.] | Temperature-limit fuel demand        | LV gate input
$V$                             | [p.u.] | LV gate output                       |
$P_\mathrm{m}$                  | [p.u.] | Mechanical power to generator        | System base; `pmech` port

### External Variables

#### Differential

None.

#### Algebraic

Symbol                  | Units  | Type    | Description                 | Note
------------------------|--------|---------|-----------------------------|------
$\omega$                | [p.u.] | Known   | Machine speed deviation     | Optional signal port `speed`; defaults to zero
$P^\mathrm{ref}$        | [p.u.] | Unknown | Active-power/load reference | Optional signal port `pref`; system base

## Model Equations

### Differential Equations

```math
\begin{aligned}
  0 &=
    -\dot{x}_V
    + \dfrac{s^\text{fix}}{T_1}
      \text{antiwindup}\!\left(
        x_V,
        V - x_V,
        V^{\min},
        V^{\max}
      \right) \\
  0 &=
    -\dot{x}_F
    + \dfrac{s^\text{fix}}{T_2}\left(-x_F + x_V\right) \\
  0 &=
    -\dot{x}_T
    + \dfrac{s^\text{fix}}{T_3}\left(-x_T + x_F\right)
\end{aligned}
```

CommonMath defines the [anti-windup](../../../../CommonMath.md#antiwindup)
function.

### Algebraic Equations

```math
\begin{aligned}
  0 &= -\omega + R(k^\mathrm{base}P^\mathrm{ref} - V_D) \\
  0 &= -V_T
       + A_T
       + K_T(A_T - x_T) \\
  0 &=
    -V
    + \text{min}\left(V_D, V_T\right) \\
  0 &= -k^\mathrm{base}P_\mathrm{m} + x_F - D^\mathrm{turb}\omega
\end{aligned}
```

CommonMath defines the smooth [minimum](../../../../CommonMath.md#min).

## Initialization

### Input Initialization

```math
\begin{aligned}
  \omega
    &\leftarrow \text{machine speed deviation, or }0\text{ if unattached} \\
  P_\mathrm{m}
    &\leftarrow \text{machine mechanical-power}
\end{aligned}
```

### Internal Initialization

```math
\begin{aligned}
  x_F
    &\leftarrow k^\mathrm{base}P_\mathrm{m} + D^\mathrm{turb}\omega \\
  V
    &\leftarrow x_F \\
  x_V
    &\leftarrow x_F \\
  x_T
    &\leftarrow x_F \\
  V_T
    &\leftarrow A_T + K_T(A_T - x_T) \\
  V_D
    &\leftarrow V_T - \rho^{-1}(V_T - x_F) \\
  \dot{x}_V, \dot{x}_F, \dot{x}_T
    &\leftarrow 0
\end{aligned}
```

Here $\rho^{-1}$ is GridKit's inverse smooth
[ramp](../../../../CommonMath.md#primitives) (`iramp`). Initialization requires
$V_T-x_F>0$. Valve-limit violations in non-Fixed modes produce a warning and
initialize at the dispatched value.

### Output Initialization

```math
\begin{aligned}
  P^\mathrm{ref}
    &\leftarrow
      \dfrac{1}{k^\mathrm{base}}
      \left(V_D + \dfrac{\omega}{R}\right)
\end{aligned}
```

The initialized reference is written to `pref` when attached; otherwise, it is
held internally.

## Monitorable Outputs

Output    | Units  | Description                              | Note
----------|--------|------------------------------------------|------
`pmech`   | [p.u.] | Mechanical power output $P_\mathrm{m}$   | System base
`xvalve`  | [p.u.] | Fuel-valve state $x_V$                   |
`xflow`   | [p.u.] | Fuel-flow state $x_F$                    |
`xtemp`   | [p.u.] | Exhaust-temperature feedback state $x_T$ |
`vload`   | [p.u.] | Speed/load fuel demand $V_D$             |
`vtemp`   | [p.u.] | Temperature-limit fuel demand $V_T$      |
