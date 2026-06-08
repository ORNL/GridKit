# **Gas Turbine-Governor Model (GASTPTI)**

GASTPTI is a gas turbine-governor model for thermal generating units. In
GridKit it is represented as a governor model that reads machine speed
deviation and supplies mechanical power to the machine through a fuel-valve,
fuel-flow, and exhaust-temperature limiting chain.

Notes:
- GASTPTI uses `Trate` as its turbine-power component base.
- The shared `pmech` signal is treated as system-base power. GridKit converts
  between system base and `Trate` at the GASTPTI signal boundary.
- The diagram shows the GASTD speed deadband block (`dbL`/`dbH`). That
  block is only used by GASTD; GASTPTI uses the speed-deviation input
  $\omega$ directly.

## Block Diagram

Standard GASTPTI block diagram.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/PhasorDynamics/GASTPTI_diagram.png">

  Figure 1: Governor GASTPTI model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters

Symbol                          | Units    | JSON     | Description                                  | Typical Value | Note
--------------------------------|----------|----------|----------------------------------------------|---------------|------
$R$                             | [p.u.]   | `R`      | Permanent droop                              | 0.05          | Required nonzero
$T_1$                           | [sec]    | `T1`     | Fuel-valve time constant                     | 0.4           | If zero, $x_{\mathrm{valve}}$ is algebraic
$T_2$                           | [sec]    | `T2`     | Fuel-flow time constant                      | 0.1           | If zero, $x_{\mathrm{flow}}$ is algebraic
$T_3$                           | [sec]    | `T3`     | Exhaust-temperature time constant            | 3.0           | If zero, $x_{\mathrm{temp}}$ is algebraic
$A_{\mathrm{T}}$                | [p.u.]   | `At`     | Ambient-temperature load limit               | 1.0           |
$K_{\mathrm{T}}$                | [p.u.]   | `Kt`     | Exhaust-temperature feedback gain            | 2.0           |
$V^{\max}$                      | [p.u.]   | `Vmax`   | Maximum fuel-valve/turbine-power limit       | 1.0           |
$V^{\min}$                      | [p.u.]   | `Vmin`   | Minimum fuel-valve/turbine-power limit       | 0.0           |
$D_{\mathrm{turb}}$             | [p.u.]   | `Dturb`  | Turbine damping coefficient                  | 0.0           | Multiplies speed deviation
$T_{\mathrm{rate}}$             | [MW]     | `Trate`  | Turbine-rating power base                    | 100.0         | Required positive value

### Parameter Validation

Invalid GASTPTI parameter sets are rejected by the following checks.

The required checks are:

```math
\begin{aligned}
  &R > 0 \\
  &T_1, T_2, T_3 \ge 0 \\
  &A_{\mathrm{T}} \ge 0 \\
  &K_{\mathrm{T}} \ge 0 \\
  &V^{\min} \le V^{\max} \\
  &D_{\mathrm{turb}} \ge 0 \\
  &T_{\mathrm{rate}} > 0
\end{aligned}
```

### Model Derived Parameters

The system power base $S_{\mathrm{sys}}$ is derived from the system model.

## Model Variables

### Internal Variables

#### Differential

Symbol                  | Units  | Description                        | Note
------------------------|--------|------------------------------------|------
$x_{\mathrm{valve}}$    | [p.u.] | Fuel-valve state                   | State 1 in Fig. 1; source label: `Fuel Valve`; algebraic when $T_1 = 0$
$x_{\mathrm{flow}}$     | [p.u.] | Fuel-flow state                    | State 2 in Fig. 1; source label: `Fuel Flow`; algebraic when $T_2 = 0$
$x_{\mathrm{temp}}$     | [p.u.] | Exhaust-temperature feedback state | State 3 in Fig. 1; source label: `Exhaust Temperature`; algebraic when $T_3 = 0$

#### Algebraic

Symbol                          | Units  | Description                          | Note
--------------------------------|--------|--------------------------------------|------
$V_{\mathrm{load}}$             | [p.u.] | Speed/load fuel demand               | $P_{\mathrm{ref}}$ less droop feedback; LV gate input
$V_{\mathrm{temp}}$             | [p.u.] | Temperature-limit fuel demand        | Exhaust-temperature branch; LV gate input
$V_{\mathrm{LV}}$               | [p.u.] | LV gate output                       | Lesser of $V_{\mathrm{load}}$ and $V_{\mathrm{temp}}$; drives the fuel-valve lag
$P_{\text{m}}$                  | [p.u.] | Mechanical power to generator        | On system base

### External Variables

#### Differential

None.

#### Algebraic

Symbol                  | Units  | Description                 | Note
------------------------|--------|-----------------------------|------
$\omega$                | [p.u.] | Machine speed deviation     | Read from a machine model
$P_{\mathrm{ref}}$      | [p.u.] | Active-power/load reference | On `Trate` component base; source label: `Pref`

## Model Equations

### Differential Equations

The lag residuals are written in descriptor form. When $T_1$, $T_2$, or
$T_3$ is zero, the corresponding residual becomes algebraic.

```math
\begin{aligned}
  0 &=
    -T_1 \dot x_{\mathrm{valve}}
    + \text{antiwindup}(
        x_{\mathrm{valve}},
        V_{\mathrm{LV}} - x_{\mathrm{valve}},
        V^{\min},
        V^{\max}
      ) \\
  0 &= -T_2 \dot x_{\mathrm{flow}} - x_{\mathrm{flow}} + x_{\mathrm{valve}} \\
  0 &= -T_3 \dot x_{\mathrm{temp}} - x_{\mathrm{temp}} + x_{\mathrm{flow}}
\end{aligned}
```

CommonMath defines the [Anti-Windup](../../../../CommonMath.md#anti-windup-indicator)
target and smooth approximation.

### Algebraic Equations

The algebraic targets form the two LV gate inputs, select between them with the
low-value gate, and assemble the mechanical power output.

```math
\begin{aligned}
  0 &= -\omega + R(P_{\mathrm{ref}} - V_{\mathrm{load}}) \\
  0 &= -V_{\mathrm{temp}}
       + A_{\mathrm{T}}
       + K_{\mathrm{T}}(A_{\mathrm{T}} - x_{\mathrm{temp}}) \\
  0 &= -V_{\mathrm{LV}} + \min(V_{\mathrm{load}}, V_{\mathrm{temp}}) \\
  0 &= -S_{\mathrm{sys}} P_{\text{m}}
       + T_{\mathrm{rate}}(x_{\mathrm{flow}} - D_{\mathrm{turb}}\omega)
\end{aligned}
```

CommonMath defines the helper targets and smooth approximations for
[min](../../../../CommonMath.md#derived-functions).

## Initialization

Initialization is performed by evaluating the steady-state residuals in
dependency order. Let subscript $0$ denote initial values and set all internal
derivatives to zero. A standard power-flow start uses synchronous machine
speed:

```math
\begin{aligned}
  \omega_0 &= 0
\end{aligned}
```

The machine initializes mechanical power first; GASTPTI reads it as
$P_{\text{m},0}$ on the system base and converts it to the `Trate`
component base:

```math
\begin{aligned}
  S_{\mathrm{sys}} P_{\text{m},0}
    &= T_{\mathrm{rate}}P^{\mathrm{gastpti}}_{\text{m},0}
\end{aligned}
```

The output residual fixes the fuel-flow state, and the two plain lag residuals
pass that value through to the fuel valve and exhaust-temperature states:

```math
\begin{aligned}
  x_{\mathrm{flow},0}  &=
    P^{\mathrm{gastpti}}_{\mathrm{mech},0}
    + D_{\mathrm{turb}}\omega_0 \\
  x_{\mathrm{valve},0} &= x_{\mathrm{flow},0} \\
  x_{\mathrm{temp},0}  &= x_{\mathrm{flow},0}
\end{aligned}
```

Then evaluate the temperature-limit branch and the LV-gate target:

```math
\begin{aligned}
  V_{\mathrm{temp},0}
    &= A_{\mathrm{T}}
       + K_{\mathrm{T}}(A_{\mathrm{T}} - x_{\mathrm{temp},0}) \\
  V_{\mathrm{LV},0} &= x_{\mathrm{valve},0}
\end{aligned}
```

For an unsaturated start where the speed/load demand is selected by the LV
gate, choose or verify:

```math
\begin{aligned}
  V_{\mathrm{load},0} &= x_{\mathrm{valve},0} \\
  P_{\mathrm{ref},0}
    &= V_{\mathrm{load},0} + \dfrac{\omega_0}{R}
\end{aligned}
```

This start requires $V^{\min} \le x_{\mathrm{valve},0} \le V^{\max}$ and
$V_{\mathrm{temp},0} \ge V_{\mathrm{load},0}$ so the LV gate selects the
speed/load demand and the fuel-valve anti-windup residual holds. If the
supplied $P_{\mathrm{ref}}$ is connected, use its initial value and verify the
same residuals.

Starts that bind the fuel-valve limit or the temperature gate are outside these
closed-form initialization formulas. The zero-derivative lag relations still
require $x_{\mathrm{temp},0}=x_{\mathrm{flow},0}=x_{\mathrm{valve},0}$.

## Model Outputs

Output         | Units  | Description                        | Note
---------------|--------|------------------------------------|------
`pmech`        | [p.u.] | Mechanical power output            | On system base; oriented as mechanical input to the machine
`fuelvalve`    | [p.u.] | Fuel-valve state                   | State 1
`fuelflow`     | [p.u.] | Fuel-flow state                    | State 2
`exhausttemp`  | [p.u.] | Exhaust-temperature feedback state | State 3
`vload`        | [p.u.] | Speed/load fuel demand             | LV gate input
`vtemp`        | [p.u.] | Temperature-limit fuel demand      | LV gate input
`vlv`          | [p.u.] | LV gate output                     | Selected fuel demand
