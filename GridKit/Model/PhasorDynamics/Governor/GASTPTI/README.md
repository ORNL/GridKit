# **Gas Turbine-Governor Model (GASTPTI)**

GASTPTI is a gas turbine-governor model for thermal generating units. In
GridKit it is represented as a governor model that reads machine speed
deviation and supplies mechanical power to the machine through a fuel-valve,
fuel-flow, and exhaust-temperature limiting chain.

Notes:
- Internal fuel-demand and mechanical-power quantities are on the per-unit base
  supplied to this model unless otherwise stated.
- The shared GAST_PTI/GASTD diagram shows the optional `Trate`
  turbine-rating base. GridKit parses and validates `Trate` for input
  compatibility, but this model does not rescale power quantities internally.
- The diagram also shows the GASTD speed deadband block (`dbL`/`dbH`). That
  block is only used by GASTD; GASTPTI uses the speed-deviation input
  $\omega$ directly.
- The speed/load demand and the temperature-limit demand meet at the LV gate;
  the fuel command is the lesser of the two.

## Block Diagram

Standard GASTPTI block diagram.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/PhasorDynamics/GASTPTI_diagram.png">

  Figure 1: Governor GASTPTI model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters

Symbol                          | Units    | Description                                  | Typical Value | Note
--------------------------------|----------|----------------------------------------------|---------------|------
$R$                             | [p.u.]   | Permanent droop                              | 0.05          | Block name: `R`; required nonzero
$T_1$                           | [sec]    | Fuel-valve time constant                     | 0.4           | Block name: `T1`; if zero, $x_{\mathrm{valve}}$ is algebraic
$T_2$                           | [sec]    | Fuel-flow time constant                      | 0.1           | Block name: `T2`; if zero, $x_{\mathrm{flow}}$ is algebraic
$T_3$                           | [sec]    | Exhaust-temperature time constant            | 3.0           | Block name: `T3`; if zero, $x_{\mathrm{temp}}$ is algebraic
$A_{\mathrm{T}}$                | [p.u.]   | Ambient-temperature load limit               | 1.0           | Block name: `At`
$K_{\mathrm{T}}$                | [p.u.]   | Exhaust-temperature feedback gain            | 2.0           | Block name: `Kt`
$V^{\max}$                      | [p.u.]   | Maximum fuel-valve/turbine-power limit       | 1.0           | Block name: `Vmax`
$V^{\min}$                      | [p.u.]   | Minimum fuel-valve/turbine-power limit       | 0.0           | Block name: `Vmin`
$D_{\mathrm{turb}}$             | [p.u.]   | Turbine damping coefficient                  | 0.0           | Block name: `Dturb`; multiplies speed deviation
$P^{\mathrm{rate}}$             | [MW]     | Optional turbine-rating power base           | 0.0           | Block name: `Trate`; parsed for input compatibility

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
  &P^{\mathrm{rate}} \ge 0
\end{aligned}
```

### Model Derived Parameters

The droop gain $1/R$ and the temperature-feedback combination
$A_{\mathrm{T}} + K_{\mathrm{T}}(A_{\mathrm{T}} - x_{\mathrm{temp}})$ are
formed inline in the algebraic equations. No internal base-conversion factor is
formed from `Trate`; callers should supply GASTPTI power quantities on the
intended per-unit base.

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
$P_{\mathrm{mech}}$             | [p.u.] | Mechanical power to generator        | Read by a machine model

### External Variables

#### Differential

None.

#### Algebraic

Symbol                  | Units  | Description                 | Note
------------------------|--------|-----------------------------|------
$\omega$                | [p.u.] | Machine speed deviation     | Read from a machine model
$P_{\mathrm{ref}}$      | [p.u.] | Active-power/load reference | External setpoint or initialized constant; source label: `Pref`

## Model Equations

For readability, define:

```math
\begin{aligned}
  f_{\mathrm{valve}} &= V_{\mathrm{LV}} - x_{\mathrm{valve}}
\end{aligned}
```

### Differential Equations

The lag residuals are written in descriptor form. When $T_1$, $T_2$, or
$T_3$ is zero, the corresponding residual becomes algebraic.

```math
\begin{aligned}
  0 &=
    -T_1 \dot x_{\mathrm{valve}}
    + \text{antiwindup}\!\left(
        x_{\mathrm{valve}},
        f_{\mathrm{valve}},
        V^{\min},
        V^{\max}
      \right) \\
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
  0 &= -V_{\mathrm{load}} + P_{\mathrm{ref}} - \dfrac{\omega}{R} \\
  0 &= -V_{\mathrm{temp}}
       + A_{\mathrm{T}}
       + K_{\mathrm{T}}\!\left(A_{\mathrm{T}} - x_{\mathrm{temp}}\right) \\
  0 &= -V_{\mathrm{LV}} + \text{min}\!\left(V_{\mathrm{load}}, V_{\mathrm{temp}}\right) \\
  0 &= -P_{\mathrm{mech}} + x_{\mathrm{flow}} - D_{\mathrm{turb}}\omega
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
$P_{\mathrm{mech},0}$. The output residual fixes the fuel-flow state, and the
two plain lag residuals pass that value through to the fuel valve and
exhaust-temperature states:

```math
\begin{aligned}
  x_{\mathrm{flow},0}  &= P_{\mathrm{mech},0} + D_{\mathrm{turb}}\omega_0 \\
  x_{\mathrm{valve},0} &= x_{\mathrm{flow},0} \\
  x_{\mathrm{temp},0}  &= x_{\mathrm{flow},0}
\end{aligned}
```

Then evaluate the temperature-limit branch and the LV-gate target:

```math
\begin{aligned}
  V_{\mathrm{temp},0}
    &= A_{\mathrm{T}}
       + K_{\mathrm{T}}\!\left(A_{\mathrm{T}} - x_{\mathrm{temp},0}\right) \\
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
`pmech`        | [p.u.] | Mechanical power output            | Oriented as mechanical input to the machine
`fuelvalve`    | [p.u.] | Fuel-valve state                   | State 1
`fuelflow`     | [p.u.] | Fuel-flow state                    | State 2
`exhausttemp`  | [p.u.] | Exhaust-temperature feedback state | State 3
`vload`        | [p.u.] | Speed/load fuel demand             | LV gate input
`vtemp`        | [p.u.] | Temperature-limit fuel demand      | LV gate input
`vlv`          | [p.u.] | LV gate output                     | Selected fuel demand
