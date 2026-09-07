# **Simplified Excitation System Model (SEXS-PTI)**

## Block Diagram

Simplified excitation system model ported from PhasorDynamics. The controller
equations and output anti-windup are unchanged; the EMT terminal port supplies
three-phase voltages in volts. An optional voltage-measurement lag defaults to
an algebraic bypass.

![](../../../../../../docs/Figures/SEXS_PTI_DIAGRAM.png)

Figure 1: Exciter SEXS-PTI model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol          | Units  | Description                                   | Typical Value | Note
----------------|--------|-----------------------------------------------|---------------|------
$V$             | [V]    | Rated line-to-line RMS terminal voltage       |               | Required
$T_r$           | [sec]  | Terminal-voltage measurement lag              | 0             | Zero bypasses
$T_A$           | [sec]  | Numerator time constant of lag-lead block     |               |
$T_B$           | [sec]  | Denominator time constant of lag-lead block   |               |
$T_E$           | [sec]  | Exciter field time constant                   |               |
$K$             | [p.u.] | Voltage regulator gain                        |               |
$E_{fd}^{\max}$ | [p.u.] | Maximum excitation output                     |               |
$E_{fd}^{\min}$ | [p.u.] | Minimum excitation output                     |               |

PowerWorld/PSS/E SEXS_PTI data often gives $T_A/T_B$ as a ratio. GridKit stores
$T_A$ and $T_B$ separately, so convert ratio-format data with
$T_A = (T_A/T_B)T_B$ before passing parameters to the model.

The six controller parameters and `V` are required. `Tr` defaults to zero.
All parameters must be finite, `V` must be positive, and `Tr` non-negative.

### Parameter Validation

Invalid SEXS-PTI parameter sets are rejected by the following checks:

```math
\begin{aligned}
  T_A &\ge 0 \\
  T_B, T_E, K &> 0 \\
  E_{fd}^{\min} &< E_{fd}^{\max}
\end{aligned}
```

### Model Derived Parameters

None.

## Model Ports

Name   | Port   | Init    | Description
-------|--------|---------|------------
`bus`  | Bus    | Known   | Three-phase terminal voltage in volts
`vref` | Input  | Known when attached | Voltage-control reference; inferred when unattached
`vs`   | Input  | Known   | Stabilizer input signal
`vuel` | Input  | Known   | Under-excitation limiter input
`voel` | Input  | Known   | Over-excitation limiter input
`efd`  | Output | Known   | Required field-voltage output seeded by the machine

## Model Variables

### Internal Variables

#### Differential

Symbol    | Units  | Description                       | Note
----------|--------|-----------------------------------|-----
$V_R$     | [p.u.] | Lag-lead block state              |
$E_{fd}$  | [p.u.] | Exciter field voltage output      |

#### Algebraic

Symbol    | Units  | Description                       | Note
----------|--------|-----------------------------------|-----
$V_{tr}$  | [p.u.] | Terminal voltage error signal     |

### External Variables

#### Differential

None.

#### Algebraic

Symbol          | Units  | Description                                  | Note
----------------|--------|----------------------------------------------|-----
$v_a,v_b,v_c$   | [V]    | Instantaneous phase voltages                  | Bus input
$V_{ref}$       | [p.u.] | Reference voltage                            | Signal port `vref`
$V_S$           | [p.u.] | Stabilizer output                            | Signal port `vs`
$V_{OEL}$       | [p.u.] | Over-excitation limiter signal               | Signal port `voel`
$V_{UEL}$       | [p.u.] | Under-excitation limiter signal              | Signal port `vuel`

## Model Equations

Define the compensated terminal voltage magnitude for readability:

```math
V_t = \frac{\sqrt{v_a^2+v_b^2+v_c^2}}{V}, \qquad
0 = -T_r \dot E_C + V_t - E_C.
```

`E_C` is an additional measured-voltage variable: differential when `Tr > 0`,
algebraic when `Tr = 0`. Balanced phase voltages of nominal line-to-line RMS
voltage `V` give `V_t = 1` at every phase angle.

### Internal Equations

#### Differential

The SEXS-PTI differential equations, as derived from the model diagram. Define the pre-limit derivative of $E_{fd}$

```math
f = \dfrac{1}{T_E}\left[-E_{fd} + \dfrac{K}{T_B}(-V_R + T_A V_{tr})\right]
```

so that $\dot E_{fd}$ can be written in piecewise form compactly.

```math
\begin{aligned}
  \dot V_R      &= -V_{tr} + \dfrac{1}{T_B}(-V_R + T_A V_{tr}) \\
  \dot E_{fd}   &=
  \begin{cases}
     f
        &  \text{if } (E_{fd}^{\min} < E_{fd} < E_{fd}^{\max}) & \lor \\
        &  \quad (E_{fd} \leq E_{fd}^{\min} \land f > 0)       & \lor \\
        &  \quad (E_{fd} \geq E_{fd}^{\max} \land f < 0)            \\
     0  &  \text{else}
  \end{cases}
\end{aligned}
```

In simulation the piecewise form above is replaced with a smooth approximation where $\phi$ is GridKit's smooth anti-windup indicator. See [CommonMath: Anti-Windup Indicator](../../../../CommonMath.md#antiwindup) for its definition, behavior, and design rationale.

#### Algebraic

```math
\begin{aligned}
0&=-V_{tr}-E_C+V_{ref}+V_S+V_{OEL}+V_{UEL}
\end{aligned}
```

### External Equations

None.

## Initialization

The generator initializes the EFD signal first. SEXS-PTI then reads that value
and any attached $V_S$, $V_{OEL}$, and $V_{UEL}$ signals and assumes steady state:

```math
\begin{aligned}
V_{tr,0} &= \dfrac{E_{fd,0}}{K} \\
V_{R,0} &= (T_A - T_B)V_{tr,0} \\
V_{ref} &= E_C + V_{tr,0} - V_S - V_{OEL} - V_{UEL}
\end{aligned}
```

The measured terminal voltage initializes to `V_t`. All derivatives initialize
to zero. An initial field voltage outside the configured limits is rejected.
Optional attached inputs are preserved and read live. The inferred reference
above is used only when `vref` is unattached. If a supplied reference differs,
consistent initialization resolves algebraic values and derivatives while
retaining the initialized differential states.

## Monitors

Monitor | Units  | Description          | Note
--------|--------|----------------------|------
`efd`   | [p.u.] | Field-voltage output | $E_{fd}$

`vts`, `vr`, and `vtr` also expose measured terminal voltage, the lead-lag state,
and voltage error, respectively.
