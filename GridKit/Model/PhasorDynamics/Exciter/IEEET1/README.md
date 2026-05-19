# **IEEE Type 1 Excitation System Model (IEEET1)**

## Block Diagram

Standard model of the IEEET1 Exciter.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/PhasorDynamics_IEEET1_Diagram.png">


  Figure 1: Exciter IEEET1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters


Symbol                 | Units    | Description                          | Typical Value | Note
-----------------------|----------|--------------------------------------|---------------|------
$T_R$                  | [sec]    | Time constant for voltage sensing    | 0             |
$K_\mathrm{a}$         | [p.u.]   | Coefficient for voltage regulation   | 50            |
$T_\mathrm{a}$         | [sec]    | Time constant for voltage regulation | 0.04          |
$K_\mathrm{e}$         | [p.u.]   | Coefficient for excitation system    | -0.06         |
$T_\mathrm{e}$         | [sec]    | Time constant for excitation system  | 0.6           |
$K_\mathrm{f}$         | [p.u.]   | Coefficient for feedback             | 0.09          |
$T_\mathrm{f}$         | [sec]    | Time constant for feedback           | 1.46          |
$V_R^{\min}$           | [p.u.]   | Lower limit to voltage regulation    | -1            |
$V_R^{\max}$           | [p.u.]   | Upper limit to voltage regulation    | 1             |
$E_1$                  | [p.u.]   | Saturation Parameter                 | 2.8           |
$E_2$                  | [p.u.]   | Saturation Parameter                 | 3.73          |
$S_{\mathrm{e}1}$      | [p.u.]   | Saturation Parameter                 | 0.04          |
$S_{\mathrm{e}2}$      | [p.u.]   | Saturation Parameter                 | 0.33          |
$I_\mathrm{spdlm}$     | [binary] | Speed Limit flag indicator           | 0             |

### Model Derived Parameters

The relationship of the derived parameters is defined by the following quadratic model. The parameters are chosen so that the quadratic model represents
the expected saturation near the operating region.
``` math
\begin{aligned}
  S_{\mathrm{e}1} &= S_B(E_1-S_A)^2 \\
  S_{\mathrm{e}2} &= S_B(E_2-S_A)^2 \\
\end{aligned}
```

Generally, this system has two solutions. The non-extraneous solution is as follows.
``` math
\begin{aligned}
  C &=  \sqrt{
   \dfrac
   {S_{\mathrm{e}2}}
   {S_{\mathrm{e}1}}
  }
  \\
  S_A &=
   \dfrac
   {C E_1 - E_2}
   {C - 1}
  \\
  S_B &=
   \dfrac
   {S_{\mathrm{e}1}}
   {(E_1-S_A)^2}
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol             | Units  | Description                       | Note
-------------------|--------|-----------------------------------|-------
$V_\mathrm{ts}$    | [p.u.] | Sensed terminal voltage           |
$V_R$              | [p.u.] | Voltage regulator                 |
$E_\mathrm{fd}'$   | [p.u.] | Field-current pre-speed multiplier|
$V_\mathrm{fx}$    | [p.u.] | Exciter feedback internal state   |


#### Algebraic


Symbol             | Units  | Description                       | Note
-------------------|--------|-----------------------------------|-------
$V_\mathrm{tr}$    | [p.u.] | Terminal Voltage Error            |
$V_\mathrm{f}$     | [p.u.] | Feedback Voltage                  |
$V_E$              | [p.u.] | Excitation control voltage        |
$E_\mathrm{fd}$    | [p.u.] | Field winding voltage             |
$k_\mathrm{sat}$   | [p.u.] | Saturation variable               |

### External Variables

#### Differential
Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$\omega$  | [p.u.] | Machine Speed Deviation                   | Read from a Machine Model

#### Algebraic


Symbol             | Units  | Description                       | Note
-------------------|--------|-----------------------------------|-------
$E_C$              | [p.u.] | Compensated machine terminal voltage magnitude |
$V_\mathrm{ref}$   | [p.u.] | Reference terminal voltage                     |
$V_\mathrm{UEL}$   | [p.u.] | Input from under excitation limiter            |
$V_\mathrm{OEL}$   | [p.u.] | Input from over excitation limiter             |
$V_S$              | [p.u.] | Input from stabilizer controller               |


## Model Equations

### Differential Equations

The IEEET1 differential equations are written in $T_x \dot{x}$ form so that
each time constant appears as a multiplier on its state derivative. When a time
constant is zero, the corresponding residual becomes algebraic.

```math
\begin{aligned}
   0 &= -T_R \, \dot{V}_\mathrm{ts}     + E_C - V_\mathrm{ts} \\
   0 &= -T_A \, \dot{V}_R               + \mathrm{antiwindup}\bigl(V_R,\ K_\mathrm{a} V_\mathrm{tr} - V_R,\ V_R^{\min},\ V_R^{\max}\bigr) \\
   0 &= -T_E \, \dot{E}_\mathrm{fd}'    + V_R - V_E - K_E E_\mathrm{fd}' \\
   0 &= -T_F \, \dot{V}_\mathrm{fx}     + V_\mathrm{f} \\
\end{aligned}
```

The piecewise hard-limit on $\dot V_R$ is encapsulated by the `antiwindup`
helper, which takes the unlimited driver $K_\mathrm{a} V_\mathrm{tr} - V_R$ and
passes it through GridKit's smooth anti-windup indicator. See
[CommonMath: Anti-Windup Indicator](../../../../CommonMath.md#anti-windup-indicator)
for the helper's definition, behavior, and design rationale.

### Algebraic Equations

The algebraic equations of the exciter.
```math
\begin{aligned}
   0 &= -V_\mathrm{tr}    + V_\mathrm{ref} + V_\mathrm{UEL} + V_\mathrm{OEL} + V_S - V_\mathrm{f} - V_\mathrm{ts} \\
   0 &= -T_F V_\mathrm{f} + K_F E_\mathrm{fd}' - T_F V_\mathrm{fx} \\
   0 &= -V_E              + k_\mathrm{sat} E_\mathrm{fd}' \\
   0 &= -E_\mathrm{fd}    + (1 + \omega I_\mathrm{spdlm})\, E_\mathrm{fd}' \\
   0 &= -k_\mathrm{sat}   +
   \begin{cases}
      S_B(E_\mathrm{fd}' - S_A)^2 & \text{if } E_\mathrm{fd}' > S_A \\
      0                           & \text{else}
   \end{cases}
\end{aligned}
```
#### Smooth Piecewise Approximation (Algebraic)

The $E_\mathrm{fd}$ relation already collapses to a single expression by absorbing the
$I_\mathrm{spdlm}$ flag as a multiplier. The saturation $k_\mathrm{sat}$ uses GridKit's
[quadratic ramp](../../../../CommonMath.md#primitives), $q(x)$:
```math
   k_\mathrm{sat} = S_B q(E_\mathrm{fd}' - S_A)
```

The approximation approaches the exact one-sided quadratic as the sigmoid
steepness increases.

## Initialization

The machine initializes $E_\mathrm{fd}$ first. IEEET1 reads that value as
$E_{\mathrm{fd},0}$, along with attached $E_C$, $\omega$, and $V_S$ signals, and
solves the steady-state algebraic chain so all residuals vanish with
$\dot y = 0$. Saturation and the speed-limit flag are included directly;
$V_\mathrm{ref}$ is set to close the $V_\mathrm{tr}$ equation with the current
auxiliary inputs.

```math
\begin{aligned}
   E_C              &= E_{C,\mathrm{signal}} \\
   E_\mathrm{fd}'   &= \dfrac{E_{\mathrm{fd},0}}{1 + I_\mathrm{spdlm}\,\omega} \\
   k_\mathrm{sat}   &= S_B q(E_\mathrm{fd}' - S_A) \\
   V_E              &= k_\mathrm{sat}\, E_\mathrm{fd}' \\
   V_R              &= K_E\, E_\mathrm{fd}' + V_E \\
   V_\mathrm{tr}    &= \dfrac{V_R}{K_\mathrm{a}} \\
   V_\mathrm{fx}    &= \dfrac{K_F}{T_F}\, E_\mathrm{fd}' \\
   V_\mathrm{ts}    &= E_C \\
   V_\mathrm{f}     &= 0 \\
   V_\mathrm{ref}   &= E_C + V_\mathrm{tr} - V_\mathrm{UEL} - V_\mathrm{OEL} - V_S
\end{aligned}
```

All internal derivatives initialize to zero.
## Model Outputs

The field voltage, $E_\mathrm{fd}$, is an internal model variable.

The magnetic saturation coefficient $k_\mathrm{sat}$ is calculated from $E_\mathrm{fd}$ using the smooth piecewise version above.
