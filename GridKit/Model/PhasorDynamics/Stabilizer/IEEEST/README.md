# **IEEE Stabilizer Model (IEEEST)**

Standard IEEE power system stabilizer: 4th-order notch filter, two lead–lag
blocks, washout, and output limiter.

Notes:
- $V_{\mathrm{cl}}$, $V_{\mathrm{cu}}$, and $T_{\mathrm{delay}}$ are accepted for input-format
  compatibility but are not modeled.
- A zero denominator time constant bypasses its corresponding lead–lag or
  washout block.

## Block Diagram

Standard IEEEST block diagram.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/stabilizer_ieeest_diagram.png">

  Figure 1: Stabilizer IEEEST model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters

Symbol                  | Units    | JSON     | Description                          | Typical Value | Note
------------------------|----------|----------|--------------------------------------|---------------|------
$A_1$                   | [sec]    | `A1`     | Notch denominator coefficient        | 1.013         |
$A_2$                   | [sec²]   | `A2`     | Notch denominator coefficient        | 0.013         |
$A_3$                   | [sec]    | `A3`     | Notch denominator coefficient        | 0.0           |
$A_4$                   | [sec²]   | `A4`     | Notch denominator coefficient        | 0.0           |
$A_5$                   | [sec]    | `A5`     | Notch numerator coefficient          | 1.013         |
$A_6$                   | [sec²]   | `A6`     | Notch numerator coefficient          | 0.113         |
$T_1$                   | [sec]    | `T1`     | Lead–lag 1 numerator time constant   | 0.0           |
$T_2$                   | [sec]    | `T2`     | Lead–lag 1 denominator time constant | 0.02          |
$T_3$                   | [sec]    | `T3`     | Lead–lag 2 numerator time constant   | 0.0           |
$T_4$                   | [sec]    | `T4`     | Lead–lag 2 denominator time constant | 0.0           |
$T_5$                   | [sec]    | `T5`     | Washout numerator time constant      | 1.65          |
$T_6$                   | [sec]    | `T6`     | Washout denominator time constant    | 1.65          |
$K_s$                   | [p.u.]   | `Ks`     | Stabilizer gain                      | 3.0           |
$L_s^{\min}$      | [p.u.]   | `Lsmin`  | Minimum stabilizer output limit      | -0.1          |
$L_s^{\max}$      | [p.u.]   | `Lsmax`  | Maximum stabilizer output limit      | 0.1           |
$V_{\mathrm{cl}}$       | [p.u.]   | `Vcl`    | Lower input cutout threshold         | 0.0           | Accepted but not modeled
$V_{\mathrm{cu}}$       | [p.u.]   | `Vcu`    | Upper input cutout threshold         | 0.0           | Accepted but not modeled
$T_{\mathrm{delay}}$    | [sec]    | `Tdelay` | Input delay                          | 0.0           | Accepted but not modeled

### Parameter Validation

The fixed realization rejects a first-order-only notch denominator:

```math
\begin{aligned}
  a_2 \ne 0 \lor a_3 \ne 0 \lor a_4 \ne 0 \lor a_1 = 0
\end{aligned}
```

### Model Derived Parameters

The notch-filter denominator expands to:

```math
\begin{aligned}
a_1 &= A_1 + A_3 \\
a_2 &= A_2 + A_4 + A_1 A_3 \\
a_3 &= A_1 A_4 + A_2 A_3 \\
a_4 &= A_2 A_4
\end{aligned}
```

The binary DAE selectors choose the active notch-filter order:

```math
\begin{aligned}
\delta_1 &= \delta_2 =
\begin{cases}
1 & a_2 \ne 0 \lor a_3 \ne 0 \lor a_4 \ne 0 \\
0 & \text{otherwise}
\end{cases}
\\
\delta_3 &=
\begin{cases}
1 & a_3 \ne 0 \lor a_4 \ne 0 \\
0 & \text{otherwise}
\end{cases}
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol                | Units  | Description           | Note
----------------------|--------|-----------------------|------
$x_1, x_2, x_3, x_4$  | [-]    | Notch filter states   | States 1–4 in Fig. 1
$x_5$                 | [-]    | Lead–lag 1 state      | State 5 in Fig. 1
$x_6$                 | [-]    | Lead–lag 2 state      | State 6 in Fig. 1
$x_7$                 | [-]    | Washout state         | State 7 in Fig. 1

For reduced-order notch filters, unused notch states remain in the fixed
component state vector and are pinned by algebraic residuals.

#### Algebraic

Symbol               | Units  | Description                              | Note
---------------------|--------|------------------------------------------|------
$v_4$                | [p.u.] | Notch filter output                      |
$v_5$                | [p.u.] | Lead–lag 1 output                        |
$v_6$                | [p.u.] | Lead–lag 2 output                        |
$v_7$                | [p.u.] | Unlimited stabilizer signal              |
$V_{\mathrm{ss}}$    | [p.u.] | Limited stabilizer signal (model output) |

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description             | Note
-------|--------|-------------------------|------
$u$    | [p.u.] | Stabilizer input signal |

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= -\delta_1\dot{x}_1 + x_2 \\
0 &= -\delta_2\dot{x}_2 + x_3 \\
0 &= -\delta_3\dot{x}_3 + x_4 \\
0 &= -a_4\dot{x}_4 - x_1 - a_1x_2 - a_2x_3 - a_3x_4 + u \\
0 &= -T_2 \dot{x}_5 - x_5 + v_4 \\
0 &= -T_4 \dot{x}_6 - x_6 + v_5 \\
0 &= -T_6 \dot{x}_7 - x_7 + v_6
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
0 &= -v_4 + x_1 + A_5x_2 + A_6x_3 \\
0 &=
\begin{cases}
-T_2(v_5 - x_5) + T_1(v_4 - x_5) & T_2 \ne 0 \\
v_4 - v_5 & T_2 = 0
\end{cases} \\
0 &=
\begin{cases}
-T_4(v_6 - x_6) + T_3(v_5 - x_6) & T_4 \ne 0 \\
v_5 - v_6 & T_4 = 0
\end{cases} \\
0 &=
\begin{cases}
-T_6 v_7 + K_s T_5(v_6 - x_7) & T_6 \ne 0 \\
K_s v_6 - v_7 & T_6 = 0
\end{cases} \\
0 &= -V_{\mathrm{ss}} + \text{clamp}(v_7, L_s^{\min}, L_s^{\max})
\end{aligned}
```

The output limiter uses GridKit's smooth
[clamp](../../../../CommonMath.md#derived-functions).

## Initialization

Initialization is performed by evaluating the steady-state residuals in
dependency order. Let subscript $0$ denote initial values and set all internal
derivatives to zero. States and derivatives initialize to the steady state
implied by the attached input $u$:

```math
\begin{aligned}
x_1 &= v_4 = x_5 = v_5 = x_6 = v_6 = x_7 = u \\
x_2 &= x_3 = x_4 = 0 \\
v_7 &=
\begin{cases}
K_su & T_6 = 0 \\
0 & \text{otherwise}
\end{cases}
&
V_{\mathrm{ss}} &= \text{clamp}(v_7, L_s^{\min}, L_s^{\max})
\end{aligned}
```

All internal derivatives initialize to zero.

## Model Outputs

Output | Units  | Description               | Note
-------|--------|---------------------------|------
`vss`  | [p.u.] | Limited stabilizer signal | Exported through `output` when assigned
