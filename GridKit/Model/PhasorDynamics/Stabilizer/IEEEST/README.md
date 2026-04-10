# IEEEST

The **IEEEST** model is a standard IEEE power system stabilizer used in transient stability simulations.
It consists of a 4th-order notch filter, two lead–lag blocks, a washout block, and an output limiter with input cutout logic.

Notes:
- The **cutout logic uses** $V_{ct}$ (as labeled in the block diagram), not $u_d$.

## Block Diagram

<div align="center">
   <img align="center" src="../../../../../docs/Figures/stabilizer_ieeest_diagram.png">

  Figure 1: Stabilizer IEEEST model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>


## Model Parameters

Symbol | Units | Description | Typical Value | Note
------ | ----- | ----------- | ----
$A_1$ | [s] | Notch filter denominator coefficient | 1.013
$A_2$ | [s²] | Notch filter denominator coefficient | 0.013
$A_3$ | [s] | Notch filter denominator coefficient | 0.0
$A_4$ | [s²] | Notch filter denominator coefficient | 0.0
$A_5$ | [s] | Notch filter numerator coefficient | 1.013
$A_6$ | [s²] | Notch filter numerator coefficient | 0.113
$T_1$ | [s] | Lead–lag 1 numerator time constant | 0.0
$T_2$ | [s] | Lead–lag 1 denominator time constant | 0.02
$T_3$ | [s] | Lead–lag 2 numerator time constant | 0.0
$T_4$ | [s] | Lead–lag 2 denominator time constant | 0.0
$T_5$ | [s] | Washout numerator time constant | 1.65
$T_6$ | [s] | Washout denominator time constant | 1.65
$K_s$ | [p.u.] | Stabilizer gain | 3.0
$L_{s\min}$ | [p.u.] | Minimum stabilizer output limit | 0.1
$L_{s\max}$ | [p.u.] | Maximum stabilizer output limit | -0.1
$V_{cl}$ | [p.u.] | Lower input cutout threshold | 0.0
$V_{cu}$ | [p.u.] | Upper input cutout threshold | 0.0
$T_{delay}$ | [s] | Input time delay | 0.0

### Model Derived Parameters
```math
\begin{aligned}
a_0 &= 1 \\
a_1 &= A_1 + A_3 \\
a_2 &= A_2 + A_4 + A_1 A_3 \\
a_3 &= A_1 A_4 + A_2 A_3 \\
a_4 &= A_2 A_4
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$x_1$ | [-] | Notch filter state |
$x_2$ | [-] | Notch filter state |
$x_3$ | [-] | Notch filter state |
$x_4$ | [-] | Notch filter state |
$x_5$ | [-] | Lead–lag 1 state |
$x_6$ | [-] | Lead–lag 2 state |
$x_7$ | [-] | Washout state |

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$v_4$ | [p.u.] | Notch filter output |
$v_5$ | [p.u.] | Lead–lag 1 output |
$v_6$ | [p.u.] | Lead–lag 2 output |
$v_7$ | [p.u.] | Unlimited stabilizer signal |
$V_{ss}$ | [p.u.] | Limited stabilizer signal |
$V_s$ | [p.u.] | Stabilizer output |

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$u$      | [p.u.] | Stabilizer input signal |
$V_{ct}$ | [p.u.] | Cutout signal (compared to $V_{cl},V_{cu}$) | from the block diagram

## Model Equations

### Differential Equations
```math
\begin{aligned}
\dot{x}_1 &= x_2 \\
\dot{x}_2 &= x_3 \\
\dot{x}_3 &= x_4 \\
\dot{x}_4 &= -\dfrac{a_0}{a_4}x_1
            -\dfrac{a_1}{a_4}x_2
            -\dfrac{a_2}{a_4}x_3
            -\dfrac{a_3}{a_4}x_4
            +\dfrac{1}{a_4}u\\
\dot{x}_5 &= \dfrac{1}{T_2}(v_4 - x_5) \\
\dot{x}_6 &= \dfrac{1}{T_4}(v_5 - x_6)\\
\dot{x}_7 &= \dfrac{1}{T_6}(v_6 - x_7)
\end{aligned}
```

### Algebraic Equations
```math
\begin{aligned}
0 &= -v_4 + x_1 + A_5 x_2 + A_6 x_3 \\
0 &= -v_5 + x_5 + \dfrac{T_1}{T_2}(v_4 - x_5) \\
0 &= -v_6 + x_6 + \dfrac{T_3}{T_4}(v_5 - x_6) \\
0 &= -v_7 + K_s \dfrac{T_5}{T_6}(v_6 - x_7) \\
0 &= -V_{ss} + \min\!\big(\max(v_7, L_{s\min}), L_{s\max}\big) \\
0 &= -V_s +
\begin{cases}
V_{ss}, & V_{cl} < V_{ct} < V_{cu} \\
0, & \text{otherwise}
\end{cases}
\end{aligned}
```
