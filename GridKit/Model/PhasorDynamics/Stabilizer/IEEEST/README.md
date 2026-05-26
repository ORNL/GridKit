# **IEEE Stabilizer Model (IEEEST)**

Standard IEEE power system stabilizer: 4th-order notch filter, two lead–lag
blocks, washout, and output limiter.

## Block Diagram

<div align="center">
   <img align="center" src="../../../../../docs/Figures/stabilizer_ieeest_diagram.png">

  Figure 1: Stabilizer IEEEST model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters

Symbol      | Units  | Description                          | Typical Value
------------|--------|--------------------------------------|--------------
$A_1$       | [s]    | Notch denominator coefficient        | 1.013
$A_2$       | [s²]   | Notch denominator coefficient        | 0.013
$A_3$       | [s]    | Notch denominator coefficient        | 0.0
$A_4$       | [s²]   | Notch denominator coefficient        | 0.0
$A_5$       | [s]    | Notch numerator coefficient          | 1.013
$A_6$       | [s²]   | Notch numerator coefficient          | 0.113
$T_1$       | [s]    | Lead–lag 1 numerator time constant   | 0.0
$T_2$       | [s]    | Lead–lag 1 denominator time constant | 0.02
$T_3$       | [s]    | Lead–lag 2 numerator time constant   | 0.0
$T_4$       | [s]    | Lead–lag 2 denominator time constant | 0.0
$T_5$       | [s]    | Washout numerator time constant      | 1.65
$T_6$       | [s]    | Washout denominator time constant    | 1.65
$K_s$       | [p.u.] | Stabilizer gain                      | 3.0
$L_s^{\min}$ | [p.u.] | Minimum stabilizer output limit      | -0.1
$L_s^{\max}$ | [p.u.] | Maximum stabilizer output limit      | 0.1

The IEEE 421.5 IEEEST also defines a cutout window ($V_{cl}$, $V_{cu}$) and an
input delay ($T_{delay}$). These parameters are accepted for input-format
compatibility but are not modeled here.

### Derived Parameters

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

Symbol                | Units  | Description
----------------------|--------|------------
$x_1, x_2, x_3, x_4$  | [-]    | Notch filter states
$x_5$                 | [-]    | Lead–lag 1 state
$x_6$                 | [-]    | Lead–lag 2 state
$x_7$                 | [-]    | Washout state

#### Algebraic

Symbol     | Units  | Description
-----------|--------|------------
$v_4$      | [p.u.] | Notch filter output
$v_5$      | [p.u.] | Lead–lag 1 output
$v_6$      | [p.u.] | Lead–lag 2 output
$v_7$      | [p.u.] | Unlimited stabilizer signal
$V_{ss}$   | [p.u.] | Limited stabilizer signal (model output)

### External Variables

#### Algebraic

Symbol | Units  | Description
-------|--------|------------
$u$    | [p.u.] | Stabilizer input signal

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= -\dot{x}_1 + x_2 \\
0 &= -\dot{x}_2 + x_3 \\
0 &= -\dot{x}_3 + x_4 \\
0 &= -\dot{x}_4 - \dfrac{a_0}{a_4}x_1 - \dfrac{a_1}{a_4}x_2 - \dfrac{a_2}{a_4}x_3 - \dfrac{a_3}{a_4}x_4 + \dfrac{1}{a_4}u \\
0 &= -T_2 \dot{x}_5 - x_5 + v_4 \\
0 &= -T_4 \dot{x}_6 - x_6 + v_5 \\
0 &= -T_6 \dot{x}_7 - x_7 + v_6
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
0 &= -v_4 + x_1 + A_5 x_2 + A_6 x_3 \\
0 &= -T_2(v_5 - x_5) + T_1(v_4 - x_5) \\
0 &= -T_4(v_6 - x_6) + T_3(v_5 - x_6) \\
0 &= -T_6 v_7 + K_s T_5(v_6 - x_7) \\
0 &= -V_{ss} + \text{clamp}(v_7, L_s^{\min}, L_s^{\max})
\end{aligned}
```

The output limiter uses GridKit's smooth
[Clamp](../../../../CommonMath.md#derived-functions).

## Initialization

All states and their derivatives initialize to zero. The stabilizer comes
online at rest and produces signal only in response to deviations in the input
$u$.
