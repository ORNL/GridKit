# **IEEE Stabilizer Model (IEEEST)**

Standard IEEE power system stabilizer: a selectable-order notch filter, two
lead–lag blocks, washout, and output limiter.

## Notes

- The notch order $n$ is the degree of the expanded denominator, so $n$ is the
  largest index with $a_n \ne 0$. Order $0$ bypasses the filter.
- The numerator is truncated to the selected order. $A_5$ and $A_6$ are
  inactive for $n=0$, and $A_6$ is inactive for $n=1$.
- A zero denominator time constant bypasses its block. $T_2=0$ gives
  $v_5=v_4$, $T_4=0$ gives $v_6=v_5$, and $T_6=0$ gives $v_7=K_s v_6$.
- $V_{cl}$, $V_{cu}$, and $T_{delay}$ are accepted for input-format
  compatibility but are not modeled.

## Block Diagram

![](../../../../../docs/Figures/stabilizer_ieeest_diagram.png)

Figure 1: Stabilizer IEEEST model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

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

Symbol | Units       | Description                          | Note
-------|-------------|--------------------------------------|------
$x_1$  | [p.u.]      | Notch filter signal state            | Active for $n\ge1$
$x_2$  | [p.u./sec]  | First derivative of filtered signal  | Active for $n\ge2$
$x_3$  | [p.u./sec²] | Second derivative of filtered signal | Active for $n\ge3$
$x_4$  | [p.u./sec³] | Third derivative of filtered signal  | Active for $n=4$
$x_5$  | [p.u.]      | Lead–lag 1 state                     |
$x_6$  | [p.u.]      | Lead–lag 2 state                     |
$x_7$  | [p.u.]      | Washout state                        |

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

Only the notch states $x_1,\ldots,x_n$ are active. The remaining chain
equations reduce to $\dot{x}_i = 0$.

```math
\begin{aligned}
0 &= -\dot{x}_1 + \begin{cases}
      \dfrac{1}{a_1}\left(-a_0 x_1 + u\right) & n = 1 \\
      x_2 & n \in \{2,3,4\}
    \end{cases} \\
0 &= -\dot{x}_2 + \begin{cases}
      \dfrac{1}{a_2}\left(-a_0 x_1 - a_1 x_2 + u\right) & n = 2 \\
      x_3 & n \in \{3,4\}
    \end{cases} \\
0 &= -\dot{x}_3 + \begin{cases}
      \dfrac{1}{a_3}\left(-a_0 x_1 - a_1 x_2 - a_2 x_3 + u\right) & n = 3 \\
      x_4 & n = 4
    \end{cases} \\
0 &= -\dot{x}_4 + \dfrac{1}{a_4}\left(-a_0 x_1 - a_1 x_2 - a_2 x_3 - a_3 x_4 + u\right),
    \quad n = 4 \\
0 &= -T_2 \dot{x}_5 - x_5 + v_4 \\
0 &= -T_4 \dot{x}_6 - x_6 + v_5 \\
0 &= -T_6 \dot{x}_7 - x_7 + v_6
\end{aligned}
```

### Algebraic Equations

The numerator $1 + A_5 s + A_6 s^2$ is applied to the filtered signal. Its
derivatives are chain states where those exist, and the rate of the highest
active state otherwise.

```math
\begin{aligned}
0 &= -v_4 + \begin{cases}
      u & n = 0 \\
      x_1 + \dfrac{A_5}{a_1}\left(-a_0 x_1 + u\right) & n = 1 \\
      x_1 + A_5 x_2 + \dfrac{A_6}{a_2}\left(-a_0 x_1 - a_1 x_2 + u\right) & n = 2 \\
      x_1 + A_5 x_2 + A_6 x_3 & n \in \{3,4\}
    \end{cases} \\
0 &= -T_2(v_5 - x_5) + T_1(v_4 - x_5) \\
0 &= -T_4(v_6 - x_6) + T_3(v_5 - x_6) \\
0 &= -T_6 v_7 + K_s T_5(v_6 - x_7) \\
0 &= -V_{ss} + \text{clamp}(v_7, L_s^{\min}, L_s^{\max})
\end{aligned}
```

The output limiter uses GridKit's smooth
[Clamp](../../../../CommonMath.md#derived-functions).

## Initialization

### Input Initialization

```math
\begin{aligned}
  u_0 &\leftarrow \text{stabilizer input signal}
\end{aligned}
```

### Internal Initialization

The notch filter and both lead–lag blocks have unity DC gain, so every block
output settles at $u_0$ while the washout removes the DC component. The
initial residual therefore vanishes for any $u_0$, and the stabilizer comes
online at rest.

```math
\begin{aligned}
  x_1 &\leftarrow u_0
    \quad n\ge1 \\
  x_i &\leftarrow 0
    \quad i=2,\ldots,4 \\
  v_4, x_5, v_5, x_6, v_6, x_7 &\leftarrow u_0 \\
  v_7 &\leftarrow 0 \\
  V_{ss} &\leftarrow \text{clamp}(v_7, L_s^{\min}, L_s^{\max}) \\
  \dot{x}_i &\leftarrow 0
\end{aligned}
```

### Output Initialization

None.

## Monitorable Outputs

Output | Units  | Description
-------|--------|------------
`vss`  | [p.u.] | Limited stabilizer signal $V_{ss}$
