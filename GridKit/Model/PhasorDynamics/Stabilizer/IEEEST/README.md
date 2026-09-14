# IEEEST

Standard IEEE power system stabilizer: notch filter of order up to four, two lead–lag
blocks, washout, and output limiter.

## Block Diagram

![](../../../../../docs/Figures/stabilizer_ieeest_diagram.png)

Figure 1: Stabilizer IEEEST model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol       | Units  | JSON    | Description                          | Typical Value | Note
-------------|--------|---------|--------------------------------------|---------------|-----
$A_1$        | [s]    | `A1`    | Notch denominator coefficient        | 1.013         |
$A_2$        | [s²]   | `A2`    | Notch denominator coefficient        | 0.013         |
$A_3$        | [s]    | `A3`    | Notch denominator coefficient        | 0.0           |
$A_4$        | [s²]   | `A4`    | Notch denominator coefficient        | 0.0           |
$A_5$        | [s]    | `A5`    | Notch numerator coefficient          | 1.013         |
$A_6$        | [s²]   | `A6`    | Notch numerator coefficient          | 0.113         |
$T_1$        | [s]    | `T1`    | Lead–lag 1 numerator time constant   | 0.0           |
$T_2$        | [s]    | `T2`    | Lead–lag 1 denominator time constant | 0.02          |
$T_3$        | [s]    | `T3`    | Lead–lag 2 numerator time constant   | 0.0           |
$T_4$        | [s]    | `T4`    | Lead–lag 2 denominator time constant | 0.0           |
$T_5$        | [s]    | `T5`    | Washout numerator time constant      | 1.65          |
$T_6$        | [s]    | `T6`    | Washout denominator time constant    | 1.65          |
$K_\mathrm{s}$        | [p.u.] | `Ks`    | Stabilizer gain                      | 3.0           |
$L_\mathrm{s}^{\min}$ | [p.u.] | `Lsmin` | Minimum stabilizer output limit      | -0.1          |
$L_\mathrm{s}^{\max}$ | [p.u.] | `Lsmax` | Maximum stabilizer output limit      | 0.1           |

The IEEE 421.5 IEEEST also defines a cutout window ($V_\mathrm{cl}$, $V_\mathrm{cu}$) and an
input delay ($T_\mathrm{delay}$). These parameters are accepted for input-format
compatibility but are not modeled here.

### Parameter Validation

A valid IEEEST parameter set must satisfy the following conditions:

```math
a_2=a_3=a_4=0 \quad\Longrightarrow\quad a_1=0
```

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

The notch order is

```math
n = \begin{cases}
4 & a_4\ne0 \\
3 & a_4=0,\ a_3\ne0 \\
2 & a_4=a_3=0,\ a_2\ne0 \\
0 & a_4=a_3=a_2=a_1=0
\end{cases}
```

## Model Ports

Name     | Port   | Init  | Description
---------|--------|-------|---------------------------------
`input`  | Input  | Known | Required stabilizer input signal
`output` | Output | Known | Limited stabilizer output signal

## Model Variables

### Internal Variables

#### Differential

Symbol | Units     | Description        | Note
-------|-----------|--------------------|------------------------
$x_1$  | [p.u.]    | Notch filter state | Held at zero when $n=0$
$x_2$  | [p.u./s]  | Notch filter state | Held at zero when $n=0$
$x_3$  | [p.u./s²] | Notch filter state | Held at zero when $n<3$
$x_4$  | [p.u./s³] | Notch filter state | Held at zero when $n<4$
$x_5$  | [p.u.]    | Lead–lag 1 state   | Algebraic when $T_2=0$
$x_6$  | [p.u.]    | Lead–lag 2 state   | Algebraic when $T_4=0$
$x_7$  | [p.u.]    | Washout state      | Algebraic when $T_6=0$

#### Algebraic

Symbol            | Units  | Description                              | Note
------------------|--------|------------------------------------------|-----
$v_4$             | [p.u.] | Notch filter output                      |
$v_5$             | [p.u.] | Lead–lag 1 output                        |
$v_6$             | [p.u.] | Lead–lag 2 output                        |
$v_7$             | [p.u.] | Unlimited stabilizer signal              |
$V_{\mathrm{ss}}$ | [p.u.] | Limited stabilizer signal (model output) |

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description             | Note
-------|--------|-------------------------|-----
$u$    | [p.u.] | Stabilizer input signal |

## Model Equations

Smooth functions: [`clamp`](../../../../CommonMath.md#clamp).

### Internal Equations

#### Differential

```math
\begin{aligned}
0 &= -\dot{x}_k+x_{k+1}, && 1\le k<n \\
0 &= -a_n\dot{x}_n-\sum_{j=0}^{n-1}a_jx_{j+1}+u, && n>0 \\
0 &= -\dot{x}_k, && n<k\le4 \\
0 &= -T_2\dot{x}_5-x_5+v_4 \\
0 &= -T_4\dot{x}_6-x_6+v_5 \\
0 &= -T_6\dot{x}_7-x_7+v_6
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
0 &= -v_4+\begin{cases}
x_1+A_5x_2+A_6\dot{x}_2 & n>0 \\
u & n=0
\end{cases} \\
0 &= \begin{cases}
-T_2(v_5-x_5)+T_1(v_4-x_5) & T_2\ne0 \\
-v_5+v_4 & T_2=0
\end{cases} \\
0 &= \begin{cases}
-T_4(v_6-x_6)+T_3(v_5-x_6) & T_4\ne0 \\
-v_6+v_5 & T_4=0
\end{cases} \\
0 &= \begin{cases}
-T_6v_7+K_\mathrm{s}T_5(v_6-x_7) & T_6\ne0 \\
-v_7+K_\mathrm{s}v_6 & T_6=0
\end{cases} \\
0 &= -V_{\mathrm{ss}}+\text{clamp}(v_7;L_\mathrm{s}^{\min},L_\mathrm{s}^{\max})
\end{aligned}
```

### External Equations

None.

## Initialization

The initial input determines the steady state; all derivatives initialize to zero.

```math
\begin{aligned}
x_1 &\leftarrow \begin{cases}u & n>0 \\ 0 & n=0\end{cases} \\
x_2,x_3,x_4 &\leftarrow 0 \\
x_5,x_6,x_7,v_4,v_5,v_6 &\leftarrow u \\
v_7 &\leftarrow \begin{cases}0 & T_6\ne0 \\ K_\mathrm{s}u & T_6=0\end{cases} \\
V_{\mathrm{ss}} &\leftarrow \text{clamp}(v_7;L_\mathrm{s}^{\min},L_\mathrm{s}^{\max})
\end{aligned}
```

## Monitors

Monitor | Units  | Description               | Note
--------|--------|---------------------------|--------------------------------
`vss`   | [p.u.] | Limited stabilizer signal | $V_{\mathrm{ss}}$; model output
