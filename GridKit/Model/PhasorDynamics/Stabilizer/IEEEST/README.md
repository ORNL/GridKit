# **IEEE Stabilizer Model (IEEEST)**

IEEEST is a standard IEEE power system stabilizer with a derived-order notch
filter, two lead-lag blocks, washout, and an output limiter.

## Block Diagram

![](../../../../../docs/Figures/stabilizer_ieeest_diagram.png)

Figure 1: Stabilizer IEEEST model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol               | Units    | JSON     | Description                          | Typical Value
---------------------|----------|----------|--------------------------------------|--------------
$A_1$                | [sec]    | `A1`     | Notch denominator coefficient        | 1.013
$A_2$                | [sec²]   | `A2`     | Notch denominator coefficient        | 0.013
$A_3$                | [sec]    | `A3`     | Notch denominator coefficient        | 0.0
$A_4$                | [sec²]   | `A4`     | Notch denominator coefficient        | 0.0
$A_5$                | [sec]    | `A5`     | Notch numerator coefficient          | 1.013
$A_6$                | [sec²]   | `A6`     | Notch numerator coefficient          | 0.113
$T_1$                | [sec]    | `T1`     | Lead-lag 1 numerator time constant   | 0.0
$T_2$                | [sec]    | `T2`     | Lead-lag 1 denominator time constant | 0.02
$T_3$                | [sec]    | `T3`     | Lead-lag 2 numerator time constant   | 0.0
$T_4$                | [sec]    | `T4`     | Lead-lag 2 denominator time constant | 0.0
$T_5$                | [sec]    | `T5`     | Washout numerator time constant      | 1.65
$T_6$                | [sec]    | `T6`     | Washout denominator time constant    | 1.65
$K_s$                | [p.u.]   | `Ks`     | Stabilizer gain                      | 3.0
$L_s^{\min}$         | [p.u.]   | `Lsmin`  | Minimum stabilizer output limit      | -0.1
$L_s^{\max}$         | [p.u.]   | `Lsmax`  | Maximum stabilizer output limit      | 0.1
$V_{\mathrm{cl}}$    | [p.u.]   | `Vcl`    | Lower input cutout threshold         | 0.0
$V_{\mathrm{cu}}$    | [p.u.]   | `Vcu`    | Upper input cutout threshold         | 0.0
$T_{\mathrm{delay}}$ | [sec]    | `Tdelay` | Input delay                          | 0.0

$V_{\mathrm{cl}}$, $V_{\mathrm{cu}}$, and $T_{\mathrm{delay}}$ are accepted for
input-format compatibility but are not modeled. A nonzero value logs a warning
and is ignored.

`IeeestFactory` derives the notch-filter order from the denominator coefficients
and selects `Ieeest<ScalarT, IdxT, order>`. Each order has its own residual and
Jacobian, with only its active notch-filter states allocated. Order is not an
input parameter.

### Parameter Validation

```math
\begin{aligned}
  n &\in \{0,1,2,3,4\} \\
  m &\le n \\
  a_n &\ne 0, \qquad n>0 \\
  T_2,T_4,T_6 &\ge 0 \\
  L_s^{\min} &< L_s^{\max}.
\end{aligned}
```

All numeric parameters and the expanded denominator coefficients must be finite.
The reciprocal of the active leading coefficient, normalized denominator
coefficients, and notch, lead-lag, and washout output coefficients must also be
finite. An underflowed leading coefficient is rejected without changing the
derived order. Denominator reciprocals are computed before residual evaluation
to avoid overflow in differentiated quotients.
The numerator order $m$ cannot exceed the denominator order $n$. Thus $A_5$ and
$A_6$ must be zero when $n=0$, and $A_6$ must be zero when $n=1$. Direct
construction of a template specialization must match the derived order.

### Model Derived Parameters

For one polynomial factor, define its degree using exact zero comparisons:

```math
d(p,q)=
\begin{cases}
  2 & q\ne0 \\
  1 & q=0,\ p\ne0 \\
  0 & p=q=0.
\end{cases}
```

The factory derives the order from the two denominator factors. The orders and
expanded denominator coefficients are

```math
\begin{aligned}
  n   &= d(A_1,A_2)+d(A_3,A_4) \\
  m   &= d(A_5,A_6) \\
  a_1 &= A_1 + A_3 \\
  a_2 &= A_2 + A_4 + A_1 A_3 \\
  a_3 &= A_1 A_4 + A_2 A_3 \\
  a_4 &= A_2 A_4.
\end{aligned}
```

The denominator is $1+a_1s+\cdots+a_ns^n$, and the numerator is
$1+A_5s+A_6s^2$. Order selection does not apply a numerical zero tolerance.
At order zero the notch filter is bypassed.

Let $\epsilon_T=10^{-3}\ \mathrm{s}$. Accepted denominator time constants below
$\epsilon_T$, including zero, are raised to that floor. A warning is logged once
per template specialization:

```math
T \leftarrow \max(T,\epsilon_T),
\qquad T\in\{T_2,T_4,T_6\}.
```

This minimum time constant approximates a fast response; it does not implement
an exact zero-time-constant bypass.

## Model Ports

Name     | Port   | Init  | Description
---------|--------|-------|------------
`input`  | Input  | Known | Stabilizer input signal
`output` | Output | Known | Stabilizer output signal

The input must be attached to a linked signal. Assigning the output signal is
optional.

## Model Variables

### Internal Variables

#### Differential

There are $n+3$ differential variables and five algebraic variables. Only the
first $n$ notch-filter states are present; $x_5$, $x_6$, and $x_7$ retain their
block-diagram names for every order.

Symbol | Units       | Description                         | Note
-------|-------------|-------------------------------------|-----
$x_1$  | [p.u.]      | Notch-filter signal state           | Present for $n\ge1$
$x_2$  | [p.u./sec]  | First derivative of filtered signal | Present for $n\ge2$
$x_3$  | [p.u./sec²] | Second derivative of filtered signal | Present for $n\ge3$
$x_4$  | [p.u./sec³] | Third derivative of filtered signal | Present for $n=4$
$x_5$  | [p.u.]      | Lead-lag 1 state                    |
$x_6$  | [p.u.]      | Lead-lag 2 state                    |
$x_7$  | [p.u.]      | Washout state                       |

#### Algebraic

Symbol            | Units  | Description
------------------|--------|------------
$v_4$             | [p.u.] | Notch filter output
$v_5$             | [p.u.] | Lead-lag 1 output
$v_6$             | [p.u.] | Lead-lag 2 output
$v_7$             | [p.u.] | Unlimited stabilizer signal
$V_{\mathrm{ss}}$ | [p.u.] | Limited stabilizer signal (model output)

### External Variables

#### Algebraic

Symbol | Units  | Description
-------|--------|------------
$u$    | [p.u.] | Stabilizer input signal

## Model Equations

### Differential Equations

For $n>0$, the notch filter satisfies

```math
\begin{aligned}
  0 &= -\dot{x}_i + x_{i+1}, \qquad 1\le i<n \\
  0 &= -\dot{x}_n
    + \dfrac{u-x_1-\sum_{i=1}^{n-1}a_i x_{i+1}}{a_n}.
\end{aligned}
```

There are no notch-filter differential equations at order zero. For every
order, the lead-lag and washout states satisfy

```math
\begin{aligned}
  0 &= -\dot{x}_5 + \dfrac{v_4-x_5}{T_2} \\
  0 &= -\dot{x}_6 + \dfrac{v_5-x_6}{T_4} \\
  0 &= -\dot{x}_7 + \dfrac{v_6-x_7}{T_6}.
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
  0 &= -v_4
    + \begin{cases}
      u & n=0 \\
      x_1 + \dfrac{A_5}{a_1}(u-x_1) & n=1 \\
      x_1 + A_5x_2 + \dfrac{A_6}{a_2}(u-x_1-a_1x_2) & n=2 \\
      x_1 + A_5x_2 + A_6x_3 & n\in\{3,4\}
    \end{cases} \\
  0 &= -v_5 + x_5 + \dfrac{T_1}{T_2}(v_4-x_5) \\
  0 &= -v_6 + x_6 + \dfrac{T_3}{T_4}(v_5-x_6) \\
  0 &= -v_7 + K_s\dfrac{T_5}{T_6}(v_6-x_7) \\
  0 &= -V_{\mathrm{ss}}
    + \text{clamp}(v_7,L_s^{\min},L_s^{\max}).
\end{aligned}
```

The output limiter uses GridKit's smooth
[clamp](../../../../CommonMath.md#clamp).

## Initialization

### Input Initialization

Read the finite initial value $u$ from the linked input signal.

### Internal Initialization

For $n>0$, initialize $x_1=u$ and the remaining active notch-filter states to
zero. For every order,

```math
\begin{aligned}
  x_5,x_6,x_7,v_4,v_5,v_6 &\leftarrow u \\
  v_7 &\leftarrow 0 \\
  V_{\mathrm{ss}} &\leftarrow \text{clamp}(0,L_s^{\min},L_s^{\max}).
\end{aligned}
```

All state derivatives initialize to zero. This is a steady state for a constant
input, including a nonzero input.

### Output Initialization

The optional `output` signal is linked to $V_{\mathrm{ss}}$.

## Monitorable Outputs

Output | Units  | Description
-------|--------|------------
`vss`  | [p.u.] | Limited stabilizer signal $V_{\mathrm{ss}}$
