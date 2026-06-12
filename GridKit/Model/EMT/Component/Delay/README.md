# Delay Model

`Delay` represents a smooth approximation of a transport delay on a scalar input
signal. The approximation uses a chain of $n$ identical first-order lag stages.

The Laplace domain representation is:

```math
e^{-s\tau}U(s) =
\lim_{n \to \infty}\left(\dfrac{1}{1 + s\tau/n}\right)^n U(s)
```

The time domain convolutional form is:

```math
u(t-\tau) = \delta(t-\tau) * u(t)
```

The model owns $n$ differential states and keeps no input history.

## Block Diagram

<div align="center">
   <img align="center" src="../../../../../docs/Figures/EMT/Delay/diagram.png">

  Figure 1: Lag-chain approximation of a pure delay
</div>

## Model Parameters

Symbol            | Units | JSON     | Description                                          | Typical Value
------------------|-------|----------|------------------------------------------------------|--------------
$\tau$            | s     | `delay`  | Delay to approximate (required, must be positive)    | --
$\Delta t_{\min}$ | s     | `dt_min` | Block resolution (required, must be positive)        | --

### Parameter Validation

```math
\begin{aligned}
\tau &> 0 \\
\Delta t_{\min} &> 0
\end{aligned}
```

## Model Derived Parameters

```math
n = \operatorname{floor}\left(\dfrac{\tau}{\Delta t_{\min}}\right)
```

## Model Variables

### Internal Variables

#### Differential

Symbol          | Units | Qty | Description
----------------|-------|-----|------------
$x_1,\dots,x_n$ | [-]   | $n$ | Lag-block states; the last state $x_n$ is the model output

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Qty | Description
-------|-------|-----|------------
$u$    | [-]   | 1   | Input signal read from the `input` port

## Model Equations

### Differential Equations

The lag-chain residuals are:

```math
\begin{aligned}
0 &= -\tau\,\dot{x}_1 + n\,(u - x_1) \\
0 &= -\tau\,\dot{x}_2 + n\,(x_1 - x_2) \\
&\vdots \\
0 &= -\tau\,\dot{x}_n + n\,(x_{n-1} - x_n)
\end{aligned}
```

### Algebraic Equations

None.

## Initialization

For a constant input $u_0$ at $t_0$, the chain is at rest:

```math
\begin{aligned}
x_1(t_0) = x_2(t_0) = \cdots = x_n(t_0) &= u_0 \\
\dot{x}_1(t_0) = \dot{x}_2(t_0) = \cdots = \dot{x}_n(t_0) &= 0
\end{aligned}
```

A steady input therefore passes through unchanged at $t_0$ and downstream
consumers initialize consistently.

## Model Outputs

Output | Units | Description
-------|-------|------------
`out`  | [-]   | Input signal smoothly delayed by $\tau$
