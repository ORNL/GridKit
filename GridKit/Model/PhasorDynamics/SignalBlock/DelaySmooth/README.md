# **DelaySmooth Model**

The DelaySmooth model approximates a transport delay $\tau$ on its input signal
by a string of $N$ first-order lag blocks, each with time constant $T = \tau/N$.
As $N$ grows, the cascade transfer function $\left[1/(1+sT)\right]^N$
approaches the ideal delay $e^{-s\tau}$. Unlike the
[Delay](../Delay/README.md) model, the delayed value is carried in solver
states, so the model owns a residual and contributes to the DAE.

Notes:
- The delay is realized entirely in differential states; the model keeps no
  input history and imposes no maximum step size.
- The approximation is smooth and strictly causal. Its step response is
  monotone with no overshoot, and its group delay at low frequency is exactly
  $\tau$.
- Accuracy increases with $N$ (smaller $\Delta t_{\min}$) at the cost of more
  states.

## Model Parameters

Symbol            | Units | JSON     | Description                                          | Typical Value
------------------|-------|----------|------------------------------------------------------|--------------
$\tau$            | s     | `delay`  | Delay to approximate (required, must be positive)    | 0.25
$\Delta t_{\min}$ | s     | `dt_min` | Block resolution (required, must be positive)        | 0.03125

## Model Derived Parameters

The matrix $\mathbf{A}$ is the incidence matrix of a directed path: $-1$ on the
diagonal and $+1$ on the first sub-diagonal.

```math
\begin{aligned}
N &= \text{floor}\left(\dfrac{\tau}{\Delta t_{\min}}\right) \\
T &= \dfrac{\tau}{N} \\
G &= \dfrac{1}{T} = \dfrac{N}{\tau}
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol          | Units | Description
----------------|-------|------------
$x_1,\dots,x_N$ | [-]   | Lag-block states; the last state $x_N$ is the model output

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description
-------|-------|------------
$u$    | [-]   | Input signal read from the `input` port

## Model Equations

### Differential Equations

With $x_0 \equiv u$, every block obeys the same first-order lag:

```math
\begin{aligned}
0 &= -\dot{x}_1 + G\,(u - x_1) \\
0 &= -\dot{x}_k + G\,(x_{k-1} - x_k), \quad k = 2,\dots,N
\end{aligned}
```

### Algebraic Equations

None.

## Initialization

For a constant input $u_0$ at $t_0$, the chain is at rest:

```math
x_k(t_0) = u_0, \qquad \dot{x}_k(t_0) = 0, \qquad k = 1,\dots,N
```

A steady input therefore passes through unchanged at $t_0$ and downstream
consumers initialize consistently.

## Model Outputs

Output | Units | Description
-------|-------|------------
`out`  | [-]   | Input signal smoothly delayed by $\tau$
