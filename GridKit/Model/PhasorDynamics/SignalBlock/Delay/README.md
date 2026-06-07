# **Delay Model**

The Delay model writes a time-shifted copy of its input signal to its output
signal: the output at time $t$ is the input as it was a fixed lag $\tau$ earlier.
The model has no solver-owned state or residual equations; it records its input
at accepted solver steps and emits the interpolated past value.

Notes:
- Input history is linearly interpolated between accepted steps, so a path
  through a `Delay` is first-order accurate in the step size.
- The solver caps its step at $\tau$ and records history each accepted step, so
  the delayed lookup always falls inside committed history.
- A `Delay` placed in a feedback path breaks the algebraic loop: its output is
  explicit forcing within a step, not coupled to the current unknowns.

## Model Parameters

Symbol      | Units | JSON          | Description                                   | Typical Value
------------|-------|---------------|-----------------------------------------------|--------------
$\tau$      | s     | `delay`       | Transport delay (required, must be positive)  | 0.25
$\varphi$   | [-]   | `prehistory`  | Constant output for $t < t_0 + \tau$ (optional) | input at $t_0$

If `prehistory` is omitted, the output before enough history exists holds the
input value sampled at the initial condition.

## Model Variables

### Internal Variables

#### Differential

None.

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

None.

### Algebraic Equations

None. The output is computed directly from recorded input history:

```math
o(t) =
\begin{cases}
  u(t-\tau), & t-\tau \ge t_0 \\
  \varphi,   & t-\tau <   t_0
\end{cases}
```

where $\varphi$ is the `prehistory` value (defaulting to $u(t_0)$).

## Initialization

The output initializes to the prehistory value, so a delayed steady input
passes through unchanged at $t_0$ and downstream consumers initialize
consistently.

## Model Outputs

Output | Units | Description
-------|-------|------------
`out`  | [-]   | Input signal delayed by $\tau$
