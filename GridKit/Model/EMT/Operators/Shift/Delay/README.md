# Delay Model

For input units $[u]$, `Delay` maps scalar input $u$ to delayed output
$y_{\mathrm{out}}$ through a chain of first-order lag sections.

> [!WARNING]
> The lag chain is an exact sampled $J$-step delay under forward Euler when
> $h = T$.
> Otherwise it is an approximation. The `fmax` parameter controls section
> density; it does not guarantee accuracy over a signal bandwidth.

## Block Diagram

![Delay operator block diagram](../../../../../../docs/Figures/EMT/Delay/diagram.png)

Figure 1: Delay model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\tau$ | [s] | `tau` | Total delay | Required, positive
$f_{\max}$ | [Hz] | `fmax` | Lag-chain section rate | Required, positive

### Parameter Validation

```math
\begin{aligned}
\tau &> 0 \\
f_{\max} &> 0
\end{aligned}
```

### Derived Parameters

```math
\begin{aligned}
J &= \left\lceil f_{\max}\tau \right\rceil \\
T &= \dfrac{\tau}{J}
\end{aligned}
```

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{y}$ | $[u]$ | Section differential states | $\mathbf{y} \in \mathbb{R}^J$

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$u$ | $[u]$ | Input signal | $u \in \mathbb{R}$

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$u$ | `input` | Input | $[u]$ | Input signal port | $u \in \mathbb{R}$
$y_{\mathrm{out}}$ | `out` | Output | $[u]$ | Delayed output port | $y_J$

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= -T\dfrac{\mathrm{d}y_1}{\mathrm{d}t} - y_1 + u \\
0 &= -T\dfrac{\mathrm{d}y_j}{\mathrm{d}t} - y_j + y_{j-1},
     \quad j \in \{2,\ldots,J\}
\end{aligned}
```

### Algebraic Equations

None.

### Wiring

```math
y_{\mathrm{out}} = y_J
```

## Initialization

Initialization assumes an affine input with zero second derivative.

### Input Initialization

```math
\begin{aligned}
u,\dfrac{\mathrm{d}u}{\mathrm{d}t}
  &\leftarrow \text{input value and derivative}
\end{aligned}
```

### Internal Initialization

```math
\begin{aligned}
y_j
  &\leftarrow
  u - jT\dfrac{\mathrm{d}u}{\mathrm{d}t},
     \quad j \in \{1,\ldots,J\} \\
\dfrac{\mathrm{d}y_j}{\mathrm{d}t}
  &\leftarrow
  \dfrac{\mathrm{d}u}{\mathrm{d}t},
     \quad j \in \{1,\ldots,J\}
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
y_{\mathrm{out}}
  &\leftarrow u - \tau\dfrac{\mathrm{d}u}{\mathrm{d}t} \\
\dfrac{\mathrm{d}y_{\mathrm{out}}}{\mathrm{d}t}
  &\leftarrow \dfrac{\mathrm{d}u}{\mathrm{d}t}
\end{aligned}
```

## Monitors

None.
