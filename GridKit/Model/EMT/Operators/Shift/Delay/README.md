# Delay Model

For input units $[u]$, `Delay` represents a scalar EMT delay operator using a
lag-block chain.
The model maps input signal $u$ to delayed output $y_{\mathrm{out}}$ and
preserves the input signal units.

> [!WARNING]
> The lag chain is an exact sampled $J$-step delay under forward Euler only
> when the integration step satisfies $h = T$. Otherwise it is a smooth
> approximation. The `fmax` parameter controls section density and does not
> guarantee delay accuracy over a signal bandwidth.

## Block Diagram

![Delay operator block diagram](../../../../../../docs/Figures/EMT/Delay/diagram.png)

Figure 1: Delay model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\tau$ | [s] | `tau` | Total delay | Required, positive
$f_{\max}$ | [Hz] | `fmax` | Lag-chain section rate | Required, positive. Not a bandwidth guarantee

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
J &= \mathrm{ceil}(f_{\max}\tau) \\
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

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$u$ | $[u]$ | Input signal | $u \in \mathbb{R}$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$u$ | `input` | Input | $[u]$ | Input signal port | $u \in \mathbb{R}$
$y_{\mathrm{out}}$ | `out` | Output | $[u]$ | Delayed output port | Final section state

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= -T\dfrac{\mathrm{d}y_0}{\mathrm{d}t} - y_0 + u \\
0 &= -T\dfrac{\mathrm{d}y_n}{\mathrm{d}t} - y_n + y_{n-1},
     \quad n \in \{1,\ldots,J-1\}
\end{aligned}
```

### Algebraic Equations

None.

### Wiring

```math
y_{\mathrm{out}} = y_{J-1}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
u,\dfrac{\mathrm{d}u}{\mathrm{d}t}
  &\leftarrow \text{affine input trajectory start}
\end{aligned}
```

### Internal Initialization

The section states use the zero-transient particular trajectory for the affine
input:

```math
\begin{aligned}
y_n
  &\leftarrow
  u - (n+1)T\dfrac{\mathrm{d}u}{\mathrm{d}t},
     \quad n \in \{0,\ldots,J-1\} \\
\dfrac{\mathrm{d}y_n}{\mathrm{d}t}
  &\leftarrow
  \dfrac{\mathrm{d}u}{\mathrm{d}t},
     \quad n \in \{0,\ldots,J-1\}
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
