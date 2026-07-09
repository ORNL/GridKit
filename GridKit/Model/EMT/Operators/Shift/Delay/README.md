# Delay Model

`Delay` represents a scalar EMT delay operator using a lag-block chain.
The model maps input signal $u$ to delayed output $y_{\mathrm{out}}$.

Note:
- This is exact with forward Euler when the integration step satisfies $\Delta t=T$.
- For other integration methods or time steps, this is a smooth approximation only.

## Block Diagram

![](../../../../../../docs/Figures/EMT/Delay/diagram.png)

Figure 1: Delay model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\tau$ | [s] | `tau` | Total delay | Required, positive
$f_{\max}$ | [Hz] | `fmax` | Highest frequency of interest | Required, positive

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
$\mathbf{y}$ | [-] | Section differential states | $\mathbf{y}\in\mathbb{R}^J$

#### Algebraic

None.

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$u$ | [-] | Input signal | $u\in\mathbb{R}$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$u$ | `input` | Input | [-] | Input signal port | $u\in\mathbb{R}$
$y_{\mathrm{out}}$ | `out` | Output | [-] | Delayed output port | $y_{\mathrm{out}} = y_{J-1}$

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= -T\dot{y}_0 - y_0 + u \\
0 &= -T\dot{y}_n - y_n + y_{n-1},
     \quad n\in\{1,\ldots,J-1\}
\end{aligned}
```

### Algebraic Equations

None.

### Wiring

None.

## Initialization

### Input Initialization

```math
\begin{aligned}
u,\dot{u}
  &\leftarrow \text{affine input trajectory start}
\end{aligned}
```

### Internal Initialization

Initialization is performed by evaluating the affine-input residuals in
dependency order:

```math
\begin{aligned}
y_n
  &\leftarrow
  u - (n+1)T\dot{u},
     \quad n\in\{0,\ldots,J-1\} \\
\dot{y}_n
  &\leftarrow
  \dot{u},
     \quad n\in\{0,\ldots,J-1\}
\end{aligned}
```

### Output Initialization

```math
y_{\mathrm{out}} \leftarrow u - \tau\dot{u}
```

## Monitors

None.
