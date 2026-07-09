# Propagation Model

`Propagation` represents the current-form EMT propagation operator used by
`LineDistributed`. It composes fitted input and output factors with one scalar
delay per mode.

```math
\mathbf{y}
=
f^{\mathbf{h}_{\mathrm{out}}}\left(
  \boldsymbol{\delta}_{\tau}\left(
    f^{\mathbf{h}_{\mathrm{in}}}(\mathbf{u})
  \right)
\right)
```

where, in the frequency domain,

```math
\begin{aligned}
f^{\mathbf{h}_{\mathrm{in}}}(s) &\approx
\widehat{\mathbf{H}}^{\mathrm{mps}}(s)\mathbf{T}_v(s)^{\mathsf H} \\
\boldsymbol{\delta}_{\tau}(s) &=
\mathrm{diag}\left(e^{-s\tau_1},\ldots,e^{-s\tau_M}\right) \\
f^{\mathbf{h}_{\mathrm{out}}}(s) &\approx
\mathbf{T}_i(s).
\end{aligned}
```

## Block Diagram

![](../../../../../../docs/Figures/EMT/Propagation/diagram.png)

Figure 1: Propagation model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$K$ | [-] | `K` | Number of conductors | Required, positive integer
$\boldsymbol{\tau}$ | [s] | `tau` | Modal propagation delays | $\boldsymbol{\tau}\in\mathbb{R}^M$
$f_{\max}$ | [Hz] | `fmax` | Delay highest frequency of interest | Required, positive

### Parameter Validation

```math
\begin{aligned}
K &> 0 \\
\boldsymbol{\tau} &\in \mathbb{R}^M \\
\tau_m &> 0,\quad m\in\{1,\ldots,M\} \\
f_{\max} &> 0
\end{aligned}
```

### Derived Parameters

```math
\begin{aligned}
J_m &= \mathrm{ceil}(f_{\max}\tau_m),\quad m\in\{1,\ldots,M\} \\
Q_{\mathbf{h}} &=
  Q_{\mathbf{h}_{\mathrm{in}}}
  + \sum_{m=1}^{M}J_m
  + Q_{\mathbf{h}_{\mathrm{out}}}
\end{aligned}
```

## Submodels

Submodel | Type | Order | Inputs | Outputs
-------- | ---- | ----- | ------ | -------
Input factor $f^{\mathbf{h}_{\mathrm{in}}}$ | [`VectorFit`](../../Rational/VectorFit/README.md) | $Q_{\mathbf{h}_{\mathrm{in}}}$ | $\mathbf{u}\in\mathbb{R}^K$ | $\mathbf{w}\in\mathbb{R}^M$
Modal delay $\delta_{\tau_m}$ | [`Delay`](../Delay/README.md) | $J_m$ | $w_m\in\mathbb{R}$ | $z_m\in\mathbb{R}$
Output factor $f^{\mathbf{h}_{\mathrm{out}}}$ | [`VectorFit`](../../Rational/VectorFit/README.md) | $Q_{\mathbf{h}_{\mathrm{out}}}$ | $\mathbf{z}\in\mathbb{R}^M$ | $\mathbf{y}\in\mathbb{R}^K$

### Submodel Validation

```math
\begin{aligned}
f^{\mathbf{h}_{\mathrm{in}}} &: \mathbb{R}^K \rightarrow \mathbb{R}^M \\
\delta_{\tau_m} &: \mathbb{R} \rightarrow \mathbb{R},\quad m\in\{1,\ldots,M\} \\
f^{\mathbf{h}_{\mathrm{out}}} &: \mathbb{R}^M \rightarrow \mathbb{R}^K
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

None.

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | [-] | Input vector | $\mathbf{u}\in\mathbb{R}^K$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | [-] | Input vector port | $\mathbf{u}\in\mathbb{R}^K$
$\mathbf{y}$ | `out` | Output | [-] | Output contribution port | $\mathbf{y}\in\mathbb{R}^K$

## Model Equations

### Differential Equations

None.

### Algebraic Equations

None.

### Wiring

```math
\begin{aligned}
\mathbf{w} &= f^{\mathbf{h}_{\mathrm{in}}}(\mathbf{u}) \\
z_m &= \delta_{\tau_m}(w_m),
     \quad m\in\{1,\ldots,M\} \\
\mathbf{z} &= \boldsymbol{\delta}_{\tau}(\mathbf{w}) \\
\mathbf{y} &= f^{\mathbf{h}_{\mathrm{out}}}(\mathbf{z}).
\end{aligned}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{u}
  &\leftarrow \text{input vector start}
\end{aligned}
```

### Internal Initialization

Initialization is performed by evaluating the submodels in dependency order.

```math
\begin{aligned}
\mathbf{w}
  &\leftarrow f^{\mathbf{h}_{\mathrm{in}}}(\mathbf{u}) \\
z_m
  &\leftarrow \delta_{\tau_m}(w_m),
     \quad m\in\{1,\ldots,M\} \\
\mathbf{z}
  &\leftarrow \boldsymbol{\delta}_{\tau}(\mathbf{w})
\end{aligned}
```

### Output Initialization

```math
\mathbf{y} \leftarrow f^{\mathbf{h}_{\mathrm{out}}}(\mathbf{z})
```

## Monitors

None.
