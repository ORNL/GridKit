# Propagation Model

`Propagation` represents the current-form EMT propagation operator used by
`LineDistributed`. It composes two fitted `VectorFit` factors with scalar
`Delay` blocks, one per propagation mode.

```math
\mathbf{H}_i(s)
=
\mathbf{F}_{\mathrm{out}}(s)
\boldsymbol{\Delta}_{\tau}(s)
\mathbf{F}_{\mathrm{in}}(s)
```

where

```math
\begin{aligned}
\mathbf{F}_{\mathrm{in}}(s) &\approx
\widehat{\mathbf{H}}^{\mathrm{mps}}(s)\mathbf{T}_v(s)^{\mathsf H} \\
\boldsymbol{\Delta}_{\tau}(s) &=
\operatorname{diag}\left(e^{-s\tau_1},\ldots,e^{-s\tau_M}\right) \\
\mathbf{F}_{\mathrm{out}}(s) &\approx
\mathbf{T}_i(s).
\end{aligned}
```

The fitted factors use the full-residue real-port
[`VectorFit`](../../Rational/VectorFit/README.md) form.

## Model Parameters

For conductor count $K$ and modal count $M$:

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\mathbf{F}_{\mathrm{in}}$ | [-] | `input` | Input-side fitted factor | `VectorFit`, $M\times K$
$\boldsymbol{\tau}$ | [s] | `tau` | Modal propagation delays | $\boldsymbol{\tau}\in\mathbb{R}^M$
$f_{\max}$ | [Hz] | `fmax` | Delay highest frequency of interest | passed to each scalar `Delay`
$\mathbf{F}_{\mathrm{out}}$ | [-] | `output` | Output-side fitted factor | `VectorFit`, $K\times M$

### Parameter Validation

```math
\begin{aligned}
\mathbf{F}_{\mathrm{in}} &: M \times K \\
\boldsymbol{\tau} &\in \mathbb{R}^M \\
\mathbf{F}_{\mathrm{out}} &: K \times M
\end{aligned}
```

### Model Derived Parameters

None.

### Model Submodels

Submodel | Inputs | Outputs
-------- | ------ | -------
[`VectorFit`](../../Rational/VectorFit/README.md) $\mathbf{F}_{\mathrm{in}}$ | $\mathbf{u}\in\mathbb{R}^K$ | $\mathbf{w}\in\mathbb{R}^M$
[`Delay`](../Delay/README.md) $\boldsymbol{\delta}_{\tau}$ | $\mathbf{w}\in\mathbb{R}^M$, $\boldsymbol{\tau}$ | $\mathbf{z}\in\mathbb{R}^M$
[`VectorFit`](../../Rational/VectorFit/README.md) $\mathbf{F}_{\mathrm{out}}$ | $\mathbf{z}\in\mathbb{R}^M$ | $\mathbf{y}\in\mathbb{R}^K$

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
\mathbf{w} &= \mathbf{f}_{\mathrm{in}} * \mathbf{u} \\
\mathbf{z} &= \boldsymbol{\delta}_{\tau} * \mathbf{w} \\
\mathbf{y} &= \mathbf{f}_{\mathrm{out}} * \mathbf{z}.
\end{aligned}
```

## Initialization

None.

## Monitors

None.
