# Propagation Model

For input units $[u]$, `Propagation` represents the current-form EMT propagation
operator used by `LineDistributed`. It composes fitted input and output factors
with one scalar delay per mode and preserves the input signal units.

```math
\mathbf{y} =
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
\mathbf{T}_i(s)
\end{aligned}
```

$\widehat{\mathbf{H}}^{\mathrm{mps}}$ is the modal propagation function with
the modal delays removed. The operator $\boldsymbol{\delta}_{\tau}$ carries
those delays exclusively. The matrices $\mathbf{T}_v$ and $\mathbf{T}_i$ are
the voltage and current modal transformations.

## Block Diagram

![Propagation operator block diagram](../../../../../../docs/Figures/EMT/Propagation/diagram.png)

Figure 1: Propagation model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$K$ | [-] | `K` | Signal dimension | Required, positive integer
$\boldsymbol{\tau}$ | [s] | `tau` | Modal propagation delays | $\boldsymbol{\tau} \in \mathbb{R}^M$
$f_{\max}$ | [Hz] | `fmax` | Modal lag-chain section rate | Required, positive

### Parameter Validation

```math
\begin{aligned}
K &> 0 \\
M &> 0 \\
\boldsymbol{\tau} &\in \mathbb{R}^M \\
\tau_m &> 0,\quad m \in \{1,\ldots,M\} \\
f_{\max} &> 0
\end{aligned}
```

### Derived Parameters

```math
\begin{aligned}
M &= \mathrm{dim}\left(\boldsymbol{\tau}\right) \\
J_m &= \left\lceil f_{\max}\tau_m \right\rceil,
  \quad m \in \{1,\ldots,M\}
\end{aligned}
```

## Submodels

Submodel | Type | Order | Inputs | JSON | Outputs
-------- | ---- | ----- | ------ | ---- | -------
Input factor $f^{\mathbf{h}_{\mathrm{in}}}$ | [`VectorFit`](../../Rational/VectorFit/README.md) | $Q_{\mathbf{h}_{\mathrm{in}}}$ | $\mathbb{R}^K$ $[u]$ | `input` | $\mathbb{R}^M$ $[u]$
Modal delay $\delta_{\tau_m}$, $m \in \{1,\ldots,M\}$ | [`Delay`](../Delay/README.md) | $J_m$ | $\mathbb{R}$ $[u]$ | `tau[m]`, `fmax` | $\mathbb{R}$ $[u]$
Output factor $f^{\mathbf{h}_{\mathrm{out}}}$ | [`VectorFit`](../../Rational/VectorFit/README.md) | $Q_{\mathbf{h}_{\mathrm{out}}}$ | $\mathbb{R}^M$ $[u]$ | `output` | $\mathbb{R}^K$ $[u]$

### Submodel Validation

```math
\mathbf{E}^{\mathbf{h}_{\mathrm{in}}} = \mathbf{0}_{M \times K}
```

The zero linear coefficient makes the input factor compatible with the
algebraic input $\mathbf{u}$.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{z}$ | $[u]$ | Modal delayed signals | $\mathbf{z} \in \mathbb{R}^M$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{w}$ | $[u]$ | Input-factor outputs | $\mathbf{w} \in \mathbb{R}^M$
$\mathbf{y}$ | $[u]$ | Output vector | $\mathbf{y} \in \mathbb{R}^K$

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | $[u]$ | Input vector | $\mathbf{u} \in \mathbb{R}^K$

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | $[u]$ | Input vector port | $\mathbf{u} \in \mathbb{R}^K$
$\mathbf{y}$ | `out` | Output | $[u]$ | Output vector port | $\mathbf{y} \in \mathbb{R}^K$

## Model Equations

### Differential Equations

None.

### Algebraic Equations

None.

### Wiring

```math
\begin{aligned}
\mathbf{w} &= f^{\mathbf{h}_{\mathrm{in}}}(\mathbf{u}) \\
\mathbf{z} &= \boldsymbol{\delta}_{\tau}(\mathbf{w}) \\
\mathbf{y} &= f^{\mathbf{h}_{\mathrm{out}}}(\mathbf{z})
\end{aligned}
```

## Initialization

Initialization assumes an affine input with zero second derivative.

### Input Initialization

```math
\begin{aligned}
\mathbf{u},\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
  &\leftarrow \text{input value and derivative}
\end{aligned}
```

### Internal Initialization

```math
\begin{aligned}
\mathbf{w}
  &\leftarrow f^{\mathbf{h}_{\mathrm{in}}}(\mathbf{u}) \\
\dfrac{\mathrm{d}\mathbf{w}}{\mathrm{d}t}
  &\leftarrow
  \dfrac{\mathrm{d}}{\mathrm{d}t}
  \left[f^{\mathbf{h}_{\mathrm{in}}}(\mathbf{u})\right] \\
z_m
  &\leftarrow
  w_m - \tau_m\dfrac{\mathrm{d}w_m}{\mathrm{d}t},
     \quad m \in \{1,\ldots,M\} \\
\dfrac{\mathrm{d}z_m}{\mathrm{d}t}
  &\leftarrow \dfrac{\mathrm{d}w_m}{\mathrm{d}t},
     \quad m \in \{1,\ldots,M\} \\
\mathbf{z}
  &\leftarrow \left[z_1,\ldots,z_M\right]^{\mathsf T} \\
\dfrac{\mathrm{d}\mathbf{z}}{\mathrm{d}t}
  &\leftarrow
  \left[
    \dfrac{\mathrm{d}z_1}{\mathrm{d}t},
    \ldots,
    \dfrac{\mathrm{d}z_M}{\mathrm{d}t}
  \right]^{\mathsf T}
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
\mathbf{y}
  &\leftarrow f^{\mathbf{h}_{\mathrm{out}}}(\mathbf{z}) \\
\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}
  &\leftarrow
  \dfrac{\mathrm{d}}{\mathrm{d}t}
  \left[f^{\mathbf{h}_{\mathrm{out}}}(\mathbf{z})\right]
\end{aligned}
```

## Monitors

None.
