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
$K$ | [-] | `K` | Number of conductors | Required, positive integer
$\boldsymbol{\tau}$ | [s] | `tau` | Modal propagation delays | $\boldsymbol{\tau} \in \mathbb{R}^M$
$f_{\max}$ | [Hz] | `fmax` | Modal lag-chain section rate | Required, positive. Not a bandwidth guarantee

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
J_m &= \mathrm{ceil}(f_{\max}\tau_m),\quad m \in \{1,\ldots,M\} \\
Q_{\mathbf{h}} &=
  Q_{\mathbf{h}_{\mathrm{in}}}
  + \sum_{m = 1}^{M}J_m
  + Q_{\mathbf{h}_{\mathrm{out}}}
\end{aligned}
```

$Q_{\mathbf{h}}$ is the aggregate submodel order, with one scalar `Delay`
instance for each mode.

## Submodels

Submodel | Type | Order | Inputs | JSON | Outputs
-------- | ---- | ----- | ------ | ---- | -------
Input factor $f^{\mathbf{h}_{\mathrm{in}}}$ | [`VectorFit`](../../Rational/VectorFit/README.md) | $Q_{\mathbf{h}_{\mathrm{in}}}$ | $\mathbf{u} \in \mathbb{R}^K$ | `input` | $\mathbf{w} \in \mathbb{R}^M$
Modal delay $\delta_{\tau_m}$, $m \in \{1,\ldots,M\}$ | [`Delay`](../Delay/README.md) | $J_m$ | $w_m \in \mathbb{R}$ | `tau[m]`, `fmax` | $z_m \in \mathbb{R}$
Output factor $f^{\mathbf{h}_{\mathrm{out}}}$ | [`VectorFit`](../../Rational/VectorFit/README.md) | $Q_{\mathbf{h}_{\mathrm{out}}}$ | $\mathbf{z} \in \mathbb{R}^M$ | `output` | $\mathbf{y} \in \mathbb{R}^K$

### Submodel Validation

```math
\begin{aligned}
f^{\mathbf{h}_{\mathrm{in}}} &: \mathbb{R}^K \rightarrow \mathbb{R}^M \\
\delta_{\tau_m} &: \mathbb{R} \rightarrow \mathbb{R},\quad m \in \{1,\ldots,M\} \\
f^{\mathbf{h}_{\mathrm{out}}} &: \mathbb{R}^M \rightarrow \mathbb{R}^K
\end{aligned}
```

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
$\mathbf{y}$ | $[u]$ | Output-factor contribution | $\mathbf{y} \in \mathbb{R}^K$

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | $[u]$ | Input vector | $\mathbf{u} \in \mathbb{R}^K$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | $[u]$ | Input vector port | $\mathbf{u} \in \mathbb{R}^K$
$\mathbf{y}$ | `out` | Output | $[u]$ | Output contribution port | $\mathbf{y} \in \mathbb{R}^K$

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
     \quad m \in \{1,\ldots,M\} \\
\mathbf{z} &= \boldsymbol{\delta}_{\tau}(\mathbf{w}) \\
\mathbf{y} &= f^{\mathbf{h}_{\mathrm{out}}}(\mathbf{z})
\end{aligned}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{u},\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
  &\leftarrow \text{affine input trajectory start}
\end{aligned}
```

### Internal Initialization

Each submodel is initialized from the value and derivative of its input
trajectory. The input factor is initialized first, followed by the $M$ scalar
delays:

```math
\begin{aligned}
\mathbf{w},\dfrac{\mathrm{d}\mathbf{w}}{\mathrm{d}t}
  &\leftarrow \text{initialized input-factor output trajectory} \\
z_m,\dfrac{\mathrm{d}z_m}{\mathrm{d}t}
  &\leftarrow \text{initialized modal-delay output trajectory},
     \quad m \in \{1,\ldots,M\} \\
\mathbf{z}
  &\leftarrow \left[z_1,\ldots,z_M\right]^{\mathsf T}
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
\mathbf{y}
  &\leftarrow f^{\mathbf{h}_{\mathrm{out}}}(\mathbf{z}) \\
\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}
  &\leftarrow \text{initialized output-factor derivative}
\end{aligned}
```

The output factor uses $\mathbf{z}$ and
$\mathrm{d}\mathbf{z}/\mathrm{d}t$ as its initialized input trajectory.

## Monitors

None.
