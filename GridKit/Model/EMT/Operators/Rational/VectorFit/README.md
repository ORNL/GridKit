# VectorFit Model

`VectorFit` represents a vector-fitted matrix rational approximation with
complex poles and residues.

The rational approximation is represented in pole form:

```math
\mathbf{H}(s) \approx \mathbf{D} + s\mathbf{E}
  + \sum_{q=1}^{Q} \frac{\mathbf{R}_q}{s - p_q}
```

The Laplace domain representation of this model is:
```math
\mathbf{Y}(s) = \mathbf{H}(s)\mathbf{U}(s)
```

The time domain operator form is:
```math
\mathbf{y}(t) = f^{\mathbf{h}}(\mathbf{u})(t)
```

## Block Diagram

![](../../../../../../docs/Figures/EMT/VecFit/diagram.png)

Figure 1: VectorFit model

## Model Parameters

For output dimension $N$, input dimension $K$, and pole count $Q$:

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\mathbf{D}$ | [-] | `D` | Constant coefficient | $\mathbf{D}\in\mathbb{R}^{N\times K}$
$\mathbf{E}$ | [s] | `E` | Linear coefficient | $\mathbf{E}\in\mathbb{R}^{N\times K}$
$\mathbf{p}$ | [1/s] | `poles` | Poles | $\mathbf{p}\in\mathbb{C}^Q$
$\mathbf{R}$ | [1/s] | `residues` | Residues | $\mathbf{R}\in\mathbb{C}^{N\times K\times Q}$

### Parameter Validation

Complex-valued poles and residues must be ordered as adjacent conjugate pairs.
For each pair, with $q$ the first index:

```math
\begin{aligned}
p_q &= (p_{q+1})^{\ast} &
\mathbf{R}_q &= (\mathbf{R}_{q+1})^{\ast}
\end{aligned}
```

### Derived Parameters

```math
\begin{aligned}
\mathbf{a} &= \mathrm{Re}(\mathbf{p}) \\
\boldsymbol{\omega} &= \mathrm{Im}(\mathbf{p}) \\
\mathbf{A} &= \mathrm{Re}(\mathbf{R}) \\
\mathbf{B} &= \mathrm{Im}(\mathbf{R})
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
$\mathbf{w}_q$ | [-] | Real memory states | $\mathbf{w}_q\in\mathbb{R}^K$
$\mathbf{v}_q$ | [-] | Imaginary memory states | $\mathbf{v}_q\in\mathbb{R}^K$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{y}$ | [-] | Output contribution | $\mathbf{y}\in\mathbb{R}^N$

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | [-] | Input vector | $\mathbf{u} \in \mathbb{R}^K$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | [-] | Input vector port | $\mathbf{u} \in \mathbb{R}^K$
$\mathbf{y}$ | `out` | Output | [-] | Output contribution port | $\mathbf{y} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= -\dot{\mathbf{w}}_q
     + a_q\mathbf{w}_q
     - \omega_q\mathbf{v}_q
     + \mathbf{u} \\
0 &= -\dot{\mathbf{v}}_q
     + \omega_q\mathbf{w}_q
     + a_q\mathbf{v}_q
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
0 &= -\mathbf{y}
     + \mathbf{D}\mathbf{u}
     + \mathbf{E}\dot{\mathbf{u}}
     + \sum_{q=1}^{Q}
       \left(
         \mathbf{A}_q\mathbf{w}_q
         - \mathbf{B}_q\mathbf{v}_q
       \right)
\end{aligned}
```

### Wiring

None.

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{u},\dot{\mathbf{u}}
  &\leftarrow \text{affine input trajectory start}
\end{aligned}
```

### Internal Initialization

Initialization is performed by evaluating the affine-input residuals in
dependency order:

```math
\begin{aligned}
\mathbf{w}_{q}
  &\leftarrow
  -\frac{a_q}{a_q^2 + \omega_q^2}\mathbf{u}
  - \frac{a_q^2 - \omega_q^2}{(a_q^2 + \omega_q^2)^2}\dot{\mathbf{u}} \\
\mathbf{v}_{q}
  &\leftarrow
  \frac{\omega_q}{a_q^2 + \omega_q^2}\mathbf{u}
  + \frac{2a_q\omega_q}{(a_q^2 + \omega_q^2)^2}\dot{\mathbf{u}}
\end{aligned}
```

### Output Initialization

```math
\mathbf{y}
  \leftarrow
  \mathbf{D}\mathbf{u}
  + \mathbf{E}\dot{\mathbf{u}}
  + \sum_{q\in\{1,\ldots,Q\}}
    \left(
      \mathbf{A}_q\mathbf{w}_{q}
      - \mathbf{B}_q\mathbf{v}_{q}
    \right)
```

## Monitors

None.
