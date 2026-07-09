# VectorFit Model

`VectorFit` represents a vector-fitted matrix rational approximation with real
or complex poles and general residue matrices:

```math
\mathbf{H}(s) \approx \mathbf{D} + s\mathbf{E}
  + \sum_{q \in \{1,\ldots,Q\}} \dfrac{\mathbf{R}_q}{s - p_q}
```

> [!NOTE]
> A template parameter will select whether the input $\mathbf{u}$ is a
> differential or algebraic vector. This page documents differential input;
> algebraic input requires $\mathbf{E} = \mathbf{0}$.

## Block Diagram

![VectorFit rational-operator block diagram](../../../../../../docs/Figures/EMT/VecFit/diagram.png)

Figure 1: VectorFit model

## Model Parameters

For output dimension $N$, input dimension $K$, pole count $Q$, input units
$[u]$, and output units $[y]$:

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\mathbf{D}$ | $[y]/[u]$ | `D` | Constant coefficient | $\mathbf{D} \in \mathbb{R}^{N \times K}$
$\mathbf{E}$ | $\mathrm{s}[y]/[u]$ | `E` | Linear coefficient | $\mathbf{E} \in \mathbb{R}^{N \times K}$
$\mathbf{p}$ | [1/s] | `poles` | Poles | $\mathbf{p} \in \mathbb{C}^Q$
$\mathbf{R}$ | $[y]/(\mathrm{s}[u])$ | `residues` | Residues | $\mathbf{R} \in \mathbb{C}^{N \times K \times Q}$

### Parameter Validation

The input and output dimensions are positive, and the pole count is
nonnegative. Poles are nonzero, real poles have real residues, and nonreal poles
and residues are stored as adjacent conjugate pairs, with $q$ the first index of
each pair.

```math
\begin{aligned}
N &> 0 \\
K &> 0 \\
Q &\ge 0 \\
p_q &\ne 0,
  \quad q \in \{1,\ldots,Q\} \\
p_q \in \mathbb{R}
  &\Longrightarrow \mathbf{R}_q \in \mathbb{R}^{N \times K} \\
p_{q+1} &= p_q^{\ast} \\
\mathbf{R}_{q+1} &= \mathbf{R}_q^{\ast}
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
$\mathbf{w}_q$ | $\mathrm{s}[u]$ | Real memory states | $\mathbf{w}_q \in \mathbb{R}^K$
$\mathbf{v}_q$ | $\mathrm{s}[u]$ | Imaginary memory states | $\mathbf{v}_q \in \mathbb{R}^K$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{y}$ | $[y]$ | Output vector | $\mathbf{y} \in \mathbb{R}^N$

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
$\mathbf{y}$ | `out` | Output | $[y]$ | Output vector port | $\mathbf{y} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= -\dfrac{\mathrm{d}\mathbf{w}_q}{\mathrm{d}t}
     + a_q\mathbf{w}_q
     - \omega_q\mathbf{v}_q
     + \mathbf{u} \\
0 &= -\dfrac{\mathrm{d}\mathbf{v}_q}{\mathrm{d}t}
     + \omega_q\mathbf{w}_q
     + a_q\mathbf{v}_q,
     \quad q \in \{1,\ldots,Q\}
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
0 &= -\mathbf{y}
     + \mathbf{D}\mathbf{u}
     + \mathbf{E}\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
     + \sum_{q \in \{1,\ldots,Q\}}
       \left(
         \mathbf{A}_q\mathbf{w}_q
         - \mathbf{B}_q\mathbf{v}_q
       \right)
\end{aligned}
```

### Wiring

None.

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
\mathbf{w}_{q}
  &\leftarrow
  -\dfrac{a_q}{a_q^2 + \omega_q^2}\mathbf{u}
  - \dfrac{a_q^2 - \omega_q^2}{(a_q^2 + \omega_q^2)^2}
    \dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \\
\mathbf{v}_{q}
  &\leftarrow
  \dfrac{\omega_q}{a_q^2 + \omega_q^2}\mathbf{u}
  + \dfrac{2a_q\omega_q}{(a_q^2 + \omega_q^2)^2}
    \dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t},
    \quad q \in \{1,\ldots,Q\} \\
\dfrac{\mathrm{d}\mathbf{w}_q}{\mathrm{d}t}
  &\leftarrow
  a_q\mathbf{w}_q - \omega_q\mathbf{v}_q + \mathbf{u} \\
\dfrac{\mathrm{d}\mathbf{v}_q}{\mathrm{d}t}
  &\leftarrow
  \omega_q\mathbf{w}_q + a_q\mathbf{v}_q,
  \quad q \in \{1,\ldots,Q\}
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
\mathbf{y}
  &\leftarrow
  \mathbf{D}\mathbf{u}
  + \mathbf{E}\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
  + \sum_{q \in \{1,\ldots,Q\}}
    \left(
      \mathbf{A}_q\mathbf{w}_{q}
      - \mathbf{B}_q\mathbf{v}_{q}
    \right) \\
\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}
  &\leftarrow
  \mathbf{D}\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
  + \sum_{q \in \{1,\ldots,Q\}}
    \left(
      \mathbf{A}_q\dfrac{\mathrm{d}\mathbf{w}_{q}}{\mathrm{d}t}
      - \mathbf{B}_q\dfrac{\mathrm{d}\mathbf{v}_{q}}{\mathrm{d}t}
    \right)
\end{aligned}
```

## Monitors

None.
