# StateSpace Model

`StateSpace` represents a vector-fitted matrix rational approximation with
real or complex poles and factorized residue terms:

```math
\mathbf{H}(s) \approx \mathbf{D} + s\mathbf{E}
  + \mathbf{C}(s\mathbf{I} - \mathbf{P})^{-1}\mathbf{B}
```

> [!NOTE]
> A template parameter will select whether the input $\mathbf{u}$ is a
> differential or algebraic vector. This page documents differential input;
> algebraic input requires $\mathbf{E} = \mathbf{0}$.

> [!NOTE]
> Each factorized term has residue
> $\mathbf{R}_q = \mathbf{C}_{:,q}\mathbf{B}_{q,:}$, so
> $\mathrm{rank}\left(\mathbf{R}_q\right) \le 1$. Terms with the same
> pole may be repeated; their residues sum and can have higher rank.

## Block Diagram

![StateSpace rational-operator block diagram](../../../../../../docs/Figures/EMT/StateSpace/diagram.png)

Figure 1: StateSpace model

## Model Parameters

For output dimension $N$, input dimension $K$, pole count $Q$, input units
$[u]$, and output units $[y]$:

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\mathbf{D}$ | $[y]/[u]$ | `D` | Constant coefficient | $\mathbf{D} \in \mathbb{R}^{N \times K}$
$\mathbf{E}$ | $\mathrm{s}[y]/[u]$ | `E` | Linear coefficient | $\mathbf{E} \in \mathbb{R}^{N \times K}$
$\mathbf{p}$ | [1/s] | `poles` | Poles | $\mathbf{p} \in \mathbb{C}^Q$
$\mathbf{C}$ | $[y]/[u]$ | `C` | Output matrix | $\mathbf{C} \in \mathbb{C}^{N \times Q}$
$\mathbf{B}$ | [1/s] | `B` | Input matrix | $\mathbf{B} \in \mathbb{C}^{Q \times K}$

### Parameter Validation

The input and output dimensions are positive, and the pole count is
nonnegative. Poles are nonzero, real poles have real factors, and nonreal poles
and factors are stored as adjacent conjugate pairs, with $q$ the first index of
each pair.

```math
\begin{aligned}
N &> 0 \\
K &> 0 \\
Q &\ge 0 \\
p_q &\ne 0,
  \quad q \in \{1,\ldots,Q\} \\
p_q \in \mathbb{R}
  &\Longrightarrow
  \mathbf{C}_{:,q} \in \mathbb{R}^{N},
  \quad \mathbf{B}_{q,:} \in \mathbb{R}^{K} \\
p_{q+1} &= p_q^{\ast} \\
\mathbf{C}_{:,q+1} &= \mathbf{C}_{:,q}^{\ast} \\
\mathbf{B}_{q+1,:} &= \mathbf{B}_{q,:}^{\ast}
\end{aligned}
```

### Derived Parameters

```math
\begin{aligned}
\mathbf{P} &= \mathrm{diag}(p_1,\ldots,p_Q) \\
\mathbf{a} &= \mathrm{Re}(\mathbf{p}) \\
\boldsymbol{\omega} &= \mathrm{Im}(\mathbf{p}) \\
\mathbf{A} &= \mathrm{diag}(a_1,\ldots,a_Q) \\
\boldsymbol{\Omega} &= \mathrm{diag}(\omega_1,\ldots,\omega_Q) \\
\mathbf{C}_{\mathrm{r}} &= \mathrm{Re}(\mathbf{C}) \\
\mathbf{C}_{\mathrm{i}} &= \mathrm{Im}(\mathbf{C}) \\
\mathbf{B}_{\mathrm{r}} &= \mathrm{Re}(\mathbf{B}) \\
\mathbf{B}_{\mathrm{i}} &= \mathrm{Im}(\mathbf{B})
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
$\mathbf{w}$ | $[u]$ | Real memory states | $\mathbf{w} \in \mathbb{R}^Q$
$\mathbf{v}$ | $[u]$ | Imaginary memory states | $\mathbf{v} \in \mathbb{R}^Q$

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
0 &= -\dfrac{\mathrm{d}\mathbf{w}}{\mathrm{d}t}
     + \mathbf{A}\mathbf{w}
     - \boldsymbol{\Omega}\mathbf{v}
     + \mathbf{B}_{\mathrm{r}}\mathbf{u} \\
0 &= -\dfrac{\mathrm{d}\mathbf{v}}{\mathrm{d}t}
     + \boldsymbol{\Omega}\mathbf{w}
     + \mathbf{A}\mathbf{v}
     + \mathbf{B}_{\mathrm{i}}\mathbf{u}
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
0 &= -\mathbf{y}
     + \mathbf{D}\mathbf{u}
     + \mathbf{E}\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
     + \mathbf{C}_{\mathrm{r}}\mathbf{w}
     - \mathbf{C}_{\mathrm{i}}\mathbf{v}
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

When poles are present, let $\mathbf{x}$ be the complex memory vector:

```math
\begin{aligned}
\mathbf{x}
  &\leftarrow
  -\mathbf{P}^{-1}\mathbf{B}\mathbf{u}
  - \mathbf{P}^{-2}\mathbf{B}
    \dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \\
\mathbf{w}
  &\leftarrow \mathrm{Re}(\mathbf{x}) \\
\mathbf{v}
  &\leftarrow \mathrm{Im}(\mathbf{x}) \\
\dfrac{\mathrm{d}\mathbf{w}}{\mathrm{d}t}
  &\leftarrow
  \mathbf{A}\mathbf{w}
  - \boldsymbol{\Omega}\mathbf{v}
  + \mathbf{B}_{\mathrm{r}}\mathbf{u} \\
\dfrac{\mathrm{d}\mathbf{v}}{\mathrm{d}t}
  &\leftarrow
  \boldsymbol{\Omega}\mathbf{w}
  + \mathbf{A}\mathbf{v}
  + \mathbf{B}_{\mathrm{i}}\mathbf{u}
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
\mathbf{y}
  &\leftarrow
  \mathbf{D}\mathbf{u}
  + \mathbf{E}\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
  + \mathbf{C}_{\mathrm{r}}\mathbf{w}
  - \mathbf{C}_{\mathrm{i}}\mathbf{v} \\
\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}
  &\leftarrow
  \mathbf{D}\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
  + \mathbf{C}_{\mathrm{r}}
    \dfrac{\mathrm{d}\mathbf{w}}{\mathrm{d}t}
  - \mathbf{C}_{\mathrm{i}}
    \dfrac{\mathrm{d}\mathbf{v}}{\mathrm{d}t}
\end{aligned}
```

## Monitors

None.
