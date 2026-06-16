# StateSpace Model

`StateSpace` represents a vector-fitted matrix rational approximation with
complex poles and factorized residues.

Notes:
- This cannot be used in general to replace `Rational`, as that model can support full-rank residuals, while this only supports rank-1 residuals.
- The benefit of this model is reduced number of internal states.

The rational approximation is represented in state-space form:

```math
\mathbf{H}(s) \approx \mathbf{D} + s\mathbf{E}
  + \mathbf{C}(s\mathbf{I} - \mathbf{P})^{-1}\mathbf{B}^T
```
The Laplace domain representation of this model is:
```math
\mathbf{Y}(s) = \mathbf{H}(s)\mathbf{U}(s)
```

The time domain representation of this model is:
```math
\mathbf{y}(t) = (\mathbf{h}*\mathbf{u})(t)
```

## Block Diagram

<div align="center">
   <img align="center" src="../../../../../../docs/Figures/EMT/StateSpace/diagram.png" width="90%">

  Figure 1: StateSpace rational approximation model
</div>

## Model Parameters

For output dimension $N$, input dimension $K$, and pole count $Q$:

Symbol | Units | JSON | Description | Typical Value | Note
------ | ----- | ---- | ----------- | ------------- | ----
$\mathbf{D}$ | [-] | `D` | Constant coefficient | | $\mathbf{D} \in \mathbb{R}^{N \times K}$
$\mathbf{E}$ | [s] | `E` | Linear coefficient | | $\mathbf{E} \in \mathbb{R}^{N \times K}$
$\mathbf{p}$ | [1/s] | `poles` | Poles | | $\mathbf{p} \in \mathbb{C}^Q$
$\mathbf{C}$ | [-] | `C` | Output matrix | | $\mathbf{C} \in \mathbb{C}^{N \times Q}$
$\mathbf{B}$ | [1/s] | `B` | Input matrix | | $\mathbf{B} \in \mathbb{C}^{K \times Q}$

### Parameter Validation

Complex-valued poles and factors must be ordered as adjacent conjugate pairs. For
each pair, with $q$ the first index:

```math
p_q = (p_{q+1})^*
```

The corresponding columns of $\mathbf{C}$ and $\mathbf{B}$ must follow the same
conjugate-pair ordering.

### Model Derived Parameters

Define the real-valued quantities used below. Entries indexed by $q$
correspond to poles.

```math
\begin{aligned}
\mathbf{P} &= \operatorname{diag}(p_1,\dots,p_Q) \\
\mathbf{a} &= \operatorname{Re}(\mathbf{p}) \\
\boldsymbol{\omega} &= \operatorname{Im}(\mathbf{p}) \\
\mathbf{A} &= \operatorname{diag}(a_1,\dots,a_Q) \\
\boldsymbol{\Omega} &= \operatorname{diag}(\omega_1,\dots,\omega_Q) \\
\mathbf{C}_{\mathrm{r}} &= \operatorname{Re}(\mathbf{C}) \\
\mathbf{C}_{\mathrm{i}} &= \operatorname{Im}(\mathbf{C}) \\
\mathbf{B}_{\mathrm{r}} &= \operatorname{Re}(\mathbf{B}) \\
\mathbf{B}_{\mathrm{i}} &= \operatorname{Im}(\mathbf{B})
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

In the general case, there are $2Q$ scalar internal differential states, grouped
as $\mathbf{x}_{\mathrm{r}}$ and $\mathbf{x}_{\mathrm{i}}$ vectors.
Real-valued poles do not need the imaginary state.

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{x}_{\mathrm{r}}$ | [-] | Real memory states | $\mathbf{x}_{\mathrm{r}} \in \mathbb{R}^Q$
$\mathbf{x}_{\mathrm{i}}$ | [-] | Imaginary memory states | $\mathbf{x}_{\mathrm{i}} \in \mathbb{R}^Q$

#### Algebraic

None.

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

For real-valued poles, the imaginary-memory equation is not needed.

```math
\begin{aligned}
0 &= -\dot{\mathbf{x}}_{\mathrm{r}}
     + \mathbf{A}\mathbf{x}_{\mathrm{r}}
     - \boldsymbol{\Omega}\mathbf{x}_{\mathrm{i}}
     + \mathbf{B}_{\mathrm{r}}^T\mathbf{u} \\
0 &= -\dot{\mathbf{x}}_{\mathrm{i}}
     + \boldsymbol{\Omega}\mathbf{x}_{\mathrm{r}}
     + \mathbf{A}\mathbf{x}_{\mathrm{i}}
     + \mathbf{B}_{\mathrm{i}}^T\mathbf{u}
\end{aligned}
```

### Algebraic Equations

None.

### Port Equations

```math
\mathbf{y} = \mathbf{D}\mathbf{u} + \mathbf{E}\dot{\mathbf{u}}
  + \mathbf{C}_{\mathrm{r}}\mathbf{x}_{\mathrm{r}}
  - \mathbf{C}_{\mathrm{i}}\mathbf{x}_{\mathrm{i}}
```

## Initialization

For an affine initial input trajectory, let subscript $0$ denote initial values.
The complex pole-memory state initializes to:

```math
\mathbf{x}_0 = -\mathbf{P}^{-1}\mathbf{B}^T\mathbf{u}_0 - \mathbf{P}^{-2}\mathbf{B}^T\dot{\mathbf{u}}_0
```

The real-valued state vectors and port contribution initialize to:

```math
\begin{aligned}
\mathbf{x}_{\mathrm{r},0} &= \operatorname{Re}(\mathbf{x}_0) \\
\mathbf{x}_{\mathrm{i},0} &= \operatorname{Im}(\mathbf{x}_0) \\
\mathbf{y}_0 &= \mathbf{D}\mathbf{u}_0 + \mathbf{E}\dot{\mathbf{u}}_0 + \operatorname{Re}(\mathbf{C}\mathbf{x}_0)
\end{aligned}
```

## Monitors

None.
