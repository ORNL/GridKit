# Gamma Model

`Gamma` computes the full-conductor propagation constant matrix and tracks its
modal propagation constants.

Notes:
- This model uses the V-definition: the series matrix pair is on the left of
  the shunt matrix pair.
- The real and imaginary parts of $\boldsymbol{\Gamma}$ are
  $\boldsymbol{\alpha}$ and $\boldsymbol{\beta}$.
- Eigenvalue computation is not expressed as a DAE, so modal branches are
  initialized from $\boldsymbol{\Gamma}^2$ at one value of $\omega$ and tracked
  as $\omega$ changes.
- The modal tracking equations assume each tracked modal branch is simple.

The propagation constant is decomposed as
$\boldsymbol{\Gamma}=\boldsymbol{\alpha}+j\boldsymbol{\beta}$.

The complex form typically given for this convention is:

```math
\boldsymbol{\Gamma} =
\sqrt{
(\mathbf{R}' + j\mathbf{X}')
(\mathbf{G}' + j\mathbf{B}')
}
```

Equivalently,

```math
\boldsymbol{\Gamma}^2 =
(\mathbf{R}' + j\mathbf{X}')
(\mathbf{G}' + j\mathbf{B}')
```

## Model Parameters

None. `Gamma` combines full-conductor parameter outputs for $K$ physical
conductors.

### Parameter Validation

The full-conductor inputs must be square and have the same dimension.

### Model Derived Parameters

None.

## Model Variables

### Internal Variables

#### Differential

A dot denotes derivative with respect to $\omega$.

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\boldsymbol{\alpha}$ | [1/m] | Attenuation constant matrix | $\mathbb{R}^{K\times K}$
$\boldsymbol{\beta}$ | [1/m] | Phase constant matrix | $\mathbb{R}^{K\times K}$
$a_m$ | [1/m] | Modal attenuation constant for mode $m$ | real
$b_m$ | [1/m] | Modal phase constant for mode $m$ | real
$\mathbf{t}^{\mathrm r}_{v,m}$ | [-] | Real part of voltage modal transformation vector | $\mathbf{t}^{\mathrm r}_{v,m}\in\mathbb{R}^K$
$\mathbf{t}^{\mathrm i}_{v,m}$ | [-] | Imaginary part of voltage modal transformation vector | $\mathbf{t}^{\mathrm i}_{v,m}\in\mathbb{R}^K$
$\mathbf{t}^{\mathrm r}_{i,m}$ | [-] | Real part of current modal transformation vector | $\mathbf{t}^{\mathrm r}_{i,m}\in\mathbb{R}^K$
$\mathbf{t}^{\mathrm i}_{i,m}$ | [-] | Imaginary part of current modal transformation vector | $\mathbf{t}^{\mathrm i}_{i,m}\in\mathbb{R}^K$

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{R}'$ | [$\Omega$/m] | Full-conductor resistance | from [`SeriesImpedance`](../../Effects/SeriesImpedance/README.md)
$\mathbf{X}'$ | [$\Omega$/m] | Full-conductor reactance | $\omega\mathbf{L}'$ from [`SeriesImpedance`](../../Effects/SeriesImpedance/README.md)
$\mathbf{G}'$ | [S/m] | Full-conductor shunt conductance | from [`ShuntAdmittance`](../../Effects/ShuntAdmittance/README.md)
$\mathbf{B}'$ | [S/m] | Full-conductor shunt susceptance | $\omega\mathbf{C}'$ from [`ShuntAdmittance`](../../Effects/ShuntAdmittance/README.md)

## Model Equations

### Differential Equations

For each mode $m=1,\dots,K$:

```math
\begin{aligned}
\mathbf{0} &=
  \dot{\boldsymbol{\alpha}}\mathbf{t}^{\mathrm r}_{v,m}
  - \dot{\boldsymbol{\beta}}\mathbf{t}^{\mathrm i}_{v,m}
  + \left(\boldsymbol{\alpha}-a_m\mathbf{I}\right)\dot{\mathbf{t}}^{\mathrm r}_{v,m}
  - \left(\boldsymbol{\beta}-b_m\mathbf{I}\right)\dot{\mathbf{t}}^{\mathrm i}_{v,m}
  - \dot{a}_m\mathbf{t}^{\mathrm r}_{v,m}
  + \dot{b}_m\mathbf{t}^{\mathrm i}_{v,m} \\
\mathbf{0} &=
  \dot{\boldsymbol{\alpha}}\mathbf{t}^{\mathrm i}_{v,m}
  + \dot{\boldsymbol{\beta}}\mathbf{t}^{\mathrm r}_{v,m}
  + \left(\boldsymbol{\beta}-b_m\mathbf{I}\right)\dot{\mathbf{t}}^{\mathrm r}_{v,m}
  + \left(\boldsymbol{\alpha}-a_m\mathbf{I}\right)\dot{\mathbf{t}}^{\mathrm i}_{v,m}
  - \dot{a}_m\mathbf{t}^{\mathrm i}_{v,m}
  - \dot{b}_m\mathbf{t}^{\mathrm r}_{v,m} \\
\mathbf{0} &=
  \dot{\boldsymbol{\alpha}}^{\mathsf{T}}\mathbf{t}^{\mathrm r}_{i,m}
  + \dot{\boldsymbol{\beta}}^{\mathsf{T}}\mathbf{t}^{\mathrm i}_{i,m}
  + \left(\boldsymbol{\alpha}^{\mathsf{T}}-a_m\mathbf{I}\right)\dot{\mathbf{t}}^{\mathrm r}_{i,m}
  + \left(\boldsymbol{\beta}^{\mathsf{T}}-b_m\mathbf{I}\right)\dot{\mathbf{t}}^{\mathrm i}_{i,m}
  - \dot{a}_m\mathbf{t}^{\mathrm r}_{i,m}
  - \dot{b}_m\mathbf{t}^{\mathrm i}_{i,m} \\
\mathbf{0} &=
  \dot{\boldsymbol{\alpha}}^{\mathsf{T}}\mathbf{t}^{\mathrm i}_{i,m}
  - \dot{\boldsymbol{\beta}}^{\mathsf{T}}\mathbf{t}^{\mathrm r}_{i,m}
  - \left(\boldsymbol{\beta}^{\mathsf{T}}-b_m\mathbf{I}\right)\dot{\mathbf{t}}^{\mathrm r}_{i,m}
  + \left(\boldsymbol{\alpha}^{\mathsf{T}}-a_m\mathbf{I}\right)\dot{\mathbf{t}}^{\mathrm i}_{i,m}
  - \dot{a}_m\mathbf{t}^{\mathrm i}_{i,m}
  + \dot{b}_m\mathbf{t}^{\mathrm r}_{i,m}
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
\mathbf{0} &= -\mathbf{R}'\mathbf{G}'
  + \mathbf{X}'\mathbf{B}'
  + \boldsymbol{\alpha}\boldsymbol{\alpha}
  - \boldsymbol{\beta}\boldsymbol{\beta} \\
\mathbf{0} &= -\mathbf{R}'\mathbf{B}'
  - \mathbf{X}'\mathbf{G}'
  + \boldsymbol{\alpha}\boldsymbol{\beta}
  + \boldsymbol{\beta}\boldsymbol{\alpha} \\
0 &= \left(\mathbf{t}^{\mathrm r}_{i,m}\right)^{\mathsf{T}}\mathbf{t}^{\mathrm r}_{v,m}
  + \left(\mathbf{t}^{\mathrm i}_{i,m}\right)^{\mathsf{T}}\mathbf{t}^{\mathrm i}_{v,m}
  - 1 \\
0 &= \left(\mathbf{t}^{\mathrm r}_{i,m}\right)^{\mathsf{T}}\mathbf{t}^{\mathrm i}_{v,m}
  - \left(\mathbf{t}^{\mathrm i}_{i,m}\right)^{\mathsf{T}}\mathbf{t}^{\mathrm r}_{v,m}
\end{aligned}
```

## Initialization

At the initial frequency, compute the eigenvalues and right eigenvectors of
$\boldsymbol{\Gamma}^2$:

```math
\boldsymbol{\Gamma}^2\mathbf{t}_m
  = \mu_m\mathbf{t}_m
```

Then initialize
$\lambda_m=\sqrt{\mu_m}=a_m+jb_m$ and construct
$\boldsymbol{\Gamma}$ with the same eigenvectors:

```math
\begin{aligned}
\boldsymbol{\Gamma}\mathbf{t}_m
  &= \lambda_m\mathbf{t}_m,\qquad
\mathbf{l}_m^*\mathbf{t}_m=1
\end{aligned}
```

Use $\lambda_m=a_m+jb_m$,
$\mathbf{t}_m=\mathbf{t}^{\mathrm r}_{v,m}+j\mathbf{t}^{\mathrm i}_{v,m}$, and
$\mathbf{l}_m=\mathbf{t}^{\mathrm r}_{i,m}+j\mathbf{t}^{\mathrm i}_{i,m}$.

## Monitors

Monitor | Symbol | Units | Shape | Description
------- | ------ | ----- | ----- | -----------
`Tv` | $\mathbf{T}_v$ | [-] | $K\times K$ | Voltage modal transformation matrix
`Ti` | $\mathbf{T}_i$ | [-] | $K\times K$ | Current modal transformation matrix
`Alpha` | $\mathbf{a}$ | [1/m] | $K$ | Modal attenuation constants
`Beta` | $\mathbf{b}$ | [1/m] | $K$ | Modal phase constants

The monitored transformation matrices use

```math
\mathbf{T}_v =
\mathbf{T}^{\mathrm r}_v+j\mathbf{T}^{\mathrm i}_v,\qquad
\mathbf{T}_i =
\mathbf{T}^{\mathrm r}_i+j\mathbf{T}^{\mathrm i}_i,\qquad
\mathbf{T}_i^{\mathsf H}\mathbf{T}_v=\mathbf{I}.
```

Both $\mathbf{T}_v$ and $\mathbf{T}_i$ are reported as $K\times K$
conductor-by-mode matrices. Column $m$ matches the ordering of `Alpha_m` and
`Beta_m`.

The voltage and current modal transforms are dual:

```math
\mathbf{v}=\mathbf{T}_v\hat{\mathbf{v}},\qquad
\hat{\mathbf{v}}=\mathbf{T}_i^{\mathsf H}\mathbf{v}
```

```math
\mathbf{i}=\mathbf{T}_i\hat{\mathbf{i}},\qquad
\hat{\mathbf{i}}=\mathbf{T}_v^{\mathsf H}\mathbf{i}.
```

GridKit uses the current-form propagation matrix because `LineDistributed`
propagates reflected current to incident current:

```math
\mathbf{H}^{\mathrm{mps}}_i(s)
=
\mathbf{T}_i(s)\widehat{\mathbf{H}}^{\mathrm{mps}}(s)
\mathbf{T}_v(s)^{\mathsf H}.
```

The dual voltage-form propagation matrix is

```math
\mathbf{H}^{\mathrm{mps}}_v(s)
=
\mathbf{T}_v(s)\widehat{\mathbf{H}}^{\mathrm{mps}}(s)
\mathbf{T}_i(s)^{\mathsf H}.
```

The phase-domain matrices $\boldsymbol{\alpha}$ and $\boldsymbol{\beta}$ are
internal `Gamma` states and are not monitors.
