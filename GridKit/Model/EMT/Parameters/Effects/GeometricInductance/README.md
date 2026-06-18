# GeometricInductance Model

`GeometricInductance` computes the external geometric series inductance matrix
for the physical conductors.

## Model Parameters

For $K$ physical conductors:

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
`tower` | [-] | Static conductor tower | [`Tower`](../../Geometry/Tower/README.md)
`conductor` | [-] | Static conductor data | [`Conductor`](../../Geometry/Conductor/README.md)
$\mu_0$ | [H/m] | Vacuum permeability | constant

### Parameter Validation

The tower and conductor models must pass their own validation.

### Model Derived Parameters

The tower model provides $h_i$, $D_{ij}$, and $D'_{ij}$. The conductor model
provides $r_i$.

```math
\left(\boldsymbol{\Lambda}^{\mathrm{geo}}\right)_{ij} =
\begin{cases}
\ln\dfrac{2h_i}{r_i}, & i=j \\
\ln\dfrac{D'_{ij}}{D_{ij}}, & i\ne j
\end{cases}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{L}^{\mathrm{geo}}$ | [H/m] | External geometric inductance | $\mathbb{R}^{K\times K}$

### External Variables

#### Differential

None.

#### Algebraic

None.

## Model Equations

### Differential Equations

None.

### Algebraic Equations

```math
\mathbf{0} = -2\pi\mathbf{L}^{\mathrm{geo}}
  + \mu_0 \boldsymbol{\Lambda}^{\mathrm{geo}}
```

## Initialization

None beyond the static tower and conductor initialization.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{L}^{\mathrm{geo}}$ | [H/m] | External geometric inductance | $\mathbb{R}^{K\times K}$
