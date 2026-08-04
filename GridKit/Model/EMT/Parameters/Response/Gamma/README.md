# Gamma Model

`Gamma` computes the full-conductor propagation constant matrix.

Notes:
- This model uses the V-definition: the series matrix pair is on the left of
  the shunt matrix pair.
- The real and imaginary parts of $\boldsymbol{\Gamma}$ are
  $\boldsymbol{\alpha}$ and $\boldsymbol{\beta}$.
- $\boldsymbol{\Gamma}$ is analytic in $\omega$ for every conductor
  geometry, degenerate modal spectra included, so it is the only propagation
  quantity carried by the sweep DAE. Its spectral factorization is not
  analytic and is observed outside the DAE by
  [`ModalDecomposition`](../../Modal/README.md).

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

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\boldsymbol{\alpha}$ | [1/m] | Attenuation constant matrix | $\mathbb{R}^{K\times K}$
$\boldsymbol{\beta}$ | [1/m] | Phase constant matrix | $\mathbb{R}^{K\times K}$

### External Variables

None.

## Model Equations

The states are pinned by the matrix square relation on the branch fixed at
initialization, split into real and imaginary parts:

```math
\boldsymbol{\alpha}^2 - \boldsymbol{\beta}^2
  - \mathbf{R}'\mathbf{G}' + \mathbf{X}'\mathbf{B}' = \mathbf{0},
\qquad
\boldsymbol{\alpha}\boldsymbol{\beta} + \boldsymbol{\beta}\boldsymbol{\alpha}
  - \mathbf{R}'\mathbf{B}' - \mathbf{X}'\mathbf{G}' = \mathbf{0}.
```

The Frechet derivative of the square, $\mathbf{E} \mapsto
\boldsymbol{\Gamma}\mathbf{E} + \mathbf{E}\boldsymbol{\Gamma}$, is
nonsingular whenever no two eigenvalues of $\boldsymbol{\Gamma}$ sum to
zero, so repeated modal eigenvalues never make these equations singular.

## Initialization

The states are seeded with the principal matrix square root of
$(\mathbf{R}' + j\mathbf{X}')(\mathbf{G}' + j\mathbf{B}')$ at the starting
frequency; the residual then holds the states on this branch across the
sweep.
