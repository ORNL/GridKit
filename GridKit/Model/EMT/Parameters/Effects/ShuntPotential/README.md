# ShuntPotential Model

`ShuntPotential` computes the real overhead-line potential-coefficient matrix
and the corresponding full-conductor shunt capacitance. The potential-derived
conductance is set to zero; direct shunt losses are added by separate passive
shunt models.

TODO: high frequnecy earth return correction

## Model Parameters

For $K$ physical conductors:

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
`tower` | [-] | Static conductor tower | [`Tower`](../../Geometry/Tower/README.md)
`conductor` | [-] | Static conductor data | [`Conductor`](../../Geometry/Conductor/README.md)
$\varepsilon_0$ | [F/m] | Vacuum permittivity | constant

### Parameter Validation

The tower and conductor models must pass their own validation.

### Model Derived Parameters

The tower model provides $h_i$, $D_{ij}$, and $D'_{ij}$. The conductor model
provides $r_i$.

```math
\Lambda^{\mathrm{pot}}_{ij} =
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
$\mathbf{P}$ | [m/F] | Potential-coefficient matrix | $\mathbb{R}^{K\times K}$
$\mathbf{G}^{\mathrm{pot}}$ | [S/m] | Potential-derived shunt conductance | zero matrix
$\mathbf{C}^{\mathrm{pot}}$ | [F/m] | Potential-derived shunt capacitance | $\mathbb{R}^{K\times K}$

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
\begin{aligned}
\mathbf{0} &= -2\pi\varepsilon_0\mathbf{P}
  + \boldsymbol{\Lambda}^{\mathrm{pot}} \\
\mathbf{0} &= -\mathbf{G}^{\mathrm{pot}} \\
\mathbf{0} &= \mathbf{P}\mathbf{C}^{\mathrm{pot}}
  - \mathbf{I}_{K}
\end{aligned}
```

## Initialization

Initialize $\mathbf{P}$ from the log-distance terms, set
$\mathbf{G}^{\mathrm{pot}}=\mathbf{0}$, and solve
$\mathbf{P}\mathbf{C}^{\mathrm{pot}}=\mathbf{I}_{K}$.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{G}^{\mathrm{pot}}$ | [S/m] | Potential-derived shunt conductance | zero matrix
$\mathbf{C}^{\mathrm{pot}}$ | [F/m] | Potential-derived shunt capacitance | $\mathbb{R}^{K\times K}$
