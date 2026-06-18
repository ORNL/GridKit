# Tau Model

`Tau` computes modal phase delay for a finite line length from modal
phase constants.

## Model Parameters

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\ell$ | [m] | Path length | from [`Path`](../../Geometry/Path/README.md)

### Parameter Validation

The connected `Gamma` model owns modal-vector validation. This model requires the
path length and current sweep point to satisfy

```math
\ell > 0,\qquad \omega > 0 .
```

### Model Derived Parameters

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\tau_m$ | [s] | Modal phase delay for mode $m$ | real

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$b_m$ | [1/m] | Modal phase constant for mode $m$ | from [`Gamma`](../Gamma/README.md)

## Model Equations

### Differential Equations

None.

### Algebraic Equations

For each mode $m=1,\dots,K$:

```math
0 = -\tau_m + \frac{\ell b_m}{\omega}
```

## Initialization

Initialize $\tau_m$ from $\ell$, the current $\omega$, and the current $b_m$
inputs.

## Model Outputs

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\boldsymbol{\tau}$ | [s] | Modal phase delay | $\boldsymbol{\tau}\in\mathbb{R}^K$
