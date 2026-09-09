# Norton Model

`Norton` realizes the shunt branch of a Norton terminal. The bus registers
the incident current and subtracts this branch's output in KCL.
The branch owns the algebraic shunt current $\mathbf{i}^\mathrm{sh}$ and admittance states.

## Block Diagram

None.

## Model Parameters

None.

### Parameter Validation

None.

### Derived Parameters

The terminal dimension is fixed at $K=3$.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `v` | Input | [V] | Terminal voltage | $\mathbb{R}^K$
$\mathbf{i}^\mathrm{sh}$ | `Ish` | Output | [A] | Shunt current | $\mathbb{R}^K$

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{y}$ | Shunt admittance | [VectorFit](../../../Operators/Rational/VectorFit/README.md) | $KQ$ | Device coefficients | $\mathbf{v}$ | $\mathbb{R}^K$

### Submodel Validation

The admittance has $K$ inputs and outputs and satisfies its coefficient constraints.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}^\mathrm{sh}$ | [A] | Shunt current | $\mathbb{R}^K$

### External Variables

#### Differential

When a connected equation depends on the bus-voltage derivative:

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Bus voltage | Algebraic when no equation depends on its derivative

#### Algebraic

Otherwise, the bus-voltage variables above are algebraic.

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

```math
0 = -\mathbf{i}^\mathrm{sh}+\mathbf{y}[\mathbf{v}]
```

### External Equations

None.

## Initialization

Initialize the admittance states, then set $\mathbf{i}^\mathrm{sh}$ from its output.

## Monitors

None.
