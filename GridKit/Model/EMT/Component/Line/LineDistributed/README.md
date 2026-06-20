# LineDistributed Model

`LineDistributed` represents a distributed EMT line with characteristic
admittance $\mathbf{y}_c$ and propagation model $\mathbf{h}$.

## Block Diagram

<div align="center">
   <img align="center" src="../../../../../../docs/Figures/EMT/LineDistributed/diagram.png">

  Figure 1: LineDistributed model
</div>

## Model Parameters

Symbol            | Units | JSON         | Description                                            | Note
----------------- | ----- | ------------ | ------------------------------------------------------ | ----
$\mathbf{P}_\phi$ | [-]   | `conductors` | Permutation matrix mapping each conductor to its phase | $\mathbf{P}_\phi \in \mathbb{R}^{N \times K}$

### Parameter Validation

None.

### Model Derived Parameters

The number of phases $N$ and conductors $K$ are the row and column counts of
$\mathbf{P}_\phi$, respectively.

### Submodels

Submodel | Inputs | Parameters | Outputs
-------- | ------ | ---------- | -------
[`VectorFit`](../../../Operators/Rational/VectorFit/README.md) characteristic admittance $\mathbf{y}_c$ | $\mathbf{v}_{1}\in\mathbb{R}^N$ | `Yc` | $\mathbf{i}^{\mathrm{sh}}_{1}\in\mathbb{R}^K$
[`VectorFit`](../../../Operators/Rational/VectorFit/README.md) characteristic admittance $\mathbf{y}_c$ | $\mathbf{v}_{2}\in\mathbb{R}^N$ | `Yc` | $\mathbf{i}^{\mathrm{sh}}_{2}\in\mathbb{R}^K$
[`Propagation`](../../../Operators/Shift/Propagation/README.md) $\mathbf{h}$ | $2\mathbf{i}^{\mathrm{sh}}_{2}-\mathbf{i}^{\mathrm{inc}}_{2}$ | `H` | $\mathbf{i}^{\mathrm{inc}}_{1}\in\mathbb{R}^K$
[`Propagation`](../../../Operators/Shift/Propagation/README.md) $\mathbf{h}$ | $2\mathbf{i}^{\mathrm{sh}}_{1}-\mathbf{i}^{\mathrm{inc}}_{1}$ | `H` | $\mathbf{i}^{\mathrm{inc}}_{2}\in\mathbb{R}^K$

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

None.

### External Variables

External variables are owned by the EMT bus.

#### Differential

Symbol                               | Units | Description                              | Note
-------------------------------------|-------|------------------------------------------|-----
$\mathbf{v}_{1}$         | [V]   | Terminal `1` voltage owned by EMT bus | $\mathbf{v}_{1} \in \mathbb{R}^N$
$\mathbf{v}_{2}$        | [V]   | Terminal `2` voltage owned by EMT bus | $\mathbf{v}_{2} \in \mathbb{R}^N$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}_{1}$ | `i1` | Output | [A] | Current injection at terminal `1` | $\mathbf{i}^{\mathrm{inj}}_{1} \in \mathbb{R}^N$
$\mathbf{i}^{\mathrm{inj}}_{2}$ | `i2` | Output | [A] | Current injection at terminal `2` | $\mathbf{i}^{\mathrm{inj}}_{2} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

None.

### Algebraic Equations

None.

### Wiring

```math
\begin{aligned}
\mathbf{i}^{\mathrm{sh}}_{1} &=
  \mathbf{y}_c*\mathbf{v}_{1} \\
\mathbf{i}^{\mathrm{sh}}_{2} &=
  \mathbf{y}_c*\mathbf{v}_{2} \\
\mathbf{i}^{\mathrm{inc}}_{1} &=
  \mathbf{h}*\left(
  2\mathbf{i}^{\mathrm{sh}}_{2}
  -
  \mathbf{i}^{\mathrm{inc}}_{2}
  \right) \\
\mathbf{i}^{\mathrm{inc}}_{2} &=
  \mathbf{h}*\left(
  2\mathbf{i}^{\mathrm{sh}}_{1}
  -
  \mathbf{i}^{\mathrm{inc}}_{1}
  \right) \\
\mathbf{i}^{\mathrm{inj}}_{1} &=
  \mathbf{P}_\phi\left(
  \mathbf{i}^{\mathrm{inc}}_{1}
  -
  \mathbf{i}^{\mathrm{sh}}_{1}
  \right) \\
\mathbf{i}^{\mathrm{inj}}_{2} &=
  \mathbf{P}_\phi
  \left(
  \mathbf{i}^{\mathrm{inc}}_{2}
  -
  \mathbf{i}^{\mathrm{sh}}_{2}
  \right)
\end{aligned}
```

## Initialization

Let subscript $0$ denote initial values. Initial values satisfy the wiring
equations:

```math
\begin{aligned}
\mathbf{i}^{\mathrm{sh}}_{1,0} &= \mathbf{y}_c*\mathbf{v}_{1,0} \\
\mathbf{i}^{\mathrm{sh}}_{2,0} &= \mathbf{y}_c*\mathbf{v}_{2,0} \\
\mathbf{i}^{\mathrm{inc}}_{1,0} &= \mathbf{h}*\left(2\mathbf{i}^{\mathrm{sh}}_{2,0} - \mathbf{i}^{\mathrm{inc}}_{2,0}\right) \\
\mathbf{i}^{\mathrm{inc}}_{2,0} &= \mathbf{h}*\left(2\mathbf{i}^{\mathrm{sh}}_{1,0} - \mathbf{i}^{\mathrm{inc}}_{1,0}\right)
\end{aligned}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i_sh_1` | [A] | Shunt current at terminal `1` | $\mathbf{i}^{\mathrm{sh}}_{1} \in \mathbb{R}^K$
`i_sh_2` | [A] | Shunt current at terminal `2` | $\mathbf{i}^{\mathrm{sh}}_{2} \in \mathbb{R}^K$
`i_inc_1` | [A] | Incident current at terminal `1` | $\mathbf{i}^{\mathrm{inc}}_{1} \in \mathbb{R}^K$
`i_inc_2` | [A] | Incident current at terminal `2` | $\mathbf{i}^{\mathrm{inc}}_{2} \in \mathbb{R}^K$
