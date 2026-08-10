# LoadZ Model

`LoadZ` represents an $N$-phase impedance load. Current $\mathbf{i}$ is injected
from the load into the EMT bus.

## Block Diagram

![LoadZ model block diagram](../../../../../../docs/Figures/EMT/LoadZ/diagram.png)

Figure 1: LoadZ model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer

### Parameter Validation

```math
N \in \mathbb{Z}_{>0}
```

### Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `v` | Input | [V] | Bus voltage at load port | $\mathbf{v} \in \mathbb{R}^N$
$\mathbf{i}$ | `i` | Output | [A] | Current injection at load port | $\mathbf{i} \in \mathbb{R}^N$

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{z}$ | Impedance | [VectorFit](../../../Operators/Rational/VectorFit/README.md) | $NQ_{\mathbf{z}}$ | `Z` | $\mathbb{R}^N$ | $\mathbb{R}^N$

### Submodel Validation

The current is differential when $\mathbf{E}^{\mathbf{z}}$ is nonsingular and
algebraic when $\mathbf{E}^{\mathbf{z}}=\mathbf{0}$.

```math
\mathbf{E}^{\mathbf{z}}=\mathbf{0}
\qquad \text{or} \qquad
\mathrm{rank}(\mathbf{E}^{\mathbf{z}})=N.
```

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$ | [A] | Current injection from load into EMT bus | $\mathbf{i} \in \mathbb{R}^N$, $\mathrm{rank}(\mathbf{E}^{\mathbf{z}})=N$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$ | [A] | Current injection from load into EMT bus | $\mathbf{i} \in \mathbb{R}^N$, $\mathbf{E}^{\mathbf{z}}=\mathbf{0}$

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Bus voltage vector owned by EMT bus | $\mathbf{v} \in \mathbb{R}^N$

#### Algebraic

None.

## Model Equations

### Internal Equations

#### Differential

For $\mathrm{rank}(\mathbf{E}^{\mathbf{z}})=N$,

```math
0 = \mathbf{z}[\mathbf{i}] + \mathbf{v}
```

#### Algebraic

For $\mathbf{E}^{\mathbf{z}}=\mathbf{0}$, the same residual is algebraic.

### External Equations

```math
\mathbf{f} \leftarrow \mathbf{i}
```

## Initialization

None beyond the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Load current injection | $\mathbf{i} \in \mathbb{R}^N$

## Development

The initial three-phase formulation uses resistance and inductance matrices.

### Differential Equations

```math
0 =
\mathbf{R}\mathbf{i}
+ \mathbf{L}\dfrac{\mathrm{d}\mathbf{i}}{\mathrm{d}t}
+ \mathbf{v}
```
