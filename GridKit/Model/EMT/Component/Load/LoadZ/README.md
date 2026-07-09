# LoadZ Model

`LoadZ` represents an $N$-phase impedance load. Current $\mathbf{i}$ is injected
from the load into the EMT bus.

> [!NOTE]
> A template parameter will select whether $\mathbf{i}$ is a differential or
> algebraic vector. This page documents the nondegenerate differential
> formulation.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer

### Parameter Validation

```math
N > 0
```

### Derived Parameters

None.

## Submodels

Submodel | Type | Order | Inputs | JSON | Outputs
-------- | ---- | ----- | ------ | ---- | -------
Impedance $f^{\mathbf{z}}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{z}}$ | $\mathbb{R}^N$ [A] | `Z` | $\mathbb{R}^N$ [V]

### Submodel Validation

```math
\mathrm{rank}\left(\mathbf{E}^{\mathbf{z}}\right) = N
```

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$ | [A] | Current injection from load into EMT bus | $\mathbf{i} \in \mathbb{R}^N$

#### Algebraic

None.

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Bus voltage vector owned by EMT bus | $\mathbf{v} \in \mathbb{R}^N$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}$ | `i` | Output | [A] | Current injection at load port | $\mathbf{i}^{\mathrm{inj}} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

```math
0 = f^{\mathbf{z}}(\mathbf{i}) + \mathbf{v}
```

### Algebraic Equations

None.

### Wiring

```math
\mathbf{i}^{\mathrm{inj}} = \mathbf{i}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{v}
  &\leftarrow \text{bus voltage}
\end{aligned}
```

### Internal Initialization

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\mathbf{i}_0$ | [A] | `init.i0` | Initial current-injection vector | Required, $\mathbf{i}_0 \in \mathbb{R}^N$

```math
\mathbf{i} \leftarrow \mathbf{i}_0
```

### Output Initialization

```math
\mathbf{i}^{\mathrm{inj}} \leftarrow \mathbf{i}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Load current injection | $\mathbf{i} \in \mathbb{R}^N$
