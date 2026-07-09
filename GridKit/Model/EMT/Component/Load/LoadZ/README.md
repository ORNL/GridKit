# LoadZ Model

`LoadZ` represents an $N$-phase impedance load. Current $\mathbf{i}$ is injected
from the load into the EMT bus.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer

### Parameter Validation

```math
\begin{aligned}
N &> 0
\end{aligned}
```

### Derived Parameters

None.

## Submodels

Submodel | Type | Order | Inputs | Outputs
-------- | ---- | ----- | ------ | -------
Impedance $f^{\mathbf{z}}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{z}}$ | $\mathbf{i}\in\mathbb{R}^N$ | $f^{\mathbf{z}}(\mathbf{i})\in\mathbb{R}^N$

### Submodel Validation

```math
\begin{aligned}
f^{\mathbf{z}} &: \mathbb{R}^N \rightarrow \mathbb{R}^N
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$        | [A]   | Current injection from load into EMT bus       | $\mathbf{i} \in \mathbb{R}^N$

#### Algebraic

None.

### External Variables

External variables are owned by the EMT bus.

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$     | [V]   | Port voltage vector                         | $\mathbf{v} \in \mathbb{R}^N$

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
  &\leftarrow \text{initialized bus voltage}
\end{aligned}
```

### Internal Initialization

Initialization sets the current start by enforcing the differential residual.

```math
\begin{aligned}
\mathbf{i}
  &\leftarrow \text{load-current start} \\
f^{\mathbf{z}}(\mathbf{i}) + \mathbf{v}
  &\leftarrow \mathbf{0}
\end{aligned}
```

### Output Initialization

```math
\mathbf{i}^{\mathrm{inj}} \leftarrow \mathbf{i}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Load current injection | $\mathbf{i} \in \mathbb{R}^N$
