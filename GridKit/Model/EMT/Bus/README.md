# Bus Model

`Bus` represents an $N$-phase bus in instantaneous phase coordinates. The bus
voltage is differential, and the model equations enforce current balance.
$\mathcal{E}$ denotes the set of connected devices.

> [!NOTE]
> A template parameter will select whether $\mathbf{v}$ is a differential or
> algebraic vector. This page documents the nondegenerate differential
> formulation. Initial end-to-end support is three-phase, although the
> formulation remains $N$-phase.

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

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Bus voltage vector | $\mathbf{v} \in \mathbb{R}^N$

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

None.

## Model Ports

Symbol | Type | Units | Description | Note
------ | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}_e$ | Input | [A] | Current injection from connected device $e$ | $e \in \mathcal{E}$, $\mathbf{i}^{\mathrm{inj}}_e \in \mathbb{R}^N$

## Model Equations

### Differential Equations

```math
0 = \sum_{e \in \mathcal{E}} \mathbf{i}^{\mathrm{inj}}_e
```

### Algebraic Equations

None.

### Wiring

None.

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{i}^{\mathrm{inj}}_e
  &\leftarrow \text{connected-device current},
     \quad e \in \mathcal{E}
\end{aligned}
```

### Internal Initialization

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\mathbf{v}_0$ | [V] | `init.v0` | Initial bus-voltage vector | Required, $\mathbf{v}_0 \in \mathbb{R}^N$

```math
\mathbf{v} \leftarrow \mathbf{v}_0
```

### Output Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`v` | [V] | Bus voltage | $\mathbf{v} \in \mathbb{R}^N$
