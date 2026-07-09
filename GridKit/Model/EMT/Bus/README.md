# Bus Model

`Bus` represents an $N$-phase bus in instantaneous phase coordinates. The
bus voltages are differential variables, and the model equations enforce
current balance at the bus.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$M$ | [-] | `M` | Number of connected-device ports | Required, positive integer
$N$ | [-] | `N` | Number of phases | Required, positive integer
$\mathbf{v}_0$ | [V] | `v0` | Initial bus voltage vector | $\mathbf{v}_0 \in \mathbb{R}^N$

### Parameter Validation

```math
\begin{aligned}
M &> 0 \\
N &> 0 \\
\mathbf{v}_0 &\in \mathbb{R}^N
\end{aligned}
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

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}_m$ | `i` | Input | [A] | Current injection from connected device $m$ | $m=1,\ldots,M$, $\mathbf{i}^{\mathrm{inj}}_m \in \mathbb{R}^N$

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= \sum_{m=1}^{M} \mathbf{i}^{\mathrm{inj}}_m
\end{aligned}
```

Each $\mathbf{i}^{\mathrm{inj}}_m$ may depend on the bus voltage and bus voltage derivative.

### Algebraic Equations

None.

### Wiring

None.

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{i}^{\mathrm{inj}}_m
  &\leftarrow \text{connected-device current starts},
     \quad m\in\{1,\ldots,M\}
\end{aligned}
```

### Internal Initialization

Initialization assigns the parameterized bus-voltage start:

```math
\mathbf{v} \leftarrow \text{parameterized bus-voltage start}
```

The differential residual determines the initial bus-voltage derivative.

### Output Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`v` | [V] | Bus voltage | $\mathbf{v} \in \mathbb{R}^N$
