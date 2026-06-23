# Bus Model

`Bus` represents an $N$-phase bus in instantaneous phase coordinates. The
bus voltages are differential variables, and the model equations enforce
current balance at the bus.

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

### Model Derived Parameters

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

## Initialization

The initial bus voltage is given by the parameter vector $\mathbf{v}_0$:

```math
\mathbf{v}(0)=\mathbf{v}_0
```

Only $\mathbf{v}(0)$ is initialized. The solver computes
$\dot{\mathbf{v}}(0)$ from the differential residual.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`v` | [V] | Bus voltage | $\mathbf{v} \in \mathbb{R}^N$
