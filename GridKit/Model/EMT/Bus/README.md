# Bus Model

`Bus` represents an $N$-phase bus in instantaneous phase coordinates. It owns
the bus voltage and contributes the current-balance residual to
the assembled DAE. $`\mathcal{D}`$ denotes the set of connected devices.
Bus voltage and its residual are algebraic when no connected model contributes
a voltage derivative; see [assembly](../README.md#assembly).

## Block Diagram

![Bus model block diagram](../../../../docs/Figures/EMT/Bus/diagram.png)

Figure 1: Bus model

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
$\mathbf{i}_d$ | `i` | Input | [A] | Current from connected device $d$ | One port per $d \in \mathcal{D}$, $\mathbf{i}_d \in \mathbb{R}^N$
$\mathbf{v}$ | `v` | Output | [V] | Bus voltage supplied to connected devices | $\mathbf{v} \in \mathbb{R}^N$

## Submodels

None.

### Submodel Validation

None.

### Submodel Wiring

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

## Model Equations

### Internal Equations

#### Differential

```math
0 = \sum_{d \in \mathcal{D}} \mathbf{i}_d
```

#### Algebraic

None.

### External Equations

None.

## Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`v` | [V] | Bus voltage | $\mathbf{v} \in \mathbb{R}^N$
