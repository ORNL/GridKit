# Bus Model

`Bus` represents an $N$-phase bus in instantaneous phase coordinates. It owns
the differential bus voltage and contributes the current-balance residual to
the assembled DAE. $\mathcal{E}$ denotes the set of connected devices.

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
$\mathbf{i}_e$ | `i` | Input | [A] | Current from connected device $e$ | One port per $e \in \mathcal{E}$, $\mathbf{i}_e \in \mathbb{R}^N$
$\mathbf{v}$ | `v` | Output | [V] | Bus voltage supplied to connected devices | $\mathbf{v} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

```math
0 = \sum_{e \in \mathcal{E}} \mathbf{i}_e
```

### Algebraic Equations

None.

### Wiring

None.

## Initialization

None beyond the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`v` | [V] | Bus voltage | $\mathbf{v} \in \mathbb{R}^N$
