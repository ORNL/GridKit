# Bus Model

`Bus` is a container of `KCL` and Norton sources. `KCL` owns the bus voltage;
each Norton source owns its shunt current and admittance states. KCL owns the
sum of all registered terminal currents and its Jacobian. Positive current
means injection into the bus.
$\mathcal{E}$ denotes the set of Norton sources.

## Block Diagram

![Bus model block diagram](../../../../../docs/Figures/EMT/Bus/diagram.png)

Figure 1: Bus model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Optional; defaults to 3, the supported phase count

### Parameter Validation

```math
N = 3
```

### Derived Parameters

Terminal $e$ has $K_e$ current channels and a conductor-to-phase map
$\mathbf{P}_e\in\mathbb{R}^{N\times K_e}$ supplied by the connected device.
For a direct phase connection, $K_e=N$ and $\mathbf{P}_e=\mathbf{I}_N$.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}_e^\mathrm{inc}$ | `i_inc` | Input | [A] | Incident current from device $e$ | One per terminal, $\mathbb{R}^{K_e}$
$\mathbf{i}_e^\mathrm{sh}$ | `Ish` | Output | [A] | Shunt current supplied to device $e$ | One per terminal, $\mathbb{R}^{K_e}$
$\mathbf{v}$ | `v` | Output | [V] | Bus voltage supplied to connected devices | $\mathbf{v} \in \mathbb{R}^N$

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{v}$ | Bus voltage and current balance | KCL | $N$ | — | Current contributions | $\mathbb{R}^N$
$\mathbf{n}_e$ | Norton source | [Norton](../Source/Norton/README.md) | $K_e(1+Q_e)$ | Device coefficients or `shunts.<name>` | $\mathbf{P}_e^\mathsf{T}\mathbf{v}$ | $\mathbf{i}_e^\mathrm{sh}$

### Submodel Validation

Each source satisfies its admittance coefficient constraints.

## Model Variables

### Internal Variables

#### Differential

None.

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

None.

#### Algebraic

None.

### External Equations

None.

### Submodel Equations

The Norton sources own their shunt-current rows; KCL owns the voltage rows:

```math
\begin{aligned}
0 &= -\mathbf{i}_e^\mathrm{sh}
     +\mathbf{y}_e[\mathbf{P}_e^\mathsf{T}\mathbf{v}],
     \quad e\in\mathcal{E} \\
0 &= \sum_{r\in\mathcal{R}}\sigma_r\mathbf{P}_r\mathbf{i}_r
\end{aligned}
```

$\mathcal{R}$ contains all current registrations, with sign $\sigma_r$ and
phase map $\mathbf{P}_r$. Norton incident and shunt currents have signs $+1$
and $-1$, respectively. Bus voltage is algebraic unless a connected equation
depends on its derivative.

Current signals are registered before allocation. KCL evaluates their signed
sum and signal gradients; connected devices do not stamp bus equations.

## Initialization

`KCL` initializes the voltage from the bus state entry. Each Norton source
initializes its admittance states and shunt current.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`v` | [V] | Bus voltage | $\mathbf{v} \in \mathbb{R}^N$
`i_sh` | [A] | Total shunt current | $\sum_{e\in\mathcal{E}}\mathbf{P}_e\mathbf{i}_e^\mathrm{sh}$
