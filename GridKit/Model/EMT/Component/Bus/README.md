# Bus Model

`Bus` owns its voltage, current balance, and shunt admittances. Positive
current denotes injection into the bus. Norton terminals additionally expose
their shunt currents to connected devices.

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
$\mathbf{i}_e^\mathrm{inc}$ | `i_inc` | Input | [A] | Incident current from device $e$ | Norton terminals, $\mathbb{R}^{K_e}$
$\mathbf{i}_e^\mathrm{sh}$ | `Ish` | Output | [A] | Shunt current supplied to device $e$ | Norton terminals, $\mathbb{R}^{K_e}$
$\mathbf{v}$ | `v` | Output | [V] | Bus voltage supplied to connected devices | $\mathbf{v} \in \mathbb{R}^N$

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{v}$ | Bus voltage and current balance | KCL | $N$ | — | Current contributions | $\mathbb{R}^N$
$\mathbf{n}_e$ | Norton source | [Norton](../Source/Norton/README.md) | $K_e(1+Q_e)$ | Device coefficients | $\mathbf{P}_e^\mathsf{T}\mathbf{v}$ | $\mathbf{i}_e^\mathrm{sh}$
$\mathbf{y}_e$ | Shunt admittance | [VectorFit](../../Operators/Rational/VectorFit/README.md) | $K_eQ_e$ | Device coefficients or `shunts.<name>` | $\mathbf{P}_e^\mathsf{T}\mathbf{v}$ | $\mathbb{R}^{K_e}$

### Submodel Validation

Each admittance satisfies its coefficient constraints.

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

```math
\begin{aligned}
0 &= -\mathbf{i}_e^\mathrm{sh}
     +\mathbf{y}_e[\mathbf{P}_e^\mathsf{T}\mathbf{v}],
     \quad e\in\mathcal{E} \\
0 &= \sum_{r\in\mathcal{R}}\sigma_r\mathbf{P}_r\mathbf{i}_r
     -\sum_{e\in\mathcal{S}}\mathbf{P}_e
       \mathbf{y}_e[\mathbf{P}_e^\mathsf{T}\mathbf{v}]
\end{aligned}
```

$\mathcal{E}$ contains Norton terminals; $\mathcal{S}$ contains direct shunts.
$\mathcal{R}$ contains registered current signals with sign $\sigma_r$ and
phase map $\mathbf{P}_r$. Norton incident and shunt currents have signs $+1$
and $-1$, respectively. Bus voltage is algebraic unless a connected equation
depends on its derivative.

## Initialization

Initialize the voltage from the bus state entry, then initialize the
admittance states and Norton shunt currents.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`v` | [V] | Bus voltage | $\mathbf{v} \in \mathbb{R}^N$
`i_sh` | [A] | Total shunt current | Includes direct shunts and Norton terminals
