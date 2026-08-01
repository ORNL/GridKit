# VoltageSource Model

`VoltageSource` represents an $N$-phase sinusoidal EMT voltage source connected
to the EMT bus through terminal admittance.

## Block Diagram

![VoltageSource model block diagram](diagram.png)

Figure 1: VoltageSource model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer
$\mathbf{E}$ | [V] | `E` | Source voltage magnitudes | $\mathbf{E} \in \mathbb{R}_{\ge 0}^N$, RMS
$\boldsymbol{\phi}$ | [rad] | `phi` | Source phase offsets | $\boldsymbol{\phi} \in \mathbb{R}^N$
$\omega$ | [rad/s] | `omega` | Source angular frequency | Required, positive

### Parameter Validation

```math
\begin{aligned}
N &\in \mathbb{Z}_{>0} \\
\mathbf{E} &\in \mathbb{R}_{\ge 0}^N \\
\boldsymbol{\phi} &\in \mathbb{R}^N \\
\omega &> 0
\end{aligned}
```

### Derived Parameters

Define the phase-index set

```math
\mathcal{N} = \{1,\ldots,N\}.
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `v` | Input | [V] | Bus voltage at source port | $\mathbf{v} \in \mathbb{R}^N$
$\mathbf{i}$ | `i` | Output | [A] | Current injection at source port | $\mathbf{i} \in \mathbb{R}^N$

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{y}$ | Terminal admittance | [VectorFit](../../../Operators/Rational/VectorFit/README.md) | $NQ_{\mathbf{y}}$ | `Y` | $\mathbb{R}^N$ | $\mathbb{R}^N$

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{e}$ | [V] | Source voltage vector | $\mathbf{e} \in \mathbb{R}^N$

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

None.

#### Algebraic

```math
0 = e_n - \sqrt{2}E_n\cos(\omega t + \phi_n),
\quad n \in \mathcal{N}
```

### External Equations

```math
\mathbf{i} \leftarrow \mathbf{y}[\mathbf{e} - \mathbf{v}]
```

## Initialization

None beyond the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`e` | [V] | Source voltage | $\mathbf{e} \in \mathbb{R}^N$
`i` | [A] | Source current injection | $\mathbf{i} \in \mathbb{R}^N$

## Development

The initial three-phase formulation realizes the terminal admittance as a
series resistance and inductance, with $\mathbf{i}$ as a differential variable.

### Differential Equations

```math
0 =
\mathbf{R}_\mathrm{s}\mathbf{i}
+ \mathbf{L}_\mathrm{s}\dfrac{\mathrm{d}\mathbf{i}}{\mathrm{d}t}
+ \mathbf{v}
- \mathbf{e}
```
