# VoltageSource Model

`VoltageSource` represents an $N$-phase sinusoidal EMT voltage source connected
to the EMT bus through terminal admittance.

## Block Diagram

![VoltageSource model block diagram](../../../../../../docs/Figures/EMT/VoltageSource/diagram.png)

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

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `v` | Input | [V] | Bus voltage at source port | $\mathbf{v} \in \mathbb{R}^N$
$\mathbf{i}$ | `i` | Output | [A] | Current injection at source port | $\mathbf{i} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

None.

### Algebraic Equations

```math
0 = e_n - \sqrt{2}E_n\cos(\omega t + \phi_n),
\quad n \in \mathcal{N}
```

### Wiring

```math
\mathbf{i} \leftarrow \mathbf{y}[\mathbf{e} - \mathbf{v}]
```

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{V}
  &\leftarrow \text{solved terminal RMS voltage phasor} \\
\mathbf{v}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\mathbf{V}) \\
\dfrac{\mathrm{d}\mathbf{v}}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(s_0\mathbf{V}).
\end{aligned}
```

### Internal Initialization

Sinusoidal steady-state initialization requires $\omega=\omega_0$. Define the
source-voltage phasor by

```math
V_n^\mathrm{s}=E_ne^{\mathrm{j}\phi_n},
\quad n \in \mathcal{N}.
```

The terminal-admittance submodel initializes from
$\mathbf{V}^\mathrm{s}-\mathbf{V}$. At $t=0$,

```math
\begin{aligned}
\mathbf{e}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\mathbf{V}^\mathrm{s}) \\
\dfrac{\mathrm{d}\mathbf{e}}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(s_0\mathbf{V}^\mathrm{s}).
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
\mathbf{I}
  &= \mathbf{Y}(s_0)(\mathbf{V}^\mathrm{s}-\mathbf{V}) \\
\mathbf{i}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\mathbf{I}) \\
\dfrac{\mathrm{d}\mathbf{i}}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(s_0\mathbf{I}).
\end{aligned}
```

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
