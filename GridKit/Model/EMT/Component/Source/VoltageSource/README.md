# VoltageSource Model

`VoltageSource` represents an $N$-phase sinusoidal EMT voltage source connected
to the EMT bus through terminal admittance.

## Block Diagram

None.

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
N &> 0 \\
\mathbf{E} &\in \mathbb{R}_{\ge 0}^N \\
\boldsymbol{\phi} &\in \mathbb{R}^N \\
\omega &> 0
\end{aligned}
```

### Derived Parameters

None.

## Submodels

Submodel | Type | Order | Inputs | JSON | Outputs
-------- | ---- | ----- | ------ | ---- | -------
Terminal admittance $f^{\mathbf{y}}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}}$ | $\mathbb{R}^N$ [V] | `Y` | $\mathbb{R}^N$ [A]

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
$\mathbf{i}$ | [A] | Current injection from source into EMT bus | $\mathbf{i} \in \mathbb{R}^N$

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
$\mathbf{i}^{\mathrm{inj}}$ | `i` | Output | [A] | Current injection at source port | $\mathbf{i}^{\mathrm{inj}} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

None.

### Algebraic Equations

```math
\begin{aligned}
0 &= \mathbf{i} + f^{\mathbf{y}}(\mathbf{v} - \mathbf{e}) \\
0 &= e_n - \sqrt{2}E_n\cos\left(\omega t + \phi_n\right),
  \quad n \in \{1,\ldots,N\}
\end{aligned}
```

### Wiring

```math
\mathbf{i}^{\mathrm{inj}} = \mathbf{i}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{v},\dfrac{\mathrm{d}\mathbf{v}}{\mathrm{d}t}
  &\leftarrow \text{bus voltage and derivative}
\end{aligned}
```

### Internal Initialization

At initial time $t_0$, the terminal-admittance submodel uses
$\mathbf{v} - \mathbf{e}$ and its time derivative as its input.

```math
\begin{aligned}
e_n
  &\leftarrow
  \sqrt{2}E_n\cos\left(\omega t_0 + \phi_n\right),
  \quad n \in \{1,\ldots,N\} \\
\dfrac{\mathrm{d}e_n}{\mathrm{d}t}
  &\leftarrow
  -\sqrt{2}\omega E_n\sin\left(\omega t_0 + \phi_n\right),
  \quad n \in \{1,\ldots,N\} \\
\mathbf{i}
  &\leftarrow -f^{\mathbf{y}}(\mathbf{v} - \mathbf{e})
\end{aligned}
```

### Output Initialization

```math
\mathbf{i}^{\mathrm{inj}} \leftarrow \mathbf{i}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`e` | [V] | Source voltage | $\mathbf{e} \in \mathbb{R}^N$
`i` | [A] | Source current injection | $\mathbf{i} \in \mathbb{R}^N$
