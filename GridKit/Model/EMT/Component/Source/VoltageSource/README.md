# VoltageSource Model

`VoltageSource` represents an $N$-phase sinusoidal EMT voltage source connected
to the EMT bus through terminal admittance.

> [!NOTE]
> The initial end-to-end implementation will support three-phase systems only
> to establish a proof of concept. The formulation below remains $N$-phase.

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

Submodel | Type | Order | Inputs | Parameters | Outputs
-------- | ---- | ----- | ------ | ---------- | -------
Terminal admittance $f^{\mathbf{y}}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}}$ | $\mathbf{v} - \mathbf{e} \in \mathbb{R}^N$ | `Y` | $f^{\mathbf{y}}(\mathbf{v} - \mathbf{e}) \in \mathbb{R}^N$

### Submodel Validation

```math
\begin{aligned}
f^{\mathbf{y}} &: \mathbb{R}^N \rightarrow \mathbb{R}^N
\end{aligned}
```

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

External variables are owned by the EMT bus.

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
  &\leftarrow \text{provisional bus-voltage trajectory}
\end{aligned}
```

### Internal Initialization

At the initial simulation time, the analytic source trajectory is initialized
first:

```math
\begin{aligned}
e_n
  &\leftarrow
  \sqrt{2}E_n\cos\left(\phi_n\right),
  \quad n \in \{1,\ldots,N\} \\
\dfrac{\mathrm{d}e_n}{\mathrm{d}t}
  &\leftarrow
  -\sqrt{2}\omega E_n\sin\left(\phi_n\right),
  \quad n \in \{1,\ldots,N\}
\end{aligned}
```

The terminal-admittance submodel is then initialized with input value
$\mathbf{v} - \mathbf{e}$ and input derivative
$\mathrm{d}\mathbf{v}/\mathrm{d}t - \mathrm{d}\mathbf{e}/\mathrm{d}t$. Its
initialized output determines the algebraic source current:

```math
\mathbf{i} \leftarrow -f^{\mathbf{y}}(\mathbf{v} - \mathbf{e})
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
