# VoltageSource Model

`VoltageSource` represents an $N$-phase sinusoidal EMT voltage source connected
to the EMT bus through terminal admittance.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer
$\mathbf{E}$                   | [V]     | `E`      | Source voltage magnitudes         | $\mathbf{E}\in\mathbb{R}^N$, RMS
$\boldsymbol{\phi}$            | [rad]   | `phi`    | Source phase offsets              | $\boldsymbol{\phi}\in\mathbb{R}^N$
$\omega$                       | [rad/s] | `omega`  | Source angular frequency          |

### Parameter Validation

```math
\begin{aligned}
N &> 0 \\
\mathbf{E} &\ge \mathbf{0} \\
\boldsymbol{\phi} &\in \mathbb{R}^N \\
\omega &> 0
\end{aligned}
```

### Derived Parameters

None.

## Submodels

Submodel | Type | Order | Inputs | Outputs
-------- | ---- | ----- | ------ | -------
Terminal admittance $f^{\mathbf{y}}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}}$ | $\mathbf{v}-\mathbf{e}\in\mathbb{R}^N$ | $f^{\mathbf{y}}(\mathbf{v}-\mathbf{e})\in\mathbb{R}^N$

### Submodel Validation

```math
\begin{aligned}
f^{\mathbf{y}} &: \mathbb{R}^N \rightarrow \mathbb{R}^N
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$ | [A] | Current injection from source into EMT bus | $\mathbf{i}\in\mathbb{R}^N$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{e}$ | [V] | Source voltage vector | $\mathbf{e}\in\mathbb{R}^N$

### External Variables

External variables are owned by the EMT bus.

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$     | [V]   | Bus voltage vector                          | $\mathbf{v} \in \mathbb{R}^N$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}$ | `i` | Output | [A] | Current injection at source port | $\mathbf{i}^{\mathrm{inj}} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

```math
0 = \mathbf{i} + f^{\mathbf{y}}(\mathbf{v}-\mathbf{e})
```

### Algebraic Equations

```math
0
=
e_n
-
\sqrt{2}E_n\cos\left(\omega t + \phi_n\right),
\quad n\in\{0,1,\ldots,N-1\}
```

### Wiring

```math
\mathbf{i}^{\mathrm{inj}}
=
\mathbf{i}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{v}
  &\leftarrow \text{initialized bus voltage}
\end{aligned}
```

### Internal Initialization

Initialization is performed by evaluating the source assignments in dependency
order:

```math
\begin{aligned}
e_n
  &\leftarrow
  \sqrt{2}E_n\cos\left(\phi_n\right),
  \quad n\in\{0,1,\ldots,N-1\} \\
\mathbf{i}
  &\leftarrow \text{source-current start} \\
\mathbf{i} + f^{\mathbf{y}}(\mathbf{v}-\mathbf{e})
  &\leftarrow \mathbf{0}
\end{aligned}
```

### Output Initialization

```math
\mathbf{i}^{\mathrm{inj}} \leftarrow \mathbf{i}
```

## Monitors

None.
