# LineDistributed Model

`LineDistributed` represents an $N$-phase, $K$-conductor distributed EMT line.

## Block Diagram

![](../../../../../../docs/Figures/EMT/LineDistributed/diagram.png)

Figure 1: LineDistributed model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer
$K$ | [-] | `K` | Number of conductors | Required, positive integer
$\mathbf{P}_\phi$ | [-]   | `conductors` | Permutation matrix mapping each conductor to its phase | $\mathbf{P}_\phi \in \mathbb{R}^{N \times K}$

### Parameter Validation

```math
\begin{aligned}
N &> 0 \\
K &> 0 \\
\mathbf{P}_\phi &\in \mathbb{R}^{N \times K}
\end{aligned}
```

### Derived Parameters

None.

## Submodels

Submodel | Type | Order | Inputs | Outputs
-------- | ---- | ----- | ------ | -------
Characteristic admittance $f^{\mathbf{y}_c}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}_c}$ | $\mathbf{v}_{1},\mathbf{v}_{2}\in\mathbb{R}^N$ | $\mathbf{i}^{\mathrm{sh}}_{1},\mathbf{i}^{\mathrm{sh}}_{2}\in\mathbb{R}^K$
Propagation $f^{\mathbf{h}}$ | [`Propagation`](../../../Operators/Shift/Propagation/README.md) | $Q_{\mathbf{h}}$ | $\mathbf{u}_1,\mathbf{u}_2\in\mathbb{R}^K$ | $\mathbf{i}^{\mathrm{inc}}_{1},\mathbf{i}^{\mathrm{inc}}_{2}\in\mathbb{R}^K$

### Submodel Validation

```math
\begin{aligned}
f^{\mathbf{y}_c} &: \mathbb{R}^N \rightarrow \mathbb{R}^K \\
f^{\mathbf{h}} &: \mathbb{R}^K \rightarrow \mathbb{R}^K
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{sh}}_{1}$ | [A] | Shunt current at terminal `1` | $\mathbf{i}^{\mathrm{sh}}_{1}\in\mathbb{R}^K$
$\mathbf{i}^{\mathrm{sh}}_{2}$ | [A] | Shunt current at terminal `2` | $\mathbf{i}^{\mathrm{sh}}_{2}\in\mathbb{R}^K$
$\mathbf{i}^{\mathrm{inc}}_{1}$ | [A] | Incident current at terminal `1` | $\mathbf{i}^{\mathrm{inc}}_{1}\in\mathbb{R}^K$
$\mathbf{i}^{\mathrm{inc}}_{2}$ | [A] | Incident current at terminal `2` | $\mathbf{i}^{\mathrm{inc}}_{2}\in\mathbb{R}^K$

### External Variables

External variables are owned by the EMT bus.

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}_{1}$         | [V]   | Terminal `1` voltage owned by EMT bus | $\mathbf{v}_{1} \in \mathbb{R}^N$
$\mathbf{v}_{2}$        | [V]   | Terminal `2` voltage owned by EMT bus | $\mathbf{v}_{2} \in \mathbb{R}^N$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}_{1}$ | `i1` | Output | [A] | Current injection at terminal `1` | $\mathbf{i}^{\mathrm{inj}}_{1} \in \mathbb{R}^N$
$\mathbf{i}^{\mathrm{inj}}_{2}$ | `i2` | Output | [A] | Current injection at terminal `2` | $\mathbf{i}^{\mathrm{inj}}_{2} \in \mathbb{R}^N$

## Model Equations

### Differential Equations

None.

### Algebraic Equations

None.

### Wiring

```math
\begin{aligned}
\mathbf{i}^{\mathrm{sh}}_{1} &=
  f^{\mathbf{y}_c}(\mathbf{v}_{1}) \\
\mathbf{i}^{\mathrm{sh}}_{2} &=
  f^{\mathbf{y}_c}(\mathbf{v}_{2}) \\
\mathbf{i}^{\mathrm{inc}}_{1} &=
  f^{\mathbf{h}}\left(
  2\mathbf{i}^{\mathrm{sh}}_{2}
  -
  \mathbf{i}^{\mathrm{inc}}_{2}
  \right) \\
\mathbf{i}^{\mathrm{inc}}_{2} &=
  f^{\mathbf{h}}\left(
  2\mathbf{i}^{\mathrm{sh}}_{1}
  -
  \mathbf{i}^{\mathrm{inc}}_{1}
  \right) \\
\mathbf{i}^{\mathrm{inj}}_{1} &=
  \mathbf{P}_\phi\left(
  \mathbf{i}^{\mathrm{inc}}_{1}
  -
  \mathbf{i}^{\mathrm{sh}}_{1}
  \right) \\
\mathbf{i}^{\mathrm{inj}}_{2} &=
  \mathbf{P}_\phi
  \left(
  \mathbf{i}^{\mathrm{inc}}_{2}
  -
  \mathbf{i}^{\mathrm{sh}}_{2}
  \right)
\end{aligned}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{v}_1,\mathbf{v}_2
  &\leftarrow \text{terminal-bus voltage}
\end{aligned}
```

### Internal Initialization

Initialization assigns the internal line currents in dependency order:

```math
\begin{aligned}
\mathbf{i}^{\mathrm{sh}}_{1}
  &\leftarrow f^{\mathbf{y}_c}(\mathbf{v}_{1}) \\
\mathbf{i}^{\mathrm{sh}}_{2}
  &\leftarrow f^{\mathbf{y}_c}(\mathbf{v}_{2}) \\
\mathbf{i}^{\mathrm{inc}}_{1}
  &\leftarrow
  f^{\mathbf{h}}\left(
    2\mathbf{i}^{\mathrm{sh}}_{2}
    -
    \mathbf{i}^{\mathrm{inc}}_{2}
  \right) \\
\mathbf{i}^{\mathrm{inc}}_{2}
  &\leftarrow
  f^{\mathbf{h}}\left(
    2\mathbf{i}^{\mathrm{sh}}_{1}
    -
    \mathbf{i}^{\mathrm{inc}}_{1}
  \right)
\end{aligned}
```

### Output Initialization

Current-injection ports are assigned from initialized internal line currents.

```math
\begin{aligned}
\mathbf{i}^{\mathrm{inj}}_{1}
  &\leftarrow
  \mathbf{P}_\phi\left(
  \mathbf{i}^{\mathrm{inc}}_{1}
  -
  \mathbf{i}^{\mathrm{sh}}_{1}
  \right) \\
\mathbf{i}^{\mathrm{inj}}_{2}
  &\leftarrow
  \mathbf{P}_\phi\left(
  \mathbf{i}^{\mathrm{inc}}_{2}
  -
  \mathbf{i}^{\mathrm{sh}}_{2}
  \right)
\end{aligned}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i_sh_1` | [A] | Shunt current at terminal `1` | $\mathbf{i}^{\mathrm{sh}}_{1} \in \mathbb{R}^K$
`i_sh_2` | [A] | Shunt current at terminal `2` | $\mathbf{i}^{\mathrm{sh}}_{2} \in \mathbb{R}^K$
`i_inc_1` | [A] | Incident current at terminal `1` | $\mathbf{i}^{\mathrm{inc}}_{1} \in \mathbb{R}^K$
`i_inc_2` | [A] | Incident current at terminal `2` | $\mathbf{i}^{\mathrm{inc}}_{2} \in \mathbb{R}^K$
