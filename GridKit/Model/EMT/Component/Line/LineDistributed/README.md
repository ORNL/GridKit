# LineDistributed Model

`LineDistributed` represents an $N$-phase, $K$-conductor distributed EMT line.

> [!NOTE]
> The initial end-to-end implementation will support three-phase systems only
> to establish a proof of concept. The formulation below remains $N$-phase.

## Block Diagram

![LineDistributed model block diagram](../../../../../../docs/Figures/EMT/LineDistributed/diagram.png)

Figure 1: LineDistributed model

The conductor-to-phase mappings are shown in the equations and omitted from
the diagram for clarity.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer
$K$ | [-] | `K` | Number of conductors | Required, positive integer
$\mathbf{P}_\phi$ | [-] | `conductors` | Conductor-to-phase assignment matrix | $\mathbf{P}_\phi \in \{0,1\}^{N \times K}$

### Parameter Validation

```math
\begin{aligned}
N &> 0 \\
K &> 0 \\
\mathbf{P}_\phi &\in \{0,1\}^{N \times K} \\
\sum_{n = 1}^{N}P_{\phi,nk} &= 1,
  \quad k \in \{1,\ldots,K\}
\end{aligned}
```

### Derived Parameters

None.

## Submodels

Submodel | Type | Order | Inputs | JSON | Outputs
-------- | ---- | ----- | ------ | ---- | -------
Terminal-1 characteristic admittance $f^{\mathbf{y}_c}_{1}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}_c}$ | $\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{1} \in \mathbb{R}^K$ | `Yc` | $\mathbf{i}^{\mathrm{sh}}_{1} \in \mathbb{R}^K$
Terminal-2 characteristic admittance $f^{\mathbf{y}_c}_{2}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}_c}$ | $\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{2} \in \mathbb{R}^K$ | `Yc` | $\mathbf{i}^{\mathrm{sh}}_{2} \in \mathbb{R}^K$
Propagation $2\rightarrow1$, $f^{\mathbf{h}}_{2\rightarrow1}$ | [`Propagation`](../../../Operators/Shift/Propagation/README.md) | $Q_{\mathbf{h}}$ | $2\mathbf{i}^{\mathrm{sh}}_{2} - \mathbf{i}^{\mathrm{inc}}_{2} \in \mathbb{R}^K$ | `H` | $\mathbf{i}^{\mathrm{inc}}_{1} \in \mathbb{R}^K$
Propagation $1\rightarrow2$, $f^{\mathbf{h}}_{1\rightarrow2}$ | [`Propagation`](../../../Operators/Shift/Propagation/README.md) | $Q_{\mathbf{h}}$ | $2\mathbf{i}^{\mathrm{sh}}_{1} - \mathbf{i}^{\mathrm{inc}}_{1} \in \mathbb{R}^K$ | `H` | $\mathbf{i}^{\mathrm{inc}}_{2} \in \mathbb{R}^K$

### Submodel Validation

```math
\begin{aligned}
f^{\mathbf{y}_c}_{1}, f^{\mathbf{y}_c}_{2} &: \mathbb{R}^K \rightarrow \mathbb{R}^K \\
f^{\mathbf{h}}_{2\rightarrow1}, f^{\mathbf{h}}_{1\rightarrow2}
  &: \mathbb{R}^K \rightarrow \mathbb{R}^K
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{sh}}_{1}$ | [A] | Shunt current at terminal `1` | $\mathbf{i}^{\mathrm{sh}}_{1} \in \mathbb{R}^K$
$\mathbf{i}^{\mathrm{sh}}_{2}$ | [A] | Shunt current at terminal `2` | $\mathbf{i}^{\mathrm{sh}}_{2} \in \mathbb{R}^K$
$\mathbf{i}^{\mathrm{inc}}_{1}$ | [A] | Incident current at terminal `1` | $\mathbf{i}^{\mathrm{inc}}_{1} \in \mathbb{R}^K$
$\mathbf{i}^{\mathrm{inc}}_{2}$ | [A] | Incident current at terminal `2` | $\mathbf{i}^{\mathrm{inc}}_{2} \in \mathbb{R}^K$

### External Variables

External variables are owned by the EMT bus.

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}_{1}$ | [V] | Terminal `1` voltage owned by EMT bus | $\mathbf{v}_{1} \in \mathbb{R}^N$
$\mathbf{v}_{2}$ | [V] | Terminal `2` voltage owned by EMT bus | $\mathbf{v}_{2} \in \mathbb{R}^N$

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
  f^{\mathbf{y}_c}_{1}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{1}\right) \\
\mathbf{i}^{\mathrm{sh}}_{2} &=
  f^{\mathbf{y}_c}_{2}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{2}\right) \\
\mathbf{i}^{\mathrm{inc}}_{1} &=
  f^{\mathbf{h}}_{2\rightarrow1}\left(
  2\mathbf{i}^{\mathrm{sh}}_{2}
  -
  \mathbf{i}^{\mathrm{inc}}_{2}
  \right) \\
\mathbf{i}^{\mathrm{inc}}_{2} &=
  f^{\mathbf{h}}_{1\rightarrow2}\left(
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
\mathbf{v}_{1},\mathbf{v}_{2},
\dfrac{\mathrm{d}\mathbf{v}_{1}}{\mathrm{d}t},
\dfrac{\mathrm{d}\mathbf{v}_{2}}{\mathrm{d}t}
  &\leftarrow \text{provisional terminal-bus voltage trajectories}
\end{aligned}
```

### Internal Initialization

The characteristic-admittance submodels are initialized from the terminal-bus
voltage trajectories:

```math
\begin{aligned}
\mathbf{i}^{\mathrm{sh}}_{1}
  &\leftarrow
  f^{\mathbf{y}_c}_{1}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{1}\right) \\
\dfrac{\mathrm{d}\mathbf{i}^{\mathrm{sh}}_{1}}{\mathrm{d}t}
  &\leftarrow
  \dfrac{\mathrm{d}}{\mathrm{d}t}
  \left[
    f^{\mathbf{y}_c}_{1}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{1}\right)
  \right] \\
\mathbf{i}^{\mathrm{sh}}_{2}
  &\leftarrow
  f^{\mathbf{y}_c}_{2}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{2}\right) \\
\dfrac{\mathrm{d}\mathbf{i}^{\mathrm{sh}}_{2}}{\mathrm{d}t}
  &\leftarrow
  \dfrac{\mathrm{d}}{\mathrm{d}t}
  \left[
    f^{\mathbf{y}_c}_{2}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{2}\right)
  \right]
\end{aligned}
```

The line prehistory or network operating point supplies both incident-current
trajectories. The coupled relations below enforce their consistency but do not,
in general, select a unique history from the terminal voltages alone. Neither
propagation direction is evaluated first, and each submodel receives the value
and time derivative of its displayed input.

```math
\begin{aligned}
\left(
  \mathbf{i}^{\mathrm{inc}}_{1},
  \dfrac{\mathrm{d}\mathbf{i}^{\mathrm{inc}}_{1}}{\mathrm{d}t},
  \mathbf{i}^{\mathrm{inc}}_{2},
  \dfrac{\mathrm{d}\mathbf{i}^{\mathrm{inc}}_{2}}{\mathrm{d}t}
\right)
  &\leftarrow \text{line prehistory or network operating point} \\
\begin{bmatrix}
\mathbf{i}^{\mathrm{inc}}_{1}
  - f^{\mathbf{h}}_{2\rightarrow1}\left(
      2\mathbf{i}^{\mathrm{sh}}_{2} - \mathbf{i}^{\mathrm{inc}}_{2}
    \right) \\
\mathbf{i}^{\mathrm{inc}}_{2}
  - f^{\mathbf{h}}_{1\rightarrow2}\left(
      2\mathbf{i}^{\mathrm{sh}}_{1} - \mathbf{i}^{\mathrm{inc}}_{1}
    \right)
\end{bmatrix}
  &\leftarrow \mathbf{0}
\end{aligned}
```

The differentiated propagation relations are enforced simultaneously:

```math
\begin{bmatrix}
\dfrac{\mathrm{d}\mathbf{i}^{\mathrm{inc}}_{1}}{\mathrm{d}t}
  - \dfrac{\mathrm{d}}{\mathrm{d}t}
    \left[
      f^{\mathbf{h}}_{2\rightarrow1}\left(
        2\mathbf{i}^{\mathrm{sh}}_{2} - \mathbf{i}^{\mathrm{inc}}_{2}
      \right)
    \right] \\
\dfrac{\mathrm{d}\mathbf{i}^{\mathrm{inc}}_{2}}{\mathrm{d}t}
  - \dfrac{\mathrm{d}}{\mathrm{d}t}
    \left[
      f^{\mathbf{h}}_{1\rightarrow2}\left(
        2\mathbf{i}^{\mathrm{sh}}_{1} - \mathbf{i}^{\mathrm{inc}}_{1}
      \right)
    \right]
\end{bmatrix}
  \leftarrow \mathbf{0}
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
