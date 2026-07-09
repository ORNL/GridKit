# LineDistributed Model

`LineDistributed` represents an $N$-phase, $K$-conductor distributed EMT line.

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
$\mathbf{c}$ | [-] | `conductors` | Conductor phase-index list | $\mathbf{c} \in \{1,\ldots,N\}^{K}$

### Parameter Validation

```math
\begin{aligned}
N &> 0 \\
K &> 0 \\
\mathbf{c} &\in \{1,\ldots,N\}^{K} \\
\{c_k \mid k \in \{1,\ldots,K\}\} &= \{1,\ldots,N\}
\end{aligned}
```

### Derived Parameters

```math
P_{\phi,nk} =
\begin{cases}
1, & n = c_k \\
0, & n \ne c_k
\end{cases},
\quad n \in \{1,\ldots,N\},\quad k \in \{1,\ldots,K\}
```

## Submodels

Submodel | Type | Order | Inputs | JSON | Outputs
-------- | ---- | ----- | ------ | ---- | -------
Terminal-1 characteristic admittance $f^{\mathbf{y}_c}_{1}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}_c}$ | $\mathbb{R}^K$ [V] | `Yc` | $\mathbb{R}^K$ [A]
Terminal-2 characteristic admittance $f^{\mathbf{y}_c}_{2}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}_c}$ | $\mathbb{R}^K$ [V] | `Yc` | $\mathbb{R}^K$ [A]
Propagation $2\rightarrow1$, $f^{\mathbf{h}}_{2\rightarrow1}$ | [`Propagation`](../../../Operators/Shift/Propagation/README.md) | Composite | $\mathbb{R}^K$ [A] | `H` | $\mathbb{R}^K$ [A]
Propagation $1\rightarrow2$, $f^{\mathbf{h}}_{1\rightarrow2}$ | [`Propagation`](../../../Operators/Shift/Propagation/README.md) | Composite | $\mathbb{R}^K$ [A] | `H` | $\mathbb{R}^K$ [A]

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{sh}}_{1}$ | [A] | Shunt current at terminal 1 | $\mathbf{i}^{\mathrm{sh}}_{1} \in \mathbb{R}^K$
$\mathbf{i}^{\mathrm{sh}}_{2}$ | [A] | Shunt current at terminal 2 | $\mathbf{i}^{\mathrm{sh}}_{2} \in \mathbb{R}^K$
$\mathbf{i}^{\mathrm{inc}}_{1}$ | [A] | Incident current at terminal 1 | $\mathbf{i}^{\mathrm{inc}}_{1} \in \mathbb{R}^K$
$\mathbf{i}^{\mathrm{inc}}_{2}$ | [A] | Incident current at terminal 2 | $\mathbf{i}^{\mathrm{inc}}_{2} \in \mathbb{R}^K$

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}_{1}$ | [V] | Terminal 1 voltage owned by EMT bus | $\mathbf{v}_{1} \in \mathbb{R}^N$
$\mathbf{v}_{2}$ | [V] | Terminal 2 voltage owned by EMT bus | $\mathbf{v}_{2} \in \mathbb{R}^N$

#### Algebraic

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inj}}_{1}$ | `i1` | Output | [A] | Current injection at terminal 1 | $\mathbf{i}^{\mathrm{inj}}_{1} \in \mathbb{R}^N$
$\mathbf{i}^{\mathrm{inj}}_{2}$ | `i2` | Output | [A] | Current injection at terminal 2 | $\mathbf{i}^{\mathrm{inj}}_{2} \in \mathbb{R}^N$

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
  &\leftarrow \text{terminal-bus values and derivatives}
\end{aligned}
```

### Internal Initialization

The system initializer supplies non-JSON incident-current values and derivatives
that must be consistent with both `Propagation` submodels.

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{inc,init}}_{j}$ | [A] | Initial incident-current vector | Required; $j \in \{1,2\}$, $\mathbf{i}^{\mathrm{inc,init}}_{j} \in \mathbb{R}^K$
$\left(\dfrac{\mathrm{d}\mathbf{i}^{\mathrm{inc}}_{j}}{\mathrm{d}t}\right)^{\mathrm{init}}$ | [A/s] | Initial incident-current derivative | Required; $j \in \{1,2\}$, $\left(\mathrm{d}\mathbf{i}^{\mathrm{inc}}_{j}/\mathrm{d}t\right)^{\mathrm{init}} \in \mathbb{R}^K$

```math
\begin{aligned}
\mathbf{i}^{\mathrm{sh}}_{j}
  &\leftarrow
  f^{\mathbf{y}_c}_{j}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{j}\right) \\
\dfrac{\mathrm{d}\mathbf{i}^{\mathrm{sh}}_{j}}{\mathrm{d}t}
  &\leftarrow
  \dfrac{\mathrm{d}}{\mathrm{d}t}
  \left[
    f^{\mathbf{y}_c}_{j}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{j}\right)
  \right] \\
\mathbf{i}^{\mathrm{inc}}_{j}
  &\leftarrow \mathbf{i}^{\mathrm{inc,init}}_{j} \\
\dfrac{\mathrm{d}\mathbf{i}^{\mathrm{inc}}_{j}}{\mathrm{d}t}
  &\leftarrow
  \left(\dfrac{\mathrm{d}\mathbf{i}^{\mathrm{inc}}_{j}}{\mathrm{d}t}\right)^{\mathrm{init}},
  \quad j \in \{1,2\}
\end{aligned}
```

### Output Initialization

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
`i_sh_1` | [A] | Shunt current at terminal 1 | $\mathbf{i}^{\mathrm{sh}}_{1} \in \mathbb{R}^K$
`i_sh_2` | [A] | Shunt current at terminal 2 | $\mathbf{i}^{\mathrm{sh}}_{2} \in \mathbb{R}^K$
`i_inc_1` | [A] | Incident current at terminal 1 | $\mathbf{i}^{\mathrm{inc}}_{1} \in \mathbb{R}^K$
`i_inc_2` | [A] | Incident current at terminal 2 | $\mathbf{i}^{\mathrm{inc}}_{2} \in \mathbb{R}^K$
