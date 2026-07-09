# LineLumped Model

`LineLumped` represents an $N$-phase, $K$-conductor lumped EMT line over length
$\Delta x$. Series current $\mathbf{i}$ is directed from terminal 1 to terminal
2.

> [!NOTE]
> A template parameter will select whether $\mathbf{i}$ is a differential or
> algebraic vector. This page documents the nondegenerate differential
> formulation.

## Block Diagram

![LineLumped model block diagram](../../../../../../docs/Figures/EMT/LineLumped/diagram.png)

Figure 1: LineLumped model

The conductor-to-phase mappings are shown in the equations and omitted from
the diagram for clarity.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer
$K$ | [-] | `K` | Number of conductors | Required, positive integer
$\mathbf{c}$ | [-] | `conductors` | Conductor phase-index list | $\mathbf{c} \in \{1,\ldots,N\}^{K}$
$\Delta x$ | [m] | `dx` | Line segment length | Required, positive

### Parameter Validation

```math
\begin{aligned}
N &> 0 \\
K &> 0 \\
\mathbf{c} &\in \{1,\ldots,N\}^{K} \\
\{c_k \mid k \in \{1,\ldots,K\}\} &= \{1,\ldots,N\} \\
\Delta x &> 0
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
Per-unit-length series impedance $f^{\mathbf{z}}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{z}}$ | $\mathbb{R}^K$ [A] | `Zp` | $\mathbb{R}^K$ [V/m]
Per-unit-length terminal-1 shunt admittance $f^{\mathbf{y}}_{1}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}}$ | $\mathbb{R}^K$ [V] | `Yp` | $\mathbb{R}^K$ [A/m]
Per-unit-length terminal-2 shunt admittance $f^{\mathbf{y}}_{2}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}}$ | $\mathbb{R}^K$ [V] | `Yp` | $\mathbb{R}^K$ [A/m]

For the pole-free RLGC specialization, the `VectorFit` coefficients map to
the per-unit-length matrices as

```math
\begin{aligned}
\mathbf{D}^{\mathbf{z}} &= \mathbf{R}', &
\mathbf{E}^{\mathbf{z}} &= \mathbf{L}' \\
\mathbf{D}^{\mathbf{y}} &= \mathbf{G}', &
\mathbf{E}^{\mathbf{y}} &= \mathbf{C}'
\end{aligned}
```

### Submodel Validation

```math
\mathrm{rank}\left(\mathbf{E}^{\mathbf{z}}\right) = K
```

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$ | [A] | Series current, directed from terminal 1 to terminal 2 | $\mathbf{i} \in \mathbb{R}^K$

#### Algebraic

None.

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

```math
0 = \Delta x f^{\mathbf{z}}(\mathbf{i}) +
  \mathbf{P}_\phi^{\mathsf T}\left(\mathbf{v}_{2} - \mathbf{v}_{1}\right)
```

### Algebraic Equations

None.

### Wiring

```math
\begin{aligned}
\mathbf{i}^{\mathrm{inj}}_{1} &=
  \mathbf{P}_\phi
  \left(
    - \dfrac{\Delta x}{2}f^{\mathbf{y}}_{1}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{1}\right)
    - \mathbf{i}
  \right) \\
\mathbf{i}^{\mathrm{inj}}_{2} &=
  \mathbf{P}_\phi
  \left(
    - \dfrac{\Delta x}{2}f^{\mathbf{y}}_{2}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{2}\right)
    + \mathbf{i}
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

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}_0$ | [A] | Initial series-current vector | Required, $\mathbf{i}_0 \in \mathbb{R}^K$

```math
\mathbf{i} \leftarrow \mathbf{i}_0
```

### Output Initialization

```math
\begin{aligned}
\mathbf{i}^{\mathrm{inj}}_{1}
  &\leftarrow
  \mathbf{P}_\phi
  \left(
    - \dfrac{\Delta x}{2}
      f^{\mathbf{y}}_{1}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{1}\right)
    - \mathbf{i}
  \right) \\
\mathbf{i}^{\mathrm{inj}}_{2}
  &\leftarrow
  \mathbf{P}_\phi
  \left(
    - \dfrac{\Delta x}{2}
      f^{\mathbf{y}}_{2}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{2}\right)
    + \mathbf{i}
  \right)
\end{aligned}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Series current | $\mathbf{i} \in \mathbb{R}^K$
