# LineLumped Model

`LineLumped` represents an $N$-phase, $K$-conductor lumped EMT line over length
$\Delta x$. Series current $\mathbf{i}$ is directed from terminal 1 to terminal
2.

> [!NOTE]
> The initial end-to-end implementation will support three-phase systems only
> to establish a proof of concept. The formulation below remains $N$-phase.

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
$\mathbf{P}_\phi$ | [-] | `conductors` | Conductor-to-phase assignment matrix | $\mathbf{P}_\phi \in \{0,1\}^{N \times K}$
$\Delta x$ | [m] | `dx` | Line segment length | Required, positive

### Parameter Validation

```math
\begin{aligned}
N &> 0 \\
K &> 0 \\
\mathbf{P}_\phi &\in \{0,1\}^{N \times K} \\
\sum_{n = 1}^{N}P_{\phi,nk} &= 1,
  \quad k \in \{1,\ldots,K\} \\
\Delta x &> 0
\end{aligned}
```

### Derived Parameters

None.

## Submodels

The line instantiates one series-impedance model and two terminal
shunt-admittance models from the configured `VectorFit` specifications.

Submodel | Type | Order | Inputs | JSON | Outputs
-------- | ---- | ----- | ------ | ---- | -------
Per-unit-length series impedance $f^{\mathbf{z}}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{z}}$ | $\mathbf{i} \in \mathbb{R}^K$ | `Zp` | $f^{\mathbf{z}}(\mathbf{i}) \in \mathbb{R}^K$
Per-unit-length terminal-1 shunt admittance $f^{\mathbf{y}}_{1}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}}$ | $\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{1} \in \mathbb{R}^K$ | `Yp` | $f^{\mathbf{y}}_{1}(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{1}) \in \mathbb{R}^K$
Per-unit-length terminal-2 shunt admittance $f^{\mathbf{y}}_{2}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}}$ | $\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{2} \in \mathbb{R}^K$ | `Yp` | $f^{\mathbf{y}}_{2}(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{2}) \in \mathbb{R}^K$

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
\begin{aligned}
f^{\mathbf{z}} &: \mathbb{R}^K \rightarrow \mathbb{R}^K \\
f^{\mathbf{y}}_{j} &: \mathbb{R}^K \rightarrow \mathbb{R}^K,
  \quad j \in \{1,2\} \\
\mathbf{E}^{\mathbf{z}} &\in \mathbb{R}^{K \times K} \\
\mathrm{rank}\left(\mathbf{E}^{\mathbf{z}}\right) &= K
\end{aligned}
```

$\mathbf{E}^{\mathbf{z}}$ is the `E` matrix of the series-impedance
`VectorFit` submodel. Its full rank ensures that $\mathbf{i}$ is differential.
Static or singular fits require a corresponding algebraic-current formulation.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$ | [A] | Series current, directed terminal `1` to terminal `2` | $\mathbf{i} \in \mathbb{R}^K$

#### Algebraic

None.

### External Variables

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
  &\leftarrow \text{provisional terminal-bus voltage trajectories}
\end{aligned}
```

### Internal Initialization

The series current and its provisional derivative define the
series-impedance input trajectory. The initialized terminal-bus trajectories
define the two shunt-admittance input trajectories:

```math
\begin{aligned}
\mathbf{i}
  &\leftarrow \text{line-current start} \\
\dfrac{\mathrm{d}\mathbf{i}}{\mathrm{d}t}
  &\leftarrow \mathbf{0} \\
\left(\mathbf{i},\dfrac{\mathrm{d}\mathbf{i}}{\mathrm{d}t}\right)
  &\leftarrow \text{series-impedance input trajectory} \\
\left(
  \mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{j},
  \mathbf{P}_\phi^{\mathsf T}
    \dfrac{\mathrm{d}\mathbf{v}_{j}}{\mathrm{d}t}
\right)
  &\leftarrow \text{terminal-$j$ shunt-admittance input trajectory},
    \quad j \in \{1,2\}
\end{aligned}
```

The child states and outputs are initialized from those trajectories. The
assembled line and network residuals then replace the provisional derivative
and determine consistent operator outputs while satisfying

```math
\Delta x f^{\mathbf{z}}(\mathbf{i})
  + \mathbf{P}_\phi^{\mathsf T}
    \left(\mathbf{v}_{2} - \mathbf{v}_{1}\right)
  \leftarrow \mathbf{0}
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
