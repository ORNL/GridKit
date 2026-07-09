# LineLumped Model

`LineLumped` represents an $N$-phase, $K$-conductor lumped EMT line over length
$\Delta x$. Series current $\mathbf{i}$ is directed from terminal 1 to terminal
2.

Notes:
- Initial implementation is manual $N=3$ because submodel mechanics are not yet
  designed.

## Block Diagram

![](../../../../../../docs/Figures/EMT/LineLumped/diagram.png)

Figure 1: LineLumped model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer
$K$ | [-] | `K` | Number of conductors | Required, positive integer
$\mathbf{P}_\phi$ | [-] | `conductors` | Permutation matrix mapping each conductor to its phase | $\mathbf{P}_\phi \in \mathbb{R}^{N \times K}$
$\Delta x$       | [m]          | `dx` | Line segment length                | $\mathbb{R}$

### Parameter Validation

```math
\begin{aligned}
N &> 0 \\
K &> 0 \\
\mathbf{P}_\phi &\in \mathbb{R}^{N \times K} \\
\Delta x &> 0
\end{aligned}
```

### Derived Parameters

None.

## Submodels

Submodel | Type | Order | Inputs | Outputs
-------- | ---- | ----- | ------ | -------
Series impedance $f^{\mathbf{z}}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{z}}$ | $\mathbf{i}\in\mathbb{R}^K$ | $f^{\mathbf{z}}(\mathbf{i})\in\mathbb{R}^K$
Shunt admittance $f^{\mathbf{y}}$ | [`VectorFit`](../../../Operators/Rational/VectorFit/README.md) | $Q_{\mathbf{y}}$ | $\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_1,\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_2\in\mathbb{R}^K$ | $f^{\mathbf{y}}(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_1),f^{\mathbf{y}}(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_2)\in\mathbb{R}^K$

### Submodel Validation

```math
\begin{aligned}
f^{\mathbf{z}} &: \mathbb{R}^K \rightarrow \mathbb{R}^K \\
f^{\mathbf{y}} &: \mathbb{R}^K \rightarrow \mathbb{R}^K
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$   | [A]    | Series current, directed terminal `1` to terminal `2` | $\mathbf{i} \in \mathbb{R}^K$

#### Algebraic

None.

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}_1$   | [V]    | Terminal `1` voltage owned by EMT bus | $\mathbf{v}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$   | [V]    | Terminal `2` voltage owned by EMT bus | $\mathbf{v}_2 \in \mathbb{R}^N$

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
  \mathbf{P}_\phi^{\mathsf T}\left(\mathbf{v}_2-\mathbf{v}_1\right)
```

### Algebraic Equations

None.

### Wiring

```math
\begin{aligned}
\mathbf{i}^{\mathrm{inj}}_{1} &=
  \mathbf{P}_\phi
  \left(
    - \dfrac{\Delta x}{2}f^{\mathbf{y}}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_1\right)
    - \mathbf{i}
  \right) \\
\mathbf{i}^{\mathrm{inj}}_{2} &=
  \mathbf{P}_\phi
  \left(
    - \dfrac{\Delta x}{2}f^{\mathbf{y}}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_2\right)
    + \mathbf{i}
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

Initialization sets the line-current start by enforcing the differential
residual.

```math
\begin{aligned}
\mathbf{i}
  &\leftarrow \text{line-current start} \\
  \Delta x f^{\mathbf{z}}(\mathbf{i}) +
  \mathbf{P}_\phi^{\mathsf T}
  \left(\mathbf{v}_{2}-\mathbf{v}_{1}\right)
  &\leftarrow \mathbf{0}
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
\mathbf{i}^{\mathrm{inj}}_{1}
  &\leftarrow
  \mathbf{P}_\phi
  \left(
    - \dfrac{\Delta x}{2}
      f^{\mathbf{y}}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{1}\right)
    - \mathbf{i}
  \right) \\
\mathbf{i}^{\mathrm{inj}}_{2}
  &\leftarrow
  \mathbf{P}_\phi
  \left(
    - \dfrac{\Delta x}{2}
      f^{\mathbf{y}}\left(\mathbf{P}_\phi^{\mathsf T}\mathbf{v}_{2}\right)
    + \mathbf{i}
  \right)
\end{aligned}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Series current | $\mathbf{i} \in \mathbb{R}^K$
