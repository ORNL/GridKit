# LineLumped Model

`LineLumped` represents an $N$-phase, $K$-conductor lumped EMT line over length
$\Delta x$. Series current $\mathbf{i}_{12}$ is directed from terminal 1 to
terminal 2. Each terminal Bus owns its shunt current and
admittance states.

## Block Diagram

![LineLumped model block diagram](../../../../../../docs/Figures/EMT/LineLumped/diagram.png)

Figure 1: LineLumped model

The conductor-to-phase mappings are shown in the equations and omitted from
the diagram for clarity.

## Model Parameters

Define the phase- and conductor-index sets

```math
\mathcal{N} = \{1,\ldots,N\},
\qquad
\mathcal{K} = \{1,\ldots,K\}.
```

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$N$ | [-] | `N` | Number of phases | Required, positive integer
$K$ | [-] | `K` | Number of conductors | Required, positive integer
$\mathbf{c}$ | [-] | `conductors` | Conductor phase-index list | $\mathbf{c} \in \mathcal{N}^K$
$\Delta x$ | [m] | `dx` | Line segment length | Required, positive

### Parameter Validation

```math
\begin{aligned}
N &\in \mathbb{Z}_{>0} \\
K &\in \mathbb{Z}_{>0} \\
\mathbf{c} &\in \mathcal{N}^K \\
\{c_k \mid k \in \mathcal{K}\} &= \mathcal{N} \\
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
\quad n \in \mathcal{N},\quad k \in \mathcal{K}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}_1$ | `v1` | Input | [V] | Terminal 1 bus voltage | $\mathbf{v}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$ | `v2` | Input | [V] | Terminal 2 bus voltage | $\mathbf{v}_2 \in \mathbb{R}^N$
$\mathbf{i}_{12}$ | `i12` | Output | [A] | Series current from terminal 1 to 2 | $\mathbb{R}^K$
$\mathbf{i}_{21}$ | `i21` | Output | [A] | Series current from terminal 2 to 1 | $\mathbf{i}_{21}=-\mathbf{i}_{12}$

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{z}$ | Per-unit-length series impedance | [VectorFit](../../../Operators/Rational/VectorFit/README.md) | $KQ_{\mathbf{z}}$ | `Zp` | $\mathbb{R}^K$ | $\mathbb{R}^K$
$\mathbf{y}_1$ | Per-unit-length shunt admittance at terminal 1 | [VectorFit](../../../Operators/Rational/VectorFit/README.md) | $KQ_{\mathbf{y}}$ | `Yp` | $\mathbb{R}^K$ | $\mathbb{R}^K$
$\mathbf{y}_2$ | Per-unit-length shunt admittance at terminal 2 | [VectorFit](../../../Operators/Rational/VectorFit/README.md) | $KQ_{\mathbf{y}}$ | `Yp` | $\mathbb{R}^K$ | $\mathbb{R}^K$

`Yp` provides one coefficient set; the two terminal instances maintain
independent states in their buses, with scale $\Delta x/2$.

### Submodel Validation

```math
\mathrm{rank}(\mathbf{E}^{\mathbf{z}}) = K
```

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}_{12}$ | [A] | Series current from terminal 1 to terminal 2 | $\mathbf{i}_{12} \in \mathbb{R}^K$

#### Algebraic

None.

### External Variables

#### Differential

When a connected equation depends on the bus-voltage derivative:

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}_1$ | [V] | Terminal 1 voltage owned by EMT bus | $\mathbf{v}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$ | [V] | Terminal 2 voltage owned by EMT bus | $\mathbf{v}_2 \in \mathbb{R}^N$

#### Algebraic

Otherwise, the bus-voltage variables above are algebraic.

## Model Equations

### Internal Equations

#### Differential

```math
0 = \Delta x\,\mathbf{z}[\mathbf{i}_{12}]
  + \mathbf{P}_\phi^\mathsf{T}(\mathbf{v}_2-\mathbf{v}_1)
```

#### Algebraic

None.

### External Equations

None. The buses register $-\mathbf{i}_{12}$ and $-\mathbf{i}_{21}$ as
series-current injections at terminals 1 and 2. Each bus owns its shunt branch.

## Initialization

None beyond the EMT initialization contract.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i12` | [A] | Series current from terminal 1 to terminal 2 | $\mathbf{i}_{12} \in \mathbb{R}^K$

## Development

The initial three-phase formulation is a subset of the generalized formulation
above.

### Derived Parameters

```math
\begin{aligned}
\mathbf{R} &= \Delta x\,\mathbf{R}' \\
\mathbf{L} &= \Delta x\,\mathbf{L}' \\
\mathbf{G} &= \Delta x\,\mathbf{G}' \\
\mathbf{C} &= \Delta x\,\mathbf{C}'
\end{aligned}
```

### Differential Equations

```math
0 = \mathbf{R}\mathbf{i}_{12}
  + \mathbf{L}\dfrac{\mathrm{d}\mathbf{i}_{12}}{\mathrm{d}t}
  + \mathbf{v}_2-\mathbf{v}_1
```
