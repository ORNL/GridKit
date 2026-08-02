# LineDistributed Model

`LineDistributed` represents an $N$-phase, $K$-conductor distributed EMT line.

## Block Diagram

![LineDistributed model block diagram](../../../../../../docs/Figures/EMT/LineDistributed/diagram.png)

Figure 1: LineDistributed model

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

### Parameter Validation

```math
\begin{aligned}
N &\in \mathbb{Z}_{>0} \\
K &\in \mathbb{Z}_{>0} \\
\mathbf{c} &\in \mathcal{N}^K \\
\{c_k \mid k \in \mathcal{K}\} &= \mathcal{N}
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
$\mathbf{i}_1$ | `i1` | Output | [A] | Current injection at terminal 1 | $\mathbf{i}_1 \in \mathbb{R}^N$
$\mathbf{i}_2$ | `i2` | Output | [A] | Current injection at terminal 2 | $\mathbf{i}_2 \in \mathbb{R}^N$

## Submodels

Symbol | Description | Type | Order | JSON | Inputs | Outputs
------ | ----------- | ---- | ----- | ---- | ------ | -------
$\mathbf{y}_1^\mathrm{c}$ | Characteristic admittance at terminal 1 | [VectorFit](../../../Operators/Rational/VectorFit/README.md) | $KQ_{\mathbf{y}^\mathrm{c}}$ | `Yc` | $\mathbb{R}^K$ | $\mathbb{R}^K$
$\mathbf{y}_2^\mathrm{c}$ | Characteristic admittance at terminal 2 | [VectorFit](../../../Operators/Rational/VectorFit/README.md) | $KQ_{\mathbf{y}^\mathrm{c}}$ | `Yc` | $\mathbb{R}^K$ | $\mathbb{R}^K$
$\mathbf{h}_{21}$ | Propagation from terminal 2 to terminal 1 | [Propagation](../../../Operators/Shift/Propagation/README.md) | Composite | `H` | $\mathbb{R}^K$ | $\mathbb{R}^K$
$\mathbf{h}_{12}$ | Propagation from terminal 1 to terminal 2 | [Propagation](../../../Operators/Shift/Propagation/README.md) | Composite | `H` | $\mathbb{R}^K$ | $\mathbb{R}^K$

`Yc` and `H` each provide one coefficient set. The model assumes a reciprocal
uniform line, so

```math
\mathbf{Y}_1^{\mathrm{c}}(s)=\mathbf{Y}_2^{\mathrm{c}}(s),
\qquad
\mathbf{H}_{12}(s)=\mathbf{H}_{21}(s).
```

The two terminal-admittance instances maintain independent states, and the two
directional propagation instances maintain independent states and histories.

### Submodel Validation

The characteristic-admittance fits must be stable, proper, and positive real.
Together with the propagation fits, they must produce a passive line model.

### Submodel Wiring

```math
\begin{aligned}
\mathbf{i}_1^\mathrm{c} &\leftarrow
  \mathbf{y}_1^\mathrm{c}[\mathbf{P}_\phi^\mathsf T\mathbf{v}_1] \\
\mathbf{i}_2^\mathrm{c} &\leftarrow
  \mathbf{y}_2^\mathrm{c}[\mathbf{P}_\phi^\mathsf T\mathbf{v}_2] \\
\mathbf{i}_1^\mathrm{inc} &\leftarrow
  \mathbf{h}_{21}[\mathbf{i}_2^\mathrm{ref}] \\
\mathbf{i}_2^\mathrm{inc} &\leftarrow
  \mathbf{h}_{12}[\mathbf{i}_1^\mathrm{ref}]
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}_1^\mathrm{ref}$ | [A] | Reflected current at terminal 1 | $\mathbf{i}_1^\mathrm{ref} \in \mathbb{R}^K$
$\mathbf{i}_2^\mathrm{ref}$ | [A] | Reflected current at terminal 2 | $\mathbf{i}_2^\mathrm{ref} \in \mathbb{R}^K$

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}_1$ | [V] | Terminal 1 voltage owned by EMT bus | $\mathbf{v}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$ | [V] | Terminal 2 voltage owned by EMT bus | $\mathbf{v}_2 \in \mathbb{R}^N$

#### Algebraic

None.

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

The residuals use the signals defined under Submodel Wiring.

```math
\begin{aligned}
0 &= -\mathbf{i}_1^\mathrm{ref}
  + 2\mathbf{i}_1^\mathrm{c}
  - \mathbf{i}_1^\mathrm{inc} \\
0 &= -\mathbf{i}_2^\mathrm{ref}
  + 2\mathbf{i}_2^\mathrm{c}
  - \mathbf{i}_2^\mathrm{inc}
\end{aligned}
```

### External Equations

```math
\begin{aligned}
\mathbf{i}_1 &\leftarrow
  \mathbf{P}_\phi(\mathbf{i}_1^\mathrm{inc}-\mathbf{i}_1^\mathrm{c}) \\
\mathbf{i}_2 &\leftarrow
  \mathbf{P}_\phi(\mathbf{i}_2^\mathrm{inc}-\mathbf{i}_2^\mathrm{c})
\end{aligned}
```

## Initialization

The characteristic-admittance and propagation submodels initialize according
to their own specifications. There is no additional line-level initialization
procedure.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i_c1` | [A] | Characteristic-admittance current at terminal 1 | $\mathbf{i}_1^\mathrm{c} \in \mathbb{R}^K$
`i_c2` | [A] | Characteristic-admittance current at terminal 2 | $\mathbf{i}_2^\mathrm{c} \in \mathbb{R}^K$
`i_inc1` | [A] | Incident current at terminal 1 | $\mathbf{i}_1^\mathrm{inc} \in \mathbb{R}^K$
`i_inc2` | [A] | Incident current at terminal 2 | $\mathbf{i}_2^\mathrm{inc} \in \mathbb{R}^K$

## Development

The initial three-phase formulation takes $N=K=3$ and
$\mathbf{P}_\phi=\mathbf{I}_3$.

```math
\begin{aligned}
0 &= -\mathbf{i}_1^\mathrm{ref}
  + 2\mathbf{i}_1^\mathrm{c}
  - \mathbf{i}_1^\mathrm{inc} \\
0 &= -\mathbf{i}_2^\mathrm{ref}
  + 2\mathbf{i}_2^\mathrm{c}
  - \mathbf{i}_2^\mathrm{inc}
\end{aligned}
```

For the initial three-phase formulation, the external equations reduce to

```math
\begin{aligned}
\mathbf{i}_1
  &\leftarrow \mathbf{i}_1^\mathrm{inc}-\mathbf{i}_1^\mathrm{c} \\
\mathbf{i}_2
  &\leftarrow \mathbf{i}_2^\mathrm{inc}-\mathbf{i}_2^\mathrm{c}
\end{aligned}
```

and the submodel wiring reduces to

```math
\begin{aligned}
\mathbf{i}_1^\mathrm{c}
  &\leftarrow \mathbf{y}_1^\mathrm{c}[\mathbf{v}_1] \\
\mathbf{i}_2^\mathrm{c}
  &\leftarrow \mathbf{y}_2^\mathrm{c}[\mathbf{v}_2] \\
\mathbf{i}_1^\mathrm{inc}
  &\leftarrow \mathbf{h}_{21}[\mathbf{i}_2^\mathrm{ref}] \\
\mathbf{i}_2^\mathrm{inc}
  &\leftarrow \mathbf{h}_{12}[\mathbf{i}_1^\mathrm{ref}]
\end{aligned}
```
