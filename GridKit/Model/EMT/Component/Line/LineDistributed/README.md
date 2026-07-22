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

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}_1$ | `v1` | Input | [V] | Terminal 1 bus voltage | $\mathbf{v}_1 \in \mathbb{R}^N$
$\mathbf{v}_2$ | `v2` | Input | [V] | Terminal 2 bus voltage | $\mathbf{v}_2 \in \mathbb{R}^N$
$\mathbf{i}_1$ | `i1` | Output | [A] | Current injection at terminal 1 | $\mathbf{i}_1 \in \mathbb{R}^N$
$\mathbf{i}_2$ | `i2` | Output | [A] | Current injection at terminal 2 | $\mathbf{i}_2 \in \mathbb{R}^N$

## Model Equations

### Differential Equations

None.

### Algebraic Equations

The residuals use the wiring signals defined below.

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

### Wiring

```math
\begin{aligned}
\mathbf{i}_1^\mathrm{c} &\leftarrow
  \mathbf{y}_1^\mathrm{c}[\mathbf{P}_\phi^\mathsf T\mathbf{v}_1] \\
\mathbf{i}_2^\mathrm{c} &\leftarrow
  \mathbf{y}_2^\mathrm{c}[\mathbf{P}_\phi^\mathsf T\mathbf{v}_2] \\
\mathbf{i}_1^\mathrm{inc} &\leftarrow
  \mathbf{h}_{21}[\mathbf{i}_2^\mathrm{ref}] \\
\mathbf{i}_2^\mathrm{inc} &\leftarrow
  \mathbf{h}_{12}[\mathbf{i}_1^\mathrm{ref}] \\
\mathbf{i}_1 &\leftarrow
  \mathbf{P}_\phi(\mathbf{i}_1^\mathrm{inc}-\mathbf{i}_1^\mathrm{c}) \\
\mathbf{i}_2 &\leftarrow
  \mathbf{P}_\phi(\mathbf{i}_2^\mathrm{inc}-\mathbf{i}_2^\mathrm{c})
\end{aligned}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
\mathbf{V}_r
  &\leftarrow \text{solved terminal RMS voltage phasor} \\
\mathbf{v}_r
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\mathbf{V}_r) \\
\dfrac{\mathrm{d}\mathbf{v}_r}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(s_0\mathbf{V}_r),
  \quad r \in \{1,2\}.
\end{aligned}
```

### Internal Initialization

The internal current phasors satisfy

```math
\begin{aligned}
\mathbf{I}_r^\mathrm{c}
  &= \mathbf{Y}_r^\mathrm{c}(s_0)
     \mathbf{P}_\phi^\mathsf T\mathbf{V}_r,
  \quad r \in \{1,2\} \\
\begin{bmatrix}
\mathbf{I}_K & \mathbf{H}_{21}(s_0) \\
\mathbf{H}_{12}(s_0) & \mathbf{I}_K
\end{bmatrix}
\begin{bmatrix}
\mathbf{I}_1^\mathrm{ref} \\
\mathbf{I}_2^\mathrm{ref}
\end{bmatrix}
  &= 2
\begin{bmatrix}
\mathbf{I}_1^\mathrm{c} \\
\mathbf{I}_2^\mathrm{c}
\end{bmatrix} \\
\mathbf{I}_1^\mathrm{inc}
  &= \mathbf{H}_{21}(s_0)\mathbf{I}_2^\mathrm{ref} \\
\mathbf{I}_2^\mathrm{inc}
  &= \mathbf{H}_{12}(s_0)\mathbf{I}_1^\mathrm{ref}.
\end{aligned}
```

The block matrix must be nonsingular. The characteristic-admittance submodels
initialize from $\mathbf{P}_\phi^\mathsf T\mathbf{V}_r$; $\mathbf{h}_{21}$ and
$\mathbf{h}_{12}$ initialize from $\mathbf{I}_2^\mathrm{ref}$ and
$\mathbf{I}_1^\mathrm{ref}$, respectively. At $t=0$, the reflected currents
initialize as

```math
\begin{aligned}
\mathbf{i}_r^\mathrm{ref}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\mathbf{I}_r^\mathrm{ref}) \\
\dfrac{\mathrm{d}\mathbf{i}_r^\mathrm{ref}}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(s_0\mathbf{I}_r^\mathrm{ref}),
  \quad r \in \{1,2\}.
\end{aligned}
```

### Output Initialization

```math
\begin{aligned}
\mathbf{I}_r
  &= \mathbf{P}_\phi(\mathbf{I}_r^\mathrm{inc}-\mathbf{I}_r^\mathrm{c}) \\
\mathbf{i}_r
  &\leftarrow \sqrt{2}\,\mathrm{Re}(\mathbf{I}_r) \\
\dfrac{\mathrm{d}\mathbf{i}_r}{\mathrm{d}t}
  &\leftarrow \sqrt{2}\,\mathrm{Re}(s_0\mathbf{I}_r),
  \quad r \in \{1,2\}.
\end{aligned}
```

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

### Algebraic Equations

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

### Wiring

```math
\begin{aligned}
\mathbf{i}_1^\mathrm{c}
  &\leftarrow \mathbf{y}_1^\mathrm{c}[\mathbf{v}_1] \\
\mathbf{i}_2^\mathrm{c}
  &\leftarrow \mathbf{y}_2^\mathrm{c}[\mathbf{v}_2] \\
\mathbf{i}_1^\mathrm{inc}
  &\leftarrow \mathbf{h}_{21}[\mathbf{i}_2^\mathrm{ref}] \\
\mathbf{i}_2^\mathrm{inc}
  &\leftarrow \mathbf{h}_{12}[\mathbf{i}_1^\mathrm{ref}] \\
\mathbf{i}_1
  &\leftarrow \mathbf{i}_1^\mathrm{inc}-\mathbf{i}_1^\mathrm{c} \\
\mathbf{i}_2
  &\leftarrow \mathbf{i}_2^\mathrm{inc}-\mathbf{i}_2^\mathrm{c}
\end{aligned}
```
