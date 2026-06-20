# LineLumped Model

`LineLumped` represents a lumped-parameter EMT transmission line using a
$\pi$-section equivalent over length $\Delta x$. Series current $\mathbf{i}$
is directed from terminal 1 to terminal 2.

Notes:
- Alternatively, $\mathbf{Z}'$ and $\mathbf{Y}'$ may be given directly as
  `VectorFit` models. Otherwise, they are derived from constant `R`, `L`,
  `G`, and `C` matrices.

## Block Diagram

<div align="center">
   <img align="center" src="../../../../../../docs/Figures/EMT/LineLumped/diagram.png">

  Figure 1: LineLumped model
</div>

## Model Parameters

Symbol           | Units          | JSON | Description                              | Note
---------------- | -------------- | ---- | ---------------------------------------- | ----
$\mathbf{R}'$    | [$\Omega$/m]   | `R`  | Series resistance matrix per unit length | $\mathbb{R}^{N \times N}$
$\mathbf{L}'$    | [H/m]          | `L`  | Series inductance matrix per unit length | $\mathbb{R}^{N \times N}$
$\mathbf{G}'$    | [S/m]          | `G`  | Shunt conductance matrix per unit length | $\mathbb{R}^{N \times N}$
$\mathbf{C}'$    | [F/m]          | `C`  | Shunt capacitance matrix per unit length | $\mathbb{R}^{N \times N}$
$\Delta x$       | [m]            | `dx` | Line segment length                      | $\mathbb{R}$
$\mathbf{Z}'$    | [$\Omega$/m]   | `Z`  | Series impedance per unit length         | $\mathbb{C}^{N\times N}$, optional `VectorFit` model 
$\mathbf{Y}'$    | [S/m]          | `Y`  | Shunt admittance per unit length         | $\mathbb{C}^{N\times N}$, optional `VectorFit` model 

### Parameter Validation

None.

### Model Derived Parameters

For the constant `R`, `L`, `G`, and `C` alternative:

```math
\begin{aligned}
  \mathbf{R} &= \mathbf{R}'\Delta x   & \mathbf{G} &= \dfrac{\mathbf{G}'\Delta x}{2} \\
  \mathbf{L} &= \mathbf{L}'\Delta x   & \mathbf{C} &= \dfrac{\mathbf{C}'\Delta x}{2}
\end{aligned}
```

If `Z` and `Y` are not given, they are derived as:

```math
\begin{aligned}
  \mathbf{Z} &= \mathbf{R} + s\mathbf{L} \qquad \text{or} \qquad = \Delta x \mathbf{Z}'\\
  \mathbf{Y} &= \mathbf{G} + s\mathbf{C} \qquad \text{or} \qquad = \dfrac{\Delta x \mathbf{Y}'}{2}
\end{aligned}
```

and modeled via `VectorFit`.

## Model Variables

### Internal Variables

#### Differential

Symbol           | Units  | Description           | Note
-----------------|--------|-----------------------|---------------------------------
$\mathbf{i}$   | [A]    | Series current, directed terminal `1` to terminal `2` | $\mathbf{i} \in \mathbb{R}^N$

#### Algebraic

None.

### External Variables

#### Differential

Symbol           | Units  | Description              | Note
-----------------|--------|--------------------------|------------------
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
0 = \mathbf{z}*\mathbf{i}+ \mathbf{v}_2 - \mathbf{v}_1
```

### Algebraic Equations

None.

### Wiring

```math
\begin{aligned}
\mathbf{i}^{\mathrm{inj}}_{1} &=
  - \mathbf{y}*\mathbf{v}_1
  - \mathbf{i} \\
\mathbf{i}^{\mathrm{inj}}_{2} &=
  - \mathbf{y}*\mathbf{v}_2
  + \mathbf{i}
\end{aligned}
```

## Initialization

Since this component is a member of the network, the power flow solution must initialize this model.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Series current | $\mathbf{i} \in \mathbb{R}^N$
