# LineLumped Model

`LineLumped` represents a lumped-parameter EMT transmission line using a
$\pi$-section equivalent over length $\Delta x$. Series current $\mathbf{i}$
is directed from terminal 1 to terminal 2.

Notes:
- Constant `R`, `L`, `G`, and `C` parameters remain the supported interface.
  Future support for per-unit-length $\mathbf{Z}'$ and $\mathbf{Y}'$ as
  convolutional models should preserve this case: using only the $\mathbf{D}$
  and $\mathbf{E}$ terms gives
  $\mathbf{Z}'=\mathbf{R}'+s\mathbf{L}'$ and
  $\mathbf{Y}'=\mathbf{G}'+s\mathbf{C}'$.

<div align="center">
   <img align="center" src="../../../../../../docs/Figures/EMT/LineLumped/diagram.svg">

  Figure 1: LineLumped model
</div>

## Model Parameters

For phase count $N$:

Symbol           | Units          | Description                              | Note
-----------------|----------------|------------------------------------------|---------------------------------
$\mathbf{R}'$    | [$\Omega$/m]   | Series resistance matrix per unit length | $\mathbb{R}^{N \times N}$
$\mathbf{L}'$    | [H/m]          | Series inductance matrix per unit length | $\mathbb{R}^{N \times N}$
$\mathbf{G}'$    | [S/m]          | Shunt conductance matrix per unit length | $\mathbb{R}^{N \times N}$
$\mathbf{C}'$    | [F/m]          | Shunt capacitance matrix per unit length | $\mathbb{R}^{N \times N}$
$\mathbf{Z}'$    | [$\Omega$/m]   | Series impedance per unit length         | $\mathbb{C}^{N \times N}$; future compatibility: supports passing a convolutional model
$\mathbf{Y}'$    | [S/m]          | Shunt admittance per unit length         | $\mathbb{C}^{N \times N}$; future compatibility: supports passing a convolutional model
$\Delta x$       | [m]            | Line segment length                      | $\mathbb{R}$

## Model Derived Parameters

```math
\begin{aligned}
  \mathbf{R} &= \mathbf{R}'\Delta x   & \mathbf{G} &= \dfrac{\mathbf{G}'\Delta x}{2} \\
  \mathbf{L} &= \mathbf{L}'\Delta x   & \mathbf{C} &= \dfrac{\mathbf{C}'\Delta x}{2}
\end{aligned}
```

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

For constant series parameters:

```math
0 = \mathbf{R}\,\mathbf{i} + \mathbf{L}\dot{\mathbf{i}} + \mathbf{v}_2 - \mathbf{v}_1
```

### Algebraic Equations

None.

### Port Equations

```math
\begin{aligned}
\mathbf{i}^{\mathrm{inj}}_{1} &=
  - \mathbf{G}\,\mathbf{v}_1
  - \mathbf{C}\,\dot{\mathbf{v}}_1
  - \mathbf{i} \\
\mathbf{i}^{\mathrm{inj}}_{2} &=
  - \mathbf{G}\,\mathbf{v}_2
  - \mathbf{C}\,\dot{\mathbf{v}}_2
  + \mathbf{i}
\end{aligned}
```

## Initialization

Given $\mathbf{i}_0$, the initial derivative follows from the series residual:

```math
\dot{\mathbf{i}}_0 = \mathbf{L}^{-1}\left(\mathbf{v}_{1,0} - \mathbf{v}_{2,0} - \mathbf{R}\,\mathbf{i}_0\right)
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Series current | $\mathbf{i} \in \mathbb{R}^N$
