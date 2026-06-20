# LineLumped Model

`LineLumped` represents a lumped-parameter EMT transmission line using a
$\pi$-section equivalent over length $\Delta x$. Series current $\mathbf{i}$
is directed from terminal 1 to terminal 2.

Notes:
- $\mathbf{Z}'$ and $\mathbf{Y}'$ may be supplied as per-unit-length
  `VectorFit` models or derived from constant `Rp`, `Lp`, `Gp`, and `Cp`
  matrices.

## Block Diagram

<div align="center">
   <img align="center" src="../../../../../../docs/Figures/EMT/LineLumped/diagram.png">

  Figure 1: LineLumped model
</div>

## Model Parameters

Symbol           | Units        | JSON | Description                        | Note
---------------- | ------------ | ---- | ---------------------------------- | ----
$\mathbf{Z}'$    | [$\Omega$/m] | `Zp` | Series impedance per unit length   | $\mathbf{Z}'\in\mathbb{C}^{N \times N}$
$\mathbf{Y}'$    | [S/m]        | `Yp` | Shunt admittance per unit length   | $\mathbf{Y}'\in\mathbb{C}^{N \times N}$
$\Delta x$       | [m]          | `dx` | Line segment length                | $\mathbb{R}$

### Parameter Validation

None.

### Model Derived Parameters

The per-unit-length impedance and admittance models are:

```math
\begin{aligned}
  \mathbf{Z}'(s) &= \mathbf{R}' + s\mathbf{L}' + \mathbf{Z}'_{\mathrm{rat}}(s) \\
  \mathbf{Y}'(s) &= \mathbf{G}' + s\mathbf{C}' + \mathbf{Y}'_{\mathrm{rat}}(s)
\end{aligned}
```

The segment impedance and shunt admittance operators used by the equations are:

```math
\begin{aligned}
  \mathbf{Z} &= \Delta x \mathbf{Z}' \\
  \mathbf{Y} &= \dfrac{\Delta x \mathbf{Y}'}{2}
\end{aligned}
```

### Submodels

Submodel | Inputs | Parameters | Outputs
-------- | ------ | ---------- | -------
[`VectorFit`](../../../Operators/Rational/VectorFit/README.md) series impedance $\mathbf{z}$ | $\mathbf{i}\in\mathbb{R}^N$ | `Zp`, `dx` | $\mathbf{z}*\mathbf{i}$
[`VectorFit`](../../../Operators/Rational/VectorFit/README.md) shunt admittance $\mathbf{y}_1$ | $\mathbf{v}_1\in\mathbb{R}^N$ | `Yp`, `dx` | $\mathbf{y}_1*\mathbf{v}_1$
[`VectorFit`](../../../Operators/Rational/VectorFit/README.md) shunt admittance $\mathbf{y}_2$ | $\mathbf{v}_2\in\mathbb{R}^N$ | `Yp`, `dx` | $\mathbf{y}_2*\mathbf{v}_2$

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
  - \mathbf{y}_1*\mathbf{v}_1
  - \mathbf{i} \\
\mathbf{i}^{\mathrm{inj}}_{2} &=
  - \mathbf{y}_2*\mathbf{v}_2
  + \mathbf{i}
\end{aligned}
```

## Initialization

Since this component is a member of the network, the power flow solution must initialize this model.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`i` | [A] | Series current | $\mathbf{i} \in \mathbb{R}^N$
