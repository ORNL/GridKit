# Clarke Model

`Clarke` transforms three-phase quantities from $abc$ to stationary
$\alpha\beta0$ coordinates with power-invariant normalization.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\mathrm{inverse}$ | [-] | `inverse` | Inverse transformation | Default `false`

### Parameter Validation

The `inverse` parameter must be Boolean and is fixed during simulation.

### Derived Parameters

```math
\mathbf{C} =
\sqrt{\dfrac{2}{3}}
\begin{bmatrix}
1 & -\dfrac{1}{2} & -\dfrac{1}{2} \\
0 & \dfrac{\sqrt{3}}{2} & -\dfrac{\sqrt{3}}{2} \\
\dfrac{1}{\sqrt{2}} & \dfrac{1}{\sqrt{2}} & \dfrac{1}{\sqrt{2}}
\end{bmatrix}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | $[u]$ | Input vector port | $\mathbf{u} \in \mathbb{R}^3$
$\mathbf{y}$ | `out` | Output | $[u]$ | Algebraic output vector | $\mathbf{y} \in \mathbb{R}^3$

Forward input and output orders are $(a,b,c)$ and $(\alpha,\beta,0)$.
The inverse transformation exchanges these orders.

All inputs must be connected and finite. In case JSON, `input` and `out` map
three signal IDs to `u1`, `u2`, `u3` and `y1`, `y2`, `y3`, respectively.

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{y}$ | $[u]$ | Transformed output | $\mathbf{y} \in \mathbb{R}^3$

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | $[u]$ | Input vector | Differential-input configuration, $\mathbf{u} \in \mathbb{R}^3$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | $[u]$ | Input vector | Algebraic-input configuration, $\mathbf{u} \in \mathbb{R}^3$

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

```math
0 = -\mathbf{y} +
\begin{cases}
\mathbf{C}\mathbf{u}, & \mathrm{inverse} = \mathrm{false} \\
\mathbf{C}^\top\mathbf{u}, & \mathrm{inverse} = \mathrm{true}
\end{cases}
```

### External Equations

None.

## Initialization

From the initialized input,

```math
\mathbf{y} \leftarrow
\begin{cases}
\mathbf{C}\mathbf{u}, & \mathrm{inverse} = \mathrm{false} \\
\mathbf{C}^\top\mathbf{u}, & \mathrm{inverse} = \mathrm{true}
\end{cases}
```

If initial output derivatives are required:

```math
\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t} \leftarrow
\begin{cases}
\mathbf{C}\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t},
  & \mathrm{inverse} = \mathrm{false} \\
\mathbf{C}^\top\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t},
  & \mathrm{inverse} = \mathrm{true}
\end{cases}
```

Required input derivatives must be available and finite (zero for constants).

## Monitors

None.
