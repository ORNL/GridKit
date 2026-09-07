# Rotation Model

`Rotation` transforms stationary $\alpha\beta0$ quantities to rotating $dq0$
coordinates and preserves the zero component.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$\mathrm{inverse}$ | [-] | `inverse` | Inverse transformation | Default `false`

### Parameter Validation

The `inverse` parameter must be Boolean and is fixed during simulation.

### Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `input` | Input | $[u]$ | Input vector port | $\mathbf{u} \in \mathbb{R}^3$
$\theta$ | `theta` | Input | [rad] | Electrical reference angle | d-axis relative to the $\alpha$-axis
$\mathbf{y}$ | `out` | Output | $[u]$ | Algebraic output vector | $\mathbf{y} \in \mathbb{R}^3$

Forward input and output orders are $(\alpha,\beta,0)$ and $(d,q,0)$.
The inverse transformation exchanges these orders and uses the same $\theta$.

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
$\theta$ | [rad] | Electrical reference angle | Differential-input configuration

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | $[u]$ | Input vector | Algebraic-input configuration, $\mathbf{u} \in \mathbb{R}^3$
$\theta$ | [rad] | Electrical reference angle | Algebraic-input configuration

## Model Equations

The transformation matrix is

```math
\mathbf{R}(\theta) =
\begin{bmatrix}
\cos\theta & \sin\theta & 0 \\
-\sin\theta & \cos\theta & 0 \\
0 & 0 & 1
\end{bmatrix}
```

### Internal Equations

#### Differential

None.

#### Algebraic

```math
0 = -\mathbf{y} +
\begin{cases}
\mathbf{R}(\theta)\mathbf{u}, & \mathrm{inverse} = \mathrm{false} \\
\mathbf{R}^\top(\theta)\mathbf{u}, & \mathrm{inverse} = \mathrm{true}
\end{cases}
```

### External Equations

None.

## Initialization

From the initialized $\mathbf{u}$ and $\theta$,

```math
\mathbf{y} \leftarrow
\begin{cases}
\mathbf{R}(\theta)\mathbf{u}, & \mathrm{inverse} = \mathrm{false} \\
\mathbf{R}^\top(\theta)\mathbf{u}, & \mathrm{inverse} = \mathrm{true}
\end{cases}
```

If initial output derivatives are required:

```math
\dfrac{\mathrm{d}\mathbf{y}}{\mathrm{d}t} \leftarrow
\begin{cases}
\mathbf{R}(\theta)\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
  + \dfrac{\partial\mathbf{R}(\theta)}{\partial\theta}
    \mathbf{u}\dfrac{\mathrm{d}\theta}{\mathrm{d}t},
  & \mathrm{inverse} = \mathrm{false} \\
\mathbf{R}^\top(\theta)\dfrac{\mathrm{d}\mathbf{u}}{\mathrm{d}t}
  + \dfrac{\partial\mathbf{R}^\top(\theta)}{\partial\theta}
    \mathbf{u}\dfrac{\mathrm{d}\theta}{\mathrm{d}t},
  & \mathrm{inverse} = \mathrm{true}
\end{cases}
```

Required input derivatives must be available and finite (zero for constants).

## Monitors

None.
