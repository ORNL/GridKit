# Park Model

`Park` transforms three-phase quantities from $abc$ to rotating $dq0$
coordinates with power-invariant normalization. The operator adds no DAE
variables or residual rows.

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
$\theta$ | `theta` | Input | [rad] | Electrical reference angle | d-axis relative to the phase-a axis
$\mathbf{y}$ | `out` | Output | $[u]$ | Algebraic output vector | $\mathbf{y} \in \mathbb{R}^3$

Forward input and output orders are $(a,b,c)$ and $(d,q,0)$.
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

None.

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
\mathbf{T}(\theta) =
\sqrt{\dfrac{2}{3}}
\begin{bmatrix}
\cos\theta
  & \cos\left(\theta-\dfrac{2\pi}{3}\right)
  & \cos\left(\theta+\dfrac{2\pi}{3}\right) \\
-\sin\theta
  & -\sin\left(\theta-\dfrac{2\pi}{3}\right)
  & -\sin\left(\theta+\dfrac{2\pi}{3}\right) \\
\dfrac{1}{\sqrt{2}} & \dfrac{1}{\sqrt{2}} & \dfrac{1}{\sqrt{2}}
\end{bmatrix}
```

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

```math
\mathbf{y} \leftarrow
\begin{cases}
\mathbf{T}(\theta)\mathbf{u}, & \mathrm{inverse} = \mathrm{false} \\
\mathbf{T}^\mathsf{T}(\theta)\mathbf{u}, & \mathrm{inverse} = \mathrm{true}
\end{cases}
```

## Initialization

The preparation pass publishes the transformed operating-point values from
$\mathbf{u}$ and $\theta$; prescribed outputs must agree:

```math
\mathbf{y} \leftarrow
\begin{cases}
\mathbf{T}(\theta)\mathbf{u}, & \mathrm{inverse} = \mathrm{false} \\
\mathbf{T}^\mathsf{T}(\theta)\mathbf{u}, & \mathrm{inverse} = \mathrm{true}
\end{cases}
```

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`out` | $[u]$ | Transformed output | $\mathbf{y} \in \mathbb{R}^3$
