# Park Model

`Park` gives three-phase quantities in coordinates that rotate at an angle
$\theta$, and the unit direction of the result. The transform keeps the
amplitude of a balanced set, and the q axis leads the d axis by $\pi/2$.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$U^\mathrm{base}$ | $[u]$ | `base` | Rated peak magnitude | Default $1$, positive

### Parameter Validation

```math
U^\mathrm{base} \in \mathbb{R}_{>0}
```

### Model Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{u}$ | `abc` | Input | $[u]$ | Three-phase quantities, phases a, b and c | $\mathbf{u} \in \mathbb{R}^3$
$\theta$ | `theta` | Input | [rad] | Angle of the d axis from the axis of phase a | Scalar
$\mathbf{u}_\mathrm{dq}$ | `dq` | Output | $[u]$ | Input in rotating coordinates | Internal variable
$(c, s)$ | `axis` | Output | [-] | Unit direction of the input in rotating coordinates | Internal variable

The angle may come from [Pll](../Pll/README.md).

## Submodels

None.

### Submodel Validation

None.

### Submodel Wiring

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}_\mathrm{dq}$ | $[u]$ | Input in rotating coordinates | $\mathbf{u}_\mathrm{dq} \in \mathbb{R}^2$
$(c, s)$ | [-] | Unit direction of the input in rotating coordinates |

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{u}$ | $[u]$ | Three-phase quantities | $\mathbf{u} \in \mathbb{R}^3$
$\theta$ | [rad] | Angle |

## Model Equations

```math
\begin{aligned}
u_\alpha &= \dfrac{1}{3}(2u_a - u_b - u_c), \qquad
u_\beta = \dfrac{1}{\sqrt{3}}(u_b - u_c) \\
r &= \dfrac{1}{U^\mathrm{base}}\sqrt{u_d^2 + u_q^2}
\end{aligned}
```

### Internal Equations

#### Differential

None.

#### Algebraic

```math
\begin{aligned}
0 &= -u_d + u_\alpha\cos\theta + u_\beta\sin\theta \\
0 &= -u_q + u_\beta\cos\theta - u_\alpha\sin\theta
\end{aligned}
```

For $r > 10^{-6}$,

```math
\begin{aligned}
0 &= -c + \dfrac{1}{r\,U^\mathrm{base}}u_d \\
0 &= -s + \dfrac{1}{r\,U^\mathrm{base}}u_q
\end{aligned}
```

Otherwise the direction is the d axis, $(c, s) = (1, 0)$.

### External Equations

None.

## Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`dq` | $[u]$ | Input in rotating coordinates | $\mathbf{u}_\mathrm{dq} \in \mathbb{R}^2$
`axis` | [-] | Unit direction of the input in rotating coordinates |
