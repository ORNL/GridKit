# Pll Model

`Pll` is a synchronous-reference-frame phase-locked loop on a three-phase bus.
It drives the q-axis voltage $v_q$ to zero, so that $\theta$ and $\omega$
track the angle and frequency of the bus voltage $\mathbf{v}$.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$K_P$ | $[\mathrm{rad}/(\mathrm{V\,s})]$ | `Kp` | Proportional gain | Required, positive
$K_I$ | $[\mathrm{rad}/(\mathrm{V\,s^2})]$ | `Ki` | Integral gain | Required, positive

### Parameter Validation

```math
\begin{aligned}
K_P &\in \mathbb{R}_{>0} \\
K_I &\in \mathbb{R}_{>0}
\end{aligned}
```

### Model Derived Parameters

```math
\omega^\mathrm{base} = 2\pi f^\mathrm{base}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `bus` | Bus | [V] | Bus voltage, phases a, b and c | $\mathbf{v} \in \mathbb{R}^3$

## Submodels

None.

### Submodel Validation

None.

### Submodel Wiring

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\theta$ | [rad] | Angle |
$x$ | [rad/s] | Integrator state |

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\omega$ | [rad/s] | Frequency |

### External Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Bus voltage vector owned by EMT bus | $\mathbf{v} \in \mathbb{R}^3$

#### Algebraic

None.

## Model Equations

```math
v_q = -\dfrac{2}{3}(v_a\sin\theta + v_b\sin(\theta - \tfrac{2\pi}{3}) + v_c\sin(\theta + \tfrac{2\pi}{3}))
```

### Internal Equations

#### Differential

```math
\begin{aligned}
0 &= -\dfrac{\mathrm{d}\theta}{\mathrm{d}t} + \omega \\
0 &= -\dfrac{\mathrm{d}x}{\mathrm{d}t} + K_I v_q
\end{aligned}
```

#### Algebraic

```math
0 = -\omega + \omega^\mathrm{base} + K_P v_q + x
```

### External Equations

None.

## Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`theta` | [rad] | Angle |
`x` | [rad/s] | Integrator state |
`omega` | [rad/s] | Frequency |
