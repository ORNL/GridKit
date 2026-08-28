# **Branch Model**

The Branch model represents a two-terminal phasor-domain $\pi$ branch with an
optional off-nominal tap magnitude and phase shift on bus 1. Terminal current
contributions are oriented entering the adjacent buses.

Notes:
- Setting $\tau = 1$ and $\theta = 0$ gives the ordinary symmetric
  transmission-line $\pi$ model.
- The total line shunt $G + jB$ is split equally between the two terminals,
  while the magnetizing shunt $G_\text{mag} + jB_\text{mag}$ is connected at
  bus 1; both shunts are added outside the $\mathbf{M}$ transformation.
- The branch has no solver-owned variables; it contributes current residuals
  directly to the connected buses.

## Model Parameters

Symbol               | Units  | JSON    | Description                             | Typical Value | Note
---------------------|--------|---------|-----------------------------------------|---------------|------
$R$                  | [p.u.] | `R`     | Branch series resistance                |               |
$X$                  | [p.u.] | `X`     | Branch series reactance                 |               |
$G$                  | [p.u.] | `G`     | Total line shunt conductance            | 0.0           |
$B$                  | [p.u.] | `B`     | Total line shunt susceptance            | 0.0           |
$G_\text{mag}$       | [p.u.] | `Gmag`  | Magnetizing shunt conductance at bus 1  | 0.0           |
$B_\text{mag}$       | [p.u.] | `Bmag`  | Magnetizing shunt susceptance at bus 1  | 0.0           |
$\tau$               | [p.u.] | `tap`   | Off-nominal tap magnitude on bus-1 side | 1.0           |
$\theta$             | [rad]  | `phase` | Phase-shift angle                       | 0.0           |

### Parameter Validation

A valid Branch parameter set must satisfy the following conditions:

```math
\begin{aligned}
  &R, X, G, B, G_\text{mag}, B_\text{mag}, \tau, \theta
    \in \mathbb{R}\ \text{and finite} \\
  &R^2 + X^2 > 0 \\
  &\tau > 0
\end{aligned}
```

### Model Derived Parameters

The series, magnetizing shunt, and line shunt admittances are:

```math
\begin{aligned}
  Y_{\mathrm{br}}  &= \dfrac{1}{R + jX} \\
  Y_{\mathrm{mag}} &= G_\text{mag} + jB_\text{mag} \\
  Y_{\mathrm{sh}}  &= G + jB
\end{aligned}
```

The series, magnetizing shunt, and line shunt admittance matrices are:

```math
\begin{aligned}
  \mathbf{Y}_0
    &=
    \begin{bmatrix}
      -Y_{\mathrm{br}} & Y_{\mathrm{br}} \\
      Y_{\mathrm{br}}
      & -Y_{\mathrm{br}}
    \end{bmatrix}
    \\
  \mathbf{Y}_\text{mag}
    &=
    \begin{bmatrix}
      -Y_{\mathrm{mag}} & 0 \\
      0 & 0
    \end{bmatrix}
    \\
  \mathbf{Y}_\text{sh}
    &=
    \dfrac{1}{2}
    \begin{bmatrix}
      -Y_{\mathrm{sh}} & 0 \\
      0 & -Y_{\mathrm{sh}}
    \end{bmatrix}
\end{aligned}
```

The off-nominal transformer transformation uses bus 1 as the tap side:

```math
\begin{aligned}
  \mathbf{M}
    &=
    \begin{bmatrix}
      \tau^{-1} & 0 \\
      0 & e^{j\theta}
    \end{bmatrix}
\end{aligned}
```

The magnetizing and line shunts are added outside the transformation:

```math
\begin{aligned}
  \mathbf{Y}
    &=
    \mathbf{M}^{\dagger}
    \mathbf{Y}_0
    \mathbf{M}
    +
    \mathbf{Y}_\text{mag}
    +
    \mathbf{Y}_\text{sh}
\end{aligned}
```

For the equations below, write each entry as $Y_{mn}=G_{mn}+jB_{mn}$.

## Model Ports

Name   | Port | Init  | Description
-------|------|-------|------------
`bus1` | Bus  | Known | Required bus-1 terminal; the tapped side
`bus2` | Bus  | Known | Required bus-2 terminal

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol      | Units   | Description                             | Note
------------|---------|-----------------------------------------|------
$V_{r1}$    | [p.u.]  | Terminal voltage, real component, bus 1 | Owned by bus object
$V_{i1}$    | [p.u.]  | Terminal voltage, imaginary component, bus 1 | Owned by bus object
$V_{r2}$    | [p.u.]  | Terminal voltage, real component, bus 2 | Owned by bus object
$V_{i2}$    | [p.u.]  | Terminal voltage, imaginary component, bus 2 | Owned by bus object

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

The branch current relation is $0 = -\mathbf{I} + \mathbf{Y}\mathbf{V}$.

```math
\begin{aligned}
  I_{r1} &= G_{11} V_{r1} - B_{11} V_{i1}
         + G_{12} V_{r2} - B_{12} V_{i2} \\
  I_{i1} &= B_{11} V_{r1} + G_{11} V_{i1}
         + B_{12} V_{r2} + G_{12} V_{i2} \\
  I_{r2} &= G_{21} V_{r1} - B_{21} V_{i1}
         + G_{22} V_{r2} - B_{22} V_{i2} \\
  I_{i2} &= B_{21} V_{r1} + G_{21} V_{i1}
         + B_{22} V_{r2} + G_{22} V_{i2}
\end{aligned}
```

These current contributions are added to the connected bus residuals with
positive sign because branch current is oriented entering the bus.

## Initialization

The Branch model has no internal state to initialize. During construction or
parameter updates, the component computes $\mathbf{Y}$ from the current
parameter values. Initial terminal current and power monitor values are
evaluated from the connected bus voltages. Parameter verification enforces the
conditions in [Parameter Validation](#parameter-validation).

## Monitors

Monitor | Units  | Description                                  | Note
--------|--------|----------------------------------------------|------
`ir1`  | [p.u.] | Terminal current, real component, bus 1      | Oriented entering bus 1
`ii1`  | [p.u.] | Terminal current, imaginary component, bus 1 | Oriented entering bus 1
`im1`  | [p.u.] | Terminal current magnitude, bus 1            |
`p1`   | [p.u.] | Active power at bus 1 terminal               | Positive entering bus 1
`q1`   | [p.u.] | Reactive power at bus 1 terminal             | Positive entering bus 1
`ir2`  | [p.u.] | Terminal current, real component, bus 2      | Oriented entering bus 2
`ii2`  | [p.u.] | Terminal current, imaginary component, bus 2 | Oriented entering bus 2
`im2`  | [p.u.] | Terminal current magnitude, bus 2            |
`p2`   | [p.u.] | Active power at bus 2 terminal               | Positive entering bus 2
`q2`   | [p.u.] | Reactive power at bus 2 terminal             | Positive entering bus 2

Current magnitudes are:

```math
\begin{aligned}
  I_{m1} &= \sqrt{I_{r1}^2 + I_{i1}^2} \\
  I_{m2} &= \sqrt{I_{r2}^2 + I_{i2}^2}
\end{aligned}
```

Complex power at each end is defined as $S=VI^{\ast}$:

```math
\begin{aligned}
  P_1 &= V_{r1} I_{r1} + V_{i1} I_{i1} \\
  Q_1 &= V_{i1} I_{r1} - V_{r1} I_{i1} \\
  P_2 &= V_{r2} I_{r2} + V_{i2} I_{i2} \\
  Q_2 &= V_{i2} I_{r2} - V_{r2} I_{i2}
\end{aligned}
```
