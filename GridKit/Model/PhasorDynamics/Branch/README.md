# **Branch Model**

The Branch model represents a two-terminal phasor-domain $\pi$ branch with an
optional off-nominal tap magnitude and phase shift on bus 1. Terminal current
contributions are oriented entering the adjacent buses.

Notes:
- Setting $\tau = 1$ and $\theta = 0$ gives the ordinary symmetric
  transmission-line $\pi$ model.
- $G$ and $B$ are total branch shunt values split equally between terminals.
- The branch has no solver-owned variables; it contributes current residuals
  directly to the connected buses.
- `closed=false` is static case-import status. It removes the branch admittance
  contribution but does not insert fake admittance; islanded or underconstrained
  topology can still make the system singular.
- Sparse automatic differentiation may materialize zero-valued structural
  entries for an open branch; the mathematical Jacobian contribution is zero.

## Circuit Diagram

An ideal complex tap is placed on the bus-1 side of the branch equivalent. The
ordinary transmission-line $\pi$ model is recovered with $\tau = 1$ and
$\theta = 0$.

![](../../../../docs/Figures/transformer-branch.png)

Figure 1: Branch equivalent circuit

## Model Parameters

Symbol                 | Units  | JSON     | Description                             | Typical Value | Note
-----------------------|--------|----------|-----------------------------------------|---------------|------
$R$                    | [p.u.] | `R`      | Branch series resistance                |               |
$X$                    | [p.u.] | `X`      | Branch series reactance                 |               |
$G$                    | [p.u.] | `G`      | Total branch shunt conductance          | 0             |
$B$                    | [p.u.] | `B`      | Total branch shunt susceptance          | 0             |
$\tau$                 | [p.u.] | `tap`    | Off-nominal tap magnitude on bus-1 side | 1             |
$\theta$               | [rad]  | `phase`  | Phase-shift angle                       | 0             |
$c_{\mathrm{br}}$      | [-]    | `closed` | Static branch closed status             | `true`        | JSON boolean

### Parameter Validation

Invalid Branch parameter sets are rejected by the following checks:

```math
\begin{aligned}
  &R, X, G, B, \tau, \theta \in \mathbb{R}\ \text{and finite} \\
  &c_{\mathrm{br}} \in \{\mathrm{true}, \mathrm{false}\} \\
  &R^2 + X^2 > 0 \\
  &\tau > 0
\end{aligned}
```

### Model Derived Parameters

The closed-status factor and branch admittances are:

```math
\begin{aligned}
  s_{\mathrm{br}} &=
    \begin{cases}
      1, & c_{\mathrm{br}} = \mathrm{true} \\
      0, & c_{\mathrm{br}} = \mathrm{false}
    \end{cases} \\
  Y_{\mathrm{br}} &= \dfrac{1}{R + jX} \\
  Y_{\mathrm{sh}} &= G + jB
\end{aligned}
```

The nominal $\pi$-branch admittance matrix is the sum of the series and shunt
admittance contributions:

```math
\begin{aligned}
  \mathbf{Y}_0
    &=
    \begin{bmatrix}
      -Y_{\mathrm{br}} & Y_{\mathrm{br}} \\
      Y_{\mathrm{br}}
      & -Y_{\mathrm{br}}
    \end{bmatrix}
    +
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
    \end{bmatrix} \\
  \mathbf{Y} &= s_{\mathrm{br}}\mathbf{M}^{\dagger}\mathbf{Y}_0\mathbf{M}
\end{aligned}
```

For each entry $Y_{mn}=G_{mn}+jB_{mn}$, the real-valued contribution from
terminal $n$ to current at terminal $m$ is:

```math
\begin{aligned}
  \begin{bmatrix}
    I_{rm} \\
    I_{im}
  \end{bmatrix}_{n}
  =
  \begin{bmatrix}
    G_{mn} & -B_{mn} \\
    B_{mn} &  G_{mn}
  \end{bmatrix}
  \begin{bmatrix}
    V_{rn} \\
    V_{in}
  \end{bmatrix}
\end{aligned}
```

The voltage derivative for the same block is:

```math
\begin{aligned}
  \frac{\partial [I_{rm}, I_{im}]^T}
       {\partial [V_{rn}, V_{in}]}
  =
  \begin{bmatrix}
    G_{mn} & -B_{mn} \\
    B_{mn} &  G_{mn}
  \end{bmatrix}
\end{aligned}
```

When `closed=false`, $s_{\mathrm{br}}=0$, so $\mathbf{Y}=0$ and every current
block and voltage derivative block is zero.

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

### Differential Equations

None.

### Algebraic Equations

```math
\begin{aligned}
  0 &= -I_{r1} + G_{11} V_{r1} - B_{11} V_{i1}
         + G_{12} V_{r2} - B_{12} V_{i2} \\
  0 &= -I_{i1} + B_{11} V_{r1} + G_{11} V_{i1}
         + B_{12} V_{r2} + G_{12} V_{i2} \\
  0 &= -I_{r2} + G_{21} V_{r1} - B_{21} V_{i1}
         + G_{22} V_{r2} - B_{22} V_{i2} \\
  0 &= -I_{i2} + B_{21} V_{r1} + G_{21} V_{i1}
         + B_{22} V_{r2} + G_{22} V_{i2}
\end{aligned}
```

## Initialization

The Branch model has no internal state to initialize. During construction or
parameter updates, the component computes $\mathbf{Y}$ from the current
parameter values. Initial terminal current and power monitor values are
evaluated from the connected bus voltages. Parameter verification rejects the
invalid cases listed above.

## Model Outputs

Output | Units  | Description                                  | Note
-------|--------|----------------------------------------------|------
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
