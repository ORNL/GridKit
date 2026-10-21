# InverseTimeRelay

An inverse-time overcurrent relay with optional lockout.

## Notes

- Above pickup, the trip rises at $x = 3/4$ once
  $\int t_\mathrm{op}^{-1}\,dt = 1$, the operate criterion of IEEE
  C37.112[^inverse-time].
- With lockout, the latch commits at $x = 1/2$. A committed relay holds the
  trip after the current falls; an uncommitted relay resets.
- With lockout, the trip rises at $x = 3/4$, after the latch commits, so the
  relay may measure the current it interrupts, such as a
  [BranchBreakers](../../Branch/BranchBreakers/README.md) terminal current.
- Without lockout, the relay resets after the current falls.

## Model Parameters

Symbol              | Units     | JSON      | Description       | Typical Value | Note
--------------------|-----------|-----------|-------------------|---------------|---------------------------------------
$I_\mathrm{pickup}$ | [p.u.]    | `Ipickup` | Pickup current    |               | Required
$D$                 | [-]       | `TD`      | Time dial         |               | Required
$A$                 | [s]       | `A`       | Curve coefficient |               | Required
$B$                 | [s]       | `B`       | Curve offset      |               | Required
$\alpha$            | [-]       | `alpha`   | Curve exponent    |               | Required
$t_r$               | [s]       | `tr`      | Reset constant    | 0.0           | 0 for instantaneous reset
$s_\mathrm{lock}$   | [boolean] | `Lockout` | Lockout selector  | `true`        | `true` = lockout, `false` = self-reset

IEEE C37.112 curves:

Curve              | `A`    | `B`    | `alpha` | `tr`
-------------------|--------|--------|---------|-----
Moderately inverse | 0.0515 | 0.1140 | 0.02    | 4.85
Very inverse       | 19.61  | 0.491  | 2.0     | 21.6
Extremely inverse  | 28.2   | 0.1217 | 2.0     | 29.1

### Parameter Validation

```math
\begin{aligned}
  I_\mathrm{pickup} &> 0 \\
  D &> 0 \\
  A &> 0 \\
  B &\ge 0 \\
  \alpha &> 0 \\
  t_r &\ge 0
\end{aligned}
```

All real parameters must be finite.

### Model Derived Parameters

```math
\begin{aligned}
  \epsilon_T &= 10^{-3}\ \mathrm{s} \\
  T_\mathrm{reset} &= \max(D t_r, \epsilon_T)
\end{aligned}
```

## Model Ports

Name   | Port   | Init  | Description
-------|--------|-------|-----------------------------------------------
`ir`   | Input  | Known | Required measured current, real component
`ii`   | Input  | Known | Required measured current, imaginary component
`trip` | Output | Known | Required trip command

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
-------|-------|-------------|-----
$x$    | [-]   | Relay latch |

#### Algebraic

Symbol | Units | Description  | Note
-------|-------|--------------|-------
$s$    | [-]   | Trip command | `trip`

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description                           | Note
-------|--------|---------------------------------------|-----
$I_r$  | [p.u.] | Measured current, real component      | `ir`
$I_i$  | [p.u.] | Measured current, imaginary component | `ii`

## Model Equations

Smooth functions: [`above`](../../../../CommonMath.md#above),
[`latch`](../../Branch/BranchBreakers/README.md#appendix-a-latch),
[`logistic`](../../../../CommonMath.md#logistic-function) ($\sigma$),
[`max`](../../../../CommonMath.md#maximum),
[`min`](../../../../CommonMath.md#minimum).

For readability, define:

```math
\begin{aligned}
  M^2 &= \dfrac{I_r^2 + I_i^2}{I_\mathrm{pickup}^2} \\
  p &= \text{logistic}(M^2 - 1) \\
  K &= \text{max}(M^2, 1)^{\alpha/2} - 1 \\
  \dfrac{1}{t_\mathrm{op}} &= \dfrac{K}{D(A + B K)} \\
  \dfrac{1}{t_\mathrm{rs}} &= \dfrac{1 - \text{min}(M^2, 1)}{T_\mathrm{reset}} \\
  \dfrac{1}{T_\mathrm{latch}} &= \ln 4\left(\dfrac{p}{t_\mathrm{op}} + \dfrac{1 - p}{t_\mathrm{rs}}\right)
\end{aligned}
```

### Internal Equations

#### Differential

```math
\begin{aligned}
  0 &= -\dot{x} + \dfrac{1}{T_\mathrm{latch}}\text{latch}(x, p, (1 - s_\mathrm{lock})(1 - p))
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
  0 &= -s + \text{above}(x; 3/4)
\end{aligned}
```

### External Equations

None.

## Initialization

```math
\begin{aligned}
  x &\leftarrow 0 \\
  s &\leftarrow \text{above}(x; 3/4)
\end{aligned}
```

All internal derivatives initialize to zero.

## Monitors

Monitor | Units  | Description                | Note
--------|--------|----------------------------|-----------------------
`im`    | [p.u.] | Measured current magnitude | $\sqrt{I_r^2 + I_i^2}$
`x`     | [-]    | Relay latch                | $x$
`trip`  | [-]    | Trip command               | $s$

[^inverse-time]: [IEEE Std C37.112-2018, *IEEE Standard for Inverse-Time Characteristics Equations for Overcurrent Relays*](https://standards.ieee.org/ieee/C37.112/7036/).
