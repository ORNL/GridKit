# DefiniteTimeRelay

A definite-time relay with optional lockout. Each relay models one
[element](../README.md#elements).

## Notes

- Designed to measure a [BranchBreakers](../../Branch/BranchBreakers/README.md)
  terminal current, a bus voltage, or a frequency, and to drive a BranchBreakers
  `trip1` or `trip2` input.
- The element operates above $M_\mathrm{pickup}$ if
  $M_\mathrm{pickup} > M_\mathrm{ref}$ and below it otherwise.
- With lockout, the relay commits to trip once $x$ passes $1/2$, so it may
  measure the current it interrupts.

## Model Parameters

Symbol              | Units     | JSON      | Description      | Typical Value | Note
--------------------|-----------|-----------|------------------|---------------|-----------------------------------------------------
$M_\mathrm{pickup}$ | [p.u.]    | `Mpickup` | Pickup level     |               | Required
$M_\mathrm{ref}$    | [p.u.]    | `Mref`    | Reference level  |               | Required; 0 for current, 1 for voltage and frequency
$T_\mathrm{trip}$   | [s]       | `Ttrip`   | Trip time        |               | Required; time to $x = 3/4$ for $p = 1$
$s_\mathrm{lock}$   | [boolean] | `Lockout` | Lockout selector | `true`        | `true` = lockout, `false` = self-reset

### Parameter Validation

```math
\begin{aligned}
  M_\mathrm{pickup} &> 0 \\
  M_\mathrm{ref} &\ge 0 \\
  M_\mathrm{pickup} &\ne M_\mathrm{ref} \\
  T_\mathrm{trip} &\ge 0
\end{aligned}
```

All real parameters must be finite.

### Model Derived Parameters

```math
\begin{aligned}
  \epsilon_T &= 10^{-3}\ \mathrm{s} \\
  T_\mathrm{trip} &\leftarrow \max(T_\mathrm{trip}, \epsilon_T) \\
  T_\mathrm{latch} &= \dfrac{T_\mathrm{trip}}{\ln 4}
\end{aligned}
```

## Model Ports

Name   | Port   | Init  | Description
-------|--------|-------|-------------------------------------------
`mr`   | Input  | Known | Required measured quantity, real component
`mi`   | Input  | Known | Measured quantity, imaginary component
`trip` | Output | Known | Required trip command

Unattached `mi` defaults to zero, so absolute frequency uses `mr` alone.

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

Symbol | Units  | Description                            | Note
-------|--------|----------------------------------------|-----
$m_r$  | [p.u.] | Measured quantity, real component      | `mr`
$m_i$  | [p.u.] | Measured quantity, imaginary component | `mi`

## Model Equations

Smooth functions: [`above`](../../../../CommonMath.md#above),
[`latch`](../../../../CommonMath.md#latch),
[`logistic`](../../../../CommonMath.md#logistic-function) ($\sigma$).

For readability, define:

```math
\begin{aligned}
  p &= \text{logistic}\left(\dfrac{m_r^2 + m_i^2 - M_\mathrm{pickup}^2}{M_\mathrm{pickup}^2 - M_\mathrm{ref}^2}\right)
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

Monitor | Units  | Description                 | Note
--------|--------|-----------------------------|-----------------------
`m`     | [p.u.] | Measured quantity magnitude | $\sqrt{m_r^2 + m_i^2}$
`x`     | [-]    | Relay latch                 | $x$
`trip`  | [-]    | Trip command                | $s$
