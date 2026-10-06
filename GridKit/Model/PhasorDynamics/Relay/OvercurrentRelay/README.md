# OvercurrentRelay

OvercurrentRelay is a definite-time overcurrent relay with lockout. It reads a
measured current from two signals and issues a trip command once the current
magnitude has exceeded pickup for the trip time.

## Notes

- The lockout latch commits at $x=1/2$. A committed relay trips and holds the
  trip after the current falls; an uncommitted relay resets.
- The trip command rises at $x=3/4$, after the latch commits, so the relay may
  measure the current it interrupts, such as a
  [BranchBreakers](../../Branch/BranchBreakers/README.md) terminal current.

## Model Parameters

Symbol              | Units  | JSON      | Description                        | Typical Value | Note
--------------------|--------|-----------|------------------------------------|---------------|---------
$I_\mathrm{pickup}$ | [p.u.] | `Ipickup` | Pickup current magnitude           |               | Required
$T_\mathrm{trip}$   | [s]    | `Ttrip`   | Time from sustained pickup to trip |               | Required

### Parameter Validation

A valid OvercurrentRelay parameter set must satisfy the following conditions:

```math
\begin{aligned}
  &I_\mathrm{pickup}, T_\mathrm{trip} \in \mathbb{R}\ \text{and finite} \\
  &I_\mathrm{pickup} > 0 \\
  &T_\mathrm{trip} \ge 0
\end{aligned}
```

### Model Derived Parameters

Let $\epsilon_T=10^{-3}\ \mathrm{s}$. A trip time below $\epsilon_T$ is raised
to that floor in place.

```math
\begin{aligned}
  T_\mathrm{trip} &\leftarrow \max(T_\mathrm{trip},\epsilon_T) \\
  T_\mathrm{latch} &= T_\mathrm{trip}/\ln 4
\end{aligned}
```

Full pickup moves the latch from zero to the trip level $x=3/4$ in
$T_\mathrm{trip}$.

## Model Ports

Name   | Port   | Init  | Description
-------|--------|-------|-----------------------------------------------
`ir`   | Input  | Known | Required measured current, real component
`ii`   | Input  | Known | Required measured current, imaginary component
`trip` | Output | Known | Required trip command

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description   | Note
-------|-------|---------------|-------------
$x$    | [-]   | Lockout latch | Zero at rest

#### Algebraic

Symbol | Units | Description  | Note
-------|-------|--------------|-------------
$s$    | [-]   | Trip command | `trip` output

### External Variables

#### Differential

None.

#### Algebraic

Symbol         | Units  | Description                           | Note
---------------|--------|---------------------------------------|-----------
$I_\mathrm{r}$ | [p.u.] | Measured current, real component      | `ir` input
$I_\mathrm{i}$ | [p.u.] | Measured current, imaginary component | `ii` input

## Model Equations

Smooth functions: [`above`](../../../../CommonMath.md#above), [`latch`](../../../../CommonMath.md#latch).

Define the pickup indicator:

```math
\begin{aligned}
  p &= \text{above}\big((I_\mathrm{r}^2+I_\mathrm{i}^2)/I_\mathrm{pickup}^2;1\big)
\end{aligned}
```

### Internal Equations

#### Differential

```math
\begin{aligned}
  0 &= -\dot{x} + \text{latch}(x,p,0)/T_\mathrm{latch}
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
  0 &= -s + \text{above}(x;3/4)
\end{aligned}
```

### External Equations

None.

## Initialization

### Internal Initialization

```math
\begin{aligned}
  x &\leftarrow 0 \\
  s &\leftarrow \text{above}(x;3/4)
\end{aligned}
```

The relay starts reset, and the latch derivative initializes to zero.

## Monitors

Monitor | Units  | Description                | Note
--------|--------|----------------------------|---------------------------------------
`im`    | [p.u.] | Measured current magnitude | $\sqrt{I_\mathrm{r}^2+I_\mathrm{i}^2}$
`x`     | [-]    | Lockout latch              | $x$
`trip`  | [-]    | Trip command               | $s$
