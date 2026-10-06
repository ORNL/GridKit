# BranchBreakers

BranchBreakers is a [Branch](../README.md) $\pi$ branch with a circuit breaker
at each terminal. Trip and reset commands arrive as signals, and an open side
is removed from the branch by Kron reduction of the two-port admittance.
Terminal current contributions are oriented entering the adjacent buses.

## Notes

- Each breaker is a closed-state latch, so a trip command has priority over a
  simultaneous reset command and a breaker holds its position once its
  commands return to zero.
- The terminal-current outputs publish values only; a model reading them sees
  no Jacobian coupling to the bus voltages or breaker states.
- The `enable` and `disable` solver events apply to [Branch](../README.md)
  only.

## Model Parameters

Symbol           | Units  | JSON    | Description                             | Typical Value | Note
-----------------|--------|---------|-----------------------------------------|---------------|-----------------------
$R$              | [p.u.] | `R`     | Branch series resistance                |               |
$X$              | [p.u.] | `X`     | Branch series reactance                 |               |
$G$              | [p.u.] | `G`     | Total line shunt conductance            | 0.0           |
$B$              | [p.u.] | `B`     | Total line shunt susceptance            | 0.0           |
$G_\mathrm{mag}$ | [p.u.] | `Gmag`  | Magnetizing shunt conductance at bus 1  | 0.0           |
$B_\mathrm{mag}$ | [p.u.] | `Bmag`  | Magnetizing shunt susceptance at bus 1  | 0.0           |
$\tau$           | [p.u.] | `tap`   | Off-nominal tap magnitude on bus-1 side | 1.0           |
$\theta$         | [rad]  | `phase` | Phase-shift angle                       | 0.0           |
$T_\mathrm{brk}$ | [s]    | `Tbrk`  | Breaker operating time                  | 0.05          | Command to half travel

### Parameter Validation

A valid BranchBreakers parameter set must satisfy the following conditions:

```math
\begin{aligned}
  &R, X, G, B, G_\mathrm{mag}, B_\mathrm{mag}, \tau, \theta, T_\mathrm{brk} \in \mathbb{R}\ \text{and finite} \\
  &R^2 + X^2 > 0 \\
  &\tau > 0 \\
  &T_\mathrm{brk} \ge 0 \\
  &Y_{11}, Y_{22} \ne 0
\end{aligned}
```

### Model Derived Parameters

Let $\epsilon_T=10^{-3}\ \mathrm{s}$. A breaker operating time below
$\epsilon_T$ is raised to that floor in place. The closed two-port admittance
$\mathbf{Y}$ is the [Branch](../README.md#model-derived-parameters) admittance
with $u=1$.

```math
\begin{aligned}
  T_\mathrm{brk} &\leftarrow \max(T_\mathrm{brk},\epsilon_T) \\
  T_\mathrm{latch} &= T_\mathrm{brk}/\ln 2 \\
  K_1 &= Y_{12}Y_{21}/Y_{22} \\
  K_2 &= Y_{21}Y_{12}/Y_{11}
\end{aligned}
```

A full command moves a breaker latch through half travel in $T_\mathrm{brk}$.
$K_1$ removes the bus-2 side from $Y_{11}$ when bus 2 opens, and $K_2$ removes
the bus-1 side from $Y_{22}$ when bus 1 opens.

## Model Ports

Name     | Port   | Init  | Description
---------|--------|-------|------------------------------------------------------
`bus1`   | Bus    | Known | Required bus-1 terminal; the tapped side
`bus2`   | Bus    | Known | Required bus-2 terminal
`trip1`  | Input  | Known | Optional bus-1 breaker trip command; defaults to zero
`reset1` | Input  | Known | Optional bus-1 breaker reset command; defaults to zero
`trip2`  | Input  | Known | Optional bus-2 breaker trip command; defaults to zero
`reset2` | Input  | Known | Optional bus-2 breaker reset command; defaults to zero
`ir1`    | Output | Known | Optional bus-1 terminal current, real component
`ii1`    | Output | Known | Optional bus-1 terminal current, imaginary component
`ir2`    | Output | Known | Optional bus-2 terminal current, real component
`ii2`    | Output | Known | Optional bus-2 terminal current, imaginary component

Attached inputs must be linked.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description                      | Note
-------|-------|----------------------------------|----------------
$z_1$  | [-]   | Bus-1 breaker closed-state latch | One when closed
$z_2$  | [-]   | Bus-2 breaker closed-state latch | One when closed

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol            | Units  | Description                                  | Note
------------------|--------|----------------------------------------------|------------------------------------
$V_{\mathrm{r}1}$ | [p.u.] | Terminal voltage, real component, bus 1      | Owned by bus object
$V_{\mathrm{i}1}$ | [p.u.] | Terminal voltage, imaginary component, bus 1 | Owned by bus object
$V_{\mathrm{r}2}$ | [p.u.] | Terminal voltage, real component, bus 2      | Owned by bus object
$V_{\mathrm{i}2}$ | [p.u.] | Terminal voltage, imaginary component, bus 2 | Owned by bus object
$s_1$             | [-]    | Bus-1 breaker trip command                   | Optional `trip1`; defaults to zero
$r_1$             | [-]    | Bus-1 breaker reset command                  | Optional `reset1`; defaults to zero
$s_2$             | [-]    | Bus-2 breaker trip command                   | Optional `trip2`; defaults to zero
$r_2$             | [-]    | Bus-2 breaker reset command                  | Optional `reset2`; defaults to zero

## Model Equations

Smooth functions: [`above`](../../../../CommonMath.md#above), [`latch`](../../../../CommonMath.md#latch).

Define the breaker closed fractions, the switched admittance, and the terminal
currents:

```math
\begin{aligned}
  u_1 &= \text{above}(z_1;1/2) \\
  u_2 &= \text{above}(z_2;1/2) \\
  \hat{Y}_{11} &= u_1(Y_{11}-(1-u_2)K_1) \\
  \hat{Y}_{12} &= u_1u_2Y_{12} \\
  \hat{Y}_{21} &= u_1u_2Y_{21} \\
  \hat{Y}_{22} &= u_2(Y_{22}-(1-u_1)K_2) \\
  I_{\mathrm{r}1} &= \hat{G}_{11}V_{\mathrm{r}1} - \hat{B}_{11}V_{\mathrm{i}1} + \hat{G}_{12}V_{\mathrm{r}2} - \hat{B}_{12}V_{\mathrm{i}2} \\
  I_{\mathrm{i}1} &= \hat{B}_{11}V_{\mathrm{r}1} + \hat{G}_{11}V_{\mathrm{i}1} + \hat{B}_{12}V_{\mathrm{r}2} + \hat{G}_{12}V_{\mathrm{i}2} \\
  I_{\mathrm{r}2} &= \hat{G}_{21}V_{\mathrm{r}1} - \hat{B}_{21}V_{\mathrm{i}1} + \hat{G}_{22}V_{\mathrm{r}2} - \hat{B}_{22}V_{\mathrm{i}2} \\
  I_{\mathrm{i}2} &= \hat{B}_{21}V_{\mathrm{r}1} + \hat{G}_{21}V_{\mathrm{i}1} + \hat{B}_{22}V_{\mathrm{r}2} + \hat{G}_{22}V_{\mathrm{i}2}
\end{aligned}
```

Each switched entry is written $\hat{Y}_{mn}=\hat{G}_{mn}+j\hat{B}_{mn}$. At
each corner of $(u_1,u_2)\in\{0,1\}^2$, $\hat{\mathbf{Y}}$ is the Kron
reduction of the open sides, and both sides open give $\hat{\mathbf{Y}}=0$.

### Internal Equations

#### Differential

```math
\begin{aligned}
  0 &= -\dot{z}_1 + \text{latch}(z_1,r_1,s_1)/T_\mathrm{latch} \\
  0 &= -\dot{z}_2 + \text{latch}(z_2,r_2,s_2)/T_\mathrm{latch}
\end{aligned}
```

The reset command is the latch set drive and the trip command is its reset
drive, so a trip command has priority.

#### Algebraic

None.

### External Equations

```math
\begin{aligned}
  \Delta I^{\mathrm{bus}1}_\mathrm{r} &\mathrel{+}= I_{\mathrm{r}1} \\
  \Delta I^{\mathrm{bus}1}_\mathrm{i} &\mathrel{+}= I_{\mathrm{i}1} \\
  \Delta I^{\mathrm{bus}2}_\mathrm{r} &\mathrel{+}= I_{\mathrm{r}2} \\
  \Delta I^{\mathrm{bus}2}_\mathrm{i} &\mathrel{+}= I_{\mathrm{i}2}
\end{aligned}
```

## Initialization

### Input Initialization

```math
\begin{aligned}
  V_{\mathrm{r}1}, V_{\mathrm{i}1} &\leftarrow \text{bus-1 voltage} \\
  V_{\mathrm{r}2}, V_{\mathrm{i}2} &\leftarrow \text{bus-2 voltage}
\end{aligned}
```

### Internal Initialization

```math
\begin{aligned}
  z_1, z_2 &\leftarrow 1
\end{aligned}
```

Both breakers start closed, and the latch derivatives initialize to zero.

### Output Initialization

```math
\begin{aligned}
  I_{\mathrm{r}1} &\leftarrow G_{11}V_{\mathrm{r}1} - B_{11}V_{\mathrm{i}1} + G_{12}V_{\mathrm{r}2} - B_{12}V_{\mathrm{i}2} \\
  I_{\mathrm{i}1} &\leftarrow B_{11}V_{\mathrm{r}1} + G_{11}V_{\mathrm{i}1} + B_{12}V_{\mathrm{r}2} + G_{12}V_{\mathrm{i}2} \\
  I_{\mathrm{r}2} &\leftarrow G_{21}V_{\mathrm{r}1} - B_{21}V_{\mathrm{i}1} + G_{22}V_{\mathrm{r}2} - B_{22}V_{\mathrm{i}2} \\
  I_{\mathrm{i}2} &\leftarrow B_{21}V_{\mathrm{r}1} + G_{21}V_{\mathrm{i}1} + B_{22}V_{\mathrm{r}2} + G_{22}V_{\mathrm{i}2}
\end{aligned}
```

The current outputs publish the closed-branch currents and are refreshed on
every residual evaluation.

## Monitors

Monitor | Units  | Description                                  | Note
--------|--------|----------------------------------------------|-------------------------------------------------------------------------
`ir1`   | [p.u.] | Terminal current, real component, bus 1      | $I_{\mathrm{r}1}$; oriented entering bus 1
`ii1`   | [p.u.] | Terminal current, imaginary component, bus 1 | $I_{\mathrm{i}1}$; oriented entering bus 1
`im1`   | [p.u.] | Terminal current magnitude, bus 1            | $\sqrt{I_{\mathrm{r}1}^2+I_{\mathrm{i}1}^2}$
`p1`    | [p.u.] | Active power at bus 1 terminal               | $V_{\mathrm{r}1}I_{\mathrm{r}1}+V_{\mathrm{i}1}I_{\mathrm{i}1}$; positive entering bus 1
`q1`    | [p.u.] | Reactive power at bus 1 terminal             | $V_{\mathrm{i}1}I_{\mathrm{r}1}-V_{\mathrm{r}1}I_{\mathrm{i}1}$; positive entering bus 1
`ir2`   | [p.u.] | Terminal current, real component, bus 2      | $I_{\mathrm{r}2}$; oriented entering bus 2
`ii2`   | [p.u.] | Terminal current, imaginary component, bus 2 | $I_{\mathrm{i}2}$; oriented entering bus 2
`im2`   | [p.u.] | Terminal current magnitude, bus 2            | $\sqrt{I_{\mathrm{r}2}^2+I_{\mathrm{i}2}^2}$
`p2`    | [p.u.] | Active power at bus 2 terminal               | $V_{\mathrm{r}2}I_{\mathrm{r}2}+V_{\mathrm{i}2}I_{\mathrm{i}2}$; positive entering bus 2
`q2`    | [p.u.] | Reactive power at bus 2 terminal             | $V_{\mathrm{i}2}I_{\mathrm{r}2}-V_{\mathrm{r}2}I_{\mathrm{i}2}$; positive entering bus 2
