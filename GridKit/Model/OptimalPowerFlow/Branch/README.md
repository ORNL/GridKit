# Branch

Line or off-nominal transformer between two buses.

## Model Parameters

Symbol           | Units  | Data         | Description                                           | Typical Value | Note
-----------------|--------|--------------|-------------------------------------------------------|---------------|-----
$R$              | [p.u.] | `R`          | Series resistance                                     |               | Default 0
$X$              | [p.u.] | `X`          | Series reactance                                      |               | Default 0
$G$              | [p.u.] | `G`          | Total line shunt conductance, split between the ports |               | Default 0
$B$              | [p.u.] | `B`          | Total line shunt susceptance, split between the ports |               | Default 0
$G_\mathrm{mag}$ | [p.u.] | `Gmag`       | Magnetizing shunt conductance at bus 1                |               | Default 0
$B_\mathrm{mag}$ | [p.u.] | `Bmag`       | Magnetizing shunt susceptance at bus 1                |               | Default 0
$\tau$           | [p.u.] | `tap`        | Off-nominal tap magnitude on the bus-1 side           |               | Default 1
$\theta$         | [rad]  | `phase`      | Phase-shift angle                                     |               | Default 0
$S^{\max}$       | [p.u.] | `Smax`       | Apparent power limit at each port                     |               | Default unbounded
$u$              | [-]    | `in_service` | Service status                                        |               | Default true

### Parameter Validation

- $R^2 + X^2 > 0$
- $\tau > 0$
- $S^{\max} > 0$

### Model Derived Parameters

```math
\begin{aligned}
g &= R / (R^2 + X^2) \\
b &= -X / (R^2 + X^2) \\
g_{11} &= u (-g / \tau^2 - G / 2 - G_\mathrm{mag}) \\
b_{11} &= u (-b / \tau^2 - B / 2 - B_\mathrm{mag}) \\
g_{12} &= u (g \cos\theta - b \sin\theta) / \tau \\
b_{12} &= u (b \cos\theta + g \sin\theta) / \tau \\
g_{21} &= u (g \cos\theta + b \sin\theta) / \tau \\
b_{21} &= u (b \cos\theta - g \sin\theta) / \tau \\
g_{22} &= u (-g - G / 2) \\
b_{22} &= u (-b - B / 2)
\end{aligned}
```

## Model Ports

Name   | Port | Description
-------|------|------------
`bus1` | Bus  | Bus-1 port, the tapped side
`bus2` | Bus  | Bus-2 port

## Model Variables

### Internal Variables

None.

### External Variables

Symbol            | Units  | Description             | Note
------------------|--------|-------------------------|-----
$V_{\mathrm{r}1}$ | [p.u.] | Bus-1 real voltage      | Owned by bus 1
$V_{\mathrm{i}1}$ | [p.u.] | Bus-1 imaginary voltage | Owned by bus 1
$V_{\mathrm{r}2}$ | [p.u.] | Bus-2 real voltage      | Owned by bus 2
$V_{\mathrm{i}2}$ | [p.u.] | Bus-2 imaginary voltage | Owned by bus 2

## Model Equations

```math
\begin{aligned}
I_{\mathrm{r}1} &= g_{11} V_{\mathrm{r}1} - b_{11} V_{\mathrm{i}1} + g_{12} V_{\mathrm{r}2} - b_{12} V_{\mathrm{i}2} \\
I_{\mathrm{i}1} &= b_{11} V_{\mathrm{r}1} + g_{11} V_{\mathrm{i}1} + b_{12} V_{\mathrm{r}2} + g_{12} V_{\mathrm{i}2} \\
I_{\mathrm{r}2} &= g_{21} V_{\mathrm{r}1} - b_{21} V_{\mathrm{i}1} + g_{22} V_{\mathrm{r}2} - b_{22} V_{\mathrm{i}2} \\
I_{\mathrm{i}2} &= b_{21} V_{\mathrm{r}1} + g_{21} V_{\mathrm{i}1} + b_{22} V_{\mathrm{r}2} + g_{22} V_{\mathrm{i}2} \\
P_1 &= V_{\mathrm{r}1} I_{\mathrm{r}1} + V_{\mathrm{i}1} I_{\mathrm{i}1} \\
Q_1 &= V_{\mathrm{i}1} I_{\mathrm{r}1} - V_{\mathrm{r}1} I_{\mathrm{i}1} \\
P_2 &= V_{\mathrm{r}2} I_{\mathrm{r}2} + V_{\mathrm{i}2} I_{\mathrm{i}2} \\
Q_2 &= V_{\mathrm{i}2} I_{\mathrm{r}2} - V_{\mathrm{r}2} I_{\mathrm{i}2}
\end{aligned}
```

### Objective

None.

### Internal Constraints

#### Equality

None.

#### Inequality

With a finite rating $S^{\max}$:

```math
\begin{aligned}
P_1^2 + Q_1^2 &\le (S^{\max})^2 \\
P_2^2 + Q_2^2 &\le (S^{\max})^2
\end{aligned}
```

### External Constraints

Power into each bus, added to its balance:

```math
\begin{aligned}
\Delta P_1^\mathrm{bus} &\mathrel{+}= P_1 \\
\Delta Q_1^\mathrm{bus} &\mathrel{+}= Q_1 \\
\Delta P_2^\mathrm{bus} &\mathrel{+}= P_2 \\
\Delta Q_2^\mathrm{bus} &\mathrel{+}= Q_2
\end{aligned}
```

## Initialization

None.
