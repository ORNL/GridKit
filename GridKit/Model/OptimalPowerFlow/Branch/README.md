# Branch

Line or off-nominal transformer between two buses, with the admittance of the
PhasorDynamics [Branch](../../PhasorDynamics/Branch/README.md).

## Model Parameters

Symbol                | Units  | JSON    | Description                                               | Typical Value | Note
----------------------|--------|---------|-----------------------------------------------------------|---------------|-----
$R$                   | [p.u.] | `R`     | Series resistance                                         |               | From the case
$X$                   | [p.u.] | `X`     | Series reactance                                          |               | From the case
$G$                   | [p.u.] | `G`     | Total line shunt conductance, split between the terminals |               | From the case
$B$                   | [p.u.] | `B`     | Total line shunt susceptance, split between the terminals |               | From the case
$G_\mathrm{mag}$      | [p.u.] | `Gmag`  | Magnetizing conductance at bus 1                          | 0             | From the case
$B_\mathrm{mag}$      | [p.u.] | `Bmag`  | Magnetizing susceptance at bus 1                          | 0             | From the case
$\tau$                | [p.u.] | `tap`   | Off-nominal tap magnitude on the bus-1 side               | 1             | From the case
$\phi$                | [rad]  | `phase` | Phase-shift angle                                         | 0             | From the case
$S^{\max}$            | [p.u.] | `Smax`  | Apparent power limit at each terminal                     |               | MATPOWER `RATE_A`, default unbounded

### Parameter Validation

- All parameters are finite, except $S^{\max}$
- $R^2 + X^2 > 0$
- $\tau > 0$
- $S^{\max} > 0$

### Model Derived Parameters

The state settings `tap` and `phase` replace $\tau$ and $\phi$. The status $s$
is zero if the state marks the branch `open` and one otherwise.

```math
\begin{aligned}
g &= \dfrac{R}{R^2 + X^2}, \quad b = -\dfrac{X}{R^2 + X^2} \\
g_{11} &= s \left(-g / \tau^2 - G/2 - G_\mathrm{mag}\right), \quad b_{11} = s \left(-b / \tau^2 - B/2 - B_\mathrm{mag}\right) \\
g_{12} &= s (g \cos\phi - b \sin\phi) / \tau, \quad b_{12} = s (b \cos\phi + g \sin\phi) / \tau \\
g_{21} &= s (g \cos\phi + b \sin\phi) / \tau, \quad b_{21} = s (b \cos\phi - g \sin\phi) / \tau \\
g_{22} &= s (-g - G/2), \quad b_{22} = s (-b - B/2)
\end{aligned}
```

## Model Ports

Name   | Port | Init  | Description
-------|------|-------|------------
`bus1` | Bus  | Known | Bus-1 terminal, the tapped side
`bus2` | Bus  | Known | Bus-2 terminal

## Model Variables

### Internal Variables

None.

### External Variables

Symbol              | Units  | Description           | Note
--------------------|--------|-----------------------|-----
$V_{\mathrm{r}1}$   | [p.u.] | Bus-1 real voltage      | Owned by bus 1
$V_{\mathrm{i}1}$   | [p.u.] | Bus-1 imaginary voltage | Owned by bus 1
$V_{\mathrm{r}2}$   | [p.u.] | Bus-2 real voltage      | Owned by bus 2
$V_{\mathrm{i}2}$   | [p.u.] | Bus-2 imaginary voltage | Owned by bus 2

## Model Equations

Terminal currents into the buses, and the power into bus $k$:

```math
\begin{aligned}
I_{\mathrm{r}1} &= g_{11} V_{\mathrm{r}1} - b_{11} V_{\mathrm{i}1} + g_{12} V_{\mathrm{r}2} - b_{12} V_{\mathrm{i}2} \\
I_{\mathrm{i}1} &= b_{11} V_{\mathrm{r}1} + g_{11} V_{\mathrm{i}1} + b_{12} V_{\mathrm{r}2} + g_{12} V_{\mathrm{i}2} \\
I_{\mathrm{r}2} &= g_{21} V_{\mathrm{r}1} - b_{21} V_{\mathrm{i}1} + g_{22} V_{\mathrm{r}2} - b_{22} V_{\mathrm{i}2} \\
I_{\mathrm{i}2} &= b_{21} V_{\mathrm{r}1} + g_{21} V_{\mathrm{i}1} + b_{22} V_{\mathrm{r}2} + g_{22} V_{\mathrm{i}2} \\
P_k &= V_{\mathrm{r}k} I_{\mathrm{r}k} + V_{\mathrm{i}k} I_{\mathrm{i}k}, \quad Q_k = V_{\mathrm{i}k} I_{\mathrm{r}k} - V_{\mathrm{r}k} I_{\mathrm{i}k}
\end{aligned}
```

### Objective

None.

### Internal Constraints

A row exists only when $S^{\max}$ is set:

```math
P_k^2 + Q_k^2 \le \left(S^{\max}\right)^2, \quad k = 1, 2
```

### External Constraints

```math
\begin{aligned}
\Delta P_k^\mathrm{bus} &\mathrel{+}= P_k \\
\Delta Q_k^\mathrm{bus} &\mathrel{+}= Q_k
\end{aligned}
```

## Initialization

The state settings `open`, `tap`, and `phase` set the derived admittance.
