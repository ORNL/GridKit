# GasBranch

GasBranch is a gas pipe between two gas buses, divided into $N$ equal segments
along its length. Terminal mass flows are oriented entering the adjacent buses.

## Notes

- Each terminal holds half a segment of gas storage, like the split line shunt
  of [Branch](../../Branch/README.md).
- Flow is isothermal. As in GasTranSim, convective acceleration and elevation
  change are neglected.

## Block Diagram

![GasBranch model block diagram](../../../../../docs/Figures/PhasorDynamics/GasBranch/diagram.png)

Figure 1: GasBranch model.

## Model Parameters

Symbol            | Units | JSON     | Description               | Typical Value | Note
------------------|-------|----------|---------------------------|---------------|-----
$L$               | [m]   | `L`      | Pipe length               | Required      |
$D$               | [m]   | `D`      | Pipe diameter             | Required      |
$\lambda$         | [-]   | `lambda` | Darcy friction factor     | Required      |
$a$               | [m/s] | `a`      | Isothermal speed of sound | Required      |
$\Delta x^{\max}$ | [m]   | `dxmax`  | Maximum segment length    | 1000.0        |

### Parameter Validation

A valid GasBranch parameter set must satisfy the following conditions:

```math
0 < L, D, \lambda, a, \Delta x^{\max} < \infty
```

### Model Derived Parameters

```math
\begin{aligned}
  N &= \left\lceil \dfrac{L}{\Delta x^{\max}} \right\rceil \\
  \Delta x &= \dfrac{L}{N} \\
  A &= \dfrac{\pi D^2}{4}
\end{aligned}
```

## Model Ports

Name   | Port   | Init  | Description
-------|--------|-------|------------------------
`bus1` | GasBus | Known | Required bus-1 terminal
`bus2` | GasBus | Known | Required bus-2 terminal

## Model Variables

### Internal Variables

#### Differential

Symbol   | Units       | Description       | Note
---------|-------------|-------------------|-----------------------------------------------------
$p_j$    | [Pa]        | Node pressure     | $j \in \{1,\ldots,N-1\}$
$\phi_k$ | [kg/(m²·s)] | Segment mass flux | $k \in \{1,\ldots,N\}$; positive from bus 1 to bus 2

#### Algebraic

None.

### External Variables

#### Differential

Symbol | Units | Description    | Note
-------|-------|----------------|--------------------
$p_0$  | [Pa]  | Bus-1 pressure | Owned by bus object
$p_N$  | [Pa]  | Bus-2 pressure | Owned by bus object

#### Algebraic

None.

## Model Equations

Node $j$ is at distance $j\Delta x$ from bus 1, and segment $k$ joins nodes
$k-1$ and $k$. Define the node and mean segment densities:

```math
\begin{aligned}
  \rho_j &= \dfrac{p_j}{a^2} \\
  \bar\rho_k &= \dfrac{\rho_{k-1} + \rho_k}{2}
\end{aligned}
```

### Internal Equations

#### Differential

```math
\begin{aligned}
  0 &= -\dot p_j + \dfrac{a^2}{\Delta x}(\phi_j - \phi_{j+1}), && j \in \{1,\ldots,N-1\} \\
  0 &= -\dot\phi_k + \dfrac{p_{k-1} - p_k}{\Delta x} - \dfrac{\lambda}{2D}\dfrac{\phi_k|\phi_k|}{\bar\rho_k}, && k \in \{1,\ldots,N\}
\end{aligned}
```

#### Algebraic

None.

### External Equations

```math
\begin{aligned}
  \Delta q^\mathrm{bus}_1 &\mathrel{+}= -A\phi_1 - \dfrac{A\Delta x}{2a^2}\dot p_0 \\
  \Delta q^\mathrm{bus}_2 &\mathrel{+}= A\phi_N - \dfrac{A\Delta x}{2a^2}\dot p_N
\end{aligned}
```

## Initialization

### Input Initialization

```math
p_0, p_N \leftarrow \text{terminal-bus pressures}
```

### Internal Initialization

Initialization requires $p_0, p_N > 0$. All internal derivatives initialize to
zero, and the values below are the steady flow between the terminal pressures:

```math
\begin{aligned}
  \phi_k &\leftarrow \text{sign}(p_0^2 - p_N^2)\sqrt{\dfrac{D\,|p_0^2 - p_N^2|}{\lambda a^2 L}} \\
  p_j &\leftarrow \sqrt{p_0^2 - \dfrac{j}{N}(p_0^2 - p_N^2)}
\end{aligned}
```

## Monitors

Monitor    | Units  | Description               | Note
-----------|--------|---------------------------|---------------------------------------------------------
`q1`       | [kg/s] | Terminal mass flow, bus 1 | $\Delta q^\mathrm{bus}_1$; oriented entering bus 1
`q2`       | [kg/s] | Terminal mass flow, bus 2 | $\Delta q^\mathrm{bus}_2$; oriented entering bus 2
`linepack` | [kg]   | Gas mass in the pipe      | $A\Delta x\,[(\rho_0+\rho_N)/2+\sum_{j=1}^{N-1}\rho_j]$
