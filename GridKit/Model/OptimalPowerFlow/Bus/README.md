# Bus

Bus voltage in Cartesian coordinates, power balance, and shunt admittance to
ground. A port to the bus reads its voltage and contributes to its balance.

## Model Parameters

Symbol     | Units  | Data        | Description                                  | Typical Value | Note
-----------|--------|-------------|----------------------------------------------|---------------|-----
$V^{\min}$ | [p.u.] | `Vmin`      | Voltage magnitude lower limit                | 0.95          | Default 0
$V^{\max}$ | [p.u.] | `Vmax`      | Voltage magnitude upper limit                | 1.05          | Default unbounded
$G$        | [p.u.] | `G`         | Shunt conductance to ground                  |               | Default 0
$B$        | [p.u.] | `B`         | Shunt susceptance to ground                  |               | Default 0
           | [-]    | `reference` | Angle reference, with zero imaginary voltage |               | Default false

### Parameter Validation

- $V^{\min} \ge 0$
- $V^{\max} \ge V^{\min}$
- $G$ and $B$ are finite

### Model Derived Parameters

None.

## Model Ports

None.

## Model Variables

### Internal Variables

Symbol         | Units  | Description       | Note
---------------|--------|-------------------|-----
$V_\mathrm{r}$ | [p.u.] | Real voltage      |
$V_\mathrm{i}$ | [p.u.] | Imaginary voltage |

### External Variables

None.

## Model Equations

### Objective

None.

### Internal Constraints

#### Equality

Power balance, the power into the bus from its shunt plus the power that the
devices at its ports add:

```math
\begin{aligned}
0 &= -G (V_\mathrm{r}^2 + V_\mathrm{i}^2) + \Delta P^\mathrm{bus} \\
0 &= B (V_\mathrm{r}^2 + V_\mathrm{i}^2) + \Delta Q^\mathrm{bus}
\end{aligned}
```

The reference condition is imposed as a variable bound:

```math
0 = V_\mathrm{i} \qquad \text{reference bus}
```

#### Inequality

```math
(V^{\min})^2 \le V_\mathrm{r}^2 + V_\mathrm{i}^2 \le (V^{\max})^2
```

### External Constraints

None.

## Initialization

```math
\begin{aligned}
V_\mathrm{r} &\leftarrow 1 \\
V_\mathrm{i} &\leftarrow 0
\end{aligned}
```
