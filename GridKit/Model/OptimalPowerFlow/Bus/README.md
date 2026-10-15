# Bus

Bus voltage in Cartesian coordinates and power balance.

## Model Parameters

Symbol     | Units  | Data        | Description                                  | Typical Value | Note
-----------|--------|-------------|----------------------------------------------|---------------|-----
$V^{\min}$ | [p.u.] | `Vmin`      | Voltage magnitude lower limit                | 0.95          | Default 0
$V^{\max}$ | [p.u.] | `Vmax`      | Voltage magnitude upper limit                | 1.05          | Default unbounded
           | [-]    | `reference` | Angle reference, with zero imaginary voltage |               | Default false

### Parameter Validation

- $V^{\min} \ge 0$
- $V^{\max} \ge V^{\min}$

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

```math
\begin{aligned}
0 &= \Delta P^\mathrm{bus} \\
0 &= \Delta Q^\mathrm{bus} \\
(V^{\min})^2 &\le V_\mathrm{r}^2 + V_\mathrm{i}^2 \le (V^{\max})^2 \\
0 &= V_\mathrm{i} && \text{reference bus}
\end{aligned}
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
