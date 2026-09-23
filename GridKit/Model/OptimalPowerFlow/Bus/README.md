# Bus

Bus voltage in Cartesian coordinates and the power balance of the bus. Devices
add the power they inject into the bus to its balance rows.

## Model Parameters

Symbol     | Units  | JSON   | Description                   | Typical Value | Note
-----------|--------|--------|-------------------------------|---------------|-----
$V^{\min}$ | [p.u.] | `Vmin` | Voltage magnitude lower limit | 0.95          | Default 0
$V^{\max}$ | [p.u.] | `Vmax` | Voltage magnitude upper limit | 1.05          | Default unbounded

A `BusInfinite` in the case makes the bus infinite.

### Parameter Validation

- $V^{\min} \ge 0$
- $V^{\max} \ge V^{\min}$

### Model Derived Parameters

None.

## Model Ports

None. Devices connect to the bus by its number.

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
0 &= \sum_{c \to k} P_{c,k} \\
0 &= \sum_{c \to k} Q_{c,k} \\
\left(V^{\min}\right)^2 &\le V_\mathrm{r}^2 + V_\mathrm{i}^2 \le \left(V^{\max}\right)^2 \\
0 &= V_\mathrm{r}^\mathrm{state} V_\mathrm{i} - V_\mathrm{i}^\mathrm{state} V_\mathrm{r}
\end{aligned}
```

The bus contributes zero to the balance rows. The magnitude row exists with a
limit, and the reference row only at the reference bus, where it keeps the
state angle. An infinite bus keeps its state voltage and has no rows.

## Initialization

The state entry `bus_id_N` gives the starting voltage, $V_\mathrm{r}^\mathrm{state}$
and $V_\mathrm{i}^\mathrm{state}$. Without an entry, the bus starts at
$V_\mathrm{r} = 1$ and $V_\mathrm{i} = 0$.
