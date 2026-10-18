# Generator

Dispatchable active and reactive power injection with a quadratic cost.

## Model Parameters

Symbol     | Units  | Data         | Description                                       | Typical Value | Note
-----------|--------|--------------|---------------------------------------------------|---------------|-----
$P^{\min}$ | [p.u.] | `Pmin`       | Active power lower limit                          | 0             | Default unbounded
$P^{\max}$ | [p.u.] | `Pmax`       | Active power upper limit                          | 1             | Default unbounded
$Q^{\min}$ | [p.u.] | `Qmin`       | Reactive power lower limit                        | -0.5          | Default unbounded
$Q^{\max}$ | [p.u.] | `Qmax`       | Reactive power upper limit                        | 0.5           | Default unbounded
$c_0$      | [\$/h] | `c0`         | Constant cost coefficient                         |               | Default 0
$c_1$      | [\$/h] | `c1`         | Linear cost coefficient                           | 1000          | Default 0
$c_2$      | [\$/h] | `c2`         | Quadratic cost coefficient                        | 100           | Default 0
$u$        | [-]    | `in_service` | Service status, 1 in service and 0 out of service |               | Default true

### Parameter Validation

- $P^{\min} \le P^{\max}$
- $Q^{\min} \le Q^{\max}$
- $c_0$, $c_1$, and $c_2$ are finite

### Model Derived Parameters

None.

## Model Ports

Name  | Port | Description
------|------|------------
`bus` | Bus  | Port bus

## Model Variables

### Internal Variables

Symbol | Units  | Description                 | Note
-------|--------|-----------------------------|-----
$P$    | [p.u.] | Active power into the bus   |
$Q$    | [p.u.] | Reactive power into the bus |

### External Variables

None.

## Model Equations

### Objective

```math
\begin{aligned}
f &\mathrel{+}= u (c_0 + c_1 P + c_2 P^2)
\end{aligned}
```

### Internal Constraints

#### Equality

Out of service, imposed as fixed variable bounds:

```math
\begin{aligned}
0 &= P \\
0 &= Q
\end{aligned}
```

#### Inequality

In service, imposed as variable bounds:

```math
\begin{aligned}
P^{\min} &\le P \le P^{\max} \\
Q^{\min} &\le Q \le Q^{\max}
\end{aligned}
```

### External Constraints

Power into the bus, added to its balance:

```math
\begin{aligned}
\Delta P^\mathrm{bus} &\mathrel{+}= P \\
\Delta Q^\mathrm{bus} &\mathrel{+}= Q
\end{aligned}
```

## Initialization

```math
\begin{aligned}
P &\leftarrow 0 \\
Q &\leftarrow 0
\end{aligned}
```
