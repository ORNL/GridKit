# Shunt

Constant admittance $Y = G + jB$ to ground.

## Model Parameters

Symbol | Units  | Data         | Description                                       | Typical Value | Note
-------|--------|--------------|---------------------------------------------------|---------------|-----
$G$    | [p.u.] | `G`          | Shunt conductance                                 |               | Default 0
$B$    | [p.u.] | `B`          | Shunt susceptance                                 |               | Default 0
$u$    | [-]    | `in_service` | Service status, 1 in service and 0 out of service |               | Default 1

### Parameter Validation

- $G$ and $B$ are finite

### Model Derived Parameters

None.

## Model Ports

Name  | Port | Description
------|------|------------
`bus` | Bus  | Terminal bus

## Model Variables

### Internal Variables

None.

### External Variables

Symbol         | Units  | Description                | Note
---------------|--------|----------------------------|-----
$V_\mathrm{r}$ | [p.u.] | Terminal real voltage      | Owned by the bus
$V_\mathrm{i}$ | [p.u.] | Terminal imaginary voltage | Owned by the bus

## Model Equations

### Objective

None.

### Internal Constraints

None.

### External Constraints

```math
\begin{aligned}
\Delta P^\mathrm{bus} &\mathrel{+}= -u G (V_\mathrm{r}^2 + V_\mathrm{i}^2) \\
\Delta Q^\mathrm{bus} &\mathrel{+}= u B (V_\mathrm{r}^2 + V_\mathrm{i}^2)
\end{aligned}
```

## Initialization

None.
