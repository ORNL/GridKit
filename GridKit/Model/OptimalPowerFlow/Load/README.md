# Load

Constant power demand.

## Model Parameters

Symbol | Units  | Data         | Description                                       | Typical Value | Note
-------|--------|--------------|---------------------------------------------------|---------------|-----
$P$    | [p.u.] | `P`          | Active power demand                               |               | Default 0
$Q$    | [p.u.] | `Q`          | Reactive power demand                             |               | Default 0
$u$    | [-]    | `in_service` | Service status, 1 in service and 0 out of service |               | Default 1

### Parameter Validation

- $P$ and $Q$ are finite

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
\Delta P^\mathrm{bus} &\mathrel{+}= -u P \\
\Delta Q^\mathrm{bus} &\mathrel{+}= -u Q
\end{aligned}
```

## Initialization

None.
