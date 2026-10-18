# Load

ZIP demand: constant-power, constant-current, and constant-impedance parts.

## Model Parameters

Symbol | Units  | Data         | Description                                  | Typical Value | Note
-------|--------|--------------|----------------------------------------------|---------------|-----
$P$    | [p.u.] | `P`          | Constant-power active demand                 |               | Default 0
$Q$    | [p.u.] | `Q`          | Constant-power reactive demand               |               | Default 0
$P_I$  | [p.u.] | `Pi`         | Constant-current active demand at 1 p.u.     |               | Default 0
$Q_I$  | [p.u.] | `Qi`         | Constant-current reactive demand at 1 p.u.   |               | Default 0
$P_Z$  | [p.u.] | `Pz`         | Constant-impedance active demand at 1 p.u.   |               | Default 0
$Q_Z$  | [p.u.] | `Qz`         | Constant-impedance reactive demand at 1 p.u. |               | Default 0
$u$    | [-]    | `in_service` | Service status                               |               | Default true

### Parameter Validation

None.

### Model Derived Parameters

None.

## Model Ports

Name  | Port | Description
------|------|------------
`bus` | Bus  | Port bus

## Model Variables

### Internal Variables

None.

### External Variables

Symbol         | Units  | Description       | Note
---------------|--------|-------------------|-----
$V_\mathrm{r}$ | [p.u.] | Real voltage      | Owned by the bus
$V_\mathrm{i}$ | [p.u.] | Imaginary voltage | Owned by the bus

## Model Equations

### Objective

None.

### Internal Constraints

#### Equality

None.

#### Inequality

None.

### External Constraints

With $|V| = \sqrt{V_\mathrm{r}^2 + V_\mathrm{i}^2}$, power into the bus, added
to its balance:

```math
\begin{aligned}
\Delta P^\mathrm{bus} &\mathrel{+}= -u (P + P_I |V| + P_Z |V|^2) \\
\Delta Q^\mathrm{bus} &\mathrel{+}= -u (Q + Q_I |V| + Q_Z |V|^2)
\end{aligned}
```

## Initialization

None.
