# Shunt

Constant admittance $Y = G + jB$ to ground. The `OptimalDispatch` application
makes one for every LoadZ device in the case.

## Model Parameters

Symbol | Units  | JSON | Description       | Typical Value | Note
-------|--------|------|-------------------|---------------|-----
$G$    | [p.u.] | `G`  | Shunt conductance |               | $R / (R^2 + X^2)$ of the LoadZ
$B$    | [p.u.] | `B`  | Shunt susceptance |               | $-X / (R^2 + X^2)$ of the LoadZ

### Parameter Validation

- $G$ and $B$ are finite

### Model Derived Parameters

None.

## Model Ports

Name  | Port | Init  | Description
------|------|-------|------------
`bus` | Bus  | Known | Terminal bus

## Model Variables

### Internal Variables

None.

### External Variables

Symbol           | Units  | Description                  | Note
-----------------|--------|------------------------------|-----
$V_\mathrm{r}$   | [p.u.] | Terminal real voltage        | Owned by the bus
$V_\mathrm{i}$   | [p.u.] | Terminal imaginary voltage   | Owned by the bus

## Model Equations

### Objective

None.

### External Constraints

Power into the bus is $-Y^* \lvert V \rvert^2$:

```math
\begin{aligned}
\Delta P^\mathrm{bus} &\mathrel{+}= -a G \left(V_\mathrm{r}^2 + V_\mathrm{i}^2\right) \\
\Delta Q^\mathrm{bus} &\mathrel{+}= a B \left(V_\mathrm{r}^2 + V_\mathrm{i}^2\right)
\end{aligned}
```

## Initialization

$a \leftarrow 0$ if the state marks the device `online` as false, and $a
\leftarrow 1$ otherwise.
