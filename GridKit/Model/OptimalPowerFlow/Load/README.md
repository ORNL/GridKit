# Load

Fixed power demand given by the state. The `OptimalDispatch` application makes
one for every LoadZIP device in the case.

## Model Parameters

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

```math
\begin{aligned}
\Delta P^\mathrm{bus} &\mathrel{+}= P^\mathrm{state} \\
\Delta Q^\mathrm{bus} &\mathrel{+}= Q^\mathrm{state}
\end{aligned}
```

## Initialization

The state terminal current and bus voltage give the injection:

```math
P^\mathrm{state} + jQ^\mathrm{state} \leftarrow V_\mathrm{bus} I^*
```

An offline load injects zero. A load without a state current is an error.
