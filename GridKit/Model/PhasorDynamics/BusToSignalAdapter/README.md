# Bus-to-Signal Adapter

This component enables signals to send and receive bus variables.

## Model Parameters

None.

## Model Ports

Name  | Port   | Init  | Description
------|--------|-------|------------
`bus` | Bus    | Known | Bus whose variables are managed by the adapter
`ir`  | Input  | Known | Real current contribution to the bus
`ii`  | Input  | Known | Imaginary current contribution to the bus
`vr`  | Output | Known | Bus voltage, real component
`vi`  | Output | Known | Bus voltage, imaginary component

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description                              | Note
-------|--------|------------------------------------------|-----
$V_r$  | [p.u.] | Bus-voltage real component               | Bus-owned value published through `vr`
$V_i$  | [p.u.] | Bus-voltage imaginary component          | Bus-owned value published through `vi`
$I_r$  | [p.u.] | Real current contribution to the bus      | Read from the optional `ir` input
$I_i$  | [p.u.] | Imaginary current contribution to the bus | Read from the optional `ii` input

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

Each attached current input is added to the corresponding connected-bus
residual:

```math
\begin{aligned}
I_r^{\mathrm{bus}} &\leftarrow I_r^{\mathrm{bus}} + I_r \\
I_i^{\mathrm{bus}} &\leftarrow I_i^{\mathrm{bus}} + I_i.
\end{aligned}
```

## Initialization

None.

## Monitors

None.
