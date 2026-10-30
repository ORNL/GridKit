# GasCompressor

GasCompressor holds the bus-2 pressure at a fixed ratio to the bus-1 pressure
and moves whatever mass flow that takes from bus 1 to bus 2.

## Notes

- As in the GasTranSim ratio control, the compressor has no volume and uses no
  gas.
- The ratio is enforced through its time derivative. A direct constraint
  between two buses that store gas would make the system index 2.

## Model Parameters

Symbol   | Units | JSON    | Description       | Typical Value | Note
---------|-------|---------|-------------------|---------------|-----
$\alpha$ | [-]   | `ratio` | Compression ratio | Required      |

### Parameter Validation

A valid GasCompressor parameter set must satisfy the following conditions:

```math
0 < \alpha < \infty
```

At least one terminal bus must be a GasBus.

### Model Derived Parameters

None.

## Model Ports

Name   | Port   | Init  | Description
-------|--------|-------|------------------------
`bus1` | GasBus | Known | Required suction bus
`bus2` | GasBus | Known | Required discharge bus

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description | Note
-------|--------|-------------|-----------------------------
$q$    | [kg/s] | Mass flow   | Positive from bus 1 to bus 2

### External Variables

#### Differential

Symbol | Units | Description        | Note
-------|-------|--------------------|--------------------
$p_1$  | [Pa]  | Suction pressure   | Owned by bus object
$p_2$  | [Pa]  | Discharge pressure | Owned by bus object

#### Algebraic

None.

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

```math
0 = \alpha\dot p_1 - \dot p_2
```

### External Equations

```math
\begin{aligned}
  \Delta q^\mathrm{bus}_1 &\mathrel{+}= -q \\
  \Delta q^\mathrm{bus}_2 &\mathrel{+}= q
\end{aligned}
```

## Initialization

### Input Initialization

```math
p_1, p_2 \leftarrow \text{terminal-bus pressures}
```

### Internal Initialization

Initialization requires $p_1, p_2 > 0$ and $p_2 = \alpha p_1$. The flow starts
at zero, and the consistent initial condition solve sets it to the flow that
balances the buses.

```math
q \leftarrow 0
```

## Monitors

Monitor | Units  | Description | Note
--------|--------|-------------|-----------------------------
`q`     | [kg/s] | Mass flow   | Positive from bus 1 to bus 2
