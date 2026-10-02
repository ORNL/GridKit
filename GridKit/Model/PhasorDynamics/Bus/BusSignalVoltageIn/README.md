# BusSignalVoltageIn

`BusSignalVoltageIn` is the mirror image of `BusSignalVoltageOut`. Its
voltage is read directly from input signals, and the sum of current
injections from components attached to it is published on output signals.
The bus stores no voltage of its own and has no unknowns and no equations.

## Notes

- Ports must be connected before `allocate()` is called. The output signals
  are linked to the current sums in `allocate()`.
- `verify()` returns the number of connected ports that have no linked
  signal source.
- The current sums are complete only after all attached components have
  evaluated their residuals. Consumers of `ir` and `ii` must be evaluated
  after them.
- Current entering the bus has positive sign.

## Model Parameters

Same as `Bus`.

## Model Ports

Port | Direction | Units  | Description                                           | Note
-----|-----------|--------|-------------------------------------------------------|-----
`vr` | in        | [p.u.] | Bus voltage, real component $V_r$                     | Sets $V_r$
`vi` | in        | [p.u.] | Bus voltage, imaginary component $V_i$                | Sets $V_i$
`ir` | out       | [p.u.] | Sum of real current injections $I_r$                  |
`ii` | out       | [p.u.] | Sum of imaginary current injections $I_i$             |

## Model Variables

### Internal Variables

None.

### External Variables

#### Algebraic

Symbol | Units  | Description                                   | Note
-------|--------|-----------------------------------------------|-----
$V_r$  | [p.u.] | Bus voltage, real component from port `vr`    | Optional
$V_i$  | [p.u.] | Bus voltage, imaginary component from port `vi` | Optional

## Model Equations

### Internal Equations

None. The published outputs are

```math
\begin{aligned}
I_r &= \sum_{d \in \mathcal{D}} I_{r,d} \\
I_i &= \sum_{d \in \mathcal{D}} I_{i,d}
\end{aligned}
```

where $\mathcal{D}$ is the set of devices attached directly to the bus.
`Vr()` and `Vi()` return the connected signal value; an unconnected input
port returns the initial value from bus data instead.

## Initialization

Current sums are set to zero. The voltage is owned by the signal sources and
is not initialized by the bus.

## Monitors

Same as `Bus`.

## Testing

Unit tests in `tests/UnitTests/PhasorDynamics/BusSignalVoltageInTests.hpp`
cover construction, voltage inputs, current output accumulation and reset,
verification of unlinked inputs, and dependency-tracking derivatives.
