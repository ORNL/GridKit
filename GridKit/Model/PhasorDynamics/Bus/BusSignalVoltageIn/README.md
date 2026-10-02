# BusSignalVoltageIn

`BusSignalVoltageIn` is the mirror image of `BusSignalVoltageOut`. Its
voltage is read directly from input signals, and the sum of current
injections from components attached to it is published on output signals.
The bus stores no voltage of its own and has no unknowns and no equations.

## Notes

- Ports must be connected before `allocate()` is called. The output signals
  are linked to the current sums in `allocate()`.
- Both voltage inlets are mandatory. `verify()` logs each problem and
  throws if a voltage inlet is not connected or not linked, or if a
  connected outlet is not linked. Reading the voltage through an unlinked inlet
  throws. No default voltage is ever used.
- The current sums are complete only after all attached components have
  evaluated their residuals. Consumers of `ir` and `ii` must be evaluated
  after them.
- Current entering the bus has positive sign.

## Model Parameters

Same as `Bus`.

## Model Ports

Port | Direction | Units  | Description                                           | Note
-----|-----------|--------|-------------------------------------------------------|-----
`vr` | in        | [p.u.] | Bus voltage, real component $V_r$                     | Required
`vi` | in        | [p.u.] | Bus voltage, imaginary component $V_i$                | Required
`ir` | out       | [p.u.] | Sum of real current injections $I_r$                  |
`ii` | out       | [p.u.] | Sum of imaginary current injections $I_i$             |

## Model Variables

### Internal Variables

None.

### External Variables

#### Algebraic

Symbol | Units  | Description                                   | Note
-------|--------|-----------------------------------------------|-----
$V_r$  | [p.u.] | Bus voltage, real component from port `vr`      |
$V_i$  | [p.u.] | Bus voltage, imaginary component from port `vi` |

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
`Vr()` and `Vi()` return the connected signal value by reference.

## Initialization

Current sums are set to zero. The voltage is owned by the signal sources and
is not initialized by the bus; the initial voltage in bus data is ignored.

## Monitors

Same as `Bus`.

## Testing

Unit tests in `tests/UnitTests/PhasorDynamics/BusSignalVoltageInTests.hpp`
cover construction, voltage inputs, current output accumulation and reset,
verification of unlinked inputs, and dependency-tracking derivatives.
