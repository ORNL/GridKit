# BusSignalVoltageOut

`BusSignalVoltageOut` is a bus with signal ports. It owns the same voltage variables
and current-balance residuals as `Bus`, publishes the voltage components on
output signals, and adds current injections received on input signals to its
residuals. Devices attached directly to the bus still add their currents after
the bus residual is evaluated.

## Notes

- Signal ports must be connected before `allocate()` is called. The signal outlets
  are linked to the bus voltage variables and their system indices in
  `allocate()`.
- Both current inlets are mandatory. `verify()` logs each problem and
  throws if a current inlet is not connected or not linked, or if a
  connected outlet is not linked. No default current is ever used.
- Current entering the bus has positive sign.

## Model Parameters

Same as `Bus`.

## Model Ports

Signal | Direction | Units  | Description                                        | Note
-----|-----------|--------|----------------------------------------------------|-----
`vr` | out       | [p.u.] | Bus voltage, real component $V_r$                  |
`vi` | out       | [p.u.] | Bus voltage, imaginary component $V_i$             |
`ir` | in        | [p.u.] | Current injection, real component $I_r^{s}$        | Required, sets $f_0$
`ii` | in        | [p.u.] | Current injection, imaginary component $I_i^{s}$   | Required, sets $f_1$

## Model Variables

### Internal Variables

#### Algebraic

Symbol | Units  | Description                      | Note
-------|--------|----------------------------------|-----
$V_r$  | [p.u.] | Bus voltage, real component      |
$V_i$  | [p.u.] | Bus voltage, imaginary component |

### External Variables

#### Algebraic

Symbol    | Units  | Description                            | Note
----------|--------|----------------------------------------|-----
$I_r^{s}$ | [p.u.] | Current injection on signal inlet `ir`  |
$I_i^{s}$ | [p.u.] | Current injection on signal inlet `ii`  |

## Model Equations

### Internal Equations

#### Algebraic

Let $\mathcal{D}$ denote the set of devices attached directly to the bus.

```math
\begin{aligned}
0 &= I_r^{s} + \sum_{d \in \mathcal{D}} I_{r,d} \\
0 &= I_i^{s} + \sum_{d \in \mathcal{D}} I_{i,d}
\end{aligned}
```

Both signal inlets must be connected; see `verify()`.

## Initialization

Same as `Bus`.

## Monitors

Same as `Bus`.

## Testing

Unit tests in `tests/UnitTests/PhasorDynamics/BusSignalVoltageOutTests.hpp` cover
construction, output signal linking, residual evaluation with input signals,
verification of unlinked inputs, dependency-tracking derivatives, and (when
Enzyme is enabled) the sparse Jacobian entries.
