# BusSignalVoltageOut

`BusSignalVoltageOut` is a bus with signal ports. It owns the same voltage variables
and current-balance residuals as `Bus`, publishes the voltage components on
output signals, and adds current injections received on input signals to its
residuals. Devices attached directly to the bus still add their currents after
the bus residual is evaluated.

## Notes

- Ports must be connected before `allocate()` is called. The output signals
  are linked to the bus voltage variables and their system indices in
  `allocate()`.
- `verify()` returns the number of connected ports that have no linked
  signal source.
- Current entering the bus has positive sign.

## Model Parameters

Same as `Bus`.

## Model Ports

Port | Direction | Units  | Description                                        | Note
-----|-----------|--------|----------------------------------------------------|-----
`vr` | out       | [p.u.] | Bus voltage, real component $V_r$                  |
`vi` | out       | [p.u.] | Bus voltage, imaginary component $V_i$             |
`ir` | in        | [p.u.] | Current injection, real component $I_r^{s}$        | Added to $f_0$
`ii` | in        | [p.u.] | Current injection, imaginary component $I_i^{s}$   | Added to $f_1$

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
$I_r^{s}$ | [p.u.] | Current injection on signal port `ir`  | Optional
$I_i^{s}$ | [p.u.] | Current injection on signal port `ii`  | Optional

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

An unconnected input port contributes zero.

## Initialization

Same as `Bus`.

## Monitors

Same as `Bus`.

## Testing

Unit tests in `tests/UnitTests/PhasorDynamics/BusSignalVoltageOutTests.hpp` cover
construction, output signal linking, residual evaluation with input signals,
verification of unlinked inputs, dependency-tracking derivatives, and (when
Enzyme is enabled) the sparse Jacobian entries.
