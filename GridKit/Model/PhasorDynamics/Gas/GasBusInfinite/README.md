# GasBusInfinite

A fixed-pressure gas bus supplies whatever mass flow the connected devices
draw. It has no variables.

## Notes

- A `pressure` event changes the pressure during a simulation.
- Mass flows added by connected devices are discarded.

## Model Parameters

None.

### Parameter Validation

None.

### Model Derived Parameters

None.

## Model Ports

None.

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

None.

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

None.

## Initialization

### Internal Initialization

The pressure must be positive.

```math
p \leftarrow \text{initial bus pressure}
```

## Monitors

Monitor | Units | Description  | Note
--------|-------|--------------|-----
`p`     | [Pa]  | Bus pressure |
