# ConstantSignalSource

Zero-state component that publishes constant real and imaginary scalar values
on two output signals.

## Model Parameters

Symbol | Units       | JSON | Description                     | Default
-------|-------------|------|---------------------------------|--------
$S_r$  | unspecified | `Sr` | Constant real output value      | 0.0
$S_i$  | unspecified | `Si` | Constant imaginary output value | 0.0

### Parameter Validation

None.

### Model Derived Parameters

None.

## Model Ports

Name | Port   | Init  | Description
-----|--------|-------|------------
`sr` | Output | Known | Constant real component $S_r$
`si` | Output | Known | Constant imaginary component $S_i$

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

None.

## Monitors

None.
