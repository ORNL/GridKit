# Constant signal source

This component emits a constant complex value on two output ports (real and
imaginary).

## Model Parameters

The complex-value parameter is intentionally ambiguous, because it may be
applied in different contexts (for different input variables).

Symbol      | Units   | Description                     | Note
------------|---------|---------------------------------| ------
$Sr$  | unspecified | Real component  |
$Si$  | unspecified | Imaginary component  | 

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
