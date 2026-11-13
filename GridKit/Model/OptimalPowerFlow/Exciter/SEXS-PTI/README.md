# SEXS-PTI

Steady state of the PhasorDynamics
[SEXS-PTI](../../../PhasorDynamics/Exciter/SEXS-PTI/README.md), with its parameters
and symbols: the limits on the field voltage.

## Model Parameters

`Efdmin` and `Efdmax` of SEXS-PTI.

### Parameter Validation

The conditions of SEXS-PTI on these parameters.

### Model Derived Parameters

None.

## Model Ports

Name  | Port   | Description
------|--------|--------------
`efd` | Signal | Field voltage

## Model Variables

### Internal Variables

None.

### External Variables

Symbol          | Units  | Description   | Note
----------------|--------|---------------|--------------------
$E_\mathrm{fd}$ | [p.u.] | Field voltage | Owned by the signal

## Model Equations

### Objective

None.

### Internal Constraints

#### Equality

None.

#### Inequality

```math
E_\mathrm{fd}^{\min} \le E_\mathrm{fd} \le E_\mathrm{fd}^{\max}
```

### External Constraints

None.

## Initialization

None.
