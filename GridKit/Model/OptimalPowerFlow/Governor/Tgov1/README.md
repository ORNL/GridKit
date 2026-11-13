# TGOV1

Steady state of the PhasorDynamics
[TGOV1](../../../PhasorDynamics/Governor/Tgov1/README.md), with its parameters
and symbols: the valve limits on the mechanical power.

## Model Parameters

`Trate`, `Pvmax`, and `Pvmin` of TGOV1.

### Parameter Validation

The conditions of TGOV1 on these parameters.

### Model Derived Parameters

None.

## Model Ports

Name    | Port   | Description
--------|--------|-----------------
`pmech` | Signal | Mechanical power

## Model Variables

### Internal Variables

None.

### External Variables

Symbol         | Units  | Description      | Note
---------------|--------|------------------|--------------------
$P_\mathrm{m}$ | [p.u.] | Mechanical power | Owned by the signal

## Model Equations

```math
P_\mathrm{v} = \dfrac{S_\mathrm{sys}}{T_\mathrm{rate}} P_\mathrm{m}
```

### Objective

None.

### Internal Constraints

#### Equality

None.

#### Inequality

```math
P_\mathrm{v}^{\min} \le P_\mathrm{v} \le P_\mathrm{v}^{\max}
```

### External Constraints

None.

## Initialization

None.
