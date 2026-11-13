# GASTPTI

Steady state of the PhasorDynamics
[GASTPTI](../../../PhasorDynamics/Governor/GASTPTI/README.md), with its parameters
and symbols: the valve and temperature limits on the mechanical power.

## Model Parameters

`Trate`, `Vmax`, `Vmin`, and `At` of GASTPTI.

### Parameter Validation

The conditions of GASTPTI on these parameters.

### Model Derived Parameters

$k_{\mathrm{base}}$ of GASTPTI.

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
x_F = k_{\mathrm{base}} P_\mathrm{m}
```

### Objective

None.

### Internal Constraints

#### Equality

None.

#### Inequality

```math
V^{\min} \le x_F \le \min(V^{\max}, A_T)
```

### External Constraints

None.

## Initialization

None.
