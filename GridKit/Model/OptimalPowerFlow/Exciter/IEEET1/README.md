# IEEET1

Steady state of the PhasorDynamics
[IEEET1](../../../PhasorDynamics/Exciter/IEEET1/README.md), with its parameters
and symbols: the regulator limits on the field voltage.

## Model Parameters

`Ke`, `Vrmin`, `Vrmax`, `E1`, `E2`, `Se1`, and `Se2` of IEEET1.

### Parameter Validation

The conditions of IEEET1 on these parameters.

### Model Derived Parameters

$S_A$, $S_B$, and $K_E^{\mathrm{eff}}$ of IEEET1.

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

Smooth functions: [`qramp`](../../../../CommonMath.md#quadratic-ramp).

```math
V_R = K_E^{\mathrm{eff}} E_\mathrm{fd} + S_B \text{qramp}(E_\mathrm{fd} - S_A)
```

### Objective

None.

### Internal Constraints

#### Equality

None.

#### Inequality

```math
V_R^{\min} \le V_R \le V_R^{\max}
```

### External Constraints

None.

## Initialization

None.
