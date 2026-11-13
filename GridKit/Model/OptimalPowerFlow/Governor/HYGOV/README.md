# HYGOV

Steady state of the PhasorDynamics
[HYGOV](../../../PhasorDynamics/Governor/HYGOV/README.md), with its parameters
and symbols: the gate limits on the mechanical power, with the identity gate curve.

## Model Parameters

`Trate`, `Gmax`, `Gmin`, `At`, `Qnl`, and `Hdam` of HYGOV.

### Parameter Validation

The conditions of HYGOV on these parameters.

### Model Derived Parameters

$k_{\mathrm{base}}$ of HYGOV.

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
g = \dfrac{1}{\sqrt{H_{\mathrm{dam}}}} \left(\dfrac{k_{\mathrm{base}}}{A_\mathrm{t} H_{\mathrm{dam}}} P_\mathrm{m} + q_{\mathrm{NL}}\right)
```

### Objective

None.

### Internal Constraints

#### Equality

None.

#### Inequality

```math
G^{\min} \le g \le G^{\max}
```

### External Constraints

None.

## Initialization

None.
