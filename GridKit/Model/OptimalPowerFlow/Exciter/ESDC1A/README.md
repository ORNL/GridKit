# ESDC1A

Steady state of the PhasorDynamics
[ESDC1A](../../../PhasorDynamics/Exciter/ESDC1A/README.md), with its parameters
and symbols: the regulator and field-voltage limits on the field voltage.

## Model Parameters

`Ke`, `Vrmax`, `Vrmin`, `E1`, `Se1`, `E2`, `Se2`, and `exclim` of ESDC1A.

### Parameter Validation

The conditions of ESDC1A on these parameters.

### Model Derived Parameters

$S_A$, $S_B$, and $K_E^{\mathrm{eff}}$ of ESDC1A.

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
\begin{aligned}
V_R^{\min} &\le V_R \le V_R^{\max} \\
0 &\le s_{\mathrm{lim}} E_\mathrm{fd}
\end{aligned}
```

### External Constraints

None.

## Initialization

None.
