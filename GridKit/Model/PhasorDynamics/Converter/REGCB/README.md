# **Renewable Energy Generator/Converter Model (REGCB)**

REGCB is a WECC renewable energy generator/converter model for inverter-coupled
resources. This document is a skeleton for the model specification; parameters,
equations, initialization details, and default values must be validated against
the REGCB source standard before implementation.

## Block Diagram

Standard model diagram for the REGCB converter interface.

![](../../../../../docs/Figures/PhasorDynamics_REGCB_Diagram.png)

Figure 1: Generator/Converter REGCB model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

Detailed REGCB parameters, variables, equations, initialization details, and
outputs will be added after validation against the REGCB source standard.

<!--
## Model Parameters

TODO: Populate from the validated REGCB source standard.

Symbol | Units | Description | Typical Value | Note
------ | ----- | ----------- | ------------- | ----
TODO | TODO | TODO | TODO | TODO

### Model Derived Parameters

TODO: Add derived parameters and their mathematical definitions after the model
parameter set and limiter behavior are validated.

## Model Variables

### Internal Variables

#### Differential

TODO: Populate from the validated REGCB state definitions.

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
TODO | TODO | TODO | TODO

#### Algebraic

TODO: Populate from the validated REGCB algebraic variables.

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
TODO | TODO | TODO | TODO

### External Variables

#### Differential

None.

#### Algebraic

TODO: Populate from the validated REGCB external ports and signal inputs.

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
TODO | TODO | TODO | TODO

## Model Equations

### Differential Equations

TODO: Add differential equations after the state definitions, limiter behavior,
and smooth approximations are validated.

```math
\begin{aligned}
  \text{TODO}
\end{aligned}
```

### Algebraic Equations

TODO: Add algebraic equations after the network-interface convention and current
limit logic are validated.

```math
\begin{aligned}
  \text{TODO}
\end{aligned}
```

## Initialization

TODO: Add the initialization procedure in the order computations are performed,
including steady-state assumptions and initial values for external commands.

## Model Outputs

TODO: Populate the monitored output list after the algebraic variables are
validated.

Real and imaginary injected currents, when defined by the model, should use
$I_r$ and $I_i$ and be oriented as leaving the converter (i.e. entering the bus).

Active and reactive power outputs, when defined, should follow the GridKit
phasor dynamics convention:

```math
\begin{aligned}
  P &= V_r I_r + V_i I_i \\
  Q &= V_i I_r - V_r I_i
\end{aligned}
```

Power outputs should be oriented as leaving the converter (i.e. entering the bus).
-->
