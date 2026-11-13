# REPCA

Steady state of the PhasorDynamics
[REPCA](../../../PhasorDynamics/Controller/REPCA/README.md), with its parameters
and symbols: the plant limits on the branch power.

## Model Parameters

`mva`, `Freqflag`, `Qmax`, `Qmin`, `Pmax`, and `Pmin` of REPCA.

### Parameter Validation

The conditions of REPCA on these parameters.

### Model Derived Parameters

$k_\mathrm{base}$ of REPCA.

## Model Ports

Name | Port   | Description
-----|--------|----------------------
`p`  | Signal | Branch active power
`q`  | Signal | Branch reactive power

## Model Variables

### Internal Variables

None.

### External Variables

Symbol | Units  | Description           | Note
-------|--------|-----------------------|--------------------
$P$    | [p.u.] | Branch active power   | Owned by the signal
$Q$    | [p.u.] | Branch reactive power | Owned by the signal

## Model Equations

### Objective

None.

### Internal Constraints

#### Equality

None.

#### Inequality

```math
\begin{aligned}
Q^{\min} &\le k_\mathrm{base} Q \le Q^{\max} \\
P^{\min} &\le s_\mathrm{freq} k_\mathrm{base} P \le P^{\max}
\end{aligned}
```

### External Constraints

None.

## Initialization

None.
