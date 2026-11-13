# REECB

Steady state of the PhasorDynamics
[REECB](../../../PhasorDynamics/Controller/REECB/README.md), with its parameters
and symbols: the power and current limits on the converter commands.

## Model Parameters

`mva`, `Pmax`, `Pmin`, `Qmax`, `Qmin`, and `Imax` of REECB.

### Parameter Validation

The conditions of REECB on these parameters.

### Model Derived Parameters

$k_\mathrm{base}$ of REECB.

## Model Ports

Name    | Port   | Description
--------|--------|-------------------------
`ipcmd` | Signal | Active-current command
`iqcmd` | Signal | Reactive-current command
`pe`    | Signal | Active-power feedback
`qgen`  | Signal | Reactive-power feedback

## Model Variables

### Internal Variables

None.

### External Variables

Symbol             | Units  | Description              | Note
-------------------|--------|--------------------------|--------------------
$I_p^\mathrm{cmd}$ | [p.u.] | Active-current command   | Owned by the signal
$I_q^\mathrm{cmd}$ | [p.u.] | Reactive-current command | Owned by the signal
$P_\mathrm{e}$     | [p.u.] | Active-power feedback    | Owned by the signal
$Q^\mathrm{gen}$   | [p.u.] | Reactive-power feedback  | Owned by the signal

## Model Equations

### Objective

None.

### Internal Constraints

#### Equality

None.

#### Inequality

```math
\begin{aligned}
P^{\min} &\le k_\mathrm{base} P_\mathrm{e} \le P^{\max} \\
Q^{\min} &\le k_\mathrm{base} Q^\mathrm{gen} \le Q^{\max} \\
(k_\mathrm{base} I_p^\mathrm{cmd})^2 + (k_\mathrm{base} I_q^\mathrm{cmd})^2 &\le (I^{\max})^2
\end{aligned}
```

### External Constraints

None.

## Initialization

None.
