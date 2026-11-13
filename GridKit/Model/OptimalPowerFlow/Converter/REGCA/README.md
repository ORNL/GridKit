# REGCA

Steady state of the PhasorDynamics
[REGCA](../../../PhasorDynamics/Converter/REGCA/README.md), with its parameters
and symbols: the converter currents and the commands and power they require.

## Model Parameters

`mva`, `sL`, `IL1`, `VL1`, `VA1`, and `Vhvmax` of REGCA.

### Parameter Validation

The conditions of REGCA on these parameters.

### Model Derived Parameters

$k_\mathrm{base}$ of REGCA.

## Model Ports

Name      | Port   | Description
----------|--------|-------------------------
`bus`     | Bus    | Port bus
`ipcmd`   | Signal | Active-current command
`iqcmd`   | Signal | Reactive-current command
`pbranch` | Signal | Branch active power
`qbranch` | Signal | Branch reactive power

## Model Variables

### Internal Variables

Symbol | Units  | Description                | Note
-------|--------|----------------------------|---------------
$I_p$  | [p.u.] | Converter active current   | Component base
$I_q$  | [p.u.] | Converter reactive current | Component base

### External Variables

Symbol             | Units  | Description              | Note
-------------------|--------|--------------------------|--------------------
$V_\mathrm{r}$     | [p.u.] | Real voltage             | Owned by the bus
$V_\mathrm{i}$     | [p.u.] | Imaginary voltage        | Owned by the bus
$I_p^\mathrm{cmd}$ | [p.u.] | Active-current command   | Owned by the signal
$I_q^\mathrm{cmd}$ | [p.u.] | Reactive-current command | Owned by the signal
$P^\mathrm{br}$    | [p.u.] | Branch active power      | Owned by the signal
$Q^\mathrm{br}$    | [p.u.] | Branch reactive power    | Owned by the signal

## Model Equations

```math
V_T = \sqrt{V_\mathrm{r}^2 + V_\mathrm{i}^2}
```

### Objective

None.

### Internal Constraints

#### Equality

```math
\begin{aligned}
0 &= -I_p^\mathrm{cmd} + \dfrac{1}{k_\mathrm{base}} I_p \\
0 &= -I_q^\mathrm{cmd} + \dfrac{1}{k_\mathrm{base}} I_q \\
0 &= -P^\mathrm{br} + \dfrac{1}{k_\mathrm{base}} V_T I_p \\
0 &= -Q^\mathrm{br} + \dfrac{1}{k_\mathrm{base}} V_T I_q
\end{aligned}
```

#### Inequality

```math
\begin{aligned}
\max(V_{A1}, s_L V_{L1}) &\le V_T \le V_\mathrm{hv}^{\max} \\
s_L I_p &\le I_{L1}
\end{aligned}
```

### External Constraints

Power into the bus, added to its balance:

```math
\begin{aligned}
\Delta P^\mathrm{bus} &\mathrel{+}= \dfrac{1}{k_\mathrm{base}} V_T I_p \\
\Delta Q^\mathrm{bus} &\mathrel{+}= \dfrac{1}{k_\mathrm{base}} V_T I_q
\end{aligned}
```

## Initialization

```math
\begin{aligned}
I_p &\leftarrow 0 \\
I_q &\leftarrow 0
\end{aligned}
```
