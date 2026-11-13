# GENROU

Steady state of the PhasorDynamics
[GENROU](../../../PhasorDynamics/SynchronousMachine/GENROU/README.md), with its
parameters and symbols: the mechanical power and field voltage that the terminal current requires.

## Model Parameters

`Ra`, `Xd`, `Xq`, `Xdpp`, `Xl`, `S10`, `S12`, and `mva` of GENROU.

### Parameter Validation

- $X_d > X_\ell$

### Model Derived Parameters

$S_A$, $S_B$, and $X_{qd}$ of GENROU.

## Model Ports

Name    | Port   | Description
--------|--------|-----------------
`bus`   | Bus    | Port bus
`pmech` | Signal | Mechanical power
`efd`   | Signal | Field voltage

## Model Variables

### Internal Variables

Symbol         | Units  | Description                    | Note
---------------|--------|--------------------------------|-------------
$I_\mathrm{r}$ | [p.u.] | Real current into the bus      | Machine base
$I_\mathrm{i}$ | [p.u.] | Imaginary current into the bus | Machine base

### External Variables

Symbol          | Units  | Description       | Note
----------------|--------|-------------------|--------------------
$V_\mathrm{r}$  | [p.u.] | Real voltage      | Owned by the bus
$V_\mathrm{i}$  | [p.u.] | Imaginary voltage | Owned by the bus
$P_\mathrm{m}$  | [p.u.] | Mechanical power  | Owned by the signal
$E_\mathrm{fd}$ | [p.u.] | Field voltage     | Owned by the signal

## Model Equations

Smooth functions: [`qramp`](../../../../CommonMath.md#quadratic-ramp).

```math
\begin{aligned}
\psi_\mathrm{r} &= V_\mathrm{r} + R_\mathrm{a} I_\mathrm{r} - X''_d I_\mathrm{i} \\
\psi_\mathrm{i} &= V_\mathrm{i} + R_\mathrm{a} I_\mathrm{i} + X''_d I_\mathrm{r} \\
k_\mathrm{sat} &= S_B \text{qramp}\left(\sqrt{\psi_\mathrm{r}^2 + \psi_\mathrm{i}^2} - S_A\right) \\
E_\mathrm{r} &= (1 + X_{qd} k_\mathrm{sat}) \psi_\mathrm{r} - (X_q - X''_d) I_\mathrm{i} \\
E_\mathrm{i} &= (1 + X_{qd} k_\mathrm{sat}) \psi_\mathrm{i} + (X_q - X''_d) I_\mathrm{r} \\
E &= \sqrt{E_\mathrm{r}^2 + E_\mathrm{i}^2} \\
I_d &= \dfrac{1}{E} (I_\mathrm{r} E_\mathrm{i} - I_\mathrm{i} E_\mathrm{r}) \\
\psi''_d &= \dfrac{1}{E} (\psi_\mathrm{r} E_\mathrm{r} + \psi_\mathrm{i} E_\mathrm{i})
\end{aligned}
```

### Objective

None.

### Internal Constraints

#### Equality

```math
\begin{aligned}
0 &= -P_\mathrm{m} + \dfrac{S^\mathrm{base}}{S^\mathrm{sys}} (V_\mathrm{r} I_\mathrm{r} + V_\mathrm{i} I_\mathrm{i} + R_\mathrm{a} (I_\mathrm{r}^2 + I_\mathrm{i}^2)) \\
0 &= -E_\mathrm{fd} + (1 + k_\mathrm{sat}) \psi''_d + (X_d - X''_d) I_d
\end{aligned}
```

#### Inequality

```math
\begin{aligned}
-1 &\le I_\mathrm{r} \le 1 \\
-1 &\le I_\mathrm{i} \le 1
\end{aligned}
```

### External Constraints

Power into the bus, added to its balance:

```math
\begin{aligned}
\Delta P^\mathrm{bus} &\mathrel{+}= \dfrac{S^\mathrm{base}}{S^\mathrm{sys}} (V_\mathrm{r} I_\mathrm{r} + V_\mathrm{i} I_\mathrm{i}) \\
\Delta Q^\mathrm{bus} &\mathrel{+}= \dfrac{S^\mathrm{base}}{S^\mathrm{sys}} (V_\mathrm{i} I_\mathrm{r} - V_\mathrm{r} I_\mathrm{i})
\end{aligned}
```

## Initialization

```math
\begin{aligned}
I_\mathrm{r} &\leftarrow 0 \\
I_\mathrm{i} &\leftarrow 0
\end{aligned}
```
