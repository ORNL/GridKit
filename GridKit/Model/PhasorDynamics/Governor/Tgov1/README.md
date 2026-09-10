# **Steam Turbine-Governor Model (TGOV1)**

## Block Diagram

Standard model of the stream turbine

![](../../../../../docs/Figures/TGOV1.JPG)

Figure 1: Governor TGOV1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol      | Units  | Description                       | Typical Value | Note
------------|--------|-----------------------------------|---------------| ------
$T_{\mathrm{rate}}$ | [MVA] | Governor component power base     | 100.0 |
$R$         | [p.u.] | Permanent droop                   | 0.05 |
$T_1$       | [sec]  | Steam-bowl time constant          | 0.5  |
$T_2$       | [sec]  | Turbine numerator time constant   | 2.5  |
$T_3$       | [sec]  | Reheater time constant            | 7.5  |
$P_v^{\max}$ | [p.u.] | Maximum valve position           | 1    |
$P_v^{\min}$ | [p.u.] | Minimum valve position           | 0    |
$D_t$       | [p.u.] | Turbine damping coefficient       | 0    |

### Parameter Validation

The component and system power bases must be positive, $R$ must be nonzero,
and $P_v^{\min}\le P_v^{\max}$.

Set $T_{\mathrm{rate}}$ equal to the connected machine MVA base. A zero
component power base is not supported.

### Model Derived Parameters

Time constants are retained as supplied. A zero lag time constant makes its
state algebraic; positive values remain differential.
At $T_1=0$, the valve position equals $(P_{\mathrm{ref}}-\omega)/R$
clamped to its valve limits.

## Model Ports

Name    | Port   | Init    | Description
--------|--------|---------|------------
`speed` | Input  | Known   | Machine speed deviation; optional, defaults to zero
`pref`  | Input  | Unknown | Governor reference; optional
`pmech` | Output | Known   | Mechanical-power signal seeded by the machine

## Model Variables

### Internal Variables

#### Differential

Symbol    | Units  | Description                       | Note
----------|--------|-----------------------------------|-------
$P_t$     | [p.u.] | Turbine-block output              | Component base
$P_v$     | [p.u.] | Valve position                    | Component base

#### Algebraic

Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$P_m$           | [p.u.] | Mechanical-power output           | System base; read by the machine model

### External Variables

#### Differential

Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$\omega$  | [p.u.] | Machine speed deviation           | Optional `speed` input; defaults to zero

#### Algebraic

Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$P_\mathrm{ref}$ | [p.u.] | Governor reference               | Component base; optional `pref` input, otherwise held internally

## Model Equations

For readability, define:

```math
g_v=-P_v+\dfrac{P_\mathrm{ref}-\omega}{R}.
```

### Differential Equations

The TGOV1 differential equations, as derived from the model diagram, are

```math
\begin{aligned}
  0 &= -T_1\dot P_v
       + \text{antiwindup}
         \left(P_v,g_v;P_v^{\min},P_v^{\max}\right) \\
  0 &= -T_3\dot P_t-P_t+P_v+T_2\dot P_v.
\end{aligned}
```

CommonMath defines the [Antiwindup](../../../../CommonMath.md#antiwindup)
target and smooth approximation.

### Algebraic Equations

The mechanical-power output is given by

```math
0=-\dfrac{S_\mathrm{sys}}{T_\mathrm{rate}}P_m
  +P_t-D_t\omega.
```

## Initialization

TGOV1 preserves the machine-provided $P_{m,0}$ and initializes the steady
state in dependency order:

```math
\begin{aligned}
  P_{m,0}^{\mathrm{TGOV1}}
    &\leftarrow \dfrac{S_\mathrm{sys}}{T_\mathrm{rate}}P_{m,0} \\
  P_{v,0}
    &\leftarrow P_{m,0}^{\mathrm{TGOV1}}+D_t\omega_0 \\
  P_{t,0}
    &\leftarrow P_{v,0} \\
  P_{\mathrm{ref},0}
    &\leftarrow \omega_0+RP_{v,0} \\
  \dot P_{v,0},\dot P_{t,0}
    &\leftarrow 0.
\end{aligned}
```

Initialization rejects $P_{v,0}$ outside the configured valve limits.
