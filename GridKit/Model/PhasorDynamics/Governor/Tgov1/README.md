# TGOV1

Steam turbine-governor model.

## Block Diagram

![](../../../../../docs/Figures/TGOV1.JPG)

Figure 1: Governor TGOV1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol              | Units  | JSON    | Description                     | Typical Value | Note
--------------------|--------|---------|---------------------------------|---------------|-----
$T_{\mathrm{rate}}$ | [MVA]  | `Trate` | Governor component power base   | 100.0         |
$R$                 | [p.u.] | `R`     | Permanent droop                 | 0.05          |
$T_1$               | [s]    | `T1`    | Steam-bowl time constant        | 0.5           |
$T_2$               | [s]    | `T2`    | Turbine numerator time constant | 2.5           |
$T_3$               | [s]    | `T3`    | Reheater time constant          | 7.5           |
$P_\mathrm{v}^{\max}$        | [p.u.] | `Pvmax` | Maximum valve position          | 1             |
$P_\mathrm{v}^{\min}$        | [p.u.] | `Pvmin` | Minimum valve position          | 0             |
$D_\mathrm{t}$               | [p.u.] | `Dt`    | Turbine damping coefficient     | 0             |

### Parameter Validation

A valid TGOV1 parameter set must satisfy the following conditions:

```math
T_\mathrm{rate},S_\mathrm{sys}>0,\qquad R\ne0,\qquad P_\mathrm{v}^{\min}\le P_\mathrm{v}^{\max}
```

Set $T_{\mathrm{rate}}$ equal to the connected machine MVA base. A zero
component power base is not supported.

### Model Derived Parameters

Let $\epsilon_T=10^{-3}\ \mathrm{s}$. A time constant below $\epsilon_T$ is
raised to that floor in place, so every equation below uses the raised value:

```math
\begin{aligned}
  T_x
    &\leftarrow \max\!(T_x,\epsilon_T),
       \quad x\in\{1,3\}
\end{aligned}
```

## Model Ports

Name    | Port   | Init    | Description
--------|--------|---------|------------
`speed` | Input  | Known   | Machine speed deviation; optional, defaults to zero
`pref`  | Input  | Unknown | Governor reference; optional
`pmech` | Output | Known   | Mechanical-power signal seeded by the machine

## Model Variables

### Internal Variables

#### Differential

Symbol | Units  | Description          | Note
-------|--------|----------------------|---------------
$P_\mathrm{v}$  | [p.u.] | Valve position       | Component base
$P_\mathrm{t}$  | [p.u.] | Turbine-block output | Component base

#### Algebraic

Symbol          | Units  | Description                       | Note
----------------|--------|-----------------------------------|-------
$P_\mathrm{m}$           | [p.u.] | Mechanical-power output           | System base; read by the machine model

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

Smooth functions: [`antiwindup`](../../../../CommonMath.md#antiwindup).

For readability, define:

```math
g_v=-P_\mathrm{v}+\dfrac{P_\mathrm{ref}-\omega}{R}
```

### Internal Equations

#### Differential

```math
\begin{aligned}
  0 &= -\dot P_\mathrm{v}
       + \dfrac{1}{T_1}\text{antiwindup}
         (P_\mathrm{v},g_v;P_\mathrm{v}^{\min},P_\mathrm{v}^{\max}) \\
  0 &= -\dot P_\mathrm{t}-\dfrac{P_\mathrm{t}-P_\mathrm{v}-T_2\dot P_\mathrm{v}}{T_3}
\end{aligned}
```

#### Algebraic

```math
0=-\dfrac{S_\mathrm{sys}}{T_\mathrm{rate}}P_\mathrm{m}
  +P_\mathrm{t}-D_\mathrm{t}\omega
```

### External Equations

None.

## Initialization

The machine provides $P_\mathrm{m}$ on system base:

```math
\begin{aligned}
P_\mathrm{v} &\leftarrow \dfrac{S_\mathrm{sys}}{T_\mathrm{rate}}P_\mathrm{m}+D_\mathrm{t}\omega \\
P_\mathrm{t} &\leftarrow P_\mathrm{v} \\
P_\mathrm{ref} &\leftarrow \omega+RP_\mathrm{v} \\
\dot P_\mathrm{v},\dot P_\mathrm{t} &\leftarrow 0
\end{aligned}
```

Initialization rejects $P_\mathrm{v}$ outside the configured valve limits.

## Monitors

None.
