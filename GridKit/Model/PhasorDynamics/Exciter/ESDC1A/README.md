# ESDC1A

ESDC1A is an IEEE DC1A excitation-system model with a voltage transducer,
input lead–lag compensation, a limited voltage regulator, exciter feedback and
saturation, under-excitation limiter routing, and an optional speed multiplier.

## Notes

- Internal voltage signals are on component base.
- The source diagram labels the optional multiplier input as `Speed`; GridKit
  uses machine speed deviation, so the enabled multiplier is $1+\omega$.
- The UEL selector routes $V_{\mathrm{uel}}$ either through the high-value gate
  or through the voltage-error summing junction.

## Block Diagram

![ESDC1A exciter block diagram](../../../../../docs/Figures/PhasorDynamics/ESDC1A/diagram.png)

Figure 1: ESDC1A exciter model. Figure courtesy of the
[PowerWorld ESDC1A model reference](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Exciter%20ESDC1A.htm).

## Model Parameters

Symbol             | Units     | JSON     | Description                                                                                          | Typical Value | Note
-------------------|-----------|----------|------------------------------------------------------------------------------------------------------|---------------|-----
$T_R$              | [s]       | `Tr`     | Voltage transducer time constant                                                                     | 0.0           |
$K_A$              | [p.u.]    | `Ka`     | Voltage-regulator gain                                                                               | 40.0          |
$T_A$              | [s]       | `Ta`     | Voltage-regulator time constant                                                                      | 0.1           |
$T_B$              | [s]       | `Tb`     | Input lead–lag denominator time constant                                                             | 0.0           |
$T_C$              | [s]       | `Tc`     | Input lead–lag numerator time constant                                                               | 0.0           |
$V_R^{\max}$       | [p.u.]    | `Vrmax`  | Maximum voltage-regulator output                                                                     | 1.0           |
$V_R^{\min}$       | [p.u.]    | `Vrmin`  | Minimum voltage-regulator output                                                                     | -1.0          |
$K_E$              | [p.u.]    | `Ke`     | Exciter field resistance line slope margin; 0 requests automatic calculation, not a zero coefficient | 0.1           |
$T_E$              | [s]       | `Te`     | Exciter time constant                                                                                | 0.5           |
$K_F$              | [p.u.]    | `Kf`     | Stabilizing feedback gain                                                                            | 0.05          |
$T_{F1}$           | [s]       | `Tf1`    | Stabilizing feedback time constant                                                                   | 0.7           |
$s_{\mathrm{spd}}$ | [boolean] | `Spdmlt` | Field-voltage speed-multiplier flag                                                                  | `false`       |
$E_1$              | [p.u.]    | `E1`     | First saturation voltage point                                                                       | 2.8           |
$S_E(E_1)$         | [p.u.]    | `Se1`    | Saturation coefficient at $E_1$                                                                      | 0.08          |
$E_2$              | [p.u.]    | `E2`     | Second saturation voltage point                                                                      | 3.7           |
$S_E(E_2)$         | [p.u.]    | `Se2`    | Saturation coefficient at $E_2$                                                                      | 0.33          |
$I_{\mathrm{uel}}$ | [integer] | `UEL`    | Under-excitation limiter input-routing selector                                                      | 0             |
$s_{\mathrm{lim}}$ | [boolean] | `exclim` | Exciter field-voltage-state lower-limit flag                                                         | `true`        |

Every parameter is optional.
All real-valued parameters must be finite. `Spdmlt` and `exclim` must be
JSON booleans, and `UEL` must be a JSON integer.

### Parameter Validation

A valid ESDC1A parameter set must satisfy the following conditions:

```math
\begin{aligned}
  K_A
    &> 0 \\
  T_R, T_A, T_B, T_E, T_{F1}
    &\ge 0 \\
  V_R^{\min}
    &\le V_R^{\max} \\
  I_{\mathrm{uel}}
    &\in \{0,1,2,3\}
\end{aligned}
```

The saturation points are either disabled together,

```math
S_E(E_1) = S_E(E_2) = 0,
```

or define a valid two-point scaled-quadratic fit:

```math
\begin{aligned}
  E_1, E_2 &> 0 \\
  S_E(E_1), S_E(E_2) &\ge 0 \\
  (E_2-E_1)
  \left[S_E(E_2)-S_E(E_1)\right] &> 0
\end{aligned}
```

### Model Derived Parameters

Let $\epsilon_T = 10^{-3}\ \mathrm{s}$. A time constant below $\epsilon_T$ is
raised to that floor in place, so every equation below uses the raised value:

```math
\begin{aligned}
  T_x
    &\leftarrow \max\!(T_x,\epsilon_T),
       \quad x\in\{R,A,B,E,F1\} \\
  s_{\mathrm{uel}}
    &=
      \begin{cases}
        1 & I_{\mathrm{uel}} \ge 2 \\
        0 & I_{\mathrm{uel}} < 2
      \end{cases}
\end{aligned}
```
When saturation is disabled, $S_A = 0$ and $S_B = 0$. Otherwise,

```math
E S_E(E) = S_B q(E-S_A)
```

When one saturation value is zero,

```math
\begin{aligned}
  S_E(E_1)=0 &: \quad S_A=E_1,\qquad
    S_B=\dfrac{E_2S_E(E_2)}{(E_2-E_1)^2} \\
  S_E(E_2)=0 &: \quad S_A=E_2,\qquad
    S_B=\dfrac{E_1S_E(E_1)}{(E_1-E_2)^2}
\end{aligned}
```

and when both saturation values are positive,

```math
\begin{aligned}
  C &= \sqrt{\dfrac{E_2S_E(E_2)}{E_1S_E(E_1)}} \\
  S_A &= \dfrac{C E_1 - E_2}{C - 1} \\
  S_B &= \dfrac{E_1S_E(E_1)}{(E_1 - S_A)^2}
\end{aligned}
```

$K_E$ is the configured exciter field resistance line slope margin. A nonzero configured value is used directly. The
[IEEE Std 421.5-2016](https://standards.ieee.org/ieee/421.5/5356/) states:

> In some programs, if $K_E$ is entered as zero, $K_E$ is automatically calculated by the program to represent a self-excited shunt field and a trimmed rheostat as its initial condition.

GridKit implements the
[PSS/E-compatible automatic-parameter rule](https://www.powerworld.com/WebHelp/Content/MainDocumentation_HTML/Transient_Stability_Dialog_Options_Power_System_Model.htm), as specified in the
[PSS/E 33.4 Program Application Guide, Vol. II](https://pdfcoffee.com/pagv2-pdf-free.html).
The divisor $10$ in $V_R=V_R^{\max}/10=0.1V_R^{\max}$ is unitless and
sets $V_R$ to 10% of the maximum regulator output, retaining the per-unit
units of $V_R^{\max}$:

```math
\begin{aligned}
0
  &= V_R
     -K_E^{\mathrm{eff}}E_{\mathrm{fd}}'
     -s_e\\
K_E^{\mathrm{eff}}
  &=
    \begin{cases}
      \dfrac{1}{E_{\mathrm{fd}}'}
      \left(\dfrac{V_R^{\max}}{10}-s_e\right)
        & K_E=0\\
      K_E
        & K_E\ne 0
    \end{cases}
\end{aligned}
```

## Model Ports

Name    | Port   | Init    | Description
--------|--------|---------|------
`bus`   | Bus    | Known   | Terminal bus voltage
`speed` | Input  | Known   | Machine speed deviation
`vref`  | Input  | Unknown | Voltage-control reference
`vs`    | Input  | Known   | Stabilizer input signal
`vuel`  | Input  | Known   | Under-excitation limiter input
`efd`   | Output | Known   | Field-voltage output

`Known` values are seeded before initialization and preserved. `Unknown` inputs
are initialized in attached signal storage or held constant when unattached. The
`efd` output must be assigned. The `speed` input is required when
$s_{\mathrm{spd}} = 1$; every other signal input is optional. Unattached `speed`,
`vs`, and `vuel` inputs default to zero.

## Model Variables

### Internal Variables

#### Differential

Symbol             | Units  | Description                         | Note
-------------------|--------|-------------------------------------|-----------------------------------------------------------------------------------------------------------
$E_{\mathrm{fd}}'$ | [p.u.] | Exciter field-voltage state         | State 1 in Fig. 1; lower bounded at zero when $s_{\mathrm{lim}} = 1$; before the optional speed multiplier
$V_C$              | [p.u.] | Filtered terminal-voltage magnitude | State 2 in Fig. 1
$V_R$              | [p.u.] | Voltage-regulator output            | State 3 in Fig. 1
$V_F$              | [p.u.] | Stabilizing feedback state          | State 4 in Fig. 1
$x_{\mathrm{LL}}$  | [p.u.] | Input lead–lag denominator state    | State 5 in Fig. 1

#### Algebraic

Symbol            | Units  | Description                              | Note
------------------|--------|------------------------------------------|----------------------------------------
$e_V$             | [p.u.] | Voltage-error summing output             |
$V_{\mathrm{LL}}$ | [p.u.] | Input lead–lag output                    |
$V_{\mathrm{HV}}$ | [p.u.] | High-value gate output                   |
$s_e$             | [p.u.] | Scaled-quadratic saturation contribution | $E_{\mathrm{fd}}'S_E(E_{\mathrm{fd}}')$
$V_{\mathrm{FE}}$ | [p.u.] | Exciter feedback drive                   |
$E_{\mathrm{fd}}$ | [p.u.] | Field-voltage output                     | Published through `efd`

### External Variables

#### Differential

None.

#### Algebraic

Symbol             | Units  | Description                           | Note
-------------------|--------|---------------------------------------|--------------------
$V_r$              | [p.u.] | Terminal voltage, real component      | Bus input
$V_i$              | [p.u.] | Terminal voltage, imaginary component | Bus input
$\omega$           | [p.u.] | Machine speed deviation               | Signal port `speed`
$V_{\mathrm{ref}}$ | [p.u.] | Voltage-control reference             | Signal port `vref`
$V_S$              | [p.u.] | Stabilizer input signal               | Signal port `vs`
$V_{\mathrm{uel}}$ | [p.u.] | Under-excitation limiter input        | Signal port `vuel`

## Model Equations

Smooth functions: [`antiwindup`](../../../../CommonMath.md#antiwindup), [`max`](../../../../CommonMath.md#maximum), [$q$](../../../../CommonMath.md#quadratic-ramp), [$\rho$](../../../../CommonMath.md#ramp).

### Internal Equations

#### Differential

Define the pre-limit exciter field-voltage rate:

```math
f_E = \dfrac{V_R-V_{\mathrm{FE}}}{T_E}
```

```math
\begin{aligned}
  0 &=
    -\dot{E}_{\mathrm{fd}}'
    + (1-s_{\mathrm{lim}})f_E
    + s_{\mathrm{lim}}\,
      \text{awmin}(E_{\mathrm{fd}}',f_E;0) \\
  0 &=
    -\dot{V}_C
    + \dfrac{1}{T_R}
      \left(
        \sqrt{V_r^2+V_i^2}
        - V_C
      \right) \\
  0 &=
    -\dot{V}_R
    + \dfrac{1}{T_A}
      \text{antiwindup}
      (
        V_R,\;
        -V_R + K_A V_{\mathrm{HV}};\,
        V_R^{\min}, V_R^{\max}
      ) \\
  0 &=
    -\dot{V}_F
    + \dfrac{1}{T_{F1}}
      \left[
        -V_F
        + \dfrac{K_F}{T_E}
          (V_R - V_{\mathrm{FE}})
      \right] \\
  0 &=
    -\dot{x}_{\mathrm{LL}}
    + \dfrac{1}{T_B}
      (e_V - x_{\mathrm{LL}})
\end{aligned}
```

The field-voltage-state limiter uses the fixed-lower-bound anti-windup rule
of [Appendix A](#appendix-a-awmin).

#### Algebraic

```math
\begin{aligned}
  0 &=
    -e_V
    + V_{\mathrm{ref}}
    + V_S
    + s_{\mathrm{uel}}V_{\mathrm{uel}}
    - V_C
    - V_F \\
  0 &=
    -V_{\mathrm{LL}}
    + x_{\mathrm{LL}}
    + \dfrac{T_C}{T_B}
      (e_V - x_{\mathrm{LL}}) \\
  0 &=
    -V_{\mathrm{HV}}
    + \begin{cases}
        \text{max}(V_{\mathrm{LL}}, V_{\mathrm{uel}})
          & s_{\mathrm{uel}} = 0 \\
        V_{\mathrm{LL}}
          & s_{\mathrm{uel}} = 1
      \end{cases} \\
  0 &=
    -s_e
    + S_B q(E_{\mathrm{fd}}' - S_A) \\
  0 &=
    -V_{\mathrm{FE}}
    + K_E^{\mathrm{eff}} E_{\mathrm{fd}}'
    + s_e \\
  0 &=
    -E_{\mathrm{fd}}
    + (1 + s_{\mathrm{spd}}\omega)E_{\mathrm{fd}}'
\end{aligned}
```

### External Equations

None.

## Initialization

### Input Initialization

```math
\begin{aligned}
  V_r, V_i
    &\leftarrow \text{terminal-bus voltage} \\
  E_{\mathrm{fd}}
    &\leftarrow \text{machine field voltage} \\
  \omega
    &\leftarrow \text{machine speed deviation or }0 \\
  V_S
    &\leftarrow \text{stabilizer signal or }0 \\
  V_{\mathrm{uel}}
    &\leftarrow \text{under-excitation limiter input or }0
\end{aligned}
```

Initialization never replaces the seeded value held in $E_{\mathrm{fd}}$.

### Internal Initialization

All internal derivatives initialize to zero. The UEL high-value gate uses the
inverse smooth [ramp](../../../../CommonMath.md#ramp) $\rho^{-1}$:

```math
\begin{aligned}
  V_C
    &\leftarrow \sqrt{V_r^2+V_i^2} \\
  E_{\mathrm{fd}}'
    &\leftarrow
      \dfrac{E_{\mathrm{fd}}}{1 + s_{\mathrm{spd}}\omega} \\
  s_e
    &\leftarrow S_B q(E_{\mathrm{fd}}' - S_A) \\
  V_{\mathrm{FE}}
    &\leftarrow K_E^{\mathrm{eff}} E_{\mathrm{fd}}' + s_e \\
  V_R
    &\leftarrow V_{\mathrm{FE}} \\
  V_{\mathrm{HV}}
    &\leftarrow \dfrac{V_R}{K_A} \\
  V_{\mathrm{LL}}
    &\leftarrow
      \begin{cases}
        V_{\mathrm{uel}}
          + \rho^{-1}
            (V_{\mathrm{HV}}-V_{\mathrm{uel}})
          & s_{\mathrm{uel}} = 0 \\
        V_{\mathrm{HV}}
          & s_{\mathrm{uel}} = 1
      \end{cases} \\
  V_F
    &\leftarrow 0 \\
  e_V
    &\leftarrow V_{\mathrm{LL}} \\
  x_{\mathrm{LL}}
    &\leftarrow e_V
\end{aligned}
```

Initialization rejects a non-finite or zero bus-voltage magnitude, a
non-finite field-voltage seed, non-finite Known signal inputs, a nonpositive
speed multiplier $1 + s_{\mathrm{spd}}\omega$, $E_{\mathrm{fd}}'<0$ while
$s_{\mathrm{lim}}=1$, initial $V_R$ outside $[V_R^{\min},V_R^{\max}]$,
and high-value-gate active
starts with $s_{\mathrm{uel}} = 0$ and
$V_{\mathrm{HV}}\le V_{\mathrm{uel}}$.

A rejected initialization leaves states and signals unchanged.

### Output Initialization

```math
\begin{aligned}
  V_{\mathrm{ref}}
    &\leftarrow
      e_V
      + V_C
      + V_F
      - V_S
      - s_{\mathrm{uel}}V_{\mathrm{uel}}
\end{aligned}
```

ESDC1A writes the resolved voltage-control reference to an attached `vref`
signal input. If no controller is connected, that value is used as a constant
reference input.

## Monitors

Monitor         | Units  | Description                         | Note
----------------|--------|-------------------------------------|------
`efd`           | [p.u.] | Field-voltage output                | $E_{\mathrm{fd}}$
`vc`            | [p.u.] | Filtered terminal-voltage magnitude | $V_C$
`vr`            | [p.u.] | Voltage-regulator output            | $V_R$
`vf`            | [p.u.] | Stabilizing feedback state          | $V_F$
`se`            | [p.u.] | Scaled-quadratic saturation contribution | $s_e$
`vfe`           | [p.u.] | Exciter feedback drive              | $V_{\mathrm{FE}}$

## Appendix A: `awmin`

The exact anti-windup rule at a fixed lower bound $\ell$ is

```math
\text{awmin}(x,f;\ell) =
  \begin{cases}
    f & x > \ell \\
    \text{max}(f,0) & x \le \ell
  \end{cases}
```

The model evaluates this rule with the following smooth approximation:

```math
\text{awmin}(x,f;\ell)
  \approx
  \left[
    \sigma(f)
    + (1-\sigma(f))\text{above}(x;\ell)
  \right]f
```

CommonMath defines the [`above`](../../../../CommonMath.md#above)
and [`sigmoid`](../../../../CommonMath.md#logistic-function) targets and smooth
approximations.
