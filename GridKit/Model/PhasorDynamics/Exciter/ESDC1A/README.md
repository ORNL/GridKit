# **IEEE DC1A Excitation System Model (ESDC1A)**

ESDC1A is an IEEE DC1A excitation-system model with a voltage transducer,
input lead-lag compensation, a limited voltage regulator, exciter feedback and
saturation, under-excitation limiter routing, and an optional speed multiplier.

## Notes

- Internal voltage signals are on component base.
- The source diagram labels the optional multiplier input as `Speed`; GridKit
  uses machine speed deviation, so the enabled multiplier is $1+\omega$.
- The UEL selector routes $V_{\mathrm{UEL}}$ either through the high-value gate
  or through the voltage-error summing junction.

## Block Diagram

![ESDC1A exciter block diagram](diagram.png)

Figure 1: ESDC1A exciter model. Figure courtesy of the
[PowerWorld ESDC1A model reference](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Exciter%20ESDC1A.htm).

## Model Parameters

Symbol                              | Units     | JSON      | Description                                     | Typical Value
------------------------------------|-----------|-----------|-------------------------------------------------|--------------
$T_R$                               | [sec]     | `Tr`      | Voltage transducer time constant                | 0.0
$K_A$                               | [p.u.]    | `Ka`      | Voltage-regulator gain                          | 40.0
$T_A$                               | [sec]     | `Ta`      | Voltage-regulator time constant                 | 0.1
$T_B$                               | [sec]     | `Tb`      | Input lead-lag denominator time constant        | 0.0
$T_C$                               | [sec]     | `Tc`      | Input lead-lag numerator time constant          | 0.0
$V_R^{\max}$                        | [p.u.]    | `Vrmax`   | Maximum voltage-regulator output                | 1.0
$V_R^{\min}$                        | [p.u.]    | `Vrmin`   | Minimum voltage-regulator output                | -1.0
$K_E$                               | [p.u.]    | `Ke`      | Exciter constant                                | 0.1
$T_E$                               | [sec]     | `Te`      | Exciter time constant                           | 0.5
$K_F$                               | [p.u.]    | `Kf`      | Stabilizing feedback gain                       | 0.05
$T_{F1}$                            | [sec]     | `Tf1`     | Stabilizing feedback time constant              | 0.7
$s_{\mathrm{spd}}$                  | [boolean] | `Spdmlt`  | Field-voltage speed-multiplier flag             | `false`
$E_1$                               | [p.u.]    | `E1`      | First saturation voltage point                  | 2.8
$S_E(E_1)$                          | [p.u.]    | `Se1`     | Saturation coefficient at $E_1$                 | 0.08
$E_2$                               | [p.u.]    | `E2`      | Second saturation voltage point                 | 3.7
$S_E(E_2)$                          | [p.u.]    | `Se2`     | Saturation coefficient at $E_2$                 | 0.33
$I_{\mathrm{UEL}}$                  | [integer] | `UEL`     | Under-excitation limiter input-routing selector | 0
$s_{\mathrm{lim}}$                  | [boolean] | `exclim`  | Exciter field-voltage-state lower-limit flag     | `true`

Every parameter is optional.
All real-valued parameters must be finite. `Spdmlt` and `exclim` must be
JSON booleans, and `UEL` must be a JSON integer.

### Parameter Validation

Invalid ESDC1A parameter sets are rejected by the following checks:

```math
\begin{aligned}
  K_A
    &> 0 \\
  T_R, T_A, T_B, T_C, T_E, T_{F1}
    &\ge 0 \\
  V_R^{\min}
    &\le V_R^{\max} \\
  s_{\mathrm{spd}}, s_{\mathrm{lim}}
    &\in \{0,1\} \\
  I_{\mathrm{UEL}}
    &\in \{0,1,2,3\}
\end{aligned}
```

The saturation points are either disabled together,

```math
S_E(E_1) = S_E(E_2) = 0,
```

or define a valid two-point quadratic fit:

```math
\begin{aligned}
  E_1, E_2, S_E(E_1), S_E(E_2) &> 0 \\
  \left(E_2-E_1\right)
  \left[S_E(E_2)-S_E(E_1)\right] &> 0
\end{aligned}
```

### Model Derived Parameters

Let $\epsilon_T = 10^{-3}\ \mathrm{s}$. A time constant below $\epsilon_T$ is
raised to that floor in place, so every equation below uses the raised value:

```math
\begin{aligned}
  T_x
    &\leftarrow \max\!\left(T_x,\epsilon_T\right),
       \quad x\in\{R,A,B,E,F1\} \\
  s_{\mathrm{UEL}}
    &=
      \begin{cases}
        1 & I_{\mathrm{UEL}} \ge 2 \\
        0 & I_{\mathrm{UEL}} < 2
      \end{cases}
\end{aligned}
```

When saturation is disabled, $S_A = 0$ and $S_B = 0$. Otherwise,

```math
\begin{aligned}
  C &= \sqrt{\dfrac{S_E(E_2)}{S_E(E_1)}} \\
  S_A &= \dfrac{C E_1 - E_2}{C - 1} \\
  S_B &= \dfrac{S_E(E_1)}{(E_1 - S_A)^2}
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

`Known` ports hold their initial values before `initialize()` and are preserved
by it. `Unknown` inputs are resolved during initialization and written to
attached signal storage, or retained as constant inputs when unattached. The
`efd` output must be assigned. The `speed` input is required when
$s_{\mathrm{spd}} = 1$; every other signal input is optional. Unattached `speed`,
`vs`, and `vuel` inputs default to zero.

## Model Variables

### Internal Variables

#### Differential

Symbol                              | Units  | Description                                      | Note
------------------------------------|--------|--------------------------------------------------|------
$E_{\mathrm{fd}}'$                  | [p.u.] | Exciter field-voltage state                     | State 1 in Fig. 1; lower bounded at zero when $s_{\mathrm{lim}} = 1$; before the optional speed multiplier
$V_C$                               | [p.u.] | Filtered terminal-voltage magnitude              | State 2 in Fig. 1
$V_R$                               | [p.u.] | Voltage-regulator output                         | State 3 in Fig. 1
$V_F$                               | [p.u.] | Stabilizing feedback state                       | State 4 in Fig. 1
$x_{\mathrm{LL}}$                   | [p.u.] | Input lead-lag denominator state                 | State 5 in Fig. 1

#### Algebraic

Symbol                              | Units  | Description                       | Note
------------------------------------|--------|-----------------------------------|------
$e_V$                               | [p.u.] | Voltage-error summing output      |
$V_{\mathrm{LL}}$                   | [p.u.] | Input lead-lag output             |
$V_{\mathrm{HV}}$                   | [p.u.] | High-value gate output            |
$S_E$                               | [p.u.] | Exciter saturation coefficient    | Evaluated at $E_{\mathrm{fd}}'$
$V_{\mathrm{FE}}$                   | [p.u.] | Exciter feedback drive            |
$E_{\mathrm{fd}}$                   | [p.u.] | Field-voltage output              | Published through `efd`

### External Variables

#### Differential

None.

#### Algebraic

Symbol                              | Units  | Init    | Description                            | Note
------------------------------------|--------|---------|----------------------------------------|------
$V_{\mathrm{r}}$                    | [p.u.] | Known   | Terminal voltage, real component       | Bus input
$V_{\mathrm{i}}$                    | [p.u.] | Known   | Terminal voltage, imaginary component  | Bus input
$\omega$                            | [p.u.] | Known   | Machine speed deviation                | Signal port `speed`
$V_{\mathrm{ref}}$                  | [p.u.] | Unknown | Voltage-control reference              | Signal port `vref`
$V_S$                               | [p.u.] | Known   | Stabilizer input signal                | Signal port `vs`
$V_{\mathrm{UEL}}$                  | [p.u.] | Known   | Under-excitation limiter input         | Signal port `vuel`

## Model Equations

Define the pre-limit exciter field-voltage rate:

```math
f_E = \dfrac{V_R-V_{\mathrm{FE}}}{T_E}.
```

### Differential Equations

```math
\begin{aligned}
  0 &=
    -\dot{E}_{\mathrm{fd}}'
    + \left(1-s_{\mathrm{lim}}\right)f_E
    + s_{\mathrm{lim}}\,
      \text{awmin}\left(E_{\mathrm{fd}}',f_E;0\right) \\
  0 &=
    -\dot{V}_C
    + \dfrac{1}{T_R}
      \left(
        \sqrt{V_{\mathrm{r}}^2+V_{\mathrm{i}}^2}
        - V_C
      \right) \\
  0 &=
    -\dot{V}_R
    + \dfrac{1}{T_A}
      \text{antiwindup}
      \left(
        V_R,\,
        -V_R + K_A V_{\mathrm{HV}};\,
        V_R^{\min}, V_R^{\max}
      \right) \\
  0 &=
    -\dot{V}_F
    + \dfrac{1}{T_{F1}}
      \left[
        -V_F
        + \dfrac{K_F}{T_E}
          \left(V_R - V_{\mathrm{FE}}\right)
      \right] \\
  0 &=
    -\dot{x}_{\mathrm{LL}}
    + \dfrac{1}{T_B}
      \left(e_V - x_{\mathrm{LL}}\right)
\end{aligned}
```

The field-voltage-state limiter uses the fixed-lower-bound anti-windup rule
of [Appendix A](#appendix-a-awmin).

### Algebraic Equations

```math
\begin{aligned}
  0 &=
    -e_V
    + V_{\mathrm{ref}}
    + V_S
    + s_{\mathrm{UEL}}V_{\mathrm{UEL}}
    - V_C
    - V_F \\
  0 &=
    -V_{\mathrm{LL}}
    + x_{\mathrm{LL}}
    + \dfrac{T_C}{T_B}
      \left(e_V - x_{\mathrm{LL}}\right) \\
  0 &=
    -V_{\mathrm{HV}}
    + \begin{cases}
        \text{max}\left(V_{\mathrm{LL}}, V_{\mathrm{UEL}}\right)
          & s_{\mathrm{UEL}} = 0 \\
        V_{\mathrm{LL}}
          & s_{\mathrm{UEL}} = 1
      \end{cases} \\
  0 &=
    -S_E
    + S_B q\left(E_{\mathrm{fd}}' - S_A\right) \\
  0 &=
    -V_{\mathrm{FE}}
    + \left(K_E + S_E\right)E_{\mathrm{fd}}' \\
  0 &=
    -E_{\mathrm{fd}}
    + \left(1 + s_{\mathrm{spd}}\omega\right)E_{\mathrm{fd}}'
\end{aligned}
```

CommonMath defines helper targets and smooth approximations for
[max](../../../../CommonMath.md#derived-functions), the [ramp](../../../../CommonMath.md#primitives)
$\rho$, and the [quadratic ramp](../../../../CommonMath.md#primitives) $q$.

## Initialization

### Input Initialization

```math
\begin{aligned}
  V_{\mathrm{r}}, V_{\mathrm{i}}
    &\leftarrow \text{terminal-bus voltage} \\
  E_{\mathrm{fd}}
    &\leftarrow \text{machine field voltage} \\
  \omega
    &\leftarrow \text{machine speed deviation or }0 \\
  V_S
    &\leftarrow \text{stabilizer signal or }0 \\
  V_{\mathrm{UEL}}
    &\leftarrow \text{under-excitation limiter input or }0
\end{aligned}
```

Initialization never replaces the seeded value held in $E_{\mathrm{fd}}$.

### Internal Initialization

All internal derivatives are set to zero. The steady-state residuals are then
resolved in dependency order. The smooth high-value gate requires its input to
be recovered through the inverse CommonMath
[ramp](../../../../CommonMath.md#primitives) $\rho^{-1}$ when the UEL input is
routed through the gate:

```math
\begin{aligned}
  V_C
    &\leftarrow \sqrt{V_{\mathrm{r}}^2+V_{\mathrm{i}}^2} \\
  E_{\mathrm{fd}}'
    &\leftarrow
      \dfrac{E_{\mathrm{fd}}}{1 + s_{\mathrm{spd}}\omega} \\
  S_E
    &\leftarrow S_B q\left(E_{\mathrm{fd}}' - S_A\right) \\
  V_{\mathrm{FE}}
    &\leftarrow \left(K_E + S_E\right)E_{\mathrm{fd}}' \\
  V_R
    &\leftarrow V_{\mathrm{FE}} \\
  V_{\mathrm{HV}}
    &\leftarrow \dfrac{V_R}{K_A} \\
  V_{\mathrm{LL}}
    &\leftarrow
      \begin{cases}
        V_{\mathrm{UEL}}
          + \rho^{-1}
            \left(V_{\mathrm{HV}}-V_{\mathrm{UEL}}\right)
          & s_{\mathrm{UEL}} = 0 \\
        V_{\mathrm{HV}}
          & s_{\mathrm{UEL}} = 1
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
$s_{\mathrm{lim}}=1$, $V_R$ outside
$[V_R^{\min},V_R^{\max}]$, and high-value-gate active starts with
$s_{\mathrm{UEL}} = 0$ and
$V_{\mathrm{HV}}\le V_{\mathrm{UEL}}$.

Every check resolves before any storage is written, so a rejected
initialization leaves state, the `efd` seed, and external signals unchanged.

### Output Initialization

```math
\begin{aligned}
  V_{\mathrm{ref}}
    &\leftarrow
      e_V
      + V_C
      + V_F
      - V_S
      - s_{\mathrm{UEL}}V_{\mathrm{UEL}}
\end{aligned}
```

ESDC1A writes the resolved voltage-control reference to an attached `vref`
signal input. If no controller is connected, that value is used as a constant
reference input.

## Monitorable Outputs

Output          | Units  | Description                         | Note
----------------|--------|-------------------------------------|------
`efd`           | [p.u.] | Field-voltage output                | $E_{\mathrm{fd}}$
`vc`            | [p.u.] | Filtered terminal-voltage magnitude | $V_C$
`vr`            | [p.u.] | Voltage-regulator output            | $V_R$
`vf`            | [p.u.] | Stabilizing feedback state          | $V_F$
`se`            | [p.u.] | Exciter saturation coefficient      | $S_E$
`vfe`           | [p.u.] | Exciter feedback drive              | $V_{\mathrm{FE}}$

## Testing

- `validation()` checks construction, documented defaults, parameter types
  and domains, signal configuration, and minimum time-constant handling.
- `initializationAndSignals()` checks steady initialization, selector
  combinations, signal publication and latching, monitor output, and
  differentiability tags.
- `initializationDomain()` checks rejected and accepted field-voltage,
  terminal-voltage, Known-input, speed-multiplier, regulator-limit, and
  high-value-gate operating points.
- `residualEquations()` checks every model residual against a fixed
  numerical answer key.
- `voltageRegulation()` checks the transducer, summing junction, lead-lag,
  stabilizing feedback, and regulator anti-windup behavior.
- `excitationLimits()` checks high-value-gate routing, saturation,
  field-voltage-state limiting, and the optional speed multiplier.
- `jacobian()` compares the dependency-tracking and Enzyme Jacobians when
  Enzyme support is enabled.

## Appendix A: `awmin`

The exact anti-windup rule at a fixed lower bound $\ell$ is

```math
\text{awmin}(x,f;\ell) =
  \begin{cases}
    f & x > \ell \\
    \text{max}(f,0) & x \le \ell
  \end{cases}
```

Above the bound the unconstrained derivative passes. At or below the bound,
outward motion is blocked and restoring motion is admitted.

The model evaluates this rule with the following smooth approximation:

```math
\text{awmin}(x,f;\ell)
  \approx
  \left[
    \sigma(f)
    + \left(1-\sigma(f)\right)\text{above}(x;\ell)
  \right]f.
```

CommonMath defines the [`above`](../../../../CommonMath.md#derived-functions)
and [`sigmoid`](../../../../CommonMath.md#primitives) targets and smooth
approximations.
