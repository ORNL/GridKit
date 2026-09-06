# **IEEE Stabilizer Model (IEEEST)**

Standard IEEE power system stabilizer: 4th-order notch filter, two lead–lag
blocks, washout, and output limiter.

## Block Diagram

![](../../../../../../docs/Figures/stabilizer_ieeest_diagram.png)

Figure 1: Stabilizer IEEEST model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol      | Units  | Description                          | Typical Value
------------|--------|--------------------------------------|--------------
$A_1$       | [s]    | Notch denominator coefficient        | 1.013
$A_2$       | [s²]   | Notch denominator coefficient        | 0.013
$A_3$       | [s]    | Notch denominator coefficient        | 0.0
$A_4$       | [s²]   | Notch denominator coefficient        | 0.0
$A_5$       | [s]    | Notch numerator coefficient          | 1.013
$A_6$       | [s²]   | Notch numerator coefficient          | 0.113
$T_1$       | [s]    | Lead–lag 1 numerator time constant   | 0.0
$T_2$       | [s]    | Lead–lag 1 denominator time constant | 0.02
$T_3$       | [s]    | Lead–lag 2 numerator time constant   | 0.0
$T_4$       | [s]    | Lead–lag 2 denominator time constant | 0.0
$T_5$       | [s]    | Washout numerator time constant      | 1.65
$T_6$       | [s]    | Washout denominator time constant    | 1.65
$K_s$       | [p.u.] | Stabilizer gain                      | 3.0
$L_s^{\min}$ | [p.u.] | Minimum stabilizer output limit      | -0.1
$L_s^{\max}$ | [p.u.] | Maximum stabilizer output limit      | 0.1

`Vcl` and `Vcu` are the lower and upper compensated-voltage cutout
thresholds in per unit. Zero disables the corresponding threshold. An
attached `vct` input is required when either threshold is enabled.
The PSS/E model has no transport delay; nonzero `Tdelay` (the PSLF
extension shown in the [PowerWorld diagram](https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Stabilizer%20IEEEST.htm))
is rejected, rather than ignored.

The parameter names above are the JSON keys, with `Ks`, `Lsmin`, and
`Lsmax` for the gain and limits. Defaults match PhasorDynamics: all `A`
coefficients, `T1`, `T3`, and `T5` are zero; `T2`, `T4`, `T6`, and `Ks`
are one; output limits are −0.1 and 0.1; cutout thresholds and delay are
zero. Parameters must be finite, time constants nonnegative, and limits
ordered.

### Derived Parameters

```math
\begin{aligned}
a_0 &= 1 \\
a_1 &= A_1 + A_3 \\
a_2 &= A_2 + A_4 + A_1 A_3 \\
a_3 &= A_1 A_4 + A_2 A_3 \\
a_4 &= A_2 A_4
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol                | Units  | Description
----------------------|--------|------------
$x_1, x_2, x_3, x_4$  | [-]    | Notch filter states
$x_5$                 | [-]    | Lead–lag 1 state
$x_6$                 | [-]    | Lead–lag 2 state
$x_7$                 | [-]    | Washout state

#### Algebraic

Symbol     | Units  | Description
-----------|--------|------------
$v_4$      | [p.u.] | Notch filter output
$v_5$      | [p.u.] | Lead–lag 1 output
$v_6$      | [p.u.] | Lead–lag 2 output
$v_7$      | [p.u.] | Unlimited stabilizer signal
$V_{ss}$   | [p.u.] | Limited stabilizer signal (model output)

### External Variables

#### Algebraic

Symbol | Units  | Description
-------|--------|------------
$u$    | [p.u.] | `input`, or `speed` minus one
$V_{ct}$ | [p.u.] | Compensated voltage magnitude for the cutout

## Model Equations

### Differential Equations

```math
\begin{aligned}
0 &= -\dot{x}_1 + x_2 \\
0 &= -\dot{x}_2 + x_3 \\
0 &= -\dot{x}_3 + x_4 \\
0 &= -\dot{x}_4 - \dfrac{a_0}{a_4}x_1 - \dfrac{a_1}{a_4}x_2 - \dfrac{a_2}{a_4}x_3 - \dfrac{a_3}{a_4}x_4 + \dfrac{1}{a_4}u \\
0 &= -T_2 \dot{x}_5 - x_5 + v_4 \\
0 &= -T_4 \dot{x}_6 - x_6 + v_5 \\
0 &= -T_6 \dot{x}_7 - x_7 + v_6
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
0 &= -v_4 + x_1 + A_5 x_2 + A_6 x_3 \\
0 &= -T_2(v_5 - x_5) + T_1(v_4 - x_5) \\
0 &= -T_4(v_6 - x_6) + T_3(v_5 - x_6) \\
0 &= -T_6 v_7 + K_s T_5(v_6 - x_7) \\
0 &= -V_{ss} + g(V_{ct})\,\text{clamp}(v_7, L_s^{\min}, L_s^{\max})
\end{aligned}
```

The output limiter uses GridKit's smooth
[Clamp](../../../../../CommonMath.md#derived-functions). With enabled
thresholds, the voltage gate is

```math
g(V_{ct}) = \sigma(V_{ct}-V_{cl})\,\sigma(V_{cu}-V_{ct}).
```

Each disabled threshold contributes a factor of one. `sigmoid` uses the
configured CommonMath `mu`; the gate tends to the voltage-window relay
as `mu` increases and is one half at an isolated enabled threshold.

The second- and third-order notch reductions retain the complete
numerator. With all denominator coefficients zero the notch is bypassed;
a lone nonzero first-order denominator is rejected, matching the
PhasorDynamics implementation. `T2=0` or `T4=0` bypasses that lead–lag;
`T6=0` gives the direct gain `Ks`. Unused notch states retain zero
initial derivatives. No time-constant floors are applied.

## Initialization

The current input initializes the filter and lag states at equilibrium.
For a washout with `T6>0`, the stabilizer output starts at zero; with
`T6=0` it starts at the limited and voltage-gated `Ks*u`.
The machine initializes before IEEEST, and exciters initialize afterward
so their reference incorporates the initial stabilizer output.

## Inputs and Outputs

Exactly one of `input` and `speed` must be connected. `input` is a generic
per-unit stabilizing signal; `speed` is absolute per-unit rotor speed,
with one subtracted inside the model. `vct` is optional when both voltage
cutouts are disabled. Required `output` connects to an exciter's `vs`.
The `vss` monitor reports the final output including the cutout.
