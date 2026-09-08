# IEEEST Model

Standard IEEE power system stabilizer: 4th-order notch filter, two lead–lag
blocks, washout, and output limiter.

## Block Diagram

![IEEEST model block diagram](../../../../../../docs/Figures/stabilizer_ieeest_diagram.png)

Figure 1: Stabilizer IEEEST model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$A_1,A_3$ | [s] | `A1`, `A3` | Notch denominator coefficients | Default $0$
$A_2,A_4$ | [s$^2$] | `A2`, `A4` | Notch denominator coefficients | Default $0$
$A_5$ | [s] | `A5` | Notch numerator coefficient | Default $0$
$A_6$ | [s$^2$] | `A6` | Notch numerator coefficient | Default $0$
$T_1$ | [s] | `T1` | First lead-lag numerator time constant | Default $0$
$T_2$ | [s] | `T2` | First lead-lag denominator time constant | Default $1$
$T_3$ | [s] | `T3` | Second lead-lag numerator time constant | Default $0$
$T_4$ | [s] | `T4` | Second lead-lag denominator time constant | Default $1$
$T_5$ | [s] | `T5` | Washout numerator time constant | Default $0$
$T_6$ | [s] | `T6` | Washout denominator time constant | Default $1$
$K_s$ | [p.u.] | `Ks` | Stabilizer gain | Default $1$
$L_s^{\min}$ | [p.u.] | `Lsmin` | Minimum stabilizer output | Default $-0.1$
$L_s^{\max}$ | [p.u.] | `Lsmax` | Maximum stabilizer output | Default $0.1$
$V_{\mathrm{cl}}$ | [p.u.] | `Vcl` | Lower voltage-cutout threshold | Default $0$; zero disables
$V_{\mathrm{cu}}$ | [p.u.] | `Vcu` | Upper voltage-cutout threshold | Default $0$; zero disables
$T_{\mathrm{delay}}$ | [s] | `Tdelay` | PSLF transport-delay extension | Only $0$ supported

### Parameter Validation

All parameters and derived coefficients must be finite.

```math
\begin{aligned}
T_1,T_2,T_3,T_4,T_5,T_6 &\ge 0 \\
L_s^{\min} &< L_s^{\max} \\
V_{\mathrm{cl}},V_{\mathrm{cu}} &\ge 0 \\
V_{\mathrm{cl}} &< V_{\mathrm{cu}},\quad\text{when both are enabled} \\
T_{\mathrm{delay}} &= 0
\end{aligned}
```

A lone first-order notch denominator is rejected. The PSS/E model has no
transport delay; the PSLF extension shown in the diagram is unsupported.

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

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$u$ | `input` | Input | [p.u.] | Stabilizing signal | Alternative to `speed`
$\omega_r$ | `speed` | Input | [p.u.] | Machine rotor speed | Alternative to `input`; $u=\omega_r-1$
$V_{ct}$ | `vct` | Input | [p.u.] | Compensated voltage magnitude | Required when either cutout is enabled
$V_{ss}$ | `output` | Output | [p.u.] | Stabilizer signal | Required; connects to exciter `vs`

Exactly one of `input` and `speed` must be connected. All attached signals
must be linked.

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$x_1$ | [p.u.] | Notch filter state |
$x_2$ | [p.u./s] | Notch filter state |
$x_3$ | [p.u./s$^2$] | Notch filter state |
$x_4$ | [p.u./s$^3$] | Notch filter state |
$x_5$ | [p.u.] | First lead-lag state | Algebraic when $T_2=0$
$x_6$ | [p.u.] | Second lead-lag state | Algebraic when $T_4=0$
$x_7$ | [p.u.] | Washout state | Algebraic when $T_6=0$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$v_4$      | [p.u.] | Notch filter output |
$v_5$      | [p.u.] | Lead–lag 1 output |
$v_6$      | [p.u.] | Lead–lag 2 output |
$v_7$      | [p.u.] | Unlimited stabilizer signal |
$V_{ss}$   | [p.u.] | Limited stabilizer signal (model output) |

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$u$ | [p.u.] | Stabilizing input | Connected `input`
$\omega_r$ | [p.u.] | Machine rotor speed | Connected `speed`
$V_{ct}$ | [p.u.] | Compensated voltage magnitude | Optional when both cutouts are disabled

## Model Equations

### Internal Equations

#### Differential

For a fourth-order notch and positive $T_2,T_4,T_6$:


```math
\begin{aligned}
0 &= -\dfrac{\mathrm{d}x_1}{\mathrm{d}t} + x_2 \\
0 &= -\dfrac{\mathrm{d}x_2}{\mathrm{d}t} + x_3 \\
0 &= -\dfrac{\mathrm{d}x_3}{\mathrm{d}t} + x_4 \\
0 &= -\dfrac{\mathrm{d}x_4}{\mathrm{d}t} - \dfrac{a_0}{a_4}x_1 - \dfrac{a_1}{a_4}x_2 - \dfrac{a_2}{a_4}x_3 - \dfrac{a_3}{a_4}x_4 + \dfrac{1}{a_4}u \\
0 &= -T_2 \dfrac{\mathrm{d}x_5}{\mathrm{d}t} - x_5 + v_4 \\
0 &= -T_4 \dfrac{\mathrm{d}x_6}{\mathrm{d}t} - x_6 + v_5 \\
0 &= -T_6 \dfrac{\mathrm{d}x_7}{\mathrm{d}t} - x_7 + v_6
\end{aligned}
```

#### Algebraic

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
g(V_{ct}) = \sigma(V_{ct}-V_{\mathrm{cl}})\,\sigma(V_{\mathrm{cu}}-V_{ct}).
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

### External Equations

None.

## Initialization

The current input initializes the filter and lag states at equilibrium.
For a washout with `T6>0`, the stabilizer output starts at zero; with
`T6=0` it starts at the limited and voltage-gated `Ks*u`.
The machine initializes before IEEEST, and exciters initialize afterward
so their reference incorporates the initial stabilizer output.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`vss` | [p.u.] | Stabilizer output | Includes voltage cutout
