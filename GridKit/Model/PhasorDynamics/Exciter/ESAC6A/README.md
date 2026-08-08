# **IEEE Type AC6A Excitation System Model (ESAC6A)**

ESAC6A is an IEEE Type AC excitation system with sensed terminal-voltage
feedback, cascaded regulator lead-lag blocks, voltage-regulator limits, an
exciter alternator state, field-current feedback limiting, rectifier loading,
saturation, and optional speed multiplier.

Notes:
- Internal voltage and current signals are on model base unless otherwise stated.
- The rectifier loading block $F_{\mathrm{ex}}=f(I_N)$ is the source AC-exciter
  loading curve from Fig. 1; it is not a CommonMath helper.
- The source diagram labels the optional multiplier input as `Speed`; GridKit
  uses machine speed deviation, so the enabled multiplier is $1+\omega$.

## Block Diagram

Standard model of the ESAC6A Exciter.

![](diagram.png)

Figure 1: Exciter ESAC6A model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol                              | Units    | JSON      | Description                                             | Typical Value | Note
------------------------------------|----------|-----------|---------------------------------------------------------|---------------|------
$T_R$                               | [sec]    | `Tr`      | Transducer time constant                                | 0.0           | Block name: `Tr`; if zero, $V_C$ is algebraic
$K_A$                               | [p.u.]   | `Ka`      | Voltage-regulator gain                                  | 40.0          | Block name: `Ka`
$T_A$                               | [sec]    | `Ta`      | Regulator denominator time constant                     | 0.1           | Block name: `Ta`
$T_K$                               | [sec]    | `Tk`      | Regulator numerator time constant                       | 0.0           | Block name: `Tk`
$T_B$                               | [sec]    | `Tb`      | Lag time constant for second lead-lag block             | 0.0           | Block name: `Tb`
$T_C$                               | [sec]    | `Tc`      | Lead time constant for second lead-lag block            | 0.0           | Block name: `Tc`
$V_A^{\max}$                        | [p.u.]   | `VaMax`   | Maximum first regulator block output                    | 1.0           | Block name: `VAMAX`
$V_A^{\min}$                        | [p.u.]   | `VaMin`   | Minimum first regulator block output                    | -1.0          | Block name: `VAMIN`
$V_R^{\max}$                        | [p.u.]   | `Vrmax`   | Maximum voltage-regulator output                        | 1.0           | Block name: `Vrmax`
$V_R^{\min}$                        | [p.u.]   | `Vrmin`   | Minimum voltage-regulator output                        | -1.0          | Block name: `Vrmin`
$T_E$                               | [sec]    | `Te`      | Exciter alternator time constant                        | 0.5           | Block name: `Te`
$V_{\mathrm{fe}}^{\mathrm{lim}}$    | [p.u.]   | `Vfelim`  | Feedback-limiter summing-junction reference             | 0.0           | Source label: `VFELIM`
$K_H$                               | [p.u.]   | `Kh`      | Feedback-limiter gain                                   | 1.0           | Block name: `KH`
$V_H^{\max}$                        | [p.u.]   | `Vhmax`   | Maximum feedback-limiter lead-lag output                | 1.0           | Block name: `VHMAX`; lower limit is zero
$T_H$                               | [sec]    | `Th`      | Feedback-limiter denominator time constant              | 0.0           | Block name: `TH`
$T_J$                               | [sec]    | `Tj`      | Feedback-limiter numerator time constant                | 0.0           | Block name: `TJ`
$K_C$                               | [p.u.]   | `Kc`      | Rectifier loading current coefficient                   | 0.0           | Block name: `Kc`; forms $I_N$
$K_D$                               | [p.u.]   | `Kd`      | Demagnetizing factor feedback gain                      | 0.0           | Block name: `Kd`
$K_E$                               | [p.u.]   | `Ke`      | Exciter field-resistance line-slope margin              | 0.1           | Block name: `Ke`
$E_1$                               | [p.u.]   | `E1`      | First saturation voltage point                          | 2.8           | Block name: `E1`
$S_E(E_1)$                          | [p.u.]   | `SE1`     | Saturation value at $E_1$                               | 0.08          | Block name: `Se1`
$E_2$                               | [p.u.]   | `E2`      | Second saturation voltage point                         | 3.7           | Block name: `E2`
$S_E(E_2)$                          | [p.u.]   | `SE2`     | Saturation value at $E_2$                               | 0.33          | Block name: `Se2`
$s_{\mathrm{spd}}$                  | [binary] | `Spdmlt`  | Speed multiplier flag                                   | 0             | Block name: `Spdmlt`; 1 enables the speed multiplier

### Parameter Validation

Invalid ESAC6A parameter sets are rejected by the following checks.

```math
\begin{aligned}
  &K_A > 0 \\
  &T_R \ge 0,\quad T_A > 0,\quad T_K \ge 0,\quad T_B \ge 0,\quad T_C \ge 0,\quad T_E > 0 \\
  &T_H \ge 0,\quad T_J \ge 0 \\
  &T_B > 0\quad\text{or}\quad(T_B = 0\ \text{and}\ T_C = 0) \\
  &T_H > 0\quad\text{or}\quad(T_H = 0\ \text{and}\ T_J = 0) \\
  &V_A^{\min} \le V_A^{\max},\quad V_R^{\min} \le V_R^{\max},\quad V_H^{\max} \ge 0 \\
  &s_{\mathrm{spd}} \in \{0,1\}
\end{aligned}
```

The saturation points are either disabled together or define a valid positive
two-point quadratic fit.

### Model Derived Parameters

The saturation curve is fitted from the two supplied saturation points. If both
saturation factors are zero, use $S_A=0$ and $S_B=0$. Otherwise:

```math
\begin{aligned}
  C &= \sqrt{\dfrac{S_E(E_2)}{S_E(E_1)}} \\
  S_A &= \dfrac{C E_1 - E_2}{C - 1} \\
  S_B &= \dfrac{S_E(E_1)}{(E_1 - S_A)^2}
\end{aligned}
```

## Model Variables

### Internal Variables

#### Differential

Symbol                              | Units  | Description                                             | Note
------------------------------------|--------|---------------------------------------------------------|------
$V_E$                               | [p.u.] | Exciter alternator voltage state before output multipliers | State 1 in Fig. 1; source label: `VE`
$V_C$                               | [p.u.] | Sensed compensated voltage                              | State 2 in Fig. 1; source label: `Sensed Vt`; algebraic when $T_R=0$
$x_A$                               | [p.u.] | First regulator lead-lag denominator state              | State 3 in Fig. 1; source label: `TA Block`
$x_{\mathrm{ll}}$                   | [p.u.] | Second lead-lag denominator state                       | State 4 in Fig. 1; source label: `VLL`
$V_F$                               | [p.u.] | Stabilizing feedback signal                             | State 5 in Fig. 1; source label: `VF`

#### Algebraic

Symbol                              | Units  | Description                                             | Note
------------------------------------|--------|---------------------------------------------------------|------
$e_V$                               | [p.u.] | Voltage-regulator input error before first lead-lag     | Summing junction after sensed voltage
$V_A$                               | [p.u.] | Limited first regulator lead-lag output                 | Limited by $V_A^{\min}$ and $V_A^{\max}$
$V_{\mathrm{ll}}$                   | [p.u.] | Second lead-lag output                                  | Input to $V_R$ summing junction
$V_H$                               | [p.u.] | Feedback-limiter lead-lag output before $K_H$           | Limited by 0 and $V_H^{\max}$
$V_H^{\mathrm{pre}}$                | [p.u.] | Feedback-limiter lead-lag output before limits          | Bypasses to $V_F$ when $T_H=T_J=0$
$V_R$                               | [p.u.] | Voltage-regulator output                                | Limited by $V_R^{\min}$ and $V_R^{\max}$
$S_E$                               | [p.u.] | Saturation coefficient evaluated at $V_E$               | Uses derived saturation curve
$I_N$                               | [p.u.] | Normalized exciter loading current                      | Source label: `IN`; satisfies $V_E I_N=K_C I_{\mathrm{fd}}$
$F_{\mathrm{ex}}$                   | [p.u.] | Rectifier loading factor                                | Source label: `FEX`; source curve $F_{\mathrm{ex}}=f(I_N)$
$V_{\mathrm{fe}}$                   | [p.u.] | Exciter feedback signal                                 | Sum of saturation/resistance, $K_D I_{\mathrm{fd}}$, and feedback-limiter paths
$E_{\mathrm{fd}}$                   | [p.u.] | Field-voltage output                                    | Output after rectifier loading and optional speed multiplier

### External Variables

#### Differential

None.

#### Algebraic

Symbol                              | Units  | Description                                             | Note
------------------------------------|--------|---------------------------------------------------------|------
$E_C$                               | [p.u.] | Compensated terminal voltage magnitude                  | Source label: `EC`
$V_{\mathrm{ref}}$                  | [p.u.] | Voltage-control reference                               | Source label: `VREF`
$V_{\mathrm{uel}}$                  | [p.u.] | Under-excitation limiter input                          | Source label: `VUEL`; optional, defaults to zero
$V_S$                               | [p.u.] | Stabilizer input signal                                 | Source label: `VS`; optional, defaults to zero
$I_{\mathrm{fd}}$                   | [p.u.] | Machine field current                                   | Source label: `IFD`
$\omega$                            | [p.u.] | Machine speed deviation                                 | Source label: `Speed`; optional when $s_{\mathrm{spd}}=0$

## Model Equations

### Differential Equations

```math
\begin{aligned}
  0 &= -T_R\dot V_C - V_C + E_C \\
  0 &= -T_A\dot x_A - x_A + K_A e_V \\
  0 &= -T_B\dot x_{\mathrm{ll}} - x_{\mathrm{ll}} + V_A \\
  0 &= -T_E\dot V_E + V_R - V_{\mathrm{fe}} \\
  0 &= -T_H\dot V_F - V_F + V_{\mathrm{fe}}
\end{aligned}
```

### Algebraic Equations

```math
\begin{aligned}
  0 &= -e_V + V_{\mathrm{ref}} + V_{\mathrm{uel}} + V_S - V_C - V_F \\
  0 &= -V_A
       + \text{clamp}\!\left(
          x_A + \dfrac{T_K}{T_A}(K_A e_V - x_A),
          V_A^{\min},
          V_A^{\max}
       \right) \\
  0 &= -T_B(V_{\mathrm{ll}} - x_{\mathrm{ll}}) + T_C(V_A - x_{\mathrm{ll}}) \\
  0 &= -V_H^{\mathrm{pre}}
       +
       \begin{cases}
         V_F, & T_H = T_J = 0 \\
         V_F + \dfrac{T_J}{T_H}(V_{\mathrm{fe}} - V_F), & T_H > 0
       \end{cases} \\
  0 &= -V_H + \text{clamp}\!\left(V_H^{\mathrm{pre}}, 0, V_H^{\max}\right) \\
  0 &= -V_R + \text{clamp}(V_{\mathrm{ll}} - K_H V_H, V_R^{\min}, V_R^{\max}) \\
  0 &= -S_E + S_B\,q(V_E - S_A) \\
  0 &= -V_E I_N + K_C I_{\mathrm{fd}} \\
  0 &= -F_{\mathrm{ex}} + f(I_N) \\
  0 &= -V_{\mathrm{fe}}
       + (K_E + S_E)V_E + K_D I_{\mathrm{fd}} + V_{\mathrm{fe}}^{\mathrm{lim}} \\
  0 &= -E_{\mathrm{fd}}
       + \left(1+s_{\mathrm{spd}}\omega\right)F_{\mathrm{ex}}V_E
\end{aligned}
```

CommonMath defines helper targets for [clamp](../../../../CommonMath.md#derived-functions)
and the primitive [quadratic ramp](../../../../CommonMath.md#primitives) $q$.
The rectifier loading function $f(I_N)$ is the source curve shown in Fig. 1.
When $T_B=T_C=0$, the second lead-lag block is bypassed. When $T_H=T_J=0$, the
feedback-limiter lead-lag block is bypassed before the 0-to-$V_H^{\max}$ clamp.

## Initialization

The machine initializes $E_{\mathrm{fd}}$ and $I_{\mathrm{fd}}$ first. For a
standard unsaturated start, ESAC6A reads those values and sets all internal
derivatives to zero. First solve the coupled rectifier-loading equations:

```math
\begin{aligned}
  0 &= -V_{E,0}I_{N,0} + K_C I_{\mathrm{fd},0} \\
  0 &= -F_{\mathrm{ex},0} + f(I_{N,0}) \\
  0 &= -E_{\mathrm{fd},0}
       + \left(1+s_{\mathrm{spd}}\omega_0\right)F_{\mathrm{ex},0}V_{E,0}
\end{aligned}
```

Then evaluate:

```math
\begin{aligned}
  S_{E,0} &= S_B\,q(V_{E,0} - S_A) \\
  V_{\mathrm{fe},0} &= (K_E + S_{E,0})V_{E,0} + K_D I_{\mathrm{fd},0} + V_{\mathrm{fe}}^{\mathrm{lim}} \\
  V_{F,0} &= V_{\mathrm{fe},0} \\
  V_{H,0}^{\mathrm{pre}} &= V_{F,0} \\
  V_{H,0} &= \text{clamp}(V_{H,0}^{\mathrm{pre}}, 0, V_H^{\max}) \\
  V_{R,0} &= V_{\mathrm{fe},0} \\
  V_{\mathrm{ll},0} &= V_{R,0} + K_H V_{H,0} \\
  V_{A,0} &= V_{\mathrm{ll},0} \\
  x_{A,0} &= V_{A,0} \\
  x_{\mathrm{ll},0} &= V_{\mathrm{ll},0} \\
  V_{C,0} &= E_{C,0} \\
  e_{V,0} &= \dfrac{x_{A,0}}{K_A} \\
  V_{\mathrm{ref},0} &= e_{V,0} + V_{C,0} + V_{F,0} - V_{\mathrm{uel},0} - V_{S,0}
\end{aligned}
```

This standard start requires $1+s_{\mathrm{spd}}\omega_0\ne 0$,
$V_{E,0}\ne 0$, inactive $V_A$, $V_R$, and $V_H$ limits, and nonsingular
regulator gains/time constants. Starts that bind those limits are outside
these closed-form equations.

## Model Outputs

Output          | Units  | Description                         | Note
----------------|--------|-------------------------------------|------
`efd`           | [p.u.] | Field-voltage output                | $E_{\mathrm{fd}}$
`ve`            | [p.u.] | Exciter alternator voltage state    | $V_E$
`vc`            | [p.u.] | Sensed compensated voltage          | $V_C$
`va`            | [p.u.] | First regulator output              | $V_A$
`vll`           | [p.u.] | Second lead-lag output              | $V_{\mathrm{ll}}$
`vf`            | [p.u.] | Feedback-limiter state              | $V_F$
`vh`            | [p.u.] | Feedback-limiter output             | $V_H$
`vr`            | [p.u.] | Voltage-regulator output            | $V_R$
`in`            | [p.u.] | Normalized exciter loading current  | $I_N$
`fex`           | [p.u.] | Rectifier loading factor            | $F_{\mathrm{ex}}$
`se`            | [p.u.] | Saturation coefficient              | $S_E$
