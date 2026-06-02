# **IEEE Type DC2A Excitation System Model (ESDC2A)**

ESDC2A is an IEEE Type DC excitation system with a voltage transducer, lead-lag
input compensation, high-value under-excitation limiter selection, limited
voltage regulator, exciter feedback, saturation, and optional speed multiplier.

Notes:
- Internal voltage signals are on model base unless otherwise stated.
- The diagram labels the optional multiplier input as `Speed`; GridKit uses
  machine speed deviation, so the enabled multiplier is $1+\omega$.
- The PowerWorld selector `UEL` routes $V_{\mathrm{uel}}$ through the summing
  junction when `UEL >= 2`, and through the high-value gate when `UEL < 2`.
- The `exclim` flag lower-limits the exciter feedback signal at zero when
  nonzero; otherwise the feedback signal is unlimited.

## Block Diagram

Standard model of the ESDC2A Exciter.

<div align="center">
   <img align="center" src="../../../../../docs/Figures/PhasorDynamics/ESDC2A_diagram.png">

  Figure 1: Exciter ESDC2A model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)
</div>

## Model Parameters

Symbol                              | Units     | JSON      | Description                                             | Typical Value | Note
------------------------------------|-----------|-----------|---------------------------------------------------------|---------------|------
$T_R$                               | [sec]     | `Tr`      | Transducer time constant                                | 0.0           | Block name: `Tr`; if zero, $V_C$ is algebraic
$K_A$                               | [p.u.]    | `Ka`      | Voltage-regulator gain                                  | 40.0          | Block name: `Ka`
$T_A$                               | [sec]     | `Ta`      | Voltage-regulator time constant                         | 0.1           | Block name: `Ta`
$T_B$                               | [sec]     | `Tb`      | Lag time constant for voltage-regulator input lead-lag  | 0.0           | Block name: `Tb`; if $T_B=T_C=0$, the lead-lag block is bypassed
$T_C$                               | [sec]     | `Tc`      | Lead time constant for voltage-regulator input lead-lag | 0.0           | Block name: `Tc`; must be zero when $T_B=0$
$V_R^{\max}$                        | [p.u.]    | `Vrmax`   | Maximum voltage-regulator output                        | 1.0           | Block name: `Vrmax`
$V_R^{\min}$                        | [p.u.]    | `Vrmin`   | Minimum voltage-regulator output                        | -1.0          | Block name: `Vrmin`
$K_E$                               | [p.u.]    | `Ke`      | Exciter field-resistance line-slope margin              | 0.1           | Block name: `Ke`
$T_E$                               | [sec]     | `Te`      | Exciter time constant                                   | 0.5           | Block name: `Te`
$K_F$                               | [p.u.]    | `Kf`      | Stabilizing feedback gain                               | 0.05          | Block name: `Kf`
$T_{F1}$                            | [sec]     | `Tf1`     | Feedback lead time constant                             | 0.7           | Block name: `Tf1`
$s_{\mathrm{spd}}$                  | [binary]  | `Spdmlt`  | Speed multiplier flag                                   | 0             | Block name: `Spdmlt`; 1 enables the speed multiplier
$E_1$                               | [p.u.]    | `E1`      | First saturation voltage point                          | 2.8           | Block name: `E1`
$S_E(E_1)$                          | [p.u.]    | `SE1`     | Saturation value at $E_1$                               | 0.08          | Block name: `Se1`
$E_2$                               | [p.u.]    | `E2`      | Second saturation voltage point                         | 3.7           | Block name: `E2`
$S_E(E_2)$                          | [p.u.]    | `SE2`     | Saturation value at $E_2$                               | 0.33          | Block name: `Se2`
$I_{\mathrm{uel}}$                  | [integer] | `UEL`     | Under-excitation limiter input-location selector        | 0             | Block name: `UEL`; 0/1 = HV gate input, 2/3 = input-error summing junction
$s_{\mathrm{lim}}$                  | [binary]  | `exclim`  | Exciter feedback lower-limit flag                       | 1             | Block name: `exclim`; nonzero enables the zero lower limit on $V_{\mathrm{fe}}$

### Parameter Validation

Invalid ESDC2A parameter sets are rejected by the following checks. Source data
may apply PowerWorld-style autocorrections before these equations are evaluated.

```math
\begin{aligned}
  &K_A > 0 \\
  &T_R \ge 0,\quad T_A > 0,\quad T_B \ge 0,\quad T_C \ge 0,\quad T_E > 0,\quad T_{F1} \ge 0 \\
  &T_B > 0\quad\text{or}\quad(T_B = 0\ \text{and}\ T_C = 0) \\
  &V_R^{\min} \le V_R^{\max} \\
  &s_{\mathrm{spd}}, s_{\mathrm{lim}} \in \{0,1\} \\
  &I_{\mathrm{uel}} \in \{0,1,2,3\}
\end{aligned}
```

The saturation points are either disabled together,

```math
\begin{aligned}
  S_E(E_1) = 0,\quad S_E(E_2) = 0
\end{aligned}
```

or define a valid two-point quadratic saturation fit:

```math
\begin{aligned}
  &E_1 > 0,\quad E_2 > 0,\quad E_1 \ne E_2 \\
  &S_E(E_1) > 0,\quad S_E(E_2) > 0,\quad S_E(E_1) \ne S_E(E_2)
\end{aligned}
```

### Model Derived Parameters

The UEL routing flag and off-mode flag complements are:

```math
\begin{aligned}
  s_{\mathrm{uel}} &=
    \begin{cases}
      1 & I_{\mathrm{uel}} \ge 2 \\
      0 & I_{\mathrm{uel}} < 2
    \end{cases} \\
  s_{\mathrm{uel}}^{\mathrm{off}} &= 1 - s_{\mathrm{uel}} \\
  s_{\mathrm{lim}}^{\mathrm{off}} &= 1 - s_{\mathrm{lim}}
\end{aligned}
```

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
$E_{\mathrm{fd}}'$                  | [p.u.] | Field-voltage state before optional speed multiplier    | State 1 in Fig. 1; source label: `EFD`
$V_C$                               | [p.u.] | Sensed compensated voltage                              | State 2 in Fig. 1; source label: `Sensed Vt`; algebraic when $T_R=0$
$V_R$                               | [p.u.] | Voltage-regulator output                                | State 3 in Fig. 1; source label: `VR`
$V_F$                               | [p.u.] | Stabilizing feedback washout output                     | State 4 in Fig. 1; source label: `VF`; algebraic when $T_{F1}=0$
$x_{\mathrm{ll}}$                   | [p.u.] | Lead-lag block state                                    | State 5 in Fig. 1; source label: `Lead-Lag`

#### Algebraic

Symbol                              | Units  | Description                                             | Note
------------------------------------|--------|---------------------------------------------------------|------
$e_V$                               | [p.u.] | Voltage-regulator input error before lead-lag block     | Includes selected $V_{\mathrm{uel}}$ summing-junction input
$V_{\mathrm{ll}}$                   | [p.u.] | Lead-lag block output                                   | Input to high-value gate
$V_{\mathrm{hv}}$                   | [p.u.] | High-value gate output                                  | Selects $V_{\mathrm{ll}}$ or alternate $V_{\mathrm{uel}}$
$S_E$                               | [p.u.] | Saturation coefficient evaluated at $E_{\mathrm{fd}}'$  | Uses derived saturation curve
$V_{\mathrm{fe}}$                   | [p.u.] | Exciter feedback signal after optional lower limit      | Lower limited at zero when $s_{\mathrm{lim}}=1$
$E_{\mathrm{fd}}$                   | [p.u.] | Field-voltage output                                    | Output after optional speed multiplier

### External Variables

#### Differential

None.

#### Algebraic

Symbol                              | Units  | Description                                             | Note
------------------------------------|--------|---------------------------------------------------------|------
$E_C$                               | [p.u.] | Compensated terminal voltage magnitude                  | Source label: `EC`
$V_{\mathrm{ref}}$                  | [p.u.] | Voltage-control reference                               | Source label: `VREF`
$V_S$                               | [p.u.] | Stabilizer input signal                                 | Source label: `VS`; optional, defaults to zero
$V_{\mathrm{uel}}$                  | [p.u.] | Under-excitation limiter input                          | Source label: `VUEL`; optional, defaults to zero
$\omega$                            | [p.u.] | Machine speed deviation                                 | Source label: `Speed`; optional when $s_{\mathrm{spd}}=0$

## Model Equations

### Differential Equations

```math
\begin{aligned}
  0 &= -T_R\dot V_C - V_C + E_C \\
  0 &= -T_B\dot x_{\mathrm{ll}} - x_{\mathrm{ll}} + e_V \\
  0 &= -T_A\dot V_R
       + \text{antiwindup}\!\left(
          V_R,
          -V_R + K_A V_{\mathrm{hv}},
          V_R^{\min},
          V_R^{\max}
       \right) \\
  0 &= -T_E\dot E_{\mathrm{fd}}' + V_R - V_{\mathrm{fe}} \\
  0 &= -T_E T_{F1}\dot V_F - T_E V_F + K_F(V_R - V_{\mathrm{fe}})
\end{aligned}
```

CommonMath defines the [Anti-Windup](../../../../CommonMath.md#anti-windup-indicator)
target and smooth approximation.

### Algebraic Equations

```math
\begin{aligned}
  0 &= -e_V + V_{\mathrm{ref}} + V_S + s_{\mathrm{uel}}V_{\mathrm{uel}} - V_C - V_F \\
  0 &= -T_B(V_{\mathrm{ll}} - x_{\mathrm{ll}}) + T_C(e_V - x_{\mathrm{ll}}) \\
  0 &= -V_{\mathrm{hv}}
       + s_{\mathrm{uel}}V_{\mathrm{ll}}
       + s_{\mathrm{uel}}^{\mathrm{off}}\text{max}(V_{\mathrm{ll}}, V_{\mathrm{uel}}) \\
  0 &= -S_E + S_B\,q(E_{\mathrm{fd}}' - S_A) \\
  0 &= -V_{\mathrm{fe}}
       + s_{\mathrm{lim}}^{\mathrm{off}}(K_E + S_E)E_{\mathrm{fd}}'
       + s_{\mathrm{lim}}\rho\!\left((K_E + S_E)E_{\mathrm{fd}}'\right) \\
  0 &= -E_{\mathrm{fd}} + \left(1 + s_{\mathrm{spd}}\omega\right)E_{\mathrm{fd}}'
\end{aligned}
```

CommonMath defines the helper targets and smooth approximations for
[max](../../../../CommonMath.md#derived-functions) and the primitives
[ramp and quadratic ramp](../../../../CommonMath.md#primitives) $\rho$ and $q$.
When $T_B=T_C=0$, the lead-lag block is bypassed so $V_{\mathrm{ll}}=e_V$.

## Initialization

The machine initializes $E_{\mathrm{fd}}$ first. For a standard unsaturated
start, ESDC2A reads that value along with attached $\omega$, $E_C$, $V_S$, and
$V_{\mathrm{uel}}$, sets all internal derivatives to zero, and evaluates:

```math
\begin{aligned}
  E_{\mathrm{fd},0}' &= \dfrac{E_{\mathrm{fd},0}}{1 + s_{\mathrm{spd}}\omega_0} \\
  S_{E,0} &= S_B\,q(E_{\mathrm{fd},0}' - S_A) \\
  V_{\mathrm{fe},0}
    &= s_{\mathrm{lim}}^{\mathrm{off}}(K_E + S_{E,0})E_{\mathrm{fd},0}'
       + s_{\mathrm{lim}}\rho\!\left((K_E + S_{E,0})E_{\mathrm{fd},0}'\right) \\
  V_{R,0} &= V_{\mathrm{fe},0} \\
  V_{\mathrm{hv},0} &= \dfrac{V_{R,0}}{K_A} \\
  V_{C,0} &= E_{C,0} \\
  V_{F,0} &= 0 \\
  x_{\mathrm{ll},0} &= V_{\mathrm{ll},0} = e_{V,0} = V_{\mathrm{hv},0} \\
  V_{\mathrm{ref},0}
    &= e_{V,0} + V_{C,0} + V_{F,0} - V_{S,0} - s_{\mathrm{uel}}V_{\mathrm{uel},0}
\end{aligned}
```

This closed-form start requires $1 + s_{\mathrm{spd}}\omega_0 \ne 0$,
$V_R^{\min} \le V_{R,0} \le V_R^{\max}$, and, when $s_{\mathrm{uel}}=0$,
$V_{\mathrm{hv},0} \ge V_{\mathrm{uel},0}$. Saturated voltage-regulator starts
and active high-value-gate starts are outside these closed-form equations.

## Model Outputs

Output          | Units  | Description                         | Note
----------------|--------|-------------------------------------|------
`efd`           | [p.u.] | Field-voltage output                | $E_{\mathrm{fd}}$
`vc`            | [p.u.] | Sensed compensated voltage          | $V_C$
`vr`            | [p.u.] | Voltage-regulator output            | $V_R$
`vf`            | [p.u.] | Stabilizing feedback state          | $V_F$
`se`            | [p.u.] | Saturation coefficient              | $S_E$
`vfe`           | [p.u.] | Exciter feedback signal             | $V_{\mathrm{fe}}$
