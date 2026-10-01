# ESDC2A

ESDC2A is an IEEE Type DC excitation system with a voltage transducer, lead–lag
input compensation, high-value under-excitation limiter selection, limited
voltage regulator, exciter feedback, saturation, and optional speed multiplier.

> [!WARNING]
> Initialization does not yet invert the smooth gates and limits used by the
> model equations.

## Notes

- Internal voltage signals are on model base unless otherwise stated.
- The diagram labels the optional multiplier input as `Speed`; GridKit uses
  machine speed deviation, so the enabled multiplier is $1+\omega$.
- The PowerWorld selector `UEL` routes $V_{\mathrm{uel}}$ through the summing
  junction when `UEL >= 2`, and through the high-value gate when `UEL < 2`.
- The `exclim` flag lower-limits the exciter feedback signal at zero when
  nonzero; otherwise the feedback signal is unlimited.

## Block Diagram

![](../../../../../docs/Figures/PhasorDynamics/ESDC2A_diagram.png)

Figure 1: Exciter ESDC2A model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol             | Units     | JSON     | Description                                             | Typical Value | Note
-------------------|-----------|----------|---------------------------------------------------------|---------------|----------------------------------------------------------
$T_R$              | [s]       | `Tr`     | Transducer time constant                                | 0.0           | if zero, $V_C$ is algebraic
$K_A$              | [p.u.]    | `Ka`     | Voltage-regulator gain                                  | 40.0          |
$T_A$              | [s]       | `Ta`     | Voltage-regulator time constant                         | 0.1           |
$T_B$              | [s]       | `Tb`     | Lag time constant for voltage-regulator input lead–lag  | 0.0           | if $T_B=T_C=0$, the lead–lag block is bypassed
$T_C$              | [s]       | `Tc`     | Lead time constant for voltage-regulator input lead–lag | 0.0           | must be zero when $T_B=0$
$V_R^{\max}$       | [p.u.]    | `Vrmax`  | Maximum voltage-regulator output                        | 1.0           |
$V_R^{\min}$       | [p.u.]    | `Vrmin`  | Minimum voltage-regulator output                        | -1.0          |
$K_E$              | [p.u.]    | `Ke`     | Exciter field-resistance line-slope margin              | 0.1           |
$T_E$              | [s]       | `Te`     | Exciter time constant                                   | 0.5           |
$K_F$              | [p.u.]    | `Kf`     | Stabilizing feedback gain                               | 0.05          |
$T_{F1}$           | [s]       | `Tf1`    | Feedback lead time constant                             | 0.7           |
$s_{\mathrm{spd}}$ | [binary]  | `Spdmlt` | Speed multiplier flag                                   | 0             | 1 enables the speed multiplier
$E_1$              | [p.u.]    | `E1`     | First saturation voltage point                          | 2.8           |
$S_E(E_1)$         | [p.u.]    | `SE1`    | Saturation value at $E_1$                               | 0.08          | Source label: `Se1`
$E_2$              | [p.u.]    | `E2`     | Second saturation voltage point                         | 3.7           |
$S_E(E_2)$         | [p.u.]    | `SE2`    | Saturation value at $E_2$                               | 0.33          | Source label: `Se2`
$I_{\mathrm{uel}}$ | [integer] | `UEL`    | Under-excitation limiter input-location selector        | 0             | 0/1 = HV gate input, 2/3 = input-error summing junction
$s_{\mathrm{lim}}$ | [binary]  | `exclim` | Exciter feedback lower-limit flag                       | 1             | nonzero enables the zero lower limit on $V_{\mathrm{FE}}$

### Parameter Validation

A valid ESDC2A parameter set must satisfy the following conditions:

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

## Model Ports

Name    | Port   | Init | Description
--------|--------|------|------------
`ec`    | Input  | TBD  | Compensated terminal voltage magnitude $E_C$
`vref`  | Input  | TBD  | Voltage-control reference $V_{\mathrm{ref}}$
`vs`    | Input  | TBD  | Stabilizer input signal $V_S$
`vuel`  | Input  | TBD  | Under-excitation limiter input $V_{\mathrm{uel}}$
`speed` | Input  | TBD  | Machine speed deviation $\omega$
`efd`   | Output | TBD  | Field-voltage output $E_{\mathrm{fd}}$

## Model Variables

### Internal Variables

#### Differential

Symbol             | Units  | Description                                          | Note
-------------------|--------|------------------------------------------------------|---------------------------------------------------------------------
$E_{\mathrm{fd}}'$ | [p.u.] | Field-voltage state before optional speed multiplier | State 1 in Fig. 1; Source label: `EFD`
$V_C$              | [p.u.] | Sensed compensated voltage                           | State 2 in Fig. 1; Source label: `Sensed Vt`; algebraic when $T_R=0$
$V_R$              | [p.u.] | Voltage-regulator output                             | State 3 in Fig. 1; Source label: `VR`
$V_F$              | [p.u.] | Stabilizing feedback washout output                  | State 4 in Fig. 1; Source label: `VF`; algebraic when $T_{F1}=0$
$x_{\mathrm{LL}}$  | [p.u.] | Lead–lag block state                                 | State 5 in Fig. 1; Source label: `Lead-Lag`

#### Algebraic

Symbol            | Units  | Description                                            | Note
------------------|--------|--------------------------------------------------------|------------------------------------------------------------
$e_V$             | [p.u.] | Voltage-regulator input error before lead–lag block    | Includes selected $V_{\mathrm{uel}}$ summing-junction input
$V_{\mathrm{LL}}$ | [p.u.] | Lead–lag block output                                  | Input to high-value gate
$V_{\mathrm{HV}}$ | [p.u.] | High-value gate output                                 | Selects $V_{\mathrm{LL}}$ or alternate $V_{\mathrm{uel}}$
$S_E$             | [p.u.] | Saturation coefficient evaluated at $E_{\mathrm{fd}}'$ | Uses derived saturation curve
$V_{\mathrm{FE}}$ | [p.u.] | Exciter feedback signal after optional lower limit     | Lower limited at zero when $s_{\mathrm{lim}}=1$
$E_{\mathrm{fd}}$ | [p.u.] | Field-voltage output                                   | Output after optional speed multiplier

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

Smooth functions: [`antiwindup`](../../../../CommonMath.md#antiwindup), [`max`](../../../../CommonMath.md#maximum), [$q$](../../../../CommonMath.md#quadratic-ramp), [$\rho$](../../../../CommonMath.md#ramp).

### Internal Equations

#### Differential

```math
\begin{aligned}
  0 &= -T_R\dot V_C - V_C + E_C \\
  0 &= -T_B\dot x_{\mathrm{LL}} - x_{\mathrm{LL}} + e_V \\
  0 &= -T_A\dot V_R
       + \text{antiwindup}\!(
          V_R,
          -V_R + K_A V_{\mathrm{HV}};
          V_R^{\min},
          V_R^{\max}
       ) \\
  0 &= -T_E\dot E_{\mathrm{fd}}' + V_R - V_{\mathrm{FE}} \\
  0 &= -T_E T_{F1}\dot V_F - T_E V_F + K_F(V_R - V_{\mathrm{FE}})
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
  0 &= -e_V + V_{\mathrm{ref}} + V_S + s_{\mathrm{uel}}V_{\mathrm{uel}} - V_C - V_F \\
  0 &= -T_B(V_{\mathrm{LL}} - x_{\mathrm{LL}}) + T_C(e_V - x_{\mathrm{LL}}) \\
  0 &= -V_{\mathrm{HV}}
       + s_{\mathrm{uel}}V_{\mathrm{LL}}
       + s_{\mathrm{uel}}^{\mathrm{off}}\text{max}(V_{\mathrm{LL}}, V_{\mathrm{uel}}) \\
  0 &= -S_E + S_B\,q(E_{\mathrm{fd}}' - S_A) \\
  0 &= -V_{\mathrm{FE}}
       + s_{\mathrm{lim}}^{\mathrm{off}}(K_E + S_E)E_{\mathrm{fd}}'
       + s_{\mathrm{lim}}\rho\!((K_E + S_E)E_{\mathrm{fd}}') \\
  0 &= -E_{\mathrm{fd}} + (1 + s_{\mathrm{spd}}\omega)E_{\mathrm{fd}}'
\end{aligned}
```

When $T_B=T_C=0$, the lead–lag block is bypassed so $V_{\mathrm{LL}}=e_V$.

### External Equations

None.

## Initialization

The machine initializes $E_{\mathrm{fd}}$ first. For a standard unsaturated
start, ESDC2A reads that value along with attached $\omega$, $E_C$, $V_S$, and
$V_{\mathrm{uel}}$, sets all internal derivatives to zero, and evaluates:

```math
\begin{aligned}
  E_{\mathrm{fd}}' &\leftarrow \dfrac{E_{\mathrm{fd}}}{1 + s_{\mathrm{spd}}\omega} \\
  S_E &\leftarrow S_B\,q(E_{\mathrm{fd}}' - S_A) \\
  V_{\mathrm{FE}}
    &\leftarrow s_{\mathrm{lim}}^{\mathrm{off}}(K_E + S_E)E_{\mathrm{fd}}'
       + s_{\mathrm{lim}}\rho\!((K_E + S_E)E_{\mathrm{fd}}') \\
  V_R &\leftarrow V_{\mathrm{FE}} \\
  V_{\mathrm{HV}} &\leftarrow \dfrac{V_R}{K_A} \\
  V_C &\leftarrow E_C \\
  V_F &\leftarrow 0 \\
  V_{\mathrm{LL}} &\leftarrow V_{\mathrm{HV}} \\
  e_V &\leftarrow V_{\mathrm{LL}} \\
  x_{\mathrm{LL}} &\leftarrow e_V \\
  V_{\mathrm{ref}}
    &\leftarrow e_V + V_C + V_F - V_S - s_{\mathrm{uel}}V_{\mathrm{uel}}
\end{aligned}
```

This closed-form start requires $1 + s_{\mathrm{spd}}\omega \ne 0$,
$V_R^{\min} \le V_R \le V_R^{\max}$, and, when $s_{\mathrm{uel}}=0$,
$V_{\mathrm{HV}} \ge V_{\mathrm{uel}}$. Saturated voltage-regulator starts
and active high-value-gate starts are outside these closed-form equations.

## Monitors

Monitor | Units  | Description                | Note
--------|--------|----------------------------|------------------
`efd`   | [p.u.] | Field-voltage output       | $E_{\mathrm{fd}}$
`vc`    | [p.u.] | Sensed compensated voltage | $V_C$
`vr`    | [p.u.] | Voltage-regulator output   | $V_R$
`vf`    | [p.u.] | Stabilizing feedback state | $V_F$
`se`    | [p.u.] | Saturation coefficient     | $S_E$
`vfe`   | [p.u.] | Exciter feedback signal    | $V_{\mathrm{FE}}$
