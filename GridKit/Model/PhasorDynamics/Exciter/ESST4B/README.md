# ESST4B

ESST4B is a static excitation system with compensated-voltage sensing, an outer
proportional/integral voltage regulator, a lag block, an inner
proportional/integral regulator with exciter-output feedback, low-value
over-excitation limiter gating, and potential- or compound-source rectifier
scaling.

> [!WARNING]
> Initialization does not yet invert the smooth gates and limits used by the
> model equations.

## Notes

- Internal voltage and current signals are on model base unless otherwise stated.
- The rectifier loading block $F_{\mathrm{ex}}=f(I_N)$ is the source
  controlled-rectifier loading curve from Fig. 1; it is not a CommonMath helper.

## Block Diagram

![](../../../../../docs/Figures/PhasorDynamics/ESST4B_diagram.png)

Figure 1: Exciter ESST4B model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol            | Units    | JSON        | Description                                           | Typical Value | Note
------------------|----------|-------------|-------------------------------------------------------|---------------|------------------------------------------------------------------------
$T_R$             | [s]      | `Tr`        | Compensated-voltage transducer time constant          | 0.0           | if zero, sensed voltage is algebraic
$K_{\mathrm{pr}}$ | [p.u.]   | `Kpr`       | Outer regulator proportional gain                     | 1.0           | Source label: `KPR`
$K_{\mathrm{ir}}$ | [p.u./s] | `Kir`       | Outer regulator integral gain                         | 0.0           | Source label: `KIR`
$V_R^{\max}$      | [p.u.]   | `Vrmax`     | Maximum outer regulator output                        | 1.0           | Source label: `VRMAX`
$V_R^{\min}$      | [p.u.]   | `Vrmin`     | Minimum outer regulator output                        | -1.0          | Source label: `VRMIN`
$T_A$             | [s]      | `Ta`        | Regulator lag time constant                           | 0.0           | if zero, $V_A$ is algebraic
$K_{\mathrm{pm}}$ | [p.u.]   | `Kpm`       | Inner regulator proportional gain                     | 1.0           | Source label: `KPM`
$K_{\mathrm{im}}$ | [p.u./s] | `Kim`       | Inner regulator integral gain                         | 0.0           | Source label: `KIM`
$V_M^{\max}$      | [p.u.]   | `VmMax`     | Maximum inner regulator output                        | 1.0           | Source label: `VMMAX`
$V_M^{\min}$      | [p.u.]   | `VmMin`     | Minimum inner regulator output                        | 0.0           | Source label: `VMMIN`
$K_G$             | [p.u.]   | `Kg`        | Exciter-output feedback gain into inner regulator     | 0.0           | Source label: `KG`
$K_P$             | [p.u.]   | `Kp`        | Potential-source voltage coefficient magnitude        | 0.0           | Source label: `KP`
$K_I$             | [p.u.]   | `Ki`        | Potential-source current coefficient                  | 0.0           | Source label: `KI`
$V_B^{\max}$      | [p.u.]   | `VbMax`     | Maximum rectifier source multiplier                   | 999.0         | Source label: `VBMAX`
$K_C$             | [p.u.]   | `Kc`        | Rectifier loading current coefficient                 | 0.0           | forms $I_N$
$X_L$             | [p.u.]   | `Xl`        | Source reactance term in potential-source calculation | 0.0           | Source label: `XL`
$\theta_P$        | [deg]    | `ThetaPDeg` | Potential-source coefficient angle                    | 0.0           | Source label: `thetaP`; forms $K_P^r$ and $K_P^i$
$V_G^{\max}$      | [p.u.]   | `VgMax`     | Maximum exciter-output feedback signal                | 999.0         | Source label: `VGMAX`; ceiling on $K_G E_{\mathrm{fd}}$

### Parameter Validation

A valid ESST4B parameter set must satisfy the following conditions:

```math
\begin{aligned}
  &T_R \ge 0,\quad T_A \ge 0 \\
  &V_R^{\min} \le V_R^{\max},\quad V_M^{\min} \le V_M^{\max} \\
  &V_B^{\max} > 0,\quad V_G^{\max} \ge 0
\end{aligned}
```

### Model Derived Parameters

The potential-source coefficient is resolved into real scalar components:

```math
\begin{aligned}
  K_P^r &= K_P\cos\theta_P \\
  K_P^i &= K_P\sin\theta_P
\end{aligned}
```

Here $\theta_P$ is converted from degrees before evaluating the trigonometric
functions.

## Model Ports

Name    | Port   | Init | Description
--------|--------|------|--------------------------------------------------
`vcomp` | Input  | TBD  | Compensated voltage input $V_{\mathrm{comp}}$
`vref`  | Input  | TBD  | Voltage-control reference $V_{\mathrm{ref}}$
`vuel`  | Input  | TBD  | Under-excitation limiter input $V_{\mathrm{uel}}$
`vs`    | Input  | TBD  | Stabilizer input signal $V_S$
`voel`  | Input  | TBD  | Over-excitation limiter input $V_{\mathrm{oel}}$
`vr`    | Input  | TBD  | Terminal-voltage real component $V_r$
`vi`    | Input  | TBD  | Terminal-voltage imaginary component $V_i$
`ir`    | Input  | TBD  | Terminal-current real component $I_r$
`ii`    | Input  | TBD  | Terminal-current imaginary component $I_i$
`ifd`   | Input  | TBD  | Machine field current $I_{\mathrm{fd}}$
`efd`   | Output | TBD  | Field-voltage output $E_{\mathrm{fd}}$

## Model Variables

### Internal Variables

#### Differential

Symbol | Units  | Description                    | Note
-------|--------|--------------------------------|---------------------------------------------------------------------
$V_M$  | [p.u.] | Inner regulator output         | State 1 in Fig. 1
$V_C$  | [p.u.] | Sensed compensated voltage     | State 2 in Fig. 1; Source label: `Sensed Vt`; algebraic when $T_R=0$
$V_A$  | [p.u.] | Lagged outer-regulator output  | State 3 in Fig. 1; algebraic when $T_A=0$
$x_R$  | [p.u.] | Outer regulator integral state | State 4 in Fig. 1; Source label: `VR`

#### Algebraic

Symbol                              | Units  | Description                                             | Note
------------------------------------|--------|---------------------------------------------------------|------
$e_V$                               | [p.u.] | Voltage-error signal into outer regulator               | Summing junction after sensed voltage
$V_R$                               | [p.u.] | Limited outer regulator output                          | Limited by $V_R^{\min}$ and $V_R^{\max}$
$V_G$                               | [p.u.] | Limited exciter-output feedback signal                  | $K_G E_{\mathrm{fd}}$ limited by $V_G^{\max}$
$e_M$                               | [p.u.] | Inner regulator error                                   | $V_A$ minus $V_G$
$V_{\mathrm{lv}}$                   | [p.u.] | Low-value gate output                                   | Lesser of $V_M$ and $V_{\mathrm{oel}}$
$V_{\mathrm{src}}^r$     | [p.u.] | Real component of the potential-source expression       | From terminal voltage/current components
$V_{\mathrm{src}}^i$     | [p.u.] | Imaginary component of the potential-source expression  | From terminal voltage/current components
$V_E$                               | [p.u.] | Potential- or compound-source voltage magnitude         | Nonnegative source magnitude
$I_N$                               | [p.u.] | Normalized exciter loading current                      | Source label: `IN`; satisfies $V_E I_N=K_C I_{\mathrm{fd}}$
$F_{\mathrm{ex}}$                   | [p.u.] | Rectifier loading factor                                | Source label: `FEX`; source curve $F_{\mathrm{ex}}=f(I_N)$
$V_B$                               | [p.u.] | Rectifier source multiplier                             | Limited by $V_B^{\max}$
$E_{\mathrm{fd}}$                   | [p.u.] | Field-voltage output                                    | Product of low-value gate and $V_B$

### External Variables

#### Differential

None.

#### Algebraic

Symbol              | Units  | Description                          | Note
--------------------|--------|--------------------------------------|----------------------------------------------------------------------
$V_{\mathrm{comp}}$ | [p.u.] | Compensated voltage input            | Source label: `VCOMP`
$V_{\mathrm{ref}}$  | [p.u.] | Voltage-control reference            | Source label: `VREF`
$V_{\mathrm{uel}}$  | [p.u.] | Under-excitation limiter input       | Source label: `VUEL`; optional, defaults to zero
$V_S$               | [p.u.] | Stabilizer input signal              | Source label: `VS`; optional, defaults to zero
$V_{\mathrm{oel}}$  | [p.u.] | Over-excitation limiter input        | Source label: `VOEL`; optional, defaults to a high value when omitted
$V_r$               | [p.u.] | Terminal-voltage real component      | Source label: `VT`
$V_i$               | [p.u.] | Terminal-voltage imaginary component | Source label: `VT`
$I_r$               | [p.u.] | Terminal-current real component      | Source label: `IT`
$I_i$               | [p.u.] | Terminal-current imaginary component | Source label: `IT`
$I_{\mathrm{fd}}$   | [p.u.] | Machine field current                | Source label: `IFD`

## Model Equations

Smooth functions: [`antiwindup`](../../../../CommonMath.md#antiwindup), [`clamp`](../../../../CommonMath.md#clamp), [`min`](../../../../CommonMath.md#minimum).

### Internal Equations

#### Differential

```math
\begin{aligned}
  0 &= -T_R\dot V_C - V_C + V_{\mathrm{comp}} \\
  0 &=
    -\dot x_R
    + \text{antiwindup}\!(
        V_R,
        K_{\mathrm{ir}}e_V;
        V_R^{\min},
        V_R^{\max}
      ) \\
  0 &= -T_A\dot V_A - V_A + V_R \\
  0 &=
    -\dot V_M
    + \text{antiwindup}\!(
        V_M,
        K_{\mathrm{im}}e_M;
        V_M^{\min},
        V_M^{\max}
      )
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
  0 &= -e_V + V_{\mathrm{ref}} + V_{\mathrm{uel}} + V_S - V_C \\
  0 &= -V_R + \text{clamp}(K_{\mathrm{pr}}e_V + x_R; V_R^{\min}, V_R^{\max}) \\
  0 &= -V_G + \text{min}(K_G E_{\mathrm{fd}}, V_G^{\max}) \\
  0 &= -e_M + V_A - V_G \\
  0 &= -V_{\mathrm{lv}} + \text{min}(V_M, V_{\mathrm{oel}}) \\
  0 &= -V_{\mathrm{src}}^r
       + K_P V_r
       - X_L K_P^i I_r
       - (K_I + X_L K_P^r)I_i \\
  0 &= -V_{\mathrm{src}}^i
       + K_P V_i
       + (K_I + X_L K_P^r)I_r
       - X_L K_P^i I_i \\
  0 &= -V_E^2
       + (V_{\mathrm{src}}^r)^2
       + (V_{\mathrm{src}}^i)^2 \\
  0 &= -V_E I_N + K_C I_{\mathrm{fd}} \\
  0 &= -F_{\mathrm{ex}} + f(I_N) \\
  0 &= -V_B + \text{min}(V_E F_{\mathrm{ex}}, V_B^{\max}) \\
  0 &= -E_{\mathrm{fd}} + V_{\mathrm{lv}}V_B
\end{aligned}
```

### External Equations

None.

## Initialization

For a standard unsaturated start, the machine initializes
$E_{\mathrm{fd}}$ and $I_{\mathrm{fd}}$ first. ESST4B reads those values,
sets all internal derivatives to zero, and evaluates:

```math
\begin{aligned}
  V_C &\leftarrow V_{\mathrm{comp}} \\
  V_{\mathrm{src}}^r
    &\leftarrow K_P V_r
       - X_L K_P^i I_r
       - (K_I + X_L K_P^r)I_i \\
  V_{\mathrm{src}}^i
    &\leftarrow K_P V_i
       + (K_I + X_L K_P^r)I_r
       - X_L K_P^i I_i \\
  V_E &\leftarrow
    \sqrt{
      (V_{\mathrm{src}}^r)^2
      + (V_{\mathrm{src}}^i)^2
    } \\
  I_N &\leftarrow \dfrac{K_C I_{\mathrm{fd}}}{V_E} \\
  F_{\mathrm{ex}} &\leftarrow f(I_N) \\
  V_B &\leftarrow \text{min}(V_E F_{\mathrm{ex}}, V_B^{\max}) \\
  V_{\mathrm{lv}} &\leftarrow \dfrac{E_{\mathrm{fd}}}{V_B} \\
  V_M &\leftarrow V_{\mathrm{lv}} \\
  V_G &\leftarrow \text{min}(K_G E_{\mathrm{fd}}, V_G^{\max}) \\
  e_M &\leftarrow 0 \\
  V_A &\leftarrow V_G \\
  V_R &\leftarrow V_A \\
  x_R &\leftarrow V_R \\
  e_V &\leftarrow 0 \\
  V_{\mathrm{ref}} &\leftarrow V_C - V_{\mathrm{uel}} - V_S
\end{aligned}
```

This closed-form start requires $V_E\ne 0$, $V_B\ne 0$, inactive
$V_R$, $V_M$, $V_G$, and $V_B$ limits, and the low-value gate selecting $V_M$.
Starts with active low-value gate limiting or saturated PI states are outside
these closed-form equations.

## Monitors

Monitor         | Units  | Description                         | Note
----------------|--------|-------------------------------------|------
`efd`           | [p.u.] | Field-voltage output                | $E_{\mathrm{fd}}$
`vm`            | [p.u.] | Inner regulator output              | $V_M$
`vc`            | [p.u.] | Sensed compensated voltage          | $V_C$
`va`            | [p.u.] | Lagged outer-regulator output       | $V_A$
`vr`            | [p.u.] | Outer regulator output              | $V_R$
`vg`            | [p.u.] | Exciter-output feedback signal      | $V_G$
`vlv`           | [p.u.] | Low-value gate output               | $V_{\mathrm{lv}}$
`ve`            | [p.u.] | Potential-source voltage magnitude  | $V_E$
`vb`            | [p.u.] | Rectifier source multiplier         | $V_B$
`in`            | [p.u.] | Normalized exciter loading current  | $I_N$
`fex`           | [p.u.] | Rectifier loading factor            | $F_{\mathrm{ex}}$
