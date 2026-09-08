# IEEET1 Model

The EMT controller uses the same four differential and five algebraic equations
as the [PhasorDynamics IEEET1](../../../../PhasorDynamics/Exciter/IEEET1/README.md).
It senses three-phase instantaneous terminal voltage and supplies the machine
field voltage on the exciter per-unit base.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$V$ | [V] | `V` | Rated line-to-line RMS terminal voltage | Required, positive
$T_R$ | [s] | `Tr` | Voltage-sensing time constant | Default $0$
$K_A$ | [p.u.] | `Ka` | Regulator gain | Default $50$
$T_A$ | [s] | `Ta` | Regulator time constant | Default $0.04$
$K_E$ | [p.u.] | `Ke` | Exciter coefficient | Default $-0.06$; zero requests initialization
$T_E$ | [s] | `Te` | Exciter time constant | Default $0.6$
$K_F$ | [p.u.] | `Kf` | Feedback gain | Default $0.09$
$T_F$ | [s] | `Tf` | Feedback time constant | Default $1.46$
$V_R^{\min}$ | [p.u.] | `Vrmin` | Minimum regulator output | Default $-1$
$V_R^{\max}$ | [p.u.] | `Vrmax` | Maximum regulator output | Default $1$
$E_1,E_2$ | [p.u.] | `E1`, `E2` | Saturation voltages | Defaults $2.8$, $3.73$
$S_1,S_2$ | [p.u.] | `Se1`, `Se2` | Saturation factors | Defaults $0.04$, $0.33$
$I_{\mathrm{spdlim}}$ | [-] | `Ispdlim` | Enable field-voltage speed multiplier | Default $0$

### Parameter Validation

All parameters and derived coefficients must be finite. `V` and `Ka` must be
positive, `Vrmin <= Vrmax`, and `Ispdlim` must be zero or one. Finite `Tr`, `Ta`,
`Te`, and `Tf` below $10^{-3}$ s are raised to that floor with a warning.

Saturation is disabled when `Se1 = Se2 = 0`. Otherwise, the saturation voltages
must be positive, the factors nonnegative, and both pairs strictly ordered in
the same direction. Invalid configurations, unattached terminal voltages,
unlinked attached inputs, and an unassigned `efd` output are rejected.
Initialization also rejects a nonpositive speed multiplier, nonfinite initial
states, or an initial regulator output outside the configured limits.

### Derived Parameters

The saturation contribution is $k_\mathrm{sat}=S_B q(E_{fd}'-S_A)$, where $q$ is
the CommonMath [quadratic ramp](../../../../../CommonMath.md#quadratic-ramp).
For two positive saturation factors,

```math
C=\sqrt{\dfrac{E_2 S_2}{E_1 S_1}},\qquad
S_A=\dfrac{C E_1-E_2}{C-1},\qquad
S_B=\dfrac{E_1 S_1}{(E_1-S_A)^2}.
```

If one factor is zero, its voltage is the knee $S_A$, and
$S_B=E_j S_j/(E_j-S_A)^2$ uses the other point. Disabled saturation sets
$S_A=S_B=0$.

A nonzero configured $K_E$ gives $K_E^\mathrm{eff}=K_E$. When `Ke = 0`,
initialization resolves
$K_E^\mathrm{eff}=(V_R^{\max}/10-k_\mathrm{sat})/E_{fd}'$ and requires
nonzero initial $E_{fd}'$. The configured coefficient remains unchanged.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `bus` | Input | [V] | Three-phase terminal voltage | Required
$\omega_r$ | `speed` | Input | [p.u.] | Machine rotor speed | Optional, defaults to one
$V_{\mathrm{ref}}$ | `vref` | Input | [p.u.] | Voltage reference | Inferred when unattached
$V_S$ | `vs` | Input | [p.u.] | Stabilizer input | Optional, defaults to zero
$V_{\mathrm{UEL}}$ | `vuel` | Input | [p.u.] | Under-excitation limiter input | Optional, defaults to zero
$V_{\mathrm{OEL}}$ | `voel` | Input | [p.u.] | Over-excitation limiter input | Optional, defaults to zero
$E_{fd}$ | `efd` | Output | [p.u.] | Field voltage | Seeded by the machine

The EMT speed signal is rotor speed $\omega_r$, with synchronous speed equal to
one. The corresponding PhasorDynamics speed deviation is $\omega_r-1$.

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$V_{ts}$ | [p.u.] | Sensed terminal voltage |
$V_R$ | [p.u.] | Regulator output |
$E_{fd}'$ | [p.u.] | Field voltage before speed multiplier |
$V_{fx}$ | [p.u.] | Feedback state |

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$V_{tr}$ | [p.u.] | Terminal-voltage error |
$V_f$ | [p.u.] | Feedback voltage |
$V_E$ | [p.u.] | Excitation control voltage |
$E_{fd}$ | [p.u.] | Field-voltage output |
$k_{\mathrm{sat}}$ | [p.u.] | Saturation contribution |

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Terminal voltage | $\mathbf{v}\in\mathbb{R}^3$
$\omega_r$ | [p.u.] | Machine rotor speed | Defaults to one
$V_{\mathrm{ref}}$ | [p.u.] | Voltage reference | Inferred when unattached
$V_S$ | [p.u.] | Stabilizer input | Defaults to zero
$V_{\mathrm{UEL}}$ | [p.u.] | Under-excitation limiter input | Defaults to zero
$V_{\mathrm{OEL}}$ | [p.u.] | Over-excitation limiter input | Defaults to zero

## Model Equations

The voltage measurement and regulator drive are

```math
E_C=\dfrac{\sqrt{v_a^2+v_b^2+v_c^2}}{V},\qquad
f_R=\dfrac{-V_R+K_A V_{tr}}{T_A}.
```

For balanced sinusoidal voltages, $E_C$ is the terminal line-to-line RMS voltage
in per unit, independent of electrical angle. Under imbalance or harmonics it
is the instantaneous aggregate phase magnitude; the sensing lag filters its
ripple. The model does not extract a positive-sequence phasor or compensate
terminal voltage for a current-dependent impedance drop. At complete voltage
collapse the measured magnitude is zero; the Jacobian uses the zero gradient
convention at the nondifferentiable origin.

### Internal Equations

#### Differential

```math
\begin{aligned}
0 &= -\dfrac{\mathrm{d}V_{ts}}{\mathrm{d}t}+(E_C-V_{ts})/T_R \\
0 &= -\dfrac{\mathrm{d}V_R}{\mathrm{d}t}+\mathrm{antiwindup}(V_R,f_R;V_R^{\min},V_R^{\max}) \\
0 &= -\dfrac{\mathrm{d}E_{fd}'}{\mathrm{d}t}+(V_R-V_E-K_E^\mathrm{eff}E_{fd}')/T_E \\
0 &= -\dfrac{\mathrm{d}V_{fx}}{\mathrm{d}t}+V_f/T_F.
\end{aligned}
```

The regulator uses the CommonMath smooth
[antiwindup](../../../../../CommonMath.md#antiwindup) function.

#### Algebraic

```math
\begin{aligned}
0 &= -V_{ts}+V_\mathrm{ref}+V_S+V_\mathrm{UEL}+V_\mathrm{OEL}-V_{tr}-V_f \\
0 &= -T_F(V_f+V_{fx})+K_F E_{fd}' \\
0 &= -V_E+k_\mathrm{sat} \\
0 &= -E_{fd}+[1+(\omega_r-1)I_\mathrm{spdlim}]E_{fd}' \\
0 &= -k_\mathrm{sat}+S_B q(E_{fd}'-S_A).
\end{aligned}
```

### External Equations

None.

## Initialization

The machine initializes first and seeds `efd`. The controller reads the
terminal voltages and attached scalar inputs, then sets

```math
\begin{aligned}
E_{fd}' &\leftarrow \dfrac{E_{fd}}{1+(\omega_r-1)I_\mathrm{spdlim}} \\
k_\mathrm{sat} &\leftarrow S_B q(E_{fd}'-S_A),\qquad V_E\leftarrow k_\mathrm{sat} \\
V_R &\leftarrow K_E^\mathrm{eff}E_{fd}'+V_E,\qquad V_{tr}\leftarrow V_R/K_A \\
V_{fx} &\leftarrow (K_F/T_F)E_{fd}',\qquad V_{ts}\leftarrow E_C,\qquad V_f\leftarrow 0 \\
V_\mathrm{ref} &\leftarrow E_C+V_{tr}-V_S-V_\mathrm{UEL}-V_\mathrm{OEL}.
\end{aligned}
```

All internal derivatives initially default to zero. The inferred reference
above is held locally for an unattached `vref`. An attached reference is
preserved and read live during residual and Jacobian evaluations; consistent
initialization resolves the derivatives if it differs from the inferred value.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`efd` | [p.u.] | Field-voltage output | $E_{fd}$
`ksat` | [p.u.] | Saturation contribution | $k_{\mathrm{sat}}$
`vts` | [p.u.] | Sensed terminal voltage | $V_{ts}$
`vr` | [p.u.] | Regulator output | $V_R$
`vref` | [p.u.] | Active voltage reference | $V_{\mathrm{ref}}$
