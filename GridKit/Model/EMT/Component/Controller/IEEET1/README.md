# IEEE Type 1 Excitation System (IEEET1)

The EMT controller uses the same four differential and five algebraic equations
as the [PhasorDynamics IEEET1](../../../../PhasorDynamics/Exciter/IEEET1/README.md).
It senses three-phase instantaneous terminal voltage and supplies the machine
field voltage on the exciter per-unit base.

## Model Parameters

Parameter | Units | Description | Default
----------|-------|-------------|--------
`V` | V | Rated line-to-line RMS terminal voltage | Required, positive
`Tr` | s | Voltage-sensing time constant | 0
`Ka` | p.u. | Regulator gain | 50
`Ta` | s | Regulator time constant | 0.04
`Ke` | p.u. | Exciter coefficient; zero requests automatic initialization | -0.06
`Te` | s | Exciter time constant | 0.6
`Kf` | p.u. | Feedback gain | 0.09
`Tf` | s | Feedback time constant | 1.46
`Vrmin` | p.u. | Minimum regulator output | -1
`Vrmax` | p.u. | Maximum regulator output | 1
`E1`, `E2` | p.u. | Saturation voltages | 2.8, 3.73
`Se1`, `Se2` | p.u. | Saturation factors at `E1`, `E2` | 0.04, 0.33
`Ispdlim` | binary | Enable the field-voltage speed multiplier | 0

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

### Model Derived Parameters

The saturation contribution is $k_\mathrm{sat}=S_B q(E_{fd}'-S_A)$, where $q$ is
the CommonMath [quadratic ramp](../../../../../CommonMath.md#quadratic-ramp).
For two positive saturation factors,

```math
C=\sqrt{\frac{E_2 S_2}{E_1 S_1}},\qquad
S_A=\frac{C E_1-E_2}{C-1},\qquad
S_B=\frac{E_1 S_1}{(E_1-S_A)^2}.
```

If one factor is zero, its voltage is the knee $S_A$, and
$S_B=E_j S_j/(E_j-S_A)^2$ uses the other point. Disabled saturation sets
$S_A=S_B=0$.

A nonzero configured $K_E$ gives $K_E^\mathrm{eff}=K_E$. When `Ke = 0`,
initialization resolves
$K_E^\mathrm{eff}=(V_R^\mathrm{max}/10-k_\mathrm{sat})/E_{fd}'$ and requires
nonzero initial $E_{fd}'$. The configured coefficient remains unchanged.

## Model Ports

Name | Kind | Initialization | Description
-----|------|----------------|------------
`bus` | Three-phase input | Known | Instantaneous terminal voltages in volts
`speed` | Scalar input | Known | Optional machine rotor speed; 1 p.u. when unattached
`vref` | Scalar input | Solved | Optional voltage reference; initialized by the controller
`vs` | Scalar input | Known | Optional stabilizer input; zero when unattached
`vuel` | Scalar input | Known | Optional under-excitation limiter input; zero when unattached
`voel` | Scalar input | Known | Optional over-excitation limiter input; zero when unattached
`efd` | Scalar output | Known | Field voltage seeded by the machine

The EMT speed signal is rotor speed $\omega_r$, with synchronous speed equal to
one. The corresponding PhasorDynamics speed deviation is $\omega_r-1$.

## Model Variables

### Differential Variables

Symbol | Description
-------|------------
$V_{ts}$ | Sensed terminal voltage
$V_R$ | Regulator output
$E_{fd}'$ | Field voltage before the speed multiplier
$V_{fx}$ | Feedback state

### Algebraic Variables

Symbol | Description
-------|------------
$V_{tr}$ | Terminal-voltage error
$V_f$ | Feedback voltage
$V_E$ | Excitation control voltage
$E_{fd}$ | Field-voltage output
$k_\mathrm{sat}$ | Saturation contribution

## Model Equations

The voltage measurement and regulator drive are

```math
E_C=\frac{\sqrt{v_a^2+v_b^2+v_c^2}}{V},\qquad
f_R=\frac{-V_R+K_A V_{tr}}{T_A}.
```

For balanced sinusoidal voltages, $E_C$ is the terminal line-to-line RMS voltage
in per unit, independent of electrical angle. Under imbalance or harmonics it
is the instantaneous aggregate phase magnitude; the sensing lag filters its
ripple. The model does not extract a positive-sequence phasor or compensate
terminal voltage for a current-dependent impedance drop. At complete voltage
collapse the measured magnitude is zero; the Jacobian uses the zero gradient
convention at the nondifferentiable origin.

### Differential Equations

```math
\begin{aligned}
0&=-\dot V_{ts}+(E_C-V_{ts})/T_R,\\
0&=-\dot V_R+\operatorname{antiwindup}(V_R,f_R;V_R^\min,V_R^\max),\\
0&=-\dot E_{fd}'+(V_R-V_E-K_E^\mathrm{eff}E_{fd}')/T_E,\\
0&=-\dot V_{fx}+V_f/T_F.
\end{aligned}
```

The regulator uses the CommonMath smooth
[antiwindup](../../../../../CommonMath.md#antiwindup) function.

### Algebraic Equations

```math
\begin{aligned}
0&=-V_{ts}+V_\mathrm{ref}+V_S+V_\mathrm{UEL}+V_\mathrm{OEL}-V_{tr}-V_f,\\
0&=-T_F(V_f+V_{fx})+K_F E_{fd}',\\
0&=-V_E+k_\mathrm{sat},\\
0&=-E_{fd}+[1+(\omega_r-1)I_\mathrm{spdlim}]E_{fd}',\\
0&=-k_\mathrm{sat}+S_B q(E_{fd}'-S_A).
\end{aligned}
```

## Initialization

The machine initializes first and seeds `efd`. The controller reads the
terminal voltages and attached scalar inputs, then sets

```math
\begin{aligned}
E_{fd}'&=\frac{E_{fd}}{1+(\omega_r-1)I_\mathrm{spdlim}},\\
k_\mathrm{sat}&=S_B q(E_{fd}'-S_A),\qquad V_E=k_\mathrm{sat},\\
V_R&=K_E^\mathrm{eff}E_{fd}'+V_E,\qquad V_{tr}=V_R/K_A,\\
V_{fx}&=(K_F/T_F)E_{fd}',\qquad V_{ts}=E_C,\qquad V_f=0,\\
V_\mathrm{ref}&=E_C+V_{tr}-V_S-V_\mathrm{UEL}-V_\mathrm{OEL}.
\end{aligned}
```

All internal derivatives initialize to zero. An attached `vref` receives the
resolved reference; otherwise it is held locally. Attached inputs are read
live during residual and Jacobian evaluations.

## Monitors

Monitor | Description
--------|------------
`efd` | Field-voltage output
`ksat` | Saturation contribution
`vts` | Sensed terminal voltage
`vr` | Regulator output
`vref` | Active voltage reference
