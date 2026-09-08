# InnerCurrentControl Model

`InnerCurrentControl` regulates converter current in power-invariant $dq$
coordinates using PI control, capacitor-voltage feedforward, and cross-coupling
compensation.[^unifi] Direction-preserving current and voltage limits use
tracking anti-windup.[^tracking-anti-windup]

## Block Diagram

![InnerCurrentControl model block diagram](../../../../../../docs/Figures/EMT/Controller/InnerCurrentControl/diagram.png)

Figure 1: InnerCurrentControl model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$L$ | [H] | `L` | Inverter-side filter inductance | Required, positive
$K_P$ | [$\Omega$] | `Kp` | Proportional gain | Required, positive
$K_I$ | [$\Omega/\mathrm{s}$] | `Ki` | Integral gain | Required, positive
$K_{\mathrm{aw}}$ | [$\mathrm{s}^{-1}$] | `Kaw` | Tracking anti-windup gain | Required, positive
$I^{\max}$ | [A] | `Imax` | Current-command norm limit | Required, positive
$M^{\max}$ | [-] | `Mmax` | Sinusoidal modulation limit | Required, $0 < M^{\max} \le 1$

For a balanced fundamental, $I^{\max} = \sqrt{3}\,I_{\mathrm{phase,rms}}^{\max}$.

### Parameter Validation

All parameters must be finite.

```math
\begin{aligned}
L, K_P, K_I, K_{\mathrm{aw}}, I^{\max} &> 0 \\
0 &< M^{\max} \le 1
\end{aligned}
```

### Derived Parameters

The limit coefficients and matrix representation of the complex unit are

```math
\begin{aligned}
a_i &= \dfrac{1}{(I^{\max})^2} \\
a_u &= \dfrac{8}{3(M^{\max})^2} \\
\mathbf{J} &=
\begin{bmatrix}
0 & -1 \\
1 & 0
\end{bmatrix}
\end{aligned}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `v` | Input | [V] | Filter-capacitor voltage | $\mathbf{v} \in \mathbb{R}^2$
$\mathbf{i}$ | `i` | Input | [A] | Inverter-side filter current | $\mathbf{i} \in \mathbb{R}^2$
$\mathbf{i}^{\mathrm{cmd}}$ | `icmd` | Input | [A] | Total current command | $\mathbf{i}^{\mathrm{cmd}} \in \mathbb{R}^2$
$\omega$ | `omega` | Input | [rad/s] | Electrical angular frequency of the $dq$ frame | Supplied through the frequency signal port
$v_{\mathrm{dc}}$ | `vdc` | Input | [V] | DC-link voltage | $v_{\mathrm{dc}} \ge 0$
$\mathbf{i}^{\mathrm{lim}}$ | `ilim` | Output | [A] | Limited current command | $\mathbf{i}^{\mathrm{lim}} \in \mathbb{R}^2$
$\mathbf{u}$ | `u` | Output | [V] | Converter voltage command | $\mathbf{u} \in \mathbb{R}^2$

All vectors use $(d,q)$ order in the same power-invariant
[Park](../../../Operators/Reference/Park/README.md) frame, with zero-sequence
components omitted. All inputs must be connected and finite.

The `omega` signal input is read throughout the simulation. In the switching
examples, [PLL](../../../Operators/Reference/PLL/README.md) supplies this frequency
and the common Park angle through its output signals. The `ilim` output
provides outer-loop anti-windup feedback.

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\boldsymbol{\xi}$ | [V] | Integral contribution | $\boldsymbol{\xi} \in \mathbb{R}^2$

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{lim}}$ | [A] | Limited current command | $\mathbf{i}^{\mathrm{lim}} \in \mathbb{R}^2$
$\mathbf{u}$ | [V] | Converter voltage command | $\mathbf{u} \in \mathbb{R}^2$

### External Variables

#### Differential

Connected voltage, current, and angle-source variables may be differential.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Filter-capacitor voltage | $\mathbf{v} \in \mathbb{R}^2$
$\mathbf{i}$ | [A] | Inverter-side filter current | $\mathbf{i} \in \mathbb{R}^2$
$\mathbf{i}^{\mathrm{cmd}}$ | [A] | Total current command | $\mathbf{i}^{\mathrm{cmd}} \in \mathbb{R}^2$
$\omega$ | [rad/s] | Electrical angular frequency of the $dq$ frame | Supplied through the frequency signal port
$v_{\mathrm{dc}}$ | [V] | DC-link voltage | $v_{\mathrm{dc}} \ge 0$

## Model Equations

The supplied frequency and the common Park angle $\theta$ satisfy

```math
\omega = \dfrac{\mathrm{d}\theta}{\mathrm{d}t}
```

The current error, feedforward voltage, unlimited voltage command, and limiter factors are

```math
\begin{aligned}
\mathbf{e} &= \mathbf{i}^{\mathrm{lim}} - \mathbf{i} \\
\mathbf{b} &= \mathbf{v} + \mathbf{J}\omega L\mathbf{i} \\
\mathbf{z} &= \mathbf{b} + K_P\mathbf{e} + \boldsymbol{\xi} \\
\mathcal{L}_i(\mathbf{i}^{\mathrm{cmd}}) &=
  \max\left(1,a_i\|\mathbf{i}^{\mathrm{cmd}}\|_2^2\right) \\
\mathcal{L}_u(v_{\mathrm{dc}},\mathbf{z}) &=
  \max\left(v_{\mathrm{dc}}^2,a_u\|\mathbf{z}\|_2^2\right)
\end{aligned}
```

The voltage limiter assumes sinusoidal PWM without zero-sequence injection.
The filter resistance $R$ remains in the physical circuit. Nominal RL-based
tuning is $K_P=L\omega_{ci}$ and $K_I=R\omega_{ci}$, where $\omega_{ci}$ is the
current-loop bandwidth in rad/s; PWM delay and filter dynamics constrain the
usable bandwidth.

Both direction-preserving limits use the CommonMath smooth
[`max`](../../../../../CommonMath.md#maximum).
The voltage limiter smooths squared voltages in $\mathrm{V}^2$. At zero DC
voltage, the voltage command is zero.

### Internal Equations

#### Differential

```math
0 = -\dfrac{\mathrm{d}\boldsymbol{\xi}}{\mathrm{d}t}
    + K_I\mathbf{e} + K_{\mathrm{aw}}(\mathbf{u}-\mathbf{z})
```

#### Algebraic

```math
\begin{aligned}
0 &= \mathbf{i}^{\mathrm{lim}}-
  \dfrac{\mathbf{i}^{\mathrm{cmd}}}{\sqrt{\mathcal{L}_i(\mathbf{i}^{\mathrm{cmd}})}} \\
0 &= \mathbf{u}-
  \dfrac{v_{\mathrm{dc}}\mathbf{z}}{\sqrt{\mathcal{L}_u(v_{\mathrm{dc}},\mathbf{z})}}
\end{aligned}
```

### External Equations

None.

Both output vectors are owned algebraic variables. The current limit bounds
the command; instantaneous filter current can overshoot. Anti-windup tracks
the limited voltage command, not the switched bridge voltage.

## Initialization

The initialized current command defines $\mathbf{i}^{\mathrm{lim}}$ and
$\mathbf{e}$. The default voltage output is the limited value of
$\mathbf{b}+K_P\mathbf{e}$, giving zero integral contribution.

The state-file keys `ud` and `uq` replace the respective defaults with finite
output values. For positive DC voltage, prescribed voltage outputs must lie
strictly inside the modulation limit. The radial smooth limiter is inverted
to recover $\mathbf{z}$, then

```math
\boldsymbol{\xi} \leftarrow \mathbf{z}-\mathbf{b}-K_P\mathbf{e}.
```

A clipped output does not uniquely determine the integral state. At zero DC
voltage, only zero voltage outputs are admissible and the integral defaults
to zero. Optional `ilimd` and `ilimq` must agree with the current limiter;
they cannot override the connected command. The integral states `xid` and
`xiq` cannot be prescribed in the state file.

Derivatives start at zero; the consistent-initial-condition solve preserves
the integral states and obtains derivatives and algebraic outputs from the
connected inputs. In unsaturated balanced steady state with zero current
error, $\boldsymbol{\xi}=R\mathbf{i}$.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`xi` | [V] | Integral contribution | $\boldsymbol{\xi} \in \mathbb{R}^2$
`ilim` | [A] | Limited current command | $\mathbf{i}^{\mathrm{lim}} \in \mathbb{R}^2$
`u` | [V] | Converter voltage command | $\mathbf{u} \in \mathbb{R}^2$

See [case connections](../../../INPUT_FORMAT.md#case-connections) for vector
ports and monitor expansion.

[^unifi]: UNIFI Consortium, [*UNIFI's Grid-Forming (GFM) Inverter Reference Design: A Tutorial on Modeling, Control, and Experimental Implementation*](https://docs.nlr.gov/docs/fy25osti/92994.pdf),
    NREL/TP-5D00-92994, July 2025, Section 2.2, equations (7), (12), and Section 5.1.
    The DC-dependent voltage limit and its tracking anti-windup extend the
    reference current PI realization.

[^tracking-anti-windup]: K. J. Åström and L. Rundqwist,
    [*Integrator Windup and How to Avoid It*](https://doi.org/10.23919/ACC.1989.4790464),
    American Control Conference, pp. 1693–1698, 1989, Fig. 2.
    Back-calculation applied componentwise, with $K_{\mathrm{aw}} = 1/T_t$.
