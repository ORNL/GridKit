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
$L$ | [H] | `L` | Inverter-side filter inductance | Required
$K_P$ | [$\Omega$] | `Kp` | Proportional gain | Required
$K_I$ | [$\Omega/\mathrm{s}$] | `Ki` | Integral gain | Required
$K_{\mathrm{aw}}$ | [$\mathrm{s}^{-1}$] | `Kaw` | Tracking anti-windup gain | Required
$I^{\max}$ | [A] | `Imax` | Current-reference norm limit | Required
$M^{\max}$ | [-] | `Mmax` | Sinusoidal modulation limit | Required

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
\mathcal{J} &=
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
$\mathbf{i}^{\mathrm{ref}}$ | `iref` | Input | [A] | Total current reference | $\mathbf{i}^{\mathrm{ref}} \in \mathbb{R}^2$
$\omega$ | `omega` | Input | [rad/s] | Electrical angular frequency of the $dq$ frame | Supplied by the angle source
$v_{\mathrm{dc}}$ | `vdc` | Input | [V] | DC-link voltage | $v_{\mathrm{dc}} \ge 0$
$\mathbf{i}^{\mathrm{lim}}$ | `ilim` | Output | [A] | Limited current reference | $\mathbf{i}^{\mathrm{lim}} \in \mathbb{R}^2$
$\mathbf{u}$ | `u` | Output | [V] | Converter voltage command | $\mathbf{u} \in \mathbb{R}^2$

All vectors use $(d,q)$ order in the same power-invariant
[Park](../../../Operators/Reference/Park/README.md) frame, with zero-sequence
components omitted. All inputs must be connected and finite.

The `omega` input is read from the angle source throughout the simulation.
It may vary with time or be supplied as a constant. The `ilim` output provides
outer-loop anti-windup feedback.

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

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}$ | [V] | Filter-capacitor voltage | $\mathbf{v} \in \mathbb{R}^2$
$\mathbf{i}$ | [A] | Inverter-side filter current | $\mathbf{i} \in \mathbb{R}^2$
$\mathbf{i}^{\mathrm{ref}}$ | [A] | Total current reference | $\mathbf{i}^{\mathrm{ref}} \in \mathbb{R}^2$
$\omega$ | [rad/s] | Electrical angular frequency of the $dq$ frame | Supplied by the angle source
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
\mathbf{b} &= \mathbf{v} + \mathcal{J}\omega L\mathbf{i} \\
\mathbf{z} &= \mathbf{b} + K_P\mathbf{e} + \boldsymbol{\xi} \\
\mathcal{L}_i(\mathbf{i}^{\mathrm{ref}}) &=
  \max\left(1,a_i\|\mathbf{i}^{\mathrm{ref}}\|_2^2\right) \\
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

None.

### External Equations

```math
\begin{aligned}
\mathbf{i}^{\mathrm{lim}} &\leftarrow
  \dfrac{\mathbf{i}^{\mathrm{ref}}}{\sqrt{\mathcal{L}_i(\mathbf{i}^{\mathrm{ref}})}} \\
\mathbf{u} &\leftarrow
  \dfrac{v_{\mathrm{dc}}\mathbf{z}}{\sqrt{\mathcal{L}_u(v_{\mathrm{dc}},\mathbf{z})}}
\end{aligned}
```

The outputs are algebraic expressions without owned DAE variables. The
current limit bounds the reference; instantaneous filter current can overshoot.
The anti-windup correction uses the limited voltage command, not the switched
bridge voltage.

## Initialization

Initialize $\boldsymbol{\xi}$ from the finite state-file values `xid` and `xiq`
(default zero). The consistent-initial-condition solve preserves these states
and obtains their derivatives from the connected inputs. In unsaturated balanced
steady state with zero current error, $\boldsymbol{\xi}=R\mathbf{i}$.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`xi` | [V] | Integral contribution | $\boldsymbol{\xi} \in \mathbb{R}^2$
`ilim` | [A] | Limited current reference | $\mathbf{i}^{\mathrm{lim}} \in \mathbb{R}^2$
`u` | [V] | Converter voltage command | $\mathbf{u} \in \mathbb{R}^2$

In case JSON, vector ports use two signal IDs in $(d,q)$ order. The vector
monitors `xi`, `ilim`, and `u` expand to `xid`, `xiq`, `ilimd`, `ilimq`, `ud`,
and `uq`, respectively.

[^unifi]: UNIFI Consortium, [*UNIFI's Grid-Forming (GFM) Inverter Reference Design: A Tutorial on Modeling, Control, and Experimental Implementation*](https://docs.nlr.gov/docs/fy25osti/92994.pdf),
    NREL/TP-5D00-92994, July 2025, Section 2.2, equations (7), (12), and Section 5.1.
    The DC-dependent voltage limit and its tracking anti-windup extend the
    reference current PI realization.

[^tracking-anti-windup]: K. J. Åström and L. Rundqwist,
    [*Integrator Windup and How to Avoid It*](https://doi.org/10.23919/ACC.1989.4790464),
    American Control Conference, pp. 1693–1698, 1989, Fig. 2.
    Back-calculation applied componentwise, with $K_{\mathrm{aw}} = 1/T_t$.
