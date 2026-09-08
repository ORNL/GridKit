# OuterVoltageControl Model

`OuterVoltageControl` regulates filter-capacitor voltage in power-invariant
$dq$ coordinates using PI control, grid-current feedforward, and cross-coupling
compensation.[^unifi] It supplies the current reference for
[InnerCurrentControl](../InnerCurrentControl/README.md).

## Block Diagram

![OuterVoltageControl model block diagram](../../../../../../docs/Figures/EMT/Controller/OuterVoltageControl/diagram.png)

Figure 1: OuterVoltageControl model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$C$ | [F] | `C` | Filter capacitance | Required
$K_P$ | [S] | `Kp` | Proportional gain | Required
$K_I$ | [S/s] | `Ki` | Integral gain | Required
$K_{\mathrm{aw}}$ | [$\mathrm{s}^{-1}$] | `Kaw` | Tracking anti-windup gain | Required

### Parameter Validation

All parameters must be finite and positive.

### Derived Parameters

The matrix representation of the complex unit is

```math
\mathcal{J} =
\begin{bmatrix}
0 & -1 \\
1 & 0
\end{bmatrix}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}^{\mathrm{ref}}$ | `vref` | Input | [V] | Capacitor-voltage reference | $\mathbf{v}^{\mathrm{ref}} \in \mathbb{R}^2$
$\mathbf{v}$ | `v` | Input | [V] | Filter-capacitor voltage | $\mathbf{v} \in \mathbb{R}^2$
$\mathbf{i}_g$ | `ig` | Input | [A] | Grid-side filter current | $\mathbf{i}_g \in \mathbb{R}^2$
$\omega$ | `omega` | Input | [rad/s] | Electrical angular frequency of the $dq$ frame | Supplied by the angle source
$\mathbf{i}^{\mathrm{lim}}$ | `ilim` | Input | [A] | Limited current reference | From InnerCurrentControl
$\mathbf{i}^{\mathrm{ref}}$ | `iref` | Output | [A] | Total current reference | $\mathbf{i}^{\mathrm{ref}} \in \mathbb{R}^2$

All vectors use $(d,q)$ order in the same power-invariant
[Park](../../../Operators/Reference/Park/README.md) frame, with zero-sequence
components omitted. All inputs must be connected and finite. A grid-forming
primary controller supplies the common frame angle, frequency, and
$\mathbf{v}^{\mathrm{ref}}=[V^\star,0]^\mathsf{T}$. For balanced voltage,
$V^\star=V_{\mathrm{LL,rms}}^\star$.

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\boldsymbol{\eta}$ | [A] | Integral contribution | $\boldsymbol{\eta} \in \mathbb{R}^2$

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{v}^{\mathrm{ref}}$ | [V] | Capacitor-voltage reference | $\mathbf{v}^{\mathrm{ref}} \in \mathbb{R}^2$
$\mathbf{v}$ | [V] | Filter-capacitor voltage | $\mathbf{v} \in \mathbb{R}^2$
$\mathbf{i}_g$ | [A] | Grid-side filter current | $\mathbf{i}_g \in \mathbb{R}^2$
$\omega$ | [rad/s] | Electrical angular frequency of the $dq$ frame | Supplied by the angle source
$\mathbf{i}^{\mathrm{lim}}$ | [A] | Limited current reference | From InnerCurrentControl

## Model Equations

The voltage error and feedforward current are

```math
\begin{aligned}
\mathbf{e} &= \mathbf{v}^{\mathrm{ref}}-\mathbf{v} \\
\mathbf{b} &= \mathbf{i}_g+\mathcal{J}\omega C\mathbf{v}
\end{aligned}
```

### Internal Equations

#### Differential

```math
0 = -\dfrac{\mathrm{d}\boldsymbol{\eta}}{\mathrm{d}t}
    + K_I\mathbf{e}
    + K_{\mathrm{aw}}(\mathbf{i}^{\mathrm{lim}}-\mathbf{i}^{\mathrm{ref}})
```

#### Algebraic

None.

### External Equations

```math
\mathbf{i}^{\mathrm{ref}} \leftarrow \mathbf{b}+K_P\mathbf{e}+\boldsymbol{\eta}
```

The output is an algebraic expression without owned DAE variables. The limited
reference feeds back to the integrator, preventing windup while the inner
current-reference limiter is active.

## Initialization

Initialize $\boldsymbol{\eta}$ from the finite state-file values `etad` and
`etaq` (default zero). The consistent-initial-condition solve preserves these
states and obtains their derivatives from the connected inputs. In unsaturated
balanced steady state with zero voltage error, $\boldsymbol{\eta}=0$.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`eta` | [A] | Integral contribution | $\boldsymbol{\eta} \in \mathbb{R}^2$
`iref` | [A] | Total current reference | $\mathbf{i}^{\mathrm{ref}} \in \mathbb{R}^2$

In case JSON, vector ports use two signal IDs in $(d,q)$ order. The vector
monitors `eta` and `iref` expand to `etad`, `etaq`, `irefd`, and `irefq`,
respectively.

[^unifi]: UNIFI Consortium, [*UNIFI's Grid-Forming (GFM) Inverter Reference Design: A Tutorial on Modeling, Control, and Experimental Implementation*](https://docs.nlr.gov/docs/fy25osti/92994.pdf),
    NREL/TP-5D00-92994, July 2025, Section 2.2, equation (6), and Figure 2.
