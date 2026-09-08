# OuterPowerControl Model

`OuterPowerControl` implements an outer-loop PI controller in power-invariant
$dq$ coordinates. It compares derived current references with measured
currents and supplies a current command to
[InnerCurrentControl](../InnerCurrentControl/README.md).

## Block Diagram

![OuterPowerControl model block diagram](../../../../../../docs/Figures/EMT/Controller/OuterPowerControl/diagram.png)

Figure 1: OuterPowerControl model

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$V$ | [V] | `V` | Rated line-to-line RMS voltage | Required, positive
$P^{\mathrm{ref}}$ | [W] | `Pref` | Active-power setpoint | Default zero
$Q^{\mathrm{ref}}$ | [var] | `Qref` | Reactive-power setpoint | Default zero
$K_P$ | [-] | `Kp` | Proportional gain | Required, positive
$K_I$ | [$\mathrm{s}^{-1}$] | `Ki` | Integral gain | Required, positive
$K_{\mathrm{aw}}$ | [$\mathrm{s}^{-1}$] | `Kaw` | Tracking anti-windup gain | Required, positive

### Parameter Validation

All parameters must be finite. The voltage rating and gains must be positive;
power setpoints may be positive, zero, or negative.

### Derived Parameters

For the voltage-aligned, power-invariant $dq$ convention,

```math
\mathbf{i}^{\mathrm{ref}} =
\begin{bmatrix}P^{\mathrm{ref}}/V\\-Q^{\mathrm{ref}}/V\end{bmatrix}.
```

These setpoints use rated voltage and remain fixed during a run. They match
requested P/Q at $v_d=V$, $v_q=0$; they do not impose constant power when the
terminal voltage differs from its rating.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{i}$ | `i` | Input | [A] | Measured current | $\mathbf{i} \in \mathbb{R}^2$
$\mathbf{i}^{\mathrm{lim}}$ | `ilim` | Input | [A] | Limited current command | From InnerCurrentControl
$\mathbf{i}^{\mathrm{cmd}}$ | `icmd` | Output | [A] | Inverter-side current command | $\mathbf{i}^{\mathrm{cmd}} \in \mathbb{R}^2$

All vectors use $(d,q)$ order in the same power-invariant
[Park](../../../Operators/Reference/Park/README.md) frame. A forward Park
transform supplies the measured currents directly. The controller uses its
derived current references internally. All inputs must be connected and finite.

Connect `icmd` to InnerCurrentControl's `icmd` input and return its `ilim`
output. For an LCL filter, the outer loop measures grid-side current while
InnerCurrentControl measures inverter-side current.

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

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}^{\mathrm{cmd}}$ | [A] | Inverter-side current command | $\mathbf{i}^{\mathrm{cmd}} \in \mathbb{R}^2$

### External Variables

#### Differential

The connected measured-current variables may be differential.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\mathbf{i}$ | [A] | Measured current | When algebraic
$\mathbf{i}^{\mathrm{lim}}$ | [A] | Limited current command |

## Model Equations

The current error is

```math
\mathbf{e} = \mathbf{i}^{\mathrm{ref}}-\mathbf{i}.
```

### Internal Equations

#### Differential

```math
0 = -\dfrac{\mathrm{d}\boldsymbol{\eta}}{\mathrm{d}t}
    + K_I\mathbf{e}
    + K_{\mathrm{aw}}(\mathbf{i}^{\mathrm{lim}}-\mathbf{i}^{\mathrm{cmd}})
```

#### Algebraic

```math
0 = \mathbf{i}^{\mathrm{cmd}}-K_P\mathbf{e}-\boldsymbol{\eta}
```

### External Equations

None.

InnerCurrentControl owns the smooth circular current limiter. Tracking the
limited command prevents outer-loop windup. Both command components are owned
algebraic variables, so the feedback has no recursive computed-signal dependency.

## Initialization

The derived current references and initialized measurements define $\mathbf{e}$.
The default current-command outputs are

```math
\mathbf{i}^{\mathrm{cmd}} \leftarrow K_P\mathbf{e}.
```

The state-file keys `icmdd` and `icmdq` replace the respective defaults with
finite output values. The integral contribution is then derived from them:

```math
\boldsymbol{\eta} \leftarrow \mathbf{i}^{\mathrm{cmd}}-K_P\mathbf{e}.
```

Omitted outputs give zero integral contribution, up to roundoff. The integral
states `etad` and `etaq` cannot be prescribed in the state file. Derivatives
start at zero; the consistent-initial-condition solve preserves the integral
states and obtains derivatives and algebraic commands from the connected inputs.

In unsaturated steady state, $\mathbf{i}=\mathbf{i}^{\mathrm{ref}}$ and
$\boldsymbol{\eta}$ equals the inverter-side current command. For an L filter,
initialize the command from measured current. For an LCL filter in balanced
steady state, include capacitor current:

```math
\mathbf{i}^{\mathrm{cmd}} \leftarrow
\mathbf{i}+\mathbf{J}\omega C\mathbf{v},\qquad
\mathbf{J}=\begin{bmatrix}0&-1\\1&0\end{bmatrix}.
```

Here $C$ is filter capacitance, $\omega$ is frame angular frequency, and
$\mathbf{v}$ is capacitor voltage. They are operating-point quantities supplied
by the surrounding circuit, not controller parameters or inputs.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`eta` | [A] | Integral contribution | $\boldsymbol{\eta} \in \mathbb{R}^2$
`icmd` | [A] | Inverter-side current command | $\mathbf{i}^{\mathrm{cmd}} \in \mathbb{R}^2$

See [case connections](../../../INPUT_FORMAT.md#case-connections) for vector
ports and monitor expansion.
