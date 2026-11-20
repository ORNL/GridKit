# HVRCM Model

`HVRCM` is a high-voltage reactive current management. It reduces a
reactive-current order as the terminal voltage rises above a threshold. All
quantities are per unit on the device base.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$V_\mathrm{hv}^{\max}$ | [p.u.] | `Vhvmax` | Voltage above which reactive current is absorbed | Required, positive
$K_\mathrm{hv}$ | [p.u.] | `Khv` | Gain of the voltage response | Nonnegative
$Q^{\min}$ | [p.u.] | `Qmin` | Lower bound of the reduced order | Nonpositive

### Parameter Validation

```math
\begin{aligned}
V_\mathrm{hv}^{\max} &> 0 \\
K_\mathrm{hv} &\ge 0 \\
Q^{\min} &\le 0
\end{aligned}
```

### Model Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$v_t$ | `vt` | Input | [p.u.] | Terminal voltage | Scalar
$i_q^\mathrm{cmd}$ | `iqcmd` | Input | [p.u.] | Reactive-current order | Scalar
$i_q$ | `iq` | Output | [p.u.] | Reduced reactive-current order | Internal variable

The reduced order can drive input `iq` of
[CurrentControl](../CurrentControl/README.md).

## Submodels

None.

### Submodel Validation

None.

### Submodel Wiring

None.

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$i_q$ | [p.u.] | Reduced reactive-current order |

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$v_t$ | [p.u.] | Terminal voltage | Held signal unless the case connects the input
$i_q^\mathrm{cmd}$ | [p.u.] | Reactive-current order | Held signal unless the case connects the input

## Model Equations

Smooth functions: [ramp](../../../../../CommonMath.md#ramp) $\rho$.

The absorbed current and the floor of the order are

```math
\begin{aligned}
a &= K_\mathrm{hv}\,\rho(v_t - V_\mathrm{hv}^{\max}) \\
b &= \min(i_q^\mathrm{cmd}, Q^{\min})
\end{aligned}
```

An order already below $Q^{\min}$ passes through.

### Internal Equations

#### Differential

None.

#### Algebraic

```math
0 = -i_q + b + \max(i_q^\mathrm{cmd} - b - a, 0)
```

### External Equations

None.

## Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`iq` | [p.u.] | Reduced reactive-current order |
