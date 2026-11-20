# LVRT Model

`LVRT` is a low-voltage ride-through. It limits an active-current order by a
ceiling that falls with the terminal voltage and recovers at a limited rate.
All quantities are per unit on the device base.

## Block Diagram

None.

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$T_M$ | [s] | `TM` | Time constant of the current ceiling | Required, positive
$V_{L0}$ | [p.u.] | `VL0` | Voltage at which the ceiling target is zero | Nonnegative
$V_{L1}$ | [p.u.] | `VL1` | Voltage at which the ceiling target is full | Above $V_{L0}$
$I_{L1}$ | [p.u.] | `IL1` | Full ceiling target | Required, positive
$R_p^{\max}$ | [p.u./s] | `Rpmax` | Maximum recovery rate of the ceiling | Required, positive
$V_{A0}$ | [p.u.] | `VA0` | Voltage at which no current is available | Nonnegative
$V_{A1}$ | [p.u.] | `VA1` | Voltage at which the whole ceiling is available | Above $V_{A0}$

### Parameter Validation

```math
\begin{aligned}
T_M &> 0 \\
0 \le V_{L0} &< V_{L1} \\
0 \le V_{A0} &< V_{A1} \\
I_{L1}, R_p^{\max} &> 0
\end{aligned}
```

### Model Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$v_t$ | `vt` | Input | [p.u.] | Terminal voltage | Scalar
$i_p^\mathrm{cmd}$ | `ipcmd` | Input | [p.u.] | Active-current order | Scalar
$i_p$ | `ip` | Output | [p.u.] | Limited active-current order | Internal variable

The limited order can drive input `ip` of
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

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$i_l$ | [p.u.] | Recovering active-current ceiling |

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$i_p$ | [p.u.] | Limited active-current order |

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$v_t$ | [p.u.] | Terminal voltage | Held signal unless the case connects the input
$i_p^\mathrm{cmd}$ | [p.u.] | Active-current order | Held signal unless the case connects the input

## Model Equations

Smooth functions: [linseg](../../../../../CommonMath.md#linear-segment) $\ell$.

The available current is

```math
i_a = i_l\,\ell(v_t, V_{A0}, V_{A1}, 1)
```

### Internal Equations

#### Differential

```math
0 = -\dfrac{\mathrm{d}i_l}{\mathrm{d}t}
    + \min\left(\dfrac{1}{T_M}\left(\ell(v_t, V_{L0}, V_{L1}, I_{L1}) - i_l\right), R_p^{\max}\right)
```

#### Algebraic

```math
0 = -i_p + i_p^\mathrm{cmd} - \max(i_p^\mathrm{cmd} - i_a, 0)
```

### External Equations

None.

## Initialization

None.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`il` | [p.u.] | Recovering active-current ceiling |
`ip` | [p.u.] | Limited active-current order |
