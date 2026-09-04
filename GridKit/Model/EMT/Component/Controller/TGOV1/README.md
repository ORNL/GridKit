# TGOV1 Model

`Tgov1` represents the TGOV1 steam turbine governor carried into the EMT
subsystem from
[PhasorDynamics TGOV1](../../../../PhasorDynamics/Governor/Tgov1/README.md).
The equations are identical; the speed input is the machine rotor speed
$\omega_r$, the droop acts on the deviation $\Delta\omega = \omega_r - 1$,
and the mechanical power is exchanged in the machine per-unit base with no
system-base conversion.

> [!NOTE]
> Set $T_\mathrm{rate}$ equal to the connected machine MVA base; the
> mechanical-power signal is in that base.

## Block Diagram

![](../../../../../../docs/Figures/TGOV1.JPG)

Figure 1: Governor TGOV1 model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$T_\mathrm{rate}$ | [MVA] | `Trate` | Governor component power base | Required, positive
$R$ | [p.u.] | `R` | Permanent droop | Nonzero
$T_1$ | [sec] | `T1` | Steam-bowl time constant | Raised to the time-constant floor
$T_2$ | [sec] | `T2` | Turbine numerator time constant |
$T_3$ | [sec] | `T3` | Reheater time constant | Raised to the time-constant floor
$P_v^{\max}$ | [p.u.] | `Pvmax` | Maximum valve position |
$P_v^{\min}$ | [p.u.] | `Pvmin` | Minimum valve position |
$D_t$ | [p.u.] | `Dt` | Turbine damping coefficient |

### Parameter Validation

```math
\begin{aligned}
T_\mathrm{rate} &> 0 \\
R &\ne 0 \\
P_v^{\min} &\le P_v^{\max}
\end{aligned}
```

### Derived Parameters

Let $\epsilon_T = 10^{-3}\ \mathrm{s}$. A time constant below $\epsilon_T$ is
raised to that floor in place, so every equation below uses the raised value:

```math
T_x \leftarrow \max\left(T_x, \epsilon_T\right),
\quad x \in \{1, 3\}
```

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\omega_r$ | `speed` | Input | [p.u.] | Machine rotor speed | Optional signal, defaults to synchronous
$P_\mathrm{ref}$ | `pref` | Input | [p.u.] | Governor reference | Optional signal, latched when unattached
$P_m$ | `pmech` | Output | [p.u.] | Mechanical-power signal seeded by the machine | Required signal

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$P_t$ | [p.u.] | Turbine-block output |
$P_v$ | [p.u.] | Valve position |

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$P_m$ | [p.u.] | Mechanical-power output | Read by the machine model

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units | Description | Note
------ | ----- | ----------- | ----
$\omega_r$ | [p.u.] | Machine rotor speed | Defaults to synchronous when unattached
$P_\mathrm{ref}$ | [p.u.] | Governor reference | Latched setpoint when unattached

## Model Equations

For readability, define the speed deviation and the valve target

```math
\Delta\omega = \omega_r - 1,
\qquad
g_v = -P_v + \dfrac{P_\mathrm{ref} - \Delta\omega}{R}.
```

### Internal Equations

#### Differential

The valve row uses the
[antiwindup](../../../../../CommonMath.md#antiwindup) target and smooth
approximation:

```math
\begin{aligned}
0 &= -\dfrac{\mathrm{d}P_v}{\mathrm{d}t}
     + \dfrac{1}{T_1}\,\mathrm{antiwindup}
       \left(P_v, g_v; P_v^{\min}, P_v^{\max}\right) \\
0 &= -\dfrac{\mathrm{d}P_t}{\mathrm{d}t}
     - \dfrac{P_t - P_v - T_2\dfrac{\mathrm{d}P_v}{\mathrm{d}t}}{T_3}
\end{aligned}
```

#### Algebraic

```math
0 = -P_m + P_t - D_t\,\Delta\omega
```

### External Equations

None.

### Wiring

The mechanical-power output $P_m$ is exposed on the `pmech` signal, which the
connected machine seeds during initialization.

## Initialization

The machine-provided $P_m$ is preserved and the steady state follows in
dependency order:

```math
\begin{aligned}
P_v &\leftarrow P_m + D_t\,\Delta\omega \\
P_t &\leftarrow P_v \\
P_\mathrm{ref} &\leftarrow \Delta\omega + R\,P_v
\end{aligned}
```

Initialization fails when the initial valve position falls outside
$[P_v^{\min}, P_v^{\max}]$.

## Monitors

None.
