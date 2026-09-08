# SEXS-PTI Model

Simplified excitation system model ported from PhasorDynamics. The controller
equations and output anti-windup are unchanged; the EMT terminal port supplies
three-phase voltages in volts. An optional voltage-measurement lag defaults to
an algebraic bypass.

## Block Diagram

![SEXS-PTI model block diagram](../../../../../../docs/Figures/SEXS_PTI_DIAGRAM.png)

Figure 1: Exciter SEXS-PTI model. Figure courtesy of [PowerWorld](https://www.powerworld.com/WebHelp/)

## Model Parameters

Symbol | Units | JSON | Description | Note
------ | ----- | ---- | ----------- | ----
$V$ | [V] | `V` | Rated line-to-line RMS terminal voltage | Required
$T_r$ | [s] | `Tr` | Terminal-voltage measurement lag | Default $0$; zero bypasses
$T_A$ | [s] | `Ta` | Lead-lag numerator time constant | Required
$T_B$ | [s] | `Tb` | Lead-lag denominator time constant | Required
$T_E$ | [s] | `Te` | Exciter field time constant | Required
$K$ | [p.u.] | `K` | Voltage regulator gain | Required
$E_{fd}^{\max}$ | [p.u.] | `Efdmax` | Maximum excitation output | Required
$E_{fd}^{\min}$ | [p.u.] | `Efdmin` | Minimum excitation output | Required

PowerWorld/PSS/E SEXS_PTI data often gives $T_A/T_B$ as a ratio. GridKit stores
$T_A$ and $T_B$ separately, so convert ratio-format data with
$T_A = (T_A/T_B)T_B$ before passing parameters to the model.

The six controller parameters and `V` are required. `Tr` defaults to zero.
All parameters must be finite, `V` must be positive, and `Tr` non-negative.

### Parameter Validation

Invalid SEXS-PTI parameter sets are rejected by the following checks:

```math
\begin{aligned}
  T_A &\ge 0 \\
  T_B, T_E, K &> 0 \\
  E_{fd}^{\min} &< E_{fd}^{\max}
\end{aligned}
```

### Derived Parameters

None.

## Model Ports

Symbol | Port | Type | Units | Description | Note
------ | ---- | ---- | ----- | ----------- | ----
$\mathbf{v}$ | `bus` | Input | [V] | Three-phase terminal voltage | Required
$V_{\mathrm{ref}}$ | `vref` | Input | [p.u.] | Voltage reference | Inferred when unattached
$V_S$ | `vs` | Input | [p.u.] | Stabilizer input | Optional, defaults to zero
$V_{\mathrm{UEL}}$ | `vuel` | Input | [p.u.] | Under-excitation limiter input | Optional, defaults to zero
$V_{\mathrm{OEL}}$ | `voel` | Input | [p.u.] | Over-excitation limiter input | Optional, defaults to zero
$E_{fd}$ | `efd` | Output | [p.u.] | Field voltage | Seeded by the machine

## Submodels

None.

### Submodel Validation

None.

## Model Variables

### Internal Variables

#### Differential

Symbol    | Units  | Description                       | Note
----------|--------|-----------------------------------|-----
$V_R$     | [p.u.] | Lag-lead block state              |
$E_{fd}$  | [p.u.] | Exciter field voltage output      |
$E_C$ | [p.u.] | Measured terminal voltage | When $T_r>0$

#### Algebraic

Symbol    | Units  | Description                       | Note
----------|--------|-----------------------------------|-----
$V_{tr}$  | [p.u.] | Terminal voltage error signal     |
$E_C$ | [p.u.] | Measured terminal voltage | When $T_r=0$

### External Variables

#### Differential

None.

#### Algebraic

Symbol          | Units  | Description                                  | Note
----------------|--------|----------------------------------------------|-----
$v_a,v_b,v_c$   | [V]    | Instantaneous phase voltages                  | Bus input
$V_{\mathrm{ref}}$       | [p.u.] | Reference voltage                            | Signal port `vref`
$V_S$           | [p.u.] | Stabilizer output                            | Signal port `vs`
$V_{\mathrm{OEL}}$       | [p.u.] | Over-excitation limiter signal               | Signal port `voel`
$V_{\mathrm{UEL}}$       | [p.u.] | Under-excitation limiter signal              | Signal port `vuel`

## Model Equations

For readability, define the terminal-voltage magnitude and pre-limit derivative:

```math
V_t = \dfrac{\sqrt{v_a^2+v_b^2+v_c^2}}{V},
\qquad
f = \dfrac{1}{T_E}\left[-E_{fd} + \dfrac{K}{T_B}(-V_R + T_A V_{tr})\right].
```

Balanced phase voltages of rated line-to-line RMS voltage $V$ give $V_t=1$.
The field-voltage limit uses the CommonMath smooth
[antiwindup](../../../../../CommonMath.md#antiwindup) function.

### Internal Equations

#### Differential

```math
\begin{aligned}
0 &= -\dfrac{\mathrm{d}V_R}{\mathrm{d}t}
     - V_{tr} + \dfrac{-V_R + T_A V_{tr}}{T_B} \\
0 &= -\dfrac{\mathrm{d}E_{fd}}{\mathrm{d}t}
     + \mathrm{antiwindup}(E_{fd},f;E_{fd}^{\min},E_{fd}^{\max}) \\
0 &= -T_r\dfrac{\mathrm{d}E_C}{\mathrm{d}t} + V_t-E_C,
     \quad T_r>0
\end{aligned}
```

#### Algebraic

```math
\begin{aligned}
0 &= -V_{tr}-E_C+V_{\mathrm{ref}}+V_S+V_{\mathrm{OEL}}+V_{\mathrm{UEL}} \\
0 &= V_t-E_C,\quad T_r=0
\end{aligned}
```

### External Equations

None.

## Initialization

The generator initializes the EFD signal first. SEXS-PTI then reads that value
and any attached $V_S$, $V_{\mathrm{OEL}}$, and $V_{\mathrm{UEL}}$ signals and assumes steady state:

```math
\begin{aligned}
V_{tr} &\leftarrow \dfrac{E_{fd}}{K} \\
V_R &\leftarrow (T_A - T_B)V_{tr} \\
V_{\mathrm{ref}} &\leftarrow E_C + V_{tr} - V_S - V_{\mathrm{OEL}} - V_{\mathrm{UEL}}
\end{aligned}
```

The measured terminal voltage initializes to `V_t`. All derivatives initialize
to zero. An initial field voltage outside the configured limits is rejected.
Optional attached inputs are preserved and read live. The inferred reference
above is used only when `vref` is unattached. If a supplied reference differs,
consistent initialization resolves algebraic values and derivatives while
retaining the initialized differential states.

## Monitors

Monitor | Units | Description | Note
------- | ----- | ----------- | ----
`efd` | [p.u.] | Field-voltage output | $E_{fd}$
`vts` | [p.u.] | Measured terminal voltage | $E_C$
`vr` | [p.u.] | Lead-lag state | $V_R$
`vtr` | [p.u.] | Voltage error | $V_{tr}$
