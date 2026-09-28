# DependentNorton

Controlled current source with a parallel admittance. Terminal current is
positive for injection into the connected bus. The source-current inputs are
supplied by an external control law or coupling interface.

All quantities use the system power base and connected-bus voltage base.
Source-current components use the same phasor reference as the bus voltage.

## Model Parameters

Symbol | Units  | JSON | Description          | Typical Value | Note
-------|--------|------|----------------------|---------------|-----
$G$    | [p.u.] | `G`  | Parallel conductance |               | Required
$B$    | [p.u.] | `B`  | Parallel susceptance |               | Required

### Parameter Validation

Both parameters must be finite. Zero admittance is permitted.

### Model Derived Parameters

None.

## Model Ports

Name  | Port  | Init  | Description
------|-------|-------|------------
`bus` | Bus   | Known | Connected bus that owns terminal voltage variables and current-balance residuals
`inr` | Input | Known | Required Norton source current, real component
`ini` | Input | Known | Required Norton source current, imaginary component

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol           | Units  | Description                                | Note
-----------------|--------|--------------------------------------------|-----
$V_r$            | [p.u.] | Terminal voltage, real component           | Owned by connected bus
$V_i$            | [p.u.] | Terminal voltage, imaginary component      | Owned by connected bus
$I_r^\mathrm{N}$ | [p.u.] | Norton source current, real component      | Signal port `inr`
$I_i^\mathrm{N}$ | [p.u.] | Norton source current, imaginary component | Signal port `ini`

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

```math
\begin{aligned}
I_r^{\mathrm{bus}} &\leftarrow I_r^{\mathrm{bus}} + I_r^\mathrm{N} - G V_r + B V_i \\
I_i^{\mathrm{bus}} &\leftarrow I_i^{\mathrm{bus}} + I_i^\mathrm{N} - B V_r - G V_i
\end{aligned}
```

## Initialization

Both source-current inputs must be connected and supplied before initialization.
The supplied currents must be consistent with the network operating point.

There are no internal variables to initialize.

## Monitors

Monitor | Units  | Description                           | Note
--------|--------|---------------------------------------|-----
`ir`    | [p.u.] | Terminal current, real component      | $I_r^\mathrm{N} - G V_r + B V_i$
`ii`    | [p.u.] | Terminal current, imaginary component | $I_i^\mathrm{N} - B V_r - G V_i$
`p`     | [p.u.] | Terminal active power                 | $V_r I_r^\mathrm{N} + V_i I_i^\mathrm{N} - G(V_r^2 + V_i^2)$
`q`     | [p.u.] | Terminal reactive power               | $V_i I_r^\mathrm{N} - V_r I_i^\mathrm{N} + B(V_r^2 + V_i^2)$

Current and power are positive for injection into the connected bus.
