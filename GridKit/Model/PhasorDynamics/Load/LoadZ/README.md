# LoadZ

Static constant-impedance load model. `LoadZ` owns terminal current states and
adds their current contribution to the connected bus residual.

## Model Parameters

Symbol | Units  | JSON | Description
-------|--------|------|------------
$R$    | [p.u.] | `R`  | Load resistance
$X$    | [p.u.] | `X`  | Load reactance

### Parameter Validation

None.

### Model Derived Parameters

```math
\begin{aligned}
G &= \frac{R}{R^2 + X^2} \\
B &= -\frac{X}{R^2 + X^2}
\end{aligned}
```

## Model Ports

Name  | Port | Init  | Description
------|------|-------|------------
`bus` | Bus  | Known | Connected bus that owns terminal voltage variables and current-balance residuals

## Model Variables

### Internal Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description                           | Note
-------|--------|---------------------------------------|------
$I_r$  | [p.u.] | Terminal current, real component      | Added to connected bus residual
$I_i$  | [p.u.] | Terminal current, imaginary component | Added to connected bus residual

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description                           | Note
-------|--------|---------------------------------------|------
$V_r$  | [p.u.] | Terminal voltage, real component      | Owned by connected bus
$V_i$  | [p.u.] | Terminal voltage, imaginary component | Owned by connected bus

## Model Equations

### Internal Equations

#### Differential

None.

#### Algebraic

```math
\begin{aligned}
0 &= I_r + G V_r - B V_i \\
0 &= I_i + B V_r + G V_i
\end{aligned}
```

### External Equations

```math
\begin{aligned}
I_r^{\mathrm{bus}} &\leftarrow I_r^{\mathrm{bus}} + I_r \\
I_i^{\mathrm{bus}} &\leftarrow I_i^{\mathrm{bus}} + I_i.
\end{aligned}
```

## Initialization

Initialization solves the algebraic current states from the connected bus
voltage. Let $V_{r0}$ and $V_{i0}$ be the initialized bus voltage components.

```math
\begin{aligned}
I_r &= -(G V_{r0} - B V_{i0}) \\
I_i &= -(B V_{r0} + G V_{i0})
\end{aligned}
```

The derivative vector entries initialize to zero.

## Monitors

Monitor | Units  | Description                                  | Note
--------|--------|----------------------------------------------|------
`p`  | [p.u.] | Active power at the connected bus terminal   | Positive for injection into the connected bus
`q`  | [p.u.] | Reactive power at the connected bus terminal | Positive for injection into the connected bus
