# LoadZIP

Static ZIP load model with constant impedance, constant current, and constant
power fractions. `LoadZIP` has no solver-owned variables; it computes terminal
current contributions from the connected bus voltage and adds them directly to
the bus current-balance residuals.

## Model Parameters

Symbol          | Units  | JSON     | Description                    | Typical Value
----------------|--------|----------|--------------------------------|--------------
$P_\text{nom}$  | [p.u.] | `Pnom`   | Nominal consumed real power    | 0.0
$Q_\text{nom}$  | [p.u.] | `Qnom`   | Nominal consumed reactive power | 0.0
$\alpha_I$      | [-]    | `alphaI` | Constant current load fraction | 0.0
$\alpha_P$      | [-]    | `alphaP` | Constant power load fraction   | 0.0

### Parameter Validation

None.

### Model Derived Parameters

$V_\text{nom}$ is the initial voltage magnitude of the respective bus.

```math
\begin{aligned}
G &= \frac{P_\text{nom}}{V_\text{nom}^2} \\
B &= \frac{Q_\text{nom}}{V_\text{nom}^2} \\
\alpha_Z &= 1 - \alpha_I - \alpha_P
\end{aligned}
```

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

Symbol | Units  | Description                              | Note
-------|--------|------------------------------------------|------
$V_r$  | [p.u.] | Terminal voltage, real component         | Owned by connected bus
$V_i$  | [p.u.] | Terminal voltage, imaginary component    | Owned by connected bus

## Wiring

Port  | Type | Description
------|------|------------
`bus` | Bus  | Connected bus that owns terminal voltage variables and current-balance residuals

## Model Equations

Let $V = \sqrt{V_r^2 + V_i^2}$.

### Differential Equations

None.

### Algebraic Equations

None.

### Bus Current-Balance Contributions

Let $I_r^{\mathrm{LoadZIP}}$ and $I_i^{\mathrm{LoadZIP}}$ denote the model
contributions to the real and imaginary current-balance residuals of the
connected bus. Positive current is oriented entering the bus.

```math
\begin{aligned}
I_r^{\mathrm{LoadZIP}}
  &= -(G V_r + B V_i)
\left[
\alpha_Z
+ \alpha_I \frac{V_\text{nom}}{V}
+ \alpha_P \frac{V_\text{nom}^2}{V^2}
\right] \\
I_i^{\mathrm{LoadZIP}}
  &= -(G V_i - B V_r)
\left[
\alpha_Z
+ \alpha_I \frac{V_\text{nom}}{V}
+ \alpha_P \frac{V_\text{nom}^2}{V^2}
\right]
\end{aligned}
```

These contributions are accumulated directly into the bus-owned residuals.

## Initialization

```math
V_\text{nom} \leftarrow \sqrt{V_r^2 + V_i^2}
```

The nominal-voltage anchor and derived admittance parameters are recomputed
from the initialized bus voltage. Initialization fails if the voltage magnitude
is not positive and finite. The model has no internal state to initialize.

## Monitors

Name | Units  | Description                                  | Note
-----|--------|----------------------------------------------|------
`ir` | [p.u.] | Terminal current, real component             | Added to connected bus residual
`ii` | [p.u.] | Terminal current, imaginary component        | Added to connected bus residual
`im` | [p.u.] | Terminal current magnitude                   |
`p`  | [p.u.] | Active power at the connected bus terminal   | Positive for injection into the connected bus
`q`  | [p.u.] | Reactive power at the connected bus terminal | Positive for injection into the connected bus
