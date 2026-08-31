# LoadZIP

Static ZIP load model with constant impedance, constant current, and constant
power fractions. `LoadZIP` owns terminal current states and adds their current
contribution to the connected bus residual.

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

Symbol | Units  | Description                              | Note
-------|--------|------------------------------------------|------
$I_r$  | [p.u.] | Terminal current, real component         | Added to connected bus residual
$I_i$  | [p.u.] | Terminal current, imaginary component    | Added to connected bus residual

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

```math
\begin{aligned}
0 &= I_r + (G V_r + B V_i)
\left[
\alpha_Z
+ \alpha_I \frac{V_\text{nom}}{V}
+ \alpha_P \frac{V_\text{nom}^2}{V^2}
\right] \\
0 &= I_i + (G V_i - B V_r)
\left[
\alpha_Z
+ \alpha_I \frac{V_\text{nom}}{V}
+ \alpha_P \frac{V_\text{nom}^2}{V^2}
\right]
\end{aligned}
```

### Bus Current-Balance Contributions

The connected bus receives $I_r$ and $I_i$ as positive contributions to its
real and imaginary current-balance residuals, respectively.

## Initialization

```math
\begin{aligned}
V_\text{nom} &\leftarrow \sqrt{V_r^2 + V_i^2} \\
I_r &\leftarrow -G V_r - B V_i \\
I_i &\leftarrow -G V_i + B V_r
\end{aligned}
```

The nominal-voltage anchor and derived admittance parameters are recomputed
from the initialized bus voltage. Initialization fails if the voltage magnitude
is not positive and finite. The derivative vector entries initialize to zero.

## Monitors

Name | Units  | Description                                  | Note
-----|--------|----------------------------------------------|------
`ir` | [p.u.] | Terminal current, real component             | Added to connected bus residual
`ii` | [p.u.] | Terminal current, imaginary component        | Added to connected bus residual
`im` | [p.u.] | Terminal current magnitude                   |
`p`  | [p.u.] | Active power at the connected bus terminal   | Positive for injection into the connected bus
`q`  | [p.u.] | Reactive power at the connected bus terminal | Positive for injection into the connected bus
