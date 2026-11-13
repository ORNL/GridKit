# LoadZIP

Static ZIP load model with constant impedance, constant current, and constant
power fractions. `LoadZIP` owns terminal current states and adds their current
contribution to the connected bus residual.

## Model Parameters

Symbol           | Units  | JSON     | Description                     | Typical Value | Note
-----------------|--------|----------|---------------------------------|---------------|-----
$\alpha_I$       | [-]    | `alphaI` | Constant current load fraction  | 0.0           |
$\alpha_P$       | [-]    | `alphaP` | Constant power load fraction    | 0.0           |

### Parameter Validation

None.

### Model Derived Parameters

$V_\mathrm{nom}$ is the initial voltage magnitude of the respective bus. $P_0$
and $Q_0$ are the power the state terminal current consumes at that voltage.

```math
\begin{aligned}
P_0 &= -(V_\mathrm{r}I_\mathrm{r}+V_\mathrm{i}I_\mathrm{i}) \\
Q_0 &= -(V_\mathrm{i}I_\mathrm{r}-V_\mathrm{r}I_\mathrm{i}) \\
G &= \dfrac{P_0}{V_\mathrm{nom}^2} \\
B &= -\dfrac{Q_0}{V_\mathrm{nom}^2} \\
\alpha_Z &= 1 - \alpha_I - \alpha_P
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

## Model Equations

Let $V = \sqrt{V_r^2 + V_i^2}$.

### Internal Equations

#### Differential

None.

#### Algebraic

```math
\begin{aligned}
0 &= I_r + (G V_r - B V_i)
\left[
\alpha_Z
+ \alpha_I \dfrac{V_\mathrm{nom}}{V}
+ \alpha_P \dfrac{V_\mathrm{nom}^2}{V^2}
\right] \\
0 &= I_i + (G V_i + B V_r)
\left[
\alpha_Z
+ \alpha_I \dfrac{V_\mathrm{nom}}{V}
+ \alpha_P \dfrac{V_\mathrm{nom}^2}{V^2}
\right]
\end{aligned}
```

### External Equations

```math
\begin{aligned}
I_r^{\mathrm{bus}} &\leftarrow I_r^{\mathrm{bus}} + I_r \\
I_i^{\mathrm{bus}} &\leftarrow I_i^{\mathrm{bus}} + I_i
\end{aligned}
```

## Initialization

```math
\begin{aligned}
    I_r, I_i &\leftarrow \text{state terminal current on system base}
\end{aligned}
```

The derivative vector entries initialize to zero.

## Monitors

Monitor | Units  | Description                                  | Note
--------|--------|----------------------------------------------|------
`ir` | [p.u.] | Terminal current, real component             | Added to connected bus residual
`ii` | [p.u.] | Terminal current, imaginary component        | Added to connected bus residual
`im` | [p.u.] | Terminal current magnitude                   |
`p`  | [p.u.] | Active power at the connected bus terminal   | Positive for injection into the connected bus
`q`  | [p.u.] | Reactive power at the connected bus terminal | Positive for injection into the connected bus
