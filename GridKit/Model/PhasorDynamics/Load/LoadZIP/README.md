# LoadZIP

Static ZIP load model with constant impedance, constant current, and constant
power fractions. `LoadZIP` has no solver-owned variables; it computes terminal
current contributions from the connected bus voltage and adds them directly to
the bus current-balance residuals.

## Model Parameters

Symbol           | Units  | JSON     | Description                     | Typical Value | Note
-----------------|--------|----------|---------------------------------|---------------|-----
$P_\mathrm{nom}$ | [p.u.] | `Pnom`   | Nominal consumed real power     | 0.0           |
$Q_\mathrm{nom}$ | [p.u.] | `Qnom`   | Nominal consumed reactive power | 0.0           |
$\alpha_I$       | [-]    | `alphaI` | Constant current load fraction  | 0.0           |
$\alpha_P$       | [-]    | `alphaP` | Constant power load fraction    | 0.0           |

### Parameter Validation

None.

### Model Derived Parameters

$V_\mathrm{nom}$ is the initial voltage magnitude of the respective bus.

```math
\begin{aligned}
G &= \dfrac{P_\mathrm{nom}}{V_\mathrm{nom}^2} \\
B &= -\dfrac{Q_\mathrm{nom}}{V_\mathrm{nom}^2} \\
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

None.

### External Variables

#### Differential

None.

#### Algebraic

Symbol | Units  | Description                              | Note
-------|--------|------------------------------------------|------
$V_{\mathrm{r}}$ | [p.u.] | Terminal voltage, real component      | Owned by connected bus
$V_{\mathrm{i}}$ | [p.u.] | Terminal voltage, imaginary component | Owned by connected bus

## Model Equations

Let $V = \sqrt{V_{\mathrm{r}}^2 + V_{\mathrm{i}}^2}$.

### Internal Equations

#### Differential

None.

#### Algebraic

None.

### External Equations

Let $I_{\mathrm{r}}^{\mathrm{LoadZIP}}$ and
$I_{\mathrm{i}}^{\mathrm{LoadZIP}}$ denote the model contributions to the real
and imaginary current-balance residuals of the connected bus. Positive current
is oriented entering the bus.

```math
\begin{aligned}
I_{\mathrm{r}}^{\mathrm{LoadZIP}}
  &= -(G V_{\mathrm{r}} - B V_{\mathrm{i}})
\left[
\alpha_Z
+ \alpha_I \dfrac{V_\mathrm{nom}}{V}
+ \alpha_P \dfrac{V_\mathrm{nom}^2}{V^2}
\right] \\
I_{\mathrm{i}}^{\mathrm{LoadZIP}}
  &= -(G V_{\mathrm{i}} + B V_{\mathrm{r}})
\left[
\alpha_Z
+ \alpha_I \dfrac{V_\mathrm{nom}}{V}
+ \alpha_P \dfrac{V_\mathrm{nom}^2}{V^2}
\right]
\end{aligned}
```

```math
\begin{aligned}
I_r^{\mathrm{bus}} &\leftarrow I_r^{\mathrm{bus}} + I_{\mathrm{r}}^{\mathrm{LoadZIP}} \\
I_i^{\mathrm{bus}} &\leftarrow I_i^{\mathrm{bus}} + I_{\mathrm{i}}^{\mathrm{LoadZIP}}
\end{aligned}
```

## Initialization

```math
V_\mathrm{nom} \leftarrow \sqrt{V_\mathrm{r}^2 + V_\mathrm{i}^2}
```

The nominal-voltage anchor and derived admittance parameters are recomputed
from the initialized bus voltage. The model has no internal state to initialize.

## Monitors

Monitor | Units  | Description                                  | Note
--------|--------|----------------------------------------------|------
`ir` | [p.u.] | Terminal current, real component             | Added to connected bus residual
`ii` | [p.u.] | Terminal current, imaginary component        | Added to connected bus residual
`im` | [p.u.] | Terminal current magnitude                   |
`p`  | [p.u.] | Active power at the connected bus terminal   | Positive for injection into the connected bus
`q`  | [p.u.] | Reactive power at the connected bus terminal | Positive for injection into the connected bus
